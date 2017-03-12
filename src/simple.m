## Reseni N-S rovnic na oblasti tvaru ctverce s okrajovou podminkou
## [u,v] = 0 na leve, spodni a prave strane, [u,v]=[1,0] na horni stene
##
## K reseni je pouzit schema typu upwind kombinovane s metodou SIMPLE
##

global N = 50;      # Velikost site 
global h = 1./N;   # Krok site
global nu=0.05;     # Viscosity

## Pocatecni podminka
u = zeros(N,N);
v = zeros(N,N);
p = zeros(N,N); 

for iter=1:1000

  #[A,b] = stokes(u,v,p,nu,h);
  [A,b] = navierstokes(u,v,p,nu,h);
  
  ## Vypocet odhadu rychlosti
  uu = reshape( pcg(A,b(:,1)), size(u) );
  vv = reshape( pcg(A,b(:,2)), size(u) );
  
  d = div(uu,vv,h);
  
  ## Vypocet korekce tlaku
  [F,g] = pEqn(A,d,h);
  pp = reshape( pcg(F,g), size(u) );
  
  ## Vypocet korekce rychlosti
  up(N,N)=0;
  for j=1:N
    for i=2:N-1
      ip = i + (j-1)*N;
      ap = A(ip,ip);
      up(i,j) = - (pp(i+1,j)-pp(i-1,j))/2/h / ap;
    endfor
  endfor
  
  vp(N,N)=0;
  for j=2:N-1
    for i=1:N
      ip = i + (j-1)*N;
      ap = A(ip,ip);
      vp(i,j) = - (pp(i,j+1)-pp(i,j-1))/2/h / ap;
    endfor
  endfor
  
  
  beta_u = 0.7;
  beta_p = 0.3;
  
  un = uu + beta_u*up;
  vn = vv + beta_u*vp;
  pn = p + beta_p*pp;
  
  printf("Iterace %i, norma divergence %f\n", iter, h*norm(div(u,v,h)) )
  printf("\tnormy rezidui %f, %f, %f\n", 
	 norm(un-u)*h, norm(vn-v)*h, norm(pn-p)*h
	 )

  u = un;
  v = vn;
  p = pn - sum(sum(pn))/N**2;

endfor


## ulozeni vysledku do VTK souboru 
f = fopen("simple.vtk", "w")
fprintf(f,"# vtk DataFile Version 2.0\n")
fprintf(f,"Driven cavity problem\n")
fprintf(f,"ASCII\n\n")
fprintf(f,"DATASET STRUCTURED_POINTS\n")
fprintf(f,"DIMENSIONS %i %i 1\n", N+1, N+1)
fprintf(f,"ORIGIN 0 0 0\n")
fprintf(f,"SPACING %f %f %f\n\n", h, h, h)

fprintf(f,"CELL_DATA %i\n\n", N*N)

fprintf(f,"VECTORS U double\n")
for j=1:N
  for i=1:N
    fprintf(f,"%f %f 0\n", u(i,j), v(i,j))
  endfor
endfor

fprintf(f,"SCALARS p double\n")
fprintf(f,"LOOKUP_TABLE default\n")
for j=1:N
  for i=1:N
    fprintf(f,"%f\n", p(i,j))
  endfor
endfor

fprintf(f,"CELL_DATA n\n")

fclose(f)

