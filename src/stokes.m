## Diskretizace Stokesova problemu 
##
##  - nu*Laplace u = - p_x
##  - nu*Laplace v = - p_y
##
## s karteszkou siti a krokem h na obdelnikove oblasti
##
## Okrajove podminky jsou [u,v]=[0,0] na leve, prave a dolni strane a
## [u,v] = [1,0] na horni strane
## 
## Vysledkem je matice A a dvojice pravych stran, tak, ze
## A*u = b(:,1) a A*v = b(:,2) 

function [A,b] = stokes(u,v,p,nu,h)
  [N,M] = size(p);
  NM = N*M;
  A = sparse(NM,NM);
  b(NM,2) = 0;

  ## Diskretizace -nu*Laplace
  for j=1:M
    for i=1:N
      ip = i + N*(j-1);
      in = ip + N;
      ie = ip + 1;
      is = ip - N;
      iw = ip - 1;

      A(ip,ip) = 4*nu/h**2;

      if (j<M)
	A(ip,in) = - nu/h**2;
      else
	b(ip,1) = 1*nu/h**2;
      endif

      if (i<N)
	A(ip,ie) = -nu/h**2;
      endif

      if (j>1)
	A(ip,is) = -nu/h**2;
      endif

      if (i>1)
	A(ip,iw) = -nu/h**2;
      endif

    endfor
  endfor

  ## Diskretizace -p_x s nulovou derivaci na okraji
  for j=1:M
    for i=2:N-1
      ip = i + N*(j-1);
      b(ip,1) += - (p(i+1,j)-p(i-1,j))/2;
    endfor
  endfor

  ## Diskretizace -p_x s nulovou derivaci na okraji
  for j=2:M-1
    for i=1:N
      ip = i + N*(j-1);
      b(ip,2) += - (p(i,j+1)-p(i,j-1))/2;
    endfor
  endfor

endfunction