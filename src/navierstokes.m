## Diskretizace Navierovych-Stokesovych problemu 
##
## s karteszkou siti a krokem h na obdelnikove oblasti
##
## Okrajove podminky jsou [u,v]=[0,0] na leve, prave a dolni strane a
## [u,v] = [1,0] na horni strane
## 
## Vysledkem je matice A a dvojice pravych stran, tak, ze
## A*u = b(:,1) a A*v = b(:,2) 

function [A,b] = navierstokes(u,v,p,nu,h)
  [N,M] = size(p);
  NM = N*M;
  
  ## Vyuziji diskretizaci stokesova problemu
  [A, b] = stokes(u,v,p,nu,h);   

  ## Pridam konvektivni cleny pomoci upwindu
  for j=1:M
    for i=1:N
      ip = i + N*(j-1);
      in = ip + N;
      ie = ip + 1;
      is = ip - N;
      iw = ip - 1;

      if (i>1)
	uw= (u(i-1,j) + u(i,j)) / 2;
	if (uw>0)
	  A(ip,iw) -= uw/h;
	else
	  A(ip,ip) -= uw/h;
	endif
      endif

      if (i<N)
	ue = (u(i,j) + u(i+1,j)) / 2;
	if (ue>=0)
	  A(ip,ip) += ue/h;
	else
	  A(ip,ie) += ue/h;
	endif
      endif

      if (j>1)
	vs = (v(i,j-1) + v(i,j)) / 2;
	if (vs>0)
	  A(ip,is) -= vs/h;
	else
	  A(ip,ip) -= vs/h;
	endif
      endif

      if (j<M)
	vn = (v(i,j) + v(i,j+1)) / 2;
	if (vn>=0)
	  A(ip,ip) += vn/h;
	else
	  A(ip,in) += vn/h;
	endif
      endif

    endfor
  endfor

endfunction
