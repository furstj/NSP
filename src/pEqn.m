## Diskretizace "Poissonovy" rovnice pro tlak

function [F,g] = pEqn(A, d, h)
  [N,M] = size(d);
  NM = N*M;
  
  F = sparse(NM,NM);

  for j=1:M
    for i=1:N
      ip = i + N*(j-1);
      in = ip + N;
      ie = ip + 1;
      is = ip - N;
      iw = ip - 1;

      if (j<M)
	an = (A(ip,ip)+A(in,in))/2;
	F(ip,in)  = - 1/h**2/an;
	F(ip,ip) += 1/h**2/an;
      endif

      if (i<N)
	ae = (A(ip,ip)+A(ie,ie))/2;
	F(ip,ie)  = - 1/h**2/ae;
	F(ip,ip) += 1/h**2/ae;
      endif

      if (j>1)
	as = (A(ip,ip)+A(is,is))/2;
	F(ip,is)  = - 1/h**2/as;
	F(ip,ip) += 1/h**2/as;
      endif

      if (i>1)
	aw = (A(ip,ip)+A(iw,iw))/2;
	F(ip,iw)  = - 1/h**2/aw;
	F(ip,ip) += 1/h**2/aw;
      endif

    endfor
  endfor

  g = -vec(d);

  ## Fixuji tlak tak, ze prictu k soustave rovnici
  ## p(1,1) = 0
  F(1,1) += 1/h**2;

endfunction