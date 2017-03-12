function d=div(u,v,h)
  [N,M] = size(u);
  d(N,M)=0;

  for j=1:M
    for i=1:N

      if (i>1)
	ul = u(i-1,j);
      else
	ul = 0;
      endif

      if (i<N)
	ur = u(i+1,j);
      else
	ur = 0;
      endif

      if (j>1)
	vl = v(i,j-1);
      else
	vl = 0;
      endif
      
      if (j<M)
	vr = v(i,j+1);
      else
	vr = 0;
      endif

      d(i,j) = (ur-ul)/2/h + (vr-vl)/2/h;
    endfor
  endfor

endfunction