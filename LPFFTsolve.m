function y=LPFFTsolve(upb,N,fgt)
  ## integer linear programming to solve binary ifft 
  # upb :previous sum vector,example:[5,3,6,7]
  # fft length = 2*length(upb)
  # fgt: sampled fft parameter
  # N total length N/2^k is section length

  
  
  upb=real(upb);
  L = length(upb)*2;
  Len = N/L;
  F=fft(eye(L));

  LF=[real(F(2,:));imag(F(2,:))];
  sLF=LF(:,1:L/2)-LF(:,L/2+1:L);
  v=[real(fgt(Len+1));imag(fgt(Len+1))]-LF(:,L/2+1:L)*upb;

  
  ## LP part
  c = ones(L/2,1);
  A=sLF;
  b = v;
  lb = max(0,upb-Len); # binary
  ub = min(upb,Len);   #binary
  #lb=0*upb; #not binary
  #ub=upb; #not binary
  ## disp(prod(ub-lb+1))
  ctype = repmat('S',1,length(LF(:,1)));
  vartype = repmat('I',1,L/2);
  s = -1;

  param.msglev = 1;
  param.itlim = 100;


  [xmin, fmin, status, extra] = ...
   glpk (c, A, b, lb, ub, ctype, vartype, s, param);
   y=[xmin;upb-xmin];