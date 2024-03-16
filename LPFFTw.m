function y=LPFFTw(t,fgt)
  N=2^t;

  if N>2^6
    fgt1=(fgt(1:N/2)+fgt(N/2+1:N))/2;
    fgt2=(fgt(1:N/2)-fgt(N/2+1:N)).*exp(2*pi*1i/N*[0:N/2-1]')/2;
    y1=LPFFTw(t-1,fgt1);
    y2=LPFFTw(t-1,fgt2);
    y=reshape([y1,y2]',N,1);
  else
    s=LPFFTsolve(fgt(1),N,fgt);
    for k=2:t
      s=LPFFTsolve(s,N,fgt);
    end
    y=s;
  endif