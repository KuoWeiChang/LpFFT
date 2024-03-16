clear all;
close all;
clc;



t=0;
e=0;
K=1;
for m=1:K

m=6;
N=2^m;
x=double(rand(N,1)>0.5);
fx=fft(x);
mask64=zeros(64,1);
mask64([0,1,2,4,8,16,32]+1)=1;
if N>64
  mask=repmat(mask64,1,N/64)';
else
  mask=mask64;
end
fx=fx.*mask;






start = clock();
s=LPFFTw(m,fx);
elapsedTime = etime(clock(), start);

t=t+elapsedTime;
e=e+max(abs(x-s));
end
disp(sprintf("error = %9.6f, elapsed Time = %9.6f ",e/K, elapsedTime) );
