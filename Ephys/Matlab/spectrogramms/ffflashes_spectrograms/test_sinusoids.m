%test
Fs=1000;
x=linspace(1,10*pi,5*Fs+1);
npts=length(x);

nnoise=randi(1000,1,npts)/1000/5;

y1=sin(2*x)+nnoise;
y2=sin(5*x)+nnoise;


N = Fs/2*1000;
f = 0:(Fs/N):(Fs/2);

% FFT
X = fft(y1,N);
%take one side only
X = X(1:N/2+1);
% power normalized by freq resolution
psp = (1/(Fs*npts)) * abs(X).^2;
% adjust because one-sided
psp(2:end-1) = 2*psp(2:end-1);

figure, 
subplot(2,3,[1:3]); plot(x,y1);
subplot(2,3,[4:6]); loglog(f,psp);


y2=sin(5*x)+nnoise;
% FFT
X = fft(y2,N);
%take one side only
X = X(1:N/2+1);
% power normalized by freq resolution
psp = (1/(Fs*npts)) * abs(X).^2;
% adjust because one-sided
psp(2:end-1) = 2*psp(2:end-1);

figure, 
subplot(2,3,[1:3]); plot(x,y2);
subplot(2,3,[4:6]); loglog(f,psp);



% FFT
X = fft(10+y1+y2,N);
%take one side only
X = X(1:N/2+1);
% power normalized by freq resolution
psp = (1/(Fs*npts)) * abs(X).^2;
% adjust because one-sided
psp(2:end-1) = 2*psp(2:end-1);

figure, 
subplot(2,3,[1:3]); plot(x,10+y1+y2);
subplot(2,3,[4:6]); loglog(f,psp);




% FFT
X = fft(y1+y2,N);
%take one side only
X = X(1:N/2+1);
% power normalized by freq resolution
psp = (1/(Fs*npts)) * abs(X).^2;
% adjust because one-sided
psp(2:end-1) = 2*psp(2:end-1);

figure, 
subplot(2,3,[1:3]); plot(x,y1+y2);
subplot(2,3,[4:6]); loglog(f,psp);


% FFT
X = fft(100*y1+y2,N);
%take one side only
X = X(1:N/2+1);
% power normalized by freq resolution
psp = (1/(Fs*npts)) * abs(X).^2;
% adjust because one-sided
psp(2:end-1) = 2*psp(2:end-1);

figure, 
subplot(2,3,[1:3]); plot(x,100*y1+y2);
subplot(2,3,[4:6]); loglog(f,psp);





