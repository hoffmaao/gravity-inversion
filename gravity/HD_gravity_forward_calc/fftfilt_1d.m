function fdata = fftfilt_1d(x,d,hmax,hmin)

% tic
xx = x';
dx = xx(2)-xx(1);
dd = d';

Fs = 1/dx;% sampling frequency
L = length(xx);
NFFT = 2^nextpow2(L);               % Next power of 2 from length of y
Y = fft(dd,NFFT);

%% low-pass filtering

hl = hmax-hmin;

f = Fs/2*linspace(0,1,NFFT/2+1);
ff = 1./f;
c = hmax-ff;
h = 0.5*(1+cos(c*2*pi/(hl*2)));
id = c<0;
h(id) = 1;
id = c>hl;
h(id) = 0;
H = [h(1:end-1) fliplr(h(1:end-1))];
fdata = real(ifft(H.*Y));
fdata = fdata(1:size(dd,2));
% toc