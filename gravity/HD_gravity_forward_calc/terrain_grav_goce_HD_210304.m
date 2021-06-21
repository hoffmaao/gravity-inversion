%% 
% This code calculates the long-wavelength Bouguer anomaly from GOCE 
% satellite gravity data. It first calculates the 3D gravity effect of the
% ice and bed topography, applies a low-pass filter to it and subtracts it
% from the GOCE gravity field. 

%%

dx = 5000;                          % prism dimension in x-direction
dy = 5000;                          % prism dimension in y-direction

% xx1 = -680000:dx:0;
% yy1 = -300000:dy:0;
xx1 = -850000:dx:150000;
yy1 = -450000:dy:150000;

xx2 = -850000:dx:150000;
yy2 = -450000:dy:150000;

% define grids
nx = length(xx1);                   % number of data points in x-direction
ny = length(yy1);                   % number of data points in y-direction

mx = length(xx2);                   % number of prisms in x-direction
my = length(yy2);                   % number of prisms in y-direction

% 
%%

load BM_5km
A = BM_5km;

% define density
rhoi = 900;
rhoc = 2670;
rho1 = [rhoi;rhoc];                 % density contrast between bedrock (2670 kg/m^3)
rho2 = rhoi-rhoc;                   % density contrast between bedrock (2670 kg/m^3)

dn = 150000;                        % distance for neighbouring prisms for gravity calculation

Pobs = makep(nx,ny,dx,dy)+dx/2;     % data positions
Pp = makep(mx,my,dx,dy); 
Pp = [Pp(:,1) Pp(:,1)+dx Pp(:,2) Pp(:,2)+dy];
xobs = vec2grid(Pobs(:,1),ny,nx);
yobs = vec2grid(Pobs(:,2),ny,nx);
distxy = [reshape(xobs',ny*nx,1) reshape(yobs',ny*nx,1)];

% mean flight level
flevel = 3800;

[m01,m02] = msplit(-A(:,3:4),flevel);
m01 = [m01 ones(my*mx,1)*flevel+1e-3];

%%

% clear memory carried over from previous mex calls
clear mex

tic
g = plouff3(Pobs,Pp,m01,m02,rho1,rho2,dx,dy,dn);
toc

terrgrav_HD_5000 = vec2grid(g,ny,nx)*1e5;

%% low-pass filter

xp = xx1;
yp = yy1;
nx = length(xp);
ny = length(yp);

maxx = max(xp);
minx = min(xp);
maxy = max(yp);
miny = min(yp);

dx = abs(maxx-minx)/(nx-1); 
dy = abs(maxy-miny)/(ny-1);
dfx = 1/(dx*(nx-1)); 
dfy = 1/(dy*(ny-1));

clear fkx fky fkx_pom fky_pom
R = rem(nx,2);
if R==0
       fkx(1:nx/2+1) = 0:dfx:(nx/2)*dfx;
       fkx_pom = -fkx(2:nx/2+1);
       fkx(nx/2+1:nx) = fkx_pom(:,[length(fkx_pom):-1:1]);
else
       fkx(1:round(nx/2))=0:dfx:floor(nx/2)*dfx;
       fkx_pom = -fkx(2:round(nx/2));
       fkx(round(nx/2)+1:nx) = fkx_pom(:,[length(fkx_pom):-1:1]);
end

R = rem(ny,2);
if R==0
       fky(1:ny/2+1) = 0:dfy:(ny/2)*dfy;
       fky_pom = -fky(2:ny/2+1);
       fky(ny/2+1:ny) = fky_pom(:,[length(fky_pom):-1:1]);
else
       fky(1:round(ny/2))=0:dfy:floor(ny/2)*dfy;
       fky_pom = -fky(2:round(ny/2));
       fky(round(ny/2)+1:ny) = fky_pom(:,[length(fky_pom):-1:1]);
end

K = ((repmat(fkx,ny,1).^2)+(repmat(fky',1,nx).^2)).^0.5;
KK = 1./K/1000;

%%

H = nan(size(KK));
for i=1:length(KK)
    c = 110-KK(:,i);
    H(:,i) = 0.5*(1+cos(c*2*pi/(2*20)));
    id = find(c<0);
    H(id,i) = 1;
    id = find(c>20);
    H(id,i) = 0;
end   

Fdata = fft2(terrgrav_HD_5000);
terrgrav_lpf = real(ifft2(H.*Fdata));

% GOCE data
load goce_HD_5km

boug_goce5000 = goce_HD_5km(:,:,3)-terrgrav_lpf;

% save boug_goce5000 boug_goce5000