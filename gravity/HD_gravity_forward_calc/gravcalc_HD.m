function gravcalc_HD(LN,dm)

% Gravity anomaly is referenced to the geoid whereas disturbance is 
% referenced to the ellipsoid. This code uses the observed gravity 
% disturbance because the elevation data are all referenced to WGS1984
% ellipsoid. 
% 
% This code takes about two minutes to run.
% 
% Input
%   LN = a character string specifying a PolarGAP flight line, P15, P16,
%        P17 or P18.
%   dm = a value of the gravity-modeling intervals in meters, 2000 m is
%        enough and the calculation is fast
% 
% Output
%   A figure of the ice-surface and the bed topography and gravity disturbances

% Written by Atsuhiro Muto
% Dept. of Earth & Environ. Sci., Temple Univ.
% amuto@temple.edu
% Last updated March 1, 2021

% tic
%% Load data 

load polargap_HD

A = polargap_HD.(LN);
xobs = round(A(:,10));

% mean level of data acquisition, mean of flight elevation
flevel = round(mean(A(:,4)));

%% Define layer densities

rhoi = 900;             % ice
rhoback = 2670; %2670        % background (bedrock, assumes standard crystalline rock)

% Convert to density contrast, which is what the forward gravity
% calculation takes.
rho = [rhoi+rhoback;rhoback*2;rhoi]; 

%% Ice-surface and bed topography

m0 = [xobs -A(:,13) -A(:,11)];
id = isnan(m0(:,2));
m0(id,:) = [];

% Split into what is above and below the ellipsoid and covert to distance
% from the observation level
[m1,m2] = msplit(m0(:,2:end),flevel);

% Add the ellipsoid level
m1 = [m1 ones(length(m1),1)*flevel+1e-3];

% Extend the line +- 50 km and modify into polygons
m11 = mmodtwodpoly([m0(:,1) m1]);
m21 = mmodtwodpoly([m0(:,1) m2]);
mm = [m11 m21(:,2:end)];

%% Gravity calculations

% Extend the observation by 20 km on either end
a = (min(xobs)-20000:dm:-dm)';
b = (max(xobs)+dm:dm:max(xobs)+20000)';

xobsb = (xobs(1):dm:xobs(end))';
xobs2=[a;xobsb;b];
id1 = find(xobs2==min(xobsb));
id2 = find(xobs2==max(xobsb));

% clear memory carried over from previous mex calls
clear mex

% Gravity calculation
g = twodpolyA(xobs2,mm,rho,rhoback);

% Low-pass filter
hmax = 20000;
hmin = 10000;
glpf = fftfilt_1d(xobs2,g,hmax,hmin)';

%% Long-wavelength signal from GOCE

load boug_goce5000

S = size(boug_goce5000);
XI = reshape(boug_goce5000(:,:,1)',S(1)*S(2),1);
YI = reshape(boug_goce5000(:,:,2)',S(1)*S(2),1);
ZI = reshape(boug_goce5000(:,:,3)',S(1)*S(2),1);

F = scatteredInterpolant(XI,YI,ZI,'natural');
z = F(A(:,8),A(:,9));
boug_goce = interp1(xobs,z,xobsb);

gmodel = glpf(id1:id2)+boug_goce;

%%
% Plot density model
mm2 = [mm(:,1) mm(:,2:end)-flevel];
col = ['c','w','c','c','c'];
figure
subplot(211);
title(LN)

npolygon = size(mm2,2);
for ip=1:npolygon-1
    zp = mm2(:,ip+1);          % z-position
    fill(mm2(:,1)/1000,-1*zp,col(ip)); 
    hold on
end
title(LN)

xlim([min(xobs)/1000 max(xobs)/1000])
ylabel('Elevation (m)')
grid on

% Gravity anomaly
gobs = interp1(xobs,A(:,5),xobsb);
subplot(212)
plot(xobsb/1000,gobs,'ko');
hold all
plot(xobsb/1000,glpf(id1:id2),'-','LineWidth',2,'Color',[0 0.447 0.741]);
plot(xobsb/1000,boug_goce,'-','LineWidth',2,'Color',[0.466 0.674 0.188]);
plot(xobsb/1000,gmodel,'-','LineWidth',2,'Color',[0.850 0.325 0.098]);
title(LN)

xlim([min(xobs)/1000 max(xobs)/1000])
xlabel('Distance along the profile (km)')
ylabel('Gravity anomaly (miligals)')
grid on
legend('Observed gravity','Gravity effect of ice and topography',...
       'Gravity effect of the Moho','Sum of modelled gravity')

% toc