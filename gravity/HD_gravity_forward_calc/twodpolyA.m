function Dg = twodpolyA(xobs,M,rho0,rhoback)

% TWODPOLYA computes gravity anomaly Dg of polygons using Talwani line 
% integral. For details, see Talwani et al. (1959), JGR, 64(1), 49-59.
% DIRECTION OF THE LOOP IS CLOCKWISE!
% 
% Input
%   xobs = n by 1 vector of observation points (m)
%   M = L by M matrix of polygon model, N = number of polygons
%       M(:,1) = x-positions of vertices (m)
%       M(2,:) = z-positions of vertices (m)
%   rho0 = densities of polygons (kg/m^3)
%   rhoback = background density (kg/m^3)

% 
% Output
%   Dg = n by 1 vector of gravity anomaly (mgal)

% Written by Atsuhiro Muto
% Dept. of Geosciences, Penn State Univ.
% aum34@psu.edu
% Last updated May. 8, 2014

% tic

N = length(rho0);	% number of polygons
mx = length(xobs);
drho = rho0-rhoback;        % density contrast
dg = nan(mx,N);      % initiate matrix for gravity anomaly
       
% -------------------------------------------------------------------------
% loop over polygons   
for ip=1:N
    
    id = find(M(2:end-1,ip+1)==M(1:end-2,ip+1)&...
              M(2:end-1,ip+1)==M(3:end,ip+1));
    
    xp = M(:,1);                % x-positions
    xp(id+1) = [];
    
    zp = M(:,ip+1);
    zp(id+1) = []; 
    zp(zp==0)=1e-10;            % if depth is zero, substitute with 1e-10 m
    clear mex
    dg(:,ip) = talwani2(xobs,xp,zp,drho(ip));
    
end
% -------------------------------------------------------------------------

% Sum over contributions from different polygons (dg) to get total gravity
% anomaly (Dg)
Dg = sum(dg,2)*1e5;

% toc
