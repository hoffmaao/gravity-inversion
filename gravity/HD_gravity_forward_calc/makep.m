function P = makep(nx,ny,dx,dy)

% P = makep(nx,ny,dx,dy)
% MAKEP creates an N by 2 matrix of x and y positions of grids from
% specified number of points and spacings in x- and y-directions.
% 
% Input
%   nx = number of points in x-direction
%   ny = number of points in y-direction
%   dx = spacing in x-direction
%   dy = spacing in y-direction
%
% Output
%   P = N by 2 matrix of x- (1st column) and y- (2nd column) positions

% Written by Atsuhiro Muto
% Dept. of Geosciences, Penn State Univ.
% aum34@psu.edu, atsumuto@gmail.com
% Last updated Jun. 18, 2011

x = (0:dx:(nx-1)*dx)';
P = repmat(x,ny,1);

count = 0;
for i=1:nx:length(P)
    P(i:i+nx-1,2) = ones(nx,1)*count*dy;
    count = count+1;
end