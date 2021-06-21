function [m1,m2] = msplit(m,varargin)

% [m1,m2] = msplit(m,varargin)
% This function is used to split an earth model to what is above and below
% the datum (usually the WGS84 ellipsoid).
% 
% Input
%   m = matrix of the earth model
%   varargin = a value of the elevation above datum of the data collection,
%       assumed to be 0 m if left blank
%
% Output
%   m1 = matrix of model parts above the datum
%   m2 = matrix of model parts below the datum

% Written by Atsuhiro Muto
% Dept. of Geosciences, Penn State Univ.
% aum34@psu.edu
% Last updated Jan. 6, 2014

if nargin>1
    dlevel = varargin{1};
else
    dlevel = 0;
end

S = size(m);
L = S(2);
% m1 = ones(S(1),L+1)*1e-10;
m1 = ones(S(1),L)*1e-10;
% m1(:,end) = 0;
m2 = m;

for i=1:L
   m1(:,i) = m1(:,i)*(L+1-i); 
   id1 = m(:,i)<0;
   m1(id1,i) = m(id1,i);
   
   m2(id1,i) = i*1e-10;
end
  
m1 = m1+dlevel;
m2 = m2+dlevel;
