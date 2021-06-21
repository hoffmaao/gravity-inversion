function m = mmodtwodpoly(m0,varargin)

% m = mmodtwodpoly(m0,varargin)
% This function modifies an earth-model profile into polygons that are
% extended by +- 50 km of the either ends. This is usually in preparation 
% for 2D graviy-anomlay calculations.
% 
% Input
%   m0 = a matrix with the distance in the first column and elevations in
%        the rest
%   varargin: if any value is entered, the model will be flipped about the
%             datum
% 
% Output
%   m = a matrix of polygons

% Written by Atsuhiro Muto
% Dept. of Earth & Environ. Sci., Temple Univ.
% amuto@temple.edu
% Last updated July 6, 2016

a = [m0(1,1)-50000 m0(1,2:end);m0;m0(end,1)+50000 m0(end,2:end)];
m = [a(:,1);flipud(a(:,1));a(1,1)];

S = size(m0);

for i=2:S(2)-1
    b = [a(:,i);flipud(a(:,i+1));a(1,i)];
    m = [m b];
    if isempty(varargin)==0
        m(m(:,i)>0,i) = varargin{1}-m(m(:,i)>0,i);
        m(m(:,i)<0,i) = m(m(:,i)<0,i)-varargin{1};
    end
end
