
LN = 'P16';
load('polargap_HD')
data = polargap_HD.(LN);

% Add a column for sediments. Gravity calculation code cannot take 0 for a
% layer thickness so put in a small number (1e-3) if there are no
% sediments.
data = [data data(:,11)-1e-3];

%% Add sediments to desired intervals along the profile. 
% Specify intervals in meters along the profile. 

% Specify the interval
x1 = 400000;
x2 = 445000;

% Find the indices that correspond to x1 and x2
[val,id1] = min(abs(data(:,10)-x1));
[val,id2] = min(abs(data(:,10)-x2));

% Add sediments of uniform thickness in the interval
% sedthick = 2000;
% data(id1:id2,end) = data(id1:id2,end)-sedthick;
data(id1:id2,end) = -3000;

%% Gravity calculation

rhoi = 900;
rhosed = 2400;
rhoback = 2670;

gravcalc_HD2(LN,data,rhoi,rhosed,rhoback,2000)