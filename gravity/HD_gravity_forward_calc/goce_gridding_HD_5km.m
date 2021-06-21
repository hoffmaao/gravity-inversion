
% This code interpolates the GOCE Antarctic gravity field onto a 5-km grid
% around Hercules Dome area. 
% grd_goce5s_do280_60s.mat is the GOCE gravity field provided by Ingo
% Sasgen at GFZ Potsdam. See section 2.2 of Muto et al. (2016, EPSL) for
% details. Originally in grd_goce5s_do280_60s.txt, added PSX and PSY in the 
% third and fourth columns of grd_goce5s_do280_60s.mat

%%

load grd_goce5s_do280_60s.mat

id = isnan(grd_goce5s_do280_60s(:,5))==1;
grd_goce5s_do280_60s(id,:) = [];

xx = -850000:5000:150000;
yy = -450000:5000:150000;
[XI,YI] = meshgrid(xx,yy);

F = TriScatteredInterp(grd_goce5s_do280_60s(:,3),grd_goce5s_do280_60s(:,4),grd_goce5s_do280_60s(:,5),'natural');
ZI = F(XI,YI);

%%
S = size(ZI);
goce_HD_5km = nan(S(1),S(2),3);

goce_HD_5km(:,:,1) = XI;
goce_HD_5km(:,:,2) = YI;
goce_HD_5km(:,:,3) = ZI;

% save goce_HD_5km goce_HD_5km

