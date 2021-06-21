
xx = -850000:2000:150000;
yy = -450000:2000:150000;
[XX,YY] = meshgrid(xx,yy);
XY = [reshape(XX',length(xx)*length(yy),1)...
       reshape(YY',length(xx)*length(yy),1)];

xi = -850000:500:150000;
yi = -450000:500:150000;
[XI,YI] = meshgrid(xi,yi);
XYI = [reshape(XI',length(xi)*length(yi),1)...
       reshape(YI',length(xi)*length(yi),1)];
   
sf = bedmachine_data('surface',xx,yy,'datum','ell');
firn = bedmachine_data('firn',xx,yy,'datum','ell');
bed = bedmachine_data('bed',xx,yy,'datum','ell');
thick = bedmachine_data('thickness',xx,yy,'datum','ell');
errbed = bedmachine_data('errbed',xx,yy);

% True surface elevation
sftrue = sf+firn;

% Flip upside down
sftrue = flipud(sftrue);
bed = flipud(bed);

% Reshape into vectors
sfvec = reshape(sftrue',length(xi)*length(yi),1);
bedvec = reshape(bed',length(xi)*length(yi),1);


Fsf = scatteredInterpolant(XYI(:,1),XYI(:,2),sfvec,'natural');
Fbed = scatteredInterpolant(XYI(:,1),XYI(:,2),bedvec,'natural');

%% Interpolate onto 2-km grid
sfi = Fsf(XY(:,1),XY(:,2));
bedi = Fbed(XY(:,1),XY(:,2));

BM_2km = [XY sfi bedi];

save BM_2km BM_2km

%% Interpolate onto 5-km grid

xx = -850000:5000:150000;
yy = -450000:5000:150000;
[XX,YY] = meshgrid(xx,yy);
XY = [reshape(XX',length(xx)*length(yy),1)...
       reshape(YY',length(xx)*length(yy),1)];

sfi = Fsf(XY(:,1),XY(:,2));
bedi = Fbed(XY(:,1),XY(:,2));

BM_5km = [XY sfi bedi];

save BM_5km BM_5km
