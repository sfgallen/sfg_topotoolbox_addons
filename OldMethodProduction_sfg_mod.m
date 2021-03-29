function [Production] = OldMethodProduction_sfg_mod(DEM,utm_zone,tag)
% Lupker et al. (2012) 10Be Erosion Rate Calculation

% The glacier raster ICE should be structured such that areas that are
% within the DB and are ice free have a value of 1 and all other areas are
% nan

pixcelsize=DEM.cellsize; % pixcel size of DEM
[~,X,Y] = GRIDobj2mat(DEM);

% find max and min X and Y, and convert to lat lon
[lat_min,lon_min]=utm2ll(min(X),min(Y),utm_zone);
[lat_max,lon_max]=utm2ll(max(X),max(Y),utm_zone);

% Create vector of degree coordinates to interpolate 
lat_in=(lat_min:0.1:lat_max+0.1)';
lon_in=(lon_min:0.1:lon_max+0.1)';

% scaling raster, pre-calculated based on values in Stone, 2000
Lat=[0,10,20,30,40,50,60,90];
a=[31.8518,34.3699,40.3153,42.0983,56.7733,69.0720,71.8733,71.8733];
b=[250.3193,258.4759,308.9894,512.6857,649.1343,832.4566,863.1927,863.1927];
c=[-0.083393,-0.089807,-0.106248,-0.120551,-0.160859,-0.199252,-0.207069,-0.207069];
d=[7.4260e-5,7.9457e-5,9.4508e-5, 1.1752e-4,1.5463e-4,1.9391e-4,2.0127e-4,2.0127e-4];
e=[-2.2397e-8,-2.3697e-8,-2.8234e-8,-3.8809e-8,-5.0330e-8,-6.3653e-8,-6.6043e-8,-6.6043e-8];
m=[0.587,0.600,0.678,0.833,0.933,1.000,1.000,1.000];

% Get x,y coordinates of our query points
utmstruct = defaultm('utm'); 
zonetext = sprintf('%dN',utm_zone); % hard coding this to be only UTM North
utmstruct.zone = zonetext;  
utmstruct.geoid = wgs84Ellipsoid;
utmstruct = defaultm(utmstruct);
[~,y]=mfwdtran(utmstruct,lat_in,repmat(lon_in(1),length(lat_in),1));

% Generate x,y arrays for grid cells (where we want to finally interpolate things
% to)
x1=min(X):pixcelsize:max(X); % create x array for interpolation, range constrained by DEM
y1=(min(Y):pixcelsize:max(Y))'; % create y array for interpolation
if length(y1) == (length(Y) - 1)
   y1 = [y1; max(y1) + pixcelsize]; 
end


% interpolation of scaling fator - to find values at [x1(:),y1(:)]
A=GRIDobj(x1,y1,repmat(interp1(y,interp1(Lat,a,lat_in),y1),1,length(x1)));
B=GRIDobj(x1,y1,repmat(interp1(y,interp1(Lat,b,lat_in),y1),1,length(x1)));
C=GRIDobj(x1,y1,repmat(interp1(y,interp1(Lat,c,lat_in),y1),1,length(x1)));
D=GRIDobj(x1,y1,repmat(interp1(y,interp1(Lat,d,lat_in),y1),1,length(x1)));
E=GRIDobj(x1,y1,repmat(interp1(y,interp1(Lat,e,lat_in),y1),1,length(x1)));
%M=GRIDobj(x1,y1,repmat(interp1(y,interp1(Lat,m,lat_in),y1),1,length(x1)));

Pn_SLHL=4.01;   % Borchers et al., 2015 
Pms_SLHL=0.012; % Braucher et al., 2011 table 6
Pmf_SLHL=0.039; % Braucher et al., 2011


pres=1013.25*exp(((-0.03417)/6.5e-3)*(log(288.15)-log(288.15-(6.5e-3*DEM.Z)))); % Stone, 2000

Pn= Pn_SLHL.*(A.Z+B.Z.*exp(-pres/150)+C.Z.*pres+D.Z.*pres.^2+E.Z.*pres.^3); % Stone, 2000
Pms=Pms_SLHL.*exp((1013.25-pres)/260); % 1013.25: sea level pressure, 260: muon attenuation lengths in the air(g/cm2) %Braucher et al., 2011
Pmf=Pmf_SLHL.*exp((1013.25-pres)/510); % Braucher et al., 2011

% Calculate area production is computed over
was_used = ~isnan(DEM.Z);
num_used_cells = sum(sum(was_used));
area = (DEM.cellsize^2.0) * num_used_cells;

% Prepare output structure
Production = struct();
Production.Pn = Pn;
Production.Pms= Pms;
Production.Pmf= Pmf;
Production.area=area;
Production.tag=tag;

end

