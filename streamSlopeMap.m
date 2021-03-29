function [sloGrid,sloGridDeg] = streamSlopeMap(DEM, fileTag, crita, smoWin)
% streamSlopeMap.m makes a river slope map shapefile from a TopoToolBox DEM
% GRIDobj
%
% it takes six inputs:
%       1) DEM: DEM GRIDobj
%       2) filetag: string used to name output shapefile
%       3) crita: critical drainage area for channel head initiation
%       4) smoWin: size of smoothing window in map units used for ksn
%       calculation
% Outputs:
%       1) ChiGrid: chi map GRIDobj
%       2) ksnGrid: ksn map GRIDobj
%       3) A shapefile with chi and ksn attributes
%
% Code uses three additional functions:
%       1) smoothChannelZ.m: smooths channel elevation
%       2) fastsmooth.z: script used to help smooth channels by T. C. O'Haver
%       3) binnedKsn.m: calcuates ksn map using a smoothing window average
%
% Author: Sean F. Gallen
% Date modified: 12/31/2015

% set nan values for common nan numbers
DEM.Z(DEM.Z==-32767) = nan;
DEM.Z(DEM.Z==-32768) = nan;
DEM.Z(DEM.Z==-32456) = nan;
DEM.Z(DEM.Z==-9999) = nan;

% Fill sinks and declare cellsize
DEM = fillsinks(DEM);
cs = DEM.cellsize;

% Calculate flow object (FD) and distance to divide (dfd)
FD  = FLOWobj(DEM,'preprocess','carve');

S1 = STREAMobj(FD,'minarea',crita/(DEM.cellsize^2));

% Declare STREAMobj variables for faster processing through forloop
ordList = S1.orderednanlist;
strmBreaks = find(isnan(ordList));

GridID = S1.IXgrid;

Sz = double(DEM.Z(GridID));                 % elevation
SmoZ = Sz;               % dumby vector to get smoothed data

% get variables ready for chi integration
chis = zeros(size(S1.distance));
Six = S1.ix;
Sixc = S1.ixc;
Sx = S1.distance;


% plot all of the river river profile data as thin gray lines
h = waitbar(0,'Smoothing elevation data for full stream network...');
id1 = 0;
for i = 1:length(strmBreaks);
    strmInds = ordList(id1+1:strmBreaks(i)-1);
    SmoZ(strmInds) = smoothChannelZ(Sz(strmInds),smoWin,cs);
    id1 = strmBreaks(i);
    f = i/length(strmBreaks);
    waitbar(f,h);
end
close(h)

% making ksn map for streams analyzed
sloStreams = binnedSlope(S1,SmoZ,smoWin,cs);

sloGrid = DEM;
sloGrid.Z = nan(size(DEM.Z));
sloGrid.Z(GridID) = sloStreams;

sloGridDeg = DEM;
sloGridDeg.Z = nan(size(DEM.Z));
sloGridDeg.Z(GridID) = atand(sloStreams);

MS = STREAMobj2mapstruct(S1,'seglength',smoWin,'attributes',...
    {'sloDeg' sloGridDeg @mean...
    'slo' sloGrid @mean});
fileName = [fileTag, '_RivSlope_map.shp'];
shapewrite(MS,fileName);