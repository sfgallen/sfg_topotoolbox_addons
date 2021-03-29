function SmoZ = smooth_channel_elevations(S,DEM,smoWin)

% Declare STREAMobj variables for faster processing through forloop
ordList = S.orderednanlist;
strmBreaks = find(isnan(ordList));
Sz = DEM.Z(S.IXgrid);
SmoZ = Sz;
cs = DEM.cellsize;

% plot all of the river river profile data as thin gray lines
h = waitbar(0,'Smoothing data and plotting all streams...');

id1 = 0;
for i = 1:length(strmBreaks)
    strmInds = ordList(id1+1:strmBreaks(i)-1);
    SmoZ(strmInds) = smoothChannelZ(SmoZ(strmInds),smoWin,cs);
    id1 = strmBreaks(i);
    f = i/length(strmBreaks);
    waitbar(f,h);
end
close(h);