function river_slope = calc_river_slope(DEM,S)


Sz = DEM.Z(S.IXgrid);

river_slope = nan(size(S.distance));

river_slope(S.ix) = (Sz(S.ixc)-Sz(S.ix))./(S.distance(S.ixc) - S.distance(S.ix));

% outlets = streampoi(S,'outlets','logical');
% outlets = find(outlets == 1);
channelheads = streampoi(S,'channelheads','logical');
channelheads  = find(channelheads  == 1);

lia = ismember(S.ix,channelheads);
river_slope(channelheads) = (Sz(channelheads)-Sz(S.ixc(lia)))./(S.distance(channelheads) - S.distance(S.ixc(lia)));

