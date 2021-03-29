function DBdata = DrainageBasinWatershed(FD,outlet)
% DrainageBasinWatershed calculates data related to the extent and location
% of the drainage basin. For example the x and y coordinates of the
% watershed (basin outline) and the linear indicies of the basin in the
% original data matrix. 
%
% Inputs:
% 1) flow direction as a FLOWobj
% 2) Outlet location as the linear index into the grid.
%
% Outputs:
% 1) drainage basin data (DBdata) as a data structure.
%   - DBdata.ix = linear indices of basin
%   - DBdata.x = x coordinates of watershed
%   - DBdata.y = y coordinates of watershed
%   - DBdata.distance = distance along watershed
%   - DBdata.outlet = outlet location (linear index)
% 
% function modified after TopoToolbox WATERSHEDobj

    % test for arguments
    p = inputParser;
    p.FunctionName = 'DrainageBasinWatershed';
    addRequired(p,'FD', @(x) isa(x,'FLOWobj'));
    addRequired(p,'outlet', @(x) isscalar(x));
    parse(p,FD,outlet);
    
    [y0,x0] = ind2sub(FD.size,outlet);

    % extract drainge basin
    db = drainagebasins(FD,outlet);
    DBdata.ix = find(db.Z);
%    XY = bwtraceboundary(db.Z,[y0,x0],'N',8);
    XY = bwtraceboundary(db.Z,[y0,x0],'N',4);
    rows = XY(:,1);
    cols = XY(:,2);

    % get coordinate pairs
    %[DB.x,DB.y] = sub2coord(FD,DB.ix);
    xy =  [double(rows(:)) double(cols(:)) ones(numel(rows),1)] * FD.refmat;
    DBdata.x = xy(1:end-1,1);
    DBdata.y = xy(1:end-1,2);

    DBdata.distance = getdistance(DBdata.x,DBdata.y);
    
    DBdata.outlet = outlet;
end