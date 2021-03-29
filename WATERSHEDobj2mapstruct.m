function C = WATERSHEDobj2mapstruct(WS)
% WATERSHEDobj2mapstruct takes an instance of a WATERSHEDobj and converts
% it into a mapstructure that can be convered into a shapefile with the
% Matlab command "shapewrite"
%
% Inputs:
% 1) WATERSHEDobj
% Outputs:
% 1) watershed map structure
%
% Author: Sean F. Gallen
% Date Modified: 02/08/2017

    % test for arguments
    p = inputParser;
    p.FunctionName = 'WATERSHEDobj2mapstruct';
    addRequired(p,'WS', @(x) isa(x,'WATERSHEDobj'));
    parse(p,WS);
    
    %[x,y] = ind2coord(DEM,WS.ix);
    C(1).Geometry = 'Polygon';
    C(1).X = WS.x;
    C(1).Y = WS.y;
    C(1).ID = 1;
end
