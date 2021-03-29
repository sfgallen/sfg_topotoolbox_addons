function line = lineXY2lineStruct(DEM,x,y)
% line = lineXY2lineStruct(DEM,x,y)
%
% lineXY2lineStruct will take a DEM GRIDobj (required) and x and y 
% coordinates (both required) of a line or verticies of a line and
% generates a structure array with properties of x and y coordinatea along
% the line, distance and elevation all in incements of the GRIDobj
% cellsize.
%
% Inputs
% 1) DEM GRIDobj(DEM).
% 2) x-coordinate (x)
% 3) y-coordinate (x)
%       x and y can be points along the line or verticies. 
%
% Outputs
% 1) 'line' as a structure with the following data properties
%  - x --> vector of x coordinates of line
%  - y --> vector of y coordinates of line
%  - d --> distance along the line
%  - z --> elevation along line
%
% Author: Sean F. Gallen
% Date Modified: 02/16/2017
% email: sean.gallen[at]erdw.ethz.ch

% test for arguments
p = inputParser;
p.FunctionName = 'selectLine';
addRequired(p,'DEM', @(x) isa(x,'GRIDobj'));
addRequired(p,'x', @(x) isscalar(x) | isvector(x));
addRequired(p,'y', @(x) isscalar(y) | isvector(y));
parse(p,DEM,x,y);

% get rid of NaNs if they exist
inds = find(~isnan(x) & ~isnan(y));
x = x(inds);
y = y(inds);

% Use TopoToolbox function "getdistance" to get distance along line
d = getdistance(x,y);

% determine distance along line in units of cellsize using interpline
[line.x,line.y,line.d] = interpline(x,y,d,DEM.cellsize);

% get elevation
line.z = interp(DEM,line.x,line.y);
end