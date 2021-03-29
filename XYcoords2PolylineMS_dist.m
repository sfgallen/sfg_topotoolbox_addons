function MS = XYcoords2PolylineMS_dist(x,y,d)
% XYcoords2Polyline takes x and y coordinates and makes a map structure
% that can be used to make a line shapefile with shapewrite
%
% example:
% MS = XYcoords2Polyline(x,y);
% shapewrite(MS,'your_file_name.shp');
%
% Inputs:
% 1) x-coordinates (required)
% 2) y-coordinates (required)
% 2) d-distance (required)
%
% Outputs:
% 1) polyline map structure
%
% Author: Sean F. Gallen
% Date Modified: 02/10/2017
% email: sean.gallen@erdw.ethz.ch

p = inputParser;
p.FunctionName = 'XYcoords2PolylineMS_dist';
addRequired(p,'x', @(x) isvector(x) | isscalar(x));
addRequired(p,'y', @(x) isvector(x) | isscalar(x));
addRequired(p,'d', @(x) isvector(x) | isscalar(x));

parse(p, x, y, d);

if length(x) ~= length(y)
    error('x and y inputs are not the same length!')
end

MS = struct('Geometry',{},...
    'X',{},...
    'Y',{},...
    'D',{},...
    'ID',{});

for i = 1:length(x)
    MS(i).Geometry = 'Line';
    MS(i).X = x(i);
    MS(i).Y = y(i);
    MS(i).D = d(i);
    MS(i).ID = 1;
end


end