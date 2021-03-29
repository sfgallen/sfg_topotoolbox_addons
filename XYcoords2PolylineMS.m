function MS = XYcoords2PolylineMS(x,y)
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
%
% Outputs:
% 1) polyline map structure
%
% Author: Sean F. Gallen
% Date Modified: 02/10/2017
% email: sean.gallen@erdw.ethz.ch

p = inputParser;
p.FunctionName = 'XYcoords2PolylineMS';
addRequired(p,'x', @(x) isvector(x) | isscalar(x));
addRequired(p,'y', @(x) isvector(x) | isscalar(x));

parse(p, x, y);

if length(x) ~= length(y)
    error('x and y inputs are not the same length!')
end

MS(1).Geometry = 'Line';
MS(1).X = x;
MS(1).Y = y;
MS(1).ID = 1;
MS(1).gridval = 1;
end