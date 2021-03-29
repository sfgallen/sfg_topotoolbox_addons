function MS = XYcoords2PointMS(x,y)
% XYcoords2PointMS takes x and y coordinates and makes a map structure
% that can be used to make a point shapefile with shapewrite
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
p.FunctionName = 'XYcoords2PointMS';
addRequired(p,'x', @(x) isvector(x) | isscalar(x));
addRequired(p,'y', @(x) isvector(x) | isscalar(x));

parse(p, x, y);

if length(x) ~= length(y)
    error('x and y inputs are not the same length!')
end

% You map also want to make a shapefile of your outlet locations
MS = struct('Geometry',{'Point'},...
    'X',num2cell(x),...
    'Y',num2cell(y),...
    'ID',num2cell((1:length(x))'));
end