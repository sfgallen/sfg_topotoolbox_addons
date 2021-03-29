function [xn, yn, dline, d2line] = snap2line(line, x, y)
% [xn, yn, dn] = snap2line(line, x, y)
% snap2line takes a line structure that has data properties of x
% coodinates, y coordinates, and distance and determines the approximate
% distance along the line of nearby x and y coordinate of points.
% The function outputs the snapped x and y coordinates and distance along
% the line for each point.
%
% Inputs:
% 1) line = a data structrue with x, y and distance vectors for the line.
% 2) x = vector or scalar of x coordinates of points
% 3) y = vector or scalar of y coordinates of points
%
% Outputs:
% 1) xn = vector or scalar of nearest x coordinate of point on line
% 2) yn = vector or scalar of nearest y coordinate ofof point on line
% 3) dline = distance along the line
% 4) d2line = shortest distance between point and line

    % test for arguments
    p = inputParser;
    p.FunctionName = 'snap2line';
    addRequired(p,'line', @(x) isstruct(x));
    addRequired(p,'x', @(x) isvector(x) | isscalar(x));
    addRequired(p,'y', @(x) isvector(x) | isscalar(x));
    parse(p,line,x,y);
    
    % some extra error handling
    if ~isfield(line,'x')
        error('line structure array missing "x" field');
    end
    if ~isfield(line,'y')
        error('line structure array missing "y" field')
    end
    if ~isfield(line,'d')
        error('line structure array missing "d" field')
    end
    if length(x) ~= length(y)
        error('x and y inputs are not the same length');
    end
    
    % declare valiables outside of loop
    xn = nan(length(x),1);
    yn = xn;
    dline = xn;
    d2line = xn;
    
    % find mimimum of least-squares residual
    for i = 1:length(x)
        res = (x(i) - line.x).^2 + (y(i) - line.y).^2;
        ind = find(res == min(res),1);
        xn(i) = line.x(ind);
        yn(i) = line.y(ind);
        dline(i) = line.d(ind);
        d2line(i) = sqrt(res(ind));
    end
end