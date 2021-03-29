function [idx,perpDist] = turingPointFinder(x,y)

% Two endpoints on the curve "data"
xend = [x(1) x(end)];
yend = [y(1) y(end)];
% The slope of the line connecting the two endpoints
m = ( yend(2) - yend(1) )/( xend(2) - xend(1) );
pm= - 1 / m;
% Point on the curve (xc,yc), point on the line (xl,yl)
perpDist = zeros(length(x),1);
for i = 1:length(x)
    xc = x(i) ; yc = y(i);
    yl = ( (m * xc) + (m^2 * yc) - (m * xend(1)) + yend(1) )/(1+ m^2);
    xl = xc - m*(yl - yc);
    % distance^2
    d2 = (xl - xc)^2 + (yl - yc)^2;
    perpDist(i) = (d2);
end
[val, idx] = max(perpDist);
end