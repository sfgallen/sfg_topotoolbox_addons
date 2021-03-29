function [x,y] = getOrderedOutlets(S)
% getOrderedOutlets.m takes a STREAMobj and gets the x and y coordinates of
% the stream outlets in the order of the STREAMobj nanorderedlist.
%
% Inputs: S - STREAMobj
%
% OUTPUTS: x - x coordinates
%          y - y corrdinates
%
% Author: Sean F. Gallen
% email: sean.gallen[at]colostate.edu
% Date modified: 02/25/2-18

% Parse Inputs
p = inputParser;         
p.FunctionName = 'getOrderedOutlets';
addRequired(p,'S',@(x) isa(x,'STREAMobj'));

parse(p,S);
S     = p.Results.S;

ol = S.orderednanlist;

sb = find(isnan(ol));

x = nan(length(sb),1);
y = x;
id1 = 0;
for j = 1:length(sb);
    strmInds = ol(id1+1:sb(j)-1);
    x(j) = S.x(strmInds(end));
    y(j) = S.y(strmInds(end));
    id1 = sb(j);
end