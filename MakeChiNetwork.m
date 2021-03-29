function [Schi] = MakeChiNetwork(S,DA,mn,Ao)
%
% MakeChiNetwork.m will make calculate chi for a TopoToolbox stream
% network.
%
% Inputs:
% S         - TopoToolbox STREAMobj.
% DA        - Drainage area grid IN MAP UNITs (e.g. m^2) as a GRIDobj.
% mn        - m/n or refence concavity (theta) value.
%
% Outputs:
% Schi      - Chi for the stream network.
%
% Author: Sean F. Gallen
% Date modified: 11/19/2018
% email: sean.gallen{at}colostate.edu

p = inputParser;
p.FunctionName = 'MakeChiNetwork';
addRequired(p,'S', @(x) isa(x,'STREAMobj'));
addRequired(p,'A', @(x) isa(x,'GRIDobj'));
addRequired(p,'mn', @(x) isscalar(x));
addRequired(p,'Ao', @(x) isscalar(x));

parse(p,S,DA,mn,Ao);

% get variables ready for chi integration
Schi = zeros(size(S.distance));
Six = S.ix;
Sixc = S.ixc;
Sx = S.distance;
Sa = (Ao./(DA.Z(S.IXgrid))).^mn;       % chi transformation variable

h = waitbar(0,'calculating \chi for full stream network...');
% calculating chi and tau_star for the entire river network
for lp = numel(Six):-1:1
    Schi(Six(lp)) = Schi(Sixc(lp)) + (Sa(Sixc(lp))+(Sa(Six(lp))-Sa(Sixc(lp)))/2) *(abs(Sx(Sixc(lp))-Sx(Six(lp))));
    f = (numel(Six)+1 - lp)/numel(Six);
    waitbar(f,h);
end
close(h);