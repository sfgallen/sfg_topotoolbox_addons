function [d, r, A, dx, outlet_nodes] = fastscape_eroder_data_prep(S, S_DA, S_K, m)

% simplified from Campforts and Schwanghart TTLEM updateDrainDir.m function
% to handel incision along STREAMobj only
%
% Inputs:
%   S               STREAMobj
%   S_DA            drainage area in meters^2 along the stream network
%   S_K             A scalar or vector of the erodibility coefficent
%   m               The drainage area exponent
%
% Outputs:
%   d               donors
%   r               recievers
%   A               the "velocity field" of the steam power incision model
%   dx              distance between nodes along network
%   outlet_nodes	locations of outlets

% get donor and recievers from STREAMobj
d = S.ix;
r = S.ixc;

% get velocity field
A = S_K.*S_DA.^m;

% Find outlet nodes in STREAMobj
outlet_nodes = streampoi(S,'outlets','ix');
for i = 1:length(outlet_nodes)
    outlet_nodes(i) = find(S.IXgrid == outlet_nodes(i));
end

% get distance between nodes along profile
Sd = S.distance;
dx = abs(Sd(d) - Sd(r));

end
