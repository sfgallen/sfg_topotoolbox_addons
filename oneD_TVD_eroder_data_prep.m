function [d, r, A, dx, outlet_nodes, dd, rr, dx_centered] = oneD_TVD_eroder_data_prep(S, S_DA, S_K, m)

% simplified from Campforts and Schwanghart TTLEM updateDrainDir.m function
% to handle incision along STREAMobj only using the TVD scheme
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

% For the second order TVD scheme, calculate upstream and second downstream
% cells
[dd,rr] =get1up1down(d,r,A);
ind_dd=dd;
ind_dd(isnan(dd))=1;
dx_dd = abs(Sd(ind_dd)-Sd(d));
dx_centered=(dx_dd+dx)/2;
dx_centered(isnan(dd))=dx(isnan(dd));



function [ii,kk] = get1up1down(i,k,A)
% get receiver of receiver and giver of giver...
% ii and kk will have nans where neighbors do not exist, 
% either because,%there is none 
% or because the upstream node has a smaller contribution area.

i = double(i);
k = double(k);

nrc = numel(A);

% Downstream neighbor
kk = nan(nrc,1);
kk(i) = k;
kk(i) = kk(k);
kk = kk(i);

% Upstream neighbor
ii = accumarray(k,i,[nrc 1],@getindexwithlargerarea,nan); ii = ii(i);



function ix = getindexwithlargerarea(ix)

[~,mix] = max(A(ix));
ix = ix(mix);
end
end
end
