%% stream power incision model and associated fucntions
function [Z0] = fastscape_incision_forward_w_SA(S,S_DA,Ui,Uf,Ki,Kf,m,n,run_time,dt,nplot,varargin)
% This function used the implicit finite difference scheme to solve the one
% dimensional stream power model (E = K*A^m*S^n) on a TopoToolbox STREAMobj
% using the implicit finite difference scheme of Braun and Willet (2013). 
% Some of the codes are modified and simplified versions of the solvers 
% available with the TTLEM. 
% 
% required inputs
% (1) S - STREAMobj from TopoToolbox
% (2) S_DA - Drainage area vector mapped to STREAMobj
% (3) Ui - the intial uplift rate in m/yr [Can be scalar or vector]
% (4) Uf - the final uplift m/yr [Can be a scalar or vector]
% (5) Ki - intial erodibility constant [Can be scalar or vector]
% (6) Kf - final erodibility constant [Can be scalar or vector]
% (7) m - drainage area exponent (scalar)
% (8) n - slope exponent (scalar)
% (9) run_time - length of model run time in years (scalar)
% (10) dt - time step used in the model in years (scalar)
% (11) nplot - number or timesteps plotted (scalar)
%
% optional input
% (12) Zi - user defined initial elevation of the river profile (e.g.
% extracted from the DEM). If no input is given an initial steady state
% river profile will be derived based on Ui and Ki.

% Parse Inputs
p = inputParser;         
p.FunctionName = 'fastscape_incision_forward';

% required inputs
addRequired(p,'S',@(x) isa(x,'STREAMobj'));
addRequired(p,'S_DA',@(x) isvector(x));
addRequired(p,'Ui',@(x) isscalar(x) | isvector(x));
addRequired(p,'Uf',@(x) isscalar(x) | isvector(x));
addRequired(p,'Ki',@(x) isscalar(x) | isvector(x));
addRequired(p,'Kf',@(x) isscalar(x) | isvector(x));
addRequired(p,'m',@(x) isscalar(x));
addRequired(p,'n',@(x) isscalar(x));
addRequired(p,'run_time',@(x) isscalar(x));
addRequired(p,'dt',@(x) isscalar(x));
addRequired(p,'nplot',@(x) isscalar(x));


% optional inputs
addOptional(p,'Zi', 0, @(x) isvector(x));

% parse inputs
parse(p,S,S_DA,Ui,Uf,Ki,Kf,m,n,run_time,dt,nplot,varargin{:});
S = p.Results.S;
S_DA = p.Results.S_DA;
Ui = p.Results.Ui;
Uf = p.Results.Uf;
Ki = p.Results.Ki;
Kf = p.Results.Kf;
m = p.Results.m;
n = p.Results.n;
run_time = p.Results.run_time;
dt = p.Results.dt;
nplot = p.Results.nplot;
Zi = p.Results.Zi;



% Author: Sean F. Gallen
% email: sean.gallen[at]colostate.edu
% date modified: 02/21/2018

% calculate the change in drainage area between cells for sed flux calc.
% fund the fastscape data preperation function
[d, r, A, dx, outlet_nodes] = fastscape_eroder_data_prep(S, S_DA, Kf, m);

% Calculate Chi for the river network
mn = m/n;
chi = calculate_chi(S,S_DA,mn);

% don't allow negative uplift rates
Uf(Uf <= 0) = 1e-7;

% set initial conditions
if Zi == 0
    % calculate initial steady-state river elevations
    S_Zi = calculate_z(S,S_DA,Ui,Ki,mn,n);
else
    S_Zi = Zi;
end

% plot initial data
Sd = S.distance;

subplot(3,1,1)
plot(Sd,S_Zi,'k.','markersize',2); hold on
xlabel('Distance (m)'); ylabel('Elevation (m)');

subplot(3,1,2)
plot(chi,S_Zi,'k.','markersize',2); hold on
xlabel('\chi'); ylabel('Elevation (m)');

slo = getslope(S_Zi,Sd);
subplot(3,1,3)
plot(S_DA,slo,'k.','markersize',3); hold on
xlabel('Drainage Area (m^2)'); ylabel('Slope');
set(gca,'xscale','log','yscale','log','xdir','reverse')

% get time loop prepped
tsteps = round(run_time/dt);
tplot = round(tsteps/nplot);
cols = jet(tsteps);

% allocate memory before running the model and calc. initial erosion rate
% and sediment flux
Z0 = S_Zi;

% run the forward model
h = waitbar(0,'running model...');
for t = 1:tsteps
    % update uplift field as needed
    Z1 = fastscape_eroder_outlets(Z0, n, dt, A, d, r, dx, Uf, outlet_nodes);
    
    % make sure river doesn't flow backwards
    Z1 = check_z(S,Z1);
    
    % update Z0 data with new elevations
    Z0 = Z1;
    
    % plot data as needed
    if rem(t,tplot)==0
        
        
        subplot(3,1,1)
        plot(Sd,Z1,'k.','color',cols(t,:),'markersize',2); hold on
        xlabel('Distance (m)'); ylabel('Elevation (m)');
        
        subplot(3,1,2)
        plot(chi,Z1,'k.','color',cols(t,:),'markersize',2); hold on
        xlabel('\chi'); ylabel('Elevation (m)');
        
        slo = getslope(Z1,Sd);
        subplot(3,1,3)
        plot(S_DA,slo,'ko','color',cols(t,:),'markersize',3); hold on
        xlabel('Drainage Area (m^2)'); ylabel('Slope');
        set(gca,'xscale','log','yscale','log','xdir','reverse')
        
        pause(2)
    end
    
    waitbar(t/tsteps,h);
end
close(h)
end


% calculate steady-state channel elevations
function z = calculate_z(S,S_DA,S_U,K,mn,n)
% Author: Sean F. Gallen
% email: sean.gallen[at]colostate.edu
% date modified: 02/21/2018

z = zeros(size(S.distance));
Six = S.ix;                         % donors
Sixc = S.ixc;                       % recievers
Sd = S.distance;                    % distance from mouth
Sa = (S_U./K).^(1/n).*(1./(S_DA)).^mn;       % chi transformation variable

%h = waitbar(0,'calculating z...');
% calculating chi for the entire river network
for lp = numel(Six):-1:1
    z(Six(lp)) = z(Sixc(lp)) + (Sa(Sixc(lp))+(Sa(Six(lp))-Sa(Sixc(lp)))/2) *(abs(Sd(Sixc(lp))-Sd(Six(lp))));
    %f = (numel(Six)+1 - lp)/numel(Six);
    %waitbar(f,h);
end
%close(h);
end

% make sure that the river doesn't flow backwards
function z = check_z(S,z)
Six = S.ix;                         % donors
Sixc = S.ixc;                       % recievers

% calculating chi for the entire river network
for lp = numel(Six):-1:1
    z_pre = z(Sixc(lp));
    z_cur = z(Six(lp));
    if z_cur <= z_pre
        z(Six(lp)) = z_pre+0.1;
    end
end
end

function [d, r, A, dx, outlet_nodes] = fastscape_eroder_data_prep(S, S_A, S_K, m)

% simplified from Campforts and Schwanghart TTLEM updateDrainDir.m function
% to handel incision along STREAMobj only
%
% Inputs:
%   S               STREAMobj
%   S_A            drainage area in meters^2 along the stream network
%   S_K             A scalar or vector of the erodibility coefficent
%   m               The drainage area exponent
%
% Outputs:
%   d               donors
%   r               recievers
%   A               the "velocity field" of the steam power incision model
%   dx              distance between nodes along network
%   outlet_nodes	locations of outlets
%
% Author: Sean F. Gallen
% email: sean.gallen[at]colostate.edu
% date modified: 02/21/2018


% get donor and recievers from STREAMobj
d = S.ix;
r = S.ixc;

% get velocity field
A = S_K.*S_A.^m;
% A = getnal(S,A);

% Find outlet nodes in STREAMobj
outlet_nodes = streampoi(S,'outlets','ix');
for i = 1:length(outlet_nodes)
    outlet_nodes(i) = find(S.IXgrid == outlet_nodes(i));
end

% get distance between nodes along profile
Sd = S.distance;
dx = abs(Sd(d) - Sd(r));

end

function S_Z = fastscape_eroder_outlets(S_Z, n, dt, A, d, r, dx, U, outlets)

% simplified from Campforts and Schwanghart TTLEM funerosion_implin.m and
% funerosion_impnlin.m to model incision along STREAMobj only using the
% detachment limited stream power incision model. There is a data
% preperation frunction, fastscape_eroder_data_prep.m that should be run
% before one goes into the time evolution forloop where this function is
% applied.
%
% Inputs:
%   S_Z     River network elevations in meters sorted as STREAMobj
%   n       Slope exponent in stream power incision model
%   dt      time step in years
%   A       the "velocity field" of the steam power incision model
%   d       STEAMobj donors
%   r       STREAMobj recievers
%   dx      distance between nodes in STREAMobj
%   U       scalar or vector of uplift rate in meters per year
%
% Outputs
%   S_Z     Updated river network elevations in meters sorted as STREAMobj
%
% Author: Sean F. Gallen
% email: sean.gallen[at]colostate.edu
% date modified: 08/30/2019

time=dt;
dte = dt;

while time>0
    time=time-dte;
    if time<0
        dte=dte+time;
        time=0;
    end
    
    S_Z=S_Z+dte.*U; % add uplift to elevation field
    if isscalar(U)
        S_Z(outlets) = S_Z(outlets)-dte.*U; %except for the outlets
    else
        S_Z(outlets) = S_Z(outlets)-dte.*U(outlets);
    end
    
    if n == 1
        for j = numel(d):-1:1;
            tt      = A(d(j))*dte/(dx(j));
            S_Z(d(j)) = (S_Z(d(j)) + S_Z(r(j))*tt)./(1+tt);
        end
    else
        for j = numel(d):-1:1;
            
            tt      = A(d(j))*dte/(dx(j));
            % z_t
            zt      = S_Z(d(j));
            % z_(t+dt) of downstream neighbor
            ztp1d   = S_Z(r(j));
            % dx
            dx_n      = dx(j);
            
            % initial value for finding root
            if ztp1d < zt;
                ztp1    = newtonraphson(zt,ztp1d,dx_n,tt,n);
            else
                ztp1    = zt;
            end
            
            if ~isreal(ztp1) || isnan(ztp1)
               % disp('Non real solutions converted to real')
                ztp1=real(ztp1);
            end
            S_Z(d(j))=ztp1;
        end
    end
end
    function ztp1 = newtonraphson(zt,ztp1d,dx,tt,n)
        
        tempz   = zt;
        tol = inf;
        
        while tol > 1e-3;
            % iteratively approximated value
            ztp1  =  tempz - (tempz-zt + ...
                (tt*dx) * ...
                ((tempz-ztp1d)./dx)^n) / ...
                (1+n*tt*((tempz-ztp1d)./dx)^(n-1));
            tol   = abs(ztp1-tempz);
            tempz = ztp1;
        end
    end
end

function chi = calculate_chi(S,S_DA,mn)
    chi = zeros(size(S.distance));
    Six = S.ix;                         % donors
    Sixc = S.ixc;                       % recievers
    Sd = S.distance;                    % distance from mouth
    Sa = (1./(S_DA)).^mn;       % chi transformation variable
    
    h = waitbar(0,'calculating \chi...');
    % calculating chi for the entire river network
    for lp = numel(Six):-1:1
        chi(Six(lp)) = chi(Sixc(lp)) + (Sa(Sixc(lp))+(Sa(Six(lp))-Sa(Sixc(lp)))/2) *(abs(Sd(Sixc(lp))-Sd(Six(lp))));
        f = (numel(Six)+1 - lp)/numel(Six);
        waitbar(f,h);
    end
    close(h);
end

function S = getslope(Z,D)
S = zeros(size(Z));
S(end) = abs(Z(end-1) - Z(end))./abs(D(end-1) - D(end));
S(1:end-1) = abs(Z(1:end-1) - Z(2:end))./abs(D(1:end-1) - D(2:end));
end
