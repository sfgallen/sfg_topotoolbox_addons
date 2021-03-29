function [KsnMod,KsnStd,A,Schi] = invert_for_ksn(DEM,GRID,varargin)


% Parse Inputs
p = inputParser;         
p.FunctionName = 'invert_for_ksn';

% required inputs
addRequired(p,'DEM', @(x) isa(x,'GRIDobj'));
addRequired(p,'GRID', @(x) isa(x,'GRIDobj'));

% optional inputs
addOptional(p,'crita', 1e6, @(x) isscalar(x));
addOptional(p,'mn', 0.5, @(x) isscalar(x));
addOptional(p,'flow_option', []);
addOptional(p,'inverse_option', []);

% declare variable names
parse(p,DEM, GRID, varargin{:});
DEM   = p.Results.DEM;
GRID     = p.Results.GRID;
crita    = p.Results.crita;
mn = p.Results.mn;

% test to make sure grids have the same cellsize
if DEM.cellsize ~= GRID.cellsize
    error('The DEM and other grid do not have the same cellsize');
end

% declare cellsize
cs = DEM.cellsize;

% make sure that grids overlap properly (also sets common nan values to
% nan)
[DEM, GRID] = largest_overlapping_extent(DEM,GRID);

% flow routing options
if isempty(p.Results.flow_option)
    FD = FLOWobj(DEM);
elseif strcmp(p.Results.flow_option, 'fill')
    DEM = fillsinks(DEM);
    FD = FLOWobj(DEM);
elseif strcmp(p.Results.flow_option, 'carve')
    FD = FLOWobj(DEM,'preprocess','carve');
    DEM = imposemin(FD,DEM);
else
    error('fillOption is not "fill" or "carve"');
end

% inverse technique options
if isempty(p.Results.inverse_option)
    in_op = 1;
elseif strcmp(p.Results.inverse_option, 'non_negative')
    in_op = 1;
elseif strcmp(p.Results.inverse_option, 'Tikhonov')
    in_op = 2;
else
    error('fillOption is not "non_negative" or "non_negative"');
end

% make flow accumulation grid
DA = flowacc(FD).*(FD.cellsize^2);

% Create stream network object (S)
S  = STREAMobj(FD,'minarea',crita/(cs)^2);

% make the forward operator
[A,Schi] = MakeChiGeoAMatrix(S,DA,GRID,mn);

[N,q] =size(A);

% get stream network elevations and set outlet to zero elevation
Sz = DEM.Z(S.IXgrid);
Sz = Sz - min(Sz);

% invert for ksn using preferred method.
if in_op == 1
    KsnMod = lsqnonneg(A,double(Sz));
elseif in_op == 2
    % make q by q identity matrix
    I = eye(q,q);
    
    % calculate KsnPri
    KsnPri = ones(q,1)*(max(Sz)./max(Schi));
    
    % find the optimal dampening (e.g. smoothing) parameter (Gamma)
    Gam = logspace(log10(1),log10(1e5),100);
    MisFit = nan(size(Gam));
    
    for i = 1:length(Gam)
        KsnMod = KsnPri + (A'*A + Gam(i)^2*I)\A'*(Sz-A*KsnPri);
        MisFit(i) = (1/(N-q))*sqrt(sum((Sz - A*KsnMod).^2));
        
    end
    
    [ind,~] = turingPointFinder(1./Gam,MisFit);
    
    figure()
    plot(1./Gam,MisFit,'k-'); hold on
    plot(1/Gam(ind),MisFit(ind),'ko','markerfacecolor',[0.5 0.8 0.5]);
    xlabel('1/\Gamma'); ylabel('normalized misfit');
    H=text(1/Gam(ind)+0.002,MisFit(ind)+0.025,['best-fit \Gamma = ',num2str(Gam(ind))]);
    set(H,'Fontsize',10);
    title('Trade off curve')
    
    % now we can invert to find U following Goren eq (21) from Tarantola, 1987
    KsnMod = KsnPri + (A'*A + Gam(ind)^2*I)\A'*(Sz-A*KsnPri);
end

% get the standard deviation of the calculation
[~,KsnStd] = lscov(A,Sz);
end


%%%%%%%%%%%%%% additional functions used in code above %%%%%%%%%%%%%%%%%%%%

%% Series of functions to clip grids to largest overlapping extent
function [GRID1, GRID2] = largest_overlapping_extent(GRID1,GRID2)
%
% largest_overlapping_extent.m take two TopoToolbox GRIDobjects and resizes
% them to their largests overlapping extent.

%% set nan values for common nan numbers
GRID1.Z(GRID1.Z<=-9999) = nan;
GRID1.Z(GRID1.Z>=9999) = nan;

% set nan values for common nan numbers
GRID2.Z(GRID1.Z<=-9999) = nan;
GRID2.Z(GRID1.Z>=9999) = nan;

%% resize grids
% make sure grids are the same size
[xMin(1,1),xMax(1,1),yMin(1,1),yMax(1,1)] = findCorners(GRID1);
[xMin(2,1),xMax(2,1),yMin(2,1),yMax(2,1)] = findCorners(GRID2);

% resizing grids to largest overlapping area
xMinP = max(xMin);
xMaxP = min(xMax);
yMinP = max(yMin);
yMaxP = min(yMax);

GRID1 = gridReSize(GRID1,xMinP,xMaxP,yMinP,yMaxP);
GRID2 = gridReSize(GRID2,xMinP,xMaxP,yMinP,yMaxP);

%% revise GRIDobj and the refmat as needed
[Ny1,Nx1] = size(GRID1.Z);
[Ny2,Nx2] = size(GRID2.Z);

% fix rows if they aren't the same length
if length(Ny1) < length(Ny2)
    l_dif = length(Ny2) - length(Ny1);
    GRID2.Z = GRID2.Z(1:end-l_dif,:);
elseif length(Ny1) > length(Ny2)
    l_dif = length(Ny1) - length(Ny2);
    GRID1.Z = GRID1.Z(1:end-l_dif,:);
end

% fix columns if they aren't the same length
if length(Nx1) < length(Nx2)
    l_dif = length(Nx2) - length(Nx1);
    GRID2.Z = GRID2.Z(:,1:end-l_dif);
elseif length(Nx1) > length(Nx2)
    l_dif = length(Nx1) - length(Nx2);
    GRID1.Z = GRID1.Z(:,1:end-l_dif);
end

GRID1.size = size(GRID1.Z);
GRID2.size = size(GRID2.Z);

% % fix refmat
% [Ny1,Nx1] = size(GRID1.Z);
% 
% GRID1.georef.SpatialRef.XWorldLimits = [0, Nx1*GRID1.cellsize];
% GRID1.georef.SpatialRef.YWorldLimits = [0, Ny1*GRID1.cellsize];

GRID2.refmat = GRID1.refmat;

end

function [xMinP,xMaxP,yMinP,yMaxP] = findCorners(inGrid)
% the findCorners function takes an acsii grid input made with the function
% 'makeGrid' and finds the UTM coorinates of the corners of the grid.
% Created: 11/06/2013
% Author: Sean F. Gallen

ncols = inGrid.size(2);
nrows = inGrid.size(1);

[yMinP,xMinP] = pix2latlon(inGrid.refmat,nrows,1);
[yMaxP,xMaxP] = pix2latlon(inGrid.refmat,1,ncols);

end

function outGrid = gridReSize(inGrid,xMinP,xMaxP,yMinP,yMaxP)
% gridReSize takes five inputs, an grid and maximum and minimum
% corrdinates you want to strink the grid to. The output is the resized
% grid.
% Created: 11/06/2013
% Author: Sean F. Gallen

ncols = inGrid.size(2);
nrows = inGrid.size(1);
R = inGrid.refmat;
[RSmaxR, RSminC] = UTMlatlon2pix(R,yMinP, xMinP);
[RSminR, RSmaxC] = UTMlatlon2pix(R,yMaxP, xMaxP);
GminR = 1;
GminC = 1;
GmaxC = ncols;
GmaxR = nrows;

%% This block of code resizes the grids according to smallest dimentions 
% and creates headers for the ascii file associated witht the resized grids

if RSminR >= GminR
    minR = RSminR;
else
    minR = GminR;
end

if RSmaxR <= GmaxR
    maxR = RSmaxR;
else
    maxR = GmaxR;
end

if RSminC >= GminC
    minC = RSminC;
else
    minC = GminC;
end

if RSmaxC <= GmaxC
    maxC = RSmaxC;
else
    maxC = GmaxC;
end
 outGrid = inGrid;
 
 minR = round(minR);
 maxR = round(maxR);
 minC = round(minC);
 maxC = round(maxC);
 
 rows = minR:maxR;
 cols = minC:maxC;
 outGrid.Z = inGrid.Z(rows,cols);
 outGrid.size = size(outGrid.Z);
 [yll,xll] = pix2UTMlatlon(inGrid.refmat,minR,minC);
 outGrid.refmat = makerefmat(xll, yll, inGrid.cellsize, -inGrid.cellsize);
 
end

function [lat, lon] = pix2UTMlatlon(R, row, col)

dx = R(2,1);
dy = R(1,2);

xll = R(3,1);
yll = R(3,2);

lat = row*dy+yll;
lon = col*dx+xll;
end

function [row, col] = UTMlatlon2pix(R, lat, lon)

dx = R(2,1);
dy = R(1,2);

xll = R(3,1);
yll = R(3,2);

row = (lat-yll)/dy;
col = (lon-xll)/dx;
end


%% This function calculates chi and sets of the forward operator matrix
function [A,Schi] = MakeChiGeoAMatrix(S,DA,GEO,mn)
%
% MakeChiGeoAMatrix.m will make an n by q matrix, A, that is n stream nodes
% long by q geologic map units. The matrix consists of the integral
% quantity chi traversed for each geologic map unit such that sum(A,2) [sum
% of all the rows] will equal chi for the stream network. This matrix can
% be used to calculate average ksn per geologic map unit using a matrix
% inversion of elevation For example, ksnPerGeo = A\Sz, where Sz is river network
% elevation and ksnPerGeo is a vector that is q long with the average ksn
% of each geologic map unit.
%
% Inputs:
% S         - TopoToolbox STREAMobj.
% DA        - Drainage area grid IN MAP UNITs (e.g. m^2) as a GRIDobj.
% GEO       - Indexed geologic map as a TopoToolbox GRIDobj.
% mn        - m/n or refence concavity (theta) value.
%
% Outputs:
% A         - A matrix as described above.
% Schi      - Chi for the stream network.
%
% Author: Sean F. Gallen
% Date modified: 06/08/2017
% email: sean.gallen{at}erdw.ethz.ch

p = inputParser;
p.FunctionName = 'MakeChiGeoAMatrix';
addRequired(p,'S', @(x) isa(x,'STREAMobj'));
addRequired(p,'A', @(x) isa(x,'GRIDobj'));
addRequired(p,'GEO', @(x) isa(x,'GRIDobj'));
addRequired(p,'mn', @(x) isscalar(x));

parse(p,S,DA,GEO,mn);

% get variables ready for chi integration
Schi = zeros(size(S.distance));
dSc = zeros(size(S.distance));
Six = S.ix;
Sixc = S.ixc;
Sx = S.distance;
Sa = DA.Z(S.IXgrid).^-mn;

h = waitbar(0,'calculating \chi for full stream network...');
% calculating chi and tau_star for the entire river network
for lp = numel(Six):-1:1
    Schi(Six(lp)) = Schi(Sixc(lp)) + (Sa(Sixc(lp))+(Sa(Six(lp))-Sa(Sixc(lp)))/2) *(abs(Sx(Sixc(lp))-Sx(Six(lp))));
    dSc(Six(lp)) =  Schi(Six(lp))- Schi(Sixc(lp));
    f = (numel(Six)+1 - lp)/numel(Six);
    waitbar(f,h);
end
close(h);

% set up variables for the A matrix
GridInd = S.IXgrid;
Sg = GEO.Z(GridInd);

geoID = unique(GEO.Z(~isnan(GEO.Z)));

N = length(Schi);
q = length(geoID);

A = zeros(N,q);

h = waitbar(0,'building the A matrix...');
for i = 1:q
    Ginds = find(Sg == geoID(i));
    A(Ginds,i) = dSc(Ginds);
    for lp = numel(Six):-1:1
        A(Six(lp),i) = A(Sixc(lp),i) + A(Six(lp),i);
    end
    waitbar(i/q,h)
end
close(h)
end

%% This function finds the turning point of the misfit function
function [idx,perpDist] = turingPointFinder(x,y)

% Two endpoints on the curve "data"
xend = [x(1) x(end)];
yend = [y(1) y(end)];
% The slope of the line connecting the two endpoints
m = ( yend(2) - yend(1) )/( xend(2) - xend(1) );
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
[~, idx] = max(perpDist);
end