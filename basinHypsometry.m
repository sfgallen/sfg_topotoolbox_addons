function [Zbin, Abin, Znorm, Acum, HI] = basinHypsometry(DEM,varargin)
% basinHypsometry will calculate the hypsometry, hypsometric curve and
% hysometric integral for a DEM that is clipped to the extent of adrainage
% basin. The function using a clipped DEM as a GRIDobj and an optional 
% binning increment (the default is 50 meters) and outputs 4 vectors. 
% Vectors 1 and 2 are the hypsometry elevation and land area data,
% respectively. Vectors 3 and 4 and the normalized elevation and normalized
% cumulative sum of land area per unit elevation, respectively, for the
% hypsometric curve. The final output it the hypsometric integral, which is
% the integral of the hypsometric curve.
%
% example:
% [Zbin, Abin, Znorm, Acum, HI] = basinHypsometry(DEM,binInc)
%
% Inputs:
% 1) DEM GRIDobj clipped to drainage basin extent(DEM) or as a vector of 
%    elevation values from a clipped DEM (required)
% 2) elevation binning increment (binInc) (option, default 50 meters)
%
% Outputs:
% 1) binned elevation data (Zbin) (vector)
% 2) binned land area per unit elevation (Abin) (vector)
% 3) normalized binned elevation (Znorm) (vector)
% 4) normalized cumulative sum of area per unit elevation (Acum) (vector)
% 5) hypsometric integral (HI) (scalar)
%
% To plot hypsometry one can use bar(Zbin,Abin) or other plotting function.
% To plot hypsometic curve plot(Acum,Znorm) will work.
%
% Author: Sean F. Gallen
% Date Modified: 02/09.2017
% email: sean.gallen[at]erdw.ethz.ch

    % test for arguments
    p = inputParser;
    p.FunctionName = 'basinHypsometry';
    addRequired(p,'DEM', @(x) isa(x,'GRIDobj') | isvector(x));
    addOptional(p,'binInc',50, @(x) isscalar(x));
    parse(p,DEM,varargin{:});

    if isa(DEM,'GRIDobj')
        % index drainage basin to get elevation data
        WSz = DEM.Z(~isnan(DEM.Z));
    else
        WSz = DEM;
    end

    % determine range and binning increment and make your binning incement vector
    maxZ = ceil(max(WSz)/10)*10;
    minZ = floor(min(WSz)/10)*10;
    binInc = p.Results.binInc;

    binVect = minZ:binInc:maxZ;

    % now declare your variables before starting your forloop.
    Zbin = zeros(length(binVect)-1,1);
    Abin = zeros(length(binVect)-1,1);
    cs = DEM.cellsize;

    % now start you forloop
    for i = 1:length(Zbin);
        Zbin(i) = (binVect(i)+binVect(i+1))/2;
        Abin(i) = length(WSz(WSz >= binVect(i) & WSz < binVect(i+1)))*cs^2;
    end


    % Normalize elevatio vector for hypsometric curve
    Znorm = (Zbin-min(Zbin))/(max(Zbin) - (min(Zbin)));

    % calculate cumulative sum of area per unit elevation and normalize
    Acum = cumsum(Abin);
    Acum = Acum/max(Acum);

    % Integrate under the hypsometric curve to get the hysometric integral
    HI = trapz(Acum,Znorm);
    
%     %plot the data
%     figure
%     subplot(1,2,1)
%     stairs(Zbin-binInc/2,Abin,'b-','LineWidth',2); axis tight
%     xlabel('elevation (m)'); ylabel('land area (m^2)');
%     title('hypsometry');
%     subplot(1,2,2);
%     plot(Acum,Znorm,'r-','LineWidth',2)
%     xlabel('normalized elevation'); ylabel('normalized cumulative area');
%     title('hypsometric curve');
%     legend(['HI = ',num2str(HI)]);
end

