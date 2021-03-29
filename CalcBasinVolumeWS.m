function [basinVolume,volDEM,interpDEM] = CalcBasinVolumeWS(DEM,WS,varargin)
% CalcBasinVolume.m is a function that will take a DEM clipped to the
% boundaries of a drainage basin (DEM), a corresponding WATERSHEDobj (WS),
% and an option interpolation method to interpolate a surface over the
% watershed boundaries and calculate a volume between that surface and the
% drainage basin DEM.
%
% Inputs:
% 1) A DEM GRIDobj (required)
% 2) A WATERSHEDobj (required)
% 3) Interpolation method as a string (optional). Interpolation methods
% include: 'nearest', 'linear', 'natural', 'cubic', and 'v4'. The Matlab
% command 'griddata' is used for the interpolation, so type "help griddata"
% for more information regarding methods. 'cubic' is the default.
%
% Outputs:
% 1) basinVolume: scalar of volume of the basin in map units cubed.
% 2) volDEM: a GRIDobj of the volume calculation.
% 3) interpDEM: a GRIDobj of the interpolation surface.
%
% Author: Sean F. Gallen
% Date Modified: 02/08/2016
% Email: sean.gallen[at]erdw.ethz.ch

    % test for arguments
    p = inputParser;
    p.FunctionName = 'CalcBasinVolumeWS';
    validMethods  = {'nearest','linear','natural','cubic','v4'};
    addRequired(p,'DEM', @(x) isa(x,'GRIDobj'));
    addRequired(p,'WS', @(x) isa(x,'WATERSHEDobj'));
    addOptional(p,'method','cubic', @(x) ischar(validatestring(x,validMethods)));
    parse(p,DEM,WS,varargin{:});
    
    % get the data organized to use griddata
    r = DEM.size(1);
    c = DEM.size(2);
    [y,x] = coord2sub(DEM,WS.x,WS.y);
    ix = coord2ind(DEM,WS.x,WS.y);
    z = double(DEM.Z(ix));
    
    % meshgrid and griddata
    [xq, yq] = meshgrid(1:c,1:r);
    zq = griddata(x,y,z,xq,yq,p.Results.method);
    
    % set nan values from original DEM
    zq(isnan(DEM.Z)) = nan;
    
    % Make the interpolation surface GRIDobj
    interpDEM = DEM;
    interpDEM.Z = zq;
    
    % Make the volume GRIDobj
    volDEM = DEM;
    volDEM.Z = interpDEM.Z - DEM.Z;
    
    % Calculate the basin volume.
    basinVolume = nansum(volDEM.Z(:)).*volDEM.cellsize^2;
end
    
