function [dataTable, attributes] = basinStatisticsFunction_plus_erosion(DEM,FD,ix,mn,critA,Ao,BeC,BeCstd,utm_zone,tag)
% basinStatisticsFunction_plus_erosion.m takes ten inputs:
%
% Inputs:
% 1) a DEM as a TopoToolbox GRIDobj projected into UTM
% 2) a flow direction grid as a TopoToolbox FLOWobj
% 3) the linear indices into the DEM of basin outlets as a vector or scalar
% 4) the reference concavity or m/n value for river profile analysis
% 5) the critical drainage area for stream channel initiation.
% 6) the reference drainage area for chi analysis (use 1)
% 7) A vector of Be10 concentrations associated with each basin outlet
% 8) A vector of the 1-sigma uncertainties associated with the Be10 data
% 9) The UTM grid zone number
% 10) A string used to define the output files.
%
% The code will make a data table with the following statistics:
% 1) mean slope
% 2) standard deviation of slope
% 3) mean local relief
% 4) standard deviation of mean local relief
% 5) total basin relief
% 6) total basin drainage area
% 7) drainage basin volume
% 8) drainage basin volume-to-area ratio
% 9) drainage basin hypsometric integral
% 10) theta from the slopearea function
% 11) ks from the slope area function
% 12) best-fit m/n from minimum variance in ChiFits
% 13) ks_var from ChiFits
% 14) ks_var uncertainty
% 15) mn_ref from ChiFits
% 16) ksn from ChiFits
% 17) ksn uncertainty. 
% 18) Production spallation
% 19) Production slow muons
% 20) Production fast muons
% 21) erosion rate mm/yr
% 22) plus 1 std erosion rate
% 23) minus 1 std erosion rate

% You will have two outputs. (1) is the data table as a matrix and (2) a
% cell array of strings that identifies each statistic in the same order as
% the columns of the statistics table.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up your data inputs properly and fix the help section so it is more
% informative. The inputs should be 'DEM', 'FD' and 'ix'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% don't worry about this part. It isn't necessary, but will stop the
% function if the inputs are not correct
% test for arguments
p = inputParser;
p.FunctionName = 'basinStatisticsFunction';
addRequired(p,'DEM', @(x) isa(x,'GRIDobj'));
addRequired(p,'FD', @(x) isa(x,'FLOWobj'));
addRequired(p,'ix', @(x) isscalar(x) | isvector(x));
addRequired(p,'mn', @(x) isscalar(x));
addRequired(p,'critA', @(x) isscalar(x));
addRequired(p,'Ao', @(x) isscalar(x));
addRequired(p,'BeC', @(x) isscalar(x) | isvector(x));
addRequired(p,'BeCstd', @(x) isscalar(x) | isvector(x));
addRequired(p,'utm_zone', @(x) isscalar(x));
addRequired(p,'tag', @(x) ischar(x));
parse(p,DEM,FD,ix,mn,critA,Ao,BeC,BeCstd,utm_zone,tag);

%% Set up a forloop
% this analysis will require a forloop and you will calculate each of the 9
% statistics above and fill in a table n-drainage basins long. This means
% you will be fill in a table with n-rows by 9 columns. This is your
% output. Cast it as a table of NaNs before you enter the loop. Use the
% command "nan" to do this, use the matlab help if you need it.
% Hint: the number of drainage basins can be determined from 'ix'.

dataTable = nan(length(ix),23);


% you will be performing the same routine n times so start the forloop
% based on the length of 'ix'

for i = 1:length(ix)
    
    % Most of the routine you will need is in
    % "Practical_4b_drainage_basins.m". Use that code, but keep in mind you
    % will need to make minor modifications to it!!!!!
    
    %% clip out drainage basin
    % You will need to first clip out each drainage basin
    % to do this first get the drainage basin statistics for the basin.
    % Use the function "DrainageBasinWatershed" to do this.
    % Hint: Make sure to index into the ix vector so that the loop knows
    % what drainage basin you are working with!!!!
    DB = drainagebasins(FD,ix(i));
    DBdata = DrainageBasinWatershed(FD,ix(i));
    
    % Look back at pratical 4b for how to make a clipped DEM of your basin
    
    % make a copy of your DEM GRIDobj
    DEMc = DEM;
    
    % in your copied DEM to nan that are not within the bounds of your 
    % drainage basin. Your drainage basin value is DBid(i).
    DEMc.Z(DB.Z ~= 1) = nan;
    
    % use the TopoToolbox "crop" function
    DEMc = crop(DEMc,DBdata.x,DBdata.y);
    
    figure
    imageschs(DEMc);
    
    %% calculate basin statistics.
    % for the next part of the function just follow the basin statstics
    % section from Practical 4b, but remember that their are a few more
    % statistics you will need to calculate here.
    %
    % The first one is coded for you and the comments will tell you what
    % other statistics to input.
    
    % 1) mean slope
    dbG = gradient8(DEMc,'degree');
    mslope = nanmean(dbG.Z(:));
    
    % add data to table
    dataTable(i,1) = mslope;
    
    % 2) standard devation of slope
    stdS = nanstd(dbG.Z(:));
    % add data to table
    dataTable(i,2) = stdS;
    
    % 3) mean local relief in a 500 m window
    lr = localtopography(DEMc,500);
    mlr = nanmean(lr.Z(:));
    % add data to table
    dataTable(i,3) = mlr;
    
    % 4) standard devation of mean  local relief
    mlrs = nanstd(lr.Z(:));
    % add data to table
    dataTable(i,4) = mlrs;
    
    % 5) total relief
    tbr = nanmax(DEMc.Z(:)) - nanmin(DEMc.Z(:));
    % add data to table
    dataTable(i,5) = tbr;
    
    % 6) drainage basin area
    dbA = length(DEMc.Z(~isnan(DEMc.Z)))*DEMc.cellsize^2;
    % add data to table
    dataTable(i,6) = dbA;
    
    % 7) drainage basin volume
    [dbV,vGrid,sGrid] = CalcBasinVolume(DEMc,DBdata.x,DBdata.y);
    % add data to table
    dataTable(i,7) = dbV;
    
    % 8) drainage basin volume-to-area ratio
    
    % add data to table
    dataTable(i,8) = dbV/dbA;
  
    % 9) drainage basin hypsometic integral
    % Note: I made a function for you called "basinHypsometry" that will
    % help you with this part. Use help for more information
    [Zbin, Abin, Znorm, Acum, HI] = basinHypsometry(DEMc,10);
    % add data to table
    dataTable(i,9) = HI;
    
    %% prep data for stream metrics
    FDc = FLOWobj(DEMc);
    dbA = flowacc(FDc);
    dbS = STREAMobj(FDc,'minarea',critA/DEMc.cellsize^2);
    
    % 10) theta from the slopearea function
    SA = slopearea(dbS,DEMc,dbA,'plot',false);
    dataTable(i,10) = -SA.theta;
    
    % 11) ks from the slope area function
    dataTable(i,11) = SA.ks;
    
    % 12) best-fit m/n from minimum variance in ChiFits
    C = ChiFits(dbS,DEMc,dbA,'a0',Ao,'mnplot',false,'mn',mn,'plot',false);
    dataTable(i,12) = C.mn_var(1);
    
    % 13) ks_var from ChiFits
    dataTable(i,13) = C.ks_var(1);
    
    % 14) ks_var uncertainty
    dataTable(i,14) = C.ks_95uc_var(1);
    
    % 15) mn_ref from ChiFits
    dataTable(i,15) = C.mnref;
    
    % 16) ksn from ChiFits
    dataTable(i,16) = C.ksn;
    
    % 17) ksn uncertainty.
    dataTable(i,17) = C.ksn_95uc;
    
    %% calculations for erosion rate
    Production = OldMethodProduction_sfg_mod(DEMc,utm_zone,tag);
    Erosion = OldMethodErosion(BeC(i),BeCstd(i),Production);
    dataTable(i,18) = nanmean(Production.Pn(:));
    dataTable(i,19) = nanmean(Production.Pms(:));
    dataTable(i,20) = nanmean(Production.Pmf(:));
    dataTable(i,21) = Erosion.Denudation_mmYr;
    dataTable(i,22) = Erosion.Denudation_UpError;
    dataTable(i,23) = Erosion.Denudation_DownError;
    
    % Create a cell array can be generated using cellarray = { }.
    % Each cell will be a string seperated by commas. 
    % Note: THESE STRINGS CANNOT CONTAIN SPACES OR '-' SIGNS.
    % Used underscores '_' of camel case.
    attributes = {'mean_slope','std_slope','mean_local_relief',...
        'std_mean_local_relief', 'total_relief', 'area', 'volume',...
        'volume_area_ratio','hypsometric_integral','theta_sa','ks_sa',...
        'mn_chi_fit','ks_chi','ks_95uc','mn_ref','ksn_chi','ksn_95uc',...
        'Pn','Pms','Pmf','ero','ero_std_p1','ero_std_m1'};
    %% That should be it. Now you are ready to test the function!
end
end



