function [dataTable, attributes] = basinStatisticsFunction(DEM,FD,ix)
% basinStatisticsFunction takes an instance of a DEM GRIDobj (DEM), an 
% instance of a FLOWobj (FD) and drainage basin outlet locations as linear 
% indicies (ix) and will make a data table with the following statistics:
% 1) mean slope
% 2) standard deviation of slope
% 3) mean local relief
% 4) standard deviation of mean local relief
% 5) total basin relief
% 6) total basin drainage area
% 7) drainage basin volume
% 8) drainage basin volume-to-area ratio
% 9) drainage basin hypsometric integral

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
parse(p,DEM,FD,ix);

%% Set up a forloop
% this analysis will require a forloop and you will calculate each of the 9
% statistics above and fill in a table n-drainage basins long. 
%
% This means you will be fill in a table with n-rows by 9 columns. 
% This is your output. 
%
% Cast it as a table of NaNs before you enter the loop. 
% Use the command "nan" to do this, use the matlab help if you need it.
%
% Hint: the number of drainage basins can be determined from 'ix'.

%% declare dataTable variable
dataTable = nan(length(ix),9);

%%
% You will also need a drainage basin GRIDobj so use "drainagebasins" to
% make one here. 
% Call your drainage basin GRIDobj "DB"
DB = drainagebasins(FD,ix);

%% you will need to get the unique id of each drainage basin
DBid = unique(DB.Z(~isnan(DB.Z) & DB.Z > 0));

%% you will be performing the same routine n times so start the forloop
% based on the length of 'ix'
for i = 1:length(ix);
    
    % Most of the routine you will need is in
    % "Practical_4b_drainage_basins.m". Use that code, but keep in mind you
    % will need to make minor modifications to it!!!!!
    
    %% clip out drainage basin
    %
    % You will need to first clip out each drainage basin
    % to do this first get the drainage basin statistics for the basin.
    % Use the function "DrainageBasinWatershed" to do this.
    %
    % Hint: Make sure to index into the ix vector so that the loop knows
    % what drainage basin you are working with!!!!
    
    DBdata = DrainageBasinWatershed(FD,ix(i));
    
    % Look back at pratical 4b for how to make a clipped DEM of your basin
    
    %% make a copy of your DEM GRIDobj
    DEMc = DEM;
    
    %% in your copied DEM to nan that are not within the bounds of your 
    % drainage basin. Your drainage basin value is DBid(i).
    DEMc.Z(DB.Z ~=DBid(i)) = nan;
    
    %% use the TopoToolbox "crop" function
    DEMc = crop(DEMc,DBdata.x,DBdata.y);
    
    %% calculate basin statistics.
    % for the next part of the function just follow the basin statstics
    % section from Practical 4b, but remember that their are a few more
    % statistics you will need to calculate here.
    %
    % The first one is coded for you and the comments will tell you what
    % other statistics to input.
    
    %% 1) mean slope
    dbG = gradient8(DEMc,'degree');
    mslope = nanmean(dbG.Z(:));
    
    % add data to table
    dataTable(i,1) = mslope;
    
    %% 2) standard devation of slope
    stdS = nanstd(dbG.Z(:));
    
    % add data to table
    dataTable(i,2) = stdS;
    
    %% 3) mean local relief in a 500 m window
    lr = localtopography(DEMc,500);
    mlr = nanmean(lr.Z(:));
    
    % add data to table
    dataTable(i,3) = mlr;
    
    %% 4) standard devation of mean  local relief
    mlrs = nanstd(lr.Z(:));
    
    % add data to table
    dataTable(i,4) = mlrs;
    
    %% 5) total relief
    tbr = nanmax(DEMc.Z(:)) - nanmin(DEMc.Z(:));
    
    % add data to table
    dataTable(i,5) = tbr;
    
    %% 6) drainage basin area
    dbA = length(DEMc.Z(~isnan(DEMc.Z)))*DEMc.cellsize^2;
    
    % add data to table
    dataTable(i,6) = dbA;
    
    %% 7) drainage basin volume
    [dbV,vGrid,sGrid] = CalcBasinVolume(DEMc,DBdata.x,DBdata.y);
    
    % add data to table
    dataTable(i,7) = dbV;
    
    %% 8) drainage basin volume-to-area ratio
    
    % add data to table
    dataTable(i,8) = dbV/dbA;
    
    %% 9) drainage basin hypsometic integral
    
    % Note: I made a function for you called "basinHypsometry" that will
    % help you with this part. Use help for more information
    [Zbin, Abin, Znorm, Acum, HI] = basinHypsometry(DEMc);
    
    % add data to table
    dataTable(i,9) = HI;
    
    %% Create a cell array can be generated using cellarray = { }.
    % Each cell will be a string seperated by commas. 
    % Note: THESE STRINGS CANNOT CONTAIN SPACES OR '-' SIGNS.
    % Used underscores '_' of camel case.
    attributes = {'mean_slope','std_slope','mean_local_relief',...
        'std_mean_local_relief', 'total_relief', 'area', 'volume',...
        'volume_area_ratio','hypsometric_integral'};
    %% That should be it. Now you are ready to test the function!
end
end



