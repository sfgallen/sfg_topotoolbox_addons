function line = selectLine(DEM,varargin)
    % line = selectLine(DEM,S)
    % selectLine will take a GRIDobj (required) and a STREAMobj
    % (optional) plot them and allow the user to interatively select a
    % line. The function outputs the x and y coordinates of the line in
    % increments of cellsize map units, distance along the profile and
    % elevation.
    %
    % Inputs
    % 1) DEM GRIDobj(DEM).
    % 2) STREAMobj (optional)
    % 3) x and y coordinates of pour points (Optional). matrix of n rows 
    %    by 2 columns where the column 1 is x and column 2 is y coordinates.
    %
    % Outputs
    % 1) 'line' as a structure with the following data properties
    %  - x --> vector of x coordinates of line
    %  - y --> vector of y coordinates of line
    %  - d --> distance along the line
    %  - z --> elevation along line
    %
    % Author: Sean F. Gallen 
    % Date Modified: 02/08/2017
    % email: sean.gallen[at]erdw.ethz.ch
    
    % test for arguments
    p = inputParser;
    p.FunctionName = 'selectLine';
    addRequired(p,'DEM', @(x) isa(x,'GRIDobj'));
    addOptional(p,'S',0, @(x) isa(x,'STREAMobj'));
    addOptional(p,'xy',[], @(x) isvector(x)| ismatrix(x));
    
    parse(p,DEM,varargin{:});
    
    % plot the data
    hfig = figure;
    imageschs(DEM); axis image
    hold on
    
    % if these is a STREAMobj plot that data on map
    if isa(p.Results.S, 'STREAMobj')
        plot(p.Results.S,'k-')
    end

    if ~isempty(p.Results.xy)
        plot(p.Results.xy(:,1),p.Results.xy(:,2),'o','MarkerFaceColor',...
            [1 1 1], 'MarkerEdgeColor', [0 0 0],'markersize', 4)
    end
    % Let the user select a line
    [x0,y0] = getline;
  
    % Use TopoToolbox function "getdistance" to get distance along line
    d0 = getdistance(x0,y0);
    
    % determine distance along line in units of cellsize using interpline
    [line.x,line.y,line.d] = interpline(x0,y0,d0,DEM.cellsize);
    plot(line.x,line.y,'k-','LineWidth',1.5);
    
    % get elevation
    line.z = interp(DEM,line.x,line.y);
end