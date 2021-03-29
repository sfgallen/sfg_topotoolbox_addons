function selectCoordinates(DEM,varargin)
    % selectCoordinates will take a GRIDobj (required) and a STREAMobj
    % (optional) plot them and allow the user to interatively select x and
    % y points on the map and will export the coordinates of these points
    % Inputs
    % GRIDobj --> for example your DEM 
    % STREAMobj
    % Output
    % x --> vector of x coordinates
    % y --> vector of y coordinates
    %
    % Author: Sean F. Gallen 
    % Date Modified: 02/07/2017
    % email: sean.gallen[at]erdw.ethz.ch
    
    % test for arguments
    p = inputParser;
    p.FunctionName = 'selectCoordinates';
    addRequired(p,'DEM', @(x) isa(x,'GRIDobj'));
    addOptional(p,'S',0, @(x) isa(x,'STREAMobj'));
    addOptional(p,'line',0, @(x) isstruct(x));
    parse(p,DEM,varargin{:});
    
    % plot the data
    hfig = figure;
    imageschs(DEM); axis image
    hold on
    
    % if these is a STREAMobj plot that data on map
    if isa(p.Results.S, 'STREAMobj')
        plot(p.Results.S,'k-')
    end
    
    % if these is a STREAMobj plot that data on map
    if isstruct(p.Results.line)
        plot(p.Results.line.x,p.Results.line.y,'-','color',[.8 .8 .8])
    end

   % run the basic GUI to get user defined x and y coordinates
    
   basicGui(hfig);
end

function basicGui(hfig)
    % basic GUI to aid in selection of XY points

    figTB = uitoolbar(hfig);
    
    % make icon
    [rr cc] = meshgrid(1:19);
    CMat = ones(19,19,3);
    C = zeros(19,19);
    cirCMat = sqrt((rr-10).^2+(cc-10).^2);
    C(cirCMat >= 7) = 1;
    C(cirCMat <= 5) = 0.5;
    CMat(:,:,1) = C;
    CMat(:,:,2) = C;
    CMat(:,:,3) = C;
    
    % assign icon location on figure and call button function
    hpt = uipushtool(figTB,'CData',CMat,'TooltipString','pick xy point',...
        'UserData',struct('x',nan,'y',nan),'ClickedCallback',...
        @pickPts_button);

    uiwait(hfig)   
end

function pickPts_button(source, eventdata)
    % this function get x and y coordinates using ginput
    
    % get input value
    [x,y] = ginput(1);
    % If this is the first round start a new vector for x and y
    if isnan(source.UserData.x)
        source.UserData.x = x;
        source.UserData.y = y;
        plot(x,y,'o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',...
            [0 0 0],'MarkerSize', 4);
    % otherwise continue adding data
    else
        source.UserData.x = [source.UserData.x;x];
        source.UserData.y = [source.UserData.y;y];
        plot(x,y,'o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',...
            [0 0 0],'MarkerSize', 4);
    end
    
    data = struct('x',source.UserData.x,'y',source.UserData.y);
    % Save data to workspace
    assignin('base','userData',data);
end

