function MS = LSGrid2ShapeFile(LS_grid, table, attributes,file_name)

% LSGrid2ShapeFile takes an instance of a integer grid where polygons in
% the grid are defined as postive integers surrounded by NaN or null values
% and will generate a mapstructure with data for each polygon based on an 
% input data table and list of attributes as strings in a cell array. The 
% map structure is then used to generate an ESRI shapefile
%
% This function is simply a modified version of the TopoToolbox function
% GRIDobj2polygon by Wolfgang Schwanghart 
% (w.schwanghart[at]geo.uni-potsdam.de)
%
% example: 
% MP = LSGrid2ShapeFile(LSgrid, table, attributes, 'your_file_name.shp')
%
% Inputs (all required): 
% 1) polygon grid as a matrix
% 2) Table with attribute data. The table should be organized such
% that each column is a different attribute or statistic and each row
% represents a drainage basin.
% 3) List of strings as a cell array to give headers to map structure
% attributes.
% 4) Name of file with '.shp' file extension
%
% Outputs:
% 1) polygon mapstructure.
% 2) shapefile saved to current directory
%
% Author: Sean F. Gallen
% Date Modified: 10/02/2019
% email: sean.gallen@colostate.edu


% test for arguments
p = inputParser;
p.FunctionName = 'LSGrid2ShapeFile';
addRequired(p,'LG_grid', @(x) isa(x,'GRIDobj'));
addRequired(p,'table', @(x) isvector(x) | ismatrix(x));
addRequired(p,'attributes', @(x) iscell(x));
addRequired(p,'file_name', @(x) ischar(x));
parse(p,LS_grid, table, attributes,file_name);

% generate basin ID numbers
dbid = unique(LS_grid.Z(~isnan(LS_grid.Z) & LS_grid.Z > 0));

% determine size of the inputs
[nbasins, nattributes] = size(table);
caSize = max(size(attributes));

% some extra error handling
if length(dbid) ~= nbasins
    error(['Number of basins in table and drainage basin GRIDobj ',...
        'do not match. Make sure the number of rows in the table ',...
        'is the same as the number of basins in the GRIDobj.']);
elseif nattributes ~= caSize
    error(['Number of attributes in the table does not match the ',...
        'number of attribute headers. Make sure that the number of'...
        ' columns in your table is equal to the number of ',...
        'attribute headers'])
end

% check underlying class of the grid
if isfloat(LS_grid.Z)
    writevalue = true;
    LS_Grid2 = GRIDobj(LS_grid,'uint32');
    I   = ~(isnan(LS_grid.Z) | LS_grid.Z == 0);
    [uniquevals,~,LS_Grid2.Z(I)] = unique(LS_grid.Z(I));
    LS_grid  = LS_Grid2;
else
    writevalue = false;
end

% identify regions and number of regions
STATS = regionprops(uint32(LS_grid.Z),'Area','PixelIdxList');
ndb = numel(STATS);  
    

for r = 1:ndb
    
    % get subscripts of basin
    [row,col] = ind2sub(LS_grid.size,STATS(r).PixelIdxList);
    % bw2poly returns the coordinates of the boundary
    C   = bw2poly([row col],0,0);
    
    % add nans at the end of coordinate vectors
    if numel(C)>1
        C = cellfun(@(x) [x;[nan nan]],C,'UniformOutput',false);
    end
    C = cell2mat(C);
    
    % write data to mapstruct
    MS(r).Geometry = 'polygon';
    [x,y] = sub2coord(LS_grid,C(:,1),C(:,2));
    MS(r).X = x;
    MS(r).Y = y;
    MS(r).ID = double(LS_grid.Z(STATS(r).PixelIdxList(1)));
    
    for j = 1:nattributes
        MS(r).(attributes{j}) = table(r,j);
    end
    
    if writevalue
        MS(r).gridval = double(uniquevals(r));
    end
    
end

% create coordinate vectors if more than one output
if nargout > 1
    for r=1:numel(MS)
        if ~isnan(MS(r).X(end))
            MS(r).X(end+1) = nan;
            MS(r).Y(end+1) = nan;
        end
    end
    
    x = {MS.X}';
    y = {MS.Y}';
    x = cell2mat(x);
    y = cell2mat(y);
end


shapewrite(MS,file_name);
end


function C = bw2poly(BW,mp,holes)

if islogical(BW)
    [r,c] = find(BW);
else
    r = BW(:,1);
    c = BW(:,2);
end
    
rc = bsxfun(@plus,r,[0 -.5 0 .5 .5 .5 0 -.5 -.5]);
cc = bsxfun(@plus,c,[0 -.5 -.5 -.5 0 .5 .5 .5 0 ]);
rc = rc(:);
cc = cc(:);

rc = rc*2;
cc = cc*2;


minrc = min(rc)-1;
maxrc = max(rc)-1;
if mp
    mincc = min(cc)-1;
else
    [mincc,ix] = min(cc);
    mincc = mincc-1;
end
maxcc = max(cc)-1;

rc = rc-minrc;
cc = cc-mincc;

siz = [maxrc-minrc+1 maxcc-mincc+1];

IX = sub2ind(siz,rc,cc);
B  = false(siz);
B(IX) = true;

if mp 
    if ~holes
        C = bwboundaries(B,4,'noholes');
        C = cellfun(@(x) modifynodelist(x),C,'UniformOutput',false);
    else
        [C,~,N] = bwboundaries(B,4,'holes');
        for iter2=1:numel(C)
            if iter2<=N
                C{iter2} = modifynodelist(C{iter2});
            else
                C{iter2} = bsxfun(@plus,C{iter2},[minrc mincc]);
                C{iter2} = C{iter2}/2;
            end
        end
            
    end
    C = C(:);
else
    C = bwtraceboundary(B,[rc(ix) cc(ix)],'S',4,inf,'counterclockwise');
    C = {modifynodelist(C)};
end
    


function nl = modifynodelist(nl)
% modify node list so that only nodes remain on pixel corners 
nl(any(mod(nl,2)==0,2),:) = [];
nl = bsxfun(@plus,nl,[minrc mincc]);
nl = nl/2;
end
end
