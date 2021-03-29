function [sloSTREAMs] = binnedSlope(S,Sz,win,cs)
% binnedSlope.m calcuates slope from a STREAMobj and elevation data using 
% a moving window of a prescribed size
%
% Inputs:
%       1) S: a TopoToolBox STREAMobj
%       2) Sz: topologically ordered vector of elevation values
%       3) win: window size in map units
%       4) cs: cellsize in map units
%
% Outputs:
%       1) topologically ordered vector of slope values [m/m]
%
% Author: Sean F. Gallen
% Date Modified: 08/07/2017

step = (win/2);
indStep = ceil(step/cs);
% declare necessary variables from stream object
ordList = S.orderednanlist;         % ordered list of streams split by nans
strmBreaks = find(isnan(ordList));  % get position of nans
Six = S.ix;                         % doners
Sixc = S.ixc;                       % recievers
Sd = S.distance;                    % distance from mouth

% cast dumby vector to catch data
sloSTREAMs = nan(size(Sd));

% from the donor position identifiy receiver nodes in stream object order
receiverID = nan(size(Sd));
donorID = nan(size(Sd));
for i = numel(Six):-1:1;
    receiverID(Six(i)) = Sixc(i);
    donorID(Six(i)) = Sixc(i);
end


id1 = 0;
h = waitbar(0,'Calculating river slope in moving window...');
for i = 1:length(strmBreaks);
    strmInds = ordList(id1+1:strmBreaks(i)-1);
    tZ = Sz(strmInds);
    tD = Sd(strmInds);
    trecID = receiverID(strmInds);
    tSlo = nan(length(strmInds),1);
    % if the stream is a trunk channel run smoothing window
    if isnan(trecID(end))
        for j = 1:length(tD)
            % find within the window size.
            dWin = tD(tD > tD(j)-win & tD <= tD(j)+win);
            zWin = tZ(tD > tD(j)-win & tD <= tD(j)+win);
            
            % make sure there are enough points for the regression
            if length(dWin) > 2
                d_segMat = [ones(size(zWin)) dWin];
                [b,bint,r,rint,stats] = regress(zWin,d_segMat,0.05);
                tSlo(j) = b(2);
            else
                tSlo(j) = nan;
            end
        end
    % if the stream is a tributary allow smoothing window to continue down
    % stream by 1/2 of moving window width
    else
        addZ = nan(1,indStep);
        addD = nan(1,indStep);
        nInd = trecID(end);
        n = 1;
        while n < length(addZ);
            addZ(n) = Sz(nInd);
            addD(n) = Sd(nInd);
            nInd = receiverID(nInd);
            if isnan(nInd) == 1 || isempty(addD(isnan(addD))) == 1
                n = length(addD) + 1;
            else
                n = n+1;
            end
        end
        addZ = addZ(~isnan(addZ));
        addD = addD(~isnan(addD));
        tZ2 = nan(length(tZ)+length(addZ),1);
        tD2 = nan(length(tD)+length(addD),1);
        
        tZ2(1:length(tZ)) = tZ;
        tZ2(length(tZ)+1:end) = addZ;
        tD2(1:length(tD)) = tD;
        tD2(length(tD)+1:end) = addD;
        
        for j = 1:length(tD)
            % find within the window size.
            dWin = tD2(tD2 > tD2(j)-win & tD2 <= tD2(j)+win);
            zWin = tZ2(tD2 > tD2(j)-win & tD2 <= tD2(j)+win);
            
            % make sure there are enough points for the regression
            if length(dWin) > 2
                d_segMat = [ones(size(zWin)) dWin];
                [b,bint,r,rint,stats] = regress(zWin,d_segMat,0.05);
                tSlo(j) = b(2);
            else
                tSlo(j) = nan;
            end
        end
    end
    tSlo(tSlo < 0) = 0;
    sloSTREAMs(strmInds) = tSlo;
    id1 = strmBreaks(i);
    f = i/length(strmBreaks);
    waitbar(f,h);
end
close(h)
end
        
        
        
        
        
