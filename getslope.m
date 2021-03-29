function S = getslope(Z,nodeSize)
S = zeros(1,length(Z));
S(end) = abs(Z(end-1) - Z(end))/nodeSize;
% S(1) = abs(Z(1) - Z(2))/nodeSize;
% S(2:end-1) = abs(Z(1:end-2) - Z(3:end))./(2*nodeSize);
S(1:end-1) = abs(Z(1:end-1) - Z(2:end))/nodeSize;
end