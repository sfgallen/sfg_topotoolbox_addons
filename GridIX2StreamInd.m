function ind = GridIX2StreamInd(S,IX)
ind =nan(length(IX),1);
gridIX = S.IXgrid;
for i = 1:length(IX)
    ind(i) = find(gridIX == IX(i));
end
end

