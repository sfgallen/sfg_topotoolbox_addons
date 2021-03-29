function [x, y] = SWATHobj2XYoutline(SW)

    IM = ~isnan(SW.Z);
    B = bwboundaries(IM,4);
    for k = 1 : length(B)
        ix = sub2ind(size(IM),B{k}(:,1),B{k}(:,2));
        x = SW.X(ix);
        y = SW.Y(ix);
    end

end