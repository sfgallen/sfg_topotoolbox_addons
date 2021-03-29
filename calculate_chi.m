function chi = calculate_chi(S,S_DA,mn)
    chi = zeros(size(S.distance));
    Six = S.ix;                         % donors
    Sixc = S.ixc;                       % recievers
    Sd = S.distance;                    % distance from mouth
    Sa = (1./(S_DA)).^mn;       % chi transformation variable
    
    h = waitbar(0,'calculating \chi...');
    % calculating chi for the entire river network
    for lp = numel(Six):-1:1
        chi(Six(lp)) = chi(Sixc(lp)) + (Sa(Sixc(lp))+(Sa(Six(lp))-Sa(Sixc(lp)))/2) *(abs(Sd(Sixc(lp))-Sd(Six(lp))));
        f = (numel(Six)+1 - lp)/numel(Six);
        waitbar(f,h);
    end
    close(h);
end