function bdry = bdry_pts(window, h_L)
    bdry = [];
    for n=1:length( window)
        pt=window(n);
        neighbors = window(abs(pt-window)<=1.1*h_L);
        if length(neighbors)<5
            bdry=unique([bdry,pt]);
        end
    end
end