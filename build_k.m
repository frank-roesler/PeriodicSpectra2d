function k = build_k(N)
%     Assembles vector of all possible combinations of frequencies between
%     -N and N such that the norms of k(:,i) are in ascending order.
    upper_length = ceil(pi*N^2);
    k = zeros(2, upper_length);
    ctr=1;
    for n = -N:N
        for m = -N:N
            k(:,ctr) = [n, m];
            ctr = ctr + 1;
        end
    end
    norms = sum(k.^2);
    k = k(:, norms < N^2);
    norms = sum(k.^2);
    [out, indices] = sort(norms);
    k = k(:, indices);
end
