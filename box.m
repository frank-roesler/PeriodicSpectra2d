function y = box(Lattice, Subset, r_outer, r_inner)
%     Computes a neighborhood of lattice points of radius r around a given
%     Subset.
    y = zeros(size(Lattice));
    for z = Subset
        y(abs(Lattice - z)<r_outer & abs(Lattice - z)>r_inner)=1;
    end
    y = logical(y);
end