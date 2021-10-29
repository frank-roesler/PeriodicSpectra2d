function val = potential(x,y)
%   Potential function
    val = 10i*cos(2*pi*x).*sin(2*pi*y) + 9*sin(2*pi*y);
end