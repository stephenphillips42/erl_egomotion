function [ hatx ] = hat( x )
%HAT Create skew-symmetric matrix from x in R^3
hatx = [ 0, -x(3), x(2); x(3), 0, -x(1); -x(2), x(1), 0 ];
end
