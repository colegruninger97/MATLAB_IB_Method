function w = BS_2x(r)
%Implementation of the BSpline1 kernel function.

%Input parameters
% r is the argument of the desired BSpline kernel
% r should be a real number in between 1 and zero

w = zeros(1,3);

if r<0 || r > 1
    msg = 'r must be between zero and 1.';
    error(msg);
end

if r < 0.5
    w(1,1) = 0.125 * (2*abs(1+r) - 3)*(2*abs(1+r) - 3);
    w(1,2) = 0.75 - r*r;
    w(1,3) = 0.125 *(2*abs(1-r) - 3)*(2*abs(1-r) - 3);
else
    w(1,1) = 0.125 * (2*abs(r) - 3)*(2*abs(r) - 3);
    w(1,2) = 0.75 - (1-r)*(1-r);
    w(1,3) = 0.125 * (2*abs(r-2) - 3)*(2*abs(r-2) - 3);
end


end
    
