function w = BS_0y(r)
%Implementation of the BSpline1 kernel function.

%Input parameters
% r is the argument of the desired BSpline kernel
% r should be a real number in between 1 and zero

w = zeros(1,1);

if r<0 || r > 1
    msg = 'r must be between zero and 1.';
    error(msg);
end

w(1,1) = 1;


end

