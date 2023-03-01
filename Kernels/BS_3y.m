function w = BS_3y(r)
%Implementation of the BSpline1 kernel function.

%Input parameters
% r is the argument of the desired BSpline kernel
% r should be a real number in between 1 and zero

w = zeros(4,1);

if r<0 || r > 1
    msg = 'r must be between zero and 1.';
    error(msg);
end

w(1,1) = (1/6)*(2 - abs(r+1)).^3;
w(2,1) = (2/3) - (r)*(r) + 0.5*(abs(r)).^3;
w(3,1) = (2/3) - (1-r)*(1-r) + 0.5*abs(1-r).^3;
w(4,1) = (1/6)*(2 - abs(2-r)).^3;


end

