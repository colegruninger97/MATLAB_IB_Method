function w = BS_1x(r)
%Implementation of the BSpline1 kernel function.

% Recall Bspline_1 looks like a hat function (1-x)Bs0(x-1/2) + (1+x)Bs0(x+1/2)

%Input parameters
% r is the argument of the desired BSpline kernel
% r should be a real number in between 1 and zero

if r<0 || r > 1
    msg = 'r must be between zero and 1.';
    error(msg);
end

w = zeros(1,2);
w(1) = 1-r;
w(2) = r;



end
