function w = sin_kerx(r)
%Implementation of the derivative of the peskin kernel 1/4(1+cos(pir/2)) to
%try and enforce neumann type boundary conditions

%Input parameters
%r needs to be an number in between zero and 1

% The kernel is supported at 4 points so we require that
w = zeros(4,1);
w(1,1) = -(pi/8)*sin(pi*(r-1)/2);
w(2,1) = -(pi/8)*sin(pi*(r)/2);
w(3,1) = -(pi/8)*sin(pi*(r+1)/2);
w(4,1) = -(pi/8)*sin(pi*(r+2)/2);

end










