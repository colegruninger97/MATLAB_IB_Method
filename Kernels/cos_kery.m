function w = cos_kery(r)
%Implementation of the derivative of the peskin kernel 1/4(1+cos(pir/2)) to
%try and enforce neumann type boundary conditions

%Input parameters
%r needs to be an number in between zero and 1

% The kernel is supported at 4 points so we require that
w = zeros(1,4);
w(1,1) = 0.25.*(1 + cos((r-1)*pi/2));
w(1,2) = 0.25.*(1 + cos((r)*pi/2));
w(1,3) = 0.25.*(1 + cos((r+1)*pi/2));
w(1,4) = 0.25.*(1 + cos((r+2)*pi/2));

end

