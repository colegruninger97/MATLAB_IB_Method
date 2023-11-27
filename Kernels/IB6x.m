function w = IB6x(r,K)
%interpolation in the x direction using the IB6 kernel family parameterized
%by the constant K. 
w = zeros(1,6);

if r<0 || r > 1
    msg = 'r must be between zero and 1.';
    error(msg);
end

%first compute the constants beta and gamma
beta = (9/4) - (3/2)*(K + r^2) + ((22/3) - 7*K)*r - (7/3)*r^3;
gamma = -(11/32)*r^2 + (3/32)*(2*K + r^2)*r^2 + (1/72)*((3*K-1)*r + r^3)^2 + (1/18)*((4-3*K)*r - r^3)^2;
w(1,6) = (1/56)*(-beta + sign(3/2 - K)*sqrt(beta*beta - 112*gamma));
w(1,5) = -3*w(1,6) - 1/16 + (K+r^2)/8 + ((3*K-1)*r)/12 + (r^3)/12;
w(1,4) = 2*w(1,6) + 1/4 + (4-3*K)*r/6 - r^3/6;
w(1,3) = 2*w(1,6) + 5/8 - (K+r^2)/4;
w(1,2) = -3*w(1,6) + 1/4 - (4-3*K)*r/6 + r^3/6;
w(1,1) = w(1,6) - 1/16 + (K+r^2)/8 - (3*K-1)*r/12 - r^3/12;

end

