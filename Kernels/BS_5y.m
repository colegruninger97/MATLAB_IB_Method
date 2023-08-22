function w = BS_5y(r)


%Input parameters
% r is the argument of the desired BSpline kernel
% r should be a real number in between 1 and zero

w = zeros(6,1);

if r<0 || r > 1
    msg = 'r must be between zero and 1.';
    error(msg);
end

w(1,1) = (189/80) - 153*abs(r+2)/40 + 99*((r+2)^2)/40 - 4*(abs(r+2)^3)/5 + (31/240)*(r+2)^4 - abs(r+2)^5/120;
w(2,1) = 119/240 + (17*abs(r+1))/24 - (77*(r+1)^2)/40 + (4*abs(r+1)^3)/3 - (31*(r+1)^4)/80 + abs(r+1)^5/24;
w(3,1) = 77/120 - (11*r^2)/20 + (31*r^4)/120 - r^5/12;
w(4,1) = 77/120 - (11*(1-r)^2)/20 + (31*(1-r)^4)/120 - ((1-r)^5)/12;
w(5,1) = 119/240 + (17*(2-r))/24 - (77*(2-r)^2)/40 + (4*(2-r)^3)/3 - (31*(2-r)^4)/80 + ((2-r)^5)/24;
w(6,1) = (189/80) - 153*(3-r)/40 + 99*((3-r)^2)/40 - 4*((3-r)^3)/5 + (31/240)*(3-r)^4 - ((3-r)^5)/120;
    
end