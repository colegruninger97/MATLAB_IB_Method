function w = BS_4x(r)


%Input parameters
% r is the argument of the desired BSpline kernel
% r should be a real number in between 1 and zero

w = zeros(1,5);

if r<0 || r > 1
    msg = 'r must be between zero and 1.';
    error(msg);
end


if r <= 0.5
    w(1,1) = (1/384)*(5 - 2*abs(2+r))^4;
    w(1,2) = 55/96 + (5/24)*abs(1+r) - (5/4)*(1+r)^2 + (5/6)*abs(1+r)^3 - (1/6)*(1+r)^4;
    w(1,3) = 115/192 - (5/8)*r*r + 0.25*r*r*r*r;
    w(1,4) = 55/96 + (5/24)*abs(1-r) - (5/4)*abs(1-r)^2 + (5/6)*abs(1-r)^3 - (1/6)*abs(1-r)^4;
    w(1,5) = (1/384)*(5 - 2*abs(2-r))^4;
else
    w(1,1) = (1/384)*(5 - 2*abs(1+r))^4;
    w(1,2) = 55/96 + (5/24)*abs(r) - (5/4)*(r)^2 + (5/6)*abs(r)^3 - (1/6)*(r)^4;
    w(1,3) = 115/192 - (5/8)*abs(1-r)^2 + 0.25*(1-r)^4;
    w(1,4) = 55/96 + (5/24)*abs(2-r) - (5/4)*abs(2-r)^2 + (5/6)*abs(2-r)^3 - (1/6)*abs(2-r)^4;
    w(1,5) = (1/384)*(5 - 2*abs(3-r))^4;
    
end
    
end



