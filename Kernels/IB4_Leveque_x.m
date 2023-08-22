function w = IB4_Leveque_x(r)

w = zeros(1,4);

if r<0 || r > 1
    msg = 'r must be between zero and 1.';
    error(msg);
end

w(1,1) = 2 - (3*(1+r)) + abs((1+r))^2;
w(1,2) = 1 - r^2;
w(1,3) = 1 - (1-r)^2;
w(1,4) = 2 - (3*(2-r)) + (2-r)^2;



end

