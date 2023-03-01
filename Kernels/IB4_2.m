function w = IB4_2(r)
%Implementation of the standard IB4 kernel as derived by peskin
w = zeros(4,1);
q=sqrt(1+4*r*(1-r));
w(4,1)=(1+2*r-q)/8;
w(3,1)=(1+2*r+q)/8;
w(2,1)=(3-2*r+q)/8;
w(1,1)=(3-2*r-q)/8;

end

