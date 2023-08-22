function w = IB4_1(r);
%Implementation of the standard IB-4 kernel as derived in Peskin's original
%work
w = zeros(1,4);
q=sqrt(1+4*r*(1-r));
w(1,4)=(1+2*r-q)/8;
w(1,3)=(1+2*r+q)/8;
w(1,2)=(3-2*r+q)/8;
w(1,1)=(3-2*r-q)/8;

end


