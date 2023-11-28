function [ffx,ffy] = spreadIB6(F,X,u,v,dx,dy,ds)

NIB = length(X(:,1));
[ru,cu] = size(u);
[rv,cv] = size(v);
ffx = zeros(ru,cu);
K = 59/60-sqrt(29)/20;
ffy = zeros(rv,cv);
C = (ds)/(dx*dy);

for k = 1:NIB
   sxx = X(k,1)/dx;
   syx = X(k,2)/dy - 1/2; %Need to subtract 1/2 to account for the staggering of the MAC grid
   sxy = X(k,1)/dx - 1/2;
   syy = X(k,2)/dy;
   sx = [sxx,syx];
   sy = [sxy,syy];
   Ixx = floor(sx);
   Iyy = floor(sy);
   rx = sx - Ixx; %compute the relevant r value defined in (0,1) in both x and y
   ry = sy - Iyy;
   %Get the relevant indices to interpolate onto the lagrangian structure
   %need to use modular arithmetic to account for periodic boundary
   %conditions
   
   %Need to include an if statement for the Bspline2 kernel implementation

   j1x = mod(Ixx(1)-2:Ixx(1)+3,cu)+1;
   i1x = mod(Ixx(2)-2:Ixx(2)+3,ru)+1;
   

   i1y = mod(Iyy(2)-2:Iyy(2)+3,rv)+1;
   j1y = mod(Iyy(1)-2:Iyy(1)+3,cv)+1;
   
   %get the weights
   wsx(:,:) = IB6x(rx(1),K).*IB6y(rx(2),K);
   wsy(:,:) = IB6x(ry(1),K).*IB6y(ry(2),K);
   
   ffx(i1x,j1x) = ffx(i1x,j1x) + C.*F(k,1).*wsx;
   ffy(i1y,j1y) = ffy(i1y,j1y) + C.*F(k,2).*wsy;
end

end

