function [ffx,ffy] = spreadIB4(F,X,u,v,dx,dy,ds)

NIB = length(X(:,1));
[ru,cu] = size(u);
[rv,cv] = size(v);
ffx = zeros(ru,cu);
ffy = zeros(rv,cv);
C = (ds)/(dx*dy);


for k = 1:NIB
   sxx = X(k,1)/dx;
   syx = X(k,2)/dy - 1/2; %Need to add 1/2 to account for the staggering of the grid
   sxy = X(k,1)/dx - 1/2;
   syy = X(k,2)/dy;
   sx = [sxx,syx];
   sy = [sxy,syy];
   Ixx = floor(sx);
   Iyy = floor(sy);
   rx = sx - Ixx;
   ry = sy - Iyy;
   %Get the relevant indices to interpolate onto the lagrangian structure
   %need to use modular arithmetic to account for periodic boundary
   %conditions
   
   j1x = mod(Ixx(1)-1:Ixx(1)+2,cu)+1;
   i1x = mod(Ixx(2)-1:Ixx(2)+2,ru)+1;
   j1y = mod(Iyy(1)-1:Iyy(1)+2,cv)+1;
   i1y = mod(Iyy(2)-1:Iyy(2)+2,rv)+1;
   
   %get the weights
   wsx(:,:) = IB4_1(rx(1)).*IB4_2(rx(2));
   wsy(:,:) = IB4_1(ry(1)).*IB4_2(ry(2));
   
   ffx(i1x,j1x) = ffx(i1x,j1x) + C.*F(k,1).*wsx;
   ffy(i1y,j1y) = ffy(i1y,j1y) + C.*F(k,2).*wsy;
end




end

