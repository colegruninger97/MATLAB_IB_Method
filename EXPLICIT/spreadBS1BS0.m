function [ffx,ffy] = spreadBS1BS0(F,X,u,v,dx,dy,ds)

%The goal of this funciton is to spread forces on the Lagrangian structure
%to the Eulerian grid
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
   
   %Need to include an if statement for the Bspline2 kernel implementation
   j1x = mod(Ixx(1):Ixx(1)+1,cu) + 1;
   if rx(2) < 0.5
       i1x = mod(Ixx(2),ru) + 1;
   else
       i1x = mod(Ixx(2)+1,ru) + 1;
   end
   
   if ry(1) < 0.5
       j1y = mod(Iyy(1),cv) + 1;
   else
       j1y = mod(Iyy(1)+1,cv) + 1;
   end
  
   i1y = mod(Iyy(2):Iyy(2)+1,rv)+1;
   
   %get the weights
   wsx(:,:) = BS_1x(rx(1)).*BS_0y(rx(2));
   wsy(:,:) = BS_0x(ry(1)).*BS_1y(ry(2));
   
   ffx(i1x,j1x) = ffx(i1x,j1x) + C.*F(k,1).*wsx;
   ffy(i1y,j1y) = ffy(i1y,j1y) + C.*F(k,2).*wsy;
end


end