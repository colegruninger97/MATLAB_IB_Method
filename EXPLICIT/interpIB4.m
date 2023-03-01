function [U,V,i1x,j1x,i1y,j1y] = interpBS2BS1(u,v,X,dx,dy)
%Right now doing this for the BSpline2,3 couple. Other kernels will be
%incorporated later

%Also, this is currently implemented for periodic boundary conditions only,
%I will update to more general physical boundary conditions later...

NIB = length(X(:,1));
[ru,cu] = size(u);
[rv,cv] = size(v);

U = zeros(NIB,1);
V = zeros(NIB,1);

for k = 1:NIB
   sxx = X(k,1)/dx;
   syx = X(k,2)/dy + 1/2; %Need to add 1/2 to account for the staggering of the grid
   sxy = X(k,1)/dx + 1/2;
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

   j1x = mod(Ixx(1)-1:Ixx(1)+2,cu)+1;
   i1x = mod(Ixx(2)-1:Ixx(2)+2,ru)+1;
   

   i1y = mod(Iyy(2)-1:Iyy(2)+2,rv)+1;
   j1y = mod(Iyy(1)-1:Iyy(1)+2,cv)+1;
   
   %get the weights
   wsx(:,:) = IB4_1(rx(1)).*IB4_2(rx(2));
   wsy(:,:) = IB4_1(ry(1)).*IB4_2(ry(2));
   U(k,1) = sum(sum(u(i1x,j1x).*wsx(:,:)));
   V(k,1) = sum(sum(v(i1y,j1y).*wsy(:,:)));
    
end




end

