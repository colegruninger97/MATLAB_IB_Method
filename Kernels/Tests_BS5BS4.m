% Tests for the BS4BS4 composite kernel

%Test scripts to ensure that kernels are indeed producing local divergence
%free interpolations
clear;
clc;
u_b = @(x,y) sin(2*pi.*x+2).*sin(2*pi.*y+4);
v_b = @(x,y) cos(2*pi.*x+2).*cos(2*pi.*y+4);

N = 64;
dx = 1/N;
dy = dx;

x = 0:dx:1-dx;
y = 0:dy:1-dy;

y_c = dy/2:dy:1-dy/2;
x_c = dx/2:dx:1-dx/2;


[xx_sidex,yy_sidex] = meshgrid(x,y_c);
[xx_sidey,yy_sidey] = meshgrid(x_c,y);


u = u_b(xx_sidex,yy_sidex);
v = v_b(xx_sidey,yy_sidey);

%compute the discrete divergence
du = (u(:,2:end)-u(:,1:end-1))./dx;
dv = (v(2:end,:)-v(1:end-1,:))./dy;
div = du(1:end-1,:)+dv(:,1:end-1);

%Now interpolate points using the composite B-spline kernel combinations
%generate a million randomly sampled points on the unit square
delta = 1e-4;

%Pick a random point in the unit square
xr = rand*(1-dx);
yr = rand*(1-dy);
% xr = 0.95;
% yr = 0.95;
scatter(xr,yr,'s')
hold on
plot(xx_sidey,yy_sidey,'ko')
plot(xx_sidex,yy_sidex,'k*')


%Using this point, create a new point that is delta away from it...

xr2 = xr + delta;
xr3 = xr - delta;

yr2 = yr + delta;
yr3 = yr - delta;
% scatter(xr,yr2,'s')
% scatter(xr,yr3,'s')
% scatter(xr2,yr,'s')
% scatter(xr3,yr,'s')
xs = [xr3 xr xr2];
ys = [yr3 yr yr2];
u_hat = zeros(3,1);
v_hat = zeros(3,1);

for i = 1:3

%interpolate at the randomly generated point
%add 1 to each point to account for the fact that MATLAB indexes by zero...
sxx = xs(i)/dx;
syx = ys(2)/dy - 1/2;
sxy = xs(2)/dx - 1/2;
syy = ys(i)/dy;
sx = [sxx,syx];
sy = [sxy,syy];

Ixx = floor([sx]);
Iyy = floor([sy]);

rx = sx-Ixx;
ry = sy-Iyy;

    
   j1x = mod(Ixx(1)-2:Ixx(1)+3,N)+1;
   if rx(2) < 0.5
       i1x = mod(Ixx(2)-2:Ixx(2)+2,N)+1;
   else
       i1x = mod(Ixx(2)-1:Ixx(2)+3,N)+1;
   end
   
   if ry(1) < 0.5
       j1y = mod(Iyy(1)-2:Iyy(1)+2,N)+1;
   else
       j1y = mod(Iyy(1)-1:Iyy(1)+3,N)+1;
   end
   
   i1y = mod(Iyy(2)-2:Iyy(2)+3,N)+1;

   %get the weights
   wsx(:,:) = BS_5x(rx(1)).*BS_4y(rx(2));
   wsy(:,:) = BS_4x(ry(1)).*BS_5y(ry(2));
   u_hat(i) = sum(sum(u(i1x,j1x).*wsx(:,:)));
   v_hat(i) = sum(sum(v(i1y,j1y).*wsy(:,:)));
end
div_interp = (u_hat(1)-u_hat(3))/(2*delta) + (v_hat(1)-v_hat(3))/(2*delta)

  plot(xx_sidey(i1y,j1y),yy_sidey(i1y,j1y),'o')
  plot(xx_sidex(i1x,j1x),yy_sidex(i1x,j1x),'*')

hold off