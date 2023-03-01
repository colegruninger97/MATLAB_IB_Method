%Main driver for the IB method solver
%In this implementation, we ignore the convective terms of the Navier Stokes equations in order 
%to linearize the problem making it tractable for analysis




%Intialize the Lagrangian structure to be immersed in the fluid
clear;
clc;
%For now we set the lagrangian structure to just be a flat plate facing
%who's normal derivative faces towards a parbolic like flow field
N = 214;
Lx = 1; %Make these variables global later...
Ly = 1;
dx = Lx/N; %Mesh spacing is uniform for simplicity
dy = Ly/N;
mu = 1.0;
rho = 1.0;
dt = 1.0;
x = 0:dx:Lx;
y = 0:dy:Ly;
Px = -0.5;

x_c = 0.5*(x(1:end-1) + x(2:end));
y_c = 0.5*(y(1:end-1) + y(2:end));
[xx_cent,yy_cent] = meshgrid(x_c,y_c);
[xx_sidex,yy_sidex] = meshgrid(x,y_c);
[xx_sidey,yy_sidey] = meshgrid(x_c,y);

mfac = 3.0;
r_cyl = 0.1;
%Specify the Lagrangian mesh spacing
%for cylinder
dX = mfac*dx/r_cyl;
jeez = ceil(2*pi/dX);
dX = 2*pi/jeez;
l = length(0:dX:2*pi-dX);
X_1 = 0.5.*ones(1,l) + r_cyl*cos(0:dX:2*pi-dX);
X_2 = 0.5.*ones(1,l) +  r_cyl*sin(0:dX:2*pi-dX);
% X_2 = Ly/4:dX:3*Ly/4;
% 
% X_1 = Lx/2.*ones(1,length(X_2));

X = [X_1;X_2];
X = X';
NIB = length(X_2);
kappa = 1;



%Setup the linear system to solve
u_exact = @(x,y,t) 0.*x + (Px/(2*mu)).*y.*(y-Ly);
v_exact = @(x,y,t) 0.*x.*y;
p_exact = @(x,y,t) 0.*y + Px.*x - Px/4;
F1 = zeros(length(X_2),1);
F2 = zeros(length(X_2),1);

u = u_exact(xx_sidex,yy_sidex,0);
u = 0.*u(:,2:end-1);
u = [u_exact(x(1),y_c,0)',u,u_exact(x(end),y_c,0)'];
v = v_exact(xx_sidey,yy_sidey,dt);
v = 0.*v(2:end-1,:);
v = [v_exact(x_c,y(1),dt);v;v_exact(x_c,y(end),0)];
p = p_exact(xx_cent,yy_cent,0);
p = 0.*p;
xx_centd = xx_cent(1:4:end,1:4:end);
yy_centd = yy_cent(1:4:end,1:4:end);

uce = 0.5.*(u(:,2:end) + u(:,1:end-1));
vce = 0.5.*(v(2:end,:) + v(1:end-1,:));
uced = uce(1:4:end,1:4:end);
vced = vce(1:4:end,1:4:end);
contourf(xx_cent,yy_cent,sqrt(uce.^2 + vce.^2))
colormap('default')
hold on
quiver(xx_centd,yy_centd,uced,vced,'w','linewidth',1.5);
axis([0+dx  Lx-dx  0+dy    Ly-dy])
plot(X_1,X_2,'ro')
hold off
pause(0)

t = 0;


%Now evovle the system 

while(t < 1*dt)
[u,v,p,F1,F2,LHS,A,B1,B2] = IMPLICIT_PRESCRIBED_MOTION(u,v,p,F1,F2,dt,dx,dy,N,NIB,mu,rho,mfac,x,y,X);
%Plot the vector field 
uce = 0.5.*(u(:,2:end) + u(:,1:end-1));
vce = 0.5.*(v(2:end,:) + v(1:end-1,:));
uced = uce(1:4:end,1:4:end);
vced = vce(1:4:end,1:4:end);
contourf(xx_cent,yy_cent,sqrt(uce.^2 + vce.^2))
colormap('default')
hold on
quiver(xx_centd,yy_centd,uced,vced,'w','linewidth',1.5);
axis([0+dx  Lx-dx  0+dy    Ly-dy])
plot(X_1,X_2,'ro')
%Plot the vector field F (of lagrange multipliers?) at each of the IB
%points
quiver(X_1,X_2,F1',F2',2,'y','linewidth',1.5)
pause(0)
% surf(x,y_c,u)
% pause(0)
t = t + dt;
end

% Compute the Schur complement ...
% Schur = -B2*inv(A)*B1;


