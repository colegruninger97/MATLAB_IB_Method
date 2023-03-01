% Main driver for an explict Immersed boundary method. In this case, the
% dirac delta functions are not fixed and the singular forces are thus
% nonlinear. I'm trying to understand how the relative spacing of the
% Lagrangian mesh affects convergence in different flow regimes...

clear;
clc;
%Define the initial geometry ...
%The following grid setup is particular to the periodic case...
N = 16;
Lx = 1.0; %Make these variables global later...
Ly = 1.0;
dy = Ly/N; %Mesh spacing is uniform for simplicity
dx = dy;
mu = 0.1;
rho = 1.0;
dt = Lx/(4*N);
eta = 3*rho*dx/dt
x = 0:dx:Lx-dx;
y = 0:dy:Ly-dy;
x_c = dx/2:dx:Lx-dx/2;
y_c = dy/2:dy:Ly-dy/2;

Px = -0.5;

[xx_cent,yy_cent] = meshgrid(x_c,y_c);
[xx_sidex,yy_sidex] = meshgrid(x,y_c);
[xx_sidey,yy_sidey] = meshgrid(x_c,y);
%Use Taylor vortices as the analytic solution
u_exact = @(x,y,t) 1+0.*2*exp(-8.*pi.*pi.*(mu/rho).*t).*sin(2*pi*(y-t)).*cos(2*pi*(x-t));
v_exact = @(x,y,t) 0.*(1-2*exp(-8.*pi.*pi.*((mu)./rho).*t).*cos(2.*pi.*(y-t)).*sin(2.*pi.*(x-t)));
p_exact = @(x,y,t) 0.*(-exp(-16*pi*pi*((mu)/rho).*t).*(cos(4*pi*(x-t))+cos(4*pi.*(y-t))));

u = u_exact(xx_sidex,yy_sidex,0);
v = v_exact(xx_sidey,yy_sidey,0);
p = p_exact(xx_cent,yy_cent,0); 
u0 = u;
v0 = v;
p0 = p;
% u = ones(size(xx_sidex));
% % u(:,2:end) = 0;
% v = zeros(size(yy_sidey));
% p = zeros(size(xx_cent));
% u = sin(2*pi.*xx_sidex).*cos(2*pi*yy_sidex);
% v = sin(2*pi.*yy_sidey).*cos(2*pi.*xx_sidey);


mfac = 1.0;
r_cyl = 0.25;
%Specify the Lagrangian mesh spacing
ds = mfac*dx/r_cyl;
jeez = ceil(2*pi/ds);
ds = 2*pi/jeez;
dX = mfac*dx;
l = length(0:ds:2*pi-ds);
X_1 = 0.5.*ones(1,l) + (1/4)*cos(0:ds:2*pi-ds);
X_2 = 0.5.*ones(1,l) +  (1/4)*sin(0:ds:2*pi-ds);
% X_1 = Ly/4:dX:3*Ly/4;
% 
% X_2 = Lx/2.*ones(1,length(X_1));
% X_1 = 0.5;
% X_2 = 0.5;
X_OG = [X_1;X_2];
X_OG = X_OG';
X = X_OG;
NIB = length(X_2);
kappa = 1;



F(:,1) = zeros(length(X_2),1);
F(:,2) = zeros(length(X_2),1);
uce = zeros(length(y_c),length(x_c));
vce = uce;
t = 0;
%Take an initial timestep using Peskin's IB method
[u,v,p,F,X] = Peskin_IB(u,v,p,F,X,X_OG,dt,dx,dy,mu,eta,dX);
t = t+dt;

%Build LHS so I don't do it everytimestep
[ru,cu] = size(u);
[rv,cv] = size(v);
% Form the matrices needed to solve the linear system ....
%Laplacian for u
e = ones(cu,1); Lapu_x = spdiags([e -2*e e], -1:1, cu, cu)/(dx*dx);
e = ones(ru,1); Lapu_y = spdiags([e -2*e e], -1:1, ru, ru)/(dy*dy);
%Append entries to each matrix to account for the periodic boundary
%conditions
Lapu_x(1,end) = 1/(dx*dx);
Lapu_x(end,1) = 1/(dx*dx);
Lapu_y(1,end) = 1/(dy*dy);
Lapu_y(end,1) = 1/(dy*dy);

Ix = speye(cu); Iy = speye(ru);
Lapu = kron(Lapu_x,Iy) + kron(Ix,Lapu_y);
Lapu = 0.5.*mu.*Lapu;

e = ones(cv,1); Lapv_x = spdiags([e -2*e e], -1:1,cv,cv)/(dx*dx);
%Edit the matrix to impose dirichlet boundart conditions
Lapv_x(1,end) = 1/(dx*dx);
Lapv_x(end,1) = 1/(dx*dx);
e = ones(rv,1); Lapv_y = spdiags([e -2*e e], -1:1, rv, rv)/(dy*dy);
Lapv_y(1,end) = 1/(dy*dy);
Lapv_y(end,1) = 1/(dy*dy);
Ix = speye(cv); Iy = speye(rv);
Lapv = kron(Lapv_x,Iy) + kron(Ix,Lapv_y);
Lapv = 0.5.*mu.*Lapv;

[rp,cp] = size(p);
G_x = zeros(cp,cp);
G_x(1:cp+1:end) = 1/dx;
G_x(2:cp+1:end) = -1/dx;
G_x(1,end) = -1/dx;
Ix = eye(rp,rp);
G_x = kron(G_x,Ix);


G_y = zeros(rp,rp);
G_y(1:rp+1:end) = 1/dy;
G_y(2:rp+1:end) = -1/dy;
G_y(1,end) = -1/dy;
Iy = eye(cp,cp);
G_y = kron(Iy,G_y);



%Build Differentiation matrices for the continuity equation
%Both are approximated via centered differences w/r/t the centerd nodes
[ru,cu] = size(u);
Dx = zeros(cu,cu);
Dx(1:cu+1:end) = -1/dx;
Dx(cu+1:cu+1:end) = 1/dx;
Dx(end,1) = 1/dx;
Ix = speye(ru,ru);
Dx = kron(Dx,Ix);


[rv,cv] = size(v);
Dy = zeros(rv,rv);
Dy(1:rv+1:end) = -1/dy;
Dy(rv+1:rv+1:end) = 1/dy;
Dy(end,1) = 1/dy;
Iy = speye(cv,cv);
Dy = kron(Iy,Dy);

%Start assemblying the linear system to be used
A1 = (rho/dt).*speye(size(Lapu));
A2 = (rho/dt).*speye(size(Lapv));
N1 = zeros(size(A1));
N2 = zeros(size(A2));
N3 = zeros(rp*cp);

LHS = [A1 - Lapu, N1, G_x; N2, A2 - Lapv, G_y; -Dx, -Dy,N3];
%I might want to instead try solving for the Schur complement

k = 1;
while t < 3000*dt
%Now use the Adams Bashforth scheme the rest of the way
[u_new,v_new,p,F,X] = IB_explicit_AB2(u,v,p,u0,v0,F,X,X_OG,dt,dx,dy,mu,eta,dX,LHS,A1,A2);
t = t+dt;
%Record the solution to use for the next time-step
u0 = u;
v0 = v;
u = u_new;
v = v_new;

% du = (u(:,2:end)-u(:,1:end-1))./dx;
% dv = (v(2:end,:)-v(1:end-1,:))./dy;
% div = du(1:end-1,:) + dv(:,1:end-1);
% surf(div)
% pause(0)
%Tracking the kernel testing
% uce(:,1:end-1) = 0.5.*(u(:,2:end) + u(:,1:end-1));
% uce(:,end) = 0.5.*(u(:,1) + u(:,end));
% vce(1:end-1,:) = 0.5.*(v(2:end,:) + v(1:end-1,:));
% vce(end,:) = 0.5.*(v(1,:) + v(end,:));
% [U,V,i1x,j1x,i1y,j1y] = interpBS3BS2(u,v,X,dx,dy);
% % contourf(xx_cent,yy_cent,sqrt(uce.^2 + vce.^2))
% plot(X(:,1),X(:,2),'*','LineWidth',2)
% hold on
% % plot(xx_sidex,yy_sidex,'ko')
% plot(xx_sidey,yy_sidey,'ko')
% % plot(xx_sidex(i1x,j1x),yy_sidex(i1x,j1x),'o','LineWidth',2)
% plot(xx_sidey(i1y,j1y),yy_sidey(i1y,j1y),'o','LineWidth',2)
% hold off
% axis([0 1 0 1])
% pause(1.0)
uce(:,1:end-1) = 0.5.*(u(:,2:end) + u(:,1:end-1));
uce(:,end) = 0.5.*(u(:,1) + u(:,end));
vce(1:end-1,:) = 0.5.*(v(2:end,:) + v(1:end-1,:));
vce(end,:) = 0.5.*(v(1,:) + v(end,:));
contourf(xx_cent,yy_cent,sqrt(uce.^2 + vce.^2))
colormap('default')
hold on
quiver(xx_cent(1:1:end),yy_cent(1:1:end),uce(1:1:end),vce(1:1:end),'w','linewidth',1.2);
axis([0+dx  Lx-dx  0+dy  Ly-dy])
plot(X(:,1),X(:,2),'o','LineWidth',2)
%quiver(X(:,1),X(:,2),F(:,1),F(:,2),'k')
hold off
pause(0)
% Area_new(k) = polyarea(X(:,1),X(:,2));
% area_error(k) = abs(Area_new(k) - polyarea(X_OG(:,1),X_OG(:,2)))/(polyarea(X_OG(:,1),X_OG(:,2)))
% surf(x_c,y_c,p)
% surf(u)
%Compute the errors
% erru = abs(u(:,:)-u_exact(xx_sidex,yy_sidex,t));
% errv = abs(v(:,:)-v_exact(xx_sidey,yy_sidey,t));
% surf(erru+errv);
% pause(0)
% max(erru + errv,[],'all')
k = k+1;
end

