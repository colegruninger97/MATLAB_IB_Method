%Main driver for the IB method solver on periodic square
%In this implementation, we ignore the convective terms of the Navier Stokes equations in order 
% %to linearize the problem making it tractable for analysis
% mfac = 10.0;
% z = pi-0.000001:0.0000000001:pi+0.000001;
% kk = length(z);
% z = complex(zeros(1,kk),z);
% f = ones(1,kk) + exp(z);
% plot(real(f),imag(f));
% hold on
% for i = 1:9

% mfac = mfac - 1.0;
%Intialize the Lagrangian structure to be immersed in the fluid

%For now we set the lagrangian structure to just be a flat plate facing
%who's normal derivative faces towards a parbolic like flow field
N = 32;
Lx = 1.0; %Make these variables global later...
Ly = 1.0;
dx = Lx/N; %Mesh spacing is uniform for simplicity
dy = Ly/N;
mu = 0.075;
rho = 1.0;
dt = 0.005;
x = 0:dx:Lx;
y = 0:dy:Ly;
Px = -0.5;

x_c = 0.5*(x(1:end-1) + x(2:end));
y_c = 0.5*(y(1:end-1) + y(2:end));
[xx_cent,yy_cent] = meshgrid(x_c,y_c);
[xx_sidex,yy_sidex] = meshgrid(x,y_c);
[xx_sidey,yy_sidey] = meshgrid(x_c,y);
xx_centd = xx_cent(1:4:end,1:4:end);
yy_centd = yy_cent(1:4:end,1:4:end);


mfac = 4.0;
%Specify the Lagrangian mesh spacing
dX = mfac*dx;
l = length(0:dX:2*pi-dX);
% X_1 = 2.5.*ones(1,l) + 0.2*cos(0:dX:2*pi-dX);
% X_2 = 0.5.*ones(1,l) +  0.2*sin(0:dX:2*pi-dX);
X_2 = Ly/4:dX:3*Ly/4;

X_1 = (Ly/4).*ones(1,length(X_2));

X = [X_1;X_2];
X = X';
NIB = length(X_2);
kappa = 1e15;



% u_exact = @(x,y,t) 1-2*exp(-8*pi*pi*(mu)/rho*t).*sin(2*pi*(y-t)).*cos(2*pi*(x-t));
% v_exact = @(x,y,t) 1+2*exp(-8*pi*pi*(mu)/rho*t).*cos(2*pi*(y-t)).*sin(2*pi*(x-t));
% p_exact = @(x,y,t) -exp(-16*pi*pi*(mu)/rho*t)*(cos(4*pi*(x-t))+cos(4*pi*(y-t)));

u_exact = @(x,y,t) 3 + 0.*x.*y;
v_exact = @(x,y,t) -0*x.*y;
p_exact = @(x,y,t) 0.*x.*y;

F1 = zeros(length(X_2),1);
F2 = zeros(length(X_2),1);

u = u_exact(xx_sidex,yy_sidex,0);
v = v_exact(xx_sidey,yy_sidey,0);
p = 0.*p_exact(xx_cent,yy_cent,0);


t = 0;


%Now evovle the system 

while(t < 30*dt)
    %tic
[u,v,p,F1,F2,LHS,D] = BE_IB_Evolve_Periodic(u,v,p,F1,F2,dt,dx,dy,N,NIB,mu,rho,kappa,mfac,x,y,X);
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
plot(X_1,X_2,'o')
hold off
pause(0)
% surf(x,y_c,u)
% pause(0)
end
% 
% plot(real(diag(D)),imag(diag(D)),'o')
% axis([-4e-10 4e-10 -2.5e-20 2.6e-20])
% rank(V)
% toc
% pause(0)
% 
% end
%Need to check the lift and drift coefficients to check accuracy....


