% Main driver for an explict Immersed boundary method. In this case, the
% dirac delta functions are not fixed and the singular forces are thus
% nonlinear. I'm trying to understand how the relative spacing of the
% Lagrangian mesh affects convergence in different flow regimes...
addpath('../Ghost_cells/','../Kernels/');
clear;
clc;
%Define the initial geometry ...
%The following grid setup is particular to the periodic case...
N = 128;
Lx = 1.0; %Make these variables global later...
Ly = 1.0;
dy = Ly/N; %Mesh spacing is uniform for simplicity
dx = dy;
mu = 0.1;
rho = 1.0;
dt = Lx/(20*N);
eta = 3*rho*dx/dt;
x = 0:dx:Lx-dx;
y = 0:dy:Ly-dy;
x_c = dx/2:dx:Lx-dx/2;
y_c = dy/2:dy:Ly-dy/2;

Px = -0.5;

[xx_cent,yy_cent] = meshgrid(x_c,y_c);
[xx_sidex,yy_sidex] = meshgrid(x,y_c);
[xx_sidey,yy_sidey] = meshgrid(x_c,y);
%Setup the intial condtions ---> Taylor vortices
u_exact = @(x,y,t) 0.*(1+2*exp(-8.*pi.*pi.*(mu/rho).*t).*sin(2*pi*(y-t)).*cos(2*pi*(x-t)));
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
%

mfac = 0.5;
r_cyl = 0.25;
%Specify the Lagrangian mesh spacing
ds = mfac*dx/r_cyl;
approx = round(2*pi/ds);
ds = 2*pi/approx;
dX = mfac*dx;
l = length(0:ds:2*pi-ds);
X_1 = Lx.*(0.5.*ones(1,l) + (0.25)*cos(0:ds:2*pi-ds));
X_2 = Lx.*(0.5.*ones(1,l) +  (0.25)*sin(0:ds:2*pi-ds));
X = [X_1;X_2]';
% X_1 = Ly/4:dX:3*Ly/4;
% 
% X_2 = Lx/2.*ones(1,length(X_1));
% X_1 = 0.5;
% X_2 = 0.5;
% X_OG = [X_1;X_2];
% X_OG = X_OG';
% X = X_OG;
NIB = length(X_2);
N_IB_tracer = 20*NIB;
ds_tracer = 2*pi/N_IB_tracer;
X_1_tracer = Lx.*(0.5.*ones(1,N_IB_tracer) + (0.25).*cos(0:ds_tracer:2*pi-ds_tracer));
X_2_tracer = Lx.*(0.5.*ones(1,N_IB_tracer) + (0.25).*sin(0:ds_tracer:2*pi-ds_tracer));
X_tracer = [X_1_tracer;X_2_tracer];
X_tracer = X_tracer';
% X_tracer_OG = X_tracer;
kappa = 1;
T_final = 0.25;
F= zeros(length(X_2),2);
uce = zeros(length(y_c),length(x_c));
vce = uce;
t = 0;
Init_Poly_area = polyarea(X(:,1),X(:,2));

%Take an initial timestep using Peskin's IB method
[u,v,p,F,X,X_tracer] = Krylov_Peskin_IB(u,v,p,F,X,X_tracer,dt,dx,dy,mu,kappa,ds);
t = t+dt;

% setup spline interpolant of tracers after initial RK2 step.
pp1 = spline([0:ds_tracer:2*pi],[X_tracer(:,1);X_tracer(1,1)]);
pp2 = spline([0:ds_tracer:2*pi],[X_tracer(:,2);X_tracer(1,2)]);
pp2_deriv = pp_deriv(pp2);
pp_mult = pp_multiply(pp1,pp2_deriv);
Spline_Area = intpp(pp_mult);
Spline_Area_error_rel(1) = abs(Spline_Area - pi*0.25*0.25)./(0.25*0.25*pi);

Rel_Poly_Area_error(1) = abs(polyarea(X(:,1),X(:,2)) - Init_Poly_area)/Init_Poly_area;


%I might want to instead try solving for the Schur complement
k = 1;
while t < T_final
%Now use the Adams Bashforth scheme the rest of the way
if t + dt > T_final
    dt = T_final - t;
end
[u_new,v_new,p,F,X,X_tracer] = Krylov_IB_explicit_AB2(u,v,p,u0,v0,F,X,X_tracer,dt,dx,dy,mu,kappa,ds);
%Record the solution to use for the next time-step
u0 = u;
v0 = v;
u = u_new;
v = v_new;

pp1 = spline([0:ds_tracer:2*pi],[X_tracer(:,1);X_tracer(1,1)]);
pp2 = spline([0:ds_tracer:2*pi],[X_tracer(:,2);X_tracer(1,2)]);
pp2_deriv = pp_deriv(pp2);
spline_mult = pp_multiply(pp1,pp2_deriv);
Spline_Area = intpp(spline_mult);
Spline_Area_error_abs(k) = abs(Spline_Area - r_cyl*r_cyl*pi);
Spline_Area_error_rel(k) = Spline_Area_error_abs(k)./r_cyl*r_cyl*pi;

Rel_Poly_Area_error(k) = abs(polyarea(X(:,1),X(:,2)) - Init_Poly_area)/Init_Poly_area;


%  du = (u(:,2:end)-u(:,1:end-1))./dx;
%  dv = (v(2:end,:)-v(1:end-1,:))./dy;
%  div = du(1:end-1,:) + dv(:,1:end-1)
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
% plot(X_tracer(:,1),X_tracer(:,2))
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
% surf(xx_cent,yy_cent,sqrt(uce.^2 + vce.^2))
colormap('default')
hold on
%quiver(xx_cent(1:1:end),yy_cent(1:1:end),uce(1:1:end),vce(1:1:end),'w','linewidth',1.2);
axis([0+dx  Lx-dx  0+dy  Ly-dy])
% plot(X_tracer(:,1),X_tracer(:,2),'LineWidth',2)
plot(X(:,1),X(:,2),'o','linewidth',2);
quiver(X(:,1),X(:,2),F(:,1),F(:,2),'k')
hold off
pause(0)
% Area_new_polygon(k) = polyarea(X_tracer(:,1),X_tracer(:,2));
% Area_spectral(k) = Spectral_Area(X_tracer,ds_tracer);
% Area_error_spectral(k) = abs(Area_spectral(k) - 0.25*0.25*pi)./(0.25*0.25*pi);
% area_error_polygon(k) = abs(Area_new_polygon(k) - pi*0.25*0.25)/(0.25*0.25*pi);
%  surf(p)
%  pause(0)
 %surf(u)
%Compute the errors
%  erru = abs(u(:,:)-u_exact(xx_sidex,yy_sidex,t));
%  errv = abs(v(:,:)-v_exact(xx_sidey,yy_sidey,t));
%  L1_error = sum(dx*dy*erru(:,:),'all');
%  L2_error = sqrt(dx*dy*sum(erru(:,:).^2,'all'));
%  L_inf_error = max(errv,[],'all')
%  surf(erru+errv);
%   pause(0)
% % max(erru + errv,[],'all')
k = k+1;
t = t + dt
end
save('Spline_Rel_Error_BS3BS2.mat','Spline_Area_error_rel');
save('Rel_Poly_Area_Error_BS3BS2.mat','Rel_Poly_Area_error');
%Take spectral derivative
% clear;
% clc;
% N=24;
% h = 2*pi/N;
% theta = 0:h:2*pi-h;
% y = 1.*sin(theta);
% y_hat = fft(y);
% y_hat_prime = 1i.*[0:N/2-1 0 -N/2+1:-1] .* y_hat;
% y_prime = real(ifft(y_hat_prime));
% y_prime_true = cos(theta);
% error = abs(y_prime-y_prime_true);
% plot(error)
% ds = 2*pi/100;
% X = 2.*cos(0:ds:2*pi-ds);
% Y = 2.*sin(0:ds:2*pi-ds);
% 
% Area = Spectral_Area([X;Y]',ds);
% True_Area = 4*pi;
% error = abs(Area - True_Area)


% clear;
% clc;
% N = 128;
% ds = 2*pi/(N);
% x = cos(0:ds:2*pi);
% y = sin(0:ds:2*pi);
% s = 0:ds:2*pi;
% pp1 = spline(s,x);
% pp2 = spline(s,y);
% 
% 
% 
% pp2_deriv = pp_deriv(pp2);
% Y_s = ppval(pp2_deriv,s);
% X = ppval(pp1,s);
% error_trap = abs(sum(X(1:end-1).*Y_s(1:end-1)*ds,'all')-pi)
% pp_times_pp = pp_multiply(pp1,pp2_deriv);
% yy = ppval(pp_times_pp,[0:0.01:2*pi]);
% q = polyint(pp1.coefs(2,:)); 
% int = intpp(pp_times_pp)
% error_exact = abs(int-pi)
% polygonal_error = abs(polyarea(x(1:end-1),y(1:end-1)) - pi)

function defint = intpp(pp)
I = 0;
for k = 1:length(pp.breaks)-1
I = I + integral(@(x)ppval(pp,x),pp.breaks(k),pp.breaks(k+1));
end
defint = I;
end


function deriv = pp_deriv(pp)


for k = 1:pp.order
    pp.coefs(:,k) = pp.coefs(:,k).*(pp.order - k);
end


deriv = mkpp(pp.breaks,pp.coefs(:,1:pp.order-1));
end



function pp_times_pp = pp_multiply(pp1,pp2);
breaks = pp1.breaks;
order = pp1.order + pp2.order-1;
coefs = zeros(length(breaks)-1,order);
for i = 1:length(breaks)-1
    coefs(i,:) = conv(pp1.coefs(i,:),pp2.coefs(i,:));
end
    pp_times_pp = mkpp(breaks,coefs);

end


function Area = Spectral_Area(X_tracer,ds_tracer)
%Use Green's theorem, spectral differentiation, and perodic trapezoidal
%rule to obtain a spectrally accurate approximation of the area of the
%structure advected by the flow map
N_tracer = length(X_tracer(:,1));
%First approximate the tangential derivative in Y
Y_hat = fft(X_tracer(:,2));
Y_hat_prime = 1i.*[0:N_tracer/2-1 0 -N_tracer/2+1:-1]' .* Y_hat;
Y_prime = real(ifft(Y_hat_prime));

Area = sum(X_tracer(:,1).*Y_prime(:).*ds_tracer,'all');



end




