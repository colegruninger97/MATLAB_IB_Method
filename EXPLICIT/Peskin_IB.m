function [u_new,v_new,p_new,F_new,X_new,X_tracer_new] = Peskin_IB(u,v,p,F,X,X_OG,dt,dx,dy,mu,kappa,ds,X_tracer)
%Explict implementation of the immersed boundary method following Peskin's
%lecture notes. 

r = 0.25;
rho = 1.0;
[ru,cu] = size(u);
[rv,cv] = size(v);
%interpolate the velocity onto the Eulerian grid
%  [U,V,i1x,j1x,i1y,j1y] = interpBS5BS4(u,v,X,dx,dy);
%  [Utracer,Vtracer,i1x,j1x,i1y,j1y] = interpBS5BS4(u,v,X_tracer,dx,dy);
 [U,V,i1x,j1x,i1y,j1y] = interpBS3BS2(u,v,X,dx,dy);
 [Utracer,Vtracer,i1x,j1x,i1y,j1y] = interpBS3BS2(u,v,X_tracer,dx,dy);
% [U,V,i1x,j1x,i1y,j1y] = interpBS1BS0(u,v,X,dx,dy);
% [U,V,i1x,j1x,i1y,j1y] = interpBS2BS1(u,v,X,dx,dy);
% [Utracer,Vtracer,i1x,j1x,i1y,j1y] = interpBS2BS1(u,v,X_tracer,dx,dy);
% [U,V,i1x,j1x,i1y,j1y] = interpIB4(u,v,X,dx,dy);
% [Utracer,Vtracer,i1x,j1x,i1y,j1y] = interpIB4(u,v,X_tracer,dx,dy);
%Update the positions of the Lagrangian structure at the half time-step
X_half(:,1) = X(:,1) + 0.5.*dt.*U(:);
X_half(:,2) = X(:,2) + 0.5.*dt.*V(:);
X_half(:,1) = mod(X_half(:,1),cu*dx);
X_half(:,2) = mod(X_half(:,2),rv*dy);
X_half_tracer(:,1) = X_tracer(:,1) + 0.5.*dt.*Utracer(:);
X_half_tracer(:,2) = X_tracer(:,2) + 0.5.*dt.*Vtracer(:);
X_half_tracer(:,1) = mod(X_half_tracer(:,1),cu*dx);
X_half_tracer(:,2) = mod(X_half_tracer(:,2),rv*dy);

%Update the forces on the Lagrangian grid
% F_half(:,1) = -kappa*(X_half(:,1) - X_OG(:,1));
% F_half(:,2) = -kappa*(X_half(:,2) - X_OG(:,2));
F_half = Elastic_force_update(X_half,kappa,ds);
% F_half = Elastic_Laplacian(X_half,kappa,ds);
% F_half = uniform_normal_force_circle(X_half,kappa,ds,r);
% F_half(:,1) = -eta*U;
% F_half(:,2) = -eta*V;
%Spread these forces onto the Eulerian grid
 [ffx_half,ffy_half] = spreadBS3BS2(F_half,X_half,u,v,dx,dy,ds);
% [ffx_half,ffy_half] = spreadBS2BS1(F_half,X_half,u,v,dx,dy,ds);
% [ffx_half,ffy_half] = spreadIB4(F_half,X_half,u,v,dx,dy,ds);

% Form the matrices needed to solve the linear system ....
%Laplacian for u
e = ones(cu,1); Lapu_x = spdiags([e -2*e e], -1:1,cu,cu)/(dx*dx);
e = ones(ru,1); Lapu_y = spdiags([e -2*e e], -1:1, ru, ru)/(dy*dy);
%Append entries to each matrix to account for the periodic boundary
%conditions
Lapu_x(1,end) = 1/(dx*dx);
Lapu_x(end,1) = 1/(dx*dx);
Lapu_y(1,end) = 1/(dy*dy);
Lapu_y(end,1) = 1/(dy*dy);

Ix = speye(cu); Iy = speye(ru);
Lapu = kron(Lapu_x,Iy) + kron(Ix,Lapu_y);
Lapu = mu.*Lapu;

%Create ghost cell entries for each term
[u_g,v_g] = Ghostnodesside_periodic(u,v);
%These terms are only used for the last solve
Lapu_rhs = (u_g(3:end,2:end-1) - 2.*u_g(2:end-1,2:end-1) + u_g(1:end-2,2:end-1))/(dy*dy) + ...
    (u_g(2:end-1,3:end) - 2.*u_g(2:end-1,2:end-1) + u_g(2:end-1,1:end-2))/(dx*dx);
Lapv_rhs = (v_g(3:end,2:end-1) - 2.*v_g(2:end-1,2:end-1) + v_g(1:end-2,2:end-1))/(dy*dy) + ...
    (v_g(2:end-1,3:end) - 2.*v_g(2:end-1,2:end-1) + v_g(2:end-1,1:end-2))/(dx*dx);

Lapu_rhs = (0.5)*(mu).*Lapu_rhs;
Lapv_rhs = (0.5)*(mu).*Lapv_rhs;



%Build the Laplacian for v

e = ones(cv,1); Lapv_x = spdiags([e -2*e e], -1:1,cv,cv)/(dx*dx);
%Edit the matrix to impose dirichlet boundart conditions
Lapv_x(1,end) = 1/(dx*dx);
Lapv_x(end,1) = 1/(dx*dx);
e = ones(rv,1); Lapv_y = spdiags([e -2*e e], -1:1, rv, rv)/(dy*dy);
Lapv_y(1,end) = 1/(dy*dy);
Lapv_y(end,1) = 1/(dy*dy);
Ix = speye(cv); Iy = speye(rv);
Lapv = kron(Lapv_x,Iy) + kron(Ix,Lapv_y);
Lapv = mu.*Lapv;


%Build Differentiaion matrix for the pressure in the X and Y directions
[rp,cp] = size(p);
G_x = zeros(cp,cp);
G_x(1:cp+1:end) = 1/dx;
G_x(2:cp+1:end) = -1/dx;
G_x(1,end) = -1/dx;
Ix = speye(rp,rp);
G_x = kron(G_x,Ix);


G_y = zeros(rp,rp);
G_y(1:rp+1:end) = 1/dy;
G_y(2:rp+1:end) = -1/dy;
G_y(1,end) = -1/dy;
Iy = speye(cp,cp);
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
A1 = 2.*(rho/dt).*speye(size(Lapu));
A2 = 2.*(rho/dt).*speye(size(Lapv));
N1 = zeros(size(A1));
N2 = zeros(size(A2));
N3 = zeros(rp*cp);
N7 = zeros(rp*cp,1);
%Setup up the explicit convective terms
[u_g,v_g] = Ghostnodesside_periodic(u,v);

uce = 0.5.*(u_g(:,1:end-1) + u_g(:,2:end));
vce = 0.5.*(v_g(1:end-1,:) + v_g(2:end,:));

u_corn = 0.25.*(uce(1:end-1,1:end-1) + uce(2:end,1:end-1) + uce(1:end-1,2:end)+ uce(2:end,2:end));
v_corn = 0.25.*(vce(1:end-1,1:end-1) + vce(2:end,1:end-1) + vce(1:end-1,2:end)+ vce(2:end,2:end));
%interpolate v onto side centered components
v_sidex = 0.5.*(vce(2:end,1:end-2)+vce(2:end,2:end-1));
u_sidey = 0.5.*(uce(1:end-2,2:end) + uce(2:end-1,2:end));
Conv_u = u.*(uce(2:end-1,2:end)-uce(2:end-1,1:end-1))./dx + v_sidex.*(u_corn(2:end,:)-u_corn(1:end-1,:))./dy;
Conv_v = u_sidey.*(v_corn(:,2:end) - v_corn(:,1:end-1))./dx + v.*(vce(2:end,2:end-1) - vce(1:end-1,2:end-1))./dy;
Conv_u = rho.*Conv_u;
Conv_v = rho.*Conv_v;

LHS = [A1 - Lapu, N1, G_x; N2, A2 - Lapv, G_y; -Dx, -Dy,N3];
RHS = [-Conv_u(:) + A1*u(:) + ffx_half(:) ;...
    -Conv_v(:) +  A2*v(:) + ffy_half(:); ...
   N7(:)];
 
vals = LHS\RHS;
u_half = vals(1:ru*(cu));
u_half = reshape(u_half,ru,cu);
v_half = vals(ru*(cu)+1:ru*(cu)+(rv)*cv);
v_half = reshape(v_half,rv,cv);





% Now again update the the lagrange points 
[U,V,i1x,j1x,i1y,j1y] = interpBS3BS2(u_half,v_half,X_half,dx,dy);
[Utracer,Vtracer,i1x,j1x,i1y,j1y] = interpBS3BS2(u_half,v_half,X_half_tracer,dx,dy);
% [U,V,i1x,j1x,i1y,j1y] = interpBS5BS4(u_half,v_half,X_half,dx,dy);
% [Utracer,Vtracer,i1x,j1x,i1y,j1y] = interpBS5BS4(u_half,v_half,X_half_tracer,dx,dy);
% [U,V,i1x,j1x,i1y,j1y] = interpBS1BS0(u_half,v_half,X_half,dx,dy);
% [U,V,i1x,j1x,i1y,j1y] = interpBS2BS1(u_half,v_half,X_half,dx,dy);
% [Utracer,Vtracer,i1x,j1x,i1y,j1y] = interpBS2BS1(u_half,v_half,X_half_tracer,dx,dy);
% [U,V,i1x,j1x,i1y,j1y] = interpIB4(u_half,v_half,X_half,dx,dy);
% [Utracer,Vtracer,i1x,j1x,i1y,j1y] = interpIB4(u_half,v_half,X_half_tracer,dx,dy);
X_new(:,1) = X(:,1) + dt*U(:);
X_new(:,2) = X(:,2) + dt*V(:);
X_new(:,1) = mod(X_new(:,1),cu*dx);
X_new(:,2) = mod(X_new(:,2),rv*dy);
X_tracer_new(:,1) = X_tracer(:,1) + dt*Utracer(:);
X_tracer_new(:,2) = X_tracer(:,2) + dt*Vtracer(:);
X_tracer_new(:,1) = mod(X_tracer_new(:,1),cu*dx);
X_tracer_new(:,2) = mod(X_tracer_new(:,2),rv*dy);
%Update the forces on the Lagrangian grid
% F_new(:,1) = -kappa*(X_new(:,1) - X_OG(:,1));
% F_new(:,2) = -kappa*(X_new(:,2) - X_OG(:,2));
F_new = Elastic_force_update(X_new,kappa,ds);
% F_new = Elastic_Laplacian(X_new,kappa,ds);
% F_new = uniform_normal_force_circle(X_new,kappa,ds,r);
% F_new(:,1) = -eta*U;
% F_new(:,2) = -eta*V;

%Spread these forces to the Eulerian grid
[ffx_new,ffy_new] = spreadBS3BS2(F_new,X_new,u_half,v_half,dx,dy,ds);
% [ffx_new,ffy_new] = spreadBS2BS1(F_new,X_new,u_half,v_half,dx,dy,ds);
% [ffx_new,ffy_new] = spreadIB4(F_new,X_new,u_half,v_half,dx,dy,ds);
 %Do another solve for the new velocity and pressure
 
[u_halfg,v_halfg] = Ghostnodesside_periodic(u_half,v_half);


uce_half = 0.5.*(u_halfg(:,1:end-1) + u_halfg(:,2:end));
vce_half = 0.5.*(v_halfg(1:end-1,:) + v_halfg(2:end,:));

u_corn_half = 0.25.*(uce_half(1:end-1,1:end-1) + uce_half(2:end,1:end-1) + uce_half(1:end-1,2:end)+ uce_half(2:end,2:end));
v_corn_half = 0.25.*(vce_half(1:end-1,1:end-1) + vce_half(2:end,1:end-1) + vce_half(1:end-1,2:end)+ vce_half(2:end,2:end));
%interpolate v onto side centered components
v_sidex = 0.5.*(vce_half(2:end,1:end-2) + vce_half(2:end,2:end-1));
u_sidey = 0.5.*(uce_half(1:end-2,2:end) + uce_half(2:end-1,2:end));
Conv_u_half = u_half.*(uce_half(2:end-1,2:end)-uce_half(2:end-1,1:end-1))./dx + v_sidex.*(u_corn_half(2:end,:)-u_corn_half(1:end-1,:))./dy;
Conv_v_half = u_sidey.*(v_corn_half(:,2:end) - v_corn_half(:,1:end-1))./dx + v_half.*(vce_half(2:end,2:end-1) - vce_half(1:end-1,2:end-1))./dy;
Conv_u_half = rho.*Conv_u_half;
Conv_v_half = rho.*Conv_v_half;


A1 = 0.5.*A1;
A2 = 0.5.*A2;
LHS = [A1 - 0.5.*Lapu, N1, G_x; N2, A2 - 0.5.*Lapv, G_y; -Dx, -Dy,N3];
RHS = [-Conv_u_half(:) + A1*u(:) + Lapu_rhs(:) + ffx_new(:);...
    -Conv_v_half(:) +  A2*v(:) + Lapv_rhs(:) + ffy_new(:);...
     N7(:)];
vals = LHS\RHS;
u_new = vals(1:ru*(cu));
u_new = reshape(u_new,ru,cu);
v_new = vals(ru*(cu)+1:ru*(cu)+(rv)*cv);
v_new = reshape(v_new,rv,cv);
p_new = vals(ru*(cu) + (rv)*cv + 1:end);
p_new = p_new - dx*dy*sum(p_new,'all');
%Normalize the pressure: the gradient is really all that matters
p_new = reshape(p_new,rp,cp);


end
 
 






