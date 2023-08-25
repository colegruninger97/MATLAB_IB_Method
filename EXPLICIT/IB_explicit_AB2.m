function [u_new,v_new,p_new,F_new,X_new,X_tracer_new] = IB_explicit_AB2(u,v,p,u0,v0,F,X,X_OG,dt,dx,dy,mu,kappa,ds,LHS,A1,A2,X_tracer)
%Explict implementation of the immersed boundary method following Peskin's
%lecture notes. 

persistent N7
rho = 1.0;
r = 0.25;
[ru,cu] = size(u);
[rv,cv] = size(v);
%interpolate the velocity onto the Eulerian grid
[U,V,i1x,j1x,i1y,j1y] = interpBS3BS2(u,v,X,dx,dy);
[Utracer,Vtracer,i1x,j1x,i1y,j1y] = interpBS3BS2(u,v,X_tracer,dx,dy);
% [U,V,i1x,j1x,i1y,j1y] = interpBS5BS4(u,v,X,dx,dy);
% [Utracer,Vtracer,i1x,j1x,i1y,j1y] = interpBS5BS4(u,v,X_tracer,dx,dy);
% [U,V,i1x,j1x,i1y,j1y] = interpBS1BS0(u,v,X,dx,dy);
% [U,V,i1x,j1x,i1y,j1y] = interpBS2BS1(u,v,X,dx,dy);
% [Utracer,Vtracer,i1x,j1x,i1y,j1y] = interpBS2BS1(u,v,X_tracer,dx,dy);
% [U,V,i1x,j1x,i1y,j1y] = interpIB4(u,v,X,dx,dy);
% [Utracer,Vtracer,i1x,j1x,i1y,j1y] = interpIB4(u,v,X_tracer,dx,dy);
%Update the positions of the Lagrangian structure at the half time-step
X_half(:,1) = X(:,1) + 0.5*dt*U(:);
X_half(:,2) = X(:,2) + 0.5*dt*V(:);
X_half(:,1) = mod(X_half(:,1),cu*dx);
X_half(:,2) = mod(X_half(:,2),rv*dy);
X_half_tracer(:,1) = X_tracer(:,1) + 0.5.*dt.*Utracer(:);
X_half_tracer(:,2) = X_tracer(:,2) + 0.5.*dt.*Vtracer(:);
X_half_tracer(:,1) = mod(X_half_tracer(:,1),cu*dx);
X_half_tracer(:,2) = mod(X_half_tracer(:,2),rv*dy);
%Update the forces on the Lagrangian grid
% F_half(:,1) = -kappa*(X_half(:,1) - X_OG(:,1));
% F_half(:,2) = -kappa*(X_half(:,2) - X_OG(:,2));
% F_half(:,1) = -eta*U;
% F_half(:,2) = -eta*V;
 F_half = Elastic_force_update(X_half,kappa,ds);
% F_half = Elastic_Laplacian(X_half,kappa,ds);
% F_half = uniform_normal_force_circle(X_half,kappa,ds,r);
%Spread these forces onto the Eulerian grid
[ffx_half,ffy_half] = spreadBS3BS2(F_half,X_half,u,v,dx,dy,ds);
% [ffx_half,ffy_half] = spreadBS2BS1(F_half,X_half,u,v,dx,dy,ds);
% [ffx_half,ffy_half] = spreadIB4(F_half,X_half,u,v,dx,dy,ds);
% Form the matrices needed to solve the linear system ....
%Laplacian for u
% e = ones(cu,1); Lapu_x = spdiags([e -2*e e], -1:1,cu,cu)/(dx*dx);
% e = ones(ru,1); Lapu_y = spdiags([e -2*e e], -1:1, ru, ru)/(dy*dy);
% %Append entries to each matrix to account for the periodic boundary
% %conditions
% Lapu_x(1,end) = 1/(dx*dx);
% Lapu_x(end,1) = 1/(dx*dx);
% Lapu_y(1,end) = 1/(dy*dy);
% Lapu_y(end,1) = 1/(dy*dy);

% Ix = speye(cu); Iy = speye(ru);
% Lapu = kron(Lapu_x,Iy) + kron(Ix,Lapu_y);
% Lapu = 0.5.*mu.*Lapu;

%Create ghost cell entries for each term
[u_g,v_g] = Ghostnodesside_periodic(u,v);
%These terms are only used for the last solve
Lapu_rhs = (u_g(3:end,2:end-1) - 2.*u_g(2:end-1,2:end-1) + u_g(1:end-2,2:end-1))/(dy*dy) + ...
    (u_g(2:end-1,3:end) - 2.*u_g(2:end-1,2:end-1) + u_g(2:end-1,1:end-2))/(dx*dx);
Lapv_rhs = (v_g(3:end,2:end-1) - 2.*v_g(2:end-1,2:end-1) + v_g(1:end-2,2:end-1))/(dy*dy) + ...
    (v_g(2:end-1,3:end) - 2.*v_g(2:end-1,2:end-1) + v_g(2:end-1,1:end-2))/(dx*dx);

Lapu_rhs = (0.5)*(mu).*Lapu_rhs;
Lapv_rhs = (0.5)*(mu).*Lapv_rhs;



% %Build the Laplacian for v
% 
% e = ones(cv,1); Lapv_x = spdiags([e -2*e e], -1:1,cv,cv)/(dx*dx);
% %Edit the matrix to impose dirichlet boundart conditions
% Lapv_x(1,end) = 1/(dx*dx);
% Lapv_x(end,1) = 1/(dx*dx);
% e = ones(rv,1); Lapv_y = spdiags([e -2*e e], -1:1, rv, rv)/(dy*dy);
% Lapv_y(1,end) = 1/(dy*dy);
% Lapv_y(end,1) = 1/(dy*dy);
% Ix = speye(cv); Iy = speye(rv);
% Lapv = kron(Lapv_x,Iy) + kron(Ix,Lapv_y);
% Lapv = 0.5.*mu.*Lapv;


%Build Differentiaion matrix for the pressure in the X and Y directions
[rp,cp] = size(p);
% G_x = zeros(cp,cp);
% G_x(1:cp+1:end) = 1/dx;
% G_x(2:cp+1:end) = -1/dx;
% G_x(1,end) = -1/dx;
% Ix = eye(rp,rp);
% G_x = kron(G_x,Ix);
% 
% 
% G_y = zeros(rp,rp);
% G_y(1:rp+1:end) = 1/dy;
% G_y(2:rp+1:end) = -1/dy;
% G_y(1,end) = -1/dy;
% Iy = eye(cp,cp);
% G_y = kron(Iy,G_y);



%Build Differentiation matrices for the continuity equation
%Both are approximated via centered differences w/r/t the centerd nodes
[ru,cu] = size(u);
% Dx = zeros(cu,cu);
% Dx(1:cu+1:end) = -1/dx;
% Dx(cu+1:cu+1:end) = 1/dx;
% Dx(end,1) = 1/dx;
% Ix = speye(ru,ru);
% Dx = kron(Dx,Ix);
% 
% 
% [rv,cv] = size(v);
% Dy = zeros(rv,rv);
% Dy(1:rv+1:end) = -1/dy;
% Dy(rv+1:rv+1:end) = 1/dy;
% Dy(end,1) = 1/dy;
% Iy = speye(cv,cv);
% Dy = kron(Iy,Dy);

%Start assemblying the linear system to be used
% A1 = (rho/dt).*speye(size(Lapu));
% A2 = (rho/dt).*speye(size(Lapv));
% N1 = zeros(size(A1));
% N2 = zeros(size(A2));
% N3 = zeros(rp*cp);
N7 = zeros(rp*cp,1);
%Setup up the explicit convective terms
[u_g,v_g] = Ghostnodesside_periodic(u,v);
[u0_g,v0_g] = Ghostnodesside_periodic(u0,v0);
uce = 0.5.*(u_g(:,1:end-1) + u_g(:,2:end));
vce = 0.5.*(v_g(1:end-1,:) + v_g(2:end,:));

u_corn = 0.25.*(uce(1:end-1,1:end-1) + uce(2:end,1:end-1) + uce(1:end-1,2:end)+ uce(2:end,2:end));
v_corn = 0.25.*(vce(1:end-1,1:end-1) + vce(2:end,1:end-1) + vce(1:end-1,2:end)+ vce(2:end,2:end));
%interpolate v onto side centered components
v_sidex = 0.5.*(vce(2:end,1:end-2)+vce(2:end,2:end-1));
u_sidey = 0.5.*(uce(1:end-2,2:end) + uce(2:end-1,2:end));
Conv_u = u.*(uce(2:end-1,2:end)-uce(2:end-1,1:end-1))./dx + v_sidex.*(u_corn(2:end,:)-u_corn(1:end-1,:))./dy;
Conv_v = u_sidey.*(v_corn(:,2:end) - v_corn(:,1:end-1))./dx + v.*(vce(2:end,2:end-1) - vce(1:end-1,2:end-1))./dy;
Conv_u = -1.5.*rho.*Conv_u;
Conv_v = -1.5.*rho.*Conv_v;

u0ce = 0.5.*(u0_g(:,1:end-1) + u0_g(:,2:end));
v0ce = 0.5.*(v0_g(1:end-1,:) + v0_g(2:end,:));

u0_corn = 0.25.*(u0ce(1:end-1,1:end-1) + u0ce(2:end,1:end-1) + u0ce(1:end-1,2:end)+ u0ce(2:end,2:end));
v0_corn = 0.25.*(v0ce(1:end-1,1:end-1) + v0ce(2:end,1:end-1) + v0ce(1:end-1,2:end)+ v0ce(2:end,2:end));
%interpolate v0 onto side centered components
v0_sidex = 0.5.*(v0ce(2:end,1:end-2)+ v0ce(2:end,2:end-1));
u0_sidey = 0.5.*(u0ce(1:end-2,2:end) + u0ce(2:end-1,2:end));
Conv_u0 = u0.*(u0ce(2:end-1,2:end)-u0ce(2:end-1,1:end-1))./dx + v0_sidex.*(u0_corn(2:end,:)-u0_corn(1:end-1,:))./dy;
Conv_v0 = u0_sidey.*(v0_corn(:,2:end) - v0_corn(:,1:end-1))./dx + v0.*(v0ce(2:end,2:end-1) - v0ce(1:end-1,2:end-1))./dy;
Conv_u0 = 0.5.*rho.*Conv_u0;
Conv_v0 = 0.5.*rho.*Conv_v0;

% LHS = [A1 - Lapu, N1, G_x; N2, A2 - Lapv, G_y; -Dx, -Dy,N3];
RHS = [Conv_u(:) + Conv_u0(:) + Lapu_rhs(:) + A1*u(:) + ffx_half(:) ;...
    Conv_v(:) + Conv_v0(:) + A2*v(:) + Lapv_rhs(:) + ffy_half(:); ...
   N7(:)];
 
vals = LHS\RHS;
u_new = vals(1:ru*(cu));
u_new = reshape(u_new,ru,cu);
v_new = vals(ru*(cu)+1:ru*(cu)+(rv)*cv);
v_new = reshape(v_new,rv,cv);
p_new = vals(ru*cu + rv*cv+1:end);
p_new = reshape(p_new,rp,cp);
%Remove the average of the pressure since it should belong to L^2_0...
p_new = p_new - dx*dy*sum(p_new,'all');


% Now again update the the lagrange points 
[U,V,i1x,j1x,i1y,j1y] = interpBS3BS2(0.5.*(u_new + u),0.5.*(v_new + v),X_half,dx,dy);
[Utracer,Vtracer,i1x,j1x,i1y,j1y] = interpBS3BS2(0.5.*(u_new + u),0.5.*(v_new + v),X_half_tracer,dx,dy);
% [U,V,i1x,j1x,i1y,j1y] = interpBS5BS4(0.5.*(u_new + u),0.5.*(v_new + v),X_half,dx,dy);
% [Utracer,Vtracer,i1x,j1x,i1y,j1y] = interpBS5BS4(0.5.*(u_new + u),0.5.*(v_new + v),X_half_tracer,dx,dy);
% [U,V,i1x,j1x,i1y,j1y] = interpBS1BS0(0.5.*(u_new + u),0.5.*(v_new + v),X_half,dx,dy);
% [U,V,i1x,j1x,i1y,j1y] = interpBS2BS1(0.5.*(u_new + u),0.5.*(v_new + v),X_half,dx,dy);
% [Utracer,Vtracer,i1x,j1x,i1y,j1y] = interpBS2BS1(0.5.*(u_new + u),0.5.*(v_new + v),X_half_tracer,dx,dy);
% [U,V,i1x,j1x,i1y,j1y] = interpIB4(0.5.*(u_new + u),0.5.*(v_new + v),X_half,dx,dy);
% [Utracer,Vtracer,i1x,j1x,i1y,j1y] = interpIB4(0.5.*(u_new + u),0.5.*(v_new + v),X_half_tracer,dx,dy);
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

end
 
 







