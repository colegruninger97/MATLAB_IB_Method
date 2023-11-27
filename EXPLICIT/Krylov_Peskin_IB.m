function [u_new,v_new,p_new,F_new,X_new] = Krylov_Peskin_IB(u,v,p,F,X,dt,dx,dy,mu,kappa,ds)
%Explict implementation of the immersed boundary method following Peskin's
%lecture notes. 
rho = 1.0;
[ru,cu] = size(u);
[rv,cv] = size(v);
[rp,cp] = size(p);
u0 = mean(u(:));
v0 = mean(v(:));
%DFIB method:
% phi = solve_for_potential(u,v,dx,dy);
% [U,V] = DFIB_interp(phi,u0,v0,X,dx,dy);
% [Utracer,Vtracer] = DFIB_interp(phi,u0,v0,X_tracer,dx,dy);

%interpolate the velocity onto the Lag grid
%  [U,V,i1x,j1x,i1y,j1y] = interpBS5BS4(u,v,X,dx,dy);
%  [Utracer,Vtracer,i1x,j1x,i1y,j1y] = interpBS5BS4(u,v,X_tracer,dx,dy);
 [U,V,i1x,j1x,i1y,j1y] = interpBS3BS2(u,v,X,dx,dy);
%  [Utracer,Vtracer,i1x,j1x,i1y,j1y] = interpBS3BS2(u,v,X_tracer,dx,dy);
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
% X_half_tracer(:,1) = X_tracer(:,1) + 0.5.*dt.*Utracer(:);
% X_half_tracer(:,2) = X_tracer(:,2) + 0.5.*dt.*Vtracer(:);
% X_half_tracer(:,1) = mod(X_half_tracer(:,1),cu*dx);
% X_half_tracer(:,2) = mod(X_half_tracer(:,2),rv*dy);

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

Vol = 1.0; %is the volume of the unit torus. need to change is domain is changed.
%DFIB spread:
% [ffx_half,ffy_half] = DFIB_spread(phi,F_half,X_half,Vol,dx,dy,ds);

%Create ghost cell entries for each term
[u_g,v_g] = Ghostnodesside_periodic(u,v);
%Compute the laplacian of the previous solution
Lapu_rhs = (u_g(3:end,2:end-1) - 2.*u_g(2:end-1,2:end-1) + u_g(1:end-2,2:end-1))/(dy*dy) + ...
    (u_g(2:end-1,3:end) - 2.*u_g(2:end-1,2:end-1) + u_g(2:end-1,1:end-2))/(dx*dx);
Lapv_rhs = (v_g(3:end,2:end-1) - 2.*v_g(2:end-1,2:end-1) + v_g(1:end-2,2:end-1))/(dy*dy) + ...
    (v_g(2:end-1,3:end) - 2.*v_g(2:end-1,2:end-1) + v_g(2:end-1,1:end-2))/(dx*dx);
%Scale
Lapu_rhs = (0.5)*(mu).*Lapu_rhs;
Lapv_rhs = (0.5)*(mu).*Lapv_rhs;

%Compute the nonlinear convective terms
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

%Setup the first solve
tol = 1e-9;
maxit = 500;
%Compute the RHS of the schur complement equation
RHS = [Conv_u(:) - ffx_half(:) - (2*rho/dt).*u(:); Conv_v(:) - ffy_half(:) - (2*rho/dt).*v(:)];
RHSp = applyInvA(RHS, cu, ru, dx, dy, tol, maxit,rho,mu,dt);
RHSp = applyDivergenceU(RHSp(1:cu*ru),cu,ru,dx) + applyDivergenceV(RHSp(cv*rv+1:end),cv,rv,dy);
tol = 1e-9;
maxit = 500;
%Solve for the updated pressure
p_half = gmres(@(p)applySchurComplement(p, rp, cp, dx, dy,rho,mu,dt),RHSp,[],tol,maxit);
[Gp_halfx,Gp_halfy] = applyGradientP(p_half,cp,rp,dx,dy);
U_half = applyInvA(-RHS - [Gp_halfx(:);Gp_halfy(:)],cu, ru, dx, dy, tol, maxit,rho,mu,dt);
u_half = U_half(1:ru*cu);
u_half = reshape(u_half,[ru,cu]);
v_half = U_half(cv*rv+1:end);
v_half = reshape(v_half,[rv,cv]);


% %DFIB method....
% phi = solve_for_potential(u_half,v_half,dx,dy);
% u0 = mean(u_half(:));
% v0 = mean(v_half(:));
% % [U,V] = DFIB_interp(phi,u0,v0,X_half,dx,dy);
% [Utracer,Vtracer] = DFIB_interp(phi,u0,v0,X_half_tracer,dx,dy);


% Now again update the the lagrange points 
[U,V,i1x,j1x,i1y,j1y] = interpBS3BS2(u_half,v_half,X_half,dx,dy);
% [Utracer,Vtracer,i1x,j1x,i1y,j1y] = interpBS3BS2(u_half,v_half,X_half_tracer,dx,dy);
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
% X_tracer_new(:,1) = X_tracer(:,1) + dt*Utracer(:);
% X_tracer_new(:,2) = X_tracer(:,2) + dt*Vtracer(:);
% X_tracer_new(:,1) = mod(X_tracer_new(:,1),cu*dx);
% X_tracer_new(:,2) = mod(X_tracer_new(:,2),rv*dy);
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

%DFIB spread
% [ffx_new,ffy_new] = DFIB_spread(phi,F_new,X_new,Vol,dx,dy,ds);

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


RHS = [-(rho/dt).*u(:) - Lapu_rhs(:) + Conv_u_half(:) - ffx_new(:); -(rho/dt).*v(:) - Lapv_rhs(:) + Conv_v_half(:) - ffy_new(:)];
RHS2 = applyInvA2(RHS, cu, ru, dx, dy, tol, maxit,rho,mu,dt);
RHS2 = applyDivergenceU(RHS2(1:cu*ru),cu,ru,dx) + applyDivergenceV(RHS2(cv*rv+1:end),cv,rv,dy);
p_new = gmres(@(p)applySchurComplement2(p, cp, rp, dx, dy,rho,mu,dt),RHS2,[],tol,maxit);
[Gpnewx,Gpnewy] = applyGradientP(p_new,cp,rp,dx,dy);
U_new = applyInvA2(-RHS - [Gpnewx(:);Gpnewy(:)], cu, ru, dx, dy, tol, maxit,rho,mu,dt);
u_new = U_new(1:ru*(cu));
u_new = reshape(u_new,ru,cu);
v_new = U_new(ru*(cu)+1:ru*(cu)+(rv)*cv);
v_new = reshape(v_new,rv,cv);
p_new = reshape(p_new,rp,cp);

end
 
 function Sp = applySchurComplement(p, nx, ny, dx, dy,rho,mu,dt)
    % Apply the Schur complement to the pressure field p
    % S = -B A^-1 B^T

    % Step 1: Apply B to p (Gradient of pressure)
    [Gpx, Gpy] = applyGradientP(p, nx, ny, dx, dy);

    % Step 2: Combine Gpx and Gpy and apply A^-1 (Inverse of Laplacian)
    Gp_combined = [Gpx; Gpy];
    tol = 1e-9;
    maxit = 1000;
    invAGp = applyInvA(Gp_combined, nx, ny, dx, dy, tol, maxit,rho,mu,dt);
    % Split the result back into u and v components
    invAGpx = invAGp(1:nx*ny);
    invAGpy = invAGp(nx*ny+1:end);

    % Step 3: Apply B^T to the result (Divergence)
    divU = applyDivergenceU(invAGpx, nx, ny, dx);
    divV = applyDivergenceV(invAGpy, nx, ny, dy);

    % Combine the divergence results
    Sp = -(divU + divV);
 end

 function Sp = applySchurComplement2(p, nx, ny, dx, dy,rho,mu,dt)
    % Apply the Schur complement to the pressure field p
    % S = -B A^-1 B^T

    % Step 1: Apply B to p (Gradient of pressure)
    [Gpx, Gpy] = applyGradientP(p, nx, ny, dx, dy);

    % Step 2: Combine Gpx and Gpy and apply A^-1 (Inverse of Laplacian)
    Gp_combined = [Gpx; Gpy];
    tol = 1e-9;
    maxit = 1000;
    invAGp = applyInvA2(Gp_combined, nx, ny, dx, dy, tol, maxit,rho,mu,dt);

    % Split the result back into u and v components
    invAGpx = invAGp(1:nx*ny);
    invAGpy = invAGp(nx*ny+1:end);

    % Step 3: Apply B^T to the result (Divergence)
    divU = applyDivergenceU(invAGpx, nx, ny, dx);
    divV = applyDivergenceV(invAGpy, nx, ny, dy);

    % Combine the divergence results
    Sp = -(divU + divV);
end
 
 
function Lu = applyLaplacianU(u, nx, ny, dx, dy)
    % Reshape u into a 2D grid
    U = reshape(u, [ny, nx]);

    % Vectorized Laplacian with periodic boundaries
    % Handle x-boundary (periodic)
    Ux = [U(:, end), U, U(:, 1)]; % Extend in x-direction for periodicity
    % Handle y-boundary (assumed to be no-slip or similar)
    Uy = [U(end, :); U; U(1, :)]; % Duplicate first and last rows in y-direction

    % Central differences in x and y directions
    Lx = (Ux(:, 3:end) - 2*U + Ux(:, 1:end-2)) / dx^2;
    Ly = (Uy(3:end, :) - 2*U + Uy(1:end-2, :)) / dy^2;

    % Combine the Laplacian components
    Lu = Lx + Ly;

    % Flatten the result
    Lu = Lu(:);
end

function [Gpx, Gpy] = applyGradientP(p, nx, ny, dx, dy)
    P = reshape(p, [ny, nx]);
    Px = [P(:,end),P];
    Py = [P(end,:);P];
    % Vectorized gradient computation with periodic boundaries
    Gpx = (Px(:,2:end) - Px(:,1:end-1)) / dx; % Gradient in x-direction
    Gpy = (Py(2:end,:) - Py(1:end-1,:)) / dy; % Gradient in y-direction
    % Flatten the output
    Gpx = Gpx(:);
    Gpy = Gpy(:);
end

function x = applyInvA(b, nx, ny, dx, dy, tol, maxit,rho,mu,dt)
    % Wrapper function for an iterative solver (CG) to approximate A^-1 * b
    % b is the input vector
    % nx, ny, dx, dy are grid parameters
    % tol is the tolerance for the solver
    % maxit is the maximum number of iterations

    % Use MATLAB's pcg function to solve A*x = b
    maxit = 100;
    tol = 1e-10;
    [x, flag, relres, iter] = minres(@(u)applyA(u, nx, ny, dx, dy,rho,dt,mu), b, tol, maxit,[],[]);
    % Check if the solver converged
    if flag ~= 0
        warning('Minres did not converge. Flag: %d, Residual: %e, Iterations: %d', flag, relres, iter);
    end
end

function x = applyInvA2(b, nx, ny, dx, dy, tol, maxit,rho,mu,dt)
    % Wrapper function for an iterative solver (CG) to approximate A^-1 * b
    % b is the input vector
    % nx, ny, dx, dy are grid parameters
    % tol is the tolerance for the solver
    % maxit is the maximum number of iterations

    % Use MATLAB's pcg function to solve A*x = b
    maxit = 100;
    tol = 1e-10;
    [x, flag, relres, iter] = minres(@(u)applyA2(u, nx, ny, dx, dy,rho,dt,mu), b, tol, maxit,[],[]);
    % Check if the solver converged
    if flag ~= 0
        warning('Minres did not converge. Flag: %d, Residual: %e, Iterations: %d', flag, relres, iter);
    end
end

function Lu = applyA2(u, nx, ny, dx, dy,rho,dt,mu)
    % This function combines the Laplacian for u and v components
    % Split u into u-velocity and v-velocity components
    uComp = u(1:nx*ny);
    vComp = u(nx*ny+1:end);

    % Apply the Laplacian to each component
    LuU = applyLaplacianU(uComp, nx, ny, dx, dy);
    LuV = applyLaplacianU(vComp, nx, ny, dx, dy);
    Ax = (rho/dt).*uComp - 0.5*mu.*LuU;
    Ay = (rho/dt).*vComp - 0.5*mu.*LuV;

    % Combine the results
    Lu = [Ax(:); Ay(:)];
end


function Lu = applyA(u, nx, ny, dx, dy,rho,dt,mu)
    % This function combines the Laplacian for u and v components
    % Split u into u-velocity and v-velocity components
    uComp = u(1:nx*ny);
    vComp = u(nx*ny+1:end);

    % Apply the Laplacian to each component
    LuU = applyLaplacianU(uComp, nx, ny, dx, dy);
    LuV = applyLaplacianU(vComp, nx, ny, dx, dy);
    Ax = (2*rho/dt).*uComp - mu.*LuU;
    Ay = (2*rho/dt).*vComp - mu.*LuV;

    % Combine the results
    Lu = [Ax; Ay];
end

function divU = applyDivergenceU(u, nx, ny, dx)
    % Reshape u into a 2D grid
    U = reshape(u, [ny, nx]);
    Ux = [U,U(:,1)];
    % Compute the divergence in x-direction
    % Account for periodic boundary conditions
    divU = (Ux(:,2:end) - Ux(:,1:end-1)) / dx;

    % Flatten the result
    divU = divU(:);
end

function divV = applyDivergenceV(v, nx, ny, dy)
    % Reshape v into a 2D grid
    V = reshape(v, [ny, nx]);
    Vy = [V;V(1,:)];

    % Compute the divergence in y-direction
    % Account for periodic boundary conditions
    divV = (Vy(2:end,:) - Vy(1:end-1,:)) / dy;

    % Flatten the result
    divV = divV(:);
end






