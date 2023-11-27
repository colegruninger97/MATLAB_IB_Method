clear;
clc;
nx = 32;
ny = 32;
dx = 1/nx;
dy = 1/ny;
x = 0:dx:1;
y = x;
x_side = x(1:end-1);
y_side = x_side;
xc = 0.5.*(x(2:end) + x(1:end-1));
yc = 0.5.*(y(2:end) + y(1:end-1));
[Xside,Y] = meshgrid(x_side,yc);
[X,Yside] = meshgrid(xc,y_side);
mu = 1;
utest = @(x,y) sin(2*pi.*x).*sin(2*pi.*y);
fx_test = @(x,y) mu*8*pi*pi.*cos(2*pi.*x).*sin(2*pi.*y) + (pi/2).*sin(2*pi.*x);
fy_test = @(x,y) -mu*8*pi*pi.*cos(2*pi.*y).*sin(2*pi.*x)  + (pi/2).*sin(2*pi.*y);
Fx = fx_test(Xside,Y);
Fy = fy_test(X,Yside);
F = [Fx(:);Fy(:)];
ux = utest(Xside,Y);
uy = utest(X,Yside);
ucent = utest(X,Y);
DU = applyDivergenceU(ux(:),nx,ny,dx);
DV = applyDivergenceV(uy(:),nx,ny,dy);
LU = applyLaplacianU(ucent(:),nx,ny,dx,dy);
Eu = abs(reshape(LU,[ny,nx]) + sin(2*pi.*Y).*sin(2*pi.*X).*8*pi*pi);
surf(Eu)
% G = applyInvLap([U;V],nx,ny,dx,dy);
% Gx = G(1:nx*ny);
% Gy = G(nx*ny+1:end);
% surf(reshape(Gx,[ny,nx]))
%Setup the right handside for the schur complement equation
tol = 1e-3; % Tolerance for the solver
maxit = 100; % Maximum iterations
g = applyInvLap(F(:), nx, ny, dx, dy, tol, maxit); % g is the velocity field
Dg = applyDivergenceU(g(1:nx*ny),nx,ny,dx) + applyDivergenceV(g(nx*ny+1:end),nx,ny,dy);
%surf(reshape(Dg,[16,16]))
[p,flag, relres, iter, resvec] = gmres(@(p)applySchurComplement(p, nx, ny, dx, dy), Dg,[],tol,maxit);
surf(reshape(p,[ny,nx]))
%Check output
if flag ~= 0
    warning('GMRES did not converge. Flag: %d, Residual: %e, Iterations: %d', flag, relres);
end
[Gpx,Gpy] = applyGradientP(p,nx,ny,dx,dy);
Gp = [Gpx;Gpy];
U = (1/mu).*applyInvLap(Gp + F(:),nx,ny,dx,dy);
u = U(1:nx*ny);
v = U(nx*ny+1:end);
Error = abs(reshape(v,[ny,nx]) - cos(2*pi.*Yside).*sin(2*pi.*X));
surf(Error)
%%

% Example parameters
nx = 64;
ny = 64;
dx = 1/nx;
dy = 1/ny;
xx = 0:dx:1;
yy = 0:dy:1;
xside = 0:dx:1-dx;
yside = xside;
xcent = dx/2:dx:1-dx/2;
ycent = xcent;
[X,Y] = meshgrid(xcent,ycent);
testb = @(x,y) sin(2*pi.*x).*sin(2*pi*y);
truesoln = @(x,y) -(1/(8*pi*pi)).*sin(2*pi.*x).*sin(2*pi*y);
b = testb(X,Y);
b = b(:);
x0 = zeros(nx*ny, 1);  % Initial guess

% Reshape for 2D operation
x0_2d = reshape(x0, nx, ny);
b_2d = reshape(b, nx, ny);

% Run the Jacobi preconditioner
x_after_jacobi = jacobiPreconditioner(x0_2d, b_2d, nx, ny, dx, dy, 1);
x_true = truesoln(X,Y);
error = abs(x_after_jacobi - x_true(:));
max(error)
surf(reshape(error,ny,nx))

function Sp = applySchurComplement(p, nx, ny, dx, dy)
    % Apply the Schur complement to the pressure field p
    % S = -B A^-1 B^T

    % Step 1: Apply B to p (Gradient of pressure)
    [Gpx, Gpy] = applyGradientP(p, nx, ny, dx, dy);

    % Step 2: Combine Gpx and Gpy and apply A^-1 (Inverse of Laplacian)
    Gp_combined = [Gpx; Gpy];
    tol = 1e-8;
    maxiter = 1000;
    invAGp = applyInvLap(Gp_combined, nx, ny, dx, dy,tol,maxiter);

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

function x = applyInvLap(b, nx, ny, dx, dy, tol, maxit)
    % Wrapper function for an iterative solver (CG) to approximate A^-1 * b
    % b is the input vector
    % nx, ny, dx, dy are grid parameters
    % tol is the tolerance for the solver
    % maxit is the maximum number of iterations

    % Use MATLAB's pcg function to solve A*x = b
    maxit = 100;
    tol = 1e-10;
    [x, flag, relres, iter] = minres(@(u)applyLaplacian(u,nx,ny,dx,dy), b, tol, maxit,[],[]);
    x = x - mean(x);
    % Check if the solver converged
    if flag ~= 0
        warning('Minres did not converge. Flag: %d, Residual: %e, Iterations: %d', flag, relres, iter);
    end
end

function Lu = applyLaplacian(u, nx, ny, dx, dy)
    % This function combines the Laplacian for u and v components
    % Split u into u-velocity and v-velocity components
    uComp = u(1:nx*ny);
    vComp = u(nx*ny+1:end);

    % Apply the Laplacian to each component
    LuU = applyLaplacianU(uComp, nx, ny, dx, dy);
    LuV = applyLaplacianU(vComp, nx, ny, dx, dy);

    % Combine the results
    Lu = [LuU; LuV];
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


function x_new = jacobiPreconditioner(x, b, nx, ny, dx, dy, iterations)
    % x is the initial guess (2D matrix)
    % b is the right-hand side vector (2D matrix)
    % nx, ny are the number of grid points in x and y directions
    % dx, dy are the grid spacings in x and y directions
    % iterations is the number of Jacobi iterations to perform
    
    % Precompute constants
    dx2 = dx * dx;
    dy2 = dy * dy;
    denom = -2 * (1/dx2 + 1/dy2);

    for iter = 1:iterations
        % Compute the update (note: using the 'old' values of x)
        x_new = -(circshift(x, [0, -1]) + circshift(x, [0, 1])) / dx2 - ...
                (circshift(x, [-1, 0]) + circshift(x, [1, 0])) / dy2 + ...
                b;
        x_new = x_new / denom;

        % Copy the updated values to x for the next iteration
        x = x_new;
    end

    % Flatten x_new to a vector for consistency with linear solvers
    x_new = x_new(:);
end



