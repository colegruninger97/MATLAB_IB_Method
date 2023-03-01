% Main driver for Poisson equation solver

% The main idea is to try and enforce Neumann Boundary condition
% constraints and Dirichlet boundary conditions constaints using the method
% of Lagrange multipliers and Peskin's Immersed boundary kernels

mfacs = 0.5:0.05:3;
conditioning = zeros(1,length(mfac));
for i = 1:length(mfacs)

L = 1.0;
N = 64;
dx = 1/N;
dy = dx;

x = 0:dx:L;
y = 0:dy:L;
xc = 0.5*(x(1:end-1) + x(2:end));
yc =  0.5*(y(1:end-1) + y(2:end));

XYc = meshgrid(x,yc);
XcY = meshgrid(xc,y);

%Define the boundary condition


%Define the curve in 2D space which can be differentiated analytically
r = 0.3;
mfac = mfacs(i);
dX = mfac*dx/r;

X_1 = 0.5*ones(1,length(0:dX:2*pi-dX)) + r.*cos(0:dX:2*pi-dX);
X_2 = 0.5*ones(1,length(0:dX:2*pi-dX)) +r.*sin(0:dX:2*pi-dX);
NIB = length(X_1);
LAMBDA=zeros(NIB,1);
X=[X_1;X_2];
X=X';
u = zeros(N+1,N+1); %This is the field we are solving for. We are only (for now) solving for the steady state equations

fun = @(x,y) 3 + 0.*x.*y;
[xx,yy] = meshgrid(x,y);
f = fun(xx,yy);%This is the forcing function
f = f(2:end-1,2:end-1);


% Use the kronecker product to form the Laplacian matrix

[ru,cu] = size(u);

e = ones(ru-2,1); Lap_y = spdiags([e -2*e e], -1:1, ru-2,ru-2)./(dy*dy);

e = ones(cu-2,1); Lap_x = spdiags([e -2*e e], -1:1, cu-2, cu-2)./(dx*dx);


Ix = speye(cu-2); Iy = speye(ru-2);
Lap = kron(Lap_x,Iy) + kron(Ix,Lap_y);

%Formulate the RHS correctly
%For now just employ homogenous Dirichlet boundary conditions
gl = @(y) 0 + 0.*y;
gr = @(y) 0 + 0.*y;
gt = @(x) 0 + 0.*x;
gb = @(x) 0 + 0.*x;

RHS1 = zeros(N-1,N-1);
RHS1(:,1) = gl(y(2:end-1))./(dy*dy);
RHS1(:,end) = gr(y(2:end-1))./(dy*dy);
RHS1(1,:) = gb(x(2:end-1))./(dx*dx);
RHS1(end,:) = gt(x(2:end-1))./(dx*dx);



%Need to augment the system now to include lagrange multipliers
%First we construct Lagrange multipliers to implement homogenous dirichlet
%boundary conditions

S = zeros((N-1)*(N-1),NIB);
S_t = S';
C = mfac/(dy);
%Now form the LHS of the equations
LHS = [-Lap + 1.*speye(length(diag(Lap))), S; S_t, zeros(NIB,NIB)];
i1 = zeros(NIB,4);
j1 = zeros(NIB,4);
for k = 1:NIB
    sx = X(k,1)/dx;
    sy = X(k,2)/dy;
    I = floor(sx);
    J = floor(sy);
    i1(k,:) = (J-1):(J+2); %These are the indices on the eulerian grid which support the discrete delta function
    j1(k,:) = (I-1):(I+2);
    rx = sx-I;
    ry = sy-J;
    ws(:,:) = IB4_1(rx).*IB4_2(ry);
    for q=1:length(i1(k,:))
        for l=1:length(j1(k,:))
            LHS((j1(k,l)-2)*(N-1) + i1(k,q), (N-1)*(N-1)+k) = LHS((j1(k,l)-2)*(N-1) + i1(k,q), (N-1)*(N-1)+k) - C.*ws(q,l);
        end
    end
    
    for q=1:length(i1(k,:))
        for l=1:length(j1(k,:))
            LHS((N-1)*(N-1) + k, (j1(k,l)-2)*(N-1) + i1(k,q)) = LHS((N-1)*(N-1) + k, (j1(k,l)-2)*(N-1) + i1(k,q)) + ws(q,l);
        end
    end


%Now reconstruct u everyone using the boundary conditions
end
%RHS = [RHS(:);zeros(NIB,1)];

%Form the schur complement
S = LHS(1:(N-1)*(N-1),(N-1)*(N-1)+1:end);
S_t = LHS((N-1)*(N-1)+1:end,1:(N-1)*(N-1));
A = -Lap + 1.*speye(length(diag(Lap)));
RHS = -S_t*inv(A)*f(:) - S_t*inv(A)*RHS1(:);
Schur = S_t*inv(A)*S;
Schur = full(Schur);
conditioning(1,i) = cond(Schur);
% LAMBDA = Schur\RHS;

% u_int = vals(1:(N-1)*(N-1));
% LAMBDA = vals((N-1)*(N-1)+1:end);
% u_int = reshape(u_int,N-1,N-1);

% u_int = inv(A)*S*LAMBDA + inv(A)*f(:);
% u_int = reshape(u_int,N-1,N-1);

u(2:end-1,2:end-1) = u_int;
u(:,1) = gl(y);
u(:,end) = gr(y);
u(1,:) = gb(x);
u(end,:) = gt(x);
Schur = full(Schur);
%Form the schur complement
S = LHS(1:(N-1)*(N-1),(N-1)*(N-1)+1:end);
S_t = LHS((N-1)*(N-1)+1:end,1:(N-1)*(N-1));
S = full(S);


%Take the solution u and perform centered differences
u_x = (-u(:,1:end-1) + u(:,2:end))/dx;
u_y = (-u(1:end-1,:) + u(2:end,:))/dy;

Ux = 0.5.*(u_x(1:end-1,:) + u_x(2:end,:));
Uy = 0.5.*(u_y(:,1:end-1) + u_y(:,2:end));
surf(u)
% quiver(xc,yc,Ux,Uy)
% 
 end




%% Now try to implement the Homogenous Neumann boundary conditions

% Main driver for Poisson equation solver

% The main idea is to try and enforce Neumann Boundary condition
% constraints and Dirichlet boundary conditions constaints using the method
% of Lagrange multipliers and Peskin's Immersed boundary kernels
L = 1.0;
N = 128;
dx = 1/N;
dy = dx;

x = 0:dx:L;
y = 0:dy:L;

%Define the boundary condition 


%Define the curve in 2D space which can be differentiated analytically
r = 0.25;
mfac = 0.2;
dX = mfac*dx/r;

X_1 = 0.5*ones(1,length(0:dX:2*pi-dX)) + r.*cos(0:dX:2*pi-dX);
X_2 = 0.5*ones(1,length(0:dX:2*pi-dX)) + r.*sin(0:dX:2*pi-dX);
NIB = length(X_1);
LAMBDA=zeros(NIB,1);
X=[X_1;X_2];
X=X';
u = zeros(N+1,N+1); %This is the field we are solving for. We are only (for now) solving for the steady state equations

fun = @(x,y)  1+0.*x.*y;
[xx,yy] = meshgrid(x,y);
f = fun(xx,yy);%This is the forcing function
f = f(2:end-1,2:end-1);


% Use the kronecker product to form the Laplacian matrix

[ru,cu] = size(u);

e = ones(ru-2,1); Lap_y = spdiags([e -2*e e], -1:1, ru-2,ru-2)./(dy*dy);

e = ones(cu-2,1); Lap_x = spdiags([e -2*e e], -1:1, cu-2, cu-2)./(dx*dx);


Ix = speye(cu-2); Iy = speye(ru-2);
Lap = kron(Lap_x,Iy) + kron(Ix,Lap_y);

%Formulate the RHS 
%For now just employ homogenous Dirichlet boundary conditions
gl = @(y) 0 + 0.*y;
gr = @(y) 0 + 0.*y;
gt = @(x) 0 + 0.*x;
gb = @(x) 0 + 0.*x;

RHS = zeros(N-1,N-1);
RHS(:,1) = gl(y(2:end-1))./(dy*dy);
RHS(:,end) = gr(y(2:end-1))./(dy*dy);
RHS(1,:) = gb(x(2:end-1))./(dx*dx);
RHS(end,:) = gt(x(2:end-1))./(dx*dx);
RHS = RHS(:,:) + f(:,:); %add in the forcing terms


%Need to augment the system now to include lagrange multipliers
%First we construct Lagrange multipliers to implement homogenous dirichlet
%boundary conditions

S = zeros((N-1)*(N-1),NIB);
S_t = S';
C = dX/(dx*dx*dy);
%Now form the LHS of the equations
LHS = [-Lap, S; S_t, zeros(NIB,NIB)];

for k = 1:NIB
    %compute the normal derivative for a particular value of k
    nx = cos(2*pi*dX*(k-1));
    ny = sin(2*pi*dX*(k-1));
    sx = X(k,1)/dx;
    sy = X(k,2)/dy;
    I = floor(sx);
    J = floor(sy);
    i1(k,:) = (J-1):(J+2); %These are the indices on the eulerian grid which support the discrete delta function
    j1(k,:) = (I-1):(I+2);
    rx = sx-I;
    ry = sy-J;
    ws(:,:) = nx.*sin_kerx(rx).*cos_kery(ry) + ny.*cos_kerx(rx).*sin_kery(ry);
    for q=1:length(i1(k,:))
        for l=1:length(j1(k,:))
            LHS((j1(k,l)-2)*(N-1) + i1(k,q), (N-1)*(N-1)+k) = LHS((j1(k,l)-2)*(N-1) + i1(k,q), (N-1)*(N-1)+k) - C.*ws(q,l);
        end
    end
    
    for q=1:length(i1(k,:))
        for l=1:length(j1(k,:))
            LHS((N-1)*(N-1) + k, (j1(k,l)-2)*(N-1) + i1(k,q)) = LHS((N-1)*(N-1) + k, (j1(k,l)-2)*(N-1) + i1(k,q)) + (1/dx).*ws(q,l);
        end
    end



%Now reconstruct u everyone using the boundary conditions
end
RHS = [RHS(:);0.*ones(NIB,1)];
vals = LHS\RHS;

u_int = vals(1:(N-1)*(N-1));
LAMBDA = vals((N-1)*(N-1)+1:end);
u_int = reshape(u_int,N-1,N-1);

u(2:end-1,2:end-1) = u_int;
u(:,1) = gl(y);
u(:,end) = gr(y);
u(1,:) = gb(x);
u(end,:) = gt(x);

u_x = (-u(:,1:end-1) + u(:,2:end))/dx;
u_y = (-u(1:end-1,:) + u(2:end,:))/dy;

Ux = 0.5.*(u_x(1:end-1,:) + u_x(2:end,:));
Uy = 0.5.*(u_y(:,1:end-1) + u_y(:,2:end));
u_xx = (u(:,1:end-2) + u(:,3:end) - 2.*u(:,2:end-1))./(dx.*dx);

surf(u)















