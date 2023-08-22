function [u_new,v_new,p_new,F1_new,F2_new,LHS] = BE_IB_Evolve_Periodic(u,v,p,F1,F2,dt,dx,dy,N,NIB,mu,rho,kappa,mfac,x,y,X)

%Build matrices used for the left hand side of the equation
[u_g,v_g] = Ghostnodesside_periodic(u,v);
%Laplacian for u
[ru,cu] = size(u);
e = ones(cu,1); Lapu_x = spdiags([e -2*e e], -1:1,cu,cu)/(dx*dx);
%Edit the matrix to impose the periodic boundary conditions
Lapu_x(1,end-1) = 1/(dx*dx);
Lapu_x(end,2) = 1/(dx*dx);

e = ones(ru,1); Lapu_y = spdiags([e -2*e e], -1:1, ru, ru)/(dy*dy);
Lapu_y(1,end) = 1/(dy*dy);
Lapu_y(end,1) = 1/(dy*dy);

Ix = speye(cu); Iy = speye(ru);
Lapu = kron(Lapu_x,Iy) + kron(Ix,Lapu_y);
Lapu= mu.*Lapu; %Matrix for Laplacian defined on x-edges


%Build the Laplacian for v
[rv,cv] = size(v);
e = ones(cv,1); Lapv_x = spdiags([e -2*e e], -1:1,cv,cv)/(dx*dx);
%Edit the matrix to impose the periodic boundary conditions
Lapv_x(1,end) = 1/(dx*dx);
Lapv_x(end,1) = 1/(dx*dx);

e = ones(rv,1); Lapv_y = spdiags([e -2*e e], -1:1, rv, rv)/(dy*dy);
Lapv_y(1,end-1) = 1/(dy*dy);
Lapv_y(end,2) = 1/(dy*dy);

Ix = speye(cv); Iy = speye(rv);
Lapv = kron(Lapv_x,Iy) + kron(Ix,Lapv_y);
Lapv= mu.*Lapv;

%Build Differentiaion matrix for the pressure in the X and Y directions
[rp,cp] = size(p);
%Initialize
%this is for the x-derivative of the pressure located at the x-edges
G_x = zeros(cp+1,cp);
G_x(1:cp+2:end) = 1;
G_x(2:cp+2:end) = -1;
G_x(1,end) = -1;
G_x(end,1) = 1;
G_x = (1/dx).*G_x; %Scale appropriately
Ix = speye(cp,cp);
G_x = kron(G_x,Ix);
G_x = sparse(G_x);


%Now for the y-derivative
G_y = zeros(rp+1,rp);
G_y(1:rp+2:end) = 1;
G_y(2:rp+2:end) = -1;
G_y(1,end) = -1;
G_y(end,1) = 1;
G_y = (1/dy).*G_y; %Scale appropriately
Iy = speye(rp,rp);
G_y = kron(Iy,G_y);
G_y = sparse(G_y);


%Build Differentiation matrices for the continuity equation
%Both are approximated via centered differences w/r/t the centerd nodes
[ru,cu] = size(u);
Dx = zeros(size(u));
Dx(1:ru+1:end) = -1;
Dx(ru+1:ru+1:end) = 1;
Ix = speye(ru);
Dx = kron(Dx,Ix);
Dx = (1/dx).*Dx;

[rv,cv] = size(v);
Dy = zeros(cv,rv);
Dy(1:cv+1:end) = -1;
Dy(cv+1:cv+1:end) = 1;
Iy = speye(cv);
Dy = kron(Iy,Dy);
Dy = (1/dy).*Dy;




A1 = (rho/dt).*speye(size(Lapu));
A2 = (rho/dt).*speye(size(Lapv));
N1 = zeros(size(A1));
N2 = zeros(size(A2));
N3 = zeros(length(p)^2);
N7 = zeros(length(p)^2,1);


S_x = zeros(N*(N+1),2*NIB);
S_y = zeros(N*(N+1),2*NIB);
S_x1 = S_x';
S_y1 = S_y';

N4 = zeros((N*N),2*NIB);
N5 = zeros(2*NIB,N*N);
N6 = zeros(2*NIB,2*NIB);

dX = mfac*dx;
C = dX/(dx*dy);

uce = 0.5.*(u_g(:,1:end-1) + u_g(:,2:end));
vce = 0.5.*(v_g(1:end-1,:) + v_g(1:end-1,:));

u_corn = 0.25.*(uce(1:end-1,1:end-1) + uce(2:end,1:end-1) + uce(1:end-1,2:end)+ uce(2:end,2:end));
v_corn = 0.25.*(vce(1:end-1,1:end-1) + vce(2:end,1:end-1) + vce(1:end-1,2:end)+ vce(2:end,2:end));

uuce = uce.*uce;
vvce = vce.*vce;
uv_corn = u_corn.*v_corn;
%Compute the nonlinear convective terms (central dofferences for now) then
%I can optionally add them in ...

Conv_u = (uuce(2:end-1,2:end) - uuce(2:end-1,1:end-1))/dx;
Conv_u = Conv_u + (uv_corn(2:end,:) - uv_corn(1:end-1,:))/dy;

Conv_v = (vvce(2:end,2:end-1) - vvce(1:end-1,2:end-1))/dy;
Conv_v = Conv_v + (uv_corn(:,2:end)-uv_corn(:,1:end-1))/dx;


%Setup the linear system to solve
LHS = [A1 - Lapu, N1, G_x, S_x; N2, A2 - Lapv, G_y, S_y; -Dx, -Dy,N3,N4;S_x1,S_y1,N5,N6];
%Now fill the relevant entries of this matrix with the interpolation and
%spreading operators....
Ixx = zeros(NIB,2);
rx = Ixx;
ry = Ixx;
Iyy = Ixx;
j1x = zeros(NIB,4); %here 4 is the width of the support of the current IB kernel
j1y = j1x;
i1x = j1x;
i1y = j1y;
wsx = zeros(4,4,NIB);
wsy = wsx;

for k = 1:NIB
    sxx = X(k,1)/dx;
    sxy = X(k,2)/dy - 1/2;
    sx = [sxx,sxy];
    syy = X(k,2)/dy;
    syx = X(k,1)/dx - 1/2;
    sy = [syx,syy];
    Ixx(k,:) = ceil(sx); %Revelant indices for the forcing on the eulerian grid
    Iyy(k,:) = ceil(sy);
    j1x(k,:) = (Ixx(k,1)-1):(Ixx(k,1)+2); %These are the indices on the eulerian grid which support the discrete delta function
    i1x(k,:) = (Ixx(k,2)-1):(Ixx(k,2)+2);
    j1y(k,:) = (Iyy(k,1)-1):(Iyy(k,1)+2);
    i1y(k,:) = (Iyy(k,2)-1):(Iyy(k,2)+2);
    rx(k,:) = (Ixx(k,:)-sx);
    ry(k,:) = (Iyy(k,:)-sy);
    wsx(:,:,k) = IB4_1(rx(k,1)).*IB4_2(rx(k,2));
    wsy(:,:,k) = IB4_1(ry(k,1)).*IB4_2(ry(k,2));
    for q = 1:length(j1x(k,:))
        for n = 1:length(i1x(k,:))
            %Fill in for the u terms
            LHS((j1x(k,q)-1)*N + i1x(k,n),N*(N+1) + N*(N+1) + N*N +k) ...
                = LHS((j1x(k,q)-1)*N + i1x(k,n),N*(N+1) + N*(N+1) + N*N +k) - C.*wsx(q,n,k);
            
            %Fill in for the v terms
            
            LHS(N*(N+1) + (j1y(k,q)-1)*(N+1) + (i1y(k,n)),N*(N+1) + N*(N+1) + N*N + NIB + k) ...
                = LHS(N*(N+1) + (j1y(k,q)-1)*(N+1) + (i1y(k,n)),N*(N+1) + N*(N+1) + N*N + NIB + k) - C.*wsy(q,n,k);
        end
    end
    % Now compute the interpolation operator which projects U onto the
    % immersed structure
    for q = 1:length(j1x(k,:))
        for n = 1:length(i1x(k,:))
    LHS(N*(N+1) + N*(N+1) + N*N + k,(j1x(k,q)-1)*N + i1x(k,n)) = ...
        LHS(N*(N+1) + N*(N+1) + N*N + k,(j1x(k,q)-1)*N + i1x(k,n)) + dt*kappa.*wsx(q,n,k);
    
    LHS(N*(N+1) + N*(N+1) + N*N + NIB + k, N*(N+1) + (j1y(k,q)-1)*(N+1) + (i1y(k,n))) = ...
        LHS(N*(N+1) + N*(N+1) + N*N + NIB + k, N*(N+1) + (j1y(k,q)-1)*(N+1) + (i1y(k,n))) + dt*kappa.*wsy(q,n,k);
    
    
        end
    end
   
    
    
    
    
end
% spy(LHS)
% Now create the RHS to invert the LHS on ...
RHS = [A1*u(:) - Conv_u(:);...
    A2*v(:) - Conv_v(:);...
    N7(:);
    F1(:);
    F2(:)];
% [V,D] = eig(LHS);
% rank(LHS)
% PP = null(LHS);
%Use the black box... 
vals = LHS\RHS;
u_new = vals(1:ru*(cu));
u_new = reshape(u_new,ru,cu);
v_new = vals(ru*(cu)+1:ru*(cu)+(rv)*cv);
v_new = reshape(v_new,rv,cv);
p_new = vals(ru*(cu) + (rv)*cv + 1:ru*(cu) + (rv)*cv + N*N);
p_new = p_new - dx*dy*sum(p_new);
F1_new = vals(ru*(cu) + (rv)*cv + N*N + 1: ru*(cu) + (rv)*cv + N*N + NIB);
F2_new = vals(ru*(cu) + (rv)*cv + N*N + NIB + 1:ru*(cu) + (rv)*cv + N*N + 2*NIB);
%Normalize the pressure: the gradient is really all that matters
p_new = reshape(p_new,rp,cp);

end
