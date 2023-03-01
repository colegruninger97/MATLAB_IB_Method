function [u_new,v_new,p_new,F1_new,F2_new,LHS,A,B1,B2] = BE_IB_Evolve_Phys(u,v,p,F1,F2,dt,dx,dy,N,NIB,mu,rho,kappa,mfac,x,y,X)
persistent Lapu_y Lapu_x Lapu Lapv_y Lapv_x Lapv G_x G_y Dx Dy A1 A2 N1 N2 N3 N7 N5 N6 S_x S_y S_x1 N4 S_y1
Px = -0.5;
y_c = 0.5.*(y(1:end-1) + y(2:end));
x_c = 0.5.*(x(1:end-1) + x(2:end));
u_exact = @(x,y,t) 0.*x + (Px/(2*mu)).*y.*(y-1);
v_exact = @(x,y,t) 0.*x.*y;
p_exact = @(x,y,t) 0.*y + Px.*x - Px/4;
u_int = u(:,2:end-1);
v_int = v(2:end-1,:);
%Laplacian for u
[ru,cu] = size(u);
e = ones(cu-2,1); Lapu_x = spdiags([e -2*e e], -1:1,cu-2,cu-2)/(dx*dx);
e = ones(ru,1); Lapu_y = spdiags([e -2*e e], -1:1, ru, ru)/(dy*dy);
%implement dirichlet boundary conditions at the walls via quadtratic
%extrapolation
 Lapu_y(1,1) = -4/(dy*dy);
 Lapu_y(end,end) = -4/(dy*dy);
 Lapu_y(1,2) = (4/3)/(dy*dy);
 Lapu_y(end,end-1) = (4/3)/(dy*dy);
%instead throw boundary terms on the RHS... hopefully this will reduce the
%error
Ix = speye(cu-2); Iy = speye(ru);
Lapu = kron(Lapu_x,Iy) + kron(Ix,Lapu_y);
Lapu = mu.*Lapu;



%Build the Laplacian for v
[rv,cv] = size(v);
e = ones(cv,1); Lapv_x = spdiags([e -2*e e], -1:1,cv,cv)/(dx*dx);
%Edit the matrix to impose dirichlet boundart conditions
Lapv_x(1,1) = -4/(dx*dx);
Lapv_x(end,end) = -4/(dx*dx);
Lapv_x(1,2) = (4/3)/(dx*dx);
Lapv_x(end,end-1) = (4/3)/(dx*dx);
e = ones(rv-2,1); Lapv_y = spdiags([e -2*e e], -1:1, rv-2, rv-2)/(dy*dy);
Ix = speye(cv); Iy = speye(rv-2);
Lapv = kron(Lapv_x,Iy) + kron(Ix,Lapv_y);
Lapv = mu.*Lapv;


%Build Differentiaion matrix for the pressure in the X and Y directions
%Build Differentiaion matrix for the pressure in the X and Y directions
[rp,cp] = size(p);
G_x = zeros(cp-1,cp);
G_x(1:cp:end) = -1;
G_x(cp:cp:end) = 1;
G_x = (1/dx).*G_x; %Scale appropriately
Ix = speye(rp,rp);
G_x = kron(G_x,Ix);


%Now for the y-derivative
G_y = zeros(rp-1,rp);
G_y(1:rp:end) = -1;
G_y(rp:rp:end) = 1;
G_y = (1/dy).*G_y; %Scale appropriately
Iy = speye(cp,cp);
G_y = kron(Iy,G_y);



%Build Differentiation matrices for the continuity equation
%Both are approximated via centered differences w/r/t the centerd nodes
[ru,cu] = size(u);
Dx = zeros(cu-1,cu-2);
Dx(1:cu:end) = 1;
Dx(2:cu:end) = -1;
Ix = speye(rp);
Dx = kron(Dx,Ix);
Dx = (1/dx).*Dx;

[rv,cv] = size(v);
Dy = zeros(rv-1,rv-2);
Dy(1:rv:end) = 1;
Dy(2:rv:end) = -1;
Iy = speye(cp);
Dy = kron(Iy,Dy);
Dy = (1/dy).*Dy;

%Start assemblying the linear system to be used
A1 = (rho/dt).*speye(size(Lapu));
A2 = (rho/dt).*speye(size(Lapv));
N1 = zeros(size(A1));
N2 = zeros(size(A2));
N3 = zeros(length(p(:,1))*length(p(1,:)));
N7 = zeros(length(p(:,1))*length(p(1,:)),1);



%Vectors for RHS to implement the inflow and outflow conditions as well as
%the dirichlet conditions at the wall.
bcu_in = zeros(size(u(:,2:end-1)));
bcu_in = bcu_in(:);
bcu_in(1:ru) = mu.*u_exact(x(2),y_c,0);
bcu_out = zeros(size(u(:,2:end-1)));
bcu_out = bcu_out(:);
bcu_out(ru*(cu-2)+1 - ru:end) = mu.*u_exact(x(end),y_c,0);


bcv_bot = zeros(size(v(2:end-1,:)));
bcv_bot = bcv_bot(:);
bcv_bot(1:rv-2:end) = mu.*v_exact(x_c,y(1),0);

bcv_top = zeros(size(v(2:end-1,:)));
bcv_top = bcv_top(:);
bcv_top(rv-2:rv-2:end) = mu.*v_exact(x_c,y(1),0);

Dxu_bclf = zeros(size(p));
Dxu_bclf = Dxu_bclf(:);
Dxu_bclf(1:ru) = u_exact(x(2), y_c);

Dxu_bcrt = zeros(size(p));
Dxu_bcrt = Dxu_bcrt(:);
Dxu_bcrt((cu-1)*ru+1 - ru:end) = u_exact(x(end),y_c);

Dxv_bcbot = zeros(size(p));
Dxv_bcbot = Dxv_bcbot(:);
Dxv_bcbot(1:rv-1:end) = v_exact(x_c,y(1),0);


Dxv_bctop = zeros(size(p));
Dxv_bctop = Dxv_bctop(:);
Dxv_bctop(rv-2:rv-1:end) = v_exact(x_c,y(end),0);

%Setup up the explicit convective terms
[u_g,v_g] = Ghostnodespoiseuille(u,v,dx,dy,mu);
u_g = u_g(:,2:end-1);
v_g = v_g(2:end-1,:);

uceg = 0.5.*(u_g(:,1:end-1) + u_g(:,2:end));
vceg = 0.5.*(v_g(1:end-1,:) + v_g(2:end,:));

u_corn = 0.25.*(uceg(1:end-1,1:end-1) + uceg(2:end,1:end-1) + uceg(1:end-1,2:end) + uceg(2:end,2:end));

vxside = 0.25.*(v_g(1:end-1,2:end-2) + v_g(2:end,2:end-2) + v_g(1:end-1,3:end-1) + v_g(2:end,3:end-1));

Dxu = (1/dx).*(uceg(2:end-1,2:end) - uceg(2:end-1,1:end-1));
Dxuu= Dxu.*u_g(2:end-1,2:end-1);

Dyu = (1/dy).*(u_corn(2:end,:) - u_corn(1:end-1,:));
Dyuv = Dyu.*vxside;

Conv_u = Dxuu + Dyuv;
Conv_u_rhs = -1.0*rho.*Conv_u;


%Now for the terms for v
v_corn = 0.25.*(vceg(1:end-1,1:end-1) + vceg(2:end,1:end-1) + vceg(1:end-1,2:end) + vceg(2:end,2:end));
uyside = 0.25.*(u_g(2:end-2,1:end-1) + u_g(2:end-2,2:end) + u_g(3:end-1,1:end-1) + u_g(3:end-1,2:end));

Dxv = (1/dx).*(v_corn(:,2:end) - v_corn(:,1:end-1));
Dxvu = Dxv.*uyside;

Dyv = (1/dy).*(vceg(2:end,2:end-1) - vceg(1:end-1,2:end-1));
Dyvv = Dyv.*v_g(2:end-1,2:end-1);
Conv_v = Dxvu + Dyvv;
Conv_v_rhs = -1.0*rho.*Conv_v;

%Construct the submatrices which represent the spreading and interpolation
%operators ... 

%Compute the relevant weights and matrix entries needed before hand
%assuming the structure is fixed
S_x = zeros(N*(N-1),2*NIB);
S_y = zeros(N*(N-1),2*NIB);
S_x1 = zeros(2*NIB,N*(N-1));
S_y1 = zeros(2*NIB,N*(N-1));
N4 = zeros((N*N),2*NIB);
N5 = zeros(2*NIB,N*N);
N6 = zeros(2*NIB,2*NIB);

dX = mfac*dx;
C = dX/(dx*dy);
%Initialize data to store IB kernel params
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


%Setup the linear system to solve
LHS = [A1 - Lapu, N1, G_x, S_x; N2, A2 - Lapv, G_y, S_y; -Dx, -Dy,N3,N4;S_x1,S_y1,N5,N6];
%Now fill the relevant entries of this matrix with the interpolation and
%spreading operators....

for k = 1:NIB
    sxx = X(k,1)/dx;
    sxy = X(k,2)/dy + 1/2; %subtract off 1/2 because of the mac grid
    sx = [sxx,sxy];
    syy = X(k,2)/dy;
    syx = X(k,1)/dx + 1/2;
    sy = [syx,syy];
    Ixx(k,:) = floor(sx); %Revelant indices for the forcing on the eulerian grid
    Iyy(k,:) = floor(sy);
    j1x(k,:) = (Ixx(k,1)-1):(Ixx(k,1)+2); %These are the indices on the eulerian grid which support the discrete delta function
    i1x(k,:) = (Ixx(k,2)-1):(Ixx(k,2)+2);
    j1y(k,:) = (Iyy(k,1)-1):(Iyy(k,1)+2);
    i1y(k,:) = (Iyy(k,2)-1):(Iyy(k,2)+2);
    rx(k,:) = (-Ixx(k,:)+sx);
    ry(k,:) = (-Iyy(k,:)+sy);
    wsx(:,:,k) = IB4_1(rx(k,1)).*IB4_2(rx(k,2));
    wsy(:,:,k) = IB4_1(ry(k,1)).*IB4_2(ry(k,2));
    for q = 1:length(j1x(k,:))
        for n = 1:length(i1x(k,:))
            %Fill in for the u terms
            LHS((j1x(k,q)-2)*N + i1x(k,n),N*(N-1) + N*(N-1) + N*N +k) ...
                = LHS((j1x(k,q)-2)*N + i1x(k,n),N*(N-1) + N*(N-1) + N*N +k) - C.*wsx(q,n,k);
            
            %Fill in for the v terms
            
            LHS(N*(N-1) + (j1y(k,q)-1)*(N-1) + (i1y(k,n)),N*(N-1) + N*(N-1) + N*N + NIB + k) ...
                = LHS(N*(N-1) + (j1y(k,q)-1)*(N-1) + (i1y(k,n)),N*(N-1) + N*(N-1) + N*N + NIB + k) - C.*wsy(q,n,k);
        end
    end
    % Now compute the interpolation operator which projects U onto the
    % immersed structure
    for q = 1:length(j1x(k,:))
        for n = 1:length(i1x(k,:))
    LHS(N*(N-1) + N*(N-1) + N*N + k,(j1x(k,q)-2)*N + i1x(k,n)) = ...
        LHS(N*(N-1) + N*(N-1) + N*N + k,(j1x(k,q)-2)*N + i1x(k,n)) + dt*kappa.*wsx(q,n,k);
    
    LHS(N*(N-1) + N*(N-1) + N*N + NIB + k, N*(N-1) + (j1y(k,q)-1)*(N-1) + (i1y(k,n)-1)) = ...
        LHS(N*(N-1) + N*(N-1) + N*N + NIB + k, N*(N-1) + (j1y(k,q)-1)*(N-1) + (i1y(k,n)-1)) + dt*kappa.*wsy(q,n,k);
    
    
        end
    end
   
    
    
    
    
end

% Now create the RHS to invert the LHS on ...

RHS = [A1*u_int(:) + 0.*Conv_u_rhs(:) +  bcu_in/(dx*dx) + bcu_out/(dx*dx);...
    A2*v_int(:) + 0.*Conv_v_rhs(:) + bcv_bot/(dy*dy) + bcv_top/(dy*dy);...
    N7(:) - Dxu_bclf/dx + Dxu_bcrt/dx - Dxv_bcbot/dy + Dxv_bctop/dy;
    0.*F1(:);
    0.*F2(:)];

%Use the black box... 
vals = LHS\RHS;
u_new = vals(1:ru*(cu-2));
u_new = reshape(u_new,ru,cu-2);
v_new = vals(ru*(cu-2)+1:ru*(cu-2)+(rv-2)*cv);
v_new = reshape(v_new,rv-2,cv);
p_new = vals(ru*(cu-2) + (rv-2)*cv + 1:ru*(cu-2) + (rv-2)*cv + N*N);
p_new = p_new - dx*dy*sum(p_new);
F1_new = vals(ru*(cu-2) + (rv-2)*cv + N*N + 1: ru*(cu-2) + (rv-2)*cv + N*N + NIB);
F2_new = vals(ru*(cu-2) + (rv-2)*cv + N*N + NIB + 1:ru*(cu-2) + (rv-2)*cv + N*N + 2*NIB);
%Normalize the pressure: the gradient is really all that matters
p_new = reshape(p_new,rp,cp);

%Also spit out the corresponding block matrices ... the naming is according
%to Benzi et als saddle point paper
A = LHS(1:2*(N-1)*N + N*N,1:2*N*(N-1) + N*N);
B1 = LHS(1:(2*(N-1)*N + N*N),2*N*(N-1)+N*N + 1:end);
B2 = LHS(2*N*(N-1) + N*N +1:end,1:(2*(N-1)*N + N*N));





%Add in the analytic solutions for the boundary
u_new = [u_exact(x(1),y_c,0)', u_new, u_exact(x(end),y_c,0)'];
v_new = [v_exact(x_c,y(1),0); v_new; v_exact(x_c,y(end),0)];
end













