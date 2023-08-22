%Script used to construct the stokes operator, spreading, and interpolation
%operators for flow past a cylinder in a periodic domain

%Setup the cartesian grid (needed to shrink the domain to use larger N)
N = 16;
dx = 1/N;
dy = 1/N;
x = 0:dx:8-dx; %Don't include the end points since the domain is periodic
y = 0:dy:4-dx;
xc = 0.5*dx:dx:8-0.5*dx;
yc = 0.5*dy:dy:4-0.5*dy;

%Setup dummy velocity and pressure values
u = zeros(length(yc),length(x));
v = zeros(length(y),length(xc));
p = zeros(length(yc),length(xc));

%Setup up some test parameters
%viscosity
mu = 1.0;
%density
rho = 1.0;
%timestep size
dt = 0.1;
%spring constant
kappa = 100;




%Setup the Cylinder Geometry
r_cyl = 1.0;
m_fac = 1.0;
dX = m_fac*dx;
N_IB = round(2*pi/dX);
ds = 2*pi/N_IB;
s = 0:ds:2*pi-ds;
X_1 = 4.0 + r_cyl.*cos(s);
X_2 = 2.0 + r_cyl.*sin(s); 


%Build the extended Saddle point system 
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
G_x = zeros(cp,cp);
G_x(1:cp+1:end) = 1/dx;
G_x(2:cp+1:end) = -1/dx;
G_x(1,end) = -1/dx;
G_x = kron(G_x,eye(rp));
 
G_y = zeros(rp,rp);
G_y(1:rp+1:end) = 1/dy;
G_y(2:rp+1:end) = -1/dy;
G_y(1,end) = -1/dy;
G_y = kron(eye(cp),G_y);
 
%Build the differentiation matrices for the incompressibility constraint
Dx = zeros(cu,cu);
Dx(1:cu+1:end) = -1/dx;
Dx(cu+1:cu+1:end) = 1/dx;
Dx(end,1) = 1/dx;
Dx = kron(Dx,eye(ru));
 
Dy = zeros(rv,rv);
Dy(1:rv+1:end) = -1/dy;
Dy(rv+1:rv+1:end) = 1/dy;
Dy(end,1) = 1/dy;
Dy = kron(eye(cv),Dy);

A1 = (rho/dt).*speye(size(Lapu));
A2 = (rho/dt).*speye(size(Lapv));


%Initialize interp and spread matrices...
S_x = zeros(length(x)*length(yc),2*N_IB);
S_y = zeros(length(y)*length(xc),2*N_IB);
S_x1 = S_x';
S_y1 = S_y';
N3 = zeros(length(p(:)),length(p(:)));
N4 = zeros(length(xc)*length(yc),2*N_IB);
N5 = zeros(2*N_IB,length(xc)*length(yc));
N6 = zeros(2*N_IB,2*N_IB);

C = dX/(dx*dy);

Extended_Saddle = [A1 - Lapu, zeros(length(A1(:,1)),length(A2(1,:))), G_x, S_x; zeros(length(A2(:,1)),length(A1(1,:))), A2 - Lapv, G_y, S_y; -Dx, -Dy,N3,N4;S_x1,S_y1,N5,N6];

% %Fill in for the interpolation and spreading operators
for k = 1:N_IB
    sxx = X_1(k)/dx; %This part of the code assumes that the lefthand corner of the domain is x = 0 and y = 0
    sxy = X_2(k)/dy - 1/2;
    sx = [sxx,sxy];
    syy = X_2(k)/dy;
    syx = X_1(k)/dx - 1/2;
    sy = [syx,syy];
    Ixx(k,:) = ceil(sx); %Revelant indices for the forcing on the eulerian grid
    Iyy(k,:) = ceil(sy);
    j1x(k,:) = (Ixx(k,1)-1):(Ixx(k,1)+2); %These are the indices on the eulerian grid which support the discrete delta function
    i1x(k,:) = (Ixx(k,2)-1):(Ixx(k,2)+2);
    j1y(k,:) = (Iyy(k,1)-1):(Iyy(k,1)+2);
    i1y(k,:) = (Iyy(k,2)-1):(Iyy(k,2)+2);
    rx(k,:) = (Ixx(k,:)-sx);
    ry(k,:) = (Iyy(k,:)-sy);
    wsx(:,:,k) = IB4_1(rx(k,1)).*IB4_2(rx(k,2)); %This code is particular to the standard IB4 kernel, the above code will need to change for different kernels
    wsy(:,:,k) = IB4_1(ry(k,1)).*IB4_2(ry(k,2));
    %Fill in the spreading operator
    for q = 1:length(j1x(k,:))
        for n = 1:length(i1x(k,:))
            %Fill in for the u terms
            Extended_Saddle((j1x(k,q)-1)*(length(yc)) + i1x(k,n),length(yc)*length(x) + length(xc)*length(y) + length(yc)*length(xc) +k) ...
                = Extended_Saddle((j1x(k,q)-1)*(length(yc)) + i1x(k,n),length(yc)*length(x) + length(xc)*length(y) + length(yc)*length(xc) +k) - C.*wsx(q,n,k);
            
            %Fill in for the v terms
            
            Extended_Saddle(length(yc)*length(x) + (j1y(k,q)-1)*length(y) + (i1y(k,n)),length(yc)*length(x) + length(xc)*length(y) + length(yc)*length(xc) + N_IB + k) ...
                = Extended_Saddle(length(yc)*length(x) + (j1y(k,q)-1)*length(y) + (i1y(k,n)),length(yc)*length(x) + length(xc)*length(y) + length(yc)*length(xc) + N_IB + k) - C.*wsy(q,n,k);
        end
    end
    % Now compute the interpolation operator which projects U onto the
    % immersed structure
    for q = 1:length(j1x(k,:))
        for n = 1:length(i1x(k,:))
    Extended_Saddle(length(yc)*length(x) + length(xc)*length(y) + length(yc)*length(xc) + k,(j1x(k,q)-1)*length(yc) + i1x(k,n)) = ...
        Extended_Saddle(length(yc)*length(x) + length(xc)*length(y) + length(yc)*length(xc) + k,(j1x(k,q)-1)*length(yc) + i1x(k,n)) + dt*kappa.*wsx(q,n,k);
    
    Extended_Saddle(length(yc)*length(x) + length(xc)*length(y) + length(yc)*length(xc) + N_IB + k, length(yc)*length(x) + (j1y(k,q)-1)*length(y) + (i1y(k,n))) = ...
        Extended_Saddle(length(yc)*length(x) + length(xc)*length(y) + length(yc)*length(xc)+ N_IB + k, length(yc)*length(x) + (j1y(k,q)-1)*length(y) + (i1y(k,n))) + dt*kappa.*wsy(q,n,k);
    
    
        end
    end
   
    
end

%Extract the stokes operator, the interpolation operator, and the spread
%operator from the extended saddle point matrix 

Stokes_Op = Extended_Saddle(1:length(yc)*length(x) + length(xc)*length(y) + length(yc)*length(xc),1:length(yc)*length(x) + length(xc)*length(y) + length(yc)*length(xc));
Spread_Op = Extended_Saddle(1:length(yc)*length(x)+length(xc)*length(y),length(yc)*length(x) + length(xc)*length(y) + length(yc)*length(xc)+1:end);
Interp_Op = Extended_Saddle(length(yc)*length(x) + length(xc)*length(y) + length(yc)*length(xc)+1:end,1:length(yc)*length(x) + length(y)*length(xc));
S_star_S = Interp_Op*Spread_Op;


function w = IB4_1(r)
%Implementation of the standard IB-4 kernel as derived in Peskin's original
%work
w = zeros(1,4);
q=sqrt(1+4*r*(1-r));
w(1,4)=(1+2*r-q)/8;
w(1,3)=(1+2*r+q)/8;
w(1,2)=(3-2*r+q)/8;
w(1,1)=(3-2*r-q)/8;

end


function w = IB4_2(r)
%Implementation of the standard IB4 kernel as derived by peskin
w = zeros(4,1);
q=sqrt(1+4*r*(1-r));
w(4,1)=(1+2*r-q)/8;
w(3,1)=(1+2*r+q)/8;
w(2,1)=(3-2*r+q)/8;
w(1,1)=(3-2*r-q)/8;

end








