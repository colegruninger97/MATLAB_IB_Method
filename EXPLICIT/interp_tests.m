%Test for interpolation functions
N = 256;
dx = 1/N;
x = 0:dx:1-dx;
y = x;
xc = dx/2:dx:1-dx/2;
yc = xc;

f = @(x,y) sin(2*pi.*x).*sin(2*pi.*y);
u_exact = @(x,y) cos(2*pi.*x).*sin(2*pi.*y);
v_exact = @(x,y) -cos(2*pi.*y).*sin(2*pi.*x);
[Xnode,Ynode] = meshgrid(x,y);
[Xside,Ycent] = meshgrid(x,yc);
[Xcent,Yside] = meshgrid(xc,y);

fx = f(Xside,Ycent);
fy = f(Xcent,Yside);
fn = f(Xnode,Ynode);
a = 0.11;
b = 0.36;
% [U,V] = interpIB6(fx,fy,[a,b],dx,dx);
% Phi = interpIB6nodes(fn,[a,b],dx,dx);
% Error = abs(U - f(a,b)) + abs(V-f(a,b));
% Error_p = abs(Phi-f(a,b));

u_trial = u_exact(Xside,Ycent);
v_trial = v_exact(Xcent,Yside);

phi = solve_for_potential(u_trial,v_trial,dx,dx);
[phixg,phiyg] = Ghostnodesside_periodic(phi,phi);
uu = -(1/dx).*(phixg(3:end,2:end-1) - phixg(2:end-1,2:end-1));
vv = (1/dx).*(phiyg(2:end-1,3:end)-phiyg(2:end-1,2:end-1));
[U,V] = DFIB_interp(phi,0,0,[a,b],dx,dx);
%Compare these values to the original velocity
Error_u = abs(u_exact(a,b)-U);
Error_v = abs(v_exact(a,b)-V);
Error_u
Error_v

%%

mfac = 0.01;
r_cyl = 0.25;
%Specify the Lagrangian mesh spacing
ds = mfac*dx/r_cyl;
approx = round(2*pi/ds);
ds = 2*pi/approx;
dX = mfac*dx;
l = length(0:ds:2*pi-ds);
X_1 = 1.*(0.5.*ones(1,l) + (0.25)*cos(0:ds:2*pi-ds));
X_2 = 1.*(0.5.*ones(1,l) +  (0.25)*sin(0:ds:2*pi-ds));
X = [X_1;X_2];
X = X';
% F = Elastic_Laplacian(X,1,dx);
F = uniform_normal_force_circle(X,1,ds,0.25);
[ffx,ffy] = DFIB_spread(zeros(N,N),F,X,1,dx,dx,ds);
[ffx_bspline,ffy_bspline] = spreadIB4(F,X,zeros(N,N),zeros(N,N),dx,dx,ds);
[ffxg,ffyg] = Ghostnodesside_periodic(ffx,ffy);
Div = (ffxg(2:end-1,3:end)-ffxg(2:end-1,2:end-1))./dx + (ffyg(3:end,2:end-1)-ffyg(2:end-1,2:end-1))./dx;


function x = fftSolveLaplacian2D(b, nx, ny, dx, dy)
    % Transform the right-hand side to the frequency domain
    b_fft = fft2(reshape(b,ny,nx));

    % Construct the wave numbers in both dimensions
    kx = (2 * pi / (nx * dx)) * [0:nx/2, -nx/2+1:-1];
    ky = (2 * pi / (ny * dy)) * [0:ny/2, -ny/2+1:-1];

    % Construct a grid of wave numbers
    [KX, KY] = meshgrid(kx, ky);

    % Eigenvalues of the periodic Laplacian in the frequency domain
    laplacian_eigenvalues = -(KX.^2 + KY.^2);

    % Solve in the frequency domain
    x_fft = b_fft ./ laplacian_eigenvalues;

    % Handle the singularity at k = 0 (DC component)
    x_fft(1,1) = 0;

    % Transform the solution back to the spatial domain
    x = ifft2(x_fft, 'symmetric'); % 'symmetric' ensures real output
    x = x(:);
end