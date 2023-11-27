function [ffx,ffy] = DFIB_spread_spectral(phi,F,X,V,dx,dy,ds)
%V is the volume of the domain
%First solve the poisson equation for each f on the grid. 
%Construct the RHS of the equation
NIB = length(X(:,1));
[rphi,cphi] = size(phi);
K = 59/60-sqrt(29)/20; 
RHS = zeros(rphi,cphi);
C = (ds)/(dx*dy);

fx0 = (1/V).*sum(F(:,1)).*ds;
fy0 = (1/V).*sum(F(:,2)).*ds;

for k = 1:NIB
   sxx = (X(k,1) + pi)/dx;
   syy = (X(k,2) + pi)/dy;
   Ixx = floor(sxx);
   Iyy = floor(syy);
   rx = sxx - Ixx; %compute the relevant r value defined in (0,1) in both x and y
   ry = syy - Iyy;
   
   %Get the relevant indices to interpolate onto the lagrangian structure
   %need to use modular arithmetic to account for periodic boundary
   %conditions
   
   j1 = mod(Ixx-2:Ixx+3,cphi)+1;
   i1 = mod(Iyy-2:Iyy+3,rphi)+1;
   
   %weights for x multiplier
   wsx(:,:) = -(1/dy).*IB6x(rx,K).*IB6y_prime(ry,K);
   wsy(:,:) = (1/dx).*IB6x_prime(rx,K).*IB6y(ry,K);
   
   RHS(i1,j1) = RHS(i1,j1) + C*wsx(:,:)*F(k,1) + C*wsy(:,:)*F(k,2);

   
    
end


%get ghostnodes
[RHS1g,RHS2g] = Ghostnodesside_periodic(RHS,RHS);
%Use centered differences to get the RHSs for the poisson equations
RHSx = (RHS1g(3:end,2:end-1) - RHS1g(1:end-2,2:end-1))./(2*dy);
RHSy = (RHS2g(2:end-1,3:end) - RHS2g(2:end-1,1:end-2))./(2*dx);

%Now that both RHS have been generated solve the poisson equation for the
%forcing

ffx = fftSolveLaplacian2D(RHSx(:),rphi,cphi,dx,dy);
ffy = fftSolveLaplacian2D(-RHSy(:),rphi,cphi,dx,dy);

%Add back the mean value of the forcing
ffx = ffx + fx0;
ffx = reshape(ffx,rphi,cphi);
ffy = ffy + fy0;
ffy = reshape(ffy,rphi,cphi);

end

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
