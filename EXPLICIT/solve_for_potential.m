function phi = solve_for_potential(u,v,dx,dy)
[ru,cu] = size(u);
%compute the discrete curl of the velocity field
[ug,vg] = Ghostnodesside_periodic(u,v);

dvdx = (vg(2:end-1,2:end-1) - vg(2:end-1,1:end-2)) ./ dx;
dudy = (ug(2:end-1,2:end-1) - ug(1:end-2,2:end-1)) ./ dy;

curl = dvdx - dudy; %Second order accurate quantity defined at the nodes of the grid

phi = fftPreconditionedSolveLaplacian2D(curl(:),cu,ru,dx,dy);
phi = reshape(phi,ru,cu);



end

function x = fftPreconditionedSolveLaplacian2D(b, nx, ny, dx, dy)
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

