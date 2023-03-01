% 1-D CIB solver to use to test out ill-conditioning theories...
clear;
clc;
alp = 3.0;
N = 128;
NIB = alp*N;

ds = 1/NIB;
h = 1/N;

% The goal is to solve the set of equations using the FFT approach in order
% to use the analysis given by Barret....

x_s = 0:ds:1-ds;
x_h = 0:h:1-h;

%Define the locations of the Lagrange multipliers
X0 = 0.5;

% The goal is to solve u''(x) + sigma*u = f(x) + lambda0*d_h(x-X0) +
% lambda1*d_h(x-X1)
sigma = 1.0;
u_hat = zeros(1,length(x_h));
f = @(x) sin(2*pi*x);
f = f(x_h);
phi_h = zeros(1,length(x_h));
Ix = X0/h;
Ix = floor(Ix);
rx = X0/h - Ix;
phi_h(Ix:Ix+3) = BS_3x(rx);
%Now try to compute the FFT of phi_h
Phi_xi = fft(phi_h);
% plot(abs(Phi_xi))

   
f_hat = fft(f);
for i = 1:length(x_h)
u_hat(i) = 0.*(h*h).*f_hat(i)/(2.*cos((2*pi*(i-1))/N) - 2 + h*h*sigma) + (h*h).*Phi_xi(i)./(2.*cos((2*pi*(i-1))/N) - 2 + h*h*sigma);
end
u = ifft(u_hat);
u = real(u);
u_true = (1/(1-4*pi*pi)).*sin(2*pi.*x_h);

plot(u)






