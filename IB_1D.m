%IB_1D 

%The following is a 1-dimensional implementation of the immersed boundary
%method used to test the frictional penalty parameter scaling derived

%rho(u_t + u_x) = /mu u_xx + friction penalty term
%X_t = U
C = 6;
%We will simulate the dynamics on a periodic grid

N = 128;
h = 1/N;
x = 0:h:1-h;
dt = 0.5*h*h;
rho = 1; %density
eta = (rho*h/dt)*C; %eta < rho*h/dt
mu = 1;


%Setup the initial conditiond
u = 0.*sin(1*pi*x);

%Get initial location for the IB point 
X = rand;
%Add in ghost cells
u_g = zeros(1,length(u)+2);
u_g(2:end-1) = u;
u_g(1) = u(end);
u_g(end) = u(1);

%Create a vector of the penalty forcing terms 
phi = IB4(mod((x-X)./h,N));
for i = 1:10000
U = One_d_interp(X,u,h);
u = u - dt.*((u_g(3:end) - u_g(1:end-2))./(2*h)) + (mu.*dt/(rho*h*h)).*(u_g(3:end) + u_g(1:end-2) - 2.*u_g(2:end-1)) - (eta.*dt/(h*rho)).*phi.*U + phi.*(eta.*dt)/(rho*h)*1;
u_g(2:end-1) = u;
u_g(1) = u(end);
u_g(end) = u(1);
X = X + dt*U;
phi = IB4(mod((x-X)./h,N));
plot(x,u)
axis([0 1 0 3])
pause(0)
end









function U = One_d_interp(X,u,h)
cu = length(u);
sxx = X/h;
Ixx = floor(sxx);
rx = sxx - Ixx;
j1x = mod(Ixx - 1:Ixx + 2,cu)+1;
ws = IB4_1(rx);
U = sum(u(j1x).*ws(:)','all');
end

function val=IB4(r)
val = zeros(1,length(r));
for i = 1:length(r)
val(i) = (1/8)*(...
    (abs(r(i))<1).*(3-2*abs(r(i))+sqrt(1+4*abs(r(i))-4*abs(r(i)).^2)) +...
    (1<=abs(r(i)) & abs(r(i))<2).*(5-2*abs(r(i))-sqrt(-7+12*abs(r(i))-4*abs(r(i)).^2)));

end
end





