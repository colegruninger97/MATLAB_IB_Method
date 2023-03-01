function [u_g,v_g] = Ghostnodespipeflow(u,v,dx,dy,mu)
Lx = 1.0;
Ly = 1.0;
Px = -0.5;
x_side = -dx:dx:Lx+dx; y_side = -dy:dy:Ly+dy;
x_cent = 0.5*(x_side(2:end)+x_side(1:end-1));
y_cent = 0.5*(y_side(2:end)+y_side(1:end-1));


%For pipeflow, need to utilize the exact solutions to employ the bcs
%employ the ghost values for the components of the stress tensor
u_exact = @(x,y,t) 0.*x + (Px/(2*(mu))).*y.*(y-Ly);
v_exact = @(x,y,t) 0.*x.*y;

[r,c] = size(u);
u_g = zeros(r+2,c+2);
u_g(2:end-1,2:end-1) = u;
%inflow ghosts
u_g(2:end-1,1) = -u_g(2:end-1,3) + 2*u_g(2:end-1,2);
%outflow ghosts
u_g(2:end-1,end) = -u_g(2:end-1,end-2) + 2*u_g(2:end-1,end-1);
%Ghost values above and below the wall (use quad extrapolation to enforce
%dirichlet condition)
u_g(1,2:end-1) = -2*u_g(2,2:end-1) + (1/3)*u_g(3,2:end-1);
u_g(end,2:end-1) = -2*u_g(end-1,2:end-1) + (1/3)*u_g(end-2,2:end-1);

[r,c] = size(v);
% v_g = zeros(r+2,c+2);
% v_g(2:end-1,2:end-1) = v;
% %inflow ghosts
% v_g(2:end-1,1) = -v_g(2:end-1,2);
% %outflow ghosts
% v_g(2:end-1,end) = -v_g(2:end-1,end-1);
% %Ghost values above and below the wall (all zero for now...)
% v_g(1,:) = -v_g(3,:);
% v_g(end,:) = -v_g(end-2,:);

v_g = zeros(r+2,c+2);
v_g(2:end-1,2:end-1) = v;
%inflow ---> analytic solution
v_g(2:end-1,1) = -v_g(2:end-1,2); %Dirichlet Boundary Condition
%same for outflow
v_g(2:end-1,end) = -v_g(2:end-1,end-1); %Dirichlet Boundary Condition
%top and bottom of the walls use quad extrap
v_g(1,2:end-1) = -2*v_g(2,2:end-1) + (1/3)*v_g(3,2:end-1);
v_g(end,2:end-1) = -2*v_g(end-1,2:end-1) + (1/3)*v_g(3,2:end-1);
end
