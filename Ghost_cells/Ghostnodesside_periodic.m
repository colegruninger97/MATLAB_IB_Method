function [u_new, v_new] = Ghostnodesside_periodic(u, v)
% Fills ghost cells using periodic boundary conditions
% for side centered vector field (u,v)

[ru,cu] = size(u);
[rv,cv] = size(v);
u_new = zeros(ru+2,cu+2);
v_new = zeros(rv+2,cv+2);
u_new(2:end-1,2:end-1) = u;
v_new(2:end-1,2:end-1) = v;

%Fill em

u_new(2:end-1,end) = u(:,1);
u_new(2:end-1,1) = u(:,end);
u_new(1,2:end-1) = u(end,:);
u_new(end,2:end-1) = u(1,:);

v_new(2:end-1,end) = v(:,1);
v_new(2:end-1,1) = v(:,end);
v_new(1,2:end-1) = v(end,:);
v_new(end,2:end-1) = v(1,:);

u_new(1,1) = u_new(1,end-1);
u_new(end,1) = u_new(2,1);
u_new(end,end) = u_new(2,end);
u_new(1,end) = u_new(end-1,end);

v_new(1,1) = v_new(end-1,1);
v_new(end,end) = v_new(2,end);
v_new(end,1) = v_new(2,1);
v_new(1,end) = v_new(end-1,end);
%fill the corners needed for convective terms
% u_new(1,1) = (1/3).*(u_new(2,1) + u(1,1) + u_new(1,2));
% u_new(end,1) = (1/3).*(u_new(end,2) + u(end,end)+u_new(end-1,1));
% u_new(end,end) = (1/3).*(u_new(end-1,end)+u(end,end)+u_new(end,end-1));
% u_new(1,end) = (1/3).*(u_new(2,end) + u_new(1,end-1) + u(1,end));
% 
% v_new(1,1) = (1/3).*(v_new(2,1) + v(1,1) + v_new(1,2));
% v_new(end,1) = (1/3).*(v_new(end,2) + v(end,end)+v_new(end-1,1));
% v_new(end,end) = (1/3).*(v_new(end-1,end)+v(end,end)+v_new(end,end-1));
% v_new(1,end) = (1/3).*(v_new(2,end) + v_new(1,end-1) + v(1,end));
end
