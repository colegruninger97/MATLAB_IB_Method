%Error tests in 2D. Here we use the analytic solution to stokes flow past a
%cylinder to test the accuracy of the IB kernels when spreading and
%interpolating 

%Construct the analytic, free space, solution to stokes flow 
h = 1/256;
x = -2:h:2;
y = -2:h:2;
K =  59/60-sqrt(29)/20;
a = 1;
%K = -0.815559903284;
u_stokes_x_deriv = @(x,y) 3*a.*x.*(-2*a.*x.^2 + 2.*x.^4 + 3*a.*y.^2 + (x.^2).*y.^2 - y.^4)./(4.*(x.^2 + y.^2).^(7/2));
u_stokes_y_deriv = @(x,y) 3*a.*y.*(-4*a.*x.^2 + 4.*x.^4 + a.*y.^2 + 5.*(x.^2).*y.^2 - y.^4)./(4.*(x.^2 + y.^2).^(7/2));
[xx,yy] = meshgrid(x,y);
a = 1.0; %radius of the cylinder
mfac = 2.0;
N_IB = round(2*pi/(mfac*h));
ds = 2*pi/N_IB;
mu = 0.1; %viscosity
U = 1; %velocity in the farfield
x_cyl = a.*cos(0:ds:2*pi-ds);
y_cyl = a.*sin(0:ds:2*pi-ds);
ux_jump = u_stokes_x_deriv(x_cyl,y_cyl);
uy_jump = u_stokes_y_deriv(x_cyl,y_cyl);
jump_normal = ux_jump.*cos(0:ds:2*pi-ds) + uy_jump.*sin(0:ds:2*pi-ds);
gamma_x = @(x) sqrt(1-x.^2);
delta_BS3_2d = @(x,y) (1/(h^2)).*BS3(x./h).*BS3(y./h)';
delta_BS2_2d = @(x,y) (1/(h^2)).*BS2(x./h).*BS2(y./h)';
delta_IB6_2d = @(x,y) (1/(h^2)).*IB6(x./h,K).*IB6(y./h,K)';
% delta_IB6_composed = @(x,y) (1/h).*IB6((x)./h,K).*IB6_kernel_composed(x,y,K,h);
% %u_interp_conform = compute_2d_integral(Stokes_x_vel(x,y,1,1),delta_IB6_composed(x,y),x,y)
u_interp_BS3 = zeros(1,length(x_cyl));
u_interp_IB6 = zeros(1,length(x_cyl));
u_interp_BS2 = zeros(1,length(x_cyl));

% u_interp_BS3 = compute_2d_integral(Stokes_x_vel(x,y,1,1),delta_BS3_2d(x - sqrt(2)/2,y - sqrt(2)/2),x,y);%[0.0029 0.0013 5.9476e-4 3.0681e-4 1.6191e-4 8.745e-5 4.5109e-5]
% u_interp_BS2 = compute_2d_integral(Stokes_x_vel(x,y,1,1),delta_BS2_2d(x - sqrt(2)/2,y - sqrt(2)/2),x,y);%[0.0026 0.0011 4.7637e-4 2.4945e-4 1.3511e-4 7.5763e-5 4.0346e-5]
% u_interp_IB6 = compute_2d_integral(Stokes_x_vel(x,y,1,1),delta_IB6_2d(x - sqrt(2)/2,y - sqrt(2)/2),x,y);%[0.0041 0.0019 9.4913e-4 4.8064e-4 2.4556e-4 1.265e-4 6.3772e-5]

Error_est_BS3 = integral(@(x)BS3(x).*x,0,10);
Error_est_BS2 = integral(@(x)BS2(x).*x,0,10);
Error_est_IB6 = integral(@(x)IB6(x,K).*x,0,10);


% vec = cos(5*pi/6):0.1:cos(pi/6);
% u_interp_LB = zeros(1,length(cos(5*pi/6):0.1:cos(pi/6)));
% for i = 1:length(u_interp_LB)
%     u_interp_LB(i) = compute_L_B_composed_interp(x,vec(i),y,h,uu);
% end
% L1_error = 0.1*sum(abs(u_interp_LB))
% max(abs(u_interp_LB))
% %surf(x,y,delta_IB6_composed(x,y))
% %contourf(x,y,IB6((y-sqrt(2)/2)./h,K)'.*IB6_kernel_composed(x,y,K,h))
% %Attempt at a second order kernel function
% %surf(x,y,IB6(x,K).*IB6_kernel_composed(x,y,K,h))
% %Sample force at the boundary
% compute_IB6_composed_interpolation(x,y,h,K,uu)
% 
% %contourf(x,y,IB6_y)
for i = 1:length(x_cyl)
    u_interp_BS3(i) = compute_2d_integral(Stokes_x_vel(x,y,1,1),delta_BS3_2d(x - x_cyl(i),y - y_cyl(i)),x,y);
    u_interp_IB6(i) = compute_2d_integral(Stokes_x_vel(x,y,1,1),delta_IB6_2d(x - x_cyl(i),y - y_cyl(i)),x,y);
    u_interp_BS2(i) = compute_2d_integral(Stokes_x_vel(x,y,1,1),delta_BS2_2d(x - x_cyl(i),y - y_cyl(i)),x,y);
end
%Need to set the velocity equal to zero inside the structure...
%scatter(x_cyl,y_cyl,[],abs(u_interp_BS3))
uu = Stokes_x_vel(x,y,1,1);
%plot(y,uu(:,128))



function u_stokes = Stokes_x_vel(x,y,a,mu)
[xx,yy] = meshgrid(x,y);
R = sqrt(xx.^2 + yy.^2);
u_stokes = zeros(length(y),length(x));
for i = 1:length(x)
    for j = 1:length(y)
        if R(j,i) < 1
            u_stokes(j,i) = 0;
        else
            u_stokes(j,i) = 1 - (6*pi*mu*a/(8*pi*mu)).*((xx(j,i).*xx(j,i))./(R(j,i).^3) +1./R(j,i)) - 0.25*a*a.*(1./(R(j,i).^3) - 3.*(xx(j,i).*xx(j,i))./(R(j,i).^5));
        end
    end
end

end




function force_x = Stokes_x_force(x,y,a,mu)
[xx,yy] = meshgrid(x,y);
R = sqrt(xx.^2 + yy.^2);
force_x = zeros(length(y),length(x));
for i = 1:length(x)
    for j = 1:length(y)
    force_x(j,i) = 3*mu*a.*(-4.*a.*xx(j,i).^2 + 4.*xx(j,i).^4 + a.*yy(j,i).^2 + 5.*(xx(j,i).^2).*yy(j,i).^2 + yy(j,i).^4)./(4.*(xx(j,i).^2 + yy(j,i).^2).^(7/2));

    end
end
end


function v_stokes = Stokes_y_vel(xx,yy,R,a,mu)
v_stokes = zeros(length(yy),length(xx));
for i = 1:length(xx)
    for j = 1:length(yy)
        if R(i,j) < 1
            v_stokes(i,j) = 0;
        else
            v_stokes(i,j) = - 6*pi*mu*a/(8*pi*mu).*((xx(i,j).*yy(i,j))./(R(i,j).^3)) + 0.25*a*a.*(3.*xx(i,j).*yy(i,j)./(R(i,j).^5));
        end
    end
end

end

function p_stokes = Stokes_pressure(xx,yy,R,a,mu)
p_stokes = zeros(length(yy),length(xx));
for i = 1:length(yy)
    for j = 1:length(xx)
        if R(i,j) < 1
            p_stokes(i,j) = 0;
        else
            p_stokes(i,j) = -1.5*mu*a*xx(i,j)/(R(i,j)^3);
        end
    end
end
end

function Force_x = Stokes_force_x_dir(xx,yy,R,a,mu)
Force_x = zeros(length(yy),length(xx));

for i = 1:length(yy)
    for j = 1:length(xx)
        if R(i,j) < 1
            Force_x(i,j) = 0;
        else
            Force_x(i,j) = 3*a*mu*(-4*a*xx(i,j)*xx(i,j) + 4*xx(i,j)^4 + a*yy(i,j)^2 + 5*(xx(i,j)^2)*(yy(i,j)^2) + yy(i,j)^4)/(4*(xx(i,j)^2 + yy(i,j)^2)^(7/2));
        end
    end
end


end

function Force_y = Stokes_force_y_dir(xx,yy,R,a,mu)
Force_y = zeros(length(yy),length(xx));

for i = 1:length(yy)
    for j = 1:length(xx)
        if R(i,j) < 1
            Force_y(i,j) = 0;
        else
            Force_y(i,j) = (3*a*xx(i,j)*yy(i,j)*mu*(-5*a + 3*(xx(i,j)^2 + yy(i,j)^2)))/(4*(xx(i,j)^2 + yy(i,j)^2)^(7/2));
        end
    end
end

end

function val = cos_kernel(r)
    %Width of the kernel is 4
    val = zeros(1,length(r));
    
    for i = 1:length(r)
        if r(i)>=2 || r(i) <= -2
            val(i) = 0;
        else
            val(i) = 0.25*(1 + cos(pi*r(i)/2));
        end
    end
    
   

end


function val = piecewise_linear_kernel(r)

    val = zeros(1,length(r));
    for i = 1:length(r)
        if r(i)>=1 || r(i)<= -1
        val(i) = 0;
        elseif r(i)<1 && r(i)>=0
        val(i) = 1-r(i);
        elseif r(i)<0 && r(i)> -1
        val(i) = r(i)+1;
        end
    end
end

function val = sinc(r)
    val = zeros(1,length(r));
    for i = 1:length(r)
        val(i) = sin(r(i))/r(i);
    end
end

function val = BS3(r)
    val = zeros(1,length(r));
    for i = 1:length(r)
        if abs(r(i)) > 2
            val(i) = 0;
        elseif abs(r(i)) <= 2 && abs(r(i))> 1
            val(i) = (1/6)*(2-abs(r(i)))^3;
        elseif abs(r(i)) <= 1
            val(i) = 2/3 - r(i)^2 + 0.5*abs(r(i))^3;
        end
    end
end

function val = BS2(r)
    val = zeros(1,length(r));
    for i = 1:length(r)
        if abs(r(i)) > 3/2
            val(i) = 0;
        elseif abs(r(i)) > 1/2 && abs(r(i))<= 3/2
            val(i) = (1/8)*(2*abs(r(i))-3)^2;
        elseif abs(r(i)) <= 1/2
            val(i) = 3/4 - r(i)^2;
        end
    end
            
end

function val = delta_5_smooth(x)

KK = (38 - sqrt(69))/60; 
% KK = (38 + sqrt(69))/60; 
val = zeros(1,length(x));

phi = @(r) (136 - 40*KK - 40*r.^2 + sqrt(2)*sqrt(3123 - 6840*KK + 3600*KK.^2 - 12440*r.^2 + 25680*KK*r.^2 - 12600*KK.^2*r.^2 + 8080*r.^4 - 8400*KK*r.^4 - 1400*r.^6))/280; 
  
for i = 1:length(x)
    if abs(x(i)) < 0.5
        r = abs(x(i)); 
        val(i) = phi(r); 
    elseif abs(x(i)) < 1.5
        r = abs(x(i)) - 1; 
        val(i) = (4 - 4*phi(r) - KK - 4*r + 3*KK*r - r.^2 + r.^3)/6; 
    elseif abs(x(i)) < 2.5
        r = abs(x(i)) - 2; 
        val(i) = (-2 + 2*phi(r) + 2*KK + r - 3*KK*r + 2*r.^2 - r.^3)/12; 
    else 
        val(i) = 0; 
    end 
end
end

function val=IB3(r)

val = zeros(1,length(r));
for i = 1:length(r)
val(i) =...
    (abs(r(i))<1/2).*((1/3)*(1+sqrt(1-3*abs(r(i)).^2))) +...
    (1/2<=abs(r(i)) & abs(r(i))<3/2).*((1/6)*(5-3*abs(r(i))-sqrt(-3*(1-abs(r(i))).^2+1)));
end

end

function val=IB4(r)
val = zeros(1,length(r));
for i = 1:length(r)
val(i) = (1/8)*(...
    (abs(r(i))<1).*(3-2*abs(r(i))+sqrt(1+4*abs(r(i))-4*abs(r(i)).^2)) +...
    (1<=abs(r(i)) & abs(r(i))<2).*(5-2*abs(r(i))-sqrt(-7+12*abs(r(i))-4*abs(r(i)).^2)));

end
end

function val=IB6(r,K)
val = zeros(1,length(r));
for i = 1:length(r)
    R = r(i);

% compute r in [0,1]
r(i) = (-3<r(i) & r(i)<=-2).*(r(i)+3) + (-2<r(i) & r(i)<=-1).*(r(i)+2) + ...
    (-1<r(i) & r(i)<=0) .*(r(i)+1) + (0<r(i) & r(i)<=1)  .*(r(i)+0) + ...
    (1<r(i) & r(i)<=2)  .*(r(i)-1) + (2<r(i) & r(i)<=3)  .*(r(i)-2);
alpha=28;
beta=(9/4)-(3/2)*(K+r(i).^2)+((22/3)-7*K)*r(i)-(7/3)*r(i).^3;
gamma=(1/4)*( ((161/36)-(59/6)*K+5*K^2)*(1/2)*r(i).^2 + (-(109/24)+5*K)*(1/3)*r(i).^4 + (5/18)*r(i).^6 );

discr=beta.^2-4*alpha*gamma;
% flag=0;
% if(any(discr<0))
%    flag=1
% end

pm3=(-beta+sign((3/2)-K)*sqrt(discr))/(2*alpha);
pm2= -3*pm3 - (1/16) + (1/8)*(K+r(i).^2) + (1/12)*(3*K-1)*r(i) + (1/12)*r(i).^3;
pm1=  2*pm3 + (1/4)                   +  (1/6)*(4-3*K)*r(i) -  (1/6)*r(i).^3;
p  =  2*pm3 + (5/8)  - (1/4)*(K+r(i).^2);
pp1= -3*pm3 + (1/4)                   -  (1/6)*(4-3*K)*r(i) +  (1/6)*r(i).^3;
pp2=    pm3 - (1/16) + (1/8)*(K+r(i).^2) - (1/12)*(3*K-1)*r(i) - (1/12)*r(i).^3;

val(i) = (-3<R & R<=-2) .* pm3 + ...
      (-2<R & R<=-1) .* pm2 + ...
      (-1<R & R<=0) .* pm1 + ...
      (0<R & R<=1)  .* p   + ...
      (1<R & R<=2)  .* pp1 + ...
      (2<R & R<=3)  .* pp2;

end
end


function Integral = compute_2d_integral(f,delta,x,y)
%f is a square matrix we wish to interpolate onto a point
%delta1 and delta2 are the two delta functions we use to interpolate with
%in the x and y directions, respectively. 
dx = x(2)-x(1);
dy = y(2)-y(1); %uniform grid...

Integral = 0;
for i = 1:length(x)
    for j = 1:length(y)
        Integral = Integral + dx*dy*f(j,i)*delta(j,i);
    end
end


end


function IB6_y = IB6_kernel_composed(x,y,K,h)
IB6_y = zeros(length(y),length(x));
for i = 1:length(x)
    for j = 1:length(y)
        IB6_y(j,i) = (1/h).*IB6((y(j) - sqrt(1 - x(i)^2))./h,K);
    end
end
end


function val = compute_IB6_composed_interpolation(x,y,h,K,u)
val = 0;
for i = 1:length(x)
    for j = 1:length(y)
        val = val + IB6(x(i)/h,K).*IB6((y(j)-sqrt(1-x(i)*x(i)))/h,K).*u(j,i);
    end
end

end


function val = Leveque_Beyer_kernel(r)
val = zeros(1,length(r));
for i = 1:length(r)
    if abs(r(i)) > 2
        val(i) = 0;
    elseif abs(r(i)) <= 2 && abs(r(i)) > 1
        val(i) = 2 - (3*abs(r(i))) + abs(r(i))^2;
    elseif abs(r(i)) <= 1
        val(i) = 1 - abs(r(i))^2;
    end
end

end


function val = compute_L_B_composed_interp(x,x_0,y,h,u)
val = 0;
for i = 1:length(x)
    for j = 1:length(y)
        val = val + Leveque_Beyer_kernel((x(i)-x_0)/h).*Leveque_Beyer_kernel((y(j)-sqrt(1-x(i)*x(i)))/h).*u(j,i);
    end
end

end