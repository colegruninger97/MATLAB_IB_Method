%Script used to test the error analysis derived for group meeting Mar 16
%64: [2.094161180335013 2.094302742901988 2.094376782818399
%2.094388440265444 2.094393846778617 2.094394935313958 2.094395034333985]
%Test 1-D function to interpolate
h = 1/128;
x = -3:h:3;
y = -3:h:3;
%K = -0.815559903284; %%Modified IB6 kernel which should be second order accurate [0.0064 0.0016 3.9822e-4 9.9556e-5 2.4889e-05]
K = 59/60-sqrt(29)/20;
% delta_cos = @(x) (1/h).*cos_kernel(x./h);
% delta_pl = @(x) (1/h).*piecewise_linear_kernel(x./h);
% delta_BS3 = @(x) (1/h).*BS3(x./h);
% delta_BS2 = @(x) (1/h).*BS2(x./h);
% delta_IB5 = @(x) (1/h).*delta_5_smooth(x./h);
% delta_IB3 = @(x) (1/h).*IB3(x./h);
% delta_IB4 = @(x) (1/h).*IB4(x./h);
% delta_IB6 = @(x) (1/h).*IB6(x./h,K);
% integral(@(x)IB6(x,K).*x,0,10);
%delta_sinc = @(x) (1/h).*sinc(x./h);
% fun_cos = @(x) poiseuille_flow(x).*delta_cos(x);
% fun_pl = @(x)  poiseuille_flow(x).*delta_pl(x);
% fun_BS3 = @(x) poiseuille_flow(x).*delta_BS3(x);
% fun_BS2 = @(x) poiseuille_flow(x).*delta_BS2(x);
% fun_IB5 = @(x) poiseuille_flow(x).*delta_IB5(x);
% fun_IB3 = @(x) poiseuille_flow(x).*delta_IB3(x);
% fun_IB4 = @(x) poiseuille_flow(x).*delta_IB4(x);
% fun_IB6 = @(x) poiseuille_flow(x).*delta_IB6(x);
% fun = @(x,y) hemi_sphere(x,y);
% int_approx = 0;
% for i = 1:length(x)
%     for j = 1:length(y)
%     int_approx = int_approx + fun(x(i),y(j))*h*h;
%     end
% end
% int_approx
%fun_sinc = @(x) abs(x).*delta_sinc(x);
jump = 5.0;
% error_derived_cos = h*(jump)*integral(@(x)cos_kernel(x).*x,0,6,'AbsTol',1e-10);
% error_derive_PL = h*(jump)*integral(@(x)piecewise_linear_kernel(x).*x,0,6,'AbsTol',1e-10);
% error_derive_BS3 = h*(jump)*integral(@(x)BS3(x).*x,0,6,'AbsTol',1e-10);
% error_derive_BS2 = h*(jump)*integral(@(x)BS2(x).*x,0,6,'AbsTol',1e-10);
% error_derive_IB5 = h*(jump)*integral(@(x)delta_5_smooth(x).*x,0,6,'AbsTol',1e-10);
% error_derive_IB3 = h*(jump)*integral(@(x)IB3(x).*x,0,6,'AbsTol',1e-10);
% error_derive_IB4 = h*(jump)*integral(@(x)IB4(x).*x,0,6,'AbsTol',1e-10);
% error_derive_IB6 = h*(jump)*integral(@(x)IB6(x,K).*x,0,6,'AbsTol',1e-10);
%interpolate fun using each kernel function             %h = 1/8 1/16 1/32
% approx_cos_ker = integral(fun_cos,-6,6,'AbsTol',1e-10);%[0.1818 0.0919 0.0462 0.0232 0.0116]
% approx_pl_ker = integral(fun_pl,-6,6,'AbsTol',1e-10);%[0.1029 0.0518 0.0260 0.0130 0.0065]
% approx_BS3_ker = integral(fun_BS3,-6,6,'AbsTol',1e-10);%[0.1432 0.0723 0.0363 0.0182 0.0091]
% approx_BS2_ker = integral(fun_BS2,-6,6,'AbsTol',1e-10);%[0.1250 0.0630 0.0316 0.0158 0.0079]
% approx_IB5_ker = integral(fun_IB5,-6,6,'AbsTol',1e-10);%[0.1742 0.0881 0.0443 0.0222 0.0111]
% approx_IB3_ker = integral(fun_IB3,-6,6,'AbsTol',1e-10);%[0.1372 0.0692 0.0347 0.0174 0.0087]
% approx_IB4_ker = integral(fun_IB4,-6,6,'AbsTol',1e-10);%[0.1819 0.0920 0.0463 0.0232 0.0116]
% approx_IB6_ker = integral(fun_IB6,-6,6,'AbsTol',1e-10);%[0.2088 0.1058 0.053 0.0267 0.0134]
%approx_sinc_ker = integral(fun_sinc,-6,6);



%Now I would like to try the same process but in 2D

xx = -5:0.01:5;

plot(xx,BS2(xx),'LineWidth',1.5)
xlabel('$r$','interpreter','latex','fontsize',20)
ylabel('$\phi^{\mathrm{BS}}_3(r)$','interpreter','latex','fontsize',20);
ax = gca;
ax.FontSize = 18;
axis([-3.5,3.5,-0.1,1.1])



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


function val = poiseuille_flow(x)
val = zeros(1,length(x));
for i = 1:length(x)
    if x(i) < -5 || x(i) > 0
        val(i) = 0;
    else
        val(i) = -x(i)*(x(i)+5);
    end
end
end

function val = hemi_sphere(x,y)
val = zeros(length(y),length(x));
for i =1:length(x)
    for j = 1:length(y)
        if x(i)^2 + y(j)^2 <= 1
            val(j,i) = sqrt(1 - x(i)^2 - y(j)^2);
        else
            val(j,i) = 0;
        end
    end
end
end



