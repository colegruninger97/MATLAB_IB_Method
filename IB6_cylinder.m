h = 1/64;
x = -1:h:1;
y = -3:h:3;
[xx,yy] = meshgrid(x,y);
K = -0.815559903284;
IB6_kernel = @(x,y) (1/1).*(1/1).*IB6((x)./1,K).*IB6((y + sqrt(1-x.^2))./1,K);
IB6_kernel_0 = @(x,y) (1/h).*(1/h).*IB6((x)./h,59/60-sqrt(29)/20).*IB6((y + 1)./h,59/60-sqrt(29)/20);
IB3_kernel = @(x,y) (1/h)*(1/h).*IB3(x./h).*IB3((y+1)./h);
BS3_kernel = @(x,y) (1/h)*(1/h).*BS3(x./h).*BS3((y+1)./h);
ymax = @(x)-sqrt(1-x.^2)+3*h;
ymin = @(x)-sqrt(1-x.^2)-3*h;
u = @(x,y) u_stokes(x,y,1,1).*IB6_kernel(x,y);
v = @(x,y) u_stokes(x,y,1,1).*IB3_kernel(x,y);
vv =  @(x,y) u_stokes(x,y,1,1).*BS3_kernel(x,y);
v_IB6_0 = @(x,y) u_stokes(x,y,1,1).*IB6_kernel_0(x,y);
g = integral2(u,-3*h,3*h,ymin,ymax)%[0.0029 7.98e-4 2.11e-4 5.43e-5]
ff = integral2(v,-1.5*h,1.5*h,-1+-1.5*h,-1+1.5*h)%[0.0192 0.0099 0.0052 0.0026]
ff_BS3 = integral2(vv,-1.5*h,1.5*h,-1+-1.5*h,-1+1.5*h)%[0.0136 0.0096 0.0056 0.0028]
ff_IB6 = integral2(v_IB6_0,-3*h,3*h,-1-3*h,-1+3*h)%[0.0216 0.0109 0.0065 0.0040]


function val=IB6(r,K)

R = r;

% compute r in [0,1]
r = (-3<r & r<=-2).*(r+3) + (-2<r & r<=-1).*(r+2) + ...
    (-1<r & r<=0) .*(r+1) + (0<r & r<=1)  .*(r+0) + ...
    (1<r & r<=2)  .*(r-1) + (2<r & r<=3)  .*(r-2);
alpha=28;
beta=(9/4)-(3/2)*(K+r.^2)+((22/3)-7*K)*r-(7/3)*r.^3;
gamma=(1/4)*( ((161/36)-(59/6)*K+5*K^2)*(1/2)*r.^2 + (-(109/24)+5*K)*(1/3)*r.^4 + (5/18)*r.^6 );

discr=beta.^2-4*alpha*gamma;
% flag=0;
% if(any(discr<0))
%    flag=1
% end

pm3=(-beta+sign((3/2)-K)*sqrt(discr))/(2*alpha);
pm2= -3*pm3 - (1/16) + (1/8)*(K+r.^2) + (1/12)*(3*K-1)*r + (1/12)*r.^3;
pm1=  2*pm3 + (1/4)                   +  (1/6)*(4-3*K)*r -  (1/6)*r.^3;
p  =  2*pm3 + (5/8)  - (1/4)*(K+r.^2);
pp1= -3*pm3 + (1/4)                   -  (1/6)*(4-3*K)*r +  (1/6)*r.^3;
pp2=    pm3 - (1/16) + (1/8)*(K+r.^2) - (1/12)*(3*K-1)*r - (1/12)*r.^3;

val = (-3<R & R<=-2) .* pm3 + ...
      (-2<R & R<=-1) .* pm2 + ...
      (-1<R & R<=0) .* pm1 + ...
      (0<R & R<=1)  .* p   + ...
      (1<R & R<=2)  .* pp1 + ...
      (2<R & R<=3)  .* pp2;

end

function val = u_stokes(x,y,a,mu)
    R = sqrt(x.^2 + y.^2);
    if R <= 1
        val = 0;
    else
        val = 1 - (6*pi*mu*a/(8*pi*mu)).*((x.*x)./(R.^3) +1./R) - 0.25*a*a.*(1./(R.^3) - 3.*(x.*x)./(R.^5));
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



function val=IB3(r)

val =...
    (abs(r)<1/2).*((1/3)*(1+sqrt(1-3*abs(r).^2))) +...
    (1/2<=abs(r) & abs(r)<3/2).*((1/6)*(5-3*abs(r)-sqrt(-3*(1-abs(r)).^2+1)));


end