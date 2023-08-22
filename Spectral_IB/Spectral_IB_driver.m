%Spectral IB driver. 


%Goals:
%Compute the periodic Green's function. Compare solution of periodic Greens
%function to the same Green's function convovled with a regularized delta
%function. Use this comparison to understand a bit more about the accuracy
%of the spreading operator. 




%Set up a periodic grid for the Green's function to be evaluated on
N = 64;
dx = 2*pi/N;
dy = 2*pi/N;
x = -pi:dx:pi-dx;
y = -pi:dy:pi-dy;
k = -pi/dx:1:round((pi-dx)/dx);
l = -pi/dx:1:round((pi-dy)/dy);
u = zeros(length(y),length(x));
ds = 0.5*dx;
theta = 0:ds:2*pi-ds;
X = cos(theta);
Y = sin(theta); %Place in the unit circle

%Construct a force which points in the normal direction along the circle
F1 = cos(theta);
F2 = sin(theta);


%Compute the periodic green's function on this domain
%Define the wave numbers on the grid

 %K = 59/60-sqrt(29)/20;
K = 0;

delta_h_cosine = @(x,y,X,Y) (1/(dy*dx)).*cos_kernel(((x-X)./dx)).*cos_kernel((y-Y)./dy)';
delta_h_BS2 = @(x,y,X,Y) (1/(dy*dx)).*BS2((x-X)./dx).*BS2((y-Y)./dy)';
delta_h_BS3 = @(x,y,X,Y) (1/(dy*dx)).*BS3((x-X)./dx).*BS3((y-Y)./dy)';
delta_h_IB4 = @(x,y,X,Y) (1/(dy*dx)).*IB4((x-X)./dx).*IB4((y-Y)./dy)';
delta_h_IB6 = @(x,y,X,Y) (1/(dy*dx)).*IB6(((x-X))./dx,K).*IB6((y-Y)./dy,K)';
delta_h_BS2BS3 = @(x,y,X,Y) (1/(dy*dx)).*BS2((x-X)./dx).*BS3((y-Y)./dy)';
delta_h_BS3BS2 = @(x,y,X,Y) (1/(dy*dx)).*BS3((x-X)./dx).*BS2((y-Y)./dy)';
delta_h_LB = @(x,y,X,Y) (1/(dy*dx)).*Leveque_Beyer_kernel((x-X)./dx).*Leveque_Beyer_kernel((y-Y)./dy)';

%Make the right hand side forcing term to convolve the green's function
%with
RHSx = zeros(length(y),length(x));
RHSy = zeros(length(y),length(x));

graph = @(x) sqrt(1-x.^2); %graph for the circle
for i = 1:length(theta) %Implicitly uses trapzoidal rule on the circle...
    RHSx(:,:) = RHSx(:,:) + delta_h_IB6(x,y,X(i),Y(i)).*F1(i)*ds;
    RHSy(:,:) = RHSy(:,:) + delta_h_IB6(x,y,X(i),Y(i)).*F2(i)*ds;
end

% for i = 1:length(theta) %Implicitly uses trapzoidal rule on the circle...
%     for n = 1:length(y)
%         for m = 1:length(x)
%             if theta(i) <= pi/4
%                 RHSx(n,m) = RHSx(n,m) + (1/(dy*dx)).*Leveque_Beyer_kernel((x(m)-graph(y(n)))./dx).*Leveque_Beyer_kernel((y(n)-Y(i))./dy)' .* (F1(i).*ds);
%                 RHSy(n,m) = RHSy(n,m) +  (1/(dy*dx)).*Leveque_Beyer_kernel((x(m)-graph(y(n)))./dx).*Leveque_Beyer_kernel((y(n)-Y(i))./dy)' .* (F2(i).*ds);
%             elseif theta(i) > pi/4 && theta(i) <= 3*pi/4
%                 RHSx(n,m) = RHSx(n,m) + (1/(dy*dx)).*Leveque_Beyer_kernel((x(m)-X(i))./dx).*Leveque_Beyer_kernel((y(n)-graph(x(m)))./dy)' .* (F1(i).*ds);
%                 RHSy(n,m) = RHSy(n,m) + (1/(dy*dx)).*Leveque_Beyer_kernel((x(m)-X(i))./dx).*Leveque_Beyer_kernel((y(n)-graph(x(m)))./dy)' .* (F2(i).*ds);
%             elseif theta(i) > 3*pi/4 && theta(i) <= 5*pi/4
%                 RHSx(n,m) = RHSx(n,m) + (1/(dy*dx)).*Leveque_Beyer_kernel((x(m)+graph(y(n)))./dx).*Leveque_Beyer_kernel((y(n)-Y(i))./dy)' .* (F1(i).*ds);
%                 RHSy(n,m) = RHSy(n,m) + (1/(dy*dx)).*Leveque_Beyer_kernel((x(m)+graph(y(n)))./dx).*Leveque_Beyer_kernel((y(n)-Y(i))./dy)' .* (F2(i).*ds);
%             elseif theta(i) > 5*pi/4 && theta(i) <= 7*pi/4
%                 RHSx(n,m) = RHSx(n,m) + (1/(dy*dx)).*Leveque_Beyer_kernel((x(m)-X(i))./dx).*Leveque_Beyer_kernel((y(n)+graph(x(m)))./dy)' .* (F1(i).*ds);
%                 RHSy(n,m) = RHSy(n,m) + (1/(dy*dx)).*Leveque_Beyer_kernel((x(m)-X(i))./dx).*Leveque_Beyer_kernel((y(n)+graph(x(m)))./dy)' .* (F2(i).*ds);
%             elseif theta(i) > 7*pi/4
%                 RHSx(n,m) = RHSx(n,m) + (1/(dy*dx)).*Leveque_Beyer_kernel((x(m)-graph(y(n)))./dx).*Leveque_Beyer_kernel((y(n)-Y(i))./dy)' .* (F1(i).*ds);
%                 RHSy(n,m) = RHSy(n,m) +  (1/(dy*dx)).*Leveque_Beyer_kernel((x(m)-graph(y(n)))./dx).*Leveque_Beyer_kernel((y(n)-Y(i))./dy)' .* (F2(i).*ds);
%             end
%         end
%     end
% end

%Take the fft of this forcing term
RHSx_hat = fft2(RHSx);
RHSy_hat = fft2(RHSy);
RHSx_hat = fftshift(RHSx_hat);
RHSy_hat = fftshift(RHSy_hat);


%Construct the periodic Green's function component by component
Gxx_hat = zeros(length(l),length(k));
Gxy_hat = zeros(length(l),length(k));
Gyy_hat = zeros(length(l),length(k));
for i = 1:length(k)
    for j = 1:length(l)
        if (k(i) == 0 && l(j) == 0)
            Gxx_hat(j,i) = 0;
            Gxy_hat(j,i) = 0;
            Gyy_hat(j,i) = 0;
        else
            Gxx_hat(j,i) = (1/(2*pi*2*pi)).*(1/((k(i)^2 + l(j)^2)^2)).*(l(j)^2);
            Gxy_hat(j,i) = (1/(2*pi*2*pi)).*(1/((k(i)^2 + l(j)^2)^2)).*(-k(i)*l(j));
            Gyy_hat(j,i) = (1/(2*pi*2*pi)).*(1/((k(i)^2 + l(j)^2)^2)).*(k(i)^2);
        end
    end
end

U_hat = zeros(length(l),length(k));
V_hat = zeros(length(l),length(k));

% construct the fourier transform of the true solution
for q = 1:length(X)
    
    for j = 1:length(l)
        for i = 1:length(k)
          U_hat(j,i) = U_hat(j,i) + F1(q).*Gxx_hat(j,i).*exp(-1i*k(i)*X(q)).*exp(-1i*l(j)*Y(q)).*ds + F2(q).*Gxy_hat(j,i).*exp(-1i*k(i)*X(q)).*exp(-1i*l(j)*Y(q)).*ds;
          V_hat(j,i) = V_hat(j,i) + F1(q).*Gxy_hat(j,i).*exp(-1i*k(i)*X(q)).*exp(-1i*l(j)*Y(q)).*ds + F2(q).*Gyy_hat(j,i).*exp(-1i*k(i)*X(q)).*exp(-1i*l(j)*Y(q)).*ds;
        end
        
    end
       
end

%compute the convolution of the regularized delta function with the
%stokeslet

U = (1/dx).*(1/dy).*real(ifftshift(ifft2(fftshift(U_hat))));
V = (1/dx).*(1/dy).*real(ifftshift(ifft2(fftshift(V_hat))));

convx_hat = Gxx_hat.*RHSx_hat + Gxy_hat.*RHSy_hat;
convy_hat = Gxy_hat.*RHSx_hat + Gyy_hat.*RHSy_hat;

convx = real(ifft2(fftshift(convx_hat)));
convy = real(ifft2(fftshift(convy_hat)));

conv_mag = sqrt(convx.^2 + convy.^2);
U_mag = sqrt(U.^2 + V.^2);
err = abs(U_mag);
L1_err = dx*dy*sum(err,'all');
L2_err = sqrt(dx*dy*sum(err.^2,'all'));
L_inf_err = max(err,[],'all');


%Perform convolution of the Green's function with the regularized delta
%function. In fourier space, this is simply multiplacation. 

% for i = 1:length(y)
%     for j = 1:length(x)
%         for q = 0:length(k)-1
%             for m = 0:length(l)-1
%                 u(i,j) = u(i,j) + Gxx_hat(m+1,q+1)*exp(2*pi*1i*q*x(j) + 2*pi*1i*m*y(i));
%             end
%         end
%     end
% end

% U = (1/dx).*(1/dy).*real(ifftshift(ifft2(fftshift(Gxx_hat))));
% V = (1/dx).*(1/dy).*real(ifftshift(ifft2(fftshift(Gxy_hat))));
% U_mag = sqrt(U.^2 + V.^2);

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



