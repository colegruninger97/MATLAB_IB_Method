%Normal modes IB kernel finding. For each IB kernel we construct the matrix
%S^*S 

%Discretize the cylinder
mfac = 0.5;
h = 1/16;
x = -2:h:2;
y = -2:h:2;
DX = mfac*h;
r_cyl = 1.0;
kappa = 1000;
N_IB = round(2*pi/DX);
ds = 2*pi/(N_IB);
K = 59/60 - sqrt(29)/20;
s = 0:ds:2*pi-ds;
X = r_cyl.*cos(s);
Y = r_cyl.*sin(s);
%Compute the number of IB points used to discretize the cylinder
N_IB = length(X);
%Form the matrix for the cosine kernel first (it is easiest to implement)
M1 = zeros(N_IB,N_IB);

%Also compute the interpolation operator 
Interp1 = zeros(N_IB,length(x)*length(y));
Interp2 = zeros(N_IB,length(x)*length(y));
Interp3 = zeros(N_IB,length(x)*length(y));


for j = 1:N_IB
    for k = 1:length(y)
        for i = 1:length(x)
        Interp1(j,length(x)*(k-1) + i) = cos_kernel((x(i)-X(j))/h)*cos_kernel((y(k)-Y(j))/h);
        end
    end
end

for j = 1:N_IB
    for k = 1:length(y)
        for i = 1:length(x)
        Interp2(j,length(x)*(k-1) + i) = piecewise_linear_kernel((x(i)-X(j))/h)*piecewise_linear_kernel((y(k)-Y(j))/h);
        end
    end
end
    
for j = 1:N_IB
    for k = 1:length(y)
        for i = 1:length(x)
        Interp3(j,length(x)*(k-1) + i) = IB6((x(i)-X(j))/h,K)*IB6((y(k)-Y(j))/h,K);
        end
    end
end
    
    

% for q = 1:N_IB
%     for w = 1:N_IB
%         for i = 1:length(x)
%             for j = 1:length(y)
%                 M1(q,w) = M1(q,w) + (kappa*DX/(h*h))*cos_kernel((x(i)-X(q))/h)*cos_kernel((y(j)-Y(q))/h)*cos_kernel((x(i)-X(w))/h)*cos_kernel((y(j)-Y(w))/h);
%             end
%         end
%         
%         
%     end
% end
% 
% 
% M2 = zeros(N_IB,N_IB);
% for q = 1:N_IB
%     for w = 1:N_IB
%         for i = 1:length(x)
%             for j = 1:length(y)
%                 M2(q,w) = M2(q,w) + (kappa*DX/(h*h))*piecewise_linear_kernel((x(i)-X(q))/h)*piecewise_linear_kernel((y(j)-Y(q))/h)*piecewise_linear_kernel((x(i)-X(w))/h)*piecewise_linear_kernel((y(j)-Y(w))/h);
%             end
%         end
%         
%         
%     end
% end
% 
% M3 = zeros(N_IB,N_IB);
% for q = 1:N_IB
%     for w = 1:N_IB
%         for i = 1:length(x)
%             for j = 1:length(y)
%                 M3(q,w) = M3(q,w) + (kappa*DX/(h*h))*IB6((x(i)-X(q))/h,K)*IB6((y(j)-Y(q))/h,K)*IB6((x(i)-X(w))/h,K)*IB6((y(j)-Y(w))/h,K);
%             end
%         end
%         
%         
%     end
% end
% 
% [V1,D1] = eig(M1);
% [V2,D2] = eig(M2);
% [V3,D3] = eig(M3);

% %% 
% 
% 
% %Discretize the cylinder
% mfac2 = 2;
% h = 1/16;
% x = -2:h:2;
% y = -2:h:2;
% DX = mfac2*h;
% r_cyl = 1.0;
% kappa = 100;
% K = 59/60 - sqrt(29)/20;
% ds2 = DX/(r_cyl);
% s2 = 0:ds2:2*pi-ds2;
% X2 = r_cyl.*cos(s2);
% Y2 = r_cyl.*sin(s2);
% %Compute the number of IB points used to discretize the cylinder
% N_IB2 = length(X2);
% %Form the matrix for the cosine kernel first (it is easiest to implement)
% M12 = zeros(N_IB2,N_IB2);
% 
% for q = 1:N_IB2
%     for w = 1:N_IB2
%         for i = 1:length(x)
%             for j = 1:length(y)
%                 M12(q,w) = M12(q,w) + (kappa*DX/(h*h))*cos_kernel((x(i)-X2(q))/h)*cos_kernel((y(j)-Y2(q))/h)*cos_kernel((x(i)-X2(w))/h)*cos_kernel((y(j)-Y2(w))/h);
%             end
%         end
%         
%         
%     end
% end
% 
% 
% M22 = zeros(N_IB2,N_IB2);
% for q = 1:N_IB2
%     for w = 1:N_IB2
%         for i = 1:length(x)
%             for j = 1:length(y)
%                 M22(q,w) = M22(q,w) + (kappa*DX/(h*h))*piecewise_linear_kernel((x(i)-X2(q))/h)*piecewise_linear_kernel((y(j)-Y2(q))/h)*piecewise_linear_kernel((x(i)-X2(w))/h)*piecewise_linear_kernel((y(j)-Y2(w))/h);
%             end
%         end
%         
%         
%     end
% end
% 
% M32 = zeros(N_IB2,N_IB2);
% for q = 1:N_IB2
%     for w = 1:N_IB2
%         for i = 1:length(x)
%             for j = 1:length(y)
%                 M32(q,w) = M32(q,w) + (kappa*DX/(h*h))*IB6((x(i)-X2(q))/h,K)*IB6((y(j)-Y2(q))/h,K)*IB6((x(i)-X2(w))/h,K)*IB6((y(j)-Y2(w))/h,K);
%             end
%         end
%         
%         
%     end
% end
% 
% [V12,D12] = eig(M12);
% [V22,D22] = eig(M22);
% [V32,D32] = eig(M32);
% 
% 
% 
% 
% %Discretize the cylinder
% mfac3 = 0.125;
% h = 1/16;
% x = -2:h:2;
% y = -2:h:2;
% DX = mfac3*h;
% r_cyl = 1.0;
% kappa = 100;
% K = 59/60 - sqrt(29)/20;
% ds3 = DX/(r_cyl);
% s3 = 0:ds3:2*pi-ds3;
% X3 = r_cyl.*cos(s3);
% Y3 = r_cyl.*sin(s3);
% %Compute the number of IB points used to discretize the cylinder
% N_IB3 = length(X3);
% %Form the matrix for the cosine kernel first (it is easiest to implement)
% M13 = zeros(N_IB3,N_IB3);
% 
% for q = 1:N_IB3
%     for w = 1:N_IB3
%         for i = 1:length(x)
%             for j = 1:length(y)
%                 M13(q,w) = M13(q,w) + (kappa*DX/(h*h))*cos_kernel((x(i)-X3(q))/h)*cos_kernel((y(j)-Y3(q))/h)*cos_kernel((x(i)-X3(w))/h)*cos_kernel((y(j)-Y3(w))/h);
%             end
%         end
%         
%         
%     end
% end
% 
% 
% M23 = zeros(N_IB3,N_IB3);
% for q = 1:N_IB3
%     for w = 1:N_IB3
%         for i = 1:length(x)
%             for j = 1:length(y)
%                 M23(q,w) = M23(q,w) + (kappa*DX/(h*h))*piecewise_linear_kernel((x(i)-X3(q))/h)*piecewise_linear_kernel((y(j)-Y3(q))/h)*piecewise_linear_kernel((x(i)-X3(w))/h)*piecewise_linear_kernel((y(j)-Y3(w))/h);
%             end
%         end
%         
%         
%     end
% end
% 
% M33 = zeros(N_IB3,N_IB3);
% for q = 1:N_IB3
%     for w = 1:N_IB3
%         for i = 1:length(x)
%             for j = 1:length(y)
%                 M33(q,w) = M33(q,w) + (kappa*DX/(h*h))*IB6((x(i)-X3(q))/h,K)*IB6((y(j)-Y3(q))/h,K)*IB6((x(i)-X3(w))/h,K)*IB6((y(j)-Y3(w))/h,K);
%             end
%         end
%         
%         
%     end
% end
% 
% [V13,D13] = eig(M13);
% [V23,D23] = eig(M23);
% [V33,D33] = eig(M33);


function val = cos_kernel(r)
    %Width of the kernel is 4
    
    if r>=2 || r<= -2
        val = 0;
    else
        val = 0.25*(1 + cos(pi*r/2));
    end
   

end


function val = piecewise_linear_kernel(r)

    if r>=1 || r<= -1
        val = 0;
    elseif r<1 && r>=0
        val = 1-r;
    elseif r<0 && r> -1
        val = r+1;
    end
end

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







