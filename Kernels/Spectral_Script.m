%Normal modes IB kernel finding. For each IB kernel we construct the matrix
%S^*S 

%Discretize the cylinder
mfac = 1;
h = 1/16;
x = -2:h:2;
y = -2:h:2;
DX = mfac*h;
r_cyl = 1.0;
N_IB = round(2*pi/DX);
ds = 2*pi/(N_IB);
K = 59/60 - sqrt(29)/20;
s = 0:ds:2*pi-ds;
X = r_cyl.*cos(s);
Y = r_cyl.*sin(s);
%Compute the number of IB points used to discretize the cylinder
N_IB = length(X);
%Form the matrix for the cosine kernel first (it is easiest to implement)
SstarS = zeros(N_IB,N_IB);

    

for q = 1:N_IB
    for w = 1:N_IB
        for i = 1:length(x)
            for j = 1:length(y)
                SstarS(q,w) = SstarS(q,w) + (DX/(h*h))*cos_kernel((x(i)-X(q))/h)*cos_kernel((y(j)-Y(q))/h)*cos_kernel((x(i)-X(w))/h)*cos_kernel((y(j)-Y(w))/h);
            end
        end
    end
end

[V,D] = eig(SstarS);
plot(real(diag(D)),imag(diag(D)),'o')


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
