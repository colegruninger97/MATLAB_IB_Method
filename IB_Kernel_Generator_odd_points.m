clear;
clc;

%Set the number of moment conditions needed to be satisfied
m = 31;
Phi =  sym('phi_%d',[1 m+2],'real')';
r = sym('r','real');
A = zeros(m+1,m+2);
if rem(m,2) == 1
    for j = 1:m+2
        A(1,j) = 1;
    end
    for i = 2:m+1
        for j = 1:m+2
            A(i,j) = (ceil(m/2) - (j-1))^(i-1);
        end
    end
else
    
end

A = sym(A);
b = zeros(m+1,1);
b(1) = 1;
b = sym(b);
for j = 2:m+1
    b(j) = r^(j-1);
end

soln1 = solve(A*Phi == b,Phi(setdiff(1:m+2, ceil(m/2)+1)));
%solve the equation for phi(0)
 phi_of_half = subs(soln1.phi_1,r,-1/2);
 phi_of_half = solve(phi_of_half == 0);
 Phi(setdiff(1:m+2,ceil(m/2)+1)) = subs(Phi(setdiff(1:m+2,ceil(m/2)+1)),soln1);
 Phi_temp = subs(Phi,Phi(ceil(m/2)+1),phi_of_half);
 Phi_temp = subs(Phi_temp,r,-1/2);
 C = sum(Phi_temp.^2);
% 
 soln2 = solve(sum(Phi(:).^2) == C,Phi(ceil(m/2)+1));
 phi_r = soln2(2);
 Phi = subs(Phi,Phi(ceil(m/2)+1),phi_r);
 x = sym('x');
 for j = 1:m+2
     %r = x + sym(ceil(m/2)+1 - j  +1 );
     Phi(j) = subs(Phi(j),r,x + sym(ceil(m/2) - j  +1));
 end
 
 syms IB_ker(x)
 IB_ker(x) = 0;
 for j = 1:m+2
     IB_ker(x) = IB_ker(x) + piecewise((x >= - ceil(m/2)-1.5 + j) & (x < - ceil(m/2)-0.5 + j), Phi(j),0);
 end
 