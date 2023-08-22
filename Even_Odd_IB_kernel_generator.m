%Script for generating the IB kernel family of regularized delta functions.
%
clear;
clc;

%Set the number of moment conditions needed to be satisfied
m = 29;
Phi =  sym('phi_%d',[1 m+3],'real')';
r = sym('r','real');
A = zeros(m+2,m+3);
if rem(m,2) == 1
    for j = 1:m+3
        A(1,j) = rem(j,2);
        A(2,j) = rem(j+1,2);
    end
    for i = 3:m+2
        for j = 1:m+3
            A(i,j) = (ceil(m/2)+1 - (j-1))^(i-2);
        end
    end
else
    
end

A = sym(A);
b = zeros(m+2,1);
b(1) = 1/2;
b(2) = 1/2;
b = sym(b);
for j = 3:m+2
    b(j) = r^(j-2);
end

soln1 = solve(A*Phi == b,Phi(setdiff(1:m+3, ceil(m/2)+2)));
%solve the equation for phi(0)
 phi_of_zero = subs(soln1.phi_1,r,0);
 phi_of_zero = solve(phi_of_zero == 0);
 Phi(setdiff(1:m+3,ceil(m/2)+2)) = subs(Phi(setdiff(1:m+3,ceil(m/2)+2)),soln1);
 Phi_temp = subs(Phi,Phi(ceil(m/2)+2),phi_of_zero);
 Phi_temp = subs(Phi_temp,r,0);
 C = sum(Phi_temp.^2);
% 
 soln2 = solve(sum(Phi(:).^2) == C,Phi(ceil(m/2)+2));
 phi_r = soln2(2);
 Phi = subs(Phi,Phi(ceil(m/2)+2),phi_r);
 x = sym('x');
 for j = 1:m+3
     %r = x + sym(ceil(m/2)+1 - j  +1 );
     Phi(j) = subs(Phi(j),r,x + sym(ceil(m/2)+1 - j  +1));
 end
 
 syms IB_ker(x)
 IB_ker(x) = 0;
 for j = 1:m+3
     IB_ker(x) = IB_ker(x) + piecewise((x >= - ceil(m/2)-2 + j) & (x < - ceil(m/2)-1 + j), Phi(j),0);
 end
%solve the quadratic equation for phi(r) for r\in[0,1].


%Compute the sum of squares coefficient
%Setup the system matrix to solve