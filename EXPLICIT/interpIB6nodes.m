function Phi = interpIB6nodes(phi,X,dx,dy)
%Right now doing this for the BSpline2,3 couple. Other kernels will be
%incorporated later

%Also, this is currently implemented for periodic boundary conditions only,
%I will update to more general physical boundary conditions later...in
%which case I will need seperate function calls for different bcs

NIB = length(X(:,1));
[rphi,cphi] = size(phi);
K = 59/60-sqrt(29)/20; 
Phi = zeros(NIB,1);


for k = 1:NIB
   sxx = X(k,1)/dx;
   syy = X(k,2)/dy;
   Ixx = floor(sxx);
   Iyy = floor(syy);
   rx = sxx - Ixx; %compute the relevant r value defined in (0,1) in both x and y
   ry = syy - Iyy;
   %Get the relevant indices to interpolate onto the lagrangian structure
   %need to use modular arithmetic to account for periodic boundary
   %conditions
   
   %Need to include an if statement for the Bspline2 kernel implementation

   j1 = mod(Ixx-2:Ixx+3,cphi)+1;
   i1 = mod(Iyy-2:Iyy+3,rphi)+1;
   
   %get the weights
   ws(:,:) = IB6x(rx,K).*IB6y(ry,K);
   Phi(k,1) = sum(sum(phi(i1,j1).*ws(:,:)));
    
end




end

