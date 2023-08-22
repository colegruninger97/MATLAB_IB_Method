function F_new = Elastic_force_update(X,kappa,dX)
F_new = zeros(length(X(:,1)),2);
ds = (dX/(1/4));
%First handle the x coordinate
F_new(2:end-1,1) = (kappa/ds).*((X(3:end,1)-X(2:end-1,1))./(sqrt((X(3:end,1)-X(2:end-1,1)).^2 + (X(3:end,2)-X(2:end-1,2)).^2)) - (X(2:end-1,1)-X(1:end-2,1))./(sqrt((X(2:end-1,1)-X(1:end-2,1)).^2 + (X(2:end-1,2)-X(1:end-2,2)).^2)));
F_new(1,1) = (kappa/ds).*((X(2,1)-X(1,1))./(sqrt((X(2,1)-X(1,1)).^2 + (X(2,2)-X(1,2)).^2)) - (X(1,1)-X(end,1))./(sqrt((X(1,1)-X(end,1)).^2 + (X(1,2)-X(end,2)).^2))); 
F_new(end,1) = (kappa/ds).*((X(1,1)-X(end,1))./(sqrt((X(1,1)-X(end,1)).^2 + (X(1,2)-X(end,2)).^2)) - (X(end,1)-X(end-1,1))./(sqrt((X(end,1)-X(end-1,1)).^2 + (X(end,2)-X(end-1,2)).^2))); 

F_new(2:end-1,2) = (kappa/ds).*((X(3:end,2)-X(2:end-1,2))./(sqrt((X(3:end,1)-X(2:end-1,1)).^2 + (X(3:end,2)-X(2:end-1,2)).^2)) - (X(2:end-1,2)-X(1:end-2,2))./(sqrt((X(2:end-1,1)-X(1:end-2,1)).^2 + (X(2:end-1,2)-X(1:end-2,2)).^2)));
F_new(1,2) = (kappa/ds).*((X(2,2)-X(1,2))./(sqrt((X(2,1)-X(1,1)).^2 + (X(2,2)-X(1,2)).^2)) - (X(1,2)-X(end,2))./(sqrt((X(end,1)-X(1,1)).^2 + (X(end,2)-X(1,2)).^2))); 
F_new(end,2) = (kappa/ds).*((X(1,2)-X(end,2))./(sqrt((X(1,1)-X(end,1)).^2 + (X(1,2)-X(end,2)).^2)) - (X(end,2)-X(end-1,2))./(sqrt((X(end,1)-X(end-1,1)).^2 + (X(end,2)-X(end-1,2)).^2))); 
end

