function F_new = Elastic_Laplacian(X,kappa,ds)
F_new = zeros(length(X(:,1)),2);
%First handle the x coordinate
F_new(2:end-1,1) = (kappa/(ds*ds)).*(X(3:end,1)+X(1:end-2,1)-2.*X(2:end-1,1));
F_new(1,1) = (kappa/(ds*ds)).*(X(2,1)+X(end,1)-2.*X(1,1));
F_new(end,1) = (kappa/(ds*ds)).*(X(1,1)+X(end-1,1)-2.*X(end,1));

F_new(2:end-1,2) = (kappa/(ds*ds)).*(X(3:end,2)+X(1:end-2,2)-2.*X(2:end-1,2));
F_new(1,2) = (kappa/(ds*ds)).*(X(2,2)+X(end,2)-2.*X(1,2));
F_new(end,2) = (kappa/(ds*ds)).*(X(1,2)+X(end-1,2)-2.*X(end,2));

end
