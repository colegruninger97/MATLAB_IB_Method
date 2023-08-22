function F_new = uniform_normal_force_circle(X,kappa,ds,r)
F_new = zeros(length(X(:,1)),2);
for k = 1:length(X(:,1))
    F_new(k,1) = -kappa.*r.*cos((k-1)*ds);
    F_new(k,2) = -kappa.*r.*sin((k-1)*ds);
end
end

