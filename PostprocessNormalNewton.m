function solution = PostprocessNormalNewton(sensitivities, criterion, initial, cost)

N = size(sensitivities, 1);

maxiter = 500;
iter = 0;
iterls = 0;
ssn_c = 1;

iterate = full(initial);
Id = eye(N);
Prox = @(x) max(x, 0);
Cost = @(x) cost*norm(x, 1);

% the solution will be Prox(iterate)

[F, AuxF] = criterion.value(Prox(iterate), sensitivities);
j = F + Cost(Prox(iterate));

tol = (1 + abs(j)) * 1e-9;
converged = false;

while ~converged && iter <= maxiter && iterls < 10*maxiter
    iter = iter+1;

    %ComputeGradient and inactiveset
    [gradF, AuxF] = criterion.eval_grad(Prox(iterate), sensitivities, AuxF);
    
    residual = ssn_c*(iterate - Prox(iterate)) + gradF + cost;

    supp = find(Prox(iterate));
    DProx = sparse(supp, supp, ones(size(supp)), N, N);
    
    %% Use L^2 norm gradient as convergence crit
    %fprintf('\tSSN niter: %i, ls: %i, res: %1.1e, supp: %d\n', iter, iterls, norm(residual, 2), length(supp));
    if norm(residual, 2) <= tol
        converged = true;
        % do one more iteration to hopefully get machine precision
    end
   
    HessF = criterion.eval_Hess(Prox(iterate), sensitivities, AuxF);
    Dres = ssn_c*(Id - DProx) + HessF*DProx;
    
%     %% gradient and Hessian check
%     dir = -residual/norm(residual);
%     deltas = 10.^(-10:5);
%     gradres = zeros(size(deltas));
%     hessres = zeros(size(deltas));
%     for cni = 1:length(deltas)
%         Fplus = criterion.value(Prox(iterate) + deltas(cni)*dir, sensitivities);
%         Fminus = criterion.value(Prox(iterate) - deltas(cni)*dir, sensitivities);
%         
%         gradres(cni) = (Fplus - Fminus)/(2*deltas(cni)) - gradF'*dir;
%         hessres(cni) = (Fplus - 2*F + Fminus)/(deltas(cni)^2) - dir'*HessF*dir;
%     end
%     figure(10000)
%     loglog(deltas, abs(gradres)/abs(gradF'*dir), 'k', deltas, abs(hessres)/abs(dir'*HessF*dir), 'g');
%     drawnow;

    stab = 1e-12 * norm(Dres, Inf);
    
    %% Compute step
    direction = - (Dres + stab*Id) \ residual;
    
    newiterate = iterate + direction;
    
    [Fnew, AuxFnew] = criterion.value(Prox(newiterate), sensitivities);
    jnew = Fnew + Cost(Prox(newiterate));

    %% Damping
    theta0 = norm(direction) / norm(residual);
    theta = 1/theta0;
    while jnew > (1 + 10*eps) * j
        newiterate = iterate - (Dres + theta*Id) \ residual;
        
        [Fnew, AuxFnew] = criterion.value(Prox(newiterate), sensitivities);
        jnew = Fnew + Cost(Prox(newiterate));

        theta = theta*2;
        iterls = iterls + 1;
    end
    iterate = newiterate;
    AuxF = AuxFnew;
    F = Fnew;
    j = jnew;

    
end

%ComputeGradient and inactiveset
gradF = criterion.eval_grad(Prox(iterate), sensitivities, AuxF);
residual = ssn_c*(iterate - Prox(iterate)) + cost + gradF;

fprintf('\tSSN niter: %i, ls: %i, res: %1.1e, supp: %d\n', iter, iterls, norm(residual,2), nnz(Prox(iterate)));

solution = max(iterate,0);

end

