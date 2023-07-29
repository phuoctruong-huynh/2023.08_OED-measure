function [gradF, trialpoint, Gap] = FirstOrder(sensitivities, criterion, AuxF, iterate, M, cost)

gradF = criterion.eval_grad(iterate, sensitivities, AuxF);

[minimum, point] = min(gradF);
if minimum >= -cost
    newcoeff = 0;
else
    newcoeff = M;
end

trialpoint = sparse(point, 1, newcoeff, size(iterate,1), 1);

Gap = gradF'*(iterate - trialpoint) + cost*(sum(iterate) - sum(trialpoint));

end