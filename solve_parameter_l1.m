function [zk, SI_dual] = solve_parameter_l1(kernel, xx, SI, p_d, beta, zinit)
global print_result

Nq = length(zinit)/2;
zk = zinit;

qk = zk(1:Nq);
yk = zk(Nq+1:end);

%% change to Robinson variables
uk = sign(qk) + qk;
P = @(u) max(0, abs(u) - 1) .* sign(u);
%% q = P(u)
assert(norm(qk - P(uk), inf) < 1e-14);
qk = P(uk);
zk = [uk; yk];

[K, dK] = kernel.matrix(xx, yk);
Gp = [K, dK .* qk'];

II = (1/beta) * (Gp' * SI * Gp);

pk = K * qk;
misfit = pk - p_d;
j = 1/(2*beta) * (misfit' * (SI * misfit)) + norm(qk, 1);

for k = 1:200
  Gp = [K, dK .* qk'];
  R = (1/beta) * (Gp' * (SI * misfit)) + [(uk - qk); zeros(Nq,1)];

  II = Gp' * SI * Gp;
  kp = dK' * (SI * misfit);
  sigma = 0.2;
  kpp = 0.1 * norm(SI * misfit, 1) * abs(qk) / (sqrt(2*pi)*sigma^3); 
  Icor = [ zeros(Nq,Nq), diag(kp);
           diag(kp),     diag(kpp) ];
  HH = (1/beta) * ( II + Icor );
  
  %% SSN correction  
  %DP = diag([abs(uk) > 1; ones(Nq,1)]);
  DP = diag([abs(uk) >= 1; abs(qk) > 0]);
  DR = HH * DP + (eye(2*Nq) - DP);
  
  %Id_R = norm(DR, 1) * eye(2*Nq);
  %DR = DR + 1e-10*Id_R;
  
  dz = - DR \ R;

  %eig(DP(1:Nq,1:Nq) * DR(1:Nq,1:Nq))

  zold = zk;
  jold = j;
  theta = 1 - 1e-14;
  has_descent = false;
  while ~has_descent && theta > 1e-6
    %DRt = theta*DR + (1-theta)*Id_R;
    %dz = - (DRt \ R);
    zk = zold + theta * dz;
    uk = zk(1:Nq);
    qk = P(uk);
    yk = zk(Nq+1:end);
    [K, dK] = kernel.matrix(xx, yk);

    pk = K * qk;
    misfit = pk - p_d;
    j = 1/(2*beta) * (misfit' * (SI * misfit)) + norm(qk, 1);
    descent = j - jold;
    pred = theta * (R' * (DP * dz));
    has_descent = (descent <= pred / 3);
    if ~has_descent
      theta = theta / 1.5;
    end
  end
  if print_result == 1
    fprintf('SSGN iter: %i, j=%f, supp=%d, desc=%1.1e, dz: %1.1e, theta: %1.1e\n', ...
           k, j, sum(abs(uk)>1), descent, norm(dz), theta);
  end
    
  if abs(pred) / theta < 1e-9 && norm(dz) < 1e-4
    %% tolerance reached
    break
  end
  
%   if ~has_descent
%     %% linesearch failed, should not happen
%     keyboard
%   end
end

%% undo Robinson variables
zk = [qk; yk];
%% dual variable
SI_dual = - 1/beta * (SI * misfit);

end
