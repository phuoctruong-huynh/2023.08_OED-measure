function [q_opt, SI_dual] = solve_TV(kernel, xx, SI, p_d, beta, q_init, x_h)
global print_result
Nq = length(q_init.u);

qk = q_init.u(:);
yk = q_init.x(:);

K_h = kernel.matrix(xx, x_h);

%% change to Robinson variables
uk = sign(qk) + qk;
Prox = @(u) max(0, abs(u) - 1) .* sign(u);
%% q = Proxu)
assert(norm(qk - Prox(uk), inf) < 1e-14);
qk = Prox(uk);
zk = [uk; yk];

[K, dK] = kernel.matrix(xx, yk);
Gp = [K, dK .* qk'];

II = (Gp' * SI * Gp);

pk = K * qk;
misfit = pk - p_d;
j = 1/(2*beta) * (misfit' * (SI * misfit)) + norm(qk, 1);

for k = 1:500
  Gp = [K, dK .* qk'];
  R = (1/beta) * (Gp' * (SI * misfit)) + [(uk - qk); zeros(Nq,1)];
  
  II = Gp' * SI * Gp;
  
  kp = dK' * (SI * misfit);
  sigma = 0.2;
  kpp = .1 * norm(SI * misfit, 1) * abs(qk) / (sqrt(2*pi)*sigma^3); 
  Icor = [ zeros(Nq,Nq), 0*diag(kp);
           0*diag(kp),     diag(kpp) ];
  
  %II = (Gp' * SI * Gp);
  %Icor = beta * diag([zeros(Nq,1); ones(Nq,1)]);
  
  HH = (1/beta) * (II + Icor);

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
  while ~has_descent && theta > 1e-9
    %DRt = theta*DR + (1-theta)*Id_R;
    %dz = - (DRt \ R);
    zk = zold + theta * dz;
    uk = zk(1:Nq);
    qk = Prox(uk);
    yk = zk(Nq+1:end);
    [K, dK] = kernel.matrix(xx, yk);

    pk = K * qk;
    misfit = pk - p_d;
    j = 1/(2*beta) * (misfit' * (SI * misfit)) + norm(qk, 1);
    descent = j - jold;
    pred = theta * (R' * (DP * dz));
    has_descent = (descent <= pred / 3 + 1e-15);
    if ~has_descent
      theta = theta / 1.5;
    end
  end

  %% active set
  suppc = (abs(uk) <= 1);

  %% constraint violation
  eta = 1/beta * K_h'*(SI*misfit);
  sh_eta = abs(Prox(eta));
  [max_sh_eta, ind] = max(sh_eta);
  if print_result == 1
    fprintf('GNAP iter: %i, j=%f, supp=(%d->%d), desc=%1.1e, dz: %1.1e, viol=%1.1e, theta: %1.1e\n', ...
           k, j, Nq, sum(~suppc), descent, norm(dz), max_sh_eta, theta);
  end

  if ~has_descent
    %% linesearch failed, should not happen
    keyboard
  end

  %% prune zero coefficient Diracs
  if any(suppc)
    Nq = sum(~suppc);
    uk(suppc) = [];
    yk(suppc) = [];
    K(:,suppc) = [];
    dK(:,suppc) = [];
    zk = [uk; yk];
    qk = Prox(uk);
  end

  %% try adding promising new zero coeffs
  grad_supp_q = 1/beta * (K' * (SI * misfit)) + (uk - qk);
  tresh_q = sum(abs(grad_supp_q));
  grad_supp_y = 1/beta * (qk .* (dK' * (SI * misfit)));
  tresh_y = sum(abs(grad_supp_y));
  
  if max_sh_eta > tresh_q + tresh_y
    Nq = Nq+1;
    uk = [uk; -sign(eta(ind))];
    qk = [qk; 0];
    yk = [yk; x_h(ind)];
    [K, dK] = kernel.matrix(xx, yk);
    zk = [uk; yk];
    fprintf('  insert: viol=%1.2e, |g_q|+|g_y|=%1.1e+%1.1e, supp:(%d->%d)\n', ...
              max_sh_eta, tresh_q, tresh_y, sum(~suppc), Nq);
    % disp(yk(:)');
  end
   
%  figure(13)
%  h1 = plot(x_h, -eta, 'b--', 'LineWidth', 1);
%  hold on
%  plot([-1;1], [[1;1],[-1;-1]], 'k:', 'LineWidth', 1);
%  h2 = plot_measure(x_h, struct('u', qk, 'x', yk), 's', 'b');
%  drawnow
%  hold off
%  pause(.1)
%  %keyboard

  if abs(pred) / theta < 1e-9 && norm(dz) + max_sh_eta < 1e-4
    %% tolerance reached
    break
  end
end

%% undo Robinson variables
q_opt = struct('u', qk, 'x', yk);
%% dual variable
SI_dual = - 1/beta * (SI * misfit);

end
