function [SI_dual, mu_min_norm] = calculate_certificate(sensor, mesh, q_dagger, y_dagger, kernel)

  %% calculate certificate by solving the noise free problem with small reg parameter

  if isstruct(sensor)
    xx = sensor.x;
    uu = sensor.u;
  else
    %% assume this is a discretized measure
    %% get the parameters of the measure
    supp = (sensor > 0);
    xx = mesh.points(supp);
    uu = sensor(supp);
  end

  [K, dK] = kernel.matrix(xx, y_dagger);
  sens = [K, dK .* q_dagger'];

  sig_vec = [sign(q_dagger); zeros(size(y_dagger))];
  
  %% calculate pre-certificate
  SI = diag(uu);
  SqrtSI = diag(sqrt(uu)); % square root of Sigma^{-1}
  SI_dual_pre = SqrtSI * ((sens' * SqrtSI) \ sig_vec);

  [K_h, dK_h] = kernel.matrix(xx, mesh.points);
  eta = K_h' * SI_dual_pre;
  linf_eta = max(max(eta), max(-eta));

  %% certificate (approximate)
  p_dagger = K * q_dagger;

  beta = 1e-5;  % chose small reg parameter
  mu_init = struct('x', y_dagger(q_dagger ~= 0), 'u', q_dagger(q_dagger ~= 0));
  [mu_min_norm, SI_dual] = solve_TV(kernel, xx, SI, p_dagger, beta, mu_init, mesh.points);

  gap = sum(abs(q_dagger)) - p_dagger'*SI_dual;
  if abs(gap) < 1e-13
    gap = 0;
  end
  fprintf("violation of pre-cert: %e; dual gap: %e\n", max(linf_eta-1-1e-13,0), gap);

end
