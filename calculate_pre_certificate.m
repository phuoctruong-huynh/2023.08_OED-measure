function SI_dual = calculate_pre_certificate(sensor, mesh, q_dagger, y_dagger, kernel)

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
  
  SqrtSI = diag(sqrt(uu));
  isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
  if isOctave
    SI_dual = SqrtSI * ((sens' * SqrtSI) \ sig_vec);
  else
    SI_dual = SqrtSI * lsqminnorm(sens' * SqrtSI, sig_vec);
  end
end

