function u_pp = postprocess_measure(ux, uu, pp_radius)

  %% take a given discrete measure u and lump together Dirac delta functions
  %% that are less than pp_radius apart
  
  [x, perm] = sort(ux');
  u = uu(perm)';

  cut = [0, find(diff(x) > pp_radius), length(x)];
  if length(x) == 0;
    cut = 0;
  end
  U = zeros(1, length(cut)-1);
  X = zeros(1, length(cut)-1);
  for cp = 1:length(cut)-1
    range = cut(cp)+1:cut(cp+1);
    %% lumped coefficient
    U(1,cp) = sum(u(1,range));
    %% center of gravity point
    if sum(abs(u(1,range))) > 0
      X(1,cp) = sum(x(1,range).*abs(u(1,range))) / sum(abs(u(1,range)));
    else
      X(1,cp) = mean(x(1,range))
    end
  end
  
  u_pp = struct('x', X, 'u', U);
  
end
