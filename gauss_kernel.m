
function kernel = gauss_kernel(param)
kernel = struct( ...
    'eval', @(x, xp) kernel_eval(x, xp, param), ...
    'matrix', @(x, xp) kernel_matrix(x, xp, param));

end

function [k, dk] = kernel_eval(x, xp, param)
  lk = - (x - xp).^2 / (2 * param.s2) - log(2*pi*param.s2) / 2;
  k = exp(lk);
  dk = k .* (x - xp) / param.s2;
end

function [K, dK] = kernel_matrix(x, xp, param)
  [Xp_h, X_h] = meshgrid(xp, x);
  [K, dK] = kernel_eval(X_h, Xp_h, param);
end
% Note: Only param.s2 is needed.