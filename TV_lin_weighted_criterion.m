function criterion = TV_lin_weighted_criterion(sig_vec, weight)

criterion = struct('value', @(x, sensitivities) value(x, sensitivities, sig_vec, weight), ...
    'eval_grad', @eval_grad, ...
    'eval_Hess', @eval_Hess, ...
    'eval_Hess_vec', @eval_Hess_vec);

end

function [psi, Aux] = value(measure, sensitivities, sig_vec, weight)

Aux = struct();
Aux.sig_vec = sig_vec;
Aux.weight = weight;

Aux.Wsqrt = diag([1./(2*weight(:)); weight(:)]);

% indices = find(measure);
% V = sensitivities(indices, :);
V = sensitivities;
if isempty(measure.x)
    Aux.C = zeros(size(sensitivities, 2));
else
    Aux.C = full(V' * bsxfun(@times, measure.u, V));
end

Aux.C = (Aux.C + Aux.C')/2;

if isempty(Aux.C) || any(isnan(Aux.C(:))) || any(isinf(Aux.C(:)))
    Aux.Cinv = Inf*ones(size(Aux.C));
    Aux.Ci_s = Inf*ones(size(Aux.sig_vec));
    psi = Inf;
else
    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    if isOctave
        [L, D, ~] = svd(Aux.C);
    else
        %% LDL is faster, but not in octave
        [L, D] = ldl(Aux.C);
    end
    
    if any(diag(D) <= 0)
        Aux.Cinv = Inf*ones(size(Aux.C));
        Aux.Ci_s = Inf*ones(size(Aux.sig_vec));
        psi = Inf;
    else
        Linv = inv(L);
        Aux.Cinv = Linv' * diag(1./diag(D)) * Linv;
        Aux.Ci_s = Aux.Cinv * Aux.sig_vec;
        psi = trace(Aux.Wsqrt * Aux.Cinv * Aux.Wsqrt) + sum((Aux.Wsqrt * Aux.Ci_s).^2);
    end
end

end

function [grad, Aux] = eval_grad(~, sensitivities, Aux)

% gradient of the trace criterion
%Aux is computed in value

%d/dtau tr( Wsqrt * C(mu + tau dmu)^{-1} * Wsqrt )
% = - tr( Wsqrt * C(mu)^{-1} * C(dmu) * C(mu)^{-1} * Wsqrt )
% = - tr( Wsqrt * C(mu)^{-1} * S' * Dmu * S * C(mu)^{-1} * Wsqrt )
% = - norms ( Wsqrt * C(mu)^{-1} * S' )^2 * dmu

%d/dtau norm( Wsqrt * C(mu + tau dmu)^{-1} * sig_vec ).^2
% = - 2 * [ Wsqrt * C(mu)^{-1} * sig_vec ]' * [ Wsqrt * C(mu)^{-1} * C(dmu) * C(mu)^{-1} * sig_vec ]
% = - 2 * sig_vec' * C(mu)^{-1} * W * C(mu)^{-1} * S' * Dmu * S * C(mu)^{-1} * sig_vec 
% = - 2 * tr( S * C(mu)^{-1} * sig_vec * sig_vec' * C(mu)^{-1} * W * C(mu)^{-1} * S'  *  Dmu )
  
Aux.K = sensitivities * Aux.Cinv;

grad_trCi = - sum((Aux.K * Aux.Wsqrt).^2, 2);
grad_Ci_s = - 2 * (Aux.K * Aux.sig_vec) .* (Aux.K * Aux.Wsqrt.^2 * Aux.Ci_s);
grad = grad_trCi + grad_Ci_s;

end

function Hess = eval_Hess(~, sensitivities, Aux)

% Hessian of the trace criterion
%Aux is computed in value and eval_grad


%d/dtau - norms ( Wsqrt * C(mu + tau tmu)^{-1} * S' )^2 * dmu
% = 2 * tr ( [ Wsqrt * C(mu)^{-1} * S' ]' * [ Wsqrt * C(mu)^{-1} * C(tmu) * C(mu)^{-1} * S' ] * Dmu )
% = 2 * tr ( K * W * K' * Tmu * S * K' * Dmu )

H1 = (sensitivities * Aux.K');   % NB: this is symmetric
H2 = (Aux.K * Aux.Wsqrt.^2 * Aux.K');  
Hess_trCi = 2 * H1 .* H2;

%d/dtau - 2 * tr( S * C(mu + tau tmu)^{-1} * sig_vec * sig_vec' * C(mu + tau tmu)^{-1} * W * C(mu + tau tmu)^{-1} * S'  *  Dmu )
%= 2 * tr( K * C(tmu) * C(mu)^{-1} * s_v*s_v' * C(mu)^{-1}                       * W * K'                        *  Dmu )
%+ 2 * tr( K *                       s_v*s_v' * C(mu)^{-1} * C(tmu) * C(mu)^{-1} * W * K'                        *  Dmu )
%+ 2 * tr( K *                       s_v*s_v' * C(mu)^{-1}                       * W * C(mu)^{-1} * C(tmu) * K'  *  Dmu )
%= 2 * tr( K * S' * Tmu * K * s_v*s_v' * C(mu)^{-1}            * W * K'                        *  Dmu )
%+ 2 * tr( K *                s_v*s_v' * K' * Tmu * C(mu)^{-1} * W * K'                        *  Dmu )
%+ 2 * tr( K *                s_v*s_v' * C(mu)^{-1}            * W * K' * Tmu * S * K'  *  Dmu )

vec_sig_1 = Aux.K * Aux.sig_vec;
vec_sig_2 = Aux.K * Aux.Wsqrt.^2 * Aux.Ci_s;
Hess_Ci_s = 2 * ( H2 .* ( vec_sig_1 * vec_sig_1' ) ...
                + H1 .* ( vec_sig_2 * vec_sig_1' + vec_sig_1 * vec_sig_2' ) );

Hess = Hess_trCi + Hess_Ci_s;

end

function Hess_vec = eval_Hess_vec(direction, sensitivities, Aux)

  error('TV_lin_criterion:eval_Hess_vec not implemented');
  
%Aux is computed in value
C_of_dir = sensitivities'*(bsxfun(@times, direction, sensitivities));
help = Aux.Cinv * (Aux.Cinv * (C_of_dir * Aux.Cinv));
Hess_vec = 2*sum(sensitivities .* (sensitivities * help'), 2);

end
