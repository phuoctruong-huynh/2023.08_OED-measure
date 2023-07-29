%% Example 1B. Reconstruction results with noisy data and uniform 9, 11 sensors.
function result = example_1B_plot(N_sensors, m, beta_0, seed)
global print_result
print_result = 0;
show_plot = 1;                  % Show plot: 1.
quantitative = 0;               % Compute quantitative estimates: 1.
rng(seed)
%% Setup of the problem:

% Number of sources and sensors:
N_sources = 3;
N_gridref = 8;

% Grid for solve_TV:
x_h = linspace(-1, 1, (2^N_gridref + 1)*N_sensors)'; 
mesh = struct('points', x_h);
Nh = size(x_h, 1);

% Initialize reference measure:
y_dagger = [-.7, -.3, .3]';
q_dagger = [.4, .3, -.2]';
z_dagger = [q_dagger; y_dagger];
mu_dagger = struct('x', y_dagger, 'u', q_dagger);

% Initialize sensor placement:
xx = linspace(-1, 1, N_sensors);
% xx = [-.8, -.6, -.4, -.1, .1, .4];
uu = 1/length(xx) * ones(length(xx), 1);
% Artificial init_sensor:

sensor = struct('x', xx, 'u', uu);

SI = diag(sensor.u);
sqrtSI = diag(sqrt(sensor.u));

% Parameters: s^2, sign_vector, beta_0, m.
param = struct();
T = 1/2*(0.2).^2;
sigma = sqrt(2*T);
param.s2 = sigma.^2;
param.sig_vec = [sign(q_dagger); zeros(N_sources, 1)];
param.beta_0 = beta_0;
param.m = m;
m_log = log10(m);

weight = sqrt(abs(q_dagger));
criterion = TV_lin_weighted_criterion(param.beta_0 * param.sig_vec, weight);

%% Calculate estimator error: Choose uniform sensor placement setup:

kernel = gauss_kernel(param);
[K_d, dK_d] = kernel.matrix(xx, y_dagger);
[K, dK] = kernel.matrix(xx, y_dagger);
sensitivities = [K, dK .* q_dagger'];
pd = K_d * q_dagger;

Gp_d = [K_d, dK_d .* q_dagger'];
SI_pre = sqrtSI * ((Gp_d' * sqrtSI) \ param.sig_vec);
II = (Gp_d' * SI * Gp_d);


% Number of realizations:

% Plots:
K_plot = kernel.matrix(sensor.x, mesh.points);
KK_plot = kernel.matrix(sensor.x, y_dagger);

epsilon = param.m^(-1/2) * (sqrtSI \ randn(size(pd)));

beta = param.beta_0 * param.m^(-1/2);
dz = II \ (- beta * param.sig_vec + Gp_d' * (SI * epsilon));

z_pert = z_dagger + dz;

q_pert = z_pert(1:N_sources);
y_pert = z_pert(N_sources+1:end);

pnoise = pd + epsilon;

SI_pre_dual_eps = -1/beta * (SI*(Gp_d*dz - epsilon));
[z_hat, SI_dual_hat] = solve_parameter_l1(kernel, xx, SI, pnoise, beta, z_dagger);

q_hat = z_hat(1:N_sources);
y_hat = z_hat(N_sources+1:end);

[mu_bar, SI_dual_bar] = solve_TV(kernel, xx, SI, pnoise, beta, mu_dagger, mesh.points);

% Plot:
if show_plot == 1
h1 = plot_measure_new(mu_dagger, 'o', 'r', 2, 4, 'r');
hold on

% h2 = plot_measure(mesh.points, struct('x', y_pert, 'u', q_pert), 'd', 'k');
h3 = plot(mesh.points, K_plot' * SI_pre_dual_eps, 'k--', 'LineWidth', 1);

h4 = plot_measure_new(struct('x', y_hat, 'u', q_hat), 's', 'g', 1, 4, 'g');
h5 = plot(mesh.points, K_plot' * SI_dual_hat, 'g-.', 'LineWidth', 2);

h6 = plot_measure_new(mu_bar, 's', 'b', 1, 3, 'b');
h7 = plot(mesh.points, K_plot' * SI_dual_bar, 'b--', 'LineWidth', 1);

h8 = plot_measure(mesh.points, sensor, '*', 0.5*[1 1 1]);

plot([-1,1], [ones(2, 1), -ones(2, 1)], 'k:', 'LineWidth', 1)
axis([-1, 1, -2.5, 1.5])
hold off
set(gca,'TickLabelInterpreter','latex', 'FontName', 'Arial', 'Fontsize', 18)
lgd = legend([h1, h3, h4, h5, h6, h7, h8], 'reference measure', 'noisy pre-certificate', 'nonlinear reconstruction $\hat{\mu}$', 'certificate for $\hat{\mu}$', 'sparse reconstruction $\bar{\mu}$', 'certificate for $\bar{\mu}$', 'sensors');
set(lgd, 'Location', 'southwest')
set(lgd,'Interpreter','latex')
lgd.FontSize = 14;
xlabel([num2str(length(sensor.x)), ' sensors, $p = $ 10\textsuperscript{', num2str(m_log) '} , $\beta_0 = $ ' num2str(param.beta_0)], 'Interpreter','latex', 'FontName', 'Arial')
set(gcf, 'renderer', 'Painters');

% file_name = ['plots', base_name, num2str(k), ext];
end
end