% Copyright (C) 2023 Phuoc Truong Huynh, Konstantin Pieper and Daniel Walter
 
    % This program is free software: you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation, either version 3 of the License, or
    % (at your option) any later version.
 
    % This program is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    % GNU General Public License for more details.
 
    % You should have received a copy of the GNU General Public License
    % along with this program. If not, see <https://www.gnu.org/licenses/>. 

    % Contact information: pieperk@ornl.gov 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    
%% Example 3. Compare (uniform) 11 sensors, (uniform) 6 sensors and selected 6 sensors.
% Change N_sensors or xx for different measurement setups.


global print_result
print_result = 0;
show_plot = 1;                  % Show plot: 1.
quantitative = 0;               % Compute quantitative estimates: 1.

%% Setup of the problem:

% Number of sources and sensors:
N_sources = 3;
N_gridref = 8;
N_sensors = 2*N_sources + 3;

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
param.beta_0 = 2;
param.p = 10^8;

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
Knoise = 20;
dist_dz = zeros(Knoise, 1);
delta_z = zeros(Knoise, 1);
dist_z_hat = zeros(Knoise, 1);
dist_z_bar = zeros(Knoise, 1);
dist_W_dz = zeros(Knoise, 1);
SI_pre_dual_eps = zeros(Knoise, 1);
bad_count = 0;
% Plots:
K_plot = kernel.matrix(sensor.x, mesh.points);
KK_plot = kernel.matrix(sensor.x, y_dagger);
base_name = 'example1B_noisy_6sensors_m1000000_beta2_';
ext = '.png';

randn ('state', 1);

for k = 1:Knoise
  epsilon = param.p^(-1/2) * (sqrtSI \ randn(size(pd)));  
  beta = param.beta_0 * param.p^(-1/2);
  dz = II \ (- beta * param.sig_vec + Gp_d' * (SI * epsilon));

  z_pert = z_dagger + dz;

  q_pert = z_pert(1:N_sources);
  y_pert = z_pert(N_sources+1:end);
  
  dist_dz(k) = compHK(q_dagger, y_dagger, q_pert, y_pert);
  pnoise = pd + epsilon;
  
  SI_pre_dual_eps = -1/beta * (SI*(Gp_d*dz - epsilon));
  [z_hat, SI_dual_hat] = solve_parameter_l1(kernel, xx, SI, pnoise, beta, z_dagger);
  
  q_hat = z_hat(1:N_sources);
  y_hat = z_hat(N_sources+1:end);
  dist_z_hat(k) = compHK(q_dagger, y_dagger, q_hat, y_hat);

  [mu_bar, SI_dual_bar] = solve_TV(kernel, xx, SI, pnoise, beta, mu_dagger, mesh.points);
  if length(mu_bar.u) > length(mu_dagger.u)
      bad_count = bad_count + 1;
  end
  dist_z_bar(k) = compHK(q_dagger, y_dagger, mu_bar.u, mu_bar.x);
  dist_W_dz(k) = sum(dz.^2 .* [1./(4*weight.^2); weight.^2]);
  
  % Plot:
  if show_plot == 1
	figure(k)
    h1 = plot_measure_new(mu_dagger, 'o', 'r', 2, 4, 'r');
    hold on

    h2 = plot_measure(mesh.points, struct('x', y_pert, 'u', q_pert), 'd', 'k');
    h3 = plot(mesh.points, K_plot' * SI_pre_dual_eps, 'k--', 'LineWidth', 1);

    h4 = plot_measure_new(struct('x', y_hat, 'u', q_hat), 's', 'g', 1, 4, 'g');
    h5 = plot(mesh.points, K_plot' * SI_dual_hat, 'g-.', 'LineWidth', 2);

    h6 = plot_measure_new(mu_bar, 's', 'b', 1, 3, 'b');
    h7 = plot(mesh.points, K_plot' * SI_dual_bar, 'b--', 'LineWidth', 1);

    h8 = plot_measure_new(sensor, 's', 'b', 1, 3, 'b');
    plot([-1,1], [ones(2, 1), -ones(2, 1)], 'k:', 'LineWidth', 1)
    axis([-1, 1, -1.5, 1.5])
    hold off
    set(gca,'TickLabelInterpreter','latex', 'FontName', 'Arial', 'Fontsize', 15)
    lgd = legend([h1, h2, h3, h4, h5, h6, h7, h8], 'reference measure', 'linear reconstruction $z^{\dagger} + \delta \hat{z}$', 'noisy pre-certificate', 'nonlinear reconstruction $\hat{\mu}$', 'certificate for $\hat{\mu}$', 'sparse reconstruction $\bar{\mu}$', 'certificate for $\bar{\mu}$', 'sensors');
    set(lgd, 'Location', 'southwest')
    set(lgd,'Interpreter','latex')
    xlabel([num2str(length(sensor.x)), ' sensors, $p = $ ' num2str(param.p) ', $\beta_0 = $ ' num2str(param.beta_0)], 'Interpreter','latex', 'FontName', 'Arial')
    set(gcf, 'renderer', 'Painters');
    % file_name = ['plots', base_name, num2str(k), ext];
  end
end

fprintf('\n');
fprintf('Number of sensors: %d \n', length(sensor.x));
fprintf('p = %d \n', param.p);
fprintf('beta_0 = %.2f \n', param.beta_0);
fprintf('Number of samples: %d \n', Knoise);
fprintf('Number of bad events: %d \n', bad_count);


if quantitative == 1
    E_dist_bound = criterion.value(sensor, sensitivities) / param.p;
    fprintf('E_bound:                [%e]\n',  E_dist_bound);

    E_dist_W_dz = mean(dist_W_dz);
    dev_E_dist_W_dz = sqrt(var(dist_W_dz) / Knoise);
    fprintf('dz(W)_squared:          [%e, %e, %e]\n',  (E_dist_W_dz + dev_E_dist_W_dz * [-3,0,3]));

    E_dist_dz = mean(dist_dz);
    dev_E_dist_dz = sqrt(var(dist_dz) / Knoise);
    fprintf('dz_squared:             [%e, %e, %e]\n',  (E_dist_dz + dev_E_dist_dz * [-3,0,3]));

    E_dist_z_hat = mean(dist_z_hat);
    dev_E_dist_z_hat = sqrt(var(dist_z_hat) / Knoise);
    fprintf('distance_z_hat_squared: [%e, %e, %e]\n',  (E_dist_z_hat + dev_E_dist_z_hat * [-3,0,3]));

    E_dist_z_bar = mean(dist_z_bar);
    dev_E_dist_z_bar = sqrt(var(dist_z_bar) / Knoise);
    fprintf('distance_z_bar_squared: [%e, %e, %e]\n',  (E_dist_z_bar + dev_E_dist_z_bar * [-3,0,3]));
end