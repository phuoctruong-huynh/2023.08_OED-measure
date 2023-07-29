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

%% Example 1A. Reconstruction results with exact data and uniform 9, 11 sensors.

global print_result
print_result = 0;

%% Setup of the problem:

% Number of sources and sensors:
N_sources = 3;
N_gridref = 8;
N_sensors = 2*N_sources + 5;

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

weight = sqrt(abs(q_dagger));
criterion = TV_lin_weighted_criterion(param.beta_0 * param.sig_vec, weight);
kernel = gauss_kernel(param);
%% Calculate estimator error: Choose uniform sensor placement setup:

% Plots:
base_name = 'example1A_11sensor';
ext = '.png';
figure(67)
plot_certificates(kernel, mesh, mu_dagger, sensor)
set(gca,'TickLabelInterpreter','latex', 'FontName', 'Arial', 'Fontsize', 15)
xlabel([num2str(length(sensor.u)), ' sensors'], 'Interpreter','latex', 'FontName', 'Arial')
set(gcf, 'renderer', 'Painters');
% file_name = ['plots/example1A/', base_name, ext];
% print(gcf, '-dpng', '-r300', file_name);
