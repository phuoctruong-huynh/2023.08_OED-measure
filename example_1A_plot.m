%% Example 1A. Reconstruction results with exact data and uniform 9, 11 sensors.
function result = example_1A_plot(N_sensors)
global print_result
print_result = 0;

%% Setup of the problem:

% Number of sources and sensors:
N_sources = 3;
N_gridref = 8;
% N_sensors = 2*N_sources + 3;

% Grid for solve_TV:
x_h = linspace(-1, 1, (2^N_gridref + 1)*N_sensors)'; 
mesh = struct('points', x_h);

% Initialize reference measure:
y_dagger = [-.7, -.3, .3]';
q_dagger = [.4, .3, -.2]';
mu_dagger = struct('x', y_dagger, 'u', q_dagger);

% Initialize sensor placement:
xx = linspace(-1, 1, N_sensors);
% xx = [-.8, -.6, -.4, -.1, .1, .4];
uu = 1/length(xx) * ones(length(xx), 1);
% Artificial init_sensor:

sensor = struct('x', xx, 'u', uu);

% Parameters: s^2, sign_vector, beta_0, m.
param = struct();
T = 1/2*(0.2).^2;
sigma = sqrt(2*T);
param.s2 = sigma.^2;
param.sig_vec = [sign(q_dagger); zeros(N_sources, 1)];
param.beta_0 = 2;

weight = sqrt(abs(q_dagger));
kernel = gauss_kernel(param);
%% Calculate estimator error: Choose uniform sensor placement setup:

% Plots:
plot_certificates(kernel, mesh, mu_dagger, sensor)
set(gca,'TickLabelInterpreter','latex', 'FontName', 'Arial', 'Fontsize', 18)
xlabel([num2str(length(sensor.u)), ' sensors'], 'Interpreter','latex', 'FontName', 'Arial')
set(gcf, 'renderer', 'Painters');
end