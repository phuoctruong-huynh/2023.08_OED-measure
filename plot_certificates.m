function plot_certificates(kernel, mesh, mu_dagger, sensor)

q_dagger = mu_dagger.u;
y_dagger = mu_dagger.x;

K_initial = kernel.matrix(sensor.x, mesh.points);
KK_initial = kernel.matrix(sensor.x, y_dagger);

SI_pre_dual_0 = calculate_pre_certificate(sensor, mesh, q_dagger, y_dagger, kernel);
[SI_dual_0, mu_min_0] = calculate_certificate(sensor, mesh, q_dagger, y_dagger, kernel);

h1 = plot(mesh.points, K_initial'*SI_pre_dual_0, 'b', 'LineWidth', 1);
hold on
h2 = plot(mesh.points, K_initial'*SI_dual_0, 'k--', 'LineWidth', 1);
h3 = plot_measure_new(mu_dagger, 'o', 'r', 2, 4, 'r');
h4 = plot_measure(mesh.points, mu_min_0, 's', 'k');
h5 = plot_measure(mesh.points, sensor, '*', 0.5*[1 1 1]);
plot(y_dagger, KK_initial'*SI_pre_dual_0, 'bx', 'LineWidth', 1);
plot(y_dagger, KK_initial'*SI_dual_0, 'kx', 'LineWidth', 1);

plot([-1,1], [ones(2, 1), -ones(2, 1)], 'k:', 'LineWidth', 1);
axis([-1, 1, -2.1, 1.5]);
hold off

lgd = legend([h1, h2, h3, h4, h5], 'sensor pre-certificate', 'sensor certificate', 'reference measure', 'minimum norm solution', 'sensors');
set(lgd, 'Location', 'southwest')
set(lgd,'Interpreter','latex')
lgd.FontSize = 14;
end
