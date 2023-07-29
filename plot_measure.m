function h = plot_measure(x_h, measure, marker, color)

  if ~isstruct(measure)
    %% assume this is a discretized mesure
    %% get the parameters of the measure
    u_h = measure;
    dx = min(diff(x_h));
    supp = (u_h ~= 0);
    ux = x_h(supp);
    uu = u_h(supp);
    
    %% postprocess_measure  
    measure = postprocess_measure(ux, uu, 1.5*dx);
  end

  %% plot measure
  xi = [measure.x(:)'; measure.x(:)'];
  yi = [zeros(size(measure.u(:)')); measure.u(:)'];
  plot(xi, yi, 'Color', color, 'LineWidth', 1);
  had_hold = ishold();
  hold on;
  plot([x_h(1),x_h(end)], [0,0], 'k:', 'LineWidth', 1);
  h = plot(measure.x, measure.u, marker, 'Color', color, 'LineWidth', 1, 'MarkerSize', 3);
  %legend(name)
  set(gca, 'FontSize', 12);
  if ~had_hold
    hold off;
  end
    
end
