function h = plot_measure_new(measure, marker, color, linewidth, markersize, markercolor)

  xi = [measure.x(:)'; measure.x(:)'];
  yi = [zeros(size(measure.u(:)')); measure.u(:)'];
  plot(xi, yi, 'Color', color, 'LineWidth', linewidth);
  had_hold = ishold();
  hold on;
  plot([-1,1], [0,0], 'k:', 'LineWidth', 1);
  h = plot(measure.x, measure.u, marker, 'Color', color, 'LineWidth', 1, 'MarkerSize', markersize, 'MarkerFaceColor', markercolor);
  %legend(name)
  if ~had_hold
    hold off;
  end
end
