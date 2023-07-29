figure('Renderer', 'painters', 'Position', [1 1 1000 490])
% subplot(2,3,[1 4])
% example_1A_plot(6)

subplot(2, 3, [2 5])
example_data(20)

% subplot(2, 3, [3 6])
% example_1A_plot(11)

fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperSize=[10 17];
fig.Position = [1 1 1000 490];