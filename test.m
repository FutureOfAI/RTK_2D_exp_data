% SamplePlotFreq = 1;
% Spin = 120;
% posPlot = [xpm_Nh(100:n),ypm_Nh(100:n),zeros(n-99,1)];
% RT_PLOT(posPlot, coordinate, Rm_data_valid_inx, ...
%             'SamplePlotFreq', SamplePlotFreq, 'Trail', 'All', ...
%             'Position', [9 39 1280 768], 'View', [(100:(Spin/(length(posPlot)-1)):(100+Spin))', 10*ones(length(posPlot), 1)], ...
%             'AxisLength', 0.1, 'ShowArrowHead', false, ...
%             'Xlabel', 'X (m)', 'Ylabel', 'Y (m)', 'Zlabel', 'Z (m)', 'ShowLegend', false, ...
%             'CreateAVI', false, 'AVIfileNameEnum', false, 'AVIfps', ((1/samplePeriod) / SamplePlotFreq));

test_point = [1,3.65,-4.025,1.38;
              2,4.65,-4.025,1.38;
              3,11.14,-4.025,1.38;
              4,12.14,-4.025,1.38;
              5,21.425,-4.025,1.38;
              6,22.425,-4.025,1.38];
for i=1:6
    test_dis(i,1) = norm(coordinate(i,2:4) - test_point(6,2:4));
end

for i=1:6
    mean_meare(i,1) = mean(pos(i,400:3000));
end
