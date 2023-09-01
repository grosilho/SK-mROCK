% plot domain boundaries
clear;

run('../../results/ODECellByCellModel/tmp_geo.m');
%figure('Position', [10 10 1000 100]);
figure;
hold on;

dom=mdom{2};
G = dom{1}.G;
plot(G(:,1),G(:,2));

dom=mdom{3};
G = dom{1}.G;
plot(G(:,1),G(:,2));