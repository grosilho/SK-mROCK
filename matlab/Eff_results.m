% efficiency results 

clear;
clc;
% close all;

dtmax = 0.25;
n = 10;

k_l = 0:n;
dt_l = dtmax./2.^k_l;

sol_name = ["mRKC","RKC1","IE","RKC2"];%,"ROCK2"];
% sol_name = ["mRKC","IE"];
% sol_name = ["mRKC_ms","mRKC"];
add_descr = '_N_100';

folder = '../dist/Release_Lavoro_GNU/GNU-Linux/Tests/IntegroDifferentialEquation/';

refsol_name = 'refsol_RK4_N_100_dt_1e-5.csv';
refsol_name = [folder refsol_name];

T = readtable(refsol_name);
ref_y = T.y;
N = numel(ref_y);
h = 1/N;

err = zeros(numel(dt_l),numel(sol_name));
time = zeros(numel(dt_l),numel(sol_name));
for s=1:numel(sol_name)
for i=1:numel(dt_l)
    file = [folder sol_name{s} add_descr '_dt_' num2str(dt_l(i)) '.csv'];
    T = readtable(file);
    y = T.y;
    err(i,s) = sqrt(h*sum((y-ref_y).^2));
    file = [folder sol_name{s} add_descr '_dt_' num2str(dt_l(i)) '_statistics.csv'];
    T = readtable(file);
    time(i,s) = T.cpu_time;
end
end

figure;
subplot(2,1,1);
for s=1:numel(sol_name)
    loglog(dt_l,err(:,s),'DisplayName',sol_name{s});
    hold on;
end
loglog(dt_l,dt_l,'k--');
legend show
title("Convergence");

subplot(2,1,2);
for s=1:numel(sol_name)
    loglog(err(:,s),time(:,s),'DisplayName',sol_name{s});
    hold on;
end
legend show
title("Efficiency");

for s=1:numel(sol_name)
    T = table(dt_l',err(:,s),time(:,s),'variablenames',{'dt','err','cpu'});
    writetable(T,[sol_name{s} add_descr '.csv']);
end



