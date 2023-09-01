% Show convergence results
clc;
clear;

integrator1 = 'RKC1';
integrator2 = 'mRKC';
refsol = 'refsol';
dtmax = 1;
n = 8;

k_l = 0:n;
dt_l = dtmax./2.^k_l;

problem = 'RobertsonChemicalSystem';
results_path = ['../dist/Release/GNU-Linux/Tests/' problem '/'];

file1 = [results_path refsol '.csv'];
T1 = readtable(file1);
yref = T1.y;

err1 = zeros(numel(k_l),1);
err2 = err1;
for k = 1:numel(k_l)
    
    dt = dt_l(k);
    dt_str = num2str(dt);
    
    file1 = [results_path integrator1 '_dt_' num2str(dt) '.csv'];
    file2 = [results_path integrator2 '_dt_' num2str(dt) '.csv'];
    T1 = readtable(file1);
    T2 = readtable(file2);
    err1(k) = max(abs(T1.y-yref));
    err2(k) = max(abs(T2.y-yref));
end

loglog(dt_l,err1);
hold on;
loglog(dt_l,err2);
loglog(dt_l,dt_l,'k--');
legend('RKC','mRKC','dt');

T = table(dt_l',err1,err1,'variablenames',{'dt','err_rkc1','err_rkcsq'});
% writetable(T,'conv_test.csv');