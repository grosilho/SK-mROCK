clear;
%clc;
%close all;

segment = 'squaredwave_op';
% folder = ['../../results/ODECellByCellModel/CVt_C_' segment '_random_gj/'];
folder = ['../../results_shared/ODECellByCellModel/NUOVI/CVt_C_' segment '_random_gj/'];

N = readtable([folder 'N.csv']);
N = table2array(N);

cl = 1.0e-2;
cw = 0.10e-2;
nc = 30;
np=5;
pa = zeros(np,2);
pb = zeros(np,2);
for n=1:np
    pa(n,:) = [0,(8.5+n)*cw];
    pb(n,:) = [0,(18.5+n)*cw];
end
V_th = -20; 

CV = zeros(N,1);
for n=1:N
    solname = ['n_' num2str(n)];
    CV(n)=get_CV(folder,solname,pa,pb,V_th);
    fprintf(solname);
    fprintf(' n = %f, CV  = %f mm/ms (m/s)\n',n,CV(n));
end

figure;
histogram(CV);

T_CV=table(CV,'VariableNames',{'CV'});
writetable(T_CV,[segment '_CV.csv']);

mu = mean(CV);
sigma = std(CV);
T_stats=table(mu,sigma,'VariableNames',{'mu','sigma'});
writetable(T_stats,[segment '_stats.csv']);


