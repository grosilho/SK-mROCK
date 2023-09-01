% Show convergence results
clc;
clear;
close all;

integrator = 'mRKC';
solname = integrator;
refsol = 'refsol';
dtmax = 1;
n = 10;

k_l = 0:n;
dt_l = dtmax./2.^k_l;

problem = 'RobertsonChemicalSystem';
results_path = ['../dist/Release/GNU-Linux/Tests/' problem '/'];

params = [results_path refsol '.csv'];
T1 = readtable(params);
yref = T1.y;

err = zeros(numel(k_l),1);
s = zeros(numel(k_l),1);
m = zeros(numel(k_l),1);
eta = zeros(numel(k_l),1);
for k = 1:numel(k_l)
    
    dt = dt_l(k);
    dt_str = num2str(dt);
    
    paramsT = [results_path solname '_dt_' num2str(dt) '_params.csv'];
    solT = [results_path solname '_dt_' num2str(dt) '.csv'];
    params = readtable(paramsT);
    sol = readtable(solT);
    
    err(k) = max(abs(sol.y-yref));
    s(k) = mean(params.s);
    
    if(strcmp(integrator,'mRKC'))
        m(k) = mean(params.m);
        eta(k) = mean(params.eta);
        if(k==1)
            rhoS=params.rhoS;
            rhoF=params.rhoF;
            st = params.s;
            mt = params.m;
            t=params.t;
        end
    else
        if(k==1)
            rho=params.rho;
            st = params.s;
            t=params.t;
        end
    end
end

figure;
loglog(dt_l,err);
hold on;
loglog(dt_l,dt_l,'k--');
legend('err','dt');

% figure;
% plot(dt_l,rhoF);
% hold on;
% plot(dt_l,rhoS);
% legend('rhoF','rhoS');

figure;
plot(dt_l,s);
hold on;
plot(dt_l,m);
legend('s','m');

figure;
loglog(dt_l,eta);
hold on;
loglog(dt_l,dt_l);
legend('eta','dt');

% figure;
% plot(t,rhoF);
% hold on;
% plot(t,rhoS);
% plot(t,rho);
% legend('rhoF','rhoS','rho');

if(strcmp(integrator,'RKC1'))
    T = table(dt_l',err,s,'variablenames',{'dt','err','s'});
%     writetable(T,['res_' solname '.csv']);
    T=table(t,rho,st,'variablenames',{'t','rho','s'});
    writetable(T,['rho_sm__' solname '.csv']);
else
    T = table(dt_l',err,s,m,eta,'variablenames',{'dt','err','s','m','eta'});
%     writetable(T,['res_' solname '.csv']);
    T=table(t,rhoF,rhoS,st,mt,'variablenames',{'t','rhoF','rhoS','s','m'});
    writetable(T,['rho_sm_' solname '.csv']);
end







