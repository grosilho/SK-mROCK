% Plot error estimator data
clear;
clc;
% close all;

nu = 1.1;

% problem = 'OneDimHodgkinHuxley';
% problem = 'InfectiousDiseaseTransmission';
% problem = 'HodgkinHuxley';
% problem = 'RobertsonChemicalSystem';
% problem = 'Brusselator';
problem = 'PDEBrusselator';
% problem = 'Oregonator';
% problem = 'PopulationDynamics';
% problem = 'CUSP';
% problem = 'mReactionDiffusion2DEquations';

folder = ['../dist/Release/GNU-Linux/Tests/' problem '/'];
sol_name = 'sol_rkcsq_257_tol_5';
file_name = [folder sol_name '_errors.m'];
run(file_name);
file_name = [folder sol_name '_eta_dt.m'];
run(file_name);

ind = 1:numel(tnpu);
Iacc = errDnpu<=nu;
indAcc = ind(Iacc);
indRej = ind(logical(1-Iacc));
tAcc = tnpu(indAcc);
tRej = tnpu(indRej);
errDnpuAcc = errDnpu(indAcc);
errDnpuRej = errDnpu(indRej);
ExacterrDnpuAcc = ExacterrDnpu(indAcc);
ExacterrDnpuRej = ExacterrDnpu(indRej);


figure;
subplot(3,1,1);
plot(tn,hn);
hold on;
% plot(tn,hn);
legend('dt');

subplot(3,1,2);
plot(tnpu,errDnpu);
hold on;
plot(tnpu,ExacterrDnpu);
plot(tnpu,errDPPI,'--');
plot(tRej,errDnpuRej,'rx');
legend('$\bar e_n$','$e_n$','Predicted PPI','Rejected','interpreter','latex');
% legend('$\bar e_n$','$e_n$','Rejected','interpreter','latex');
ylim([0,2]);

subplot(3,1,3);


semilogy(tnpu,abs(errDnpu-ExacterrDnpu));
hold on
semilogy(tnpu,abs(errDnpu-errDPPI));
semilogy(tnpu,0*tnpu+1e-1,'k');
legend('$e_n-\bar e_n$','$e_n-errDPPI$','$0.1$','interpreter','latex');
ylim([1e-5,1]);

figure;
plot(tn,hn);
hold on;
plot(tn,eta);

% sgtitle(sol_name);
% 
T = table(tn',tnpu',hn',errDnpu',errDPPI',ExacterrDnpu',eta','variablenames',{'tn','tnpu','hn','errDnpu','errDPPI','ExacterrDnpu','eta'});
TRej = table(tRej',errDnpuRej','variablenames',{'tRej','errDnpuRej'});
writetable(T,['err_' sol_name(5:end) '.csv']);
writetable(TRej,['errRej_' sol_name(5:end) '.csv']);


