clear;
%clc;
%close all;

folder = '../../results/ODECellByCellModel/';
folder = '../../results_cluster/ODECellByCellModel/nx=30/';
folder = '../../results_shared/ODECellByCellModel/NUOVI/';
subfolder = 'CVl_C_freq_smooth/';

freq_list = readtable([folder subfolder 'freq_list.csv']);
amp = 1;

first_freq = 1;
last_freq= freq_list.n(end)-1;
freq_list = freq_list(first_freq:last_freq,:);

%n_states = 21; % 8 for LuoRudy, 21 for Courtemanche, 41 for Ohara Rudy
cl = 1.0e-2;
cw = 0.10e-2;
nc = 30;
ps = zeros(2,2);
% longitudinal 
np = 5;
pa = zeros(np,2);
pb = zeros(np,2);
for i=1:np
    pa(i,:) = [(7.5+i)*cl,0];
    pb(i,:) = [(17.5+i)*cl,0];
end
%pa = [9.5*cl,0];
%pb = [19.5*cl,0];
V_th = -20;

CV = zeros(size(freq_list,1),1);
for freqi=1:size(freq_list,1)
    solname = ['freq_' num2str(freq_list.n(freqi)) '_amp_' num2str(amp)];
    CV(freqi) = get_CV([folder subfolder],solname,pa,pb,V_th);
    fprintf(solname);
    fprintf(' CV  = %f mm/ms (m/s)\n',CV(freqi));
end


figure;
plot(freq_list.freq,CV);
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 14;
xlabel('$k$','interpreter','latex')
ylabel('$CV$','interpreter','latex')

T=table(freq_list.freq,CV,'VariableNames',{'freq','CV'});
writetable(T,'CV_vs_freq.csv');


