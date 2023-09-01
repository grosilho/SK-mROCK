clear;
%clc;
%close all;

folder = '../../results/ODECellByCellModel/';
folder = '../../results_shared/ODECellByCellModel/';
folder = '../../results_shared/ODECellByCellModel/NUOVI/';
subfolder = 'CVl_C_amp_smooth/';

amp_list = readtable([folder subfolder 'amp_list.csv']);

first_amp = 1;
last_amp= amp_list.n(end);
amp_list = amp_list(first_amp:last_amp,:);
freq = 2;

cl = 1.0e-2;
cw = 0.10e-2;
nc = 30;
np = 1;
pa = zeros(np,2);
pb = zeros(np,2);
for i=1:np
    pa(i,:) = [(7.5+i)*cl,0];
    pb(i,:) = [(17.5+i)*cl,0];
end
V_th = -20;

CV = zeros(size(amp_list,1),1);
for ampi=1:size(amp_list,1)
    solname = ['freq_' num2str(freq) '_amp_' num2str(amp_list.n(ampi))];
    CV(ampi) = get_CV([folder subfolder],solname,pa,pb,V_th,1);
    fprintf(solname);
    fprintf(' CV  = %f mm/ms (m/s)\n',CV(ampi));
end


figure;
plot(amp_list.amp,CV);
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 14;
xlabel('$a$','interpreter','latex')
ylabel('$CV$','interpreter','latex')

T=table(amp_list.amp,amp_list.amp/2,CV,'VariableNames',{'amp','trueamp','CV'});
writetable(T,'CV_vs_amp.csv');


