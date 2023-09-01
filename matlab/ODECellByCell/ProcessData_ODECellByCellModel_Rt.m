clear;
%clc;
%close all;

folder = '../../results/ODECellByCellModel/CVt_C_Rt/';

Rt_list = readtable([folder 'Rt_list.csv']);

first_Rt = 4;
last_Rt= Rt_list.n(end);
Rt_list = Rt_list(first_Rt:last_Rt,:);

cl = 1.0e-2;
cw = 0.10e-2;
nc = 30;
pa = [0,9.5*cw];
pb = [0,19.5*cw];  
V_th = -20;

CV = zeros(size(Rt_list,1),1);
for i=1:size(Rt_list,1)
    solname = ['Rt_' num2str(Rt_list.n(i))];
    CV(i)=get_CV(folder,solname,pa,pb,V_th);
    fprintf([folder solname]);
    fprintf(' Rt = %f, CV  = %f mm/ms (m/s)\n',Rt_list.Rt(i),CV(i));
end


figure;
loglog(Rt_list.Rt,CV);
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 14;
xlabel('$R_t$','interpreter','latex')
ylabel('$CV$','interpreter','latex')

T=table(Rt_list.Rt,1./Rt_list.Rt,CV,'VariableNames',{'Rt','k','CV'});
writetable(T,'CV_vs_Rt.csv');


