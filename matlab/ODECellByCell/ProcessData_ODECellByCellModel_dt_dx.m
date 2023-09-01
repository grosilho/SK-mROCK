function ProcessData_ODECellByCellModel_dt_dx()

folder = '../../results/ODECellByCellModel/';
folder = '../../results_cluster/ODECellByCellModel/nx=30/';
folder = '../../results_shared/ODECellByCellModel/';
subfolder = 'CVl_C_dtdx/';
filename_base = '';

dt_list = readtable([folder subfolder 'dt_list.csv']);
dG_list = readtable([folder subfolder 'dG_list.csv']);

first_dt = 1;
last_dt = dt_list.n(end);
first_dG = 1;
last_dG = dG_list.n(end);
dt_list = dt_list(first_dt:last_dt,:);
dG_list = dG_list(first_dG:last_dG,:);

cl = 1.0e-2;
np = 5;
pa = zeros(np,2);
pb = zeros(np,2);
for i=1:np
    pa(i,:) = [(7.5+i)*cl,0];
    pb(i,:) = [(17.5+i)*cl,0];
end
V_th = -20;

CV_ref = get_CV([folder subfolder],[filename_base 'dt_ref_dG_ref'],pa,pb,V_th,false);
fprintf([filename_base 'dt_ref_dG_ref']);
fprintf(' CV  = %f mm/ms (m/s)\n',CV_ref);
%CV_ref=0;

CV = zeros(size(dt_list,1),size(dG_list,1));
for dti=1:size(dt_list,1)
for dGj=1:size(dG_list,1)
    solname = [filename_base 'dt_' num2str(dt_list.n(dti)) '_dG_' num2str(dG_list.n(dGj))];
    %CV(dti,dGj) = get_CV_old(folder,subfolder,solname,n_states,ps,V_th);
    CV(dti,dGj) = get_CV([folder subfolder],solname,pa,pb,V_th,false);
    fprintf(solname);
    fprintf(' CV  = %f mm/ms (m/s)\n',CV(dti,dGj));
end
end

if CV_ref==0
    CV_ref = CV(1,1);
end

CV_err = (CV-CV_ref)/CV_ref;

% CV_err = abs(CV_err);

dG_list.dG = dG_list.dG*1e4; %from cm to um

[dG,dt]=meshgrid(dG_list.dG,dt_list.dt);
% levels = linspace(-0.25,0,11);
%levels = [-0.10,-0.05,-0.04,-0.03,-0.02,-0.01,0];
figure;
%contourf(log2(dG),log2(dt),CV_err,levels);
contourf(log2(dG),log2(dt),CV_err);
cb = colorbar;
cb.Label.String = '$Err=(CV-CV^*)/CV$';
cb.Label.Interpreter = 'latex';
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 14;
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 14;
xlabel('$\Delta x\ [\mu m]$','interpreter','latex')
ylabel('$\Delta t\ [ms]$','interpreter','latex')
xticks(log2(dG_list.dG));
yticks(log2(dt_list.dt));
ylabels = cell(1,size(dt_list,1));
xlabels = cell(1,size(dG_list,1));
for i=1:size(dt_list,1)
    ylabels{i} = ['$' num2str(dt_list.dt(i)) '$'];
end
for i=1:size(dG_list,1)
    xlabels{i} = ['$' num2str(dG_list.dG(i)) '$'];
end
xticklabels(xlabels);
yticklabels(ylabels);


CV_vs_dt=array2table([dt_list.dt,CV,CV_err,abs(CV_err)]);
col_CV = cell(1,size(dG_list,1));
col_CVerr = cell(1,size(dG_list,1));
col_absCVerr = cell(1,size(dG_list,1));
for i=1:size(dG_list,1)
    col_CV{i} = ['CV_dG_' num2str(dG_list.n(i))];
    col_CVerr{i} = ['CVerr_dG_' num2str(dG_list.n(i))];
    col_absCVerr{i} = ['absCVerr_dG_' num2str(dG_list.n(i))];
end
CV_vs_dt.Properties.VariableNames = [{'dt'} col_CV col_CVerr col_absCVerr];
writetable(CV_vs_dt,'CV_vs_dt.csv');

CV_vs_dx=array2table([dG_list.dG,CV',CV_err',abs(CV_err')]);
col_CV = cell(1,size(dt_list,1));
col_CVerr = cell(1,size(dt_list,1));
col_absCVerr = cell(1,size(dt_list,1));
for i=1:size(dt_list,1)
    col_CV{i} = ['CV_dt_' num2str(dt_list.n(i))];
    col_CVerr{i} = ['CVerr_dt_' num2str(dt_list.n(i))];
    col_absCVerr{i} = ['absCVerr_dt_' num2str(dt_list.n(i))];
end
CV_vs_dx.Properties.VariableNames = [{'dx'} col_CV col_CVerr col_absCVerr];
writetable(CV_vs_dx,'CV_vs_dx.csv');    

dtT=dt';
dGT=dG';
CV_errT=CV_err';

CV_vs_dt_vs_dx=array2table([log2(dGT(:)),log2(dtT(:)),CV_errT(:)]);
writetable(CV_vs_dt_vs_dx,'CV_vs_dt_vs_dx.csv');

CV_vs_dt_vs_dx_nll=array2table([dGT(:),dtT(:),CV_errT(:)]);
writetable(CV_vs_dt_vs_dx_nll,'CV_vs_dt_vs_dx_notloglog.csv');

end
