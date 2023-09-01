clear;
%clc;
%close all;

% folder = '../../results/ODECellByCellModel/CVl_VS_cells_size/';
folder = '../../results_shared/ODECellByCellModel/CVl_VS_cells_size/';
subfolder_list = {'cl/','cw/','clcw/'};

nc = 30;
pa = @(cl) [8.5*cl,0;9.5*cl,0;10.5*cl,0;11.5*cl,0;12.5*cl,0];
pb = @(cl) [18.5*cl,0;19.5*cl,0;20.5*cl,0;21.5*cl,0;22.5*cl,0];
V_th = -20;

for sf_ind=1:size(subfolder_list,2)
    subfolder = subfolder_list{sf_ind};
    cl_list = readtable([folder subfolder 'cl_list.csv']);
    cw_list = readtable([folder subfolder 'cw_list.csv']);
    
    first_clcw = 1;
    last_clcw = cl_list.n(end)-1;
    cl_list = cl_list(first_clcw:last_clcw,:);
    cw_list = cw_list(first_clcw:last_clcw,:);       
    
    CV = zeros(size(cl_list,1),1);
    for i=1:size(cl_list,1)
        solname = [subfolder(1:end-1) '_' num2str(cl_list.n(i))];
        CV(i)=get_CV([folder subfolder],solname,pa(cl_list.cl(i)),pb(cl_list.cl(i)),V_th,0);
        fprintf(solname);
        fprintf(' cl = %f, cw = %f, CV  = %f mm/ms (m/s)\n',cl_list.cl(i),cw_list.cw(i),CV(i));
    end
    fprintf('\n');
    
    subplot(3,1,sf_ind);
    plot(CV);
    
    T=table(cl_list.cl,cl_list.cl*1e4,cw_list.cw,cw_list.cw*1e4,cl_list.cl.*cw_list.cw,cl_list.cl.*cw_list.cw*1e8,CV,...
        'VariableNames',{'cl_cm','cl_um','cw_cm','cw_um','area_cm','area_um','CV'});
    writetable(T,['CV_vs_' subfolder(1:end-1) '.csv']);    
end

