clear;
%clc;
%close all;

segment = 'flat';
folder = ['../../results_cluster/ODECellByCellModel/CVltd_C_' segment '_randomly_deactivated_gj_nx_20_ny_20/'];
%folder = ['../../results_cluster/ODECellByCellModel/CVltd_C_' segment '_randomly_deactivated_gj_nx_20_ny_20_stim_cell_middle/'];

N = readtable([folder 'p_0/N.csv']);
N = table2array(N);
minN = 46;
maxN = 50;

cl = 1.0e-2;
cw = 0.10e-2;
ncx = 20;
ncy = 20;

write_files = 1;

% np_l=1;
% pa_l = zeros(np_l,2);
% pb_l = zeros(np_l,2);
% for n=1:np_l
%     pa_l(n,:) = [(-0.5+n)*cl,0];
%     pb_l(n,:) = [(8.5+n)*cl,0];
% end
% np_t=1;
% pa_t = zeros(np_t,2);
% pb_t = zeros(np_t,2);
% for n=1:np_t
%     pa_t(n,:) = [0,(0.5+n)*cw];
%     pb_t(n,:) = [0,(13.5+n)*cw];
% end
% pa_d = [0,0];
% pb_d = [ncx*cl,ncy*cw];

pa_l = [0,0];
pb_l = [ncx*cl,0];
pa_t = [0,0];
pb_t = [0,ncy*cw];  
pa_d = [0,0];
pb_d = [ncx*cl,ncy*cw];    

V_th = -20; 

CVl = zeros(maxN-minN+1,11);
CVt = zeros(maxN-minN+1,11);
CVd = zeros(maxN-minN+1,11);

col_names = cell(1,11);

tic;
for pi=1:11
    p = 10*(pi-1);
    subfolder = ['p_' num2str(p) '/'];    
    col_names{1,pi} = num2str(p);

    for n=minN:maxN
        solname = ['n_' num2str(n)];

        [CVl(n-minN+1,pi),transmembrane_G,y,t]=get_CV([folder subfolder],solname,pa_l,pb_l,V_th,0);
        CVt(n-minN+1,pi)=get_CV([folder subfolder],solname,pa_t,pb_t,V_th,0,transmembrane_G,y,t);
        CVd(n-minN+1,pi)=get_CV([folder subfolder],solname,pa_d,pb_d,V_th,0,transmembrane_G,y,t);

        fprintf('p = %i, n = %f, CVl  = %f mm/ms (m/s)\n',p,n,CVl(n-minN+1,pi));
        fprintf('p = %i, n = %f, CVt  = %f mm/ms (m/s)\n',p,n,CVt(n-minN+1,pi));
        fprintf('p = %i, n = %f, CVd  = %f mm/ms (m/s)\n',p,n,CVd(n-minN+1,pi));
    end
    
%     figure;
%     sgtitle(['p = ' num2str(p)])
%     subplot(1,3,1);
%     histogram(CVl(:,pi));
%     title('L');
%     subplot(1,3,2);
%     histogram(CVt(:,pi));
%     title('T');
%     subplot(1,3,3);
%     histogram(CVd(:,pi));
%     title('D');
end
toc;

T_CVl = array2table(CVl);
T_CVl.Properties.VariableNames = col_names;

T_CVt =array2table(CVt);
T_CVt.Properties.VariableNames = col_names;

T_CVd =array2table(CVd);
T_CVd.Properties.VariableNames = col_names;

mul = mean(CVl,1)';
sigmal = std(CVl,1)';
CVl_stats=table((0:10:100)',mul,sigmal,'VariableNames',{'p','mul','sigmal'});

mut = mean(CVt,1)';
sigmat = std(CVt,1)';
CVt_stats=table((0:10:100)',mut,sigmat,'VariableNames',{'p','mut','sigmat'});

mud = mean(CVd,1)';
sigmad = std(CVd,1)';
CVd_stats=table((0:10:100)',mud,sigmad,'VariableNames',{'p','mud','sigmad'});

if write_files==1
    writetable(T_CVl,['CVl_' segment '_randomly_deactivated.csv'],'WriteMode','Append','WriteVariableNames',false);
    writetable(T_CVt,['CVt_' segment '_randomly_deactivated.csv'],'WriteMode','Append','WriteVariableNames',false);
    writetable(T_CVd,['CVd_' segment '_randomly_deactivated.csv'],'WriteMode','Append','WriteVariableNames',false);
    writetable(CVl_stats,[segment '_CVl_stats.csv']);
    writetable(CVt_stats,[segment '_CVt_stats.csv']);
    writetable(CVd_stats,[segment '_CVd_stats.csv']);
end
