function ProcessData_ODECellByCellModel_CVt
% clc;
%close all;

folder = '../../results_shared/ODECellByCellModel/Try_2_CVt_Better/';
folder = '../../results_shared/ODECellByCellModel/CVt/';

writing_folder = '~/Dropbox/Applicazioni/Overleaf/BidomainBEM/images/cv_transversal/';
writing_folder = './';
%writing_folder = ''; % do not write processed data, only plots

segment = 'flat';
%segment = 'wave';
%segment = 'squaredwave';
%segment = 'squaredwave_op';

subfolder1 = ['CVt_C_' segment '/'];
subfolder2 = {'center','side','side_alt'};
subfolder3 = {'length','amplitude','both','p'};

cw = 10e-4;
V_th = -20;

pa = zeros(1,2);
pb = zeros(1,2);
pa(1,:) = [0, 9.5*cw];
pb(1,:) = [0, 19.5*cw];

np=5;
pa = zeros(np,2);
pb = zeros(np,2);
for i=1:np
    pa(i,:) = [0,(8.5+i)*cw];
    pb(i,:) = [0,(18.5+i)*cw];
end

new_figure=0;

amplitude_list = cell(size(subfolder2,2),size(subfolder3,2));
length_list = cell(size(subfolder2,2),size(subfolder3,2));
p_list = cell(size(subfolder2,2),size(subfolder3,2));
for i=1:size(subfolder2,2)
    for j=1:size(subfolder3,2)
        current_folder = [folder subfolder1 subfolder2{i} '/' subfolder3{j} '/'];
        if exist(current_folder, 'dir')
            amplitude_list{i,j} = readtable([current_folder 'amplitude_list.csv']);
            length_list{i,j} = readtable([current_folder 'length_list.csv']);
            p_list{i,j} = readtable([current_folder 'p_list.csv']);
            n_list_elems = max([size(amplitude_list{i,j},1),size(length_list{i,j},1),size(p_list{i,j},1)]);
            amplitude_list{i,j}=resize_table(amplitude_list{i,j},n_list_elems);
            length_list{i,j}=resize_table(length_list{i,j},n_list_elems);
            p_list{i,j}=resize_table(p_list{i,j},n_list_elems);            
            for k=1:n_list_elems
                amp_n = amplitude_list{i,j}.n(k);
                len_n = length_list{i,j}.n(k);
                p_n = p_list{i,j}.n(k);
                filename = ['length_' num2str(len_n) '_amplitude_' num2str(amp_n) '_p_' num2str(p_n) '_statistics.csv'];
                if ~exist([current_folder filename], 'file')
                    fprintf([subfolder1 subfolder2{i} subfolder3{j} filename ' does not exist\n']);
                end
            end
        end
    end
end


CV = cell(size(subfolder2,2),size(subfolder3,2));
area_list = cell(size(subfolder2,2),size(subfolder3,2));

if new_figure
    figure;
    sgtitle(segment); 
else
    figure(1);
end
for i=1:size(subfolder2,2)
    for j=1:size(subfolder3,2)
        current_folder = [folder subfolder1 subfolder2{i} '/' subfolder3{j} '/'];
        if exist(current_folder, 'dir')
            n_list_elems = size(amplitude_list{i,j},1);
            CV{i,j} = zeros(n_list_elems,1);
            area_list{i,j} = zeros(n_list_elems,1);
            for k=1:n_list_elems
                amp_n = amplitude_list{i,j}.n(k);
                len_n = length_list{i,j}.n(k);
                p_n = p_list{i,j}.n(k);
                filename = ['length_' num2str(len_n) '_amplitude_' num2str(amp_n) '_p_' num2str(p_n)];
                fprintf(['For ' subfolder1 subfolder2{i} subfolder3{j} filename ' ']);
                CV{i,j}(k) = get_CV(current_folder,filename,pa,pb,V_th,0);      
                area_list{i,j}(k) = contact_area(length_list{i,j}.length(k),amplitude_list{i,j}.amplitude(k),segment);
                fprintf([' CV = ' num2str(CV{i,j}(k)) '\n']);        
            end
            subplot(size(subfolder2,2),size(subfolder3,2),(i-1)*size(subfolder3,2)+j);
            semilogy(CV{i,j},'displayname',segment(1));
            hold on;
            legend show;
            title([subfolder2{i} subfolder3{j}]);
        end
    end
end

if ~strcmp(writing_folder,'')
    write_subfolder = [writing_folder segment];
    system(['mkdir -p ' write_subfolder]);
    for i=1:size(subfolder2,2)
        for j=1:size(subfolder3,2)     
            if ~isempty(p_list{i,j})
                write_file = [write_subfolder '/' subfolder2{i} '_' subfolder3{j} '.csv'];
                T = table(CV{i,j},...
                    length_list{i,j}.n,length_list{i,j}.length,1e4*length_list{i,j}.length,...
                    amplitude_list{i,j}.n,amplitude_list{i,j}.amplitude,1e4*amplitude_list{i,j}.amplitude,...
                    p_list{i,j}.n,p_list{i,j}.p,...
                    area_list{i,j},1e4*area_list{i,j},...
                    'variablenames',{'CV','len_n','len_len_cm','len_len_um','amp_n','amp_amp_cm','amp_amp_um','p_n','p_p','area_cm','area_um'});
                writetable(T,write_file);
            end
        end
    end
end

end

function [new_tab,stay_same] = resize_table(tab,n,stay_same)
    onesvec = ones(n,1);
    varnames = tab.Properties.VariableNames;
    if size(tab,1)<n        
        tab = tab.Variables;
        new_tab = table(tab(1,1)*onesvec,tab(1,2)*onesvec,...
                        'VariableNames',varnames);
    else
        if nargin==3
            stay_same = [stay_same varnames{2} ', '];
        end
        new_tab=tab;
    end
end

function A = contact_area(l,a,segment)

    if strcmp(segment,'flat')
        A = l;
    elseif strcmp(segment,'wave')
        k=1;
        norm_dg = @(t) sqrt( l*l+(2.*a*k*pi*cos(t*k*2*pi)).^2 );
        A = integral(norm_dg,0,1);
    elseif strcmp(segment,'squaredwave')
        A = l+4*a;
    elseif strcmp(segment,'squaredwave_op')
        A = 4*a;
    end

end

