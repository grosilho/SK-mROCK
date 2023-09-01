% Plot solutions Van Der Pol
clear;
clc;
%close all;

folder = '../../Results/Tests/NeuronCable/';
refsol_name = 'RKC1';
sol_name = 'RKU1_RKU1';
ref_file_name = [folder refsol_name '_evolution'];

iter = 4;
n_threads = 10;

run(ref_file_name);
neqn = numel(y(1,:));
y_ref = y;
t_ref = t;
clear y t;
x_ref = linspace(0,1,neqn);
[X_ref,T_ref]=meshgrid(x_ref,t_ref);

file_name = [folder sol_name '_out_sol_iter_' num2str(0) '_evolution.m'];
run(file_name);
y_outer = y;
t_outer = t;
clear y t;
[X_outer,T_outer]=meshgrid(x_ref,t_outer);

figure;
surf(X_ref,T_ref,y_ref,'FaceColor','r', 'FaceAlpha',1, 'EdgeColor','r');
hold on;
surf(X_outer,T_outer,y_outer,'FaceColor','b', 'FaceAlpha',0.5, 'EdgeColor','b','LineStyle','-');


for k=1:iter
    file_name = [folder sol_name '_out_sol_iter_' num2str(k) '_evolution.m'];
    run(file_name);
    y_outer = y;
    t_outer = t;
    clear y t;
    [X_outer,T_outer]=meshgrid(x_ref,t_outer);

    figure;
    surf(X_ref,T_ref,y_ref,'FaceColor','r', 'FaceAlpha',1, 'EdgeColor','r');
    hold on;
    surf(X_outer,T_outer,y_outer,'FaceColor','b', 'FaceAlpha',0.5, 'EdgeColor','b','LineStyle','-');


    for n=1:n_threads
        file_name = [folder sol_name '_in_sol_iter_' num2str(k) '_thread_' num2str(n-1) '_evolution.m'];
        run(file_name);
        y_inner = y;
        t_inner = t;
        clear y t;
        [X_inner,T_inner]=meshgrid(x_ref,t_inner);

        surf(X_inner,T_inner,y_inner,'FaceColor','g', 'FaceAlpha',0.5, 'EdgeColor','g','LineStyle','-');
    end
end



