% Plot solutions 
clear;
clc;
% close all;

iter = 4;
n_threads = 10;

folder = '../../Results/Tests/ScalarNonStiffNonLinearTest/';
sol_name = 'RKC1';
refsol_name = 'RKC1';
ref_file_name = [folder refsol_name '_evolution.m'];

run(ref_file_name);

y_ref = y;
t_ref = t;
clear y t;

% plot first coarse guess
figure;
file_name = [folder refsol_name '_out_sol_iter_' num2str(0) '_evolution.m'];
run(file_name);
y_outer = y;
t_outer = t;
clear y t;

subplot(iter+1,1,1);
plot(t_ref,y_ref,'k');
hold on;
plot(t_outer,y_outer,'ro');

for k=1:iter
    file_name = [folder refsol_name '_out_sol_iter_' num2str(k) '_evolution.m'];
    run(file_name);
    y_outer = y;
    t_outer = t;
    clear y t;

    subplot(iter+1,1,k+1);
    plot(t_ref,y_ref,'k');
    hold on;
    plot(t_outer,y_outer,'ro');

    for n=1:n_threads
        file_name = [folder refsol_name '_in_sol_iter_' num2str(k) '_thread_' num2str(n-1) '_evolution.m'];
        run(file_name);
        y_inner = y;
        t_inner = t;
        clear y t;

        subplot(iter+1,1,k+1);
        plot(t_inner,y_inner,'b.');
    end
end