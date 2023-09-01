% Plot solutions Van Der Pol
clear;
clc;
% close all;

iter = 4;
n_threads = 10;

folder = '../../Results/Tests/OneDimHodgkinHuxley/';
refsol_name = 'RKL1';
sol_name = 'RKU2_RKU2';
ref_file_name = [folder refsol_name '_evolution.m'];

run(ref_file_name);

V_ref = y(:,1);
n_ref = y(:,2);
m_ref = y(:,3);
h_ref = y(:,4);
t_ref = t;
clear y t;

% plot first coarse guess
figure;
file_name = [folder sol_name '_out_sol_iter_' num2str(0) '_evolution.m'];
run(file_name);
V_outer = y(:,1);
n_outer = y(:,2);
m_outer = y(:,3);
h_outer = y(:,4);
t_outer = t;
clear y t;

subplot(iter+1,4,1);
plot(t_ref,V_ref,'k');
hold on;
plot(t_outer,V_outer,'ro');
subplot(iter+1,4,2);
plot(t_ref,n_ref,'k');
hold on;
plot(t_outer,n_outer,'ro');
subplot(iter+1,4,3);
plot(t_ref,m_ref,'k');
hold on;
plot(t_outer,m_outer,'ro');
subplot(iter+1,4,4);
plot(t_ref,h_ref,'k');
hold on;
plot(t_outer,h_outer,'ro');

for k=1:iter
    file_name = [folder sol_name '_out_sol_iter_' num2str(k) '_evolution.m'];
    run(file_name);
    V_outer = y(:,1);
    n_outer = y(:,2);
    m_outer = y(:,3);
    h_outer = y(:,4);
    t_outer = t;
    clear y t;

    subplot(iter+1,4,4*k+1);
    plot(t_ref,V_ref,'k');
    hold on;
    plot(t_outer,V_outer,'ro');
    subplot(iter+1,4,4*k+2);
    plot(t_ref,n_ref,'k');
    hold on;
    plot(t_outer,n_outer,'ro');
    subplot(iter+1,4,4*k+3);
    plot(t_ref,m_ref,'k');
    hold on;
    plot(t_outer,m_outer,'ro');
    subplot(iter+1,4,4*k+4);
    plot(t_ref,h_ref,'k');
    hold on;
    plot(t_outer,h_outer,'ro');

    for n=1:n_threads
        file_name = [folder sol_name '_in_sol_iter_' num2str(k) '_thread_' num2str(n-1) '_evolution.m'];
        run(file_name);
        V_inner = y(:,1);
        n_inner = y(:,2);
        m_inner = y(:,3);
        h_inner = y(:,4);
        t_inner = t;
        clear y t;

        subplot(iter+1,4,4*k+1);
        plot(t_inner,V_inner,'b.');
        subplot(iter+1,4,4*k+2);
        plot(t_inner,n_inner,'b.');
        subplot(iter+1,4,4*k+3);
        plot(t_inner,m_inner,'b.');
        subplot(iter+1,4,4*k+4);
        plot(t_inner,h_inner,'b.');

    end
end



