% Plot solutions Van Der Pol
clear;
clc;
% close all;

folder = '../results/ODEIonicModel/';
solname = 'RKC1_dt_2e_2_F';

n_states = 13; % 8 for LuoRudy, 13 for Fox, 21 for Courtemanche, 41 for Ohara Rudy
fileID = fopen([folder solname '_evolution.bin']);
y = fread(fileID,'double');
y = reshape(y,[1+n_states,numel(y)/(1+n_states)]);
t = y(1,:);
y = y(2:end,:);


plot_gating = 0 ;
t_end = 2000;
t_end = min(t_end,t(end));
tol = 1e-6;
last = find(t>=t_end-tol,1);

t = t(1:last);
V = y(1,1:last);
z = y(2:end,1:last);
n_gating_vars = size(z,1);

figure;


if plot_gating
    subplot(1+n_gating_vars,1,1);
    plot(t,V);
    for i=1:n_gating_vars
        subplot(1+n_gating_vars,1,1+i);
        plot(t,z(i,:));
    end
else
    plot(t,V);
end

% reduce precision to 1ms
n=round(1/t(2));
t = t(1:n:end);
V = V(1:n:end);

T = table(t',V','variablenames',{'t','V'});
writetable(T,'cell.csv');


