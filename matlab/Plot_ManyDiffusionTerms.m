% Plot solutions Van Der Pol
clear;
clc;
% close all;

folder = '../results/ManyDiffusionTerms/';
solname = 'sol';

n_states = 7;
fileID = fopen([folder solname '_evolution.bin']);
y = fread(fileID,'double');
y = reshape(y,[1+n_states,numel(y)/(1+n_states)]);
t = y(1,:);
y = y(2:end,:);

t_end = 2000;
t_end = min(t_end,t(end));
tol = 1e-6;
last = find(t>=t_end-tol,1);

t = t(1:last);
y = y(:,1:last);

figure;
for i=1:n_states
    subplot(4,2,i);
    plot(t,y(i,:));
end



