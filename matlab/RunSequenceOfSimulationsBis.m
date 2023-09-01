% Running convergence tests
clc;
clear;

integrator = ["mRKC","IE","RKC1","RKC2"];%,"ROCK2"];
% integrator = ["RKC1","RKC2","ROCK2"];
% integrator = ["ROCK2"];
% integrator = 'RKC1';
% integrator = 'RKC2';
% integrator = 'mRKC';
integrator = ['IE'];

dtmax = 0.25;
n = 10;

k_l = 0:n;
dt_l = dtmax./2.^k_l;

dta = 'false';
eq_str = ' -ode';
ofreq = 0;
ntest = 18;
verb = 'false';
add_descr = '_N_100_T_1';

prog_path = '../dist/Release_Lavoro_GNU/GNU-Linux';
prog_name = './multirate_integrators';


for s = 1:numel(integrator)
for k = 1:numel(k_l)
    
    dt = dt_l(k);
    
    dta_str = [' -dta ' dta];
    dt_str = [' -dt ' num2str(dt)];
    rk_str = [' -rk ' integrator{s}];
    ofreq_str = [' -ofreq ' num2str(ofreq)];
    ntest_str = [' -ntest ' num2str(ntest)];
    ofile_str = [' -ofile ' integrator{s} add_descr '_dt_' num2str(dt)];
    verb_str = [' -verb ' verb];
    command = [prog_name eq_str dta_str dt_str rk_str ofreq_str ntest_str ofile_str verb_str];
    fprintf([command '\n']);
    system(['cd ' prog_path ' && ' command]); 
end
end