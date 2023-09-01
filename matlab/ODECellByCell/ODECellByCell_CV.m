clear;
clc;
close all;

si=1.7; % mS/cm
L=125*1e-4; % cm

si2 = si/10; % mS/mm
L2 = L*10; % mm

filepath = '/home/giacomo/Dropbox/Applicazioni/Overleaf/BidomainBEM/images/dep_cv_on_kappa/CV_vs_kappa.csv';
T=readtable(filepath);

T.kappa2 = T.kappa/100; % mS/mm^2
maxCV = max(T.cv); % mm/ms

De = @(k) si*k./(k+si/L); % mS/cm
CV = @(k) maxCV*sqrt(De(k))/sqrt(si); % cm/ms

De2 = @(k) si2*k./(k+si2/L2); % mS/mm
CV2 = @(k) sqrt(De2(k)); %mm/ms

kappas = 10.^linspace(-6,5,20);

semilogx(T.kappa,T.cv);
hold on;
% semilogx(T.kappa2,T.cv);
semilogx(kappas,CV(kappas));
% semilogx(T.kappa2,maxCV*CV2(T.kappa2)/max(CV2(T.kappa2)) );
legend('experimental','formula');