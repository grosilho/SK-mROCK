

N = [100,200, 400, 800, 1600,3200, 6400];
time_IE = [0.0050, 0.0177, 0.0783, 0.3567, 2.1333, 18.6667, 146];
time_mRKC = [0.01,0.03,0.1,0.3,1.24,4.8,16];

loglog(N,time_IE,'displayname','IE');
hold on;
loglog(N,time_mRKC,'displayname','mRKC');
loglog(N,N.*N/N(1)^2*time_IE(1),'k--','displayname','N^2');
loglog(N,N.*N.*N/N(1)^3*time_IE(1),'k-','displayname','N^3');
legend show;