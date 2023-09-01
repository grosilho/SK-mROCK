% plot domain boundaries

run('../../results/ODECellByCellModel/tmp_geo.m');
figure;
hold on;

domj = [1:size(mdom,1)];
%domj = [1];
%npts= 25 ;
for j=domj
    dom=mdom{j};
    for i=1:size(dom,1) 
        G = dom{i}.G;
        N = dom{i}.N;
        plot(G(:,1),G(:,2));
        %plot(G(1:npts,1),G(1:npts,2),'kx');
        quiver(G(:,1),G(:,2),N(:,1),N(:,2));
    end
end
