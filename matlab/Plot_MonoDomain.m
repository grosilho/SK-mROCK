clear;
%clc;
%close all;

solname = 'tmp';
read_gating_vars = false;

folder = '../results/MonoDomain/';
folder = '../results/ODECellByCellModel/';

Lx = 1.25;
Ly = 0.05;
h = 0.01;

Nx = ceil(1+Lx/h);
Ny = ceil(1+Ly/h);
hx = Lx/(Nx-1.);
hy = Ly/(Ny-1.);

n_pts = Nx*Ny;
n_states = 8; % 8 for LuoRudy, 13 for Fox, 21 for Courtemanche, 41 for Ohara Rudy
fileID = fopen([folder solname '_evolution.bin']);
if read_gating_vars
    A = fread(fileID,'double');
    A = reshape(A,[1+n_pts*n_states,numel(A)/(1+n_pts*n_states)]);
    t = A(1,:);
    A = A(2:end,:);
    V = A(1:n_pts,:);
    z = A(n_pts+1:end,:);
else
    precision = [num2str(1+n_pts) '*double'];
    skip = n_pts*(n_states-1)*8;
    A = fread(fileID,precision,skip);
    A = reshape(A,[1+n_pts,numel(A)/(1+n_pts)]);
    t = A(1,:);
    V = A(2:end,:);
end

V = reshape(V,[Nx,Ny,numel(t)]);

x = linspace(0,Lx,Nx);
y = linspace(0,Ly,Ny);
[Y,X] = meshgrid(y,x);

if 0    
    surf(X,Y,V(:,:,21));
end

if 1
    ps = zeros(2,2);
    ps(1,:) = [0.10*Lx,0];
    ps(2,:) = [0.5*Lx,0];
    p = zeros(2,2);
    Vp = zeros(size(ps,1),size(V,3));
    hist_j = zeros(size(ps,1),2);
    figure;
    hold on;
    for i=1:size(ps,1)
        d = (X-ps(i,1)).^2+(Y-ps(i,2)).^2;
        [~,j] = min(d(:));
        [r,c] = ind2sub(size(d),j);
        hist_j(i,:) = [r,c];
        p(i,:) = [X(r,c),Y(r,c)];
        Vp(i,:) = V(r,c,:);
        plot(t,Vp(i,:),'DisplayName',['$x = (' num2str(p(i,1)) ', ' num2str(p(i,2)) ')$'],'LineWidth',1.5);
    end
    title('$V(t)$','Interpreter','latex','FontSize',16);
    legend('FontSize',16,'Interpreter','latex');
    legend show;    

    if(size(ps,1)==2)% two points, we measure the conductive velocity
        d = norm(p(1,:)-p(2,:)); % in cm
        d = 10*d; % in mm
        % measure first moment where V>=-65 mV
        th = -65;
        t1_j = find(Vp(1,:)>=th,1);
        t2_j = find(Vp(2,:)>=th,1);
        dt = t(t2_j)-t(t1_j); % in ms
        CV = d/dt;
        fprintf(solname);
        fprintf(' CV = %f mm/ms (m/s)\n',CV);
    end
end

