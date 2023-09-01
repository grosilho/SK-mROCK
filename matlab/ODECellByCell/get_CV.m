function [CV,transmembrane_G,y,t]=get_CV(folder,solname,pa,pb,V_th,old_file,transmembrane_G,y,t)
    
    if(nargin<=6)
        run([folder solname '_geo.m']);
        transmembrane_G = [];
        dom=mdom{1};
        for i=2:size(dom,1) 
            transmembrane_G = [transmembrane_G;dom{i}.G];
        end
        n_pts = size(transmembrane_G,1);
        if ~old_file
            % new reading, since we do not write ionic model vars any more    
            fileID = fopen([folder solname '_V_evolution.bin']);
            y = fread(fileID,'double');
            y = reshape(y,[1+n_pts,numel(y)/(1+n_pts)]);
            t = y(1,:);
            y = y(2:end,:);
        else
            % use this for old files where ionic vars are still part of the sol
            fileID = fopen([folder solname '_evolution.bin']);
            n_states=21; %change if not Courtemanche
            precision = [num2str(1+n_pts) '*double'];
            skip = n_pts*(n_states-1)*8;
            y = fread(fileID,precision,skip);
            y = reshape(y,[1+n_pts,numel(y)/(1+n_pts)]);
            t = y(1,:);
            y = y(2:end,:);
        end
    end

    ps = [pa;pb];
    np = size(ps,1);
    p_ind=zeros(np,1);
    t_ind = zeros(np,1);
    t_th = zeros(np,1);
    p = zeros(np,2);
    no_propag=0;
    for i=1:np
        d = transmembrane_G-ps(i,:);
        d = sqrt(d(:,1).^2+d(:,2).^2);
        [~,j] = min(d);
        p_ind(i) = j;
        p(i,:) = transmembrane_G(j,:);

        % measure first moment where V==V_th.
        % first find first moment where V>=V_th and then interpolate
        t_tmp = find(y(p_ind(i),:)>=V_th,1);
        if ~isempty(t_tmp)
            t_ind(i) = t_tmp;
            t_th(i) = t(t_ind(i)-1)+(t(t_ind(i))-t(t_ind(i)-1))*((V_th-y(p_ind(i),t_ind(i)-1))/(y(p_ind(i),t_ind(i))-y(p_ind(i),t_ind(i)-1)));
        else
            no_propag=1;
        end
    end
    if no_propag
        CV=0;
    else
        np = np/2;
        pa = p(1:np,:);
        pb = p((np+1):end,:);
        ta = t_th(1:np);
        tb = t_th((np+1):end);
        CV = zeros(np,1);
        for i=1:np
            dp = norm(pb(i,:)-pa(i,:)); % in cm
            dp = 10*dp; % in mm    
            dt = tb(i)-ta(i); % in ms
            CV(i) = dp/dt;
        end
        CV = mean(CV);
    end
end