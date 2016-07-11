classdef ion_simulation<handle
    properties
        Natoms
        x0
        m
        ions
        ejected
        t
        r_esc
        d
        fig
        
        quadrupole
        doppler
        lattice
    end
    
    properties (Constant,Hidden)
        %% Units
        cm = 1e-2;
        mm = 1e-3;
        um = 1e-6;
        nm = 1e-9;
        kW = 1e3;
        us = 1e-6;
        ns = 1e-9;
        uK = 1e-6;
        MHz = 1e6;
        %% funcamental constants
        amu = 1.660539040e-27;
        q = 1.60217662e-19;
        ke = 8.987551787 * 1e9;
        eps0 = 8.854187817e-12;
        h = 6.62607004e-34;
        hbar = 6.62607004e-34 / (2*pi);
        kB = 1.38064852e-23;
    end
    
    methods
        function obj = ion_simulation(N)
            obj.Natoms = N;
            obj.m = 25*obj.amu;
            obj.x0 = zeros(N,6);
            obj.ions = zeros(1,N,8);
            obj.t = 0:0.1*obj.us:20*obj.us;
            obj.r_esc = 100*obj.um;
            obj.d = 10*obj.um;
            
            obj.quadrupole.V = 20;
            obj.quadrupole.r = 5*obj.mm;
            
            obj.doppler.kv = [1 1 1];
            obj.doppler.lambda = 280*obj.nm;
            obj.doppler.Gamma = 2*pi*40*obj.MHz;
            
            laser.P = 10*obj.kW;
            laser.w0 = 60*obj.um;
            laser.lambda = 532*obj.nm;
            laser.alpha = -2.2*obj.cm^2;
            
            obj.lattice{1} = laser;
            obj.lattice{1}.w0 = 60*obj.um;
            obj.lattice{1}.kv = [0 0 1];
            
%             obj.lattice{2} = laser;
%             obj.lattice{2}.P = 10*obj.kW;
%             obj.lattice{2}.kv = [0 1 0];
% % 
%             obj.lattice{3} = laser;
%             obj.lattice{3}.kv = [0 0 1];
%             obj.lattice{3}.w0 = 300*obj.um;
        end
        
        function obj = init_pancake(obj)
            % generate 2D lattice starting grid
            obj.x0 = zeros(obj.Natoms, 6);
            obj.x0(:,1:2) = normrnd(zeros(obj.Natoms,2),(20*obj.um)*ones(obj.Natoms,2));
            obj.x0(:,1:3) = obj.x0(:,1:3) + normrnd(zeros(obj.Natoms,3),(10*obj.nm)*ones(obj.Natoms,3));
            
            Tionization = 100*obj.uK;
            velUncert = sqrt(obj.kB * Tionization / obj.m);
            obj.x0(:,4:6) = normrnd(zeros(obj.Natoms,3),velUncert*ones(obj.Natoms,3));
            obj.ions = zeros(numel(obj.t),obj.Natoms,6);
            obj.ions(1,:,:) = obj.x0;
            obj.ejected = [];
        end        
        
        function obj = init_2pancakes(obj)
            % generate 2D lattice starting grid
            obj.x0 = zeros(obj.Natoms, 6);
            obj.x0(:,1:2) = normrnd(zeros(obj.Natoms,2),(20*obj.um)*ones(obj.Natoms,2));
            spacing = round(obj.d/obj.lattice{1}.lambda);
            obj.x0(1:2:end,3) = -spacing*obj.lattice{1}.lambda/2;
            obj.x0(2:2:end,3) = spacing*obj.lattice{1}.lambda/2;
            obj.x0(:,1:3) = obj.x0(:,1:3) + normrnd(zeros(obj.Natoms,3),(10*obj.nm)*ones(obj.Natoms,3));
            
            Tionization = 100*obj.uK;
            velUncert = sqrt(obj.kB * Tionization / obj.m);
            obj.x0(:,4:6) = normrnd(zeros(obj.Natoms,3),velUncert*ones(obj.Natoms,3));
            obj.ions = zeros(numel(obj.t),obj.Natoms,6);
            obj.ions(1,:,:) = obj.x0;
            obj.ejected = [];
        end
        
        function obj = init_1d(obj)
            % generate 1d array
            obj.x0 = zeros(obj.Natoms, 6);
            spacing = ceil(2*obj.d/obj.lattice{1}.lambda);
            real_d = spacing*obj.lattice{1}.lambda/2;
            xl = obj.Natoms;
            xgrid = linspace(-(xl-1)/2,(xl-1)/2,xl)*real_d;
            obj.x0(:,3) = xgrid;
            pos_uncert = 10*obj.nm;
            obj.x0(:,1:3) = obj.x0(:,1:3) + normrnd(zeros(obj.Natoms,3),pos_uncert*ones(obj.Natoms,3));
            
            Tionization = 100*obj.uK;
            vel_uncert = sqrt(obj.kB * Tionization / obj.m);
            obj.x0(:,4:6) = normrnd(zeros(obj.Natoms,3),vel_uncert*ones(obj.Natoms,3));
            obj.ions = zeros(numel(obj.t),obj.Natoms,6);
            obj.ions(1,:,:) = obj.x0;
            obj.ejected = [];
        end        
        
        function obj = init_2d(obj)
            % generate 2d array
            obj.x0 = zeros(obj.Natoms, 6);
            spacing = ceil(2*obj.d/obj.lattice{1}.lambda);
            real_d = spacing*obj.lattice{1}.lambda/2;
            xl = round(obj.Natoms/sqrt(obj.Natoms));
            yl = obj.Natoms/xl;
            xgrid = linspace(-(xl-1)/2,(xl-1)/2,xl)*real_d;
            ygrid = linspace(-(yl-1)/2,(yl-1)/2,yl)*real_d;
            [X,Y] = meshgrid(xgrid,ygrid);
            obj.x0(:,1) = X(:);
            obj.x0(:,2) = Y(:);
            pos_uncert = 10*obj.nm;
            obj.x0(:,1:3) = obj.x0(:,1:3) + normrnd(zeros(obj.Natoms,3),pos_uncert*ones(obj.Natoms,3));
            
            Tionization = 100*obj.uK;
            vel_uncert = sqrt(obj.kB * Tionization / obj.m);
            obj.x0(:,4:6) = normrnd(zeros(obj.Natoms,3),vel_uncert*ones(obj.Natoms,3));
            obj.ions = zeros(numel(obj.t),obj.Natoms,6);
            obj.ions(1,:,:) = obj.x0;
            obj.ejected = [];
        end
        
        function obj = init_3d(obj)
            % generate 2d array
            obj.x0 = zeros(obj.Natoms, 6);
            spacing = ceil(2*obj.d/obj.lattice{1}.lambda);
            real_d = spacing*obj.lattice{1}.lambda/2;
            xl = round(obj.Natoms^(1/3));
            yl = xl;
            zl = yl;
            xgrid = linspace(-(xl-1)/2,(xl-1)/2,xl)*real_d;
            ygrid = linspace(-(yl-1)/2,(yl-1)/2,yl)*real_d;
            zgrid = linspace(-(zl-1)/2,(zl-1)/2,zl)*real_d;            
            [X,Y,Z] = meshgrid(xgrid,ygrid,zgrid);
            obj.x0(:,1) = X(:);
            obj.x0(:,2) = Y(:);
            obj.x0(:,3) = Z(:);
            pos_uncert = 10*obj.nm;
            obj.x0(:,1:3) = obj.x0(:,1:3) + normrnd(zeros(obj.Natoms,3),pos_uncert*ones(obj.Natoms,3));
            
            Tionization = 100*obj.uK;
            vel_uncert = sqrt(obj.kB * Tionization / obj.m);
            obj.x0(:,4:6) = normrnd(zeros(obj.Natoms,3),vel_uncert*ones(obj.Natoms,3));
            obj.ions = zeros(numel(obj.t),obj.Natoms,6);
            obj.ions(1,:,:) = obj.x0;
            obj.ejected = [];
        end
        
        function obj = reinit(obj)
            obj.x0 = squeeze(obj.ions(end,:,:));
            obj.Natoms = size(obj.x0,1);
        end
        
        function run(obj)
            time = clock;
            fprintf('Starting simulation with %d ions and %d time steps.\nTime: %i.%i.%i, %i:%i\n',obj.Natoms,numel(obj.t),time(3),time(2),time(1),time(4),time(5));
            tic();
            odefun = @(t,x) obj.Force(x);
            obj.ions = zeros(numel(obj.t),size(obj.x0,1),6);
            obj.ejected = [];
            obj.ions(1,:,:) = obj.x0;
            for i = 2:numel(obj.t)
                if any(obj.ions)
                    vec0 = reshape(obj.ions(i-1,:,:),[size(obj.ions,2)*size(obj.ions,3),1]);
                    [~, vec1] = ode45(odefun, [obj.t(i-1) mean([obj.t(i-1),obj.t(i)]) obj.t(i)], vec0);
                    vec0 = vec1(3,:);
                    obj.ions(i, :, :) = reshape(vec0, [size(obj.ions,2),size(obj.ions,3)]);
                    % remove all escaping atoms
                    r_abs = squeeze(sqrt(sum(obj.ions(i,:,1:3).^2,3)));
                    esc = r_abs > obj.r_esc;
                    if any(esc)
                        for j=find(esc)
                            obj.ejected(:,size(obj.ejected,2)+1,:) = obj.ions(:,j,:);
                            obj.ejected(i+1:end,end,:) = repmat(obj.ejected(i,end,:),size(obj.ejected,1)-i,1,1);
                        end
                        obj.ions = obj.ions(:,~esc,:);
                        fprintf('Ion Escaped! %d ions are left.\n',size(obj.ions,2));
                    end
                end
                if i==2
                    fprintf('Simulation expected to take %.2f seconds.\n',toc()*length(3:numel(obj.t)))
                end
            end
            fprintf('Simulation finished. Time = %.2f seconds.\n',toc())
        end
        
        function F = Force(obj,X)
            N = length(X)/6;
            X = reshape(X,[N,6]);
            F = zeros(N,6);
            F(:,1:3) = X(:,4:6);
            F(:,4:6) =  (...
                obj.F_opt(X(:,1:3)) + ...
                obj.F_quad(X(:,1:3)) + ...
                obj.F_damp(X(:,4:6)) + ...
                obj.F_Coulomb(X(:,1:3)) ...
                )/obj.m;
            F = reshape(F,[N*6,1]);
        end
        
        function normal_modes(obj)
            X = squeeze(obj.ions(end,:,1:3));
            N = size(X,1);
            function f = force(x)
                f = zeros(size(x));
                f = f + obj.F_opt(x);
                f = f + obj.F_quad(x);
                f = f + obj.F_Coulomb(x);
            end
            
            function f = lin_force(x,i)
                dx = zeros(size(x));
                dx(i) = 10*eps;
                f = (force(x+dx)-force(x))/(10*eps);
            end

            % Calculate the Hessian
            H = -cell2mat(arrayfun(@(i) reshape(lin_force(X,i),N*3,1),1:(N*3),'uni',false))/obj.m;
            [v,l] = eig(H);
            % get the eigenfrequencies in kHz
            freq = sqrt(abs(diag(l)))/(2*pi*1e3);
            % sort the eigenvectors by eigenfrequency
            [freq,i] = sort(freq);
            v = v(:,i);
            
            % Plotting
            fig1 = figure(1); clf;
            fig1.KeyPressFcn = {@plot_mode};
            ax1 = subplot(6,1,6);
            ax1.YTick = [];
            zoom('xon');
            zoom('off');
            ylim([0,2]);
            xlabel('Normal Mode Frequency (kHz)');
            
            ax2 = subplot(6,1,1:4); box on;
            c = ax2.ColorOrder;
            menu = uicontextmenu;
            ax2.UIContextMenu = menu;
            uimenu(menu,'label','isometric','callback',@(~,~) view([1 1 1]));
            uimenu(menu,'label','x-y','callback',@(~,~) view([0 0 1]));
            drum = uimenu(menu,'label','drum','callback',@toggle);
            function toggle(source,~)
                if strcmp(source.Checked,'on')
                    source.Checked = 'off';
                else
                    source.Checked = 'on';
                end
                plot_mode();
            end
            
            last_mode = 1;
            plot_mode()
            rmax = 1.2*sqrt(max(max(sum(obj.ions(:,:,1:3).^2,3))))/obj.um;
            xlim([-rmax,rmax]);
            ylim([-rmax,rmax]);
            zlim([-rmax,rmax]);
            view([1 1 1]);
            grid('on');
            xlabel('$x (\mu m)$','fontsize',16,'interpreter','latex');
            ylabel('$y (\mu m)$','fontsize',16,'interpreter','latex');
            zlabel('$z (\mu m)$','fontsize',16,'interpreter','latex');
            
            function plot_mode(~,event)
                if nargin>0
                    if isa(event,'matlab.graphics.eventdata.Hit')
                        [~,mode_i] = min(abs(event.IntersectionPoint(1)-freq));
                    elseif isa(event,'matlab.ui.eventdata.KeyData')
                        if strcmp(event.Key,'leftarrow')
                            mode_i = last_mode-1;
                            mode_i = max(min(mode_i,length(freq)),1);
                        elseif strcmp(event.Key,'rightarrow')
                            mode_i = last_mode+1;
                            mode_i = max(min(mode_i,length(freq)),1);
                        end
                    end
                else
                    mode_i = last_mode;
                end
                last_mode = mode_i;
                
                cla(ax1); axes(ax1); hold on; box on;
                p1 = stem(freq,ones(size(freq)),'.','linewidth',2,'color',[1 1 1]*0.6);
                p1.ButtonDownFcn = {@plot_mode};
                p2 = stem(freq(mode_i),1,'.','linewidth',3,'color',[0 0.5 1]);
                p2.ButtonDownFcn = {@plot_mode};
                cla(ax2); axes(ax2); hold on;             
                
                u = obj.d*reshape(v(:,mode_i),N,3);
                if strcmp(drum.Checked,'on')
                    [xg,yg] = meshgrid(linspace(min(X(:,1)),max(X(:,1)),100),linspace(min(X(:,2)),max(X(:,2)),100));
                    zg = griddata(X(:,1),X(:,2),u(:,3),xg,yg);
                    pcolor(xg/obj.um,yg/obj.um,zg/obj.um);
                    shading('flat');
                    plot(X(:,1)/obj.um,X(:,2)/obj.um,'k.','markersize',30)
                else
                    if std(X(:,3))>obj.d/10
                        z = X(:,3)>0;
                    else
                        z = logical(X(:,3));
                    end
                    quiver3(X(z,1)/obj.um,X(z,2)/obj.um,X(z,3)/obj.um,u(z,1)/obj.um,u(z,2)/obj.um,u(z,3)/obj.um,'color',c(1,:),'linewidth',2,'autoscale','on','autoscalefactor',0.5,'marker','.','markersize',30,'alignvertexcenters','on');
                    quiver3(X(~z,1)/obj.um,X(~z,2)/obj.um,X(~z,3)/obj.um,u(~z,1)/obj.um,u(~z,2)/obj.um,u(~z,3)/obj.um,'color',c(2,:),'linewidth',2,'autoscale','on','autoscalefactor',0.5,'marker','.','markersize',30,'alignvertexcenters','on');
                end
            end
        end            
        
        function F = F_quad(obj,X)
            F = -obj.q*(obj.quadrupole.V/obj.quadrupole.r^2)*X;
            F(:,3) = -2*F(:,3);
        end  
        
        function F = F_damp(obj,V)
            damping = (1/4)*obj.hbar*(2*pi/obj.doppler.lambda*obj.doppler.kv/norm(obj.doppler.kv)).^2;
            F = -bsxfun(@times, damping, V);
        end
        
        function F = F_Coulomb(obj,X)
            N = size(X,1);
            Rv = zeros(3, N, N);
            for i = 1:3
                Rv(i,:,:) = bsxfun(@minus,X(:,i),X(:,i)');
            end
            % calculate distance 3/2
            R3 = (squeeze(1 ./ sum(Rv.^2, 1))).^(3/2);
            R3(1:(N+1):end) = 0;
            % calculate the Coulomb force
            F = zeros(N,3);
            for i = 1:3
                F(:,i) = obj.ke*obj.q^2*sum(R3.*squeeze(Rv(i,:,:)), 2);
            end
        end
        
        function F = F_opt(obj,X)
            F = zeros(size(X));
            for i=1:numel(obj.lattice)
                kv = obj.lattice{i}.kv;
                P = obj.lattice{i}.P;
                w0 = obj.lattice{i}.w0;
                lambda = obj.lattice{i}.lambda;
                zR = pi*w0^2/lambda;
                k = 2*pi/lambda;
                U0 = obj.lattice{i}.alpha*8*P./(pi*w0.^2)*obj.h;
                
                r = sqrt(kv(1)^2 + kv(2)^2);
                p = sqrt(kv(1)^2 + kv(2)^2 + kv(3)^2);
                if r>0
                    R = [   kv(1)*kv(2)/(r*p),      kv(2)*kv(3)/(r*p),      -r/p;
                            -kv(2)/r,               kv(1)/r,                0;
                            kv(1)/p,                kv(2)/p,                kv(3)/p];
                else
                    R = eye(3);
                end
                Xr = (R*X')';
                
                x = Xr(:,1);
                y = Xr(:,2);
                z = Xr(:,3);
                
                Theta = 1 + z.^2 / zR^2;
                rsq = x.^2 + y.^2;
                cosZ = cos(k*z);
                Pre = 2*U0*cosZ ./ Theta .* exp(-2*rsq ./ (w0^2 * Theta));
                
                vec = [ ...
                    2 * x .* cosZ ./ (w0^2 * Theta), ...
                    2 * y .* cosZ ./ (w0^2 * Theta), ...
                    -2 * z .* rsq .* cosZ ./ (w0^2 * Theta.^2 * zR^2) + ...
                    z .* cosZ ./ (Theta * zR^2) + ...
                    k * sin(k*z)];
                
                F = F + ((R^-1)*bsxfun(@times, Pre, vec)')';
            end
        end     
                
        function plot_potential(obj)
            xl = linspace(-obj.r_esc/50,obj.r_esc/50,200);
            yl = linspace(-obj.r_esc/50,obj.r_esc/50,200);
            [X,Y] = meshgrid(xl,yl);
            U = zeros(size(X));
            z = zeros(1000,1);
            for i=1:numel(U)
                x = linspace(0,X(i),size(z,1))';
                y = linspace(0,Y(i),size(z,1))';
                U(i) = sum(bsxfun(@dot,obj.F_opt([x(2:end) y(2:end) z(2:end)])',[diff(x) diff(y) diff(z)]'));
            end
            
            obj.fig = figure(1); clf; hold on;
            obj.fig.WindowStyle = 'docked';
            pcolor(X/obj.um,Y/obj.um,(U-max(U(:)))/obj.h);
            colormap('bone');
            shading('interp');
            xlabel('$x (\mu m)$','fontsize',16,'interpreter','latex');
            ylabel('$y (\mu m)$','fontsize',16,'interpreter','latex');
            title('Potential (Hz)','fontsize',16);
            colorbar();
        end
        
        function plot_ion(obj,ion_i)
            if ~exist('ion_i','var')
                ion_i = 1:size(obj.ions,2);
            end
            obj.fig = figure(1); clf; hold on;
            obj.fig.WindowStyle = 'docked';
            s(1) = subplot(3,1,1); hold on; box on;
            plot(obj.t/obj.us,obj.ions(:,ion_i,1)/obj.um);
            ylabel('$x (\mu m)$','fontsize',16,'interpreter','latex');
            s(1).XTickLabel = '';
            s(2) = subplot(3,1,2); hold on; box on;
            plot(obj.t/obj.us,obj.ions(:,ion_i,2)/obj.um);
            ylabel('$y (\mu m)$','fontsize',16,'interpreter','latex');
            s(2).XTickLabel = '';
            s(3) = subplot(3,1,3); hold on; box on;
            plot(obj.t/obj.us,obj.ions(:,ion_i,3)/obj.um);
            xlabel('$t (\mu s)$','fontsize',16,'interpreter','latex');
            ylabel('$z (\mu m)$','fontsize',16,'interpreter','latex');
            linkaxes(s,'x');
            xlim([min(obj.t)/obj.us,max(obj.t)/obj.us]);
        end
        
        function plot_ejected(obj,ion_i)
            if ~exist('ion_i','var')
                ion_i = 1:size(obj.ejected,2);
            end
            obj.fig = figure(1); clf; hold on;
            obj.fig.WindowStyle = 'docked';
            s(1) = subplot(3,1,1); hold on; box on;
            plot(obj.t/obj.us,obj.ejected(:,ion_i,1)/obj.um);
            ylabel('$x (\mu m)$','fontsize',16,'interpreter','latex');
            s(1).XTickLabel = '';
            s(2) = subplot(3,1,2); hold on; box on;
            plot(obj.t/obj.us,obj.ejected(:,ion_i,2)/obj.um);
            ylabel('$y (\mu m)$','fontsize',16,'interpreter','latex');
            s(2).XTickLabel = '';
            s(3) = subplot(3,1,3); hold on; box on;
            plot(obj.t/obj.us,obj.ejected(:,ion_i,3)/obj.um);
            xlabel('$t (\mu s)$','fontsize',16,'interpreter','latex');
            ylabel('$z (\mu m)$','fontsize',16,'interpreter','latex');
            linkaxes(s,'x');
            xlim([min(obj.t)/obj.us,max(obj.t)/obj.us]);
        end
        
        function plot_frequency(obj,ion_i)
            if ~exist('ion_i','var')
                ion_i = 1:size(obj.ions,2);
            end
            obj.fig = figure(5); clf; hold on; box on;
            obj.fig.WindowStyle = 'docked';

            ffty = fftshift(fft(obj.ions,[],1));
            freq = 1e-6*0.5*linspace(-1,1,length(obj.t))/obj.dt();
            ffty = ffty(freq>0,:,:);
            freq = freq(freq>0);
            
            s(1) = subplot(3,1,1); hold on; box on;
            plot(freq,abs(ffty(:,ion_i,1)));
            ylabel('$x (\mu m)$','fontsize',16,'interpreter','latex');
            s(1).XTickLabel = '';
            s(2) = subplot(3,1,2); hold on; box on;
            plot(freq,abs(ffty(:,ion_i,2)));
            ylabel('$y (\mu m)$','fontsize',16,'interpreter','latex');
            s(2).XTickLabel = '';
            s(3) = subplot(3,1,3); hold on; box on;
            plot(freq,abs(ffty(:,ion_i,3)));
            xlabel('$t (\mu s)$','fontsize',16,'interpreter','latex');
            ylabel('$z (\mu m)$','fontsize',16,'interpreter','latex');
            linkaxes(s,'x');
            xlim([0,max(freq)]);
        end        
        
        function scatter(obj,ion_i)
            if ~exist('ion_i','var')
                ion_i = 1:size(obj.ions,2);
            end            
            obj.fig = figure(1); clf; hold on; box on;
            obj.fig.WindowStyle = 'docked';
            a = gca;
            c = a.ColorOrder;
            plot(obj.ions(1,ion_i,1)/obj.um, obj.ions(1,ion_i,2)/obj.um,'o','color',c(1,:),'MarkerSize',10);
            plot(obj.ions(:,ion_i,1)/obj.um, obj.ions(:,ion_i,2)/obj.um, '-');
            plot(obj.ions(end,ion_i,1)/obj.um,obj.ions(end,ion_i,2)/obj.um,'.','MarkerSize',30,'color',c(1,:));
            xlabel('$x (\mu m)$','fontsize',16,'interpreter','latex');
            ylabel('$y (\mu m)$','fontsize',16,'interpreter','latex');            
        end
        
        function scatter_initial(obj,varargin)
            p = inputParser;
            addOptional(p,'dim',2,@isnumeric);
            parse(p,varargin{:});
            dim = p.Results.dim;            
            
            obj.fig = figure(1); clf; hold on; box on; grid on;
            obj.fig.WindowStyle = 'docked';
            if ~isempty(obj.ions)
                if dim==3
                    plot3(obj.ions(1,:,1)/obj.um,obj.ions(1,:,2)/obj.um,obj.ions(1,:,3)/obj.um,'.','MarkerSize',30);
                else
                    plot(obj.ions(1,:,1)/obj.um,obj.ions(1,:,2)/obj.um,'.','MarkerSize',40);
                end
            end
            if ~isempty(obj.ejected)
                if dim==3
                    plot3(obj.ejected(1,:,1)/obj.um,obj.ejected(1,:,2)/obj.um,obj.ejected(1,:,3)/obj.um,'o','MarkerSize',8);
                else
                    plot3(obj.ejected(1,:,1)/obj.um,obj.ejected(1,:,2)/obj.um,'o','MarkerSize',10);
                end
            end
            xlabel('$x (\mu m)$','fontsize',16,'interpreter','latex');
            ylabel('$y (\mu m)$','fontsize',16,'interpreter','latex');
            if dim==3
                view([2 1 1]);
                zlabel('$z (\mu m)$','fontsize',16,'interpreter','latex');
                
                xl = round(obj.Natoms^(1/3));
                yl = xl;
                zl = yl;
                spacing = ceil(2*obj.d/obj.lattice{1}.lambda);
                real_d = spacing*obj.lattice{1}.lambda/2;
                xgrid = linspace(-(xl-1)/2,(xl-1)/2,xl)*real_d/obj.um;
                ygrid = linspace(-(yl-1)/2,(yl-1)/2,yl)*real_d/obj.um;
                zgrid = linspace(-(zl-1)/2,(zl-1)/2,zl)*real_d/obj.um;
                for z=zgrid
                    for x=xgrid
                        plot3(x*[1 1],[min(ygrid) max(ygrid)],z*[1 1],'k','color',0.75*[1 1 1])
                    end
                    for y=xgrid
                        plot3([min(xgrid) max(xgrid)],y*[1 1],z*[1 1],'k','color',0.75*[1 1 1])
                    end                    
                end
                for x=xgrid
                    for y=ygrid
                        plot3(x*[1 1],y*[1 1],[min(zgrid) max(zgrid)],'k','color',0.75*[1 1 1])
                    end
                end
            else
                legend({'confined ion','escaping ion'},'fontsize',16);
            end
        end
        
        function scatter_final(obj,varargin)
            p = inputParser;
            addOptional(p,'partitioned',false,@islogical);
            parse(p,varargin{:});
            partitioned = p.Results.partitioned;
            
            obj.fig = figure(1);
            obj.fig.WindowStyle = 'docked';
            clf; hold on; box on;
            rmax = 1.2*sqrt(max(max(sum(obj.ions(:,:,1:3).^2,3))))/obj.um;
            if ~isempty(obj.ions)
                if partitioned
                    plane = squeeze(obj.ions(end,:,3))>0;
                    plot(obj.ions(end,plane,1)/obj.um, obj.ions(end,plane,2)/obj.um,'.','MarkerSize',40);
                    plot(obj.ions(end,~plane,1)/obj.um, obj.ions(end,~plane,2)/obj.um,'.','MarkerSize',40);
                else
                    plot(obj.ions(end,:,1)/obj.um, obj.ions(end,:,2)/obj.um,'.','MarkerSize',40);
                end
            end
            xlim([-rmax,rmax]);
            ylim([-rmax,rmax]);
            grid(gca,'on');
            a = gca;
            a.TickLabelInterpreter = 'latex';
            a.FontSize = 12;
            xlabel('$x (\mu m)$','fontsize',16,'interpreter','latex');
            ylabel('$y (\mu m)$','fontsize',16,'interpreter','latex');
            if partitioned
                legend({'$z > 0$','$z < 0$'},'fontsize',16,'interpreter','latex');
            end
        end
        
        function plot_chain(obj)
            obj.fig = figure(1);
            obj.fig.WindowStyle = 'docked';
            clf; hold on; box on;
            xlabel('$x (\mu m)$','fontsize',16,'interpreter','latex');
            ylabel('$y (\mu m)$','fontsize',16,'interpreter','latex');
            zlabel('$z (\mu m)$','fontsize',16,'interpreter','latex');
            rmax = 1.2*sqrt(max(max(sum(obj.ions(:,:,1:3).^2,3))))/obj.um;
            if ~isempty(obj.ions)
                plot3(obj.ions(end,:,1)/obj.um, obj.ions(end,:,2)/obj.um, obj.ions(end,:,3)/obj.um,'k-','MarkerSize',30);                
                plot3(obj.ions(end,:,1)/obj.um, obj.ions(end,:,2)/obj.um, obj.ions(end,:,3)/obj.um,'.','MarkerSize',30);
            end
            xlim([-rmax,rmax]);
            ylim([-rmax,rmax]);
            zlim([-rmax,rmax]);
            view([4 3 1]);
            grid(gca,'on');
        end
        
        function varargout = movie(obj)
            obj.fig = figure(1);
            obj.fig.WindowStyle = 'docked';
            clf; hold on; box on;
            xlabel('$x (\mu m)$','fontsize',16,'interpreter','latex');
            ylabel('$y (\mu m)$','fontsize',16,'interpreter','latex');
            zlabel('$z (\mu m)$','fontsize',16,'interpreter','latex');
            rmax = 1.2*sqrt(max(max(sum(obj.ions(:,:,1:3).^2,3))))/obj.um;
            if nargout>0
                F(1) = getframe(gca);
                F = repmat(F,1,numel(obj.t));
            end
            for i=2:numel(obj.t)
                cla();
                past = max(1,i-10);
                if ~isempty(obj.ions)
                    plot3(obj.ions(i,:,1)/obj.um, obj.ions(i,:,2)/obj.um, obj.ions(i,:,3)/obj.um,'.-','MarkerSize',30,'linestyle','none');
                    plot3(obj.ions(past:i,:,1)/obj.um, obj.ions(past:i,:,2)/obj.um, obj.ions(past:i,:,3)/obj.um,'-k');
                end
                if ~isempty(obj.ejected)
                    plot3(obj.ejected(i,:,1)/obj.um, obj.ejected(i,:,2)/obj.um, obj.ejected(i,:,3)/obj.um,'.-','MarkerSize',30,'linestyle','none');
                    plot3(obj.ejected(past:i,:,1)/obj.um, obj.ejected(past:i,:,2)/obj.um, obj.ejected(past:i,:,3)/obj.um,'-k');
                end
                xlim([-rmax,rmax]);
                ylim([-rmax,rmax]);
                zlim([-rmax,rmax]);
                view([1 1 1]);
                grid(gca,'on');
                drawnow();
                if nargout>0
                    F(i) = getframe(gca);
                end
            end
            if nargout==1
                varargout{1} = F;
            end
        end
        
        function movie2d(obj,varargin)
            p = inputParser;
            addOptional(p,'partitioned',false,@islogical);
            parse(p,varargin{:});
            partitioned = p.Results.partitioned;
            
            obj.fig = figure(1);
            obj.fig.WindowStyle = 'docked';
            clf; hold on; box on;
            xlabel('$x (\mu m)$','fontsize',16,'interpreter','latex');
            ylabel('$y (\mu m)$','fontsize',16,'interpreter','latex');
            rmax = 1.2*sqrt(max(max(sum(obj.ions(:,:,1:3).^2,3))))/obj.um;
            for i=1:numel(obj.t)
                cla();
                past = max(1,i-10);
                if ~isempty(obj.ions)
                    plane = squeeze(obj.ions(i,:,3))>0;
                    if partitioned
                        plot(obj.ions(i,plane,1)/obj.um, obj.ions(i,plane,2)/obj.um,'.','MarkerSize',30);
                        plot(obj.ions(i,~plane,1)/obj.um, obj.ions(i,~plane,2)/obj.um,'.','MarkerSize',30);
                    else
                        plot(obj.ions(i,:,1)/obj.um, obj.ions(i,:,2)/obj.um,'.','MarkerSize',30);
                    end
                    plot(obj.ions(past:i,:,1)/obj.um, obj.ions(past:i,:,2)/obj.um,'-k');
                end
                xlim([-rmax,rmax]);
                ylim([-rmax,rmax]);
                grid(gca,'on');
                drawnow();
            end
            if partitioned
                legend({'z > 0','z < 0'},'fontsize',16);
            end
        end
        
        function varargout = tf(obj,tf_)
            if ~exist('tf_','var')
                varargout{1} = max(obj.t);
            else
                t0_ = min(obj.t);
                dt_ = min(diff(obj.t));
                obj.t = t0_:dt_:tf_;
            end
        end
        
        function varargout = dt(obj,dt_)
            if ~exist('dt_','var')
                varargout{1} = min(diff(obj.t));
            else            
                t0_ = min(obj.t);
                tf_ = max(obj.t);
                obj.t = t0_:dt_:tf_;
            end
        end
        
        function min_d = min_d(obj)
            Rv = zeros(3, size(obj.ions,2), size(obj.ions,2));
            for i = 1:3
                Rv(i,:,:) = bsxfun(@minus,obj.ions(end,:,i),obj.ions(end,:,i)');
            end
            % calculate distance 3/2
            R = squeeze(sqrt(sum(Rv.^2, 1)));
            R(1:(size(R,1)+1):end) = inf;
            min_d = min(R);
        end
        
        function print_figure(obj,filename)
            if ~exist('filename','var')
                filename = 'figure.pdf';
            end
            if exist('filename','file')
                delete('filename');
            end
            temp_fig = figure('visible','off');
            copyobj(obj.fig.Children,temp_fig);
            temp_fig.PaperPosition = [0 0 8 6];
            temp_fig.PaperSize = [8 6]; %
            temp_fig.PaperPositionMode = 'manual';
            print(temp_fig,filename,'-dpdf');
            fprintf('figure printed\n');
            delete(temp_fig);
        end
        
    end
end