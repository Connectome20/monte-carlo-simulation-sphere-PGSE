% ********** Setup the directory on your computer **********
clear
restoredefaultpath
filePath = matlab.desktop.editor.getActiveFilename;
root0 = fileparts(filePath);
addpath(genpath(fullfile(root0,'lib')));
root = fullfile(root0,'data');
root_code = fullfile(root0,'lib','packing');

projname = 'sphere_PGSE';

%% Generate random sphere packing

% Sphere diameter (um) and intra-spherical volume fraction
d_mean = [5, 10, 15];      % mean diameter
f      = [0.2, 0.4, 0.6];  % volume fraction

% Membrane permeability, um/ms
kappa = [0.01, 0.03, 0.05];

% Simulation parameters
dt = 2e-4;              % time of each step, ms
TN = ceil((120+6)/dt);  % # steps
NPar = 1e5;             % # random walkers, usually larger than 1e5
Din = 1;                % ICS diffusivity, um^2/ms
Dex = 2;                % ECS diffusivity, um^2/ms
pinit = 3;              % initial position, 1=ICS, 2=ECS, 3=ICS+ECS, 4=center
threadpb = 256;         % thread per block for cuda

% Generate randomly packed spheres and save simulation parameters
% You may need to complie the packing generation code. You can setup
% complieFlag = 1 when seed = 1.
seed = 0;
for i = 1:numel(d_mean)
    di = d_mean(i);
    for j = 1:numel(f)
        fi = f(j);
        for k = 1:numel(kappa)
            seed = seed+1;
            kappai = kappa(k);
            target = fullfile(root,projname,sprintf('sphere_%04u',seed));
            mkdir(target);

            % Generate packing
            ps = packsph();     % class for randomly packed spehres
            maxdensity = fi;    % maximal spherical volume fraction
            hardwallBC = 0;     % boundary condition, 0 for periodic, 1 for hard wall
            N = 999;            % # sphere                
            rinit = (di/2)*ones(N,1); % sphere radius (um)
            compileFlag = 0;    % compile Donev's C code, 0 for not compile, 1 for compile
            n = 200;            % matrix size of the lookup table, n x n x n
            gap = 0;            % minimal distance between spheres (µm), default=0
            ps.packing(maxdensity,hardwallBC,rinit,root_code,target,compileFlag,n,gap)

            % Save simulation parameters
            fileID = fopen(fullfile(target,'simParamInput.txt'),'w');
            fprintf(fileID,sprintf('%g\n%u\n%u\n%g\n%g\n%g\n%u\n%u\n',dt,TN,NPar,Din,Dex,kappai,pinit,threadpb));
            fclose(fileID);

            % Create b-table for wide pulse PGSE
            % For PGSE, we save [diffusion time, pulse width] and
            % [b-value, gradient direction] seperately. This helps to
            % accelerate the phase calculation for diffusion signals.
            TD = 6:120;                    % diffusion time, ms
            Td = 6*ones(size(TD));         % pulse width, ms
            
            bval = 0.2:0.2:1;              % b-value, ms/um^2
            bvec = [1 0 0; 0 1 0; 0 0 1];  % gradient direction
        
            NDelta = numel(TD);
            DELdel = zeros(NDelta,2);
            for ii = 1:numel(TD)
                DELdel(ii,:) = [TD(ii), Td(ii)];
            end
            DELdel = DELdel.';
            DELdel = DELdel(:);
            fid = fopen(fullfile(target,'gradient_NDelta.txt'),'w');
            fprintf(fid,sprintf('%u\n',NDelta));
            fclose(fid);

            fid = fopen(fullfile(target,'gradient_DELdel.txt'),'w');
            fprintf(fid,sprintf('%.8f\n',DELdel));
            fclose(fid);
            
            ig = 0;
            btab = zeros(numel(bval)*size(bvec,1),4);
            for jj = 1:numel(bval)
                bvalj = bval(jj);
                for kk = 1:size(bvec,1)
                    ig = ig+1;
                    bveck = bvec(kk,:);
                    btab(ig,:) = [bvalj bveck];
                end
            end
            btab = btab.';
            btab = btab(:);
            fid = fopen(fullfile(target,'gradient_Nbtab.txt'),'w');
            fprintf(fid,sprintf('%u\n',numel(bval)*size(bvec,1)));
            fclose(fid);
        
            fid = fopen(fullfile(target,'gradient_btab.txt'),'w');
            fprintf(fid,sprintf('%.8f\n',btab));
            fclose(fid);
        end
    end
end

%% Create a shell script to run the codes

fileID = fopen(fullfile(root,projname,'job.sh'),'w');
root_cuda = fullfile(root0,'lib','rms');
fprintf(fileID,'#!/bin/bash\n');
for i = 1:numel(d_mean)*numel(f)*numel(kappa)
    target = fullfile(root,projname,sprintf('sphere_%04u',i));
    fprintf(fileID,sprintf('cd %s\n',target));
    fprintf(fileID,sprintf('cp -a %s .\n',fullfile(root_cuda,'main_PGSE_bvec_cuda')));
    fprintf(fileID,'./main_PGSE_bvec_cuda\n');
end
fclose(fileID);

% Open the terminal window in the project folder and run "sh job.sh"

% You may need to open the terminal in the root_cuda folder and compile the 
% CUDA code using "nvcc main_PGSE_bvec.cu -o main_PGSE_bvec_cuda"

%% Load data for all time points
Dtot = [];  % time-dependent diffusivity for all geometries
Ktot = [];  % time-dependent kurtosis for all geometries

dtot = [];  % mean diffusivity (um) for all geometries
ftot = [];  % intra-spherical volume fraction for all geometries
katot = []; % membrane permeability (um/ms) for all geometries
Ditot = []; % intrinsic diffusivity inside spheres (um^2/ms) for all geometries
Detot = []; % intrinsic diffusivity outside spheres (um^2/ms) for all geometries
tptot = []; % time to the kurtosis peak (ms) for all geometries
for i = 1:27
    rms = simul3Dsphere_cuda_pgse_bvec(fullfile(root,projname,sprintf('sphere_%04u',i)));
    sigi = rms.sig;
    bvali = rms.bval;
    Di = zeros(rms.NDel,1);
    Ki = zeros(rms.NDel,1);
    for j = 1:rms.NDel
        sigj = sigi(:,j);
        [C, IA, IC] = unique(bvali);
        sigk = zeros(numel(C),1);
        for k = 1:numel(C)
            Ik = IC==k;
            sigk(k) = mean(sigj(Ik));
        end
        A = [-C, 1/6*C.^2, C.^(3:5)];
        X = A\log(sigk(:));
        Di(j) = X(1);
        Ki(j) = X(2)/X(1)^2;
    end
    Dtot = cat(2,Dtot,Di);
    Ktot = cat(2,Ktot,Ki);

    [~,I] = max(Ki);
    tptot = cat(1, tptot, rms.Del(I));

    ri = load(fullfile(root,projname,sprintf('sphere_%04u',i),'phantom_rCir.txt'));
    vi = load(fullfile(root,projname,sprintf('sphere_%04u',i),'phantom_res.txt'));
    dtot = cat(1, dtot, mean(ri)*vi*2);
    ftot = cat(1, ftot, sum(4/3*pi*ri.^3));
    katot = cat(1, katot, rms.kappa);
    Ditot = cat(1, Ditot, rms.Din);
    Detot = cat(1, Detot, rms.Dex);

end

% Theoretical prediction of time to the kurtosis peak for narrow pulse PG
tpeak_np = sqrt(2/5*(dtot/2).^3./katot./Detot) + rms.del(1)/3;

% Theoretical prediction of time to the kurtosis peak for wide pulse PG
tpeak_wp = sqrt(2/5*(dtot/2).^3./katot./Detot) .* sqrt(16/35*(dtot/2).^2./Ditot/rms.del(1)) + rms.del(1)/3;

% Diffusion time (ms)
Del = rms.Del;

% intra-spherical residence time (ms)
tr = 1./(katot*3./(dtot/2));

% intra-spherical correlation time (ms)
tD = (dtot/2).^2./Ditot;

%% Plot figure: time-dependent kurtosis
close all
d_mean = [5 10 15];               % mean diameter (um)
f      = [0.2 0.4 0.6];           % intra-spherical volume fraction
kappa =  [0.01 0.03 0.05];        % membrane permeability (um/ms)

nrow = 3;
ncol = 3;
figure('unit','inch','position',[0 0 19.37/2 8.98]);
cmap = colormap('lines');
seed = 0;
for i = 1:numel(d_mean)
    for j = 1:numel(f)
        clear h lgtxt
        for k = 1:numel(kappa)
            seed = seed+1;
            Ki = Ktot(:,seed);
            [~,I] = min(abs(tptot(seed)-Del));

            subplot(nrow,ncol,(i-1)*numel(f)+j)
            hold on;
            h(k) = plot(Del, Ki, '-', 'color', cmap(k,:), 'linewidth', 1);
            plot(Del(I), Ki(I), 'v', 'color', cmap(k,:), 'markersize', 8, 'linewidth', 1);
            lgtxt{k} = sprintf('$\\kappa=%.2f$', kappa(k));
        end
        title(sprintf('($d, f$) = (%u $\\mu$m, %.1f)', d_mean(i), f(j)), 'interpreter','latex','fontsize',14)
        if i==1 && j==1
            hl = legend(h, lgtxt, 'interpreter', 'latex', 'box', 'off', 'location', 'northeast', 'fontsize', 14);
            hl.Position = hl.Position + [0.015 0.01 0 0];
            hl.ItemTokenSize = [18 18];
        end

        xticks(0:20:120)

        xlim([0 120]);
        if j==1
            ylim([0 1]);
        elseif j==2
            ylim([0 2]);
        elseif j==3
            ylim([0 4]);
        end
        box on;
        grid on;
        if i==3
            xlabel('$\Delta$, ms', 'interpreter','latex','fontsize',14)
        end
        if j==1
            ylabel('$K(\Delta)$', 'interpreter','latex','fontsize',14)
        end
    end
end

%% Plot figure: time to the kurtosis peak

nrow = 1;
ncol = 3;
figure('unit','inch','position',[0 0 19.37/4*3 8.98/2]);
cmap = colormap('lines');
seed = 0;
mk = {'o','d','s'};
for i = 1:numel(d_mean)
    for j = 1:numel(f)
        clear h lgtxt
        for k = 1:numel(kappa)
            seed = seed+1;
            subplot(nrow,ncol,j)
            hold on;
            if tr(seed) < tD(seed)
                h(k) = plot(tptot(seed), tpeak_np(seed), '-', 'color', cmap(i,:),...
                'marker', mk{k}, 'markersize', 8, 'linewidth', 1, 'MarkerFaceColor', 0.85*[1 1 1]);
            else
                h(k) = plot(tptot(seed), tpeak_np(seed), '-', 'color', cmap(i,:),...
                    'marker', mk{k}, 'markersize', 8, 'linewidth', 1);
            end
        end
        title(sprintf('$f$ = %.1f', f(j)), 'interpreter','latex','fontsize',16)
        xlim([0 100]);
        ylim([0 100]);
        box on;
        grid on;
        xlabel('$\Delta_{\rm peak}$, ms', 'interpreter','latex','fontsize',16)
        ylabel('$t_{\rm peak}+\frac{1}{3}\delta$, ms', 'interpreter','latex','fontsize',16)
        hr = refline(1,0); set(hr, 'color','k')
        pbaspect([1 1 1])
        xticks(0:20:100);
        yticks(0:20:100);
    end
end

clear h lgtxt
seed = 0;
for i = 1:numel(kappa)
    seed = seed+1;
    h(seed) = plot(-1, -1, 'k.', 'marker', mk{i}, 'markersize', 8, 'linewidth', 1);
    lgtxt{seed} = sprintf('$\\kappa$ = %.2f $\\mu$m/ms', kappa(i));
end

for i = 1:numel(d_mean)
    seed = seed+1;
    h(seed) = plot(-1, -1, 'k-', 'color', cmap(i,:), 'linewidth', 1);
    lgtxt{seed} = sprintf('$d$ = %u $\\mu$m', d_mean(i));
end

seed = seed+1;
h(seed) = fill([-1 0 0 -1],[1 1 2e2 2e2],0.85*[1 1 1],'edgecolor',0.5*[1 1 1]);
lgtxt{seed} = '$\,\, t_r<t_D$';

hl = legend(h, lgtxt, 'box', 'off', 'location', 'northwest', 'fontsize', 12, 'interpreter', 'latex');
hl.Position = hl.Position + [-0.015 0.015 0 0];
hl.ItemTokenSize = [18 18];
