% ********** Setup the directory on your computer **********
clear
restoredefaultpath
filePath = matlab.desktop.editor.getActiveFilename;
root0 = fileparts(filePath);
addpath(genpath(fullfile(root0,'lib')));
root = fullfile(root0,'data');
root_code = fullfile(root0,'lib','packing');

projname = 'sphere_1';

%% Generate random sphere packing

% Sphere diameter (um) and intra-spherical volume fraction
d_mean = [5,    7.5,   10,    15,    20];  % mean diameter, um
f      = 0.1:0.1:0.6;                      % volume fraction

% Membrane permeability, um/ms
kappa = 0:0.01:0.05;

% b-table: [b-value gx gy gz]
bval = 0.5:0.5:3;             % b-value, ms/um^2
bvec = [1 0 0; 0 1 0; 0 0 1]; % [gx gy gz]
btab = [];
for i = 1:numel(bval)
    bvali = bval(i);
    for j = 1:size(bvec,1)
        bvecj = bvec(j,:);
        btab = cat(1,btab,[bvali; bvecj(:)]);
    end
end
Nbvec = numel(bval)*size(bvec,1);

% Simulation parameters
dt = 2e-4;              % time of each step, ms
TN = 1e6;               % # steps
NPar = 1e5;             % # random walkers, usually more than 1e5
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
            maxdensity = fi;    % maximal intra-spherical volume fraction, the packing code will aim for this value.
            hardwallBC = 0;     % boundary condition, 0 for periodic, 1 for hard wall
            N = 999;            % # sphere
            rinit = (di/2)*ones(N,1); % sphere radius (um)
            compileFlag = 0;    % compile Donev's C code, 0 for not compile, 1 for compile
            n = 200;            % matrix size of the lookup table, n x n x n
            gap = 0;            % minimal distance between spheres (um), default=0
            ps.packing(maxdensity,hardwallBC,rinit,root_code,target,compileFlag,n,gap)

            % Save b-table, not really used in this study
            fileID = fopen(fullfile(target,'btable.txt'),'w');
            fprintf(fileID,sprintf('%g\n',btab));
            fclose(fileID);

            % Save simulation parameters
            fileID = fopen(fullfile(target,'simParamInput.txt'),'w');
            fprintf(fileID,sprintf('%g\n%u\n%u\n%u\n%g\n%g\n%g\n%u\n%u\n',dt,TN,NPar,Nbvec,Din,Dex,kappai,pinit,threadpb));
            fclose(fileID);
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
    fprintf(fileID,sprintf('cp -a %s .\n',fullfile(root_cuda,'main_cuda')));
    fprintf(fileID,'./main_cuda\n');
end
fclose(fileID);

% Open the terminal window in the project folder and run "sh job.sh"

% You may need to open the terminal in the root_cuda folder and compile the 
% CUDA code using "nvcc main.cu -o main_cuda"

%% Have a look for the microstructure
% The lookup table A saves two axon labels in one integer. If the first
% and the second axon labels are ax1 and ax2, ax1 = mod(A,Nmax), and ax2 =
% floor(A/Nmax).
% Other parameters:
%   NPix: matrix size of lookup table
%   fov:  field of view of the entire gemoetry and lookup table A, um
%   Nax:  # spheres
%   rCir: sphere radius, normalized by fov
%   [xCir,yCir,zCir]: sphere center position, normalized by fov
i = 1;
root_packing = fullfile(root,projname,sprintf('sphere_%04u',i));

A    = load(fullfile(root_packing,'phantom_APix.txt'));
Nmax = load(fullfile(root_packing,'phantom_Nmax.txt'));
NPix = load(fullfile(root_packing,'phantom_NPix.txt'));
fov  = load(fullfile(root_packing,'phantom_res.txt'));

Nax  = load(fullfile(root_packing,'phantom_NAx.txt'));
rCir = load(fullfile(root_packing,'phantom_rCir.txt'));
xCir = load(fullfile(root_packing,'phantom_xCir.txt'));
yCir = load(fullfile(root_packing,'phantom_yCir.txt'));
zCir = load(fullfile(root_packing,'phantom_zCir.txt'));

%% Plot the microstructure and the lookup table
% Plot the microstructure
ps = packsph();
figure; ps.plotpack(xCir,yCir,zCir,rCir);

% Plot the lookup table
Ar = reshape(A,NPix,NPix,NPix);
figure; ps.plotlookup(Ar,Nax,100);
