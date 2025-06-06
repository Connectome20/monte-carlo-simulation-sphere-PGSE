% ********** Setup the directory on your computer **********
clear
restoredefaultpath
filePath = matlab.desktop.editor.getActiveFilename;
root0 = fileparts(filePath);
addpath(genpath(fullfile(root0,'lib')));
root = fullfile(root0,'data');
root_code = fullfile(root0,'lib','packing');

projname = 'sphere_2';

%% Read Each data set

files = dir(fullfile(root,projname,'sphere_*'));
sim = struct([]);
for i = 1:numel(files)
    filename = sprintf('sphere_%04u',i);
    sim(i).data = simul3Dsphere_cuda(fullfile(root,projname,filename));
end
save(fullfile(root,projname,'simResult_allsetup.mat'),'sim');

%% Analyze packing
files = dir(fullfile(root,projname,'sphere_*'));
pck = struct([]);
for i = 1:numel(files)
    rooti = fullfile(root,projname,sprintf('sphere_%04u',i));

    % field-of-view of the whole geometry (um)
    pck(i).fov = load(fullfile(rooti,'phantom_res.txt'));

    % sphere radius (um)
    pck(i).r   = load(fullfile(rooti,'phantom_rCir.txt'))*pck(i).fov;

    pck(i).rmean = mean(pck(i).r);     % mean sphere radius (um)
    vol = sum(4/3*pi*pck(i).r.^3);     % sphere volume (um^3)
    sur = sum(4*pi*pck(i).r.^2);       % sphere surface area (um^2)
    pck(i).f_in = vol/pck(i).fov^3;    % intra-spherical volume fraction
    pck(i).f_ex = 1-pck(i).f_in;       % extra-spherical volume fraction
    pck(i).sv_in = sur/vol;            % intra-spherical SV ratio, 1/um
    pck(i).sv_ex = sur/(pck(i).fov^3-vol);  % extra-spherical SV ratio, 1/um
    pck(i).K0 = 3*pck(i).f_in/pck(i).f_ex;  % kurtosis at time=0 in Karger model
    pck(i).D_in = sim(i).data.Din;     % intrinsic diffusivity inside sphere, um^2/ms
    pck(i).D_ex = sim(i).data.Dex;     % intrinsic diffusivity outside sphere, um^2/ms
    pck(i).kappa = sim(i).data.kappa;  % membrane permeability, um/ms

    % exchange time, ms
    pck(i).tex = round((1-pck(i).f_in)/(pck(i).kappa * pck(i).sv_in));
end
save(fullfile(root,projname,'packing.mat'),'pck');

%% Calculate gradient direction for MK calculation
ndir = 64;  % # gradient direction

% You need to install mrtrix first
gdir = dirgen(ndir);

% % Alternatively, you can create many random directions to calculate MK
% ndir = 1000;
% gdir = randn(ndir,3);
% gdir = bvec./sqrt(sum(gdir.^2,2));

%% Plot figure: narrow pulse PG, mean kurtosis (MK)

nfb = numel(sim);
t = sim(1).data.TD;    % diffusion time (ms)

Din = [pck.D_in];      % intra-spherical intrinsic diffusivity (um^2/ms)
Dex = [pck.D_ex];      % extra-spherical intrinsic diffusivity (um^2/ms)
kappa = [pck.kappa];   % membrane permeability (um/ms)
kappa = round(kappa,1,'significant');

Dinu = unique(Din);
Dexu = unique(Dex);
kappau = unique(kappa);

f = [pck.f_in];        % intra-spherica volume fraction
f = round(f*100)/100;
r = [pck.rmean];       % sphere radius (um)
r = round(r*100)/100;

karger = @(x) 2./(x).*(1-1./(x).*(1-exp(-x)));

tex = [pck.tex];

% close all
figure('unit','inch','position',[0 0 19.37 8.98]);
cmap = colormap('jet'); cmap = cmap(64:64:end,:);
tpeak = zeros(numel(Din),1);
tex_fit = zeros(numel(f),1);
K0_fit = zeros(numel(f),1);
f_fit = zeros(numel(f),1);
for i = 1:numel(Din)
    It = find(Din(i)==Dinu);
    Ie = find(Dex(i)==Dexu);
    Ik = find(kappa(i)==kappau);
    
    simi = sim(i).data;
    [Ki,Di] = simi.akc_mom(gdir);
    MDi = mean(Di,2);
    MKi = mean(Ki,2);

    
    subplot(3,4,(Ik-1)*numel(Dexu)+Ie)
    hold on
    plot(t,MKi,'-','color',cmap(It,:),'linewidth',1);
    
    if Ie<=2
        if Ie==1
            [~,Ii] = min(abs(t-10));
        else
            [~,Ii] = min(abs(t-2));
        end
        listi = Ii:numel(t);
        [~,Ii] = max(MKi(listi));
        tpeak(i) = t(listi(Ii));
        plot(tpeak(i),MKi(listi(Ii)),'v','color',cmap(It,:),'linewidth',1);
    else
        [~,Ii] = max(MKi);
        tpeak(i) = t(Ii);
        plot(tpeak(i),MKi(Ii),'v','color',cmap(It,:),'linewidth',1);
    end
    title({sprintf('$\\kappa$ = %.2f $\\mu$m/ms',kappa(i)),sprintf('$D_{0,\\rm ec}$ = %.1f $\\mu$m$^2$/ms',Dex(i))},...
       'interpreter','latex','fontsize',16);
   
   if tpeak(i) < 75
        [~,Ii] = min(abs(t-tpeak(i)*2));
        X = kargerfit(t(Ii:end),MKi(Ii:end),[5,20]);
        K0_fit(i) = X(1);
        f_fit(i) = X(1)/(X(1)+3);
        tex_fit(i) = X(2);
    end
end

subplot(341)
hold on;
clear hd lgtxt
for i = 1:numel(Dinu)
    hd(i) = plot(-1,0,'-','color',cmap(i,:),'linewidth',1);
    lgtxt{i} = sprintf('$D_{0,\\rm ic}$ = %.1f $\\mu$m$^2$/ms',Dinu(i));
end
hd(i+1) = plot(-1,-1,'kv');
lgtxt{i+1} = 'peak';
hl = legend(hd,lgtxt,'interpreter','latex','fontsize',12,...
    'location','northwest','box','off','NumColumns',1);
hl.Position = hl.Position + [-0.005 0.005 0 0];
hl.ItemTokenSize = [18 18];

for i = 1:12
   subplot(3,4,i);
   xlim([0 200]); ylim([0 1.5]);
   box on; grid on;
   if i>8, xlabel('$t$, ms','interpreter','latex','fontsize',20); end
   if mod(i,4)==1, ylabel('$K(t)$','interpreter','latex','fontsize',20); end
end

%% Plot universal scaling function
figure(100);
cmap = colormap('parula');
for i = 1:numel(f)
    Ie = find(Dex(i)==Dexu);
    Ik = find(kappa(i)==kappau);
    
    simi = sim(i).data;
    [Ki,Di] = simi.akc_mom(gdir);
    MDi = mean(Di,2);
    MKi = mean(Ki,2);
    
    fi = f(i);
    ri = r(i);
    ki = kappa(i);
    K0 = 3*fi/(1-fi);
    tr = 1/(ki*3/ri);
    tc = ri^2/Dex(i);
    ci = sqrt(2/15)/(1-fi)*sqrt(tc/tr);

    hold on;
    x = t(:);
    tpeak_tilde = sqrt(6/5*tc*tr) * (1+sqrt(tc/30/tr));
    x = x/tpeak_tilde;
    y = MKi(:);

    y = (y/K0 -1 - ((1-fi)/2*(1./x-x) + (1+2*fi)/4./x.^2 + 3/4.*x.^2 + 3-2*fi)*ci^2 )/ci;
    It = round(tc/tr /3*size(cmap,1));
    It = min(It,size(cmap,1)); It = max(It,1);
    
    plot(x,y,'-','color',cmap(It,:),'linewidth',1);
    
end
hb = colorbar;
hb.Ticks = [0 1/3 2/3 1];
hb.TickLabels = {'0','1','2','3'};
hb.FontSize = 12;
x = 0.01:0.001:30;
y = -(x + 1./x);
plot(x,y,'r.');
hx = xline(1); set(hx,'color','r','linestyle',':','linewidth',2);
hs = plot(-1,-1,'-','color',0.5*[1 1 1],'linewidth',1);
ht = plot(-1,-1,'-r','linewidth',1);
hl = legend([hs,ht,hx],{'simulation','$-(\tilde{\tau}^{-1}+\tilde{\tau})$','$\tilde{\tau}=1$'},...
    'interpreter','latex','fontsize',14,'box','off','location','northwest');
pbaspect([1 1 1]); box on;
set(gca,'xscale','log');
xlim([0.01 100]); ylim([-10 0]);

xlabel('$\tilde{\tau} = t/\tilde{t}_{\rm peak,NP}$','interpreter','latex','fontsize',20);
ylabel('$(K/K_0-1-{\cal O}(c^2))/c$','interpreter','latex','fontsize',20);

%% Theory curve of K_var and h[q]
ri = 5;
fi = 0.5;
ki = 0.1;
texi = (1-fi)/(ki*3/ri);
ti = 0.1:0.5:200;

D0ei = 1;
D0ii = 2;
Dini = sphereDK(ri,D0ii,ti,10);
Dexi = D0ei;
Dmi  = fi*Dini + (1-fi)*Dexi;
Kvar = fi*(1-fi)*(Dini - D0ei).^2 ./ Dmi.^2;
hq = 2*(ti/texi-1+exp(-ti/texi))./(ti/texi).^2;
figure;
cmap = colormap('lines');
hold on;
hk = plot(ti,Kvar.*hq,'-k','color',0.7*[1 1 1],'linewidth',2);
hv = plot(ti,Kvar,'--r','color',cmap(7,:),'linewidth',1);
he = plot(ti,hq,'--b','color',cmap(1,:),'linewidth',1);
legend([hv,he,hk],{'$K_{\rm var}(t)$','$2(r_{\rm ex}t - 1 + e^{-r_{\rm ex}t})/(r_{\rm ex}t)^2$','$K(t)$'},'location','east','box','off',...
    'interpreter','latex','fontsize',14);
pbaspect([1 1 1]); box on;
xticks([]);
yticks([]);
xlabel('diffusion time','fontsize',20);
ylabel('kurtosis','fontsize',20);
title('$D_{\rm 0,ic} > D_{\rm 0,ec}$','interpreter','latex','fontsize',20);

D0ei = 2;
D0ii = 1;
Dini = sphereDK(ri,D0ii,ti,10);
Dexi = D0ei;
Dmi  = fi*Dini + (1-fi)*Dexi;
Kvar = fi*(1-fi)*(Dini - D0ei).^2 ./ Dmi.^2;
hq = 2*(ti/texi-1+exp(-ti/texi))./(ti/texi).^2;
figure;
cmap = colormap('lines');
hold on;
hk = plot(ti,Kvar.*hq,'-k','color',0.7*[1 1 1],'linewidth',2);
hv = plot(ti,Kvar,'--r','color',cmap(7,:),'linewidth',1);
he = plot(ti,hq,'--b','color',cmap(1,:),'linewidth',1);
legend([hv,he,hk],{'$K_{\rm var}(t)$','$2(r_{\rm ex}t - 1 + e^{-r_{\rm ex}t})/(r_{\rm ex}t)^2$','$K(t)$'},'location','east','box','off',...
    'interpreter','latex','fontsize',14);
pbaspect([1 1 1]); box on;
xticks([]);
yticks([]);
xlabel('diffusion time','fontsize',20);
ylabel('kurtosis','fontsize',20);
title('$D_{\rm 0,ic} \leq D_{\rm 0,ec}$','interpreter','latex','fontsize',20);

%% Plot figure, tpeak vs kappa, f, diameter
Dinu = unique(Din);
Dexu = unique(Dex);
kappau = unique(kappa);
nDin = numel(Dinu);
nDex = numel(Dexu);
nkappa = numel(kappau);

y = zeros(nDin, nDex, nkappa);
x = zeros(nDin, nDex, nkappa, 3);
z = zeros(nDin, nDex, nkappa, 3);
for i = 1:nDin
    for j = 1:nDex
        for k = 1:nkappa
            I = find(ismember([Din(:), Dex(:), kappa(:)], [Dinu(i), Dexu(j), kappau(k)],'rows'));
            y(i,j,k) = tpeak(I);
            x(i,j,k,:) = [Din(I), Dex(I), kappa(I)];
            z(i,j,k,:) = [K0_fit(I), f_fit(I), tex_fit(I)];
        end
    end
end

aDin = zeros(nkappa, nDex);
aDin_std = zeros(nkappa, nDex);
aDex = zeros(nkappa, nDin);
aDex_std = zeros(nkappa, nDin);

ri = 5;
tr = 1./(x(:,:,:,3)*3/ri);
tD = ri^2./x(:,:,:,1);

figure('unit','inch','position',[0 0 19.37/6*4 8.98/3*2]);
cmap = colormap('lines');
mk = {'o','v','d','^','s'};
% tpeak vs Din
clear h1 lgtxt1 h2 lgtxt2
for i = 1:nkappa % kappa
    subplot(2,4,i);
    hold on;
    plot([0.5 2],1e2*[1 (2/0.5)^(-0.5)],'--k','linewidth',1);
    if i==nkappa
        ht = text(sum([0.5 1.5])/2, sum(1e2*[1 (1.5/0.5)^(-0.5)])/2, 'slope=$-0.5$','interpreter','latex',...
            'VerticalAlignment','bottom','HorizontalAlignment','center',...
            'FontSize',16);
        set(ht,'Rotation',-7);
    end
    for j = 1:nDin
        tri = tr(j,:,i); tri = tri(:);
        tDi = tD(j,:,i); tDi = tDi(:);
        list = tri - tDi > -1e-10;
        
        yi = y(j,:,i); yi = yi(:);
        xi = Dexu(:);
        plot(xi(list),yi(list),mk{j},'color',cmap(j,:),'linewidth',1);
        plot(xi(~list),yi(~list),mk{j},'color',0.5*[1 1 1],'linewidth',1);
        
        if nnz(list) > 1
            xi = xi(list);
            yi = yi(list);
            A = [ones(length(xi),1), log(xi)];
            [X, STDX] = lscov(A,log(yi),1./yi);
            xi = [0.5 2].'; A = [ones(length(xi),1), log(xi)];
            Y = A*X;
            aDex(i,j) = X(2);
            aDex_std(i,j) = STDX(2);
            plot(xi,exp(Y),':','color',cmap(j,:),'linewidth',1);
        end
        xlim([0.5 2]);
        ylim([1 200]);
        set(gca,'yscale','log','xscale','log')
        box on; grid on;
        xlabel('$D_{0,\rm ec}$, $\mu$m$^2$/ms','interpreter','latex','fontsize',16);
        ylabel('$t_{\rm peak}$, ms','interpreter','latex','FontSize',16);
        title(sprintf('$\\kappa$=%.2f $\\mu$m/ms',kappau(i)),'interpreter','latex','fontsize',16);
        pbaspect([1 1 1]);
    end
    
    subplot(2,4,i+4)
    hold on;
    fill([0.5 kappau(i)*3*ri kappau(i)*3*ri 0.5],[1 1 2e2 2e2],0.92*[1 1 1],'edgecolor','w');
    set(gca, "Layer", "top")
    plot([0.5 2],1e2*[1 (2/0.5)^(-0)],'--k','linewidth',1);
    if i==nkappa
        ht = text(sum([0.5 1.5])/2, sum(1e2*[1 (1.5/0.5)^(-0)])/2, 'slope=$0$','interpreter','latex',...
            'VerticalAlignment','top','HorizontalAlignment','center',...
            'FontSize',16);
        set(ht,'Rotation',0);
    end
    for j = 1:nDex
        tri = tr(:,j,i); tri = tri(:);
        tDi = tD(:,j,i); tDi = tDi(:);
        list = tri - tDi > -1e-10;
        
        yi = y(:,j,i); yi = yi(:);
        xi = Dinu(:);
        plot(xi(list),yi(list),mk{j},'color',cmap(j,:),'linewidth',1);
        plot(xi(~list),yi(~list),mk{j},'color',0.5*[1 1 1],'linewidth',1);
        
        if nnz(list) > 1
            xi = xi(list);
            yi = yi(list);
            A = [ones(length(xi),1), log(xi)];
            [X, STDX] = lscov(A,log(yi),1./yi);
            xi = [0.5 2].'; A = [ones(length(xi),1), log(xi)];
            Y = A*X;
            aDin(i,j) = X(2);
            aDin_std(i,j) = STDX(2);
            plot(xi,exp(Y),':','color',cmap(j,:),'linewidth',1);
        end
        xlim([0.5 2]);
        ylim([1 200]);
        set(gca,'yscale','log','xscale','log')
        box on; grid on;
        xlabel('$D_{0,\rm ic}$, $\mu$m$^2$/ms','interpreter','latex','fontsize',16);
        ylabel('$t_{\rm peak}$, ms','interpreter','latex','FontSize',16);
        title(sprintf('$\\kappa$=%.2f $\\mu$m/ms',kappau(i)),'interpreter','latex','fontsize',16);
        pbaspect([1 1 1]);
    end
    
end

subplot(2,4,4);
hold on;
for j = 1:nDin
    h1(j) = plot(-1,-1,mk{j},'color',cmap(j,:),'linewidth',1);
    lgtxt1{j} = sprintf('$D_{0,\\rm ic}=%.1f$ $\\mu$m$^2$/ms',Dinu(j));
end
h1(nDin+1) = plot(-1,-1,'-','linewidth',1,'color',0.5*[1 1 1]);
lgtxt1{nDex+1} = '$\quad t_r<t_D$';
lg1 = legend(h1, lgtxt1,'interpreter','latex','fontsize',16,'location','west','box','off');
lg1.Position = lg1.Position + [-0.025 0 0 0];
xlim([0 1]); ylim([0 1]);
axis off

subplot(2,4,8);
hold on;
for j = 1:nDex
    h2(j) = plot(-1,-1,mk{j},'color',cmap(j,:),'linewidth',1);
    lgtxt2{j} = sprintf('$D_{0,\\rm ec}=%.1f$ $\\mu$m$^2$/ms',Dexu(j));
end
h2(nDex+1) = fill([-1 0 0 -1],[1 1 2e2 2e2],0.92*[1 1 1],'edgecolor',0.5*[1 1 1]);
lgtxt2{nDex+1} = '$\quad t_r<t_D$';
lg2 = legend(h2, lgtxt2,'interpreter','latex','fontsize',16,'location','west','box','off');
lg2.Position = lg2.Position + [-0.025 0 0 0];
xlim([0 1]); ylim([0 1]);
axis off

%% report power law exponent
clc
fprintf('alpha_Dex = %.2f +- %.2f\n',mean(aDex(aDex~=0),'all'),rms(aDex_std(aDex~=0),'all'));
fprintf('alpha_Din = %.2f +- %.2f\n',mean(aDin(:)),rms(aDin_std(:)));

%% Simulation vs theory for tpeak

ri = 5;
f = 0.4;
tr = 1./(x(:,:,:,3)*3/ri);
tD = ri^2./x(:,:,:,1);
tc = ri^2./x(:,:,:,2);

t_theory = sqrt(6/5.*tr.*tc);
figure('unit','inch','position',[0 0 19.37/4*3 8.98/2]);
cmap = colormap('lines');
for i = 1:nkappa
    subplot(1,nkappa,i);
    hold on;
    clear h lgtxt
    for j = 1:nDin
        for k = 1:nDex
            yi = y(j,k,i);
            xi = t_theory(j,k,i);
            if tr(j,k,i) - tD(j,k,i) <= -1e-10
                plot(yi,xi,mk{j},'color',cmap(k,:),'markerfacecolor',0.85*[1 1 1],'markersize',8,...
                'linewidth',1);
            else
                plot(yi,xi,mk{j},'color',cmap(k,:),'markersize',8,...
                    'linewidth',1);
            end
        end
    end
    xlim([0 100]);
    ylim([0 100]);
    xticks(0:20:100);
    yticks(0:20:100);
    box on; grid on;
    ylabel('$\sqrt{\frac{6}{5} \, t_c \, t_r}$, ms','interpreter','latex','fontsize',16);
    xlabel('$t_{\rm peak}$, ms','interpreter','latex','FontSize',16);
    title(sprintf('$\\kappa$ = %.2f $\\mu$m/ms',kappau(i)),'interpreter','latex','fontsize',16);
    hr = refline(1,0); set(hr,'color','k');
    pbaspect([1 1 1]);
    if i==3
        l = 0;
        for j = 1:nDin
            l = l+1;
            h(l) = plot(-1,-1,mk{j},'color','k','markersize',8,'linewidth',1);
            lgtxt{l} = sprintf('$D_{0,\\rm ic}$=%.1f',Dinu(j));
        end
        for j = 1:nDex
            l = l+1;
            h(l) = plot(-1,-1,'-','color',cmap(j,:),'markersize',8,'linewidth',1);
            lgtxt{l} = sprintf('$D_{0,\\rm ec}$=%.1f',Dexu(j));
        end
        l = l+1;
        h(l) = fill([-1 0 0 -1],[1 1 2e2 2e2],0.85*[1 1 1],'edgecolor',0.5*[1 1 1]);
        lgtxt{l} = '$\,\, t_r<t_D$';
        lg = legend(h, lgtxt, 'interpreter','latex','fontsize',12,'location','northwest','box','off');
        lg.Position = lg.Position + [-0.01 0.015 0 0];
        lg.ItemTokenSize = [18 18];
    end
end

%% Karger model fitting vs ground truth for f, tex, tc

% y: tpeak, nDin x nDex x nkappa
% x: Din, Dex, kappa
% z: K0_fit, f_fit, tex_fit

fi = 0.4;
ri = 5;

Dini = x(:,:,:,1);
Dexi = x(:,:,:,2);
ki   = x(:,:,:,3);

Kfi = z(:,:,:,1);
ffi = z(:,:,:,2);
tfi = z(:,:,:,3);

tc = ri.^2./Dexi;
tr = 1./(ki*3./ri);

t_theory = sqrt(6/5.*tr.*tc);

figure('unit','inch','position',[0 0 19.37/6*5 8.98]);
cmap = colormap('lines');
for i = 1:nDin % Din
    subplot(3,nDin+1,i)
    hold on;
    for j = 1:nDex % Dex
        for k = 1:nkappa % kappa
            yi = squeeze(ffi(i,j,k));
            xi = fi;
            
            if yi
                if tr(i,j,k) - tD(i,j,k) <= -1e-10
                    plot(xi, yi, mk{j}, 'color', cmap(k,:), 'markerfacecolor',0.85*[1 1 1],...
                        'markersize', 8, 'linewidth', 1);
                else
                    plot(xi, yi, mk{j}, 'color', cmap(k,:), 'markersize', 8, 'linewidth', 1);
                end
            end
        end
    end
    xlim([0 0.8]);
    ylim([0 0.8]);
    xticks(0:0.2:1);
    yticks(0:0.2:1);
    box on; grid on;
    xlabel('$f$, ground truth','interpreter','latex','fontsize',14);
    ylabel('$f$, fitted','interpreter','latex','FontSize',14);
    title(sprintf('$D_{\\rm 0,ic}=%.2g$ $\\mu$m$^2$/ms',Dinu(i)),'interpreter','latex','fontsize',14)
    hr = refline(1,0); set(hr,'color','k');
    pbaspect([1 1 1]);
    
    subplot(3,nDin+1,i+nDin+1)
    hold on;
    for j = 1:nDex % Dex
        for k = 1:nkappa % kappa
            yi = squeeze(tfi(i,j,k));
            xi = squeeze((1-fi).*tr(i,j,k));
            if yi
                if tr(i,j,k) - tD(i,j,k) <= -1e-10
                    plot(xi, yi, mk{j}, 'color', cmap(k,:), 'markerfacecolor',0.85*[1 1 1],...
                        'markersize', 8, 'linewidth', 1);
                else
                    plot(xi, yi, mk{j}, 'color', cmap(k,:), 'markersize', 8, 'linewidth', 1);
                end
            end
        end
    end
    xlim([0 100]);
    ylim([0 100]);
    xticks(0:20:100)
    box on; grid on;
    xlabel('$t_{\rm ex}$, ms, ground truth','interpreter','latex','fontsize',14);
    ylabel('$t_{\rm ex}$, ms, fitted','interpreter','latex','FontSize',14);
    title(sprintf('$D_{\\rm 0,ic}=%.2g$ $\\mu$m$^2$/ms',Dinu(i)),'interpreter','latex','fontsize',14)
    hr = refline(1,0); set(hr,'color','k');
    pbaspect([1 1 1]);
    
    subplot(3,nDin+1,i+(nDin+1)*2)
    hold on;
    for j = 1:nDex % Dex
        for k = 1:nkappa % kappa
            yi = squeeze(y(i,j,k)^2/( 6/5*tfi(i,j,k) ./ (1-ffi(i,j,k)) ) );
            xi = squeeze(tc(i,j,k));
            if yi
                if tr(i,j,k) - tD(i,j,k) <= -1e-10
                    plot(xi, yi, mk{j}, 'color', cmap(k,:), 'markerfacecolor',0.85*[1 1 1],...
                        'markersize', 8, 'linewidth', 1);
                else
                    plot(xi, yi, mk{j}, 'color', cmap(k,:), 'markersize', 8, 'linewidth', 1);
                end
            end
        end
    end
    xlim([0 60]);
    ylim([0 60]);
    yticks(0:20:200);
    box on; grid on;
    xlabel('$t_c$, ms, ground truth','interpreter','latex','fontsize',14);
    ylabel('$t_c$, ms, fitted','interpreter','latex','FontSize',14);
    title(sprintf('$D_{\\rm 0,ic}=%.2g$ $\\mu$m$^2$/ms',Dinu(i)),'interpreter','latex','fontsize',14)
    hr = refline(1,0); set(hr,'color','k');
    pbaspect([1 1 1]);
    
    
end

clear h lgtxt
for i = 1:3
    subplot(3,nDin+1,(nDin+1)*i);
    hold on;
    for j = 1:nDex
        h(j) = plot(-1,-1,mk{j},'color','k','markersize',8,'linewidth',1);
        lgtxt{j} = sprintf('$D_{\\rm 0,ec}=%.1f$ $\\mu$m$^2$/ms',Dexu(j));
    end
    for k = 1:nkappa
        h(k+nDex) = plot(-1,-1,'-','color',cmap(k,:),'linewidth',1);
        lgtxt{k+nDex} = sprintf('$\\kappa=%.2f$ $\\mu$m/ms',kappau(k));
    end
    h(nDex+nkappa+1) = fill([-1 0 0 -1],[1 1 2e2 2e2],0.85*[1 1 1],'edgecolor',0.5*[1 1 1]);
    lgtxt{nDex+nkappa+1} = '$\,\, t_r<t_D$';
    xlim([0 1]);
    ylim([0 1]);
    
    hl = legend(h,lgtxt,'interpreter','latex','fontsize',14,'location','west',...
        'box','off');
    hl.Position = hl.Position + [-0.02/6*5 0.01 0 0];
    hl.ItemTokenSize = [18 18];
    axis off
end
