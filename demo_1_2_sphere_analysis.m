% ********** Setup the directory on your computer **********
clear
restoredefaultpath
filePath = matlab.desktop.editor.getActiveFilename;
root0 = fileparts(filePath);
addpath(genpath(fullfile(root0,'lib')));
root = fullfile(root0,'data');
root_code = fullfile(root0,'lib','packing');

projname = 'sphere_1';

%% Read data

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

f = [pck.f_in];        % intra-spherical volume fraction
f = round(f*100)/100;
r = [pck.rmean];       % sphere radius (um)
r = round(r*100)/100;

[C, IA, IC] = unique([r(:),f(:)],'row');
ru = C(:,1); fu = C(:,2);
tex = [pck.tex];       % exchange time (ms)
kappa = [pck.kappa];   % membrane permeability (um/ms)
ku = unique(kappa); nt = numel(ku);
Din = [pck.D_in];      % intra-spherical intrinsic diffusivity (um^2/ms)
Dex = [pck.D_ex];      % extra-spherical intrinsic diffusivity (um^2/ms)

f = f(:); r = r(:); kappa = kappa(:);

karger = @(x) 2./(x).*(1-1./(x).*(1-exp(-x)));

nrow = 3;
ncol = 10;

figure('unit','inch','position',[0 0 19.37 8.98]);
cmap = colormap('jet'); cmap = cmap(23:23:end,:);
tpeak = zeros(numel(f),1);
tex_fit = zeros(numel(f),1);
K0_fit = zeros(numel(f),1);
f_fit = zeros(numel(f),1);
for i = 1:numel(f)
    Ir = find((r(i)==ru).*(f(i)==fu));
    It = find(kappa(i)==ku);
    simi = sim(i).data;
    [Ki,Di] = simi.akc_mom(gdir);
    MDi = mean(Di,2);
    MKi = mean(Ki,2);

    [~,Ii] = max(MKi);
    tpeak(i) = t(Ii);
    
    subplot(nrow,ncol,Ir)
    hold on;
    plot(t(1:10:end),MKi(1:10:end),'-','color',cmap(It,:),'linewidth',1);
    if sim(i).data.kappa > 0
        plot(tpeak(i),MKi(Ii),'v','color',cmap(It,:),'linewidth',1);
    end
    
    if tpeak(i) < 75
        [~,Ii] = min(abs(t-tpeak(i)*2));
        X = kargerfit(t(Ii:end),MKi(Ii:end),[5,20]);
        K0_fit(i) = X(1);
        f_fit(i) = X(1)/(X(1)+3);
        tex_fit(i) = X(2);
    end
end

subplot(nrow,ncol,1)
hold on;
clear hd lgtxt
for i = 1:nt
    hd(i) = plot(-1,-1,'-','color',cmap(i,:),'linewidth',1);
    lgtxt{i} = sprintf('$\\kappa$ = %.2f',kappa(i));
end
lg = legend(hd,lgtxt,'interpreter','latex','fontsize',12,...
    'location','north','box','off','NumColumns',1);
lg.Position = lg.Position + [0 1e-2 0 0];
lg.ItemTokenSize = [18 18];

for i = 1:numel(IA)
   subplot(nrow,ncol,i);
   xlim([0 200]);
   if i <= 10
       ylim([0 1]);
   elseif i <= 20
       ylim([0 2]);
   else
       ylim([0 3]);
   end
   box on; grid on;
   set(gca,'fontsize',10,'xtick',0:50:200);
   if i == 1,  yticks(0:0.2:0.8); end
   if i == 11, yticks(0:0.5:1.5); end
   if i == 21, yticks(0:0.5:2.5); end
   if i > 20, xlabel('$t$, ms','interpreter','latex','fontsize',16); end
   if mod(i,10)==1, ylabel('$K(t)$','interpreter','latex','fontsize',16); end
   ht = title(sprintf('($%.2g$ $\\mu$m, $%.1f$)',2*ru(i),fu(i)),...
       'interpreter','latex','fontsize',12);
   if mod(i,10) == 1
       xi = ht.Position(1);
       yi = ht.Position(2);
       text(xi-80,yi,'($d$, $f$)=','interpreter','latex','fontsize',12,...
           'HorizontalAlignment','right','VerticalAlignment','bottom');
   end
end

%% Plot universal scaling function
figure(100);
cmap = colormap('parula');
for i = 1:numel(f)
    Ir = find((r(i)==ru).*(f(i)==fu));
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
pbaspect([1 1 1]);
box on;
set(gca,'xscale','log');
xlim([0.01 100]);
ylim([-10 0]);
xlabel('$\tilde{\tau} = t/\tilde{t}_{\rm peak,NP}$','interpreter','latex','fontsize',20);
ylabel('$\frac{c}{1+{\cal O}(c^2)-K(\tau)/K_0}$','interpreter','latex','fontsize',24);
ylabel('$(K/K_0-1-{\cal O}(c^2))/c$','interpreter','latex','fontsize',20);

%% Plot figure, tpeak vs kappa, f, diameter
ku = unique(kappa);
fu = unique(f);
ru = unique(r);

nk = numel(ku);
nf = numel(fu);
nr = numel(ru);

y = zeros(nk, nf, nr);
x = zeros(nk, nf, nr, 3);
z = zeros(nk, nf, nr, 3);
for i = 1:nk
    for j = 1:nf
        for k = 1:nr
            I = find( (kappa==ku(i)) .* (f==fu(j)) .* (r==ru(k)) );
            y(i,j,k) = tpeak(I);
            x(i,j,k,:) = [kappa(I), f(I), r(I)];
            z(i,j,k,:) = [K0_fit(I), f_fit(I), tex_fit(I)];
        end
    end
end

Din = sim(1).data.Din;
tr = 1./(x(:,:,:,1)*3./x(:,:,:,3));
tD = x(:,:,:,3).^2./Din;

af = zeros(nk,nr);
ak = zeros(nf,nr);
ad = zeros(nk,nf);
af_std = zeros(nk,nr);
ak_std = zeros(nf,nr);
ad_std = zeros(nk,nf);

figure('unit','inch','position',[0 0 19.37*4/6 8.98]);
cmap = colormap('lines');
mk = {'o','v','d','^','s','>'};
% tpeak vs f
for i = 1:2:nr % radius
    subplot(3,4,(i+1)/2)
    plot([0.1 0.6],10*[1 (0.25/0.21)^(0)],'--k','linewidth',1);
    if i==nr
        ht = text(sum([0.1 0.4])/2, sum(10*[1 (0.4/0.1)^(0)])/2, 'slope=$0$','interpreter','latex',...
            'VerticalAlignment','bottom','HorizontalAlignment','center',...
            'FontSize',16);
        set(ht,'Rotation',1);
    end
    hold on;
    k = 0;
    for j = 2:nk % kappa
        k = k+1;
        yi = squeeze(y(j,:,i)); yi = yi(:);
        xi = fu;%.*(1-fu);
        
        tri = tr(j,:,i); tri = tri(:);
        tDi = tD(j,:,i); tDi = tDi(:);
        list = tri - tDi > -1e-10;
        
        plot(xi(list),yi(list),mk{k},'color',cmap(k,:),'linewidth',1);
        plot(xi(~list),yi(~list),mk{k},'color',0.5*[1 1 1],'linewidth',1);
        
        if nnz(list) > 1
            xi = xi(list);
            yi = yi(list);
            A = [ones(length(xi),1), log(xi)];
            [X, STDX] = lscov(A,log(yi),1./yi);
            xi = [0.1 0.5].'; A = [ones(length(xi),1), log(xi)];
            Y = A*X;
            af(j,i) = X(2);
            af_std(j,i) = STDX(2);
            plot(xi,exp(Y),':','color',cmap(k,:),'linewidth',1);
        end
    end
    xlim([0.1 0.6]);
    ylim([5 200]);
    xticks(0.1:0.1:1);
    yticks([10 100]);
    set(gca,'yscale','log','xscale','log');
    box on; grid on;
    xlabel('$f$','interpreter','latex','fontsize',16);
    ylabel('$t_{\rm peak}$, ms','interpreter','latex','FontSize',16);
    title(sprintf('$d=%.2g$ $\\mu$m',2*ru(i)),'interpreter','latex','fontsize',16);
    pbaspect([1 1 1]);
end
subplot(3,4,4);
hold on;
clear h lgtxt
k = 0;
for j = 2:nk
    k = k+1;
    h(k) = plot(-1,-1,mk{k},'color',cmap(k,:),'linewidth',1);
    lgtxt{k} = sprintf('$\\kappa = %.2f$ $\\mu$m/ms',ku(j));
end
k = k+1;
h(k) = plot(-1,-1,'-','linewidth',1,'color',0.5*[1 1 1]);
lgtxt{k} = '$\quad t_r<t_D$';
xlim([0 1]); ylim([0 1]); axis off
hl = legend(h, lgtxt,'interpreter','latex','fontsize',16,'location','west','box','off');
hl.Position = hl.Position + [-0.025 0 0 0];

% tpeak vs kappa
for i = 1:2:nr % radius
    subplot(3,4,(i+1)/2+4)
    plot([0.01 0.05],100*[1 5^(-0.5)],'--k','linewidth',1);
    if i==nr
        ht = text(sum([0.01 0.035])/2, sum(40*[1 (0.035/0.01)^(-0.5)])/2, 'slope=$-0.5$','interpreter','latex',...
            'VerticalAlignment','bottom','HorizontalAlignment','center',...
            'FontSize',16);
        set(ht,'Rotation',-10);
    end
    hold on;
    k = 0;
    for j = 1:nf % f
        k = k+1;
        yi = squeeze(y(2:6,j,i));
        xi = ku(2:6);
        
        tri = tr(2:6,j,i); tri = tri(:);
        tDi = tD(2:6,j,i); tDi = tDi(:);
        list = tri - tDi > -1e-10;
        
        plot(xi(list),yi(list),mk{k},'color',cmap(k,:),'linewidth',1);
        plot(xi(~list),yi(~list),mk{k},'color',0.5*[1 1 1],'linewidth',1);
        
        if nnz(list) > 1
            xi = xi(list);
            yi = yi(list);
            A = [ones(length(xi),1), log(xi)];
            [X, STDX] = lscov(A,log(yi),1./yi);
            xi = [0.01 0.05].'; A = [ones(length(xi),1), log(xi)];
            Y = A*X;
            ak(j,i) = X(2);
            ak_std(j,i) = STDX(2);
            plot(xi,exp(Y),':','color',cmap(k,:),'linewidth',1);
        end
    end
    ylim([5 200]);
    xtickangle(45);
    set(gca,'yscale','log','xscale','log')
    box on; grid on;
    xlabel('$\kappa$, $\mu$m/ms','interpreter','latex','fontsize',16);
    ylabel('$t_{\rm peak}$, ms','interpreter','latex','FontSize',16);
    title(sprintf('$d=%.2g$ $\\mu$m',2*ru(i)),'interpreter','latex','fontsize',16);
    pbaspect([1 1 1]);
end
subplot(3,4,8);
hold on;
clear h lgtxt
k = 0;
for j = 1:nf
    k = k+1;
    h(k) = plot(-1,-1,mk{k},'linewidth',1);
    lgtxt{k} = sprintf('$f = %.1f$',fu(j));
end
k = k+1;
h(k) = plot(-1,-1,'-','linewidth',1,'color',0.5*[1 1 1]);
lgtxt{k} = '$\quad t_r<t_D$';
xlim([0 1]); ylim([0 1]); axis off
hl = legend(h, lgtxt,'interpreter','latex','fontsize',16,'location','west','box','off');
hl.Position = hl.Position + [-0.025 0 0 0];

% tpeak vs diameter
for i = 2:2:nk % f
    subplot(3,4,i/2+8)
    plot([5 20],10*[1 (20/5)^1.5],'--k','linewidth',1);
    if i==nk
        ht = text(sum([5 15])/2, sum(10*[1 (15/5)^(1.5)])/2, 'slope=$1.5$','interpreter','latex',...
            'VerticalAlignment','bottom','HorizontalAlignment','center',...
            'FontSize',16);
        set(ht,'Rotation',30);
    end
    hold on;
    k = 0;
    for j = 1:nf % kappa
        k = k+1;
        yi = squeeze(y(i,j,:));
        xi = 2*ru;
        
        tri = tr(i,j,:); tri = tri(:);
        tDi = tD(i,j,:); tDi = tDi(:);
        list = tri - tDi > -1e-10;
        
        plot(xi(list),yi(list),mk{k},'color',cmap(k,:),'linewidth',1);
        plot(xi(~list),yi(~list),mk{k},'color',0.5*[1 1 1],'linewidth',1);
        
        if nnz(list) > 1
            xi = xi(list);
            yi = yi(list);
            A = [ones(length(xi),1), log(xi)];
            [X, STDX] = lscov(A,log(yi),1./yi);
            xi = [5 20].'; A = [ones(length(xi),1), log(xi)];
            Y = A*X;
            ad(i,j) = X(2);
            ad_std(i,j) = STDX(2);
            plot(xi,exp(Y),':','color',cmap(k,:),'linewidth',1);
        end
    end
    ylim([5 200]);
    set(gca,'yscale','log','xscale','log')
    box on; grid on;
    xlabel('$d$, $\mu$m','interpreter','latex','fontsize',16);
    ylabel('$t_{\rm peak}$, ms','interpreter','latex','FontSize',16);
    title(sprintf('$\\kappa = %.2f$ $\\mu$m/ms',ku(i)),'interpreter','latex','fontsize',16)
    pbaspect([1 1 1]);
end
subplot(3,4,12);
hold on;
clear h lgtxt
k = 0;
for j = 1:nf
    k = k+1;
    h(k) = plot(-1,-1,mk{k},'linewidth',1);
    lgtxt{k} = sprintf('$f = %.1f$',fu(j));
end
k = k+1;
h(k) = plot(-1,-1,'-','linewidth',1,'color',0.5*[1 1 1]);
lgtxt{k} = '$\quad t_r<t_D$';
xlim([0 1]); ylim([0 1]); axis off
hl = legend(h, lgtxt,'interpreter','latex','fontsize',16,'location','west','box','off');
hl.Position = hl.Position + [-0.025 0 0 0];

%% report power law exponent
clc
fprintf('alpha_f = %.2f +- %.2f\n',mean(af(af~=0),'all'),rms(af_std(af~=0),'all'));
fprintf('alpha_k = %.2f +- %.2f\n',mean(ak(:)),rms(ak_std(:)));
fprintf('alpha_d = %.2f +- %.2f\n',mean(ad(ad~=0),'all'),rms(ad_std(ad~=0),'all'));

%% Simulation vs theory for tpeak
Dex = sim(1).data.Dex;

ki = x(:,:,:,1);
ri = x(:,:,:,3);

tc = ri.^2/Dex;
tr = 1./(ki*3./ri);

t_theory = sqrt(6/5.*tr.*tc);

figure('unit','inch','position',[0 0 19.37/4*3 8.98]);
cmap = colormap('lines');
for i = 1:nf % f
    subplot(2,ceil(nf/2),i)
    hold on;
    for j = 2:nk % kappa
        for k = 1:nr % diameter
            yi = squeeze(y(j,i,k));
            xi = squeeze(t_theory(j,i,k));
            
            if tr(j,i,k) - tD(j,i,k) <= -1e-10
                plot(yi,xi,mk{j-1},'color',cmap(k,:),'markerfacecolor',0.85*[1 1 1],'markersize',8,...
                'linewidth',1);
            else
                plot(yi,xi,mk{j-1},'color',cmap(k,:),'markersize',8,...
                    'linewidth',1);
            end
%             plot(xi,yi,mk{j-1},'color',cmap(k,:),'markersize',8,...
%                 'linewidth',1);
        end
    end
    xlim([0 120]);
    ylim([0 120]);
    box on; grid on;
    ylabel('$\sqrt{\frac{6}{5}\, t_c \, t_r}$, ms','interpreter','latex','fontsize',16);
    xlabel('$t_{\rm peak}$, ms','interpreter','latex','FontSize',16);
    title(sprintf('$f=%.1f$',fu(i)),'interpreter','latex','fontsize',16)
    hr = refline(1,0); set(hr,'color','k');
    pbaspect([1 1 1]);
end

clear h lgtxt
for i = 6
    k = 0;
    subplot(2,ceil(nf/2),i)
    for j = 2:nk
        k = k+1;
        h(k) = plot(-1,-1,mk{j-1},'color','k','markersize',8);
        lgtxt{k} = sprintf('$\\kappa=%.2f$ $\\mu$m/ms',ku(j));
    end
    
    for j = 1:nr
        k = k+1;
        h(k) = plot(-1,-1,'-','color',cmap(j,:),'linewidth',1);
        if j==2
            lgtxt{k} = sprintf('$d=%.1f$ $\\mu$m',2*ru(j));
        else
            lgtxt{k} = sprintf('$d=%u$ $\\mu$m',2*ru(j));
        end
    end
    k = k+1;
    h(k) = fill([-1 0 0 -1],[1 1 2e2 2e2],0.85*[1 1 1],'edgecolor',0.5*[1 1 1]);
    lgtxt{k} = '$\,\, t_r<t_D$';
    
    
    hl = legend(h,lgtxt,'interpreter','latex','fontsize',12,'location','northwest',...
        'box','off');
    hl.Position = hl.Position + [-0.015 0.015 0 0];
    hl.ItemTokenSize = [18 18];
    
end

%% Karger model fitting vs ground truth for f, tex, tc

% y: tpeak, nk x nf x nr
% x: kappa, f, r
% z: K0_fit, f_fit, tex_fit

Din = sim(1).data.Din;
Dex = sim(1).data.Dex;

ki = x(:,:,:,1);
fi = x(:,:,:,2);
ri = x(:,:,:,3);

Kfi = z(:,:,:,1);
ffi = z(:,:,:,2);
tfi = z(:,:,:,3);

tc = ri.^2/Dex;
tr = 1./(ki*3./ri);

t_theory = sqrt(6/5.*tr.*tc);

figure('unit','inch','position',[0 0 19.37*4/6 8.98]);
cmap = colormap('lines');
for i = 2:2:nk % kappa
    subplot(3,4,i/2)
    hold on;
    for j = 1:nr % diameter
        for k = 1:nf % f
            yi = squeeze(ffi(i,k,j));
            xi = squeeze(fi(i,k,j));
            
            if yi
                if tr(i,k,j) - tD(i,k,j) <= -1e-10
                    plot(xi, yi, mk{j}, 'color', cmap(j,:), 'markerfacecolor',0.85*[1 1 1],...
                        'markersize', 8, 'linewidth', 1);
                else
                    plot(xi, yi, mk{j}, 'color', cmap(j,:), 'markersize', 8, 'linewidth', 1);
                end
            end
        end
    end
    xlim([0 0.8]);
    ylim([0 0.8]);
    box on; grid on;
    xlabel('$f$, ground truth','interpreter','latex','fontsize',14);
    ylabel('$f$, fitted','interpreter','latex','FontSize',14);
    title(sprintf('$\\kappa=%.1g$ $\\mu$m/ms',ku(i)),'interpreter','latex','fontsize',14)
    hr = refline(1,0); set(hr,'color','k');
    pbaspect([1 1 1]);
    
    subplot(3,4,i/2+4)
    hold on;
    for j = 1:nr % diameter
        for k = 1:nf % f
            yi = squeeze(tfi(i,k,j));
            xi = squeeze((1-fi(i,k,j)).*tr(i,k,j));
            if yi
                if tr(i,k,j) - tD(i,k,j) <= -1e-10
                    plot(xi, yi, mk{j}, 'color', cmap(k,:), 'markerfacecolor',0.85*[1 1 1],...
                        'markersize', 8, 'linewidth', 1);
                else
                    plot(xi, yi, mk{j}, 'color', cmap(k,:), 'markersize', 8, 'linewidth', 1);
                end
            end
        end
    end
    xlim([0 200]); ylim([0 200]);
    xticks(0:50:200); yticks(0:50:200);
    box on; grid on;
    xlabel('$t_{\rm ex}$, ms, ground truth','interpreter','latex','fontsize',14);
    ylabel('$t_{\rm ex}$, ms, fitted','interpreter','latex','FontSize',14);
    title(sprintf('$\\kappa=%.1g$ $\\mu$m/ms',ku(i)),'interpreter','latex','fontsize',14)
    hr = refline(1,0); set(hr,'color','k');
    pbaspect([1 1 1]);
    
    subplot(3,4,i/2+8)
    hold on;
    for j = 1:nr % diameter
        for k = 1:nf % f
            yi = squeeze( y(i,k,j)^2/( 6/5*tfi(i,k,j)./( 1-ffi(i,k,j) ) ) );
            xi = squeeze(tc(i,k,j));
            if yi
                if tr(i,k,j) - tD(i,k,j) <= -1e-10
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
    yticks(0:20:100);
    box on; grid on;
    xlabel('$t_c$, ms, ground truth','interpreter','latex','fontsize',14);
    ylabel('$t_c$, ms, fitted','interpreter','latex','FontSize',14);
    title(sprintf('$\\kappa=%.1g$ $\\mu$m/ms',ku(i)),'interpreter','latex','fontsize',14)
    hr = refline(1,0); set(hr,'color','k');
    pbaspect([1 1 1]);
    
    
end

clear h lgtxt
for i = 1
    subplot(3,4,4);
    hold on;
    for j = 1:nr
        h(j) = plot(-1,-1,mk{j},'color',cmap(j,:),'markersize',8,'linewidth',1);
        lgtxt{j} = sprintf('$d=%.2g$ $\\mu$m',2*ru(j));
    end
    h(nr+1) = fill([-1 0 0 -1],[1 1 2e2 2e2],0.85*[1 1 1],'edgecolor',0.5*[1 1 1]);
    lgtxt{nr+1} = '$\,\, t_r<t_D$';
    xlim([0 1]);
    ylim([0 1]);
    
    hl = legend(h,lgtxt,'interpreter','latex','fontsize',14,'location','west',...
        'box','off');
    hl.Position = hl.Position + [-0.02 0.0 0 0];
    hl.ItemTokenSize = [18 18];
    axis off
    
    subplot(3,4,8);
    hold on;
    for j = 1:nr
        h(j) = plot(-1,-1,mk{j},'color','k','markersize',8,'linewidth',1);
        lgtxt{j} = sprintf('$d=%.2g$ $\\mu$m',2*ru(j));
    end
    h(nr+1) = fill([-1 0 0 -1],[1 1 2e2 2e2],0.85*[1 1 1],'edgecolor',0.5*[1 1 1]);
    lgtxt{nr+1} = '$\,\, t_r<t_D$';
    for k = 1:nf
        h(k+nr+1) = plot(-1,-1,'-','color',cmap(k,:),'linewidth',1);
        lgtxt{k+nr+1} = sprintf('$f=%.1f$',fu(k));
    end
    xlim([0 1]);
    ylim([0 1]);
    
    hl = legend(h,lgtxt,'interpreter','latex','fontsize',14,'location','west',...
        'box','off','NumColumns',2);
    hl.Position = hl.Position + [-0.02 0.0 0 0];
    hl.ItemTokenSize = [18 18];
    axis off
    
    subplot(3,4,12);
    hold on;
    for j = 1:nr
        h(j) = plot(-1,-1,mk{j},'color','k','markersize',8,'linewidth',1);
        lgtxt{j} = sprintf('$d=%.2g$ $\\mu$m',2*ru(j));
    end
    h(nr+1) = fill([-1 0 0 -1],[1 1 2e2 2e2],0.85*[1 1 1],'edgecolor',0.5*[1 1 1]);
    lgtxt{nr+1} = '$\,\, t_r<t_D$';
    for k = 1:nf
        h(k+nr+1) = plot(-1,-1,'-','color',cmap(k,:),'linewidth',1);
        lgtxt{k+nr+1} = sprintf('$f=%.1f$',fu(k));
    end
    xlim([0 1]);
    ylim([0 1]);
    
    hl = legend(h,lgtxt,'interpreter','latex','fontsize',14,'location','west',...
        'box','off','NumColumns',2);
    hl.Position = hl.Position + [-0.02 0.0 0 0];
    hl.ItemTokenSize = [18 18];
    axis off
    
    
end


