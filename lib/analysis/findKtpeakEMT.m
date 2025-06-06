function [tpeak, tinit] = findKtpeakEMT(f,kappa,R,Din,Dex)
p1  = f;
p2  = 1-f;
r   = kappa*(3/R);
tex = p2/r;
rex = 1/tex;
tc = R^2/Dex;
% a = @(t) sphereD(R,Din,t).*t*2;
% % a = @(t) R.^2/5*2;
tinit = sqrt(6/5*tex/p2*tc);

a = @(t) sphereD(R,Din,t).*t*2;
b = 2*Dex;
K = @(t) ( 2*p1*p2*(a(t).^2 - 2*a(t).*(a(t)*rex + b).*t + 2/rex^2*(a(t)*rex + b).^2 .* (exp(-rex*t) - 1 + rex*t)) ) ./ ( p1*a(t) + p2*b*t ).^2;
ti = 0.5*tinit : 0.01 : 4*tinit;
[~, I] = max(K(ti));
tpeak = ti(I);

% syms t
% a = sphereD(R,Din,t).*t*2;
% b = 2*Dex;
% K = ( 2*p1*p2*(a.^2 - 2*a*(a*rex + b)*t + 2/rex^2*(a*rex + b)^2 * (exp(-rex*t) - 1 + rex*t)) ) / ( p1*a + p2*b*t )^2;
% fun = diff(K,t)==0;
% tpeak = vpasolve(fun, t, 2*tinit);
% tpeak = double(tpeak);

end

function D = sphereD(r,D0,t)
%SPHEREDK    diffusivity and kurtosis inside a sphere
%   [D,K] = SPHEREDK(r,D0,t,N) returns diffusivity D and kurtosis K inside
%   a sphere in radius r at diffusion time t. The intrinsic diffusivity
%   inside the sphere is D0, and the number of summation terms is N.
%
%   Author: Hong-Hsi Lee, September, 2018 (orcid.org/0000-0002-3663-6559)

% sbesselj = @(n,x) sqrt(pi/2./(x+eps)).*besselj(n+1/2,x);

% dj0 = @(x) -sbesselj(1,x);
% alpha0k = @(k) fzero(dj0, [(k-1)*pi, k*pi]); 
% alpha0 = zeros(1,N); for k=1:N, alpha0(k) = alpha0k(k+1); end
% % alpha0 = [4.4934    7.7253   10.9041   14.0662   17.2208   20.3713   23.5195   26.6661   29.8116   32.9564];
% alpha0 = [NaN alpha0(1:N-1)];

% dj1 = @(x) sbesselj(1,x)./(x + eps) - sbesselj(2,x);
% alpha1k = @(k) fzero(dj1, [(k-1)*pi+1e-10, k*pi-1e-10]); 
% alpha1 = zeros(1,N); for k=1:N, alpha1(k) = alpha1k(k); end
alpha1 = [2.0816    5.9404    9.2058   12.4044   15.5792   18.7426   21.8997   25.0528   28.2034   31.3521];

% dj2 = @(x) sbesselj(2,x)*2./(x + eps) - sbesselj(3,x);
% alpha2k = @(k) fzero(dj2,[floor((k-1)*pi), k*pi]);
% alpha2 = zeros(1,N); for k=1:N, alpha2(k) = alpha2k(k+1); end
% % alpha2 = [3.3421    7.2899   10.6139   13.8461   17.0429   20.2219   23.3905   26.5526   29.7103   32.8649];
% alpha2 = [NaN alpha2(1:N-1)];

D = r^2/5./t;
tc = r^2/D0;
for i = 1:numel(alpha1)
    D = D - 2*r^2./t.*exp(-alpha1(i)^2*t/tc)/alpha1(i)^2/(alpha1(i)^2-2);
end

% numer = 3/25/7 + exp(-alpha1(1)^2*t/tc)/alpha1(1)^2/(alpha1(1)^2-2)*(4/alpha1(1)^2-6/5);
% denom = D./(r^2./t);
% for i = 2:numel(alpha0)
%     numer = numer + exp(-alpha0(i)^2*t/tc)/alpha0(i)^4*2/3 +...
%         exp(-alpha1(i)^2*t/tc)/alpha1(i)^2/(alpha1(i)^2-2)*(4/alpha1(i)^2-6/5) +...
%         exp(-alpha2(i)^2*t/tc)/alpha2(i)^2/(alpha2(i)^2-6)*8/15;
% end
% numer = 6*numer;
% denom = (denom).^2;
% K = numer./denom-3;

end