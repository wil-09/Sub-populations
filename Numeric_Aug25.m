function Numeric_Aug25
% Numerical integration of the two populations chemotaxis model with
% cortical flows: application of the Euler method. 11th April 2025
tic
clear all; close all; clc
set(0, 'defaultaxesfontsize', 20, 'defaultaxesfontWeight', 'bold', 'defaultaxesLineWidth', 1)
i = sqrt(-1); 

% Gridspace
xlim = 10; x0 = -xlim; xf = xlim; dx = 0.1; t0 = 0; tf = 1; 
dt = 1e-4; x = x0:dx:xf; t = t0:dt:tf; N = numel(x); M = numel(t);

%%%%%%% parameters fo Fig 3.
% r1 = 0.1; r2 = 0.9; alpha1 = 2; alpha2 = 5; mu = 0.035; beta = 100;
% r1 = 0.1; r2 = 0.1; alpha1 = 0.2; alpha2 = 0.5; mu = 0.035; beta = 100;

D1 = 1.76; D2 = 3.29; D3 = 8; chi1 = 0.0649; chi2 = 0.288; u0 = 0.5; 
r1 = 0.1; r2 = 0.9; alpha1 = 0.5; alpha2 = 0.5; mu = 0.035;
h1 = 1; h2 = 1; beta1 = 8; beta2 = 10; k = pi; epsi = 0.01; 

for beta = -5
    sigma1 = chi1*(1 + beta); sigma2 = chi2*(1 - beta);
    n11 = r1 - 2.*D1./dx.^2; n22 = D1./dx.^2 - 0.5.*u0./dx; n66 = - sigma2./dx.^2;
    n33 = D1./dx.^2 + 0.5.*u0./dx; n44 = -chi1./dx.^2; n55 = - sigma1./dx.^2;
    m111 = r2 - 2.*D2./dx.^2; m221 = D2./dx.^2 - 0.5.*u0./dx; 
    m331 = D2./dx.^2 + 0.5.*u0./dx; m441 = -chi2./dx.^2;
    c11 = - mu - 2.*D3./dx.^2; c22 = D3./dx.^2 - 0.5.*u0./dx;
    c33 = D3./dx.^2 + 0.5.*u0./dx;
    
%     Parameters of the steady state
    n0 = (alpha1 - 1)./(alpha1.*alpha2 - 1); m0 = (alpha2 - 1)./(alpha1.*alpha2 - 1);
    c0 = (h1.*(alpha1 - 1) + h2.*(alpha2 - 1))./(mu.*(alpha1.*alpha2 - 1) + beta1.*(alpha1 - 1) + beta2.*(alpha2 - 1));
    eps = [n0, m0, c0]; eps = min(eps)*epsi; x00 = 0; x00n = -1; x00m = 1; 
    nn = zeros(1,N); mm = zeros(1,N); cc = zeros(1,N); Er0 = zeros(1,N);
for l = 1:M
    for j = 1:N
        nn(l,j) = n0 + 1.*eps.*sech(x(j) - x00n).*cos(k*(x(j)));
        mm(l,j) = m0 + 1.*eps.*sech(x(j) - x00m).*cos(k*(x(j)));
        cc(l,j) = c0 + 1.*eps.*sech(x(j) - x00).*cos(k*(x(j)));
        Er0(l,j) = 0*nn(l,j);
    end
end
% First row of of the above defined matries as initial conditions
n(1:N) = nn(1,1:N); c(1:N) = cc(1,1:N); m(1:N) = mm(1,1:N);
Nf1 = []; Cf1 = []; Mf1 = []; Er1 = [];
% the central loop
for l = 1:M
    for j=1:N
        if(j==1)
          jp = j + 1;
          jm = N;
        elseif(j==N)
          jp = 1;
          jm = j-1;
        else
          jp = j + 1;
          jm = j - 1;
        end
        n(j) = n(j) + dt.*(n11.*n(j) + n22.*n(jp) + n33.*n(jm) - r1.*n(j).*(n(j) + alpha1.*m(j)) + ...
            0.5.*n44.*((n(jp) + n(j)).*(c(jp) - c(j)) - (n(j) + n(jm)).*(c(j) - c(jm))) + ...
            0.5.*n55.*((n(jp) + n(j)).*(n(jp) - n(j)) - (n(j) + n(jm)).*(n(j) - n(jm))) + ...
            0.5.*n66.*((n(jp) + n(j)).*(m(jp) - m(j)) - (n(j) + n(jm)).*(m(j) - m(jm))));
        
        m(j) = m(j) + dt.*(m111.*m(j) + m221.*m(jp) + m331.*m(jm) - r2.*m(j).*(m(j) + alpha2.*n(j)) + ...
            0.5.*m441.*((m(jp) + m(j)).*(c(jp) - c(j)) - (m(j) + m(jm)).*(c(j) - c(jm))) + ...
            0.5.*n55.*((m(jp) + m(j)).*(n(jp) - n(j)) - (m(j) + m(jm)).*(n(j) - n(jm))) + ...
            0.5.*n66.*((m(jp) + m(j)).*(m(jp) - m(j)) - (m(j) + m(jm)).*(m(j) - m(jm))));

        c(j) = c(j) + dt.*(c11.*c(j) + c22.*c(jp) + c33.*c(jm) + n(j).*(h1 - beta1.*c(j)) + m(j).*(h2 - beta2.*c(j)) + ...
            0.5.*n55.*((c(jp) + c(j)).*(n(jp) - n(j)) - (c(j) + c(jm)).*(n(j) - n(jm))) + ...
            0.5.*n66.*((c(jp) + c(j)).*(m(jp) - m(j)) - (c(j) + c(jm)).*(m(j) - m(jm))));

        Nf1(l,j) = n(j); Mf1(l,j) = m(j); Cf1(l,j) = c(j);
        Er1(l,j) = (Nf1(l,j) - nn(l,j));
    end
end
plot(x, Mf1(round(0.75*M),:), 'linewidth', 3); hold on
end

for beta = -15
    sigma1 = chi1*(1 + beta); sigma2 = chi2*(1 - beta);
    n11 = r1 - 2.*D1./dx.^2; n22 = D1./dx.^2 - 0.5.*u0./dx; n66 = - sigma2./dx.^2;
    n33 = D1./dx.^2 + 0.5.*u0./dx; n44 = -chi1./dx.^2; n55 = - sigma1./dx.^2;
    m111 = r2 - 2.*D2./dx.^2; m221 = D2./dx.^2 - 0.5.*u0./dx; 
    m331 = D2./dx.^2 + 0.5.*u0./dx; m441 = -chi2./dx.^2;
    c11 = - mu - 2.*D3./dx.^2; c22 = D3./dx.^2 - 0.5.*u0./dx;
    c33 = D3./dx.^2 + 0.5.*u0./dx;
    
%     Parameters of the steady state
    n0 = (alpha1 - 1)./(alpha1.*alpha2 - 1); m0 = (alpha2 - 1)./(alpha1.*alpha2 - 1);
    c0 = (h1.*(alpha1 - 1) + h2.*(alpha2 - 1))./(mu.*(alpha1.*alpha2 - 1) + beta1.*(alpha1 - 1) + beta2.*(alpha2 - 1));
    eps = [n0, m0, c0]; eps = min(eps)*epsi; x00 = 0; x00n = -1; x00m = 1; 
    nn = zeros(1,N); mm = zeros(1,N); cc = zeros(1,N); Er0 = zeros(1,N);
for l = 1:M
    for j = 1:N
        nn(l,j) = n0 + 1.*eps.*sech(x(j) - x00n).*cos(k*(x(j)));
        mm(l,j) = m0 + 1.*eps.*sech(x(j) - x00m).*cos(k*(x(j)));
        cc(l,j) = c0 + 1.*eps.*sech(x(j) - x00).*cos(k*(x(j)));
        Er0(l,j) = 0*nn(l,j);
    end
end
% First row of of the above defined matries as initial conditions
n(1:N) = nn(1,1:N); c(1:N) = cc(1,1:N); m(1:N) = mm(1,1:N);
Nf1 = []; Cf1 = []; Mf1 = []; Er1 = [];
% the central loop
for l = 1:M
    for j=1:N
        if(j==1)
          jp = j + 1;
          jm = N;
        elseif(j==N)
          jp = 1;
          jm = j-1;
        else
          jp = j + 1;
          jm = j - 1;
        end
        n(j) = n(j) + dt.*(n11.*n(j) + n22.*n(jp) + n33.*n(jm) - r1.*n(j).*(n(j) + alpha1.*m(j)) + ...
            0.5.*n44.*((n(jp) + n(j)).*(c(jp) - c(j)) - (n(j) + n(jm)).*(c(j) - c(jm))) + ...
            0.5.*n55.*((n(jp) + n(j)).*(n(jp) - n(j)) - (n(j) + n(jm)).*(n(j) - n(jm))) + ...
            0.5.*n66.*((n(jp) + n(j)).*(m(jp) - m(j)) - (n(j) + n(jm)).*(m(j) - m(jm))));
        
        m(j) = m(j) + dt.*(m111.*m(j) + m221.*m(jp) + m331.*m(jm) - r2.*m(j).*(m(j) + alpha2.*n(j)) + ...
            0.5.*m441.*((m(jp) + m(j)).*(c(jp) - c(j)) - (m(j) + m(jm)).*(c(j) - c(jm))) + ...
            0.5.*n55.*((m(jp) + m(j)).*(n(jp) - n(j)) - (m(j) + m(jm)).*(n(j) - n(jm))) + ...
            0.5.*n66.*((m(jp) + m(j)).*(m(jp) - m(j)) - (m(j) + m(jm)).*(m(j) - m(jm))));

        c(j) = c(j) + dt.*(c11.*c(j) + c22.*c(jp) + c33.*c(jm) + n(j).*(h1 - beta1.*c(j)) + m(j).*(h2 - beta2.*c(j)) + ...
            0.5.*n55.*((c(jp) + c(j)).*(n(jp) - n(j)) - (c(j) + c(jm)).*(n(j) - n(jm))) + ...
            0.5.*n66.*((c(jp) + c(j)).*(m(jp) - m(j)) - (c(j) + c(jm)).*(m(j) - m(jm))));

        Nf1(l,j) = n(j); Mf1(l,j) = m(j); Cf1(l,j) = c(j);
        Er1(l,j) = (Nf1(l,j) - nn(l,j));
    end
end
plot(x, Mf1(round(0.75*M),:), 'linewidth', 3); hold on
end

for beta = -25
    sigma1 = chi1*(1 + beta); sigma2 = chi2*(1 - beta);
    n11 = r1 - 2.*D1./dx.^2; n22 = D1./dx.^2 - 0.5.*u0./dx; n66 = - sigma2./dx.^2;
    n33 = D1./dx.^2 + 0.5.*u0./dx; n44 = -chi1./dx.^2; n55 = - sigma1./dx.^2;
    m111 = r2 - 2.*D2./dx.^2; m221 = D2./dx.^2 - 0.5.*u0./dx; 
    m331 = D2./dx.^2 + 0.5.*u0./dx; m441 = -chi2./dx.^2;
    c11 = - mu - 2.*D3./dx.^2; c22 = D3./dx.^2 - 0.5.*u0./dx;
    c33 = D3./dx.^2 + 0.5.*u0./dx;
    
%     Parameters of the steady state
    n0 = (alpha1 - 1)./(alpha1.*alpha2 - 1); m0 = (alpha2 - 1)./(alpha1.*alpha2 - 1);
    c0 = (h1.*(alpha1 - 1) + h2.*(alpha2 - 1))./(mu.*(alpha1.*alpha2 - 1) + beta1.*(alpha1 - 1) + beta2.*(alpha2 - 1));
    eps = [n0, m0, c0]; eps = min(eps)*epsi; x00 = 0; x00n = -1; x00m = 1; 
    nn = zeros(1,N); mm = zeros(1,N); cc = zeros(1,N); Er0 = zeros(1,N);
for l = 1:M
    for j = 1:N
        nn(l,j) = n0 + 1.*eps.*sech(x(j) - x00n).*cos(k*(x(j)));
        mm(l,j) = m0 + 1.*eps.*sech(x(j) - x00m).*cos(k*(x(j)));
        cc(l,j) = c0 + 1.*eps.*sech(x(j) - x00).*cos(k*(x(j)));
        Er0(l,j) = 0*nn(l,j);
    end
end
% First row of of the above defined matries as initial conditions
n(1:N) = nn(1,1:N); c(1:N) = cc(1,1:N); m(1:N) = mm(1,1:N);
Nf1 = []; Cf1 = []; Mf1 = []; Er1 = [];
% the central loop
for l = 1:M
    for j=1:N
        if(j==1)
          jp = j + 1;
          jm = N;
        elseif(j==N)
          jp = 1;
          jm = j-1;
        else
          jp = j + 1;
          jm = j - 1;
        end
        n(j) = n(j) + dt.*(n11.*n(j) + n22.*n(jp) + n33.*n(jm) - r1.*n(j).*(n(j) + alpha1.*m(j)) + ...
            0.5.*n44.*((n(jp) + n(j)).*(c(jp) - c(j)) - (n(j) + n(jm)).*(c(j) - c(jm))) + ...
            0.5.*n55.*((n(jp) + n(j)).*(n(jp) - n(j)) - (n(j) + n(jm)).*(n(j) - n(jm))) + ...
            0.5.*n66.*((n(jp) + n(j)).*(m(jp) - m(j)) - (n(j) + n(jm)).*(m(j) - m(jm))));
        
        m(j) = m(j) + dt.*(m111.*m(j) + m221.*m(jp) + m331.*m(jm) - r2.*m(j).*(m(j) + alpha2.*n(j)) + ...
            0.5.*m441.*((m(jp) + m(j)).*(c(jp) - c(j)) - (m(j) + m(jm)).*(c(j) - c(jm))) + ...
            0.5.*n55.*((m(jp) + m(j)).*(n(jp) - n(j)) - (m(j) + m(jm)).*(n(j) - n(jm))) + ...
            0.5.*n66.*((m(jp) + m(j)).*(m(jp) - m(j)) - (m(j) + m(jm)).*(m(j) - m(jm))));

        c(j) = c(j) + dt.*(c11.*c(j) + c22.*c(jp) + c33.*c(jm) + n(j).*(h1 - beta1.*c(j)) + m(j).*(h2 - beta2.*c(j)) + ...
            0.5.*n55.*((c(jp) + c(j)).*(n(jp) - n(j)) - (c(j) + c(jm)).*(n(j) - n(jm))) + ...
            0.5.*n66.*((c(jp) + c(j)).*(m(jp) - m(j)) - (c(j) + c(jm)).*(m(j) - m(jm))));

        Nf1(l,j) = n(j); Mf1(l,j) = m(j); Cf1(l,j) = c(j);
        Er1(l,j) = (Nf1(l,j) - nn(l,j));
    end
end
plot(x, Mf1(round(0.75*M),:), 'linewidth', 3)
end
xlabel 'x', ylabel 'm(r,t)'
legend('\beta = -5', '\beta = -15', '\beta = -25')
toc

