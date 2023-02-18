function TwoChemPop_Dim_Good
% Numerical integration of the two populations chemotaxis model

clear all; close all; clc
set(0, 'defaultaxesfontsize', 20, 'defaultaxesfontWeight', 'bold', 'defaultaxesLineWidth', 1)
i = sqrt(-1); 
D1 = 1.76; D2 = 3.29; D3 = 8; chi1 = 0.649; chi2 = 0.288; delta = 0.9;
r1 = 0.1; r2 = 0.9; alpha1 = 5; alpha2 = 1; beta1 = 0.1; beta2 = 0.3; 
h1 = 0.2; h2 = 0.5; mu = 5*10^4; 

% steepness of the analytical solution
lambda = 0.9;

% Gridspace
xlim = 20; x0 = -xlim; xf = xlim; dx = 0.01; t0 = 0; tf = 0.05; 
dt = dx*(tf - t0)/(xf - x0), x = x0:dx:xf; t = t0:dt:tf; N = numel(x); M = numel(t)

% [x,t] = meshgrid(x,t);
nn = zeros(1,N); mm = zeros(1,N); cc = zeros(1,N);

n11 = r1 - 2.*D1./dx.^2; n22 = D1./dx.^2 - 0.5.*delta./dx; 
n33 = D1./dx.^2 + 0.5.*delta./dx; n44 =  -chi1./dx.^2;
m11 = r2 - 2.*D2./dx.^2; m22 = D2./dx.^2 - 0.5.*delta./dx; 
m33 = D2./dx.^2 + 0.5.*delta./dx; m44 =  -chi2./dx.^2;
c11 = -mu - 2.*D3./dx.^2; c22 = D3./dx.^2 - 0.5.*delta./dx;
c33 = D3./dx.^2 + 0.5.*delta./dx;

% initial conditions of the model
for k = 1:M
    for j = 1:N
        nn(j,k) = 0.1.*exp(-(x(j) - 0./3).^2./0.5).*sech(x(j) - lambda.*t(k));
%         mm(j,k) = 1 - nn(j,k);
        mm(j,k) = 0.3.*sech(-(x(j) - x0./115).^2./0.9).*exp(x(j) - lambda.*t(k));
        cc(j,k) = 0.0005;%.*cos(0.11.*pi.*x(j)).*tanh(x(j) - lambda.*t(k));
    end
end
n(1:N) = nn(1:N,1); c(1:N) = cc(1:N,1); m(1:N) = mm(1:N,1);
% % % % figure
% % % % subplot(1,3,1)
% % % % plot(x, n, 'linewidth',3)
% % % % %   plot(n, v(:,1), n, v(:,20), n, v(:,50),'linewidth',3)
% % % % xlabel 'x'; ylabel 'n(x,0)'
% % % % 
% % % % subplot(1,3,2)
% % % % plot(x, m, 'linewidth',3)
% % % % %   plot(n, v(:,1), n, v(:,20), n, v(:,50),'linewidth',3)
% % % % xlabel 'x'; ylabel 'm(x,0)'
% % % % 
% % % % subplot(1,3,3)
% % % % plot(x, c, 'linewidth',3)
% % % % %   plot(n, v(:,1), n, v(:,20), n, v(:,50),'linewidth',3)
% % % % xlabel 'x'; ylabel 'c(x,0)'

k1n = zeros(1,N); k2n = zeros(1,N); k3n = zeros(1,N); k4n = zeros(1,N);
k1c = zeros(1,N); k2c = zeros(1,N); k3c = zeros(1,N); k4c = zeros(1,N);
k1m = zeros(1,N); k2m = zeros(1,N); k3m = zeros(1,N); k4m = zeros(1,N);
n1 = zeros(1,N); n2 = zeros(1,N); n3 = zeros(1,N);
c1 = zeros(1,N); c2 = zeros(1,N); c3 = zeros(1,N);
m1 = zeros(1,N); m2 = zeros(1,N); m3 = zeros(1,N);

Nf = []; Cf = []; Mf = [];
% Periodic boundary conditions
for k = 1:M
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
   
%  RK4_M1
%  now the central loop
%-----------------------------------------------------------------------
  k1n(j) = n11.*n(j) + n22.*n(jp) + n33.*n(jm) - r1.*n(j).*(n(j) + alpha1.*m(j)) + ...
      0.5.*n44.*((n(jp) + n(j)).*(c(jp) - c(j)) - (n(j) + n(jm)).*(c(j) - c(jm)));
  k1m(j) = m11.*m(j) + m22.*m(jp) + m33.*m(jm) - r2.*m(j).*(m(j) + alpha2.*n(j)) + ...
      0.5.*m44.*((m(jp) + m(j)).*(c(jp) - c(j)) - (m(j) + m(jm)).*(c(j) - c(jm)));
  k1c(j) = c11.*c(j) + c22.*c(jp) + c33.*c(jm) + n(j).*(h1 - beta1.*c(j)) + m(j).*(h2 - beta2.*c(j));
  
  n1(j) = n(j) + dt.*k1n(j)./2;
  m1(j) = m(j) + dt.*k1m(j)./2;
  c1(j) = c(j) + dt.*k1c(j)./2;
   %---------------------------------------------------------------------
   
   %----------------------------------------------------------------------
   
  k2n(j) = n11.*n1(j) + n22.*n1(jp) + n33.*n1(jm) - r1.*n1(j).*(n1(j) + alpha1.*m1(j)) + ...
      0.5.*n44.*((n1(jp) + n1(j)).*(c1(jp) - c1(j)) - (n1(j) + n1(jm)).*(c1(j) - c1(jm)));
  k2m(j) = m11.*m1(j) + m22.*m1(jp) + m33.*m1(jm) - r2.*m1(j).*(m1(j) + alpha2.*n1(j)) + ...
      0.5.*m44.*((m1(jp) + m1(j)).*(c1(jp) - c1(j)) - (m1(j) + m1(jm)).*(c1(j) - c1(jm)));
  k2c(j) = c11.*c1(j) + c22.*c1(jp) + c33.*c1(jm) + n1(j).*(h1 - beta1.*c1(j)) + m1(j).*(h2 - beta2.*c1(j));
  
  n2(j) = n(j) + dt.*k1n(j)./2;
  m2(j) = m(j) + dt.*k1m(j)./2;
  c2(j) = c(j) + dt.*k1c(j)./2;
        %----------------------------------------------------------------
        
        %----------------------------------------------------------------------
  k3n(j) = n11.*n2(j) + n22.*n2(jp) + n33.*n2(jm) - r1.*n2(j).*(n2(j) + alpha1.*m2(j)) + ...
      0.5.*n44.*((n2(jp) + n2(j)).*(c2(jp) - c2(j)) - (n2(j) + n2(jm)).*(c2(j) - c2(jm)));
  k3m(j) = m11.*m2(j) + m22.*m2(jp) + m33.*m2(jm) - r2.*m2(j).*(m2(j) + alpha2.*n2(j)) + ...
      0.5.*m44.*((m2(jp) + m2(j)).*(c2(jp) - c2(j)) - (m2(j) + m2(jm)).*(c2(j) - c2(jm)));
  k3c(j) = c11.*c2(j) + c22.*c2(jp) + c33.*c2(jm) + n2(j).*(h1 - beta1.*c2(j)) + m2(j).*(h2 - beta2.*c2(j));
  
  n3(j) = n(j) + dt.*k2n(j)./2;
  m3(j) = m(j) + dt.*k2m(j)./2;
  c3(j) = c(j) + dt.*k2c(j)./2;      
        %----------------------------------------------------------------
        
        %----------------------------------------------------------------------
  k4n(j) = n11.*n3(j) + n22.*n3(jp) + n33.*n3(jm) - r1.*n3(j).*(n3(j) + alpha1.*m3(j)) + ...
      0.5.*n44.*((n3(jp) + n3(j)).*(c3(jp) - c3(j)) - (n3(j) + n3(jm)).*(c3(j) - c3(jm)));
  k4m(j) = m11.*m3(j) + m22.*m3(jp) + m33.*m3(jm) - r2.*m3(j).*(m3(j) + alpha2.*n3(j)) + ...
      0.5.*m44.*((m3(jp) + m3(j)).*(c3(jp) - c3(j)) - (m3(j) + m3(jm)).*(c3(j) - c3(jm)));
  k4c(j) = c11.*c3(j) + c22.*c3(jp) + c33.*c3(jm) + n3(j).*(h1 - beta1.*c3(j)) + m3(j).*(h2 - beta2.*c3(j));
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
  n(j) = n(j) + dt.*(k1n(j) + 2.*k2n(j) + 2.*k3n(j) + k4n(j))./6;
  m(j) = m(j) + dt.*(k1m(j) + 2.*k2m(j) + 2.*k3m(j) + k4m(j))./6;
  c(j) = c(j) + dt.*(k1c(j) + 2.*k2c(j) + 2.*k3c(j) + k4c(j))./6;
  
  Nf(j,k) = n(j); Mf(j,k) = m(j); Cf(j,k) = c(j);
    end
end

% % % figure
% % % mesh(x,t,Nf); colormap(jet); 
% % % xlabel('x'); ylabel('t'); zlabel('n(x,t)')
% % % 
% % % figure
% % % mesh(x,t,Mf); colormap(jet); 
% % % xlabel('x'); ylabel('t'); zlabel('m(x,t)')
% % % 
% % % figure
% % % mesh(x,t,Cf); colormap(jet); 
% % % xlabel('x'); ylabel('t'); zlabel('c(x,t)')

figure
subplot(2,2,1)
plot(x, nn(:,1), 'g', x, Nf(:,1), 'r', 'linewidth', 3)
xlabel 'x'; ylabel 'n(x,t)'
legend('analytical', 'numerical')

subplot(2,2,2)
plot(x, nn(:,20), 'g', x, Nf(:,20), 'r', 'linewidth', 3)
xlabel 'x'; ylabel 'n(x,t)'
legend('analytical', 'numerical')

subplot(2,2,3)
plot(x, nn(:,80), 'g', x, Nf(:,80), 'r', 'linewidth', 3)
xlabel 'x'; ylabel 'n(x,t)'
legend('analytical', 'numerical')

subplot(2,2,4)
plot(x, nn(:,99), 'g', x, Nf(:,99), 'r', 'linewidth', 3)
xlabel 'x'; ylabel 'n(x,t)'
legend('analytical', 'numerical')

figure
plot(x, Nf(:,1), 'g', x, Nf(:,100), 'r', x, Nf(:,500), 'b', x, Nf(:,1500), 'c', x, Nf(:,2000), 'y', 'linewidth', 3); hold on
xlabel 'x'; ylabel 'n(x,t)'
legend('t = 0', 't = 100', 't = 500', 't = 1500', 't = 2000')

figure
plot(x, Mf(:,1), 'g', x, Mf(:,100), 'r', x, Mf(:,500), 'b', x, Mf(:,1500), 'c', x, Mf(:,2000), 'y', 'linewidth', 3)
xlabel 'x'; ylabel 'm(x,t)'
legend('t = 0', 't = 100', 't = 500', 't = 1500', 't = 2000')

figure
plot(x, Cf(:,1000), 'r', x, Cf(:,2000), 'b', x, Cf(:,3000), 'c', x, Cf(:,4000), 'y', 'linewidth', 3)
xlabel 'x'; ylabel 'c(x,t)'
legend('t = 100', 't = 500', 't = 1500', 't = 2000')