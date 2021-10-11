% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Reproducible file for the paper
% APPROXIMABILITY MODELS AND OPTIMAL SYSTEM IDENTIFICATION
% by M. Ettehad and S. Foucart
% Written by M. Ettehad and S. Foucart in October 2018
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %


%% Experiment 1: Optimal Identification in the Hardy Space H_2

clear all; clc;
load('seeds.mat');
rng(Seed_Exp1);   % comment out to produce a different figure
% Define random and equispaced points on an inner circle
m = 12;
r = 0.95;
zeta_rand = r*exp(1i*2*pi*rand(m,1));
zeta_equi = r*exp(1i*2*pi/m*(0:m-1)');
% Define the function to be identified, its H_2-norm, 
% its Taylor coefficients, and its point evaluations
F = @(z) 2*z./(2-z.^2);
norm_H2_F = 2/sqrt(3);
V_coeffs_F = zeros(m,1);
V_coeffs_F(2:2:m) = 2*1./2.^(1:m/2); 
y_rand = F(zeta_rand);
y_equi = F(zeta_equi);

% For each type of points, obtain the indicator mu(L_zeta,P_n) 
% when the dimension n of the polynomial space P_n varies 
mu_rand = zeros(1,m);
mu_equi = zeros(1,m);
for n = 1:m
    % define an orthonornal basis for the polynomial space
    clear V;
    V = cell(n,1);
    for j=1:n
        V{j} = @(z) z.^(j-1); 
    end
    % compute the values of the indicator mu(L_zeta,P_n)
    mu_rand(n) = indicator_H2(zeta_rand,V);
    mu_equi(n) = indicator_H2(zeta_equi,V); 
end

% For each type of points, obtain the function recovered from F(zeta),
% its H_2-norm, its coefficients on V, and its coefficients on 
% the representers of the point evaluation functionals
% when the dimension n of the polynomial space P_n varies
F_star_rand   = cell(m,1);   F_star_equi   = cell(m,1);
norm_H2_rand  = cell(m,1);   norm_H2_equi  = cell(m,1);
V_coeffs_rand = cell(m,1);   V_coeffs_equi = cell(m,1);
L_coeffs_rand = cell(m,1);   L_coeffs_equi = cell(m,1);
for n=1:m
    % define an orthonornal basis for the polynomial space
    clear V;
    V = cell(n,1);
    for j=1:n
        V{j} = @(z) z.^(j-1); 
    end
    % compute the desired quantitites
    [F_star_rand{n},norm_H2_rand{n},V_coeffs_rand{n},L_coeffs_rand{n}] ...
        = opt_algo_H2(y_rand,zeta_rand,V);
    [F_star_equi{n},norm_H2_equi{n},V_coeffs_equi{n},L_coeffs_equi{n}] ...
        = opt_algo_H2(y_equi,zeta_equi,V);
end

% Visualization of the results
% the indicator mu
figure(1)  
plot(1:m,mu_rand,'-rx',1:m,mu_equi,'-bo','LineWidth',2)
xlabel('dimension n','FontSize',20)
ylabel('indicator \mu','FontSize',20)
legend({'Random points','Equispaced points'},'FontSize',16,...
    'Location','northwest')
title(strcat('m=',num2str(m),', r=',num2str(r)),'FontSize',20)
% the identification error
figure(2)  
Err_rand = zeros(1,m);
Err_equi = zeros(1,m);
for n=1:m
    Err_rand(n) = norm_H2_F^2 + norm_H2_rand{n}^2 - 2*...
        real( V_coeffs_F(1:n)'*V_coeffs_rand{n} + L_coeffs_rand{n}'*y_rand );
    Err_equi(n) = norm_H2_F^2 + norm_H2_equi{n}^2 - 2*...
        real( V_coeffs_F(1:n)'*V_coeffs_equi{n} + L_coeffs_equi{n}'*y_equi );
end
plot(1:m,Err_rand,'-rx',1:m,Err_equi,'-bo','LineWidth',2)
xlabel('dimension n','FontSize',20)
ylabel('identification error','FontSize',20)
legend({'Random points','Equispaced points'},'FontSize',16,...
    'Location','northeast')
title(strcat('m=',num2str(m),', r=',num2str(r),', F(z) = 2z/(2-z^2)'),'FontSize',20)


%% Experiment 2: Optimal Estimation in the Disc Algebra

clear all; clc;
load('seeds.mat');
rng(Seed_Exp2);   % comment out to produce a different figure
% Define random and equispaced points on the torus
m = 12;
zeta_rand = exp(1i*2*pi*rand(m,1));
zeta_equi = exp(1i*2*pi/m*(0:m-1)');
% Define the target function and its point evaluations
F = @(z) 2*z./(2-z.^2);
y_rand = F(zeta_rand);
y_equi = F(zeta_equi);
% Define the quantity of interest --- point evaluation at zeta0
zeta0 = exp(1i*2*pi*rand);
% select solver used by CVX: gurobi is recommended, default is sdpt3
cvx_solver gurobi;

% For each type of points, obtain the indicator mu(L_zeta,P_n,e_zeta0)
% as well as an optimal estimation of F(zeta0) based on F(zeta)
mu_rand = zeros(1,m);
mu_equi = zeros(1,m);
Err_rand = zeros(1,m);
Err_equi = zeros(1,m);
for n = 1:m
    % define a basis for the polynomial space
    clear V;
    V = cell(n,1);
    for j=1:n
        V{j} = @(z) z.^(j-1); 
    end
    % compute the desired quantities
    [a_star,mu] = opt_algo_DA(zeta_rand,V,zeta0);
    mu_rand(n) = mu;
    Err_rand(n) = abs( F(zeta0)-a_star*y_rand ); 
    [a_star,mu] = opt_algo_DA(zeta_equi,V,zeta0);
    mu_equi(n) = mu; 
    Err_equi(n) = abs( F(zeta0)-a_star*y_equi );
end

% visualization of the results
% the indicator mu
figure(3)  
plot(2:m,mu_rand(2:m),'-rx',2:m,mu_equi(2:m),'-bo', 'LineWidth', 2)
xlabel('dimension n','FontSize',20)
ylabel('indicator \mu','FontSize',20)
legend({'Random points','Equispaced points'},'FontSize',16,...
    'Location','northwest')
title(strcat('m=',num2str(m)),'FontSize',20)
% the estimation error
figure(4)  
plot(2:m,Err_rand(2:m),'-rx',2:m,Err_equi(2:m),'-bo', 'LineWidth', 2)
xlabel('dimension n','FontSize',20)
ylabel('Estimation error','FontSize',20)
legend({'Random points','Equispaced points'},'FontSize',16,...
    'Location','northeast')
title(strcat('m=',num2str(m),', F(z) = 2z/(2-z^2)'),'FontSize',20)

%% 
