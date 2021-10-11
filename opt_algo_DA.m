%%
% opt_algo_DA.m
% Construct an optimal estimation algorithm in A(D)
% given the acquisition process and the approximability model
%
% Implements the optimal algorithm described in the paper
% "Approximability models and optimal system identification"
% by M. Ettehad and S. Foucart
% Note: CVX is needed to perform the L1-minimization
% 
% Usage: [a_star,mu] = opt_algo_DA(zeta,V,zeta0)
%
% zeta: vector containing points where the target function is evaluated
%       these points must belong to the torus (i.e., unit circle)
% V: cell containing a basis for the space V
% zeta0: the point where one wishes to estimate the target function
%        this point must also belong to the torus 
%
% a_star: a row vector describing the optimal algorithm
% mu: numerical value of the indicator mu(L_zeta,V,e_zeta0)

% Written by Mahmood Ettehad and Simon Foucart in October 2018
% Send comments to simon.foucart@centraliens.net


function [a_star,mu] = opt_algo_DA(zeta,V,zeta0)

m = length(zeta);
n = length(V);
% the cross-Gramian matrix G and the RHS
G = zeros(m,n);
b = zeros(n,1);
for j=1:n
    Vj = V{j};
    G(:,j) = Vj(zeta);
    b(j) = Vj(zeta0);
end

cvx_quiet true;
% formulate the optimization program
cvx_begin
variable a(m) complex
minimize sum(abs(a))
subject to
transpose(G)*a == b;
cvx_end

% return the outputs
mu = 1 + cvx_optval;
a_star = transpose(a);

end