%%
% opt_algo_H2.m
% Computes the output of an optimal recovery algorithm in H_2
% given the acquisition process and the approximability model
%
% Implements the optimal algorithm described in the paper
% "Approximability models and optimal system identification"
% by M. Ettehad and S. Foucart
% 
% Usage: [F_star,norm_H2,V_coeffs,L_coeffs] = opt_algo_H2(y,zeta,V)
%
% y: vector containing values of the target function 
% zeta: vector containing points where the target function is evaluated
% V: cell containing an orthonormal basis for the space V
%
% F_star: function handle representing the recovered function
% norm_H2: the H2-norm of the recovered function
% V_coeffs: the coefficients of the recovered function on the basis for V
% L_coeffs: the coefficients of the recovered function on the system made 
% of the Riesz representers of point evaluations at the zeta's

% Written by Mahmood Ettehad and Simon Foucart in October 2018
% Send comments to simon.foucart@centraliens.net

function [F_star,norm_H2,V_coeffs,L_coeffs] = opt_algo_H2(y,zeta,V)

m = length(zeta);
n = length(V);
% the Riesz representers of point evaluations at the zeta's
L = cell(m);
for k=1:m
   L{k} = @(z) 1./(1-zeta(k)'*z); 
end
% the cross-Gramian matrix G
G = zeros(m,n);
for j=1:n
    Vj = V{j};
    G(:,j) = Vj(zeta);
end
% the Gramian matrix H and its inverse
H = zeros(m,m);
for j=1:m
    H(:,j) = L{j}(zeta);
end
Hinv = inv(H);
% the coefficients on V and L
V_coeffs = (G'*Hinv*G)\(G'*Hinv*y);
L_coeffs = H\(y-G*V_coeffs);
% the recovered function
F_star = @(z) 0;
for j=1:n
   F_star = @(z) F_star(z) + V_coeffs(j)*V{j}(z); 
end
for k=1:m
   F_star = @(z) F_star(z) + L_coeffs(k)*L{k}(z);  
end
% the H_2 norm of the recovered function
norm_H2 = sqrt( real(V_coeffs'*V_coeffs + L_coeffs'*H*L_coeffs ...
    + 2*L_coeffs'*G*V_coeffs) ); 

end