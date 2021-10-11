%%
% indicator_H2.m
% Computes the indicator of compatibility in H_2 
% of the acquisition process and the approximability model
%
% Implements the eigendecomposition formula presented in the paper
% "Approximability models and optimal system identification"
% by M. Ettehad and S. Foucart
% i.e., computes 1/lambda_min(G'*inv(H)*G)^(1/2)
% with cross-Gramian and Gramian matrices G and H defined in the paper
% 
% Usage: mu = indicator_H2(zeta,V)
%
% zeta: vector containing points where the target function is evaluated
% V: cell containing an orthonormal basis for the space V
%
% mu: numerical value of the indicator mu(L_zeta,V)

% Written by Mahmood Ettehad and Simon Foucart in October 2018
% Send comments to simon.foucart@centraliens.net

function mu = indicator_H2(zeta,V)

m = length(zeta);
n = length(V);
% the Riesz representers of the point evaluations
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
% the Gramian matrix H
H = zeros(m,m);
for j=1:m
    H(:,j) = L{j}(zeta);
end
% the value of mu
mu = 1/sqrt(real(eigs(G'*inv(H)*G,1,'SM')));

end