function [y, e, psi] = adapt_nlms(N ,M, d_i, u_i, mu, psi_prev )
%ADAPT_NLMS Performs the adaptation of N nlms filters.
%
%   It used Normalized least-mean-squares (NLMS) algorithm to update in
%   N estimations for using in the adaptation step of Diffusion networks.
%
%
% INPUT N: number of nodes in the network.
%       M: number of coefficients of w to be estimated.
%       d_i: Vector outputs (N, one per node) a iteration i. 
%       u_i: Matrix of inputs (NxM) a iteration i.
%       mu: Unormalized step size for the NMLS filter
%       psi_prev: matrix with previous combined estimations of all the nodes.
%
%
%
% OUTPUT y: vector with the predicted output of each node.
%        e: vector with the a priori errors of each node.
%        psi: vector with the local estimations of each node.
%
%
%
% created with MATLAB ver.: 8.0.0.783 (R2012b) on Mac OS X  Version: 10.9.4
%
% created by: Jesus Fernandez-Bes (<a href="http://www.tsc.uc3m.es/~jesusfbes">web</a>)
% DATE: Oct-2014


psi = zeros(M, N);
y = zeros(N, 1);
e = zeros(N, 1);

for k = 1:N % for each node k
    
    
    % output of each local filter
    y(k) = u_i(k, :) * psi_prev(:, k);
    
    % Error a priori
    e(k) = d_i(k) - y(k);
    
    % Compute normalized step size
    mu_norm = (mu(k) ./ sum(u_i(k, :) .* u_i(k, :)));
    
    % Compute each estimation
    psi(:, k) = psi_prev(:, k) + mu_norm * (u_i(k, :)' *  e(k));
    
end


end