function [msd_comb, e_comb, c_out, w_out]= ...
    atc_nlms_acw(d, u, w0, mu, A, N, Tmax, nu, error_params, is_returned_w)
%ATC_NLMS_ACW Simulates a Diffusion network with ATC and Adaptive Combiners
%
%
% It computes the adaptive combination weights according to:
%   [1]	S.-Y. Tu, A.H. Sayed
%	Optimal combination rules for adaptation and learning over networks
%   Proc. of 4th IEEE Intl. Workshop on Computational Advances in Multi-Sensor
%   Adaptive Process, CAMSAP, San Juan, Puerto Rico (2011), pp. 317?320
%
%
% For details in the algorithm, please see [1].
%
% INPUT d: desired ouput signals. Matrix N x Tmax (one row per node).
%       u: input signals. Matrix N x Tmax (one row per node).
%       w0: vector to estimate. Array MxNxTmax (vector Mx1 per node per iteration)
%       mu: unnormalized stepsize for the NLMS filter
%       A: adjacency matrix of the network.
%       N: number of nodes in the network.
%       Tmax: maximum number of iterations.
%       nu: learning parameter for the combination.
%       error_params: Parameters of disconnection model.
%       is_returned_w: Flag, if 1 the full w estimated matrix is returned.
%
%
% OUTPUT msd_comb:
%        e_comb:
%        c_out:
%        w_out:
%
%
% See Also UPDATE_COMBINE_ACW, ADAPT_NLMS
%
% created with MATLAB ver.: 8.0.0.783 (R2012b) on Mac OS X  Version: 10.9.4
%
% created by: Jesus Fernandez-Bes (<a href="http://www.tsc.uc3m.es/~jesusfbes">web</a>)
% DATE: Oct-2014

%% Parameters



% number of coefficients
M = size(w0, 1);


% The dead nodes
dead_nodes = zeros(N, 1);


% Initialize filter coefficients

w = randn(M, N);

w_prev = w;


% Initialize combination coefficientes and intermediate variables

Nk = sum(A); % Number of neighbours per node

b = cell(N, 1);
ones_Nk = cell(N, 1);

gamma2 = zeros(N, N);

for k = 1:N
    
    ones_Nk{k} = ones(Nk(k), 1);
    b{k} = ones(Nk(k), 1) / Nk(k);
    
end

% Preallocate results
msd_comb = zeros(N,Tmax);


c_out = zeros(N, N, Tmax);

e_comb = zeros(N,Tmax);


if is_returned_w == 1
    w_out = zeros(M,N,Tmax);
else
    w_out = [];
end




for i = M:Tmax-1
    
    dead_nodes = get_dead_nodes(dead_nodes, N, i, error_params);
    
    w_prev(:, dead_nodes == 0) = w(:, dead_nodes == 0);
    
    %% 1- FILTER and ADAPT
    u_k = u(:, i:-1:i-M+1) ;
    
    [dummy, e_k, psi] = adapt_nlms(N , M , d(:, i), u_k, mu, w_prev ); %#ok<*ASGLU>
    
    %% 2 .- COMBINATION WEIGHTS UPDATE
    
    [w, c, gamma2] = update_combine_acw(N,M, A, psi, w_prev, gamma2, nu);
    
    % save results
    msd_comb(:, i) = sum((w0(:, :, i) - w_prev).^2);
    
    c_out(:, :, i) = c;
    
    e_comb(:, i) = e_k;
    
    
    if is_returned_w == 1
        w_out(:, :, i) = w;
    end
    
end
