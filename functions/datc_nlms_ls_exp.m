function [msd_comb, e_comb, c_out, w_out] = ...
    datc_nlms_ls_exp(d, u, w0, mu, A, N, gamma, regul, Tmax, error_params, is_returned_w)
%LE_ADN_NLMS_LS Simulates a Diffusion network with local estimation ATC and Adaptive Combiners
%
%
% It computes the adaptive combination weights according to:
%   [1]  Fernandez-Bes, J., Azpicueta-Ruiz, L. A., Arenas-García, J., & Silva, M. T. (2014). 
%   Distributed estimation in diffusion networks using affine least-squares combiners. 
%   Digital Signal Processing.
% 
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
%       error_params: Parameters of disconnection model.
%       is_returned_w: Flag, if 1 the full w estimated matrix is returned.
%
%
% OUTPUT msd_comb, msd_sep, e_comb, c_out, w_out
%
% See Also UPDATE_COMBINE_LS, ADAPT_NLMS, GET_DEAD_NODES
%
% created with MATLAB ver.: 8.0.0.783 (R2012b) on Mac OS X  Version: 10.9.4
%
% created by: Jesus Fernandez-Bes (<a href="http://www.tsc.uc3m.es/~jesusfbes">web</a>)
% DATE: Oct-2014


%% Parameters

% number of coefficients
M = size(w0, 1);

delta = regul;


dead_nodes = zeros(N, 1);

% Initialize filter coefficients
psi = zeros(M,N);
w = zeros(M,N);


P_sum = cell(N, Tmax);
z_sum = cell(N, Tmax);

%Y_tilde_list = zeros(N, N, L);
%z_list = zeros(N, N, L);

% Preallocate results
msd_comb = zeros(N, Tmax);

c_out = zeros(N, N, Tmax);

e_comb = zeros(N, Tmax);

ecomb_i = zeros(N,1);

Nk = sum(A) - 1;

if is_returned_w == 1
    w_out = zeros(M, N, Tmax);
end

%% Initialize estimation structures
for k = 1:N
    
    % We update based in the estimation of each l-neighbour
    % Neighbours are the nodes connected to k
    
    l_k = (A(:, k) == 1);
    
    l_k(k) = 0;
    
    y_aux = zeros(N);
    y_aux(l_k, l_k ) = sqrt(delta) * eye(Nk(k));
    
    
    P_sum{k, M-1} = y_aux(l_k, l_k).' * y_aux(l_k, l_k);
    
    
    z_sum{k, M-1} = zeros(1, Nk(k));
end

%% Simulation of Diff-RLS N-LMS.ATC Algorithm

for i = M:Tmax-1
    
    % Save psi(n-1)         
    psi_prev = psi;
    w_prev = w;
    
    % get what nodes will be disconnectes in this iteration
    dead_nodes = get_dead_nodes(dead_nodes, N, i, error_params);
    
    
    %update the w only of the connected nodes
    w_prev(:, dead_nodes == 0 ) = w(:, dead_nodes == 0);
    
    %% 0 - FILTER
    u_k = u(:, i:-1:i-M+1);
    
    % Filter each input u_k with all the available weights w(n-1) for n_bar_k 
    % and for k : psi (n-1)
    Y_kl = (u_k * w_prev) .* (A - eye(N)).' + (u_k * psi_prev) .* eye(N);
   
    % Error computation
    for k = 1:N % for each node k
        ecomb_i(k) = d(k, i) - u_k(k,:) * w(:,k);
    end
   
    %% 1- ADAPT
    [y_k , e_k,  psi] = adapt_nlms(N , M , d(:,i), u_k, mu, psi_prev);
    
    
    %% UPDATE WEIGHTS AND COMBINE
    [w, c, P_sum, z_sum] = update_combine_ls_exp(N, A, psi,...
        w_prev, e_k, Y_kl, y_k, P_sum, z_sum, gamma, regul, i);
    
      
    
    % save results  
    msd_comb(:, i) = sum((w0(:, :, i) - w_prev).^2);
    c_out(:, :, i) = c;
    e_comb(:, i) = ecomb_i;
 
    if is_returned_w == 1
        w_out(:,:,i) = w;
    else
        w_out = [];
    end
    
end


