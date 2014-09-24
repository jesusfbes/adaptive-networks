function [msd_comb, e_comb, c_out, w_out] = ...
    atc_nlms_nocoop(d, u, w0, mu, A, N, Tmax, error_params, is_returned_w)
%ATC_NLMS_NOCOOP Simulates a Diffusion network with ATC but without cooperation
%
%
% 
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
% OUTPUT msd_comb:
%        e_comb:
%        c_out:
%        w_out:
%
%
%
% created with MATLAB ver.: 8.0.0.783 (R2012b) on Mac OS X  Version: 10.9.4
%
% created by: Jesus Fernandez-Bes (<a href="http://www.tsc.uc3m.es/~jesusfbes">web</a>)
% DATE: Oct-2014


%% Parameters


% number of coefficients
M = size(w0,1);

% Initialize filter coefficients

w = zeros(M,N) ;

w_prev = w;

% Combination matrix
c = eye(N) ; 

% Initialize combination coefficientes and intermediate variables

dead_nodes = zeros(N, 1);

Nk = sum(A) ; % Number of neighbours per node

b = cell(N,1) ;
ones_Nk = cell(N,1) ;


for k=1:N
    
    ones_Nk{k}=ones(Nk(k),1);
    b{k}= ones(Nk(k),1)/Nk(k);
    
end

% Preallocate results
msd_comb = zeros(N, Tmax);
c_out = zeros(N, N, Tmax);
e_comb = zeros(N, Tmax);


if(is_returned_w == 1)
    w_out = zeros(M, N, Tmax);
else
    w_out = [];
end


%% Simulation of N-LMS non-cooperative ATC Algorithm

for i=M:Tmax-1
    
    dead_nodes = get_dead_nodes(dead_nodes, N, i, error_params);
    
    w_prev(:, dead_nodes == 0) = w(:, dead_nodes == 0);
    
    
    %% 1- FILTER and ADAPT
    u_k =u(:, i:-1:i-M+1) ;
    
   
    [dummy , e_k , psi] = adapt_nlms(N , M , d(:, i), u_k, mu, w_prev ) ; %#ok<*ASGLU>
    
    
    %% 2 .- COMBINATION WEIGHTS UPDATE
    

    w = psi ;

    
    % save results
    msd_comb(:, i) = sum((w0(:, :, i) - w_prev).^2);
    
    c_out(:, :, i) = c;
    e_comb(:, i) = e_k;
       
    if(is_returned_w == 1)
        w_out(:, :, i) = w ;
    end
    
end
