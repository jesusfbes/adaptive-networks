function [msd, errors, c, w0, u, v, d, w_out] = sim_an(algorithms,...
    Tmax, n_nodes, A, sigma2_u, SNR, w0_1, w0_2, mu, params, error_params, ...
    is_returned_w)
%SIM_AN Main function that simulates the diffusion network in stationary environments.
%
%
% It generates the signals, execute the algorithms and return errors and input signals
%
% INPUT algorithms: Cell array with the names of the algorithms to
% simulate.
%       Tmax: Number of iterations.
%       n_nodes: Number of nodes of the network.
%       A: Adjacency matrix of the network.
%       sigma2_u: Variance of the regressor input u.
%       SNR: in db of the input (x) vs. the noise v.
%       w0_1: vector of unknown parameters to estimate (Mx1)
%       w0_2: vector of unknown parameters to estimate for the second part of the simulation, if no change
%       you must leave it empty ([]).
%       mu: unormalized stepsize for the NLMS filters in all the
%       algorithms.
%       params: Parameters for the algorithm, different meaning for different
%       algorithm.
%       error_params. Parameters for node errors models.
%       is_returned_w. FLAG to indicate if w_out should be returned (big data)
%
%
% OUTPUT msd_comb. MSD using combined estimations (not used in ACW )%       
%        errors. error of the estimation (use outside to compute EMSE) 
%        c. combination coefficients evolution 
%        w0. unknown parameter to estimate in each time step.  M x N_nodes x TMax
%        u. input regressors
%        v. additive noise signal 
%        d. input measurement signal 
%        w_out. estimation of w_0 (MxNxTmax). Only return it if neaded
%
%
%
% See Also ATC_NLMS_NOCOOP, SIM_AN_TRACK, ACW, ADN_LS, ADN_HYBRID
%
% created with MATLAB ver.: 8.0.0.783 (R2012b) on Mac OS X  Version: 10.9.4
%
% created by: Jesus Fernandez-Bes (<a href="http://www.tsc.uc3m.es/~jesusfbes">web</a>)
% DATE: Oct-2014


% prealocation for filter inputs and output
y = zeros(n_nodes, Tmax);
u = zeros(n_nodes, Tmax);
v = zeros(n_nodes, Tmax);

% Prealocation of results

msd = struct;
errors = struct;
c = struct;



for a = 1:length(algorithms)
    
    algorithm = algorithms{a};
    
    msd.(algorithm) = []; 
  
    errors.(algorithm) = []; 

    c.(algorithm) = [];     
    
    w_out.(algorithm) = []; 
    
end


%%% Variance per each node

if(length(sigma2_u) == 1)
    
    sigma2_u = ones(1, n_nodes) * sigma2_u;
    
elseif(length(sigma2_u) ~= n_nodes)
    error('var_u has an incorrect size');
    
end
 
%%% SNR per each node

if(length(SNR) == 1)
    SNR = ones(1, n_nodes) * SNR;
    
elseif(length(SNR) ~= n_nodes)
    error('SNR has an incorrect size');
    
end

%% Simulation


%% Generate the signals of each simulation
for k = 1:n_nodes
    
    u(k, :) = normrnd(0, sqrt(sigma2_u(k)), 1, Tmax);
    
    sigma2_v_k = (power(10, (-SNR(k) / 10))) * sigma2_u(k);
    v(k, :) = normrnd(0, sqrt(sigma2_v_k), 1, Tmax);
    
    
    if(isempty(w0_2)) % no filter change
        
        y(k, :) = filter(w0_1, 1, u(k, :));
        w0 = repmat(w0_1, [1 n_nodes Tmax]);
        
    else % filter change at the middle of the simulation
        y(k, :) = filter_with_change(w0_1, w0_2, u(k, :));
        w0(:, :, 1:floor(Tmax/2)) = repmat(w0_1, [1 n_nodes floor(Tmax/2)]);
        w0(:, :, floor(Tmax/2)+1:Tmax) = repmat(w0_2, [1 n_nodes floor(Tmax/2)]);
    end
    
end

d = y + v;

%% Simulate selected algorithms
for a = 1:length(algorithms)
    
    algorithm = algorithms{a};
    disp(algorithm);
    
    switch lower(algorithm)
       
        case 'atc_nlms_nocoop' 
            [msd_aux, e_aux,c_aux, w_aux] = ...
                atc_nlms_nocoop(d,u,w0,mu,A,n_nodes,Tmax,error_params, is_returned_w);
            
         case 'atc_nlms_metropolis' 
            [msd_aux, e_aux,c_aux, w_aux] = ...
                atc_nlms_metropolis(d,u,w0,mu,A,n_nodes,Tmax,error_params, is_returned_w);
            
        case 'atc_nlms_acw' 
           
            nu = params.atc_nlms_acw.nu ;
            
            [msd_aux, e_aux, c_aux, w_aux] = ...
                atc_nlms_acw(d,u,w0,mu,A,n_nodes,Tmax,nu,error_params, is_returned_w);
        
%         case 'datc_nlms_ls_rect'
%             L = params.datc_nlms_ls_rect.L ;
%             regul = params.datc_nlms_ls_rect.regul;
%             
%             [msd_aux, e_aux, c_aux, w_out] = ...
%                 datc_nlms_ls_rect(d, u, w0, mu, A, n_nodes, L, regul, Tmax, error_params, is_returned_w);
            
         case 'datc_nlms_ls_exp'
            gamma = params.datc_nlms_ls_exp.gamma ;
            regul = params.datc_nlms_ls_exp.regul;
            
            [msd_aux, e_aux, c_aux, w_out] = ...
                datc_nlms_ls_exp(d, u, w0, mu, A, n_nodes, gamma, regul, Tmax, error_params, is_returned_w);
            
         
        otherwise
            error('Unknown diffusion algorithm.')
            
    end
    
    msd.(algorithm) = msd_aux;    
    errors.(algorithm) = e_aux;
    c.(algorithm) = c_aux;

    if is_returned_w
        w_out.(algorithm) = w_aux;
    end
    
end



