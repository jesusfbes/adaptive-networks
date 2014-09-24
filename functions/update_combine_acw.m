function [w, c, gamma2] = update_combine_acw(N, M, A, psi, w_prev, gamma2, nu)
%UPDATE_COMBINE_ACW. Update the combiners and use them in Diffusion
%Networks.
%
%
% It computes the adaptive combination weights of the ACW algorithm
% according to paper:
%   [1]	S.-Y. Tu, A.H. Sayed
%	Optimal combination rules for adaptation and learning over networks
%   Proc. of 4th IEEE Intl. Workshop on Computational Advances in Multi-Sensor
%   Adaptive Process, CAMSAP, San Juan, Puerto Rico (2011), pp. 317?320
%
% 
% For details in the algorithm, please see [1].
%
% INPUT N: number of nodes in the network.
%       M: number of coefficients in w.
%       A: adjacency matrix of the network.
%       psi: vector with local estimations of the nodes.
%       w_prev: matrix with previous combined estimations of all the nodes.
%       gamma_2: vector of estimated sigma_2 of each node.
%       nu: learning parameter for the combiners.
%           
%
% OUTPUT w: matrix with new combined estimations of all the nodes.
%        c: matrix with new combination weights of all the nodes.
%        gamma_2: new vector of estimated sigma_2 of each node.
%
%
% See Also ACW
%
% created with MATLAB ver.: 8.0.0.783 (R2012b) on Mac OS X  Version: 10.9.4
%
% created by: Jesus Fernandez-Bes (<a href="http://www.tsc.uc3m.es/~jesusfbes">web</a>)
% DATE: Oct-2014

% Initialize outputs
c = zeros(N, N);
w = zeros(M, N);

for k = 1:N % for each node k
        
        % list of neighbours of k
        l_k = (A(:, k) == 1); 
             
                     
        %eq 34 of [1]       
        error_aux = psi(:, l_k) - w_prev(:, k) * ones(1, sum(l_k));
       
       
        gamma2(l_k, k) = (1 - nu) * gamma2(l_k, k) ...
                        + nu * diag((error_aux)' * (error_aux));
        
        %eq 36 of [1]
        sum_gammainv2 = sum(1 ./ gamma2(l_k, k));
        
        c(l_k, k) = 1 ./ gamma2(l_k, k) ./ sum_gammainv2;
        
        
        %% 3. - Combination
        w(:, k) = psi * c(:, k);
        
end
    

end

