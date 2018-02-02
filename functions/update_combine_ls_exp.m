function [w, c, P_sum, z_sum] = ...
    update_combine_ls_exp(N, A, psi, w_prev, e_k, Y_kl, y_k, P_sum, z_sum, ...
    gamma, regul, iter)
%UPDATE_COMBINE_LS. Update the combiners and use them in Diffusion
%Networks.
%
%
% It computes the adaptive combination weights of the ADN-LS algorithm
% according to paper:
%   [1]  Fernandez-Bes, J., Azpicueta-Ruiz, L. A., Arenas-García, J., & Silva, M. T. (2014). 
%   Distributed estimation in diffusion networks using affine least-squares combiners. 
%   Digital Signal Processing.
% 
% 
% For details in the algorithm, please see [1].
%
% INPUT N: number of nodes in the network.
%       A: adjacency matrix of the network.
%       psi: vector with local estimations of the nodes.
%       w_prev: matrix with previous combined estimations of all the nodes.
%       e_k: prediction error of this iteration, computed with previus psi.
%       Y_kl: predicted output using the estimations received by the
%       neighbours.
%       y_k: predicted output using previous local estimation psi.
%       Y_tilde_list: list of the previous predicted outputs in the
%       estimation window.
%       z_list: list of the previous z's in the estimation window.
%       P_sum: Acummulated matrix P to compute the combiners.
%       z_sum: Acummulated vector z to compute the combiners.
%       iter: index of the iteration.
%           
%
% OUTPUT w: matrix with new combined estimations of all the nodes.
%        c: matrix with new combination weights of all the nodes.
%        P_sum: New acummulated matrix P to compute the combiners.
%        z_sum: New acummulated vector z to compute the combiners.
%        Y_tilde_list: new list of the previous predicted outputs in the
%       estimation window.
%        z_list: new list of the previous z's in the estimation window.
%
%
% See Also ADN_LS
%
% created with MATLAB ver.: 8.0.0.783 (R2012b) on Mac OS X  Version: 10.9.4
%
% created by: Jesus Fernandez-Bes (<a href="http://jesusfbes.es">web</a>)
% DATE: Oct-2014



% Initialize variables

% y_bk-y_k (filter with psi for all nodes.
Y_tilde = ( Y_kl - y_k*ones(1,N) ) .* A.' ;

c = zeros(N, N);
w = w_prev ;

% y_bk-y_k (filter with psi for all nodes.
Y_tilde = (Y_kl - y_k*ones(1, N)) .* A.' ;


for k = 1:N
    
    % We update based in the estimation of each l-neighbour
    % Neighbours are the nodes connected to k
    l_k = (A(:, k)==1);
    l_k(k) = 0;
    
    
    %%% P udpate
    
    
    %calculate the last term in the window
    P_aux = Y_tilde(k, l_k).' * Y_tilde(k, l_k);
  
    
    % Calculate the current P value
    P = gamma*(P_sum{ k, iter - 1 } ) +(1-gamma)* P_aux  ;
    
    
    P_sum{k, iter} = P;
    
    
    
    %%% z update
    
    %calculate the last term in the window
    z_aux = e_k( k ) .* Y_tilde( k , l_k ) ;
    
      
    z = gamma*(z_sum{ k, iter - 1  } ) + (1-gamma)* z_aux ;

    z_sum{k, iter} = z;
    
    
    
    % Calculate the new combination weights
    
    R = chol(P + regul*eye(sum(l_k)) ) ;  
    c(l_k, k) = R \ (R'\ (z).') ;
    
    c(k, k) = 1 - sum(c(l_k, k));
    
    c(c(:,k) < 0,k ) = 0 ;
    
    
    c(:,k)  = c(:,k) ./  sum(c(:,k)) ;
    
    %% 3. - COMBINE
        
    % Combine weights   
    w(:, k) =  w_prev(:, l_k) * c(l_k, k) + c(k, k) * psi(:, k);
    
end

end