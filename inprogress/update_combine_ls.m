function [w, c, P_sum, z_sum, Y_tilde_list, z_list] = ...
    update_combine_ls(N, A, psi, w_prev, e_k, Y_kl, y_k, Y_tilde_list, z_list,...
    P_sum, z_sum, iter)
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
% created by: Jesus Fernandez-Bes (<a href="http://www.tsc.uc3m.es/~jesusfbes">web</a>)
% DATE: Oct-2014



% Initialize variables

old_Y_tilde = Y_tilde_list(:, :, 1);

old_Z = z_list(:, :, 1);

z_new = zeros(N, N);

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
    
    
    % Add it to the other L
    old_P = old_Y_tilde(k, l_k).' * old_Y_tilde(k, l_k);
    
    
    % Calculate the current P value
    P = P_sum{k, iter - 1} - old_P + P_aux;
    
    
    P_sum{k, iter} = P;
    
    
    
    %%% z update
    
    %calculate the last term in the window
    z_aux = e_k( k ) .* Y_tilde( k , l_k ) ;
    
    % Add it to the other L
    old_z = old_Z( k , l_k ) ;
    
    
    % Calculate the current z value
        
    z = z_sum{k, iter - 1} - old_z + z_aux;
   
    z_sum{k, iter} = z;
    
    z_new (k, l_k) = z_aux;
    
    
    % Calculate the new combination weights
    
    
    c(l_k, k) = P \ (z).';
    
    c(k, k) = 1 - sum(c(l_k, k));
    
    
    %% 3. - COMBINE
        
    % Combine weights   
    w(:, k) =  w_prev(:, l_k) * c(l_k, k) + c(k, k) * psi(:, k);
    
end


% Saving the output in the ouput structures

Y_tilde_list(:, :, 1:end-1) = Y_tilde_list(:, :, 2:end);

Y_tilde_list(:, :, end) = Y_tilde;

z_list(:, :, 1:end-1) = z_list(:, :, 2:end);

z_list(:, :, end) = z_new;




end