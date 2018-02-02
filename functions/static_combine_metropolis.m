function [w , c] = static_combine_metropolis(N, M, A, psi, w_prev )
%UPDATE_COMBINE_ACW. Update the combiners and use them in Diffusion
%Networks.
%
%
% It computes the combination weights according to metropolis rule and
% combine accordingly:
%Cattivelli, F. S., & Sayed, A. H. (2010). Diffusion LMS strategies 
%for distributed estimation. IEEE Transactions on Signal Processing, 
%58(3), 1035-1048.
% 
% For details in the algorithm, please see [1].
%
% INPUT N: number of nodes in the network.
%       M: number of coefficients in w.
%       A: adjacency matrix of the network.
%       psi: vector with local estimations of the nodes.
%       w_prev: matrix with previous combined estimations of all the nodes.
%           
%
% OUTPUT w: matrix with new combined estimations of all the nodes.
%        c: matrix with new combination weights of all the nodes.
%
%
% See Also ATC_NLMS_METROPOLIS
%
% created with MATLAB ver.: 8.0.0.783 (R2012b) on Mac OS X  Version: 10.9.4
%
% created by: Jesus Fernandez-Bes (<a href="http://jesusfbes.es">web</a>)
% DATE: Feb-2018

w = zeros(M,N) ;

n_k = sum(A == 1) ;

c_aux = (A-eye(N)).*(1./max(repmat(n_k,N,1),repmat(n_k.',1,N))) ;

c = c_aux + diag( 1 - sum(c_aux) ) ;


for k = 1:N
    
    c_bar_k = c(:,k ) ;
    c_bar_k(k) = 0 ;
    
    if (~isempty(w_prev) && ~isempty(psi))
        w(:,k) =  w_prev*c_bar_k + c(k,k)*psi(:,k) ;
    end
    
end

end