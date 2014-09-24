function [msd_comb msd_sep e_comb c_out w_out] = ...
    adn_hybrid(d,u,w0,mu,A,N,L,Tmax,mu_a,error_params,is_returned_w, debug)
%ADN_HYBRID Function that simulates the Diffusion NLMS filter
% with adaptive coefficients (calculated using Least-Squares ). The adaptation
% is done over the loal estimation but we shared the combined weights
%
%INPUTS
%d. desired signal NxT
%u. Input vector NxT
%w0. Vector to estimate MxNxT
%mu. Adaptation step
%N. Number of nodes
%L. Square Window size
%Tmax. Number of iterations of the algorithm
%
%OUTPUTS
%msd_comb. MSD using the combination  of each node NxT
%msd_sep. MSD of each individual filter NxT
%e_comb. combination errors NxT
%e_partial. Partial errors NxNxT
%c_out. Mixing parameters NxNxT

%% Parameters

% number of coefficients
M = size(w0,1) ;



delta = 1e-2 ;


dead_nodes = zeros(N,1) ;

% Initialize filter coefficients
psi_ls = zeros(M,N) ;
w_ls = zeros(M,N) ;
w_prev_ls = w_ls ;


% Parameters ACW

nu = 0.01 ;

%psi_acw = zeros(M,N) ;
w_acw = zeros(M,N) ;
w_prev_acw = w_acw ;

w = w_ls ;

gamma2 = zeros(N,N) ;

% parameters of output combination
beta = 0.9 ;

a_plus = 4 ;
p = ones( N , 1 ) ;
a = a_plus*ones( N , 1 ) ;
alpha = zeros( N , 1 ) ;

cte_a = 1 ./ ( logsig( a_plus )- logsig(- a_plus ) ) ;


P_sum = cell( N , Tmax ) ;
z_sum = cell( N , Tmax ) ;

Y_tilde_list = zeros( N , N , L ) ;
z_list = zeros(N , N , L ) ;

% Preallocate results
msd_comb=zeros(N,Tmax) ;
msd_sep=zeros(N,Tmax) ;
c_out=zeros(N,Tmax);

e_comb=zeros(N,Tmax);
e_partial=zeros(N,Tmax);


Nk = sum(A) - 1 ;

if is_returned_w ==1
    w_out = zeros(M,N,Tmax) ;
end

for k=1:N
    % We update based in the estimation of each l-neighbour
    % Neighbours are the nodes connected to k
    
    l_k = ( A(:,k)==1 ) ;
    
    l_k(k)=0;
    
    y_aux = zeros(N) ;
    y_aux( l_k , l_k ) = sqrt( delta ) * eye( Nk(k) ) ;
    
    
    Y_tilde_list( : , : , 1 ) =  y_aux ;
    
    
    
    P_sum{ k , M-1 } = y_aux( l_k , l_k ).' * y_aux( l_k , l_k )  ;
    
    
    z_sum{ k , M-1  } = zeros( 1 , Nk(k) ) ;
end

%% Simulation of Diff-RLS N-LMS.ATC Algorithm

for i = M:Tmax-1
    
    % Save psi(n-1) and w(n-1)
    
    
    psi_prev_ls = psi_ls ;
    
    dead_nodes = get_dead_nodes(dead_nodes,N,i, error_params) ;
    
    
    %update the W of the alive nodes
    w_prev_ls(:,dead_nodes == 0 ) = w_ls(:,dead_nodes == 0 ) ;
    
    w_prev_acw(:,dead_nodes == 0 ) = w_acw(:,dead_nodes == 0 ) ;
    
    %% 0 - FILTER
    u_k = u( : , i:-1:i-M+1 ) ;
    
    % Filter each input u_k with all the available weights w(n-1) for n_bar_k
    % and for k : psi (n-1)
    Y_kl_ls = ( u_k*w_prev_ls ) .*( A ).'  ;
   % Y_kl_ls = ( u_k*w_prev_ls ) .*(A - eye(N)).' + (u_k*psi_prev_ls).*eye(N) ;
 
    %% COMBINATION
    
    y_ls = diag(Y_kl_ls) ;
    
    Y_kl_acw = ( u_k*w_prev_acw ) .*A.' ;
    y_acw = diag( Y_kl_acw ) ;
    
    
    y_out = alpha.*y_ls + (1-alpha).* y_acw ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Combination Coefficient learning
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y = y_ls - y_acw ;
    
    e_alpha = d(:, i ) - y_out ;
    
    p = beta * p + (1-beta)*e_alpha.^2 ;
    
    der_a = ( logsig(a).*(1- logsig(a)) ) .* cte_a  ;
    
    a  = a + (mu_a./p).*e_alpha.*y.*der_a ;
    
    a = max( min( a , a_plus ) , -a_plus );
    
    alpha = ( logsig(a) - logsig(-a_plus) ) .* cte_a ;
    

    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Hybrid output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % combination coefficient learning
    alpha_aux = alpha ;
    alpha_aux(alpha > 0.99 ) = 1 ;
    
    for k = 1:N
       
        w(:,k) = alpha_aux(k) .* w_ls(:,k) + ( 1- alpha_aux(k) ).*w_acw(:,k) ;
        
    end
    
    
    
    ecomb_i = d(:, i ) - diag( Y_kl_ls ) ;
    
    %% 1- ADAPT
    [y_k , e_k_ls,  psi_ls] = adapt_nlms(N , M , d(:,i), u_k, mu, psi_prev_ls ) ;
    
    
    %% UPDATE WEIGHTS AND COMBINE
    [w_ls ,c ,P_sum, z_sum ,Y_tilde_list, z_list] = update_combine_ls( N , A , psi_ls, w_prev_ls, e_k_ls, ...
        Y_kl_ls, y_k, Y_tilde_list, z_list , P_sum , z_sum, i ) ;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ACW
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save psi(n-1) and w(n-1)
    %psi_prev = psi ;
   
    
    
    %% 1- ADAPT
    [ dummy , e_k_acw , psi_acw] = adapt_nlms(N , M , d(:,i), u_k, mu, w_prev_acw ) ;
    
    
    %% 2 .- COMBINATION WEIGHTS UPDATE
    
    [ w_acw c_acw gamma2] = update_combine_acw( N , M, A, psi_acw ,w_prev_acw, gamma2 , nu ) ;
    
    
    
    
    % save results
    msd_sep(:,i) = sum( (w0(:,:,i)-psi_ls).^2 ) ;
    msd_comb(:,i) = sum( (w0(:,:,i)-w).^2 ) ;
    c_out(:,i) = alpha ;
    e_comb(:,i) = ecomb_i ;
    e_partial(:,i) = e_k_ls ;
    
    if is_returned_w == 1
        w_out(:,:,i) = w ;
    else
        w_out = [] ;
    end
    
end


if debug
    
    %for k=1:N
    figure;
    subplot(311);
    plot(1./reshape([P_sum{4,M:Tmax-1}],(Nk(4)*Nk(4)),Tmax-M).');
    xlim([M Tmax-1]);
    title('P^{-1}') ;
    
    subplot(312);
    plot(reshape([z_sum{4,M:Tmax-1}],( Nk(4) ),Tmax-M).');
    xlim([M Tmax-1]);
    title('z') ;
    
    subplot(313);
    plot(  reshape( c_out( :, 4 , M:Tmax-1) , N , Tmax-M ).'  ) ;
    axis([M Tmax-1 -0.5 1.5]);
    
    title('c') ;
    legend
    %end
end
