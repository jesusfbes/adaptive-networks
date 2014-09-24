clc; clear ; close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Simulation Parameters

DEBUG = 0 ; % FLAG if 1 shows debug figures and messages
FLAG_RETURN_W = 0 ;% FLAG if 1 returns all the evolution of w


algorithms = {'adn_ls','adn_hybrid','acw'} ;

Tmax = 60000 ; % Number of iterations

n_sim  = 1 ; % Number of simulations

%% Load Network
% We load the network
load Atak ;

n_nodes = size(A,1) ; % Number of nodes

display(A);  % display Adjacency Matrix

%% Signal Parameters
params.adn_ls = 100 ; %LS square window size
params.adn_hybrid = [100 1] ;
error_param.mode = 'random' ;
error_param.P = 0.5 ;
%error_param.mode = 'none' ;

%error_param.mode = 'fixed' ;
%error_param.T = floor( Tmax / 4 ) ;
%error_param.percent = 0.1 ;

var_u = 1 ; % variance of input signals u



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization of output structures
for a=1:length(algorithms)
    
   algorithm=algorithms{a};
    
    Msd_comb_array.(algorithm)=0;
    Msd_sep_array.(algorithm)=0;
    Emse_comb.(algorithm)=0;    
    c.(algorithm)=0;
    
end

for iter = 1:n_sim
    
    
    [msd_comb msd_sep e_comb c_aux w0 u v d w_out] = simAN( ...
        algorithms,Tmax,n_nodes,A,var_u,snr,w0_1,w0_2,mu_filter,params,error_param,FLAG_RETURN_W,DEBUG) ;
    
    
%     %Accumulate results
    for a=1:length(algorithms)
        
        algorithm=algorithms{a};
        Msd_comb_array.(algorithm)=msd_comb.(algorithm) + Msd_comb_array.(algorithm);
        Msd_sep_array.(algorithm)=msd_sep.(algorithm) + Msd_sep_array.(algorithm);
        Emse_comb.(algorithm)=( e_comb.(algorithm)  - v ).^2 + Emse_comb.(algorithm);       
        c.(algorithm)=c_aux.(algorithm) + c.(algorithm);
        
    end
    
end


% Divide to obtain average
for a=1:length(algorithms)
    
    algorithm=algorithms{a};
    Msd_comb_array.(algorithm)= Msd_comb_array.(algorithm) ./ n_sim ;
    Msd_sep_array.(algorithm)= Msd_sep_array.(algorithm)./ n_sim ;
    Emse_comb.(algorithm) = Emse_comb.(algorithm)./n_sim ;
    c.(algorithm) = c.(algorithm) ./ n_sim ;
end

h_nmsd = show_nmsd_figure( Msd_comb_array , Msd_sep_array , algorithms , Tmax ) ;

%h_emse = show_nmsd_figure( Emse_comb , Emse_comb , algorithms , Tmax ) ;

%save( [ 'paraAzpi_basica' num2str(iter) '.mat'] ) ;
save('results/Perror0_5_mu_a_1.mat') ;


