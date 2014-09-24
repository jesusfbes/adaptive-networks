function [msd_comb msd_sep e_comb c w0 u v d w_out] = ...
    simAN_track(algorithms,Tmax,Nnodes,A,var_u,SNR,w0_1,Tr_Q,mu,params,is_returned_w,debug)
%SIM_AN_TRACK It is the main function that simulates the algorithms of
%distributed estimation on diffusion adaptive networks in tracking enviroments. 
%It generates the signals, execute the algorithms and return errors and input signals
%
%INPUTS:
%algorithms. Cell array with the names of the algorithms
%Tmax. Number of iterations
%Nnodes. Number of nodes
%A. Adjacency matrix
%var_u. Variance of the regressor input
%SNR. in db
%w0_1. unknow parameter to estimate M x 1
%Tr_Q. Trace of covariance matrix of perturbations in w_0 to track
%mu. stepsize for the NLMS filters in all the algorithms
%params. PArameters for the algorithm, different meaning for different
%       users.
%is_returned_w. FLAG to indicate if w_out should be returned (big data)
%debug. FLAG to show debug messages
%
%OUTPUTS:
%msd_comb.MSD using combined estimations (not used in ACW )
%msd_sep. MSD using not combined estimations (not used in ACW )
%e_comb. error of the estimation (use outside to compute EMSE) 
%c. combination coefficients evolution 
%w0. unknown parameter to estimate in each time step.  M x N_nodes x TMax
%u. input regressors
%v. additive noise signal 
%d. input measurement signal 
%w_out. estimation of w_0 (MxNxTmax). Only return it if neaded
%
%author: Jesus Fernandez-Bes

%% PARAMETERS 
% Gamma in the dynamic model
gamma = 0.99 ;

% number of coefficients
M = size(w0_1,1) ;

% prealocation for filter inputs and output
y = zeros(Nnodes,Tmax) ;
u = zeros(Nnodes,Tmax) ;
v = zeros(Nnodes,Tmax) ;
w0 = zeros(M,Nnodes,Tmax) ;


% Prealocation of results

msd_comb = struct ;
msd_sep = struct ;
e_comb = struct;
%e_partial = struct ;
c = struct ;

% Initialization of output structures
for a=1:length(algorithms)
    
    algorithm=algorithms{a};
    
    msd_comb.(algorithm) = [] ;
    msd_sep.(algorithm) = [] ;
    e_comb.(algorithm) = [] ;
    %   e_partial.(algorithm)=0;
    c.(algorithm) = [] ;
    
    w_out.(algorithm) = [] ;
    
    
end


%%% Variance per each node

if length(var_u) == 1
    
    var_u=ones(1,Nnodes)*var_u;
    
elseif length(var_u) ~= Nnodes
    error('var_u has an incorrect size');
    
end

%%% SNR per each node

if length(SNR) == 1
    SNR=ones(1,Nnodes)*SNR;
    
elseif length(SNR) ~= Nnodes
    error('SNR has an incorrect size');
    
end

%% Simulation


%init buffer
ubuffer = zeros(1,M);

%% Generate the signals of each simulation
for k=1:Nnodes
    
    u(k,:) = normrnd(0,sqrt(var_u(k)),1,Tmax);
    v(k,:)= normrnd(0,sqrt((power(10,(-SNR(k)/10)))*var_u(k)),1,Tmax);
    
    for i=1:Tmax-1
        ubuffer=[u(k,i) ubuffer];
        ubuffer(end)=[];
        
        %  if(isempty(w0_2))
        q = normrnd(0,sqrt(Tr_Q/M),M,1);
        
        if(i==1)
            w0(:,k,i)=q+w0_1;
        else
            w0(:,k,i) = w0_1+ gamma*(w0(:,k,i-1)-w0_1)+q;
        end
        
        
        y(k,i)=ubuffer*w0(:,k,i);
        
    end
end

d = y + v ;

for a=1:length(algorithms)
    
    algorithm=algorithms{a};
    disp(algorithm);
    switch lower(algorithm)
        
        case 'adn_ls'
            
            L = params.adn_ls ;
            
            %% RLS Algorithm (combination weights shared)
            [msd_comb_aux msd_sep_aux e_comb_aux c_aux w_aux] =...
                adn_ls(d,u,w0,mu,A,Nnodes,L,Tmax,is_returned_w,debug);
            
            % case 'adn_ls_transfer'
            %
            %     %% RLS Algorithm (combination weights shared)
            %     [msd_comb_aux msd_sep_aux e_comb_aux e_partial_aux c_aux] =...
            %         adn_ls_transfer(d,u,w0,mu,A,Nnodes,L,Tmax,debug);
            
            
        case 'acw'
            %% ADAPTIVE COMBINATION WEIGHTS Zhao and Sayed
            [msd_comb_aux msd_sep_aux e_comb_aux c_aux w_aux] =...
                acw(d,u,w0,mu,A,Nnodes,Tmax,is_returned_w,debug);
            
            
            %  case 'acw_notransfer'
            %      %% ATC Algorithm (shared, original of Takahashi)
            %      [msd_comb_aux msd_sep_aux e_comb_aux e_partial_aux c_aux] =...
            %          acw_notransfer(d,u,w0,mu,A,Nnodes,Tmax);
            
        case 'adn_hybrid'
            
            L = params.adn_ls(1) ;
            mu_a = params.adn_hybrid(2) ;
            [msd_comb_aux msd_sep_aux e_comb_aux c_aux w_aux] =...
                adn_hybrid(d,u,w0,mu,A,Nnodes,L,Tmax,mu_a,is_returned_w,debug);
            
        otherwise
            error('Unknown diffusion algorithm.')
    end
    
    msd_comb.(algorithm) = msd_comb_aux ;
    msd_sep.(algorithm) = msd_sep_aux ;
    e_comb.(algorithm) = e_comb_aux ;
    c.(algorithm) = c_aux ;
    
    if is_returned_w
        w_out.(algorithm) = w_aux ;
    end
    
end



