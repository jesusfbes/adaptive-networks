%% Example of execution of Adaptive Diffusion Networks
%
%

clc; clear; close all;


%% Location of code
addpath('functions');

%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FLAG_RETURN_W = 0;% FLAG if 1 returns all the evolution of w


Tmax = 20000; % Number of iterations

n_sim = 2; % Number of simulations to average


% We load the network

load('inputs/example_complex.mat');
%load('inputs/example_basic.mat');


% Algorithms to execute
algorithms = { 'atc_nlms_nocoop', 'atc_nlms_acw', 'le_atc_nlms_ls' };

% Algorithm parameters
params.atc_nlms_acw.nu = 0.01; % learning parameter for the combination

params.le_atc_ls.L = 100; % Window size for combination estimation


% Error model parameters (no errors)
error_param.mode = 'none';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_nodes = size(A,1); % Number of nodes

display(A);  % display Adjacency Matrix




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization of output structures
for a = 1:length(algorithms)
    
    algorithm = algorithms{a};
    
    msd_array.(algorithm) = 0;
    emse_array.(algorithm) = 0;
    c.(algorithm) = 0;
    
end

for iter = 1:n_sim
    
    disp(['Simulation ' num2str(iter)]);
    
    [msd, errors, c_aux, w0, u, v, d] = sim_an( algorithms, Tmax, ...
        n_nodes, A, sigma2_u, snr, w0_1, w0_2, mu_filter, params, error_param,...
        FLAG_RETURN_W);
    
    
    % Accumulate results
    for a = 1:length(algorithms)
        
        algorithm = algorithms{a};
        msd_array.(algorithm) = msd.(algorithm) + msd_array.(algorithm);
        emse_array.(algorithm) = ...
            (errors.(algorithm)  - v).^2 + emse_array.(algorithm);        
        c.(algorithm) = c_aux.(algorithm) + c.(algorithm);
        
    end
    
    disp(' ');
    
end


% Divide to obtain averages
for a = 1:length(algorithms)
    
    algorithm = algorithms{a};
    msd_array.(algorithm)= msd_array.(algorithm) ./ n_sim;    
    emse_array.(algorithm) = emse_array.(algorithm)./n_sim;
    c.(algorithm) = c.(algorithm) ./ n_sim;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

algorithms_plot = ...
    cellfun(@get_algorithm_name_plot, algorithms, 'UniformOutput', 0);

colors = 'rgbmcky';
figure(1)
for a = 1:length(algorithms)
   algorithm = algorithms{a};
   
   network_msd = mean(msd_array.(algorithm)) ; 
   
   plot(10*log10(network_msd),colors(a)); hold on;
   
    
end
title('NETWORK MSD (dB)');
legend(algorithms_plot);

figure(2)
for a = 1:length(algorithms)
   algorithm = algorithms{a};
   
   network_emse = mean(emse_array.(algorithm)) ; 
   plot(10*log10(network_emse),colors(a)); hold on;
    
end
title('NETWORK EMSE (dB)');
legend(algorithms_plot);

%% Combiners
figure(3);

node_to_show = 1;
for a = 1:length(algorithms)
   algorithm = algorithms{a};
   c_to_plot = c.(algorithm);
   
   subplot(1,length(algorithms),a);
   plot(squeeze(c_to_plot(:, node_to_show,:)).'); hold on;
   axis([1 Tmax 0 1]);
   title(['Combiners of ' algorithms_plot{a} ' for node' num2str(node_to_show)]);
    
end


