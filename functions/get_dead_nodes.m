function dead_nodes = get_dead_nodes(dead_nodes,N,iter, param)
%GET_DEAD_NODES Generates the vector of nodes that do not respond in this
%iteration.
%
%
% It computes which nodes will be able to perform an estimation and
% comunicate in this iteration. There are now 3 nodes: 
%   'none': no errors at all.
% 	'random': At every iteration nodes get disconnect with probability P.% 	
%	'fixed': At iteration T some nodes are disconnected forever.
%
% INPUT dead_nodes: binary vector with the flags 1 if dead, 0 if alive in the
% previous iteration. 
%       N: Number of nodes
%       iter. iteration index
%       param: struct with params of the node error model. it always
%       constains a param `mode' and maybe other parameters depending on
%       the mode.
%           
%
% OUTPUT dead_nodes: binary vector with the flags 1 if dead, 0 if alive.
%
% created with MATLAB ver.: 8.0.0.783 (R2012b) on Mac OS X  Version: 10.9.4
%
% created by: Jesus Fernandez-Bes (<a href="http://www.tsc.uc3m.es/~jesusfbes">web</a>)
% DATE: Oct-2014

switch param.mode
    
    case 'fixed' 
    % From iteration param.T a fixed percent of nodes are disconnected forever
        if(iter == param.T) 
            P = randperm(N);
            
            number_nodes = floor(param.percent * N);
            dead_nodes(P(1:number_nodes)) = 1;
        end
        
    case 'random'
    % At every iteration nodes are disconnected with probability param.p
        
        dead_nodes = rand(N, 1) <= param.P;
        
    case 'none'
    % do nothing
    
    otherwise 
        
        disp('error in dead node method')         
    

end

end