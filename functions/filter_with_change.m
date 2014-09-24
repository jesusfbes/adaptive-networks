function [y] = filter_with_change(h1, h2, x)
%FILTER_WITH_CHANGE filter one input signal with two different systems.
%
%
% It computes the output y corresponding to filtering signal x with system h1
% for half of the simulation and h2 the rest. 
%   
%
% INPUT h1: System to filter in the first half of the simulation. 
%       h2: System to filter in the second half of the simulation. 
%       x: Signal to be filtered
%           
%
% OUTPUT y: output signal
%
%
% created by: Luis A. Azpicueta-Ruiz
% DATE: Oct-2014



%-----------------
% i. PREPARATIONS:
%-----------------


%determine additional information
N1 = length(h1); 
N2 = length(h2); 
Lx = length(x);   %length of signal

%init buffers, outputs & adaptive filter
xbuffer = zeros(N1, 1);
y = zeros(Lx, 1);

%-------------------------------------------------
% ii. TIME-DOMAIN FILTERING:
%-------------------------------------------------

for k = 1:1:Lx/2 
  
  %shift & refill generic buffer
  xbuffer = [x(k); xbuffer(1:N1-1)];
  
  %filter
  for n = 1:1:N1
    y(k) = y(k) + h1(n) * xbuffer(n);
  end
  
  
end

for k = (Lx/2+1):1:Lx 
  
  %shift & refill generic buffer
  xbuffer = [x(k); xbuffer(1:N1-1)];
  
  %filter
  for n = 1:1:N2
    y(k) = y(k) + h2(n) * xbuffer(n);
  end
  
  
end

y = y';