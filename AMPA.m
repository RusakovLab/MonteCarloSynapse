function dy = AMPA(t, y, f, ft, TimeControl)
% AMPA - Detailed model of glutamate AMPA receptors
%
% 6-state gating model:
%      O
%      |
% C0 - C1 - C2
%      |    |
%      D1   D2
%
% Based on the kinetic model by Patneau and Mayer, Neuron 6:785 (1991)
% and Destexhe et al. (1996)
%
% Input:
%   t - current time point
%   y - state variables [C0, C1, C2, D1, D2, O]
%   f - glutamate concentration time series
%   ft - time points for glutamate concentration
%   TimeControl - time control parameter (optional)
%
% Output:
%   dy - derivatives of state variables

    % Interpolate glutamate concentration at current time
    glutamate_conc = interp1(ft, f, t, 'linear', 0);
    
    % Rate constants from the original MOD file
    % Binding and unbinding rates
    Rb  = 13;       % (/mM /ms): binding (diffusion limited)
    Ru1 = 0.0059;   % (/ms)    : unbinding (1st site)
    Ru2 = 86;       % (/ms)    : unbinding (2nd site)
    
    % Desensitization and resensitization rates
    Rd  = 0.9;      % (/ms)    : desensitization
    Rr  = 0.064;    % (/ms)    : resensitization 
    
    % Opening and closing rates
    Ro  = 2.7;      % (/ms)    : opening
    Rc  = 0.2;      % (/ms)    : closing
    
    % Compute binding rate
    rb = Rb * glutamate_conc;
    
    % Initialize derivatives vector
    dy = zeros(6,1);
    
    % Normalization parameter
    NormParam = 1;
    
    % State transition derivatives
    % dC0/dt: Unbound state
    dy(1) = NormParam * (-rb * y(1) + Ru1 * y(2));
    
    % dC1/dt: Single glutamate bound
    dy(2) = NormParam * (rb * y(1) - (Ru1 + rb + Rd) * y(2) + Rr * y(4) + Ru2 * y(3));
    
    % dC2/dt: Double glutamate bound
    dy(3) = NormParam * (rb * y(2) - (Ru2 + Rd + Ro + Rc) * y(3) + Rr * y(5));
    
    % dD1/dt: Single glutamate bound, desensitized
    dy(4) = NormParam * (Rd * y(2) - Rr * y(4));
    
    % dD2/dt: Double glutamate bound, desensitized
    dy(5) = NormParam * (Rd * y(3) - Rr * y(5));
    
    % dO/dt: Open state
    dy(6) = NormParam * (Ro * y(3) - Rc * y(6));
end