function dy = AMPA(t, y, f, ft, TimeControl)

% AMPA_RECEPTOR_KINETICS - 6-state model of glutamate AMPA receptors
%
% Description:
%   Implements a 6-state kinetic model of AMPA receptor gating based on:
%   - Patneau and Mayer, Neuron 6:785 (1991)
%   - Destexhe et al. (1996)
%
% States:
%           O
%           |
% C0 - C1 - C2
%      |    |
%      D1   D2
%
% Topology (from the equations below): C0<->C1<->C2; C1<->D1; C2<->D2; C2<->O.
% The open state O opens from the doubly-bound C2 (dy(6)=Ro*y(3)-Rc*y(6)),
% NOT from C1.
%
% State variables: y = [C0, C1, C2, D1, D2, O]
%
% Inputs:
%   t  - Current time (ms)
%   y  - State vector (6 elements, fractions must sum to 1)
%   f  - Glutamate concentration time series (mM)
%   ft - Time points for glutamate concentration (ms)
%   TimeControl - Integration horizon / time-control parameter (ms)
%
% Output:
%   dy - Derivatives of state variables (6x1 column vector)
%
% Example:
%   [t,y] = ode45(@(t,y) AMPA(t,y,f,ft,TimeControl), tspan, y0);

% Validate state vector
if length(y) ~= 6
    error('State vector y must have 6 elements corresponding to [C0, C1, C2, D1, D2, O]');
end

% Validate glutamate data
if length(f) ~= length(ft)
    error('Glutamate concentration f and time vector ft must have same length');
end
    % Interpolate glutamate concentration at current time
    glutamate_conc = interp1(ft, f, t, 'linear', 0);
    
    % Rate constants from the original MOD file
    % Binding and unbinding rates
    Rb  = 13;       % (/uM /ms): binding (diffusion limited)
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
    dy(3) = NormParam * (rb * y(2) - (Ru2 + Rd + Ro) * y(3) + Rr * y(5) + Rc * y(6));
    
    % dD1/dt: Single glutamate bound, desensitized
    dy(4) = NormParam * (Rd * y(2) - Rr * y(4));
    
    % dD2/dt: Double glutamate bound, desensitized
    dy(5) = NormParam * (Rd * y(3) - Rr * y(5));
    
    % dO/dt: Open state
    dy(6) = NormParam * (Ro * y(3) - Rc * y(6));
end