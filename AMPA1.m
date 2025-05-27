function dy = AMPA(t, y, f, ft, TimeControl)
% AMPA - Model of AMPA receptors kinetics
%
% 6-state gating model (Raman and Trussell, Neuron 9:173-186, 1992)
%      O1
%      |
% C -- C1 -- C2 -- O2
%      |
%      D
%
% Input:
%   t - time point
%   y - state variables [C0, C1, C2, D, O1, O2]
%   f - glutamate concentration time series
%   ft - time points for glutamate concentration
%   TimeControl - time control parameter
%
% Output:
%   dy - derivatives of state variables

    % Interpolate glutamate concentration at current time
    f = interp1(ft, f, t); % Interpolate the data set (ft, f) at times t
    
    % Rate constants from the MOD file
    Rb  = 13;       % (/mM /ms) : binding (diffusion limited)
    Ru1 = 0.3;      % (/ms)     : unbinding (1st site)
    Ru2 = 200;      % (/ms)     : unbinding (2nd site)      
    Rd  = 30.0;     % (/ms)     : desensitization
    Rr  = 0.02;     % (/ms)     : resensitization 
    Ro1 = 100;      % (/ms)     : opening (fast)
    Rc1 = 2;        % (/ms)     : closing
    Ro2 = 2;        % (/ms)     : opening (slow)
    Rc2 = 0.25;     % (/ms)     : closing
    
    % Convert mM to μM for consistency with the NMDA model if needed
    % (Rb is already in /mM/ms, so multiply by 1e-3 to convert from /μM/ms)
    % rb = Rb * (1e-3) * f;
    
    % Or keep in mM as per the original MOD file:
    rb = Rb * f;
    
    % Initialize derivatives vector
    dy = zeros(6,1);    % a column vector for [C0, C1, C2, D, O1, O2]
    
    % Normalization parameter (similar to NMDA model)
    NormParam = 1;
    
    % Calculate derivatives using kinetic scheme
    dy(1) = NormParam * (-rb * y(1) + Ru1 * y(2));                                  % dC0/dt
    dy(2) = NormParam * (rb * y(1) - (Ru1 + rb) * y(2) + Ru2 * y(3));               % dC1/dt
    dy(3) = NormParam * (rb * y(2) - (Ru2 + Rd + Ro1 + Ro2) * y(3) + ...
                          Rr * y(4) + Rc1 * y(5) + Rc2 * y(6));                      % dC2/dt
    dy(4) = NormParam * (Rd * y(3) - Rr * y(4));                                    % dD/dt
    dy(5) = NormParam * (Ro1 * y(3) - Rc1 * y(5));                                  % dO1/dt
    dy(6) = NormParam * (Ro2 * y(3) - Rc2 * y(6));                                  % dO2/dt
end