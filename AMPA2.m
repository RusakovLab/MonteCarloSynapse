function dy = AMPA2(t, y, f, ft, TimeControl)
% AMPA_KINETICS - Extended model of AMPA receptor gating (GluR1 parameters)
%
% This function models the kinetics of AMPA receptors based on the provided
% state diagram with multiple open, closed, and desensitized states.
% 
% Input:
%   t - current time point
%   y - state variables [R0, R1, R2, R3, R4, O2, O3, O4, D1, D2, D3, D4]
%   f - glutamate concentration time series
%   ft - time points for glutamate concentration
%   TimeControl - time control parameter (optional)
%
% Output:
%   dy - derivatives of state variables

    % Interpolate glutamate concentration at current time
    glutamate_conc = interp1(ft, f, t, 'linear', 0);
    
    % Rate constants for GluR1 (from table)
    alpha = 3100;      % (/s)
    beta = 8000;       % (/s)
    k_f = 2e7;         % (/M/s)
    k_r = 9000;        % (/s)
    gamma = 7.6; %3.3e-3;    % (/s)
    delta = 1800; % 1;         % (/s)
    
    CorrectopnCoef= 1000;
    % Rescale from s/M to ms/mM.
    %   First-order rates (/s):   divide by 1000  -> (/ms)
    %   k_f (/M/s):               divide by 1000^2 -> (/mM/ms)
    %   (i.e. /s->/ms is 1e-3, and /M->/mM is another 1e-3)
    alpha = alpha / CorrectopnCoef;         % (/s)  -> (/ms)
    beta = beta / CorrectopnCoef;           % (/s)  -> (/ms)
    k_f = k_f /CorrectopnCoef^2;            % (/M/s)-> (/mM/ms)
    k_r = k_r / CorrectopnCoef;             % (/s)  -> (/ms)

    gamma = gamma / CorrectopnCoef;         % (/s)  -> (/ms)
    delta = delta / CorrectopnCoef;         % (/s)  -> (/ms)
    
    % Compute binding rate
    rb = k_f * glutamate_conc;
    
    % Initialize derivatives vector
    dy = zeros(12,1);

    % Rate normalization parameter (as in AMPA and AMPA1)
    NormParam = 1;

    % State transition equations
    dy(1) = NormParam * (k_r*y(2) - 4*rb*y(1));  % R0
    dy(2) = NormParam * (4*rb*y(1) - (3*rb + k_r)*y(2)   + 2*k_r*y(3) -  delta*y(2) + gamma*y(9)); %R1
    dy(3) = NormParam * (3*rb*y(2) - (2*rb + 2*k_r)*y(3) + 3*k_r*y(4) -  delta*y(3) + gamma*y(10) + alpha*y(6) - 2*beta*y(3)); %R2
    dy(4) = NormParam * (2*rb*y(3) - (rb + 3*k_r)*y(4)   + 4*k_r*y(5) -  delta*y(4) + gamma*y(11) + alpha*y(7) - 3* beta * y(4)); %R3
    dy(5) = NormParam * (rb*y(4) - 4*k_r*y(5) - delta*y(5) +  gamma*y(12) + alpha*y(8) - 4*beta*y(5)); %R4

    dy(6) = NormParam * (2*beta*y(3) - alpha*y(6)); % O2
    dy(7) = NormParam * (3*beta*y(4) - alpha*y(7)); % O3
    dy(8) = NormParam * (4*beta*y(5) - alpha*y(8)); % O4

    dy(9) =  NormParam * (delta* y(2) -  gamma* y(9) - 3 * rb * y(9) + k_r * y(10)); % D1
    dy(10) =  NormParam * (delta* y(3) - gamma * y(10) - k_r * y(10) - 2 * rb * y(10) + 3*rb*y(9)  + 2 * k_r * y(11)); % D2
    dy(11) =  NormParam * (delta * y(4) - gamma * y(11) + 2 * rb * y(10) - 2 * k_r * y(11) - rb * y(11) + 3 * k_r * y(12)); % D3
    dy(12) =  NormParam * (delta* y(5) -  gamma* y(12) + rb * y(11) - 3 * k_r * y(12)); % D4

    % NOTE: The desensitized-row glutamate binding/unbinding (D1<->D2<->D3<->D4)
    % is already fully included in dy(9)..dy(12) above. A second
    % "Deep desensitization transitions" block that repeated exactly those terms
    % was removed here, because it doubled the D-row binding/unbinding rates.


    


end