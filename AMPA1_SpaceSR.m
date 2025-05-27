

TimeGlu = 8; % ms 
BinGlu=100;
timeCalculus = 200; %ms
DistanceReal=(1:1:BinGlu+1)*(1/BinGlu);
InitialDistribution = zeros(BinGlu, 101);

myFolder = pwd; % or 'C:\Users\yourUserName\Documents\My Pictures' or whatever...
% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isfolder(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, 'DistanceFree*.txt'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
for K = 1 : length(theFiles)
    baseFileName = theFiles(K).name;
    fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    % Now do whatever you want with this file name,
    % such as reading it in as an image array with imread()
    thisStructure = load(fullFileName);
    output =thisStructure;
    InitialDistribution=InitialDistribution+output;
end
InitialDistribution=InitialDistribution/length(theFiles);

%*************************Comput Glu concentration  from molecules number ***********************************
%NormParameter = 10^21/(6.023*10^23);
%InitialDistribution(:, 2:end) = NormParameter * InitialDistribution(:, 2:end) ./ Volume;
%  Volume is calculating the difference of the volumes of spheres with radii of 2 times the distances in the DistanceReal Vector.
Volume = ((4/3) * pi * ((2 * DistanceReal(2:end)').^3 - (2 * DistanceReal(1:end-1)').^3));
ConcentrationOneMolecules = 1.66 * 10^-6; % mM 1 molecules creat such concentration in mM 
Con= (1000)*ConcentrationOneMolecules; % uM
InitialDistribution(:, 2:end) = Con*InitialDistribution(:, 2:end) ./ Volume;


figure(1)
contourf(log(InitialDistribution(:, 2:end)),10, 'LineColor','none')
colorbar
clim([-8,4.5]);
% Set correct axis scaling
ax = gca; % Get current axis
ax.XTick = linspace(1, size(InitialDistribution(:, 2:end), 2), 5); % 5 tick marks along X
ax.XTickLabel = linspace(0, TimeGlu, 5); % Label as 0 to 200 ms
ax.YTick = linspace(1, size(InitialDistribution(:, 2:end), 1), 5); % 5 tick marks along Y
ax.YTickLabel = linspace(0, 2, 5); % Label as 0 to 2 mkm
xlabel('Time (ms)')
ylabel('Distance (um)')
title('Glu concentration')
%********************* Differential AMPA Analysis ***********************************
% Calculate total number of bins based on provided parameters
BinTotal = round(BinGlu * timeCalculus / TimeGlu);

% Generate time bins for interpolation (done once outside loop)
AMPATimeBin = linspace(0, timeCalculus, BinTotal);
TwoBountRaw = 1:100; % Define output time points once

% Preallocate all arrays for better performance
GluTotal = zeros(1, BinTotal);
TwoBount = zeros(100, 100);
Open = zeros(100, 100); % Preallocation was missing

% Set ODE options outside the loop (unchanged during iterations)
options = odeset('RelTol', 1e-8, 'AbsTol', [1e-7 1e-7 1e-7 1e-7 1e-7 1e-7]);
ode = @(t, y) AMPA1(t, y, zeros(BinTotal, 1), AMPATimeBin, timeCalculus); % Template function

% Create figure once, outside the loop
figure(4);
hold on;

% Main processing loop
for Jitter = 1:BinGlu
    fprintf('Processing iteration %d of %d\n', Jitter, BinGlu);
    % Extract glutamate distribution for this iteration
    TestGlu = InitialDistribution(Jitter, 2:end);
    
    % Pad with zeros if needed
    if numel(TestGlu) < BinTotal
        TestGlu(BinTotal) = 0;
    else
        TestGlu = TestGlu(1:BinTotal);
    end
    
    % Update ODE function with current TestGlu
    ode = @(t, y) AMPA1(t, y, TestGlu', AMPATimeBin, timeCalculus);
    
    % Solve ODE system
    [T, Y] = ode45(ode, [0 timeCalculus], [1 0 0 0 0 0], options);
    
    % Extract and process desensitization state
    Desensitization = Y(:, 4);
    plot(T, Desensitization); % Plot to figure created outside loop
    
    % Interpolate to standard time points (once for each state)
    Desensitization = interp1(T, Desensitization, TwoBountRaw);
    TwoBount(Jitter, :) = Desensitization;
    
    % Process open state
    OpenState = interp1(T, Y(:, 6), TwoBountRaw);
    Open(Jitter, :) = OpenState;
    
end

% Release figure hold state
hold off;
title('Desensitization States');
xlabel('Time');
ylabel('Desensitization');

axis manual
figure(2);                  % display image
%axis([0 timeCalculus 0 1])

subplot(2, 4, 1); plot(T,Y(:,1))
%axis([0 timeCalculus 0 1])

subplot(2, 4, 2); plot(T,Y(:,2))
%axis([0 timeCalculus 0 1])

subplot(2, 4, 3); plot(T,Y(:,3))
%axis([0 timeCalculus 0 1])

subplot(2, 4, 4); plot(T,Y(:,4))
%axis([0 timeCalculus 0 0.02])

subplot(2, 4, 5); plot(T,Y(:,5))
%axis([0 timeCalculus 0 0.01])

subplot(2, 4, 6); plot(T,Y(:,6))
%axis([0 timeCalculus 0 0.01])


% Figure 3: TwoBount visualization
figure(3)
contourf(TwoBount(:,1:end), 20, 'LineColor', 'none')
colorbar
clim([0, 0.4])

% Set correct axis scaling
ax = gca; % Get current axis
ax.XTick = linspace(1, size(TwoBount, 2), 5); % 5 tick marks along X
ax.XTickLabel = linspace(0, timeCalculus, 5); % Label as 0 to 200 ms
ax.YTick = linspace(1, size(TwoBount, 1), 5); % 5 tick marks along Y
ax.YTickLabel = linspace(0, 2, 5); % Label as 0 to 2 mkm
xlabel('Time (ms)')
ylabel('Distance (um)')
title('AMPA desentization state')

% Figure 6: Open visualization
figure(6)
contourf(Open(:,1:end), 20, 'LineColor', 'none')
colorbar
clim([0, 0.03])

% Set correct axis scaling
ax = gca; % Get current axis
ax.XTick = linspace(1, size(Open, 2), 5); % 5 tick marks along X
ax.XTickLabel = linspace(0, timeCalculus, 5); % Label as 0 to 200 ms
ax.YTick = linspace(1, size(Open, 1), 5); % 5 tick marks along Y
ax.YTickLabel = linspace(0, 2, 5); % Label as 0 to 2 mkm
xlabel('Time (ms)')
ylabel('Distance (um)')
title('AMPA Open state')




