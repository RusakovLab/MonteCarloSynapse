function Unified_AMPA_SpaceSR(model_type)
% Unified version of AMPA_SpaceSR, AMPA1_SpaceSR, and AMPA2_SpaceSR
% Input: model_type - 'AMPA', 'AMPA1', or 'AMPA2'

% Common parameters
TimeGlu = 8; % ms 
BinGlu = 100;
timeCalculus = 200; %ms
DistanceReal = (1:1:BinGlu+1)*(1/BinGlu);
InitialDistribution = zeros(BinGlu, 101);

% Set model-specific parameters
switch model_type
    case 'AMPA'
        ode_function = @AMPA;
        initial_conditions = [1 0 0 0 0 0];
        desens_calc = @(Y) Y(:, 4) + Y(:, 5);
        open_calc = @(Y) Y(:, 6);
        desens_clim = [0, 1];
        open_clim = [0, 1];
    case 'AMPA1'
        ode_function = @AMPA1;
        initial_conditions = [1 0 0 0 0 0];
        desens_calc = @(Y) Y(:, 4);
         open_calc = @(Y) Y(:, 5)+Y(:, 6);
        desens_clim = [0, 0.4];
        open_clim = [0, 0.5];
    case 'AMPA2'
        ode_function = @AMPA2;
        initial_conditions = [1 0 0 0 0 0 0 0 0 0 0 0];
        desens_calc = @(Y) Y(:, 9) + Y(:, 10) + Y(:, 11) + Y(:, 12);
        open_calc = @(Y) Y(:, 6) + Y(:, 7) + Y(:, 8);
        desens_clim = [0, 1]; % Adjust if needed
        open_clim = [0, 1]; % Adjust if needed
    otherwise
        error('Invalid model type. Use ''AMPA'', ''AMPA1'', or ''AMPA2''');
end

% Load and process data files
myFolder = pwd;
if ~isfolder(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
    uiwait(warndlg(errorMessage));
    return;
end

filePattern = fullfile(myFolder, 'DistanceFree*.txt');
theFiles = dir(filePattern);
for K = 1:length(theFiles)
    baseFileName = theFiles(K).name;
    fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    thisStructure = load(fullFileName);
    output = thisStructure;
    InitialDistribution = InitialDistribution + output;
end
InitialDistribution = InitialDistribution/length(theFiles);

% Calculate glutamate concentration
Volume = ((4/3) * pi * ((2 * DistanceReal(2:end)').^3 - (2 * DistanceReal(1:end-1)').^3));
ConcentrationOneMolecules = 1.66 * 10^-6; % mM
Con = (1000)*ConcentrationOneMolecules; % uM
InitialDistribution(:, 2:end) = Con*InitialDistribution(:, 2:end) ./ Volume;

% Plot glutamate concentration
figure(1)
contourf(log(InitialDistribution(:, 2:end)), 10, 'LineColor', 'none')
colorbar
clim([-8,4.5]);
ax = gca;
ax.XTick = linspace(1, size(InitialDistribution(:, 2:end), 2), 5);
ax.XTickLabel = linspace(0, TimeGlu, 5);
ax.YTick = linspace(1, size(InitialDistribution(:, 2:end), 1), 5);
ax.YTickLabel = linspace(0, 2, 5);
xlabel('Time (ms)')
ylabel('Distance (um)')
title('Glu concentration')

% Differential AMPA Analysis
BinTotal = round(BinGlu * timeCalculus / TimeGlu);
AMPATimeBin = linspace(0, timeCalculus, BinTotal);
TwoBountRaw = 1:100;

% Preallocate arrays
GluTotal = zeros(1, BinTotal);
TwoBount = zeros(100, 100);
Open = zeros(100, 100);

% Set ODE options based on model type
if strcmp(model_type, 'AMPA2')
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-20, ...
                    'NonNegative', 1:12, 'Refine', 10, 'MaxStep', 0.1);
else
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-7*ones(1, length(initial_conditions)));
end

% Create figure for desensitization states
figure(4);
hold on;

% Main processing loop
for Jitter = 1:BinGlu
    fprintf('Processing iteration %d of %d\n', Jitter, BinGlu);
    TestGlu = InitialDistribution(Jitter, 2:end);
    
    if numel(TestGlu) < BinTotal
        TestGlu(BinTotal) = 0;
    else
        TestGlu = TestGlu(1:BinTotal);
    end
    
    current_ode = @(t, y) ode_function(t, y, TestGlu', AMPATimeBin, timeCalculus);
    [T, Y] = ode45(current_ode, [0 timeCalculus], initial_conditions, options);
    
    % Calculate and plot desensitization
    Desensitization = desens_calc(Y);
    plot(T, Desensitization);
    
    % Interpolate results
    Desensitization = interp1(T, Desensitization, TwoBountRaw);
    TwoBount(Jitter, :) = Desensitization;
    
    OpenState = interp1(T, open_calc(Y), TwoBountRaw);
    Open(Jitter, :) = OpenState;
end

hold off;
title('Desensitization States');
xlabel('Time');
ylabel('Desensitization');

% Plot state variables
figure(2);
num_states = length(initial_conditions);
for i = 1:num_states
    subplot(2, ceil(num_states/2), i);
    plot(T, Y(:,i));
    title(['State ' num2str(i)]);
end

% Plot desensitization state
figure(3)
contourf(TwoBount(:,1:end), 20, 'LineColor', 'none')
colorbar
clim(desens_clim);
ax = gca;
ax.XTick = linspace(1, size(TwoBount, 2), 5);
ax.XTickLabel = linspace(0, timeCalculus, 5);
ax.YTick = linspace(1, size(TwoBount, 1), 5);
ax.YTickLabel = linspace(0, 2, 5);
xlabel('Time (ms)')
ylabel('Distance (um)')
title([model_type ' desensitization state'])

% Plot open state
figure(6)
contourf(Open(:,1:end), 20, 'LineColor', 'none')
colorbar
clim(open_clim);
ax = gca;
ax.XTick = linspace(1, size(Open, 2), 5);
ax.XTickLabel = linspace(0, timeCalculus, 5);
ax.YTick = linspace(1, size(Open, 1), 5);
ax.YTickLabel = linspace(0, 2, 5);
xlabel('Time (ms)')
ylabel('Distance (um)')
title([model_type ' Open state'])

% Save data in ORIGIN format with model_type in filename
% Create time and distance vectors for headers
time_points = linspace(0, timeCalculus, size(TwoBount, 2));
distance_points = linspace(0, 2, size(TwoBount, 1));

% Save TwoBount data
desens_filename = sprintf('%s_Desensitization_State.dat', model_type);
fid = fopen(desens_filename, 'w');
% Write header line with time points
%fprintf(fid, 'Distance (um)\t');
%fprintf(fid, '%f\t', time_points(1:end-1));
%fprintf(fid, '%f\n', time_points(end));
% Write data
for i = 1:size(TwoBount, 1)
    %fprintf(fid, '%f\t', distance_points(i));
    fprintf(fid, '%f\t', TwoBount(i, 1:end-1));
    fprintf(fid, '%f\n', TwoBount(i, end));
end
fclose(fid);

% Save Open data
open_filename = sprintf('%s_Open_State.dat', model_type);
fid = fopen(open_filename, 'w');
% Write header line with time points
% fprintf(fid, 'Distance (um)\t');
% fprintf(fid, '%f\t', time_points(1:end-1));
% fprintf(fid, '%f\n', time_points(end));
% Write data
for i = 1:size(Open, 1)
    %fprintf(fid, '%f\t', distance_points(i));
    fprintf(fid, '%f\t', Open(i, 1:end-1));
    fprintf(fid, '%f\n', Open(i, end));
end
fclose(fid);

disp('Data saved in ORIGIN format:');
disp(['  - ' desens_filename]);
disp(['  - ' open_filename]);
end


