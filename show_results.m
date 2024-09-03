InitTestBench;

lutFig = openfig('static_pump_LUT.fig','reuse');

axLUT = get(lutFig, 'CurrentAxes');
linesLUT = get(axLUT, 'Children');

%Extract LUT 
% Initialize an array to hold the original data
originalLUTData = [];

% Loop through each line object and check for a distinguishing feature of the 'Original' lines
for i = 1:length(linesLUT)
    lineObj = linesLUT(i);

    if strcmp(lineObj.LineStyle, '-') 
        xData = get(lineObj, 'XData');
        yData = get(lineObj, 'YData');

        % Extract pump speed from the DisplayName
        displayName = get(lineObj, 'DisplayName');
        pumpSpeedStr = regexp(displayName, '=\s*(\d+)', 'tokens');
        if ~isempty(pumpSpeedStr)
            pumpSpeed = str2double(pumpSpeedStr{1});
        else
            pumpSpeed = NaN; 
        end

        % Concatenate the data to the originalLUTData array
        % Each row contains pump speed, flowrate, and pressure
        originalLUTData = [originalLUTData; repmat(pumpSpeed, length(xData), 1), xData', yData'];
    end
end
close(lutFig);

originalLUTData = sortrows(originalLUTData, 1); 

% Load vertices of the operating area
data = load("op_bereich_vertices.mat");
vertices = data.vertices;

minX = min(vertices(:, 1));
maxX = max(vertices(:, 1));
minY = min(vertices(:, 2));
maxY = max(vertices(:, 2));

% Generate random points inside the operating area
N = 2000;  % Number of random points
randX = minX+0.6 + (maxX - minX + 0.6) * rand(N, 1);
randY = minY + (maxY - minY) * rand(N, 1);

% Check if points are inside the polygon
inside = inpolygon(randX, randY, vertices(:,1), vertices(:,2));
randX = randX(inside);
randY = randY(inside);

% Plotting initial operating range
figure(1);
fill([vertices(:,1); vertices(1,1)],[vertices(:,2); vertices(1,2)], [0.2196, 0.4627, 0.749]);
hold on
plot(5.08, 67, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0.9529, 0.9412, 0.7922], 'MarkerFaceColor', [0.9529, 0.9412, 0.7922]);

hold off

xlabel('Flowrate (L/min)');
ylabel('Pressure (mmHg)');
title('Operating Range of Sputnik Pump');
legend("Borderline", "Operating Point for tuning controller")
grid on;


%% Load controller and run tests for 4 different step inputs.
data = load('opt_controllers.mat');
Ks = data.opt_controller;
open_system("TestController2.slx");

% Choose the controller you want to test.
a = Ks(1).a;
b = Ks(1).b;
K_I = Ks(1).KP_KI;
K_P = K_I;

steps = {[0.5, 4.35]; [0.5, -4.35]; [-0.5, 4.35]; [-0.5, -4.35]};

performance_tot = {};

for j=1:size(steps, 1)
    step = steps{j};
    OS_pressure = [];
    OS_flowrate = [];
    ssTime_pressure = [];
    ssTime_flowrate = [];
    performance = struct();
    for i=1:length(randX)
        OS_tmp = zeros(2, 1);
        ssTime_tmp = zeros(2,1);

        % Set parameters for the v_ref step block
        set_param('TestController2/PH_soll', 'Time', '120', ...
                                              'Before', num2str(randY(i)), ...
                                              'After', num2str(randY(i)+step(2)))

        % Set parameters for the second input block
        set_param('TestController2/Q_soll', 'Time', '120', ...
                                              'Before', num2str(randX(i)), ...
                                              'After', num2str(randX(i)+step(1)))


        % Run the simulation
        simOut = sim('TestController2.slx', 'SimulationMode', 'normal');

        % Extract data from logsout
        PH_meas = simOut.logsout.get('PH_meas').Values.Data;
        Q_meas = simOut.logsout.get('Q_meas').Values.Data;
        PH_soll = simOut.logsout.get('PH_soll').Values.Data;
        Q_soll = simOut.logsout.get('Q_soll').Values.Data;
        time = simOut.tout;

        pressure = randY(i);
        flow_rate = randX(i);

        for k = 1:2
            if k == 1
                outputSignal = PH_meas;
                refInput = PH_soll;
            else
                outputSignal = Q_meas;
                refInput = Q_soll;
            end

            % Find index where the step occurs in the reference signal
            stepIndex = find(diff(refInput) ~= 0, 1, 'first');
            timeOfStepChange = time(stepIndex);

            initialValue = refInput(stepIndex);
            steadyStateValue = refInput(end);

            % Analyze the signal only after the step
            postStepSignal = outputSignal(stepIndex:end);

            if steadyStateValue > initialValue
                % Positive step
                peakValue = max(postStepSignal);
                OS_tmp(k) = abs((peakValue - steadyStateValue) / steadyStateValue);
            else
                % Negative step
                minValue = min(postStepSignal);
                OS_tmp(k) = abs((steadyStateValue - minValue) / steadyStateValue);
            end
            % Define upper and lower bounds for settling
            lowerBound = steadyStateValue * (1 - 0.01);
            upperBound = steadyStateValue * (1 + 0.01);

            % Find the index where the signal enters and stays in the tolerance band
            inBand = (postStepSignal >= lowerBound) & (postStepSignal <= upperBound);

            consecutiveInBand = strfind(inBand', ones(1, 50));

            if ~isempty(consecutiveInBand)
                firstSettlingIndex = consecutiveInBand(1) + 50 - 1;
                ssTime_tmp(k) = time(stepIndex + firstSettlingIndex - 1) - timeOfStepChange;
            else
                ssTime_tmp(k) = NaN;
            end
        end
        performance(i).OS_pressure = OS_tmp(1);
        OS_pressure = [OS_pressure; OS_tmp(1)];

        performance(i).OS_flowrate = OS_tmp(2);
        OS_flowrate = [OS_flowrate; OS_tmp(2)];

        performance(i).ssTime_pressure = ssTime_tmp(1);
        ssTime_pressure = [ssTime_pressure; ssTime_tmp(1)];

        performance(i).ssTime_flowrate = ssTime_tmp(2);
        ssTime_flowrate = [ssTime_flowrate; ssTime_tmp(2)];
        
        performance(i).pressure_ap = pressure;
        performance(i).flowrate_ap = flow_rate;
    end
    performance_tot{j} = {performance, step};

    figure;
    subplot(1,2,1);
    scatter(randX, randY, [], OS_pressure*100, 'filled');
    hold on;
    plot([vertices(:,1); vertices(1,1)],[vertices(:,2); vertices(1,2)], '-', "Color", [0.098, 0.149, 0.333]);
    plot([0.5, 0.5], [20, 120], "r-");
    hold off;
    colormap('winter'); % Apply colormap
    colorbar; % Add colorbar
    xlabel('Flowrate (L/min)');
    ylabel('Pressure (mmHg)');
    legend('Operating Point', 'Operation Limit', 'flowrate=0.5L/min');
    title(['OS Pressure (%), step: (', num2str(step(1)), 'L/min,', num2str(step(2)), 'mmHg)']);
    
    subplot(1,2,2);
    scatter(randX, randY, [], OS_flowrate*100, 'filled');
    hold on;
    plot([vertices(:,1); vertices(1,1)],[vertices(:,2); vertices(1,2)], '-', "Color", [0.098, 0.149, 0.333]);
    plot([0.5, 0.5], [20, 120], "r-");
    hold off;
    colormap('winter'); % Apply colormap
    colorbar; % Add colorbar
    xlabel('Flowrate (L/min)');
    ylabel('Pressure (mmHg)');
    legend('Operating Point', 'Operation Limit', 'flowrate=0.5L/min');
    title(['OS Flowrate (%), step: (', num2str(step(1)), 'L/min,', num2str(step(2)), 'mmHg)']);
    
    figure;
    subplot(1,2,1);
    scatter(randX, randY, [], ssTime_pressure, 'filled');
    hold on;
    plot([vertices(:,1); vertices(1,1)],[vertices(:,2); vertices(1,2)], '-', "Color", [0.098, 0.149, 0.333]);
    plot([0.5, 0.5], [20, 120], "r-");
    hold off;
    colormap('winter'); % Apply colormap
    colorbar; % Add colorbar
    xlabel('Flowrate (L/min)');
    ylabel('Pressure (mmHg)');
    legend('Operating Point', 'Operation Limit', 'flowrate=0.5L/min');
    title(['Settling Time (s) Pressure, step: (', num2str(step(1)), 'L/min,', num2str(step(2)), 'mmHg)']);
    
    subplot(1,2,2);
    scatter(randX, randY, [], ssTime_flowrate, 'filled');
    hold on;
    plot([vertices(:,1); vertices(1,1)],[vertices(:,2); vertices(1,2)], '-', "Color", [0.098, 0.149, 0.333]);
    plot([0.5, 0.5], [20, 120], "r-");
    hold off;
    colormap('winter'); % Apply colormap
    colorbar; % Add colorbar
    xlabel('Flowrate (L/min)');
    ylabel('Pressure (mmHg)');
    legend('Operating Point', 'Operation Limit', 'flowrate=0.5L/min');
    title(['Settling Time (s) Flowrate, step: (', num2str(step(1)), 'L/min,', num2str(step(2)), 'mmHg)']);

end

numPoints = length(randX); 
numJ = length(performance_tot);

% Initialize sum variables
sumOSPressure = zeros(numPoints, 1);
sumOSFlowrate = zeros(numPoints, 1);
sumSSTimePressure = zeros(numPoints, 1);
sumSSTimePressure = zeros(numPoints, 1);

% Loop over all j
for j = 1:numJ
    performance = performance_tot{j}{1};
    % Accumulate sums for each i
    for i = 1:numPoints
        sumOSPressure(i) = sumOSPressure(i) + performance(i).OS_pressure;
        sumOSFlowrate(i) = sumOSFlowrate(i) + performance(i).OS_flowrate;
        sumSSTimePressure(i) = sumSSTimePressure(i) + performance(i).ssTime_pressure;
        sumSSTimeFlowrate(i) = sumSSTimePressure(i) + performance(i).ssTime_flowrate;
    end
end

% Calculate averages
avgOSPressure = sumOSPressure / numJ;
avgOSFlowrate = sumOSFlowrate / numJ;
avgSSTimePressure = sumSSTimePressure / numJ;
avgSSTimeFlowrate = sumSSTimeFlowrate / numJ;

%plot averages
figure;
subplot(1,2,1);
scatter(randX, randY, [], avgOSPressure*100, 'filled');

hold on;
plot([vertices(:,1); vertices(1,1)],[vertices(:,2); vertices(1,2)], '-', "Color", [0.098, 0.149, 0.333]);
plot([0.5, 0.5], [20, 120], "r-");
hold off;
colormap('winter'); % Apply colormap
colorbar; % Add colorbar
xlabel('Flowrate (L/min)');
ylabel('Pressure (mmHg)');
legend('Operating Point', 'Operation Limit', 'flowrate=0.5L/min');
title('Average OS Pressure (%) Heatmap');

subplot(1,2,2);
scatter(randX, randY, [], avgOSFlowrate*100, 'filled');
hold on;
plot([vertices(:,1); vertices(1,1)],[vertices(:,2); vertices(1,2)], '-', "Color", [0.098, 0.149, 0.333]);
plot([0.5, 0.5], [20, 120], "r-");
hold off;
colormap('winter'); % Apply colormap
colorbar; % Add colorbar
xlabel('Flowrate (L/min)');
ylabel('Pressure (mmHg)');
legend('Operating Point', 'Operation Limit', 'flowrate=0.5L/min');
title('Average OS Flowrate (%)');

figure;
subplot(1,2,1);
scatter(randX, randY, [], avgSSTimePressure, 'filled');
hold on;
plot([vertices(:,1); vertices(1,1)],[vertices(:,2); vertices(1,2)], '-', "Color", [0.098, 0.149, 0.333]);
plot([0.5, 0.5], [20, 120], "r-");
hold off;
colormap('winter'); % Apply colormap
colorbar; % Add colorbar
xlabel('Flowrate (L/min)');
ylabel('Pressure (mmHg)');
legend('Operating Point', 'Operation Limit', 'flowrate=0.5L/min');
title('Average Settling Time (s) Pressure');

subplot(1,2,2);
scatter(randX, randY, [], avgSSTimeFlowrate, 'filled');
hold on;
plot([vertices(:,1); vertices(1,1)],[vertices(:,2); vertices(1,2)], '-', "Color", [0.098, 0.149, 0.333]);
plot([0.5, 0.5], [20, 120], "r-");
hold off;
colormap('winter'); % Apply colormap
colorbar; % Add colorbar
xlabel('Flowrate (L/min)');
ylabel('Pressure (mmHg)');
legend('Operating Point', 'Operation Limit', 'flowrate=0.5L/min');
title('Average Settling Time (s) Flowrate');



