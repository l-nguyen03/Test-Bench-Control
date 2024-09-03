%Initialise Test Bench 
InitTestBench;

open_system('HemoTestBench.slx')
format short g

% Define a range of step function parameters [stepTime, clamp_pos, v]
arbeits_punkte = [];
%change granularity of v_ap as well as pos_ap to your liking
for v_ap = 6000:500:10000
    for pos_ap = 12:0.125:20
        arbeits_punkte = [arbeits_punkte; 120, pos_ap, v_ap];
    end
end

% Array to store outputs
PH_meas = [];
Q_meas = [];
PH_before_step = [];
Q_before_step = [];
iteration = 1; 
outputs = [];

%% Loop over each set of step function parameters and do 2 simulation:
% first simulation: v = const, pos jumps
% second simulation: v jumps, pos = const
for i = 1:size(arbeits_punkte, 1)
    %Define the absolute height of the jump
    delta_v = 250;
    delta_pos = 0.125;
    
    %First simulation
    % Set parameters for the v_ref step block
    set_param('HemoTestBench/v_ref', 'Time', num2str(arbeits_punkte(i, 1)), ...
                                          'Before', num2str(arbeits_punkte(i, 3)), ...
                                          'After', num2str(arbeits_punkte(i, 3)))

    % Set parameters for the pos_ref step block
    set_param('HemoTestBench/pos_ref', 'Time', num2str(arbeits_punkte(i, 1)), ...
                                          'Before', num2str(arbeits_punkte(i, 2)), ...
                                          'After', num2str(arbeits_punkte(i, 2) + delta_pos))

    % Run the simulation
    simOut = sim('HemoTestBench', 'SimulationMode', 'normal');

    % Find the time index right before jump at 120 seconds.
    PH_time = simOut.logsout.get("PH_meas").Values.Time;
    Q_time = simOut.logsout.get("Q_meas").Values.Time;
    [~, index_PH] = max(PH_time(PH_time < 120));
    [~, index_Q] = max(Q_time(Q_time < 120));

    % Retrieve the values at the nearest time point before 120 seconds.
    PH_before = simOut.logsout.get("PH_meas").Values.Data(index_PH);
    Q_before = simOut.logsout.get("Q_meas").Values.Data(index_Q);

    %Extract end values of outputs
    PH_meas = simOut.logsout.get("PH_meas").Values.Data(end);
    Q_meas = simOut.logsout.get("Q_meas").Values.Data(end);

    v_ref = arbeits_punkte(i, 3);
    pos_ref = arbeits_punkte(i, 2);

    delta_PH = PH_meas - PH_before;
    delta_Q = Q_meas - Q_before;
    outputs = [outputs; pos_ref, v_ref, PH_meas, Q_meas, PH_before, Q_before, 0, delta_pos, delta_PH, delta_Q];

    disp(['Simulation run ' num2str(iteration)])
    iteration = iteration + 1;

    % Second Simulation
    % Set parameters for the v_ref step block
    set_param('HemoTestBench/v_ref', 'Time', num2str(arbeits_punkte(i, 1)), ...
                                          'Before', num2str(arbeits_punkte(i, 3)), ...
                                          'After', num2str(arbeits_punkte(i, 3) + delta_v))

    % Set parameters for the second input block
    set_param('HemoTestBench/pos_ref', 'Time', num2str(arbeits_punkte(i, 1)), ...
                                          'Before', num2str(arbeits_punkte(i, 2)), ...
                                          'After', num2str(arbeits_punkte(i, 2)))

    % Run the simulation
    simOut = sim('HemoTestBench', 'SimulationMode', 'normal');

    % Find the index of the nearest time point before 120 seconds
    PH_time = simOut.logsout.get("PH_meas").Values.Time;
    Q_time = simOut.logsout.get("Q_meas").Values.Time;
    [~, index_PH] = max(PH_time(PH_time < 120));
    [~, index_Q] = max(Q_time(Q_time < 120));

    % Retrieve the values at the nearest time point before 120 seconds
    PH_before = simOut.logsout.get("PH_meas").Values.Data(index_PH);
    Q_before = simOut.logsout.get("Q_meas").Values.Data(index_Q);


    PH_meas = simOut.logsout.get("PH_meas").Values.Data(end);
    Q_meas = simOut.logsout.get("Q_meas").Values.Data(end);

    delta_PH = PH_meas - PH_before;
    delta_Q = Q_meas - Q_before;
    outputs = [outputs; pos_ref, v_ref, PH_meas, Q_meas, PH_before, Q_before, delta_v, 0, delta_PH, delta_Q];

    disp(['Simulation run ' num2str(iteration)])
    iteration = iteration + 1;

end

simulationResults = table(outputs(:,1), outputs(:,2), outputs(:,3) ,...
 outputs(:,4), outputs(:,5), outputs(:,6), outputs(:,7), outputs(:,8), outputs(:,9), outputs(:,10), 'VariableNames', {'pos_ref', 'v_ref', 'PH_meas', 'Q_meas', 'PH_before', 'Q_before', 'delta_v', 'delta_pos', 'delta_PH', 'delta_Q'});

% plot_sim_res(simulationResults)
