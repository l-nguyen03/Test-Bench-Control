%Start Simulation for different clamp position and rotation speed
simulation;

%Calculate Static Matrix
count = 1;
Ks = struct();
for i=1:2:size(simulationResults,1)

    y1 = [simulationResults{i, "delta_PH"}; simulationResults{i, "delta_Q"}];
    y2 = [simulationResults{i+1, "delta_PH"}; simulationResults{i+1, "delta_Q"}];
    u1 = [simulationResults{i, "delta_pos"}; simulationResults{i, "delta_v"}];
    u2 = [simulationResults{i+1, "delta_pos"}; simulationResults{i+1, "delta_v"}];

    y = [y1 y2];
    u = [u1 u2];
    u_inv = inv(u);


    Ks(count).Ks = y*u_inv;
    Ks(count).y = y;
    Ks(count).u = u;
    Ks(count).pos_ap = simulationResults{i, "pos_ref"};
    Ks(count).v_ap = simulationResults{i, "v_ref"};
    Ks(count).PH_ap = simulationResults{i, "PH_before"};
    Ks(count).Q_ap = simulationResults{i, "Q_before"};
    Ks(count).KP_KI = inv(Ks(count).Ks);
    Ks(count).delta_PH = (simulationResults{i, "delta_PH"} + simulationResults{i+1, "delta_PH"})/2;
    Ks(count).delta_Q = (simulationResults{i, "delta_Q"} + simulationResults{i+1, "delta_Q"})/2;
    count = count + 1;
end

% % Plot Ks and its vs different clamp position
% % How many clamp position did we simulated?
% batch_size = length(unique(simulationResults.pos_ref));
% plot_static(Ks, batch_size)
% plot_inv_ks(Ks, batch_size)



format long g
indices = [];
%Specify point of interest
aps = [
    7500, 5.08;
    6000, 0.5;
    6000, 7.12;
    9000, 4.33;
    9000, 8.43;
];

% Fetch nearest working point for each point of interest from the calculation above
for i = 1:size(aps, 1)
    v_ap = aps(i, 1); 
    Q_ap = aps(i, 2);

    % Calculate the Q_ap error for each instance in Ks
    errors = arrayfun(@(static_matrix) abs(static_matrix.Q_ap - Q_ap), Ks);

    v_matches = arrayfun(@(static_matrix) static_matrix.v_ap == v_ap, Ks);
    min_error = min(errors(v_matches));

    % Find indices of Ks instances that meet both criteria
    matching_indices = find(v_matches & errors == min_error);

    indices = [indices; matching_indices];
end

%Tune a and b for K_I = a*inv(K_s) and K_P = b*inv(K_s)
open_system('ControlledHemoTestBench.slx')
opt_controller = struct('Ks', [], 'y', [], 'u', [], 'pos_ap', [], 'v_ap', [], 'PH_ap', [], 'Q_ap', [], ...
    'KP_KI', [], 'delta_PH', [], 'delta_Q', [], 'a', [], 'b', [], 'ssTime', [], 'OS', [], 'OS_controll', []);
for i=1:length(indices)
    index = indices(i);
    K_P = Ks(index).KP_KI;
    K_I = K_P;

    %Initial Values for a and b
    a = 1;
    b = 1.5;
    
    %Start tuning a and b
    [settlingTime, OS, a, b, OS_controll] = find_opt(a, b, Ks(index).PH_ap, Ks(index).Q_ap, Ks(index).delta_PH, Ks(index).delta_Q)
    
    Ks(index).a = a;
    Ks(index).b = b;
    Ks(index).ssTime = settlingTime;
    Ks(index).OS = OS;
    Ks(index).OS_controll = OS_controll;
    opt_controller(i) = Ks(index);
end
% save('opt_controllers.mat', 'opt_controller')
