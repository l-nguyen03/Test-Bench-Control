InitTestBench

data = load('opt_controllers.mat');
Ks = data.opt_controller;

% Test Points: 
% (7500rpm, 5.08L/min) -> (8000 rpm, 5.69 L/min)
% (6000rpm, 7.12L/min) -> (6500 rpm, 5.88 L/min)
% (6000rpm, 0.5 L/min) -> (6500 rpm, 2.29 L/min)
% (6000rpm, 3.81 L/min) -> (7000 rpm, 4.37 L/min)
% (7000rpm, 0.5 L/min) -> (6500 rpm, 2.29 L/min)
% (8000rpm, 2.06 L/min) -> (7500 rpm, 3.15 L/min)
% (8000rpm, 9.32 L/min) -> (7500 rpm, 6.91 L/min)
% (9000rpm, 6.38 L/min) -> (8500 rpm, 7.46 L/min)
% (9000rpm, 8.43 L/min) ->  (8500 rpm, 7.46 L/min)
% (9000rpm, 4.33 L/min) -> (8500 rpm, 4.61 L/min)
% Test AP = [Time, PH_before, PH_after, Q_before, Q_after]
test_ap = [
    120, 67, 5.08;
    120, 32, 7.25;
];
title = {"Tuning Point", "Limitation Point"};
a = Ks(1).a;
b = Ks(1).b;
K_I = Ks(1).KP_KI;
K_P = K_I;

%Set up PT1
K = 1;
T = 0.5; 
sys = tf(K, [T, 1]);

t = 0:0.1:180;
steps = {[0.5, 4.35]; [0.5, -4.35]; [-0.5, 4.35]; [-0.5, -4.35]};

for i=1:size(test_ap, 1)
    for j=1:size(steps,1)
        step = steps{j};
        % % PT1 input attempt
        % open_system("TestController.slx");
        % PH_soll_step = test_ap(i, 2) + step(2) * (t >= test_ap(i, 1));
        % PH_soll_PT1_response = lsim(sys, PH_soll_step, t);
        % 
        % Q_soll_step = test_ap(i, 3) + step(1) * (t >= test_ap(i, 1));
        % Q_soll_PT1_response = lsim(sys, Q_soll_step, t);
        % 
        % PH_soll_PT1 = [t', PH_soll_PT1_response];
        % Q_soll_PT1 = [t', Q_soll_PT1_response];

        % % Run the simulation
        % simOut = sim('TestController.slx', 'SimulationMode', 'normal');

        % % Extract data from logsout
        % PH_meas_data = simOut.logsout.get('PH_meas').Values;
        % Q_meas_data = simOut.logsout.get('Q_meas').Values;
        % PH_soll_data = simOut.logsout.get('PH_soll_PT1').Values;
        % Q_soll_data = simOut.logsout.get('Q_soll_PT1').Values;

           
        % Step input attempt
        open_system("TestController2.slx");
        % Set parameters for the v_ref step block
        set_param('TestController2/PH_soll', 'Time', num2str(test_ap(i, 1)), ...
                                              'Before', num2str(test_ap(i, 2)), ...
                                              'After', num2str(test_ap(i,2)+step(2)))

        % Set parameters for the second input block
        set_param('TestController2/Q_soll', 'Time', num2str(test_ap(i, 1)), ...
                                              'Before', num2str(test_ap(i, 3)), ...
                                              'After', num2str(test_ap(i, 3)+step(1)))
        
        
        % Run the simulation
        simOut = sim('TestController2.slx', 'SimulationMode', 'normal');
        
        % Extract data from logsout
        PH_meas_data = simOut.logsout.get('PH_meas').Values;
        Q_meas_data = simOut.logsout.get('Q_meas').Values;
        PH_soll_data = simOut.logsout.get('PH_soll').Values;
        Q_soll_data = simOut.logsout.get('Q_soll').Values;
        
        v_ref = simOut.logsout.get('v_ref').Values;
        pos_ref = simOut.logsout.get('pos_ref').Values;
        v_meas = simOut.logsout.get('v_meas [rpm]').Values;
        pos_meas = simOut.logsout.get('pos_meas').Values;
    
        % Plotting the data
        figure; % Open a new figure window
    
        % Plot PH
        ax1 = subplot(4,1,1); % % Create a subplot for PH, change subplot config to 4,1,1 when also plotting speed and position
        timeIndex = PH_meas_data.Time >= 118 & PH_meas_data.Time < 180;
        plot(PH_meas_data.Time(timeIndex), PH_meas_data.Data(timeIndex), 'b', PH_soll_data.Time(timeIndex), PH_soll_data.Data(timeIndex), 'r--', "LineWidth", 1.5);
        ylabel('Pressure (mmHg)');
        legend('measured', 'setpoint');
    
        % Plot Q
        ax2 = subplot(4,1,2); % Create a subplot for Q, change subplot config to 4,1,2 when also plotting speed and position
        plot(Q_meas_data.Time(timeIndex), Q_meas_data.Data(timeIndex), 'b', Q_soll_data.Time(timeIndex), Q_soll_data.Data(timeIndex), 'r--', "LineWidth", 1.5);
        ylabel('Flowrate (L/min)');
        
        
        % Plot rotation speed
        ax3 = subplot(4,1,3);
        plot(v_ref.Time(timeIndex), v_ref.Data(timeIndex), 'b', v_meas.Time(timeIndex), v_meas.Data(timeIndex), 'r--', "LineWidth", 1.5);
        ylabel('Speed (rpm)');
        

        % Plot position
        ax4 = subplot(4,1,4);
        plot(pos_ref.Time(timeIndex), pos_ref.Data(timeIndex), 'b', pos_meas.Time(timeIndex), pos_meas.Data(timeIndex), "r--", "LineWidth", 1.5);
        ylabel('Position (mm)');
        xlabel('Time (s)');

        if (i<2) % Show results at tuning point
            xlim(ax1, [118, 131]);
            xlim(ax2, [118, 131]);
            xlim(ax3, [118, 131]);
            xlim(ax4, [118, 131]);
        else
            xlim(ax1, [118, 180]);
            xlim(ax2, [118, 180]);
            xlim(ax3, [118, 180]);
            xlim(ax4, [118, 180]);
        end


        linkaxes([ax1, ax2, ax3, ax4], 'x'); %add ax3, ax4 if plotting speed and position
     
        % Ensure plots do not overlap
        % sgtitle(["Step Response of controlled system at ", title{i}])
        sgtitle(["Step Response of controlled system at ", title{i}]);
        drawnow; % Update the figure window
    end
end