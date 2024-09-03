function [settlingTime, OS, a, b, OS_controll] = find_opt(a, b, PH_ref, Q_ref, delta_PH, delta_Q)
   
    overshootDetected = false;
    tolerancePercentage = 0.01;
    OS = zeros(2, 1);
    OS_controll = zeros(2,1);
    %Tune a
    while true
        % Set parameters for the v_ref step block
        set_param('ControlledHemoTestBench/PH_soll', 'Before', num2str(PH_ref), ...
                                              'After', num2str(PH_ref+delta_PH))

        % Set parameters for the second input block
        set_param('ControlledHemoTestBench/Q_soll', 'Before', num2str(Q_ref), ...
                                              'After', num2str(Q_ref+delta_Q))

        % Run the simulation
        simOut = sim('ControlledHemoTestBench', 'SimulationMode', 'normal');
        
        % Extract the logsout object
        logsout = simOut.get('logsout');
        
        % Extract the output signal
        PH_meas = logsout.get('PH_meas').Values.Data;
        Q_meas = logsout.get('Q_meas').Values.Data;
        PH_soll = logsout.get('PH_soll').Values.Data;
        Q_soll = logsout.get('Q_soll').Values.Data;

        pos_soll = logsout.get('pos_soll').Values.Data;
        v_soll = logsout.get('v_soll').Values.Data;
        pos_meas = logsout.get('pos_meas').Values.Data;
        v_meas = logsout.get('v_meas').Values.Data;

        time = simOut.tout;
        
        % Determine step direction and calculate overshoot and settling time for each output
        % Process each signal pair (PH and Q)
        for i = 1:2
            if i == 1
                outputSignal = PH_meas;
                refInput = PH_soll;
                controll_input = v_soll;
                controll_output = v_meas;
            else
                outputSignal = Q_meas;
                refInput = Q_soll;
                controll_input = pos_soll;
                controll_output = pos_soll;
            end

            % Find index where the step occurs in the reference signal
            stepIndex = find(diff(refInput) ~= 0, 1, 'first');
            timeOfStepChange = time(stepIndex);

            initialValue = refInput(stepIndex);
            steadyStateValue = refInput(end);
            
            initialValue_controll = controll_input(stepIndex-1);
            steadyStateValue_controll = controll_input(end);

            % Analyze the signal only after the step
            postStepSignal = outputSignal(stepIndex:end);
            postStepControll = controll_input(stepIndex:end);

            if steadyStateValue > initialValue
                % Positive step
                peakValue = max(postStepSignal);
                OS(i) = abs((peakValue - steadyStateValue) / steadyStateValue);
            else
                % Negative step
                minValue = min(postStepSignal);
                OS(i) = abs((steadyStateValue - minValue) / abs(steadyStateValue));
            end
            
            if steadyStateValue_controll > initialValue_controll
                % Positive step
                peakValue = max(postStepControll);
                OS_controll(i) = abs((peakValue - steadyStateValue_controll) / steadyStateValue_controll);
            else
                % Negative step
                minValue = min(postStepControll);
                OS_controll(i) = abs((steadyStateValue_controll - minValue) / abs(steadyStateValue_controll));
            end

        end
        if any(OS > 0.0008)|| any(OS_controll > 0.095)
            a = a - 0.1;
            assignin('base', 'a', a);
            break
        else
            a = a + 0.1;
            assignin('base', 'a', a);
        end

    end
    fprintf("Done tuning a: %d \n", a);
    fprintf("Start tuning b")
    %Tune b
    tmp = [];
    settlingTime = [];
    b = b + 0.1;
    assignin('base', 'b', b);
    while true
        
        % Run the simulation
        simOut = sim('ControlledHemoTestBench', 'SimulationMode', 'normal');
        
        % Extract the logsout object
        logsout = simOut.get('logsout');
        
        % Extract the output signal
        PH_meas = logsout.get('PH_meas').Values.Data;
        Q_meas = logsout.get('Q_meas').Values.Data;
        PH_soll = logsout.get('PH_soll').Values.Data;
        Q_soll = logsout.get('Q_soll').Values.Data;

        pos_soll = logsout.get('pos_soll').Values.Data;
        v_soll = logsout.get('v_soll').Values.Data;
        pos_meas = logsout.get('pos_meas').Values.Data;
        v_meas = logsout.get('v_meas').Values.Data;

        time = simOut.tout;
        
        % Determine step direction and calculate overshoot for each output
        % Process each signal pair (PH and Q)
        for i = 1:2
            if i == 1
                outputSignal = PH_meas;
                refInput = PH_soll;
                controll_input = v_soll;
                controll_output = v_meas;
            else
                outputSignal = Q_meas;
                refInput = Q_soll;
                controll_input = pos_soll;
                controll_output = pos_soll;
            end

            % Find index where the step occurs in the reference signal
            stepIndex = find(diff(refInput) ~= 0, 1, 'first');
            timeOfStepChange = time(stepIndex);

            initialValue = refInput(stepIndex);
            steadyStateValue = refInput(end);

            initialValue_controll = controll_input(stepIndex-1);
            steadyStateValue_controll = controll_input(end);


            % Analyze the signal only after the step
            postStepSignal = outputSignal(stepIndex:end);
            postStepControll = controll_input(stepIndex:end);

            OS_before = OS;
            OS_controll_before = OS_controll;

            if steadyStateValue > initialValue
                % Positive step
                peakValue = max(postStepSignal);
                OS(i) = abs((peakValue - steadyStateValue) / steadyStateValue);
            else
                % Negative step
                minValue = min(postStepSignal);
                OS(i) = abs((steadyStateValue - minValue) / steadyStateValue);
            end

            if steadyStateValue_controll > initialValue_controll
                % Positive step
                peakValue = max(postStepControll);
                OS_controll(i) = abs((peakValue - steadyStateValue_controll) / steadyStateValue_controll);
            else
                % Negative step
                minValue = min(postStepControll);
                OS_controll(i) = abs((steadyStateValue_controll - minValue) / abs(steadyStateValue_controll));
            end


            % Define upper and lower bounds for settling
            lowerBound = steadyStateValue * (1 - tolerancePercentage);
            upperBound = steadyStateValue * (1 + tolerancePercentage);

            % Find the index where the signal enters and stays in the tolerance band
            inBand = (postStepSignal >= lowerBound) & (postStepSignal <= upperBound);
            
            consecutiveInBand = strfind(inBand', ones(1, 50));

            if ~isempty(consecutiveInBand)
                firstSettlingIndex = consecutiveInBand(1) + 50 - 1;
                tmp(i) = time(stepIndex + firstSettlingIndex - 1) - timeOfStepChange;
            else
                tmp(i) = NaN;
            end
        end

        if any(OS > 0.001) || any(OS_controll > 0.1)
            % first run of b-scaling?
            if isempty(settlingTime)
                settlingTime = tmp;
            end
            b = b - 0.1;
            OS = OS_before;
            OS_controll = OS_controll_before
            assignin('base', 'b', b);
            break
        else
            b = b + 0.1;
            assignin('base', 'b', b);
            settlingTime = tmp;
        end
    end
    

    % Return the last steady-state time before overshoot was detected
end
