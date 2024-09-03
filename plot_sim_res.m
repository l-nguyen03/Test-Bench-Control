function [] = plot_sim_res(simulationResults)

    % Get unique values of v_ref
    uniqueVRef = unique(simulationResults.v_ref);

    % Initialize figures for the plots
    figure; % For all PH_meas vs pos_ref plots
    phPlot = gcf; % Get the current figure handle for PH_meas plot
    figure; % For all Q_meas vs pos_ref plots
    qPlot = gcf; % Get the current figure handle for Q_meas plot

    % Loop through each unique value of v_ref
    for i = 1:length(uniqueVRef)
        % Extract subset of data for current v_ref value
        currentData = simulationResults(simulationResults.v_ref == uniqueVRef(i), :);

        % Plot PH_meas vs pos_ref
        figure(phPlot);
        plot(currentData.pos_ref, currentData.PH_meas, 'o-', 'DisplayName', sprintf('v_{ref} = %.2f', uniqueVRef(i)));
        hold on;
    
        % Plot Q_meas vs pos_ref
        figure(qPlot);
        plot(currentData.pos_ref, currentData.Q_meas, 'o-', 'DisplayName', sprintf('v_{ref} = %.2f', uniqueVRef(i)));
        hold on;
    end

    % Add titles, labels, and legends
    figure(phPlot);
    title('PH_{meas} vs pos_{ref}');
    xlabel('pos_{ref}');
    ylabel('PH_{meas}');
    legend('show', 'Location', 'best');

    figure(qPlot);
    title('Q_{meas} vs pos_{ref}');
    xlabel('pos_{ref}');
    ylabel('Q_{meas}');
    legend('show', 'Location', 'best');
end