function [] = plot_static(Ks, batchSize)
    n = numel(Ks);
    legendLabels = cell(1, n/batchSize); % Initialize cell array for legend labels

    % Loop over each batch (each rotation speed)
    figure; % Create a new figure
    for batch = 1:n/batchSize
        % Initialize arrays to store Ks values and pos_ref for each batch
        K_values = zeros(batchSize, 4); % 4 entries in Ks (2x2 matrix)
        pos_aps = zeros(batchSize, 1);
    
        % Extract data for this batch
        for i = 1:batchSize
            index = (batch - 1) * batchSize + i;
            K_values(i, :) = reshape(Ks(index).Ks, 1, []); % Reshape 2x2 matrix to 1x4 [a11, a12, a21, a22]
            pos_aps(i) = Ks(index).pos_ap;
        end
        
        legendLabels{batch} = sprintf('v_{ap} = %.2f', Ks(index).v_ap);

        % Plot each entry of Ks vs pos_ap
        for j = 1:4
            subplot(2, 2, j); % Create a subplot for each entry of Ks
            plot(pos_aps, K_values(:, j), 'o-', 'DisplayName', legendLabels{batch}); % Plot the j-th entry of Ks vs pos_ref
            xlabel('pos_{ap}');
            switch j
                case 1
                    title('a_{11}: Delta PH response to 1mm step in pos at each (pos,v)');
                    ylabel('Delta PH');
                case 2
                    title('a_{21}: Delta Q response to 1mm step in pos at each (pos,v)');
                    ylabel('Delta Q');
                case 3
                    title('a_{12}: Delta PH response to 200 rpm step in v at each (pos,v)');
                    ylabel('(Delta PH)/200');
                case 4
                    title('a_{22}: Delta Q response to 200 rpm step in v at each (pos,v)');
                    ylabel('(Delta Q)/200');
            end

            hold on; % Hold on for overlapping plots from different batches
        end
    end
    % Add legend to each subplot
    for j = 1:4
        subplot(2, 2, j);
        legend('show', 'Location', 'best');
    end
end
