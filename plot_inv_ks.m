function [] = plot_inv_ks(Ks, batchSize)
    n = numel(Ks);
    legendLabels = cell(1, n/batchSize); % Initialize cell array for legend labels

    % Loop over each batch
    figure; % Create a new figure
    for batch = 1:n/batchSize
        % Initialize arrays to store Ks values and pos_ref for each batch
        K_values = zeros(batchSize, 4); % 4 entries in Ks (2x2 matrix)
        pos_aps = zeros(batchSize, 1);
    
        % Extract data for this batch
        for i = 1:batchSize
            index = (batch - 1) * batchSize + i;
            K_values(i, :) = reshape(Ks(index).KP_KI, 1, []); % Reshape 2x2 matrix to 1x4
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
                    title("a11");
                    ylabel('a11');
                case 2
                    title("a21");
                    ylabel('a21');
                case 3
                    title('a12');
                    ylabel('a12');
                case 4
                    title('a22');
                    ylabel('a22');
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
