clear all;
close all;
clc;

%% load data and plot raw signals:
%data_raw = load('..\Static Characterization MCL\Data\Exp_static_02_mVisc.mat');
data_raw = load('./2023_10_23_GlycerinWater/Static_Map_Glycerin.mat');
data = load_single_dSpace_data(data_raw,2);

% figure(1)
% 
% subplot(3,1,1);
% plot(data.time,data.Q_meas);
% hold on;
% axis tight
% title("Flowrate")
% xlabel("Time")
% ylabel("Flowrate (L/min)")

% subplot(3,1,2)
% plot(data.time,data.P_ao_meas);
% hold on;
% plot(data.time,data.P_lv_meas);
% hold on;
% axis tight
% title("Pressure")
% xlabel("Time")
% ylabel("Pressure (mmHg)")
% legend("Measured Aortic Pressure", "Measured Left Ventricular Pressure")
% 
% subplot(3,1,3);
% plot(data.time,data.Pos_meas);
% hold on;
% plot(data.time,data.Pos_ref);
% hold on;
% axis tight
% title("Clamp position")
% xlabel("Time")
% ylabel("Position (mm)")
% legend("Measured Position", "Reference Position")

%% Determine Static Pump Lookup-Table:
Pos_ref_vals = unique(data.Pos_ref);
v_ref_vals = unique(data.v_ref);

LUT_Clamp_PH = zeros(length(Pos_ref_vals),length(v_ref_vals));
LUT_P_ao = zeros(length(Pos_ref_vals),length(v_ref_vals));
LUT_Clamp_Q = zeros(length(Pos_ref_vals),length(v_ref_vals));
for i = 1:1:length(Pos_ref_vals)
    for j = 1:1:length(v_ref_vals)
        indices = find(data.Pos_ref == Pos_ref_vals(i) & data.v_ref == v_ref_vals(j));
        Q_cut1 = data.Q_meas(indices);
        Q_cut2 = Q_cut1(2000:end);
        LUT_Clamp_Q(i,j) = mean(Q_cut2);
        PH_cut1 = data.P_ao_meas(indices)-data.P_lv_meas(indices);
        PH_cut2 = PH_cut1(2000:end);
        LUT_Clamp_PH(i,j) = mean(PH_cut2);
        P_ao_cut = data.P_ao_meas(indices);
        LUT_P_ao(i,j) = mean(P_ao_cut);
    end
end

indices = find(LUT_P_ao>max(LUT_P_ao)-1);
LUT_P_ao_deleted = LUT_P_ao;
LUT_P_ao_deleted(indices) = NaN;
LUT_Clamp_PH(indices) = NaN;

% figure(2)
% plot(LUT_Clamp_Q,LUT_Clamp_PH);
% title("Static Pump Lookup-Table")
% xlabel("Flowrate (L/min)")
% ylabel("Pressure Head (mmHg)")


%% Determine Lookup-Table from Speed (v_ref) and Flow (Q_meas) to Pressure Head (PH):
Q_vals = [0:0.5:13]; %desired Flow values for Lookup-Table
% legendLabels = cell(1, length(v_ref_vals)*2); % Initialize cell array for legend
LUT_Static_Pump = zeros(length(v_ref_vals),length(Q_vals)); %initialize Lookup-Table with zeros
for i = 1:1:length(v_ref_vals)

    %
    x = LUT_Clamp_Q(:,i);
    y1 = LUT_Clamp_PH(:,i);           
    
    %Delete NaN entries:
    NaN_indices = isnan(y1);
    y1 = y1(~NaN_indices);
    x = x(~NaN_indices);
    
    %Fit linear model:
    P = polyfit(x,y1,1);
    yfit = P(1)*x+P(2);
    %hold on;

    %Plot fitted signal:
    % figure(3)
    % plot(x,yfit,'r-.');hold on;
    % plot(x,y1);
    % legendLabels{i*2-1} = sprintf('Fitted v_{ref}=%.2f', v_ref_vals(i));
    % legendLabels{i*2} = sprintf('Original v_{ref}=%.2f', v_ref_vals(i));
    %Save values to Lookup-Table:
    LUT_Static_Pump(i,:) = P(1)*Q_vals+P(2);
    
end
% title("Static Pump Lookup-Table")
% xlabel("Flowrate (L/min)")
% ylabel("Pressure Head (mmHg)")
% legend(legendLabels, 'Location', 'best');

%% Determine turbulent tube (with bag) resistance:
PH_clamp_open = LUT_Clamp_PH(1,:);
Q_clamp_open = LUT_Clamp_Q(1,:);
exponents = [1.4:0.001:1.8]; %Exponents under investigation
R_means = zeros(1,length(exponents));
R_variances = zeros(1,length(exponents));
for i=1:1:length(exponents)
    R_clamp_open_vals = PH_clamp_open./Q_clamp_open.^exponents(i);
    R_variances(i) = var(R_clamp_open_vals);
    R_means(i) = mean(R_clamp_open_vals);
end
[value_min,indice_min] = min(R_variances);
exponent = exponents(indice_min);
R_clamp_open = R_means(indice_min);

% figure(4)
% subplot(3,1,1);
% plot(v_ref_vals,Q_clamp_open);hold on;
% xlabel("Reference Velocity")
% ylabel("Flowrate, clamp open")
% subplot(3,1,2);
% plot(v_ref_vals,PH_clamp_open);
% xlabel("Reference Velocity")
% ylabel("Pressure Head, clamp open")
% subplot(3,1,3);
% plot(v_ref_vals,PH_clamp_open./Q_clamp_open.^exponent);
% xlabel("Reference Velocity")
% ylabel("Resistance")

LUT_Clamp_Q(find(LUT_Clamp_Q < 0.05)) = NaN;
PH_tube = R_clamp_open.*(LUT_Clamp_Q.^exponent);


%% Determine clamp resistance:
% figure(5)
% subplot(3,1,1);
% plot(Pos_ref_vals,LUT_Clamp_Q);
% subplot(3,1,2);
% plot(Pos_ref_vals,LUT_Clamp_PH);
% subplot(3,1,3);
% plot(Pos_ref_vals,(LUT_Clamp_PH-LUT_Clamp_PH)./LUT_Clamp_Q);hold on;

PH_clamp = LUT_Clamp_PH-PH_tube;
Q_clamp = LUT_Clamp_Q;
R_clamp = PH_clamp./Q_clamp;
R_clamp_mean = mean(R_clamp,2,"omitnan");

indices = find(isnan(R_clamp_mean));
Pos_step_size = diff(Pos_ref_vals(indices(1)-2:indices(1)-1));
R_clamp_step = diff(R_clamp_mean(indices(1)-2:indices(1)-1));
for i = 1:1:length(indices)
    R_clamp_mean(indices(i)) = R_clamp_mean(indices(i)-1) + R_clamp_step;
end

% figure(6)
% plot(Pos_ref_vals,R_clamp);hold on;
% plot(Pos_ref_vals,R_clamp_mean,'--');hold on;
% xlabel("Reference Position")
% ylabel("Clamp Resistance")


