function [data_out] = load_single_dSpace_data(data_set,select)
    %Static map ==> select == 1    
    %
    
    %load relevant data:
    fname = fieldnames(data_set)';
    thisFieldName = fname{:};
    
    if select == 1
        experiment_interval = data_set.(thisFieldName).Y(1).Data;

        time = data_set.(thisFieldName).X(1).Data(experiment_interval == 1)';
        time = time - time(1);

        Pos_ref = data_set.(thisFieldName).Y(2).Data(experiment_interval == 1)';
        v_ref = data_set.(thisFieldName).Y(3).Data(experiment_interval == 1)';
        P_ao_meas = data_set.(thisFieldName).Y(6).Data(experiment_interval == 1)';
        P_lv_meas = data_set.(thisFieldName).Y(7).Data(experiment_interval == 1)'; 
        Pos_meas = data_set.(thisFieldName).Y(8).Data(experiment_interval == 1)';
        Q_meas = data_set.(thisFieldName).Y(9).Data(experiment_interval == 1)';
        %Temp = data_set.(thisFieldName).Y(10).Data(experiment_interval == 1)';
        v_meas = data_set.(thisFieldName).Y(13).Data(experiment_interval == 1)';

        %create data struct:
        data_out = struct('time',{time},'v_meas',{v_meas},'v_ref',{v_ref},'P_ao_meas',{P_ao_meas},...
            'P_lv_meas',{P_lv_meas},'Q_meas',{Q_meas},'Pos_ref',{Pos_ref},'Pos_meas',{Pos_meas});

    elseif select == 2
        experiment_interval = data_set.(thisFieldName).Y(3).Data;

        time = data_set.(thisFieldName).X(1).Data(experiment_interval == 1)';
        time = time - time(1);

        Pos_ref = data_set.(thisFieldName).Y(2).Data(experiment_interval == 1)';
        v_ref = data_set.(thisFieldName).Y(4).Data(experiment_interval == 1)';
        P_ao_meas = data_set.(thisFieldName).Y(8).Data(experiment_interval == 1)';
        P_lv_meas = data_set.(thisFieldName).Y(9).Data(experiment_interval == 1)'; 
        Pos_meas = data_set.(thisFieldName).Y(10).Data(experiment_interval == 1)';
        Q_meas = data_set.(thisFieldName).Y(11).Data(experiment_interval == 1)';
        %Temp = data_set.(thisFieldName).Y(12).Data(experiment_interval == 1)';
        v_meas = data_set.(thisFieldName).Y(15).Data(experiment_interval == 1)';

        %create data struct:
        data_out = struct('time',{time},'v_meas',{v_meas},'v_ref',{v_ref},'P_ao_meas',{P_ao_meas},...
            'P_lv_meas',{P_lv_meas},'Q_meas',{Q_meas},'Pos_ref',{Pos_ref},'Pos_meas',{Pos_meas});     
    end