%% return Model Builder
%% version 1: change the way of defining input variables.
%% version 2: can used for CV; change output structure from cell to struct; use parameter ModelType
%%IndexCV: index for testing data; 1 for testing data
%%2016-08-29 change name from CalibrateGC_FitModel_function2_2 to Analysis_FitModel_function2
function ModelBuilder = Analysis_FitModel_function2(Input_var,IndexCV,EnvironPara,ModelType)
disp('Analysis_FitModel_function2');

if(isempty(IndexCV))
   IndexCV = false(length(Input_var),1); 
end


%% 1:GCNon
%%variables that WILL be used as input variables
if(strcmp('GCNon',ModelType))
    try
        [Index_Var,~]=ismember(EnvironPara.VariableList,{...
            'ModelDataGEOSChem','ModelDataGEOSChem_Scale',...
            'ModelData_PriRD','ModelData_PriSecRD','ModelData_Elv','ModelData_Urban','ModelData_Pop','EmissionData_PM25',...
            'ModelData_PriRD_neighbor_1','ModelData_PriRD_neighbor_2','ModelData_PriRD_neighbor_3',...
            'ModelData_PriSecRD_neighbor_1','ModelData_PriSecRD_neighbor_2','ModelData_PriSecRD_neighbor_3',...
            'ModelData_Elv_neighbor_1','ModelData_Elv_neighbor_2','ModelData_Elv_neighbor_3',...
            'ModelData_Urban_neighbor_1','ModelData_Urban_neighbor_2','ModelData_Urban_neighbor_3',...
            'NCDataSite_AIRSFC','NCDataSite_APCP','NCDataSite_DSWRF','NCDataSite_EVAP','NCDataSite_HPBL','NCDataSite_LCDC','NCDataSite_PR_WTR','NCDataSite_PRATE','NCDataSite_VEG','NCDataSite_VIS','NCDataSite_PRESSFC','NCDataSite_RHUM2M','NCDataSite_WINDSPEED',...
            'NCDataSite_albedo','NCDataSite_hcdc','NCDataSite_mcdc'...
            });
        Index_temp = (1:1:EnvironPara.N_Var);
        Index_Var = Index_temp(Index_Var);
        % days we have available data in all data (including output variable and input variables)
        % used in training process
        Index_Day = isnan_matrix(Input_var(:,[1,Index_Var]))&~IndexCV;
        % days we have available data in input data
        % used in testing process
        Index_Day_1 = isnan_matrix(Input_var(:,Index_Var));
    %     ModelBuilder = cat(1,ModelBuilder,struct('NAME','GCNon','Index_Var',Index_Var,'Index_Day',Index_Day,'Index_Day_1',Index_Day_1));
        ModelBuilder = struct('NAME','GCNon','Index_Var',Index_Var,'Index_Day',Index_Day,'Index_Day_1',Index_Day_1);
    catch
        disp('constructing model builder errer, GCNon');
        ModelBuilder = struct('NAME','GCNon','Index_Var',[],'Index_Day',[],'Index_Day_1',[]);
    end
end

%% 2:GC
if(strcmp('GC',ModelType))
    try
        [Index_Var,~]=ismember(EnvironPara.VariableList,{...
            'ModelDataGEOSChem','ModelDataGEOSChem_Scale',...
            'ModelData_PriRD','ModelData_PriSecRD','ModelData_Elv','ModelData_Urban','ModelData_Pop','EmissionData_PM25',...
            'ModelData_PriRD_neighbor_1','ModelData_PriRD_neighbor_2','ModelData_PriRD_neighbor_3',...
            'ModelData_PriSecRD_neighbor_1','ModelData_PriSecRD_neighbor_2','ModelData_PriSecRD_neighbor_3',...
            'ModelData_Elv_neighbor_1','ModelData_Elv_neighbor_2','ModelData_Elv_neighbor_3',...
            'ModelData_Urban_neighbor_1','ModelData_Urban_neighbor_2','ModelData_Urban_neighbor_3',...
            'NCDataSite_AIRSFC','NCDataSite_APCP','NCDataSite_DSWRF','NCDataSite_EVAP','NCDataSite_HPBL','NCDataSite_LCDC','NCDataSite_PR_WTR','NCDataSite_PRATE','NCDataSite_VEG','NCDataSite_VIS','NCDataSite_PRESSFC','NCDataSite_RHUM2M','NCDataSite_WINDSPEED',...
            'NCDataSite_albedo','NCDataSite_hcdc','NCDataSite_mcdc',...
            'Spatial_Lagged',...
            'Spatial_Lagged_1','Spatial_Lagged_2','Spatial_Lagged_3','Temporal_Lagged_1','Temporal_Lagged_2','Temporal_Lagged_3'...
            });
        Index_temp = (1:1:EnvironPara.N_Var);
        Index_Var = Index_temp(Index_Var);
        Index_Day = isnan_matrix(Input_var(:,[1,Index_Var]))&~IndexCV;
        Index_Day_1 = isnan_matrix(Input_var(:,Index_Var));
        ModelBuilder = struct('NAME','GC','Index_Var',Index_Var,'Index_Day',Index_Day,'Index_Day_1',Index_Day_1);
    catch
        disp('constructing model builder errer, GC');
        ModelBuilder = struct('NAME','GC','Index_Var',[],'Index_Day',[],'Index_Day_1',[]);
    end
end

%% 3:AOD-Aqua
if(strcmp('AODAqua',ModelType))
    try
        [Index_Var,~]=ismember(EnvironPara.VariableList,{...
            'AODAqua_Data','AODAqua_Data_neighbor_1','AODAqua_Data_neighbor_2','AODAqua_Data_neighbor_3',...
            'ModelDataOMI',...
            'ModelDataGEOSChem',...
            'ModelData_PriRD','ModelData_PriSecRD','ModelData_Elv','ModelData_Urban','ModelData_Pop','EmissionData_PM25',...
            'ModelData_PriRD_neighbor_1','ModelData_PriRD_neighbor_2','ModelData_PriRD_neighbor_3',...
            'ModelData_PriSecRD_neighbor_1','ModelData_PriSecRD_neighbor_2','ModelData_PriSecRD_neighbor_3',...
            'ModelData_Elv_neighbor_1','ModelData_Elv_neighbor_2','ModelData_Elv_neighbor_3',...
            'ModelData_Urban_neighbor_1','ModelData_Urban_neighbor_2','ModelData_Urban_neighbor_3',...
            'NCDataSite_AIRSFC','NCDataSite_APCP','NCDataSite_DSWRF','NCDataSite_EVAP','NCDataSite_HPBL','NCDataSite_LCDC','NCDataSite_PR_WTR','NCDataSite_PRATE','NCDataSite_VEG','NCDataSite_VIS','NCDataSite_PRESSFC','NCDataSite_RHUM2M','NCDataSite_WINDSPEED',...
            'NCDataSite_albedo','NCDataSite_hcdc','NCDataSite_mcdc',...
            'Spatial_Lagged',...
            'Spatial_Lagged_1','Spatial_Lagged_2','Spatial_Lagged_3','Temporal_Lagged_1','Temporal_Lagged_2','Temporal_Lagged_3'...
            });
        Index_temp = (1:1:EnvironPara.N_Var);
        Index_Var = Index_temp(Index_Var);
        Index_Day = isnan_matrix(Input_var(:,[1,Index_Var]))&~IndexCV;
        Index_Day_1 = isnan_matrix(Input_var(:,Index_Var));
        ModelBuilder= struct('NAME','AODAqua','Index_Var',Index_Var,'Index_Day',Index_Day,'Index_Day_1',Index_Day_1);
    catch
        disp('constructing model builder errer, AODAqua');
        ModelBuilder= struct('NAME','AODAqua','Index_Var',[],'Index_Day',[],'Index_Day_1',[]);    
    end
end

%% 4:AOD-Terra
if(strcmp('AODTerra',ModelType))
    try
        [Index_Var,~]=ismember(EnvironPara.VariableList,{...
            'AODTerra_Data','AODTerra_Data_neighbor_1','AODTerra_Data_neighbor_2','AODTerra_Data_neighbor_3',...
            'ModelDataOMI',...
            'ModelDataGEOSChem_Scale',...
            'ModelData_PriRD','ModelData_PriSecRD','ModelData_Elv','ModelData_Urban','ModelData_Pop','EmissionData_PM25',...
            'ModelData_PriRD_neighbor_1','ModelData_PriRD_neighbor_2','ModelData_PriRD_neighbor_3',...
            'ModelData_PriSecRD_neighbor_1','ModelData_PriSecRD_neighbor_2','ModelData_PriSecRD_neighbor_3',...
            'ModelData_Elv_neighbor_1','ModelData_Elv_neighbor_2','ModelData_Elv_neighbor_3',...
            'ModelData_Urban_neighbor_1','ModelData_Urban_neighbor_2','ModelData_Urban_neighbor_3',...
            'NCDataSite_AIRSFC','NCDataSite_APCP','NCDataSite_DSWRF','NCDataSite_EVAP','NCDataSite_HPBL','NCDataSite_LCDC','NCDataSite_PR_WTR','NCDataSite_PRATE','NCDataSite_VEG','NCDataSite_VIS','NCDataSite_PRESSFC','NCDataSite_RHUM2M','NCDataSite_WINDSPEED',...
            'NCDataSite_albedo','NCDataSite_hcdc','NCDataSite_mcdc',...
            'Spatial_Lagged',...
            'Spatial_Lagged_1','Spatial_Lagged_2','Spatial_Lagged_3','Temporal_Lagged_1','Temporal_Lagged_2','Temporal_Lagged_3'...
            });
        Index_temp = (1:1:EnvironPara.N_Var);
        Index_Var = Index_temp(Index_Var);
        Index_Day = isnan_matrix(Input_var(:,[1,Index_Var]))&~IndexCV;
        Index_Day_1 = isnan_matrix(Input_var(:,Index_Var));
        ModelBuilder= struct('NAME','AODTerra','Index_Var',Index_Var,'Index_Day',Index_Day,'Index_Day_1',Index_Day_1);
    catch
        disp('constructing model builder errer, AODTerra');
        ModelBuilder= struct('NAME','AODTerra','Index_Var',[],'Index_Day',[],'Index_Day_1',[]);
    end
end

%% 5:Final model
if(strcmp('FINAL',ModelType))
    try
        [Index_Var,~]=ismember(EnvironPara.VariableList,{...
            'ModelDataGEOSChem','ModelDataGEOSChem_Scale',...
            'ModelData_PriRD','ModelData_PriSecRD','ModelData_Elv','ModelData_Urban','ModelData_Pop','EmissionData_PM25',...
            'ModelData_PriRD_neighbor_1','ModelData_PriRD_neighbor_2','ModelData_PriRD_neighbor_3',...
            'ModelData_PriSecRD_neighbor_1','ModelData_PriSecRD_neighbor_2','ModelData_PriSecRD_neighbor_3',...
            'ModelData_Elv_neighbor_1','ModelData_Elv_neighbor_2','ModelData_Elv_neighbor_3',...
            'ModelData_Urban_neighbor_1','ModelData_Urban_neighbor_2','ModelData_Urban_neighbor_3',...
            'NCDataSite_AIRSFC','NCDataSite_APCP','NCDataSite_DSWRF','NCDataSite_EVAP','NCDataSite_HPBL','NCDataSite_LCDC','NCDataSite_PR_WTR','NCDataSite_PRATE','NCDataSite_VEG','NCDataSite_VIS','NCDataSite_PRESSFC','NCDataSite_RHUM2M','NCDataSite_WINDSPEED',...
            'NCDataSite_albedo','NCDataSite_hcdc','NCDataSite_mcdc',...
            'MonitorPredicted_GC','MonitorPredicted_AODAqua','MonitorPredicted_AODTerra'...
            });
        Index_temp = (1:1:EnvironPara.N_Var);
        Index_Var = Index_temp(Index_Var);
        Index_Day = isnan_matrix(Input_var(:,[1,Index_Var]))&~IndexCV;
        Index_Day_1 = isnan_matrix(Input_var(:,Index_Var));
        ModelBuilder= struct('NAME','FINAL','Index_Var',Index_Var,'Index_Day',Index_Day,'Index_Day_1',Index_Day_1);
    catch
        disp('constructing model builder errer, FINAL');
        ModelBuilder= struct('NAME','FINAL','Index_Var',[],'Index_Day',[],'Index_Day_1',[]);    
    end
end

%% 6:General Model
if(strcmp('General',ModelType))
    Index_Var = 2:1:EnvironPara.N_Var;
    Index_Day = isnan_matrix(Input_var(:,[1,Index_Var]))&~IndexCV;
    Index_Day_1 = isnan_matrix(Input_var(:,Index_Var));
    ModelBuilder= struct('NAME','General','Index_Var',Index_Var,'Index_Day',Index_Day,'Index_Day_1',Index_Day_1);
end

%% 7:random General Model AODAqua
if(strcmp('RandAquaNearby',ModelType))
    try
        [Index_Var,~]=ismember(EnvironPara.VariableList,{...
            'AODAqua_Data_Nearby',...
            'ModelDataOMI_UV','ModelDataOMI_VIS'...
            'ModelDataGEOSChem','ModelDataGEOSChemSum',...
            'ModelData_PriRD','ModelData_PriSecRD','ModelData_Elv','ModelData_Urban','ModelData_Pop','EmissionData_PM25',...
            'NCDataSite_AIRSFC','NCDataSite_APCP','NCDataSite_DSWRF','NCDataSite_EVAP','NCDataSite_HPBL','NCDataSite_LCDC','NCDataSite_PR_WTR','NCDataSite_PRATE','NCDataSite_VEG','NCDataSite_VIS','NCDataSite_PRESSFC','NCDataSite_RHUM2M','NCDataSite_WINDSPEED',...
            'NCDataSite_albedo','NCDataSite_hcdc','NCDataSite_mcdc'...  
            });
        Index_temp = (1:1:EnvironPara.N_Var);
        Index_Var = Index_temp(Index_Var);
        % take a subset of training data, which is labelled by Index_Day
        Index_Day = isnan_matrix(Input_var(:,[1,Index_Var]));
        Temp = (1:length(Index_Day))';
        Index_Day = Temp(Index_Day);
        Index_Day_subset = datasample(Index_Day,min(length(Index_Day),EnvironPara.SizeTrainingData),'Replace',false);
        % testing data label
        Index_Day_1 = isnan_matrix(Input_var(:,Index_Var));
        ModelBuilder= struct('NAME',ModelType,'Index_Var',Index_Var,'Index_Day',Index_Day_subset,'Index_Day_1',Index_Day_1);
    catch
        disp(['constructing model builder errer',ModelType]);
        ModelBuilder= struct('NAME',ModelType,'Index_Var',[],'Index_Day',[],'Index_Day_1',[]);
    end
end


%% 8:random General Model TerraAqua
if(strcmp('RandTerraNearby',ModelType))
    try
        [Index_Var,~]=ismember(EnvironPara.VariableList,{...
            'AODTerra_Data_Nearby',...
            'ModelDataOMI_UV','ModelDataOMI_VIS'...
            'ModelDataGEOSChem','ModelDataGEOSChemSum',...
            'ModelData_PriRD','ModelData_PriSecRD','ModelData_Elv','ModelData_Urban','ModelData_Pop','EmissionData_PM25',...
            'NCDataSite_AIRSFC','NCDataSite_APCP','NCDataSite_DSWRF','NCDataSite_EVAP','NCDataSite_HPBL','NCDataSite_LCDC','NCDataSite_PR_WTR','NCDataSite_PRATE','NCDataSite_VEG','NCDataSite_VIS','NCDataSite_PRESSFC','NCDataSite_RHUM2M','NCDataSite_WINDSPEED',...
            'NCDataSite_albedo','NCDataSite_hcdc','NCDataSite_mcdc'...  
            });
        Index_temp = (1:1:EnvironPara.N_Var);
        Index_Var = Index_temp(Index_Var);
        % take a subset of training data, which is labelled by Index_Day
        Index_Day = isnan_matrix(Input_var(:,[1,Index_Var]));
        Temp = (1:length(Index_Day))';
        Index_Day = Temp(Index_Day);
        Index_Day_subset = datasample(Index_Day,min(length(Index_Day),EnvironPara.SizeTrainingData),'Replace',false);
        % testing data label
        Index_Day_1 = isnan_matrix(Input_var(:,Index_Var));
        ModelBuilder= struct('NAME',ModelType,'Index_Var',Index_Var,'Index_Day',Index_Day_subset,'Index_Day_1',Index_Day_1);
    catch
        disp(['constructing model builder errer',ModelType]);
        ModelBuilder= struct('NAME',ModelType,'Index_Var',[],'Index_Day',[],'Index_Day_1',[]);
    end
end

%% 9:random General Model GC
if(strcmp('RandGC',ModelType))
    try
        [Index_Var,~]=ismember(EnvironPara.VariableList,{...
            'ModelDataGEOSChemSum',...
            'ModelData_PriRD','ModelData_PriSecRD','ModelData_Elv','ModelData_Urban','ModelData_Pop','EmissionData_PM25',...
            'NCDataSite_AIRSFC','NCDataSite_APCP','NCDataSite_DSWRF','NCDataSite_EVAP','NCDataSite_HPBL','NCDataSite_LCDC','NCDataSite_PR_WTR','NCDataSite_PRATE','NCDataSite_VEG','NCDataSite_VIS','NCDataSite_PRESSFC','NCDataSite_RHUM2M','NCDataSite_WINDSPEED',...
            'NCDataSite_albedo','NCDataSite_hcdc','NCDataSite_mcdc'...  
            });
        Index_temp = (1:1:EnvironPara.N_Var);
        Index_Var = Index_temp(Index_Var);
        % take a subset of training data, which is labelled by Index_Day
        Index_Day = isnan_matrix(Input_var(:,[1,Index_Var]));
        Temp = (1:length(Index_Day))';
        Index_Day = Temp(Index_Day);
        Index_Day_subset = datasample(Index_Day,min(length(Index_Day),EnvironPara.SizeTrainingData),'Replace',false);
        % testing data label
        Index_Day_1 = isnan_matrix(Input_var(:,Index_Var));
        ModelBuilder= struct('NAME',ModelType,'Index_Var',Index_Var,'Index_Day',Index_Day_subset,'Index_Day_1',Index_Day_1);
    catch
        disp(['constructing model builder errer',ModelType]);
        ModelBuilder= struct('NAME',ModelType,'Index_Var',[],'Index_Day',[],'Index_Day_1',[]);
    end
end

%% 9:random General Model OMI
if(strcmp('RandOMI',ModelType))
    try
        [Index_Var,~]=ismember(EnvironPara.VariableList,{...
            'ModelDataOMI_UV','ModelDataOMI_VIS',...
            'ModelData_PriRD','ModelData_PriSecRD','ModelData_Elv','ModelData_Urban','ModelData_Pop','EmissionData_PM25',...
            'NCDataSite_AIRSFC','NCDataSite_APCP','NCDataSite_DSWRF','NCDataSite_EVAP','NCDataSite_HPBL','NCDataSite_LCDC','NCDataSite_PR_WTR','NCDataSite_PRATE','NCDataSite_VEG','NCDataSite_VIS','NCDataSite_PRESSFC','NCDataSite_RHUM2M','NCDataSite_WINDSPEED',...
            'NCDataSite_albedo','NCDataSite_hcdc','NCDataSite_mcdc'...  
            });
        Index_temp = (1:1:EnvironPara.N_Var);
        Index_Var = Index_temp(Index_Var);
        % take a subset of training data, which is labelled by Index_Day
        Index_Day = isnan_matrix(Input_var(:,[1,Index_Var]));
        Temp = (1:length(Index_Day))';
        Index_Day = Temp(Index_Day);
        Index_Day_subset = datasample(Index_Day,min(length(Index_Day),EnvironPara.SizeTrainingData),'Replace',false);
        % testing data label
        Index_Day_1 = isnan_matrix(Input_var(:,Index_Var));
        ModelBuilder= struct('NAME',ModelType,'Index_Var',Index_Var,'Index_Day',Index_Day_subset,'Index_Day_1',Index_Day_1);
    catch
        disp(['constructing model builder errer',ModelType]);
        ModelBuilder= struct('NAME',ModelType,'Index_Var',[],'Index_Day',[],'Index_Day_1',[]);
    end
end

%% 9:non-random General Model GC
if(strcmp('NonRandGC',ModelType))
    try
        [Index_Var,~]=ismember(EnvironPara.VariableList,{...
            'ModelDataGEOSChemSum',...
            'ModelData_PriRD','ModelData_PriSecRD','ModelData_Elv','ModelData_Urban','ModelData_Pop','EmissionData_PM25',...
            'NCDataSite_AIRSFC','NCDataSite_APCP','NCDataSite_DSWRF','NCDataSite_EVAP','NCDataSite_HPBL','NCDataSite_LCDC','NCDataSite_PR_WTR','NCDataSite_PRATE','NCDataSite_VEG','NCDataSite_VIS','NCDataSite_PRESSFC','NCDataSite_RHUM2M','NCDataSite_WINDSPEED',...
            'NCDataSite_albedo','NCDataSite_hcdc','NCDataSite_mcdc'...  
            });
        Index_temp = (1:1:EnvironPara.N_Var);
        Index_Var = Index_temp(Index_Var);
        Index_Day = isnan_matrix(Input_var(:,[1,Index_Var]))&~IndexCV;
        Index_Day_1 = isnan_matrix(Input_var(:,Index_Var));
        ModelBuilder= struct('NAME','NonRandGC','Index_Var',Index_Var,'Index_Day',Index_Day,'Index_Day_1',Index_Day_1);
    catch
        disp(['constructing model builder errer',ModelType]);
        ModelBuilder= struct('NAME',ModelType,'Index_Var',[],'Index_Day',[],'Index_Day_1',[]);
    end
end

%% 9:random General Model GC only
if(strcmp('RandGCOnly',ModelType))
    try
        [Index_Var,~]=ismember(EnvironPara.VariableList,{...
            'ModelDataGEOSChemSum'...
            });
        Index_temp = (1:1:EnvironPara.N_Var);
        Index_Var = Index_temp(Index_Var);
        % take a subset of training data, which is labelled by Index_Day
        Index_Day = isnan_matrix(Input_var(:,[1,Index_Var]));
        Temp = (1:length(Index_Day))';
        Index_Day = Temp(Index_Day);
        Index_Day_subset = datasample(Index_Day,min(length(Index_Day),EnvironPara.SizeTrainingData),'Replace',false);
        % testing data label
        Index_Day_1 = isnan_matrix(Input_var(:,Index_Var));
        ModelBuilder= struct('NAME',ModelType,'Index_Var',Index_Var,'Index_Day',Index_Day_subset,'Index_Day_1',Index_Day_1);
    catch
        disp(['constructing model builder errer',ModelType]);
        ModelBuilder= struct('NAME',ModelType,'Index_Var',[],'Index_Day',[],'Index_Day_1',[]);
    end
end

%% 6:General Model
if(strcmp('RandGeneral',ModelType))
    Index_Var = 2:1:EnvironPara.N_Var;
    % take a subset of training data, which is labelled by Index_Day
    Index_Day = isnan_matrix(Input_var(:,[1,Index_Var]));
    Temp = (1:length(Index_Day))';
    Index_Day = Temp(Index_Day);
    Index_Day_subset = datasample(Index_Day,min(length(Index_Day),EnvironPara.SizeTrainingData),'Replace',false);
    % testing data label
    Index_Day_1 = isnan_matrix(Input_var(:,Index_Var));
    ModelBuilder= struct('NAME',ModelType,'Index_Var',Index_Var,'Index_Day',Index_Day_subset,'Index_Day_1',Index_Day_1);
end

