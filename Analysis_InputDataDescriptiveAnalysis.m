%% evalute/make descriptive stat of input variables
%%version1: OPTION1, evalute the whole thing or subset?
%%OPTION2: save figures? OPTION3: save figure for land-use variables?
function Analysis_InputDataDescriptiveAnalysis(Input_var,EnvironPara,IsSaved,NOTICE,OPTION)
disp('Analysis_InputDataDescriptiveAnalysis');

N_var = size(Input_var,2);
VariableList = EnvironPara.VariableList;
if(isempty(EnvironPara.VariableList))
   VariableList = cell(N_var,1);
   for i =1:N_var
       VariableList{i} = ['Variable ',num2str(i)];
   end
end

if(strcmp(OPTION{1},'SUBSET1'))
    Locb1 = ~cellfun(@isempty,(regexp(EnvironPara.VariableList,'Spatial_Lagged_+|Temporal_Lagged_+')));
    Input_var = Input_var(:,Locb1);
    VariableList = VariableList(Locb1);
    N_var = size(Input_var,2);
elseif(strcmp(OPTION{1},'WHOLE'))
    Locb1 = ones(1,N_var);
end

%% we need to check
if(~EnvironPara.IsTest)

    if(sum(Locb1)>0)
        %% descriptive analysis
        N = length(Input_var);
        N1 = sum(isnan_matrix(Input_var));
        N2 = sum(isnan_matrix(Input_var(:,2:N_var)));
        N3 = sum(isnan_matrix(Input_var(:,1)));
        N_List = sum(~isnan(Input_var));

        fprintf('\n%s\t',NOTICE);
        fprintf('total variable:%d\n',N_var);
        fprintf('total records:%d\n',N);
        fprintf('complete records:%d\n',N1);
        fprintf('complete records(%%):%f\n',N1/N);
        fprintf('complete input:%d\n',N2);
        fprintf('complete input(%%):%f\n',N2/N);
        fprintf('monitoring data:%d\n',N3);
        fprintf('monitoring data (%%):%f\n',N3/N);
        for i = 1:length(N_List)
           fprintf('%s\t%s\tcomplete records:%d\tcomplete %%:%f\tis real value:%d\t',EnvironPara.CommonStartPoint,VariableList{i},N_List(i),N_List(i)/N,isreal(Input_var(:,i))); 
           fprintf('min:%f\t5%%:%f\t20%%:%f\tmean:%f\t80%%:%f\t95%%:%f\tmax:%f\n',nanmin(Input_var(:,i)),quantile(Input_var(:,i),0.05),quantile(Input_var(:,i),0.2),nanmean(Input_var(:,i)),quantile(Input_var(:,i),0.8),quantile(Input_var(:,i),0.95),nanmax(Input_var(:,i))); 
        end

        if(IsSaved)
            fprintf(EnvironPara.FID_CHECK,'\n%s\t',NOTICE);
            fprintf(EnvironPara.FID_CHECK,'total variable:%d\n',N_var);
            fprintf(EnvironPara.FID_CHECK,'total records:%d\n',N);
            fprintf(EnvironPara.FID_CHECK,'complete records:%d\n',N1);
            fprintf(EnvironPara.FID_CHECK,'complete records(%%):%f\n',N1/N);
            fprintf(EnvironPara.FID_CHECK,'complete input:%d\n',N2);
            fprintf(EnvironPara.FID_CHECK,'complete input(%%):%f\n',N2/N);
            fprintf(EnvironPara.FID_CHECK,'monitoring data:%d\n',N3);
            fprintf(EnvironPara.FID_CHECK,'monitoring data (%%):%f\n',N3/N);
            for i = 1:length(N_List)
                fprintf(EnvironPara.FID_CHECK,'%s\t%s\tcomplete records:%d\tcomplete %%:%f\tis real value:%d\t',EnvironPara.CommonStartPoint,VariableList{i},N_List(i),N_List(i)/N,isreal(Input_var(:,i)));
                fprintf(EnvironPara.FID_CHECK,'min:%f\t5%%:%f\t20%%:%f\tmean:%f\t80%%:%f\t95%%:%f\tmax:%f\n',nanmin(Input_var(:,i)),quantile(Input_var(:,i),0.05),quantile(Input_var(:,i),0.2),nanmean(Input_var(:,i)),quantile(Input_var(:,i),0.8),quantile(Input_var(:,i),0.95),nanmax(Input_var(:,i)));
            end
        end 
    end
    
    %% visual check 1   
    if(strcmp(OPTION{1},'WHOLE') && OPTION{2})
        %% land-use variable
        if(OPTION{3})
            TempList = {...
            'ModelData_PriRD','ModelData_PriSecRD','ModelData_Elv','ModelData_Urban','ModelData_Pop','EmissionData_PM25',...
            'ModelData_PriRD_neighbor_1','ModelData_PriRD_neighbor_2','ModelData_PriRD_neighbor_3',...
            'ModelData_PriSecRD_neighbor_1','ModelData_PriSecRD_neighbor_2','ModelData_PriSecRD_neighbor_3',...
            'ModelData_Elv_neighbor_1','ModelData_Elv_neighbor_2','ModelData_Elv_neighbor_3',...
            'ModelData_Urban_neighbor_1','ModelData_Urban_neighbor_2','ModelData_Urban_neighbor_3',...
            'ModelData_Pop_neighbor_1','ModelData_Pop_neighbor_2','ModelData_Pop_neighbor_3'...
            };
            [Index_Var,Index_Locb]=ismember(TempList,VariableList);
            for i=1:length(Index_Var)
                if(0==Index_Locb(i))
                   continue; 
                end
                
                FileName = strrep(TempList{i},'_',' ');
            end
        end

        %% AOD, meteorological       
        TempList = {...
            'AODAqua_Data_neighbor_1','AODAqua_Data_neighbor_2','AODAqua_Data_neighbor_3',...
            'AODTerra_Data_neighbor_1','AODTerra_Data_neighbor_2','AODTerra_Data_neighbor_3',...
            'ModelDataOMI',...
            'ModelDataGEOSChem','ModelDataGEOSChem_Scale',...
                'NCDataSite_AIRSFC','NCDataSite_APCP','NCDataSite_DSWRF','NCDataSite_EVAP','NCDataSite_HPBL','NCDataSite_LCDC','NCDataSite_PR_WTR','NCDataSite_PRATE','NCDataSite_VEG','NCDataSite_VIS','NCDataSite_PRESSFC','NCDataSite_RHUM2M','NCDataSite_WINDSPEED',...
            'NCDataSite_albedo','NCDataSite_hcdc','NCDataSite_mcdc',...
            'Spatial_Lagged'...
            };
        [Index_Var,Index_Locb]=ismember(TempList,VariableList);
        %% check every week
        if(rem(EnvironPara.StartDate,7)==0)
            for i=1:length(Index_Var)
                if(0==Index_Locb(i))
                   continue; 
                end
                
                FileName = [EnvironPara.CommonStartPoint,' ',strrep(TempList{i},'_',' ')];
            end
        end
    end
    if(strcmp(OPTION{1},'SUBSET1')&& OPTION{2})
        %% lagged terms        
        TempList = {...
            'Spatial_Lagged_1','Spatial_Lagged_2','Spatial_Lagged_3','Temporal_Lagged_1','Temporal_Lagged_2','Temporal_Lagged_3'...
            };
        [Index_Var,Index_Locb]=ismember(TempList,VariableList);
        %% check every week
        if(rem(EnvironPara.StartDate,1)==0)
            for i=1:length(Index_Var)
                if(0==Index_Locb(i))
                   continue; 
                end
                
                FileName = [EnvironPara.CommonStartPoint,' ',strrep(TempList{i},'_',' ')];
            end
        end
        
    end
end
