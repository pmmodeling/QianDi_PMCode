%% just a component to do the fitting model, updating lagged terms stuff.
%%avoid duplicate codes
%%Input_var: input variables; CurrentMonitorPredicted: current modeling
%%results; ModelBuilder:some parameters (just the struct variable for one model, not the whole list); MaxRun: how many times you would
%%like the model to run. Inf if running recursively. MaxFail: what
%%is the maxinum failure time that model can fit. MaxFail<MaxRun, otherwise
%%meaningless
%%count: how many time the model has run
%%2016-08-29: change name from CalibrateGC_FitModel_function1 to Analysis_FitModel_function1
function [Input_var,MonitorPredicted,NetPredicted,Count] = Analysis_FitModel_function1(Input_var,MonitorPredicted,NetPredicted,ModelBuilder,EnvironPara,MaxRun,MaxFail,Count)
disp('Analysis_FitModel_function1');
%% initilization parameter and output
Index_Var = ModelBuilder.Index_Var;
% days we have available data in all data (including output variable and input variables)
Index_Day = ModelBuilder.Index_Day;
% days we have available data in input data
Index_Day_1 = ModelBuilder.Index_Day_1;
FailureTime = 0;

%% running the model

fprintf('%s\n',repmat('.',[1,50]));
fprintf('count:%d\trunning...%s\n',Count,ModelBuilder.NAME); 
%%temp output
MonitorPredicted_temp = nan(length(Input_var),1);

%%calculate current model fit
Temp = reshape(MonitorPredicted,[EnvironPara.N_Day*EnvironPara.N_Site,1]);
[CurrentR,~,~] = CalculateRsquare(Input_var(Index_Day_1,1),Temp(Index_Day_1,:),'out-of-sample');
fprintf('Current R2 is %d\n',CurrentR); 

%%run the model!!!
fprintf('%s\n',strjoin(EnvironPara.VariableList(Index_Var),'\n'));
%     fprintf(EnvironPara.FID,'%s\n',strjoin(EnvironPara.VariableList(Index_Var),'\n'));

fprintf('%d,%d\n',size(Input_var(Index_Day,Index_Var),1),size(Input_var(Index_Day,Index_Var),2))
[MonitorPredicted_temp(Index_Day_1),NetPredicted_temp] = Analysis_NeuralNetwork(Input_var(Index_Day,1),Input_var(Index_Day,Index_Var),Input_var(Index_Day_1,Index_Var),EnvironPara.N_Layer,[],EnvironPara);
try
    TempVariable = EnvironPara.VariableList(Index_Var);
    fprintf('variables used by neural network:\n%s\n',strjoin(TempVariable(NetPredicted_temp{3,1}.keep),'\n'));
catch

end
%%calculate R2
[R,~,~] = CalculateRsquare(Input_var(Index_Day_1,1),MonitorPredicted_temp(Index_Day_1,:),'out-of-sample');

%%there is some improvement
%     if(R>CurrentR)
if(R>0)
    disp([ModelBuilder.NAME,' improved!!!']);
    %%store results
    MonitorPredicted = reshape(MonitorPredicted,[EnvironPara.N_Day*EnvironPara.N_Site,1]);
    MonitorPredicted(Index_Day_1) = MonitorPredicted_temp(Index_Day_1);
    MonitorPredicted = reshape(MonitorPredicted,[EnvironPara.N_Day,EnvironPara.N_Site]);
    NetPredicted_temp = cat(1,NetPredicted_temp,ModelBuilder.NAME);
    NetPredicted_temp = cat(1,NetPredicted_temp,MonitorPredicted);%store every step's interim results;

    % update estimation-based terms -- if needed
    [Index_temp,~]=ismember(EnvironPara.VariableList,'Spatial_Lagged_1');
    if(sum(Index_temp)>0 && MaxRun == 1)% there is estimation-based terms. we should update them
        All_lagged_terms = Analysis_TwoStep(MonitorPredicted,[],[],EnvironPara.SITENAME_PREDICT,EnvironPara.SITENAME_MODEL,EnvironPara,'UPDATE_IN_TRAINING');
        NetPredicted_temp = cat(1,NetPredicted_temp,All_lagged_terms);
        Input_var = MatchVariablebyName(Input_var,All_lagged_terms(:,1),EnvironPara.VariableList,'Spatial_Lagged_1');
        Input_var = MatchVariablebyName(Input_var,All_lagged_terms(:,2),EnvironPara.VariableList,'Spatial_Lagged_2');
        Input_var = MatchVariablebyName(Input_var,All_lagged_terms(:,3),EnvironPara.VariableList,'Spatial_Lagged_3');
        Input_var = MatchVariablebyName(Input_var,All_lagged_terms(:,4),EnvironPara.VariableList,'Temporal_Lagged_1');
        Input_var = MatchVariablebyName(Input_var,All_lagged_terms(:,5),EnvironPara.VariableList,'Temporal_Lagged_2');
        Input_var = MatchVariablebyName(Input_var,All_lagged_terms(:,6),EnvironPara.VariableList,'Temporal_Lagged_3');
    end
    %%count
    Count = Count + 1;
    FailureTime = 0;
    %%record results
    NetPredicted = cat(1,NetPredicted,{NetPredicted_temp});
    if(EnvironPara.IsRecord)
        Analysis_ResultDescriptiveAnalysis(MonitorPredicted,reshape(Input_var(:,1),[EnvironPara.N_Day,EnvironPara.N_Site]),NetPredicted,EnvironPara,[num2str(Count),ModelBuilder.NAME],{true,false,false,false});
        Analysis_InputDataDescriptiveAnalysis(Input_var,EnvironPara,true,[num2str(Count),ModelBuilder.NAME],{'SUBSET1',false,false});
    else
        Analysis_ResultDescriptiveAnalysis(MonitorPredicted,reshape(Input_var(:,1),[EnvironPara.N_Day,EnvironPara.N_Site]),NetPredicted,EnvironPara,[num2str(Count),ModelBuilder.NAME],{false,false,false,false});
        Analysis_InputDataDescriptiveAnalysis(Input_var,EnvironPara,false,[num2str(Count),ModelBuilder.NAME],{'SUBSET1',false,false});
    end
else
    FailureTime = FailureTime+1;
    fprintf('%s does NOT work...failing %d\n',ModelBuilder.NAME,FailureTime);
end

if(FailureTime>=MaxFail)
   return; 
end   

end

%% assign variable to the data set based on variable names
function Input_var = MatchVariablebyName(Input_var,Data,VariableList,VariableName)
    try
        [Index_var,~]=ismember(VariableList,VariableName);
        Input_var(:,Index_var)=Data;
    catch

    end

end

