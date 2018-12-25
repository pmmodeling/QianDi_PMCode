%% the main function of fitting model (GC,AOD) for CV cross-validation; all; prediction
%%Index: index matrix used for CV cross validation
%%Input_var_Grid,NetPredicted are data and pre-trained neural network used
%%for grid cell prediction
%%output in N*1 format
%%MonitorPredicted: prior estimation

%%version1: can take weight (convolution) as inputs, fit the model recursively
%%version2: fit GEOS-Chem first and use that to fit AOD (a hybrid model)
%%version2.1: the order of fitting GC and AOD model is different
%%version2.2 the order of fitting model is different:
%%non-GC-->GC-->AOD(aqua+terra)/2;
%%version2.3: save intermittent modeling results;
%%version2.3.1: deal with MonitorPredicted --- which can fit into the model
%%in the final run;
%% 2016-08-29: change name from CalibrateGC_FitModel_2_3_1 to Analysis_FitModel

function [Input_var,MonitorPredicted,NetPredicted] = Analysis_FitModel(Input_var,MonitorPredicted,NetPredicted,Count,IndexCV,N_Site,N_Day,EnvironPara)
disp('Analysis_FitModel');

%% number of variables
% N_All = length(Input_var);
% ModelBuilder = CalibrateGC_FitModel_function2_1(Input_var,EnvironPara);
MonitorData = reshape(Input_var(:,1),[EnvironPara.N_Day,EnvironPara.N_Site]);

if(strcmp(EnvironPara.OPTION,'CV_round2'))
    MonitorPredicted_1 = MonitorPredicted;
    MonitorPredicted = reshape(MonitorPredicted,[EnvironPara.N_Day*EnvironPara.N_Site,1]);
    
    i=1;
    try
        load([EnvironPara.OUTPUTPATH,'TempOutput_',EnvironPara.NAME,'_',EnvironPara.SITENAME_MODEL,'_',EnvironPara.CommonStartPoint,'_',EnvironPara.CommonEndPoint,'.mat']);
        disp('previous temporary results found! the code is to continue!');
    catch
        disp('NO previous temporary results found... the code starts from stratch!');
    end
    
    while i<=size(IndexCV,2)
        fprintf('%s\n',repmat('*',[1,50]));
        fprintf('CV in process!....%d\n',i);
       
        Count_temp = Count;
        MonitorPredicted_temp = MonitorPredicted_1;
        
        %%round 1
        EnvironPara.ModelName = 'Hybrid';
        ModelBuilder = Analysis_FitModel_function2(Input_var,IndexCV(:,i),EnvironPara,'General');
        [Input_var,MonitorPredicted_temp,NetPredicted,Count_temp] =...
            Analysis_FitModel_function1(Input_var,MonitorPredicted_temp,NetPredicted,ModelBuilder,EnvironPara,1,Inf,Count_temp);
        
        %%round 2
        EnvironPara.ModelName = 'NeuralNetwork';
        ModelBuilder = Analysis_FitModel_function2(Input_var,IndexCV(:,i),EnvironPara,'General');
        [Input_var,MonitorPredicted_temp,NetPredicted,Count_temp] =...
            Analysis_FitModel_function1(Input_var,MonitorPredicted_temp,NetPredicted,ModelBuilder,EnvironPara,2,Inf,Count_temp);
        
        %%record results in the testing dataset
        MonitorPredicted(IndexCV(:,i)) = MonitorPredicted_temp(IndexCV(:,i));
        
        Analysis_ResultDescriptiveAnalysis(MonitorPredicted(IndexCV(:,i)),MonitorData(IndexCV(:,i)),NetPredicted,EnvironPara,[num2str(i),ModelBuilder.NAME,'_',EnvironPara.OPTION],{true,false,false,false});
        
        i = i+1;
                
        save([EnvironPara.OUTPUTPATH,'TempOutput_',EnvironPara.NAME,'_',EnvironPara.SITENAME_MODEL,'_',EnvironPara.CommonStartPoint,'_',EnvironPara.CommonEndPoint,'.mat'],...
            'i','MonitorPredicted','NetPredicted','MonitorPredicted_1','Count','IndexCV','-v7.3');
    end
    
    MonitorPredicted = reshape(MonitorPredicted,[EnvironPara.N_Day,EnvironPara.N_Site]);
elseif(strcmp(EnvironPara.OPTION,'CV_round3'))
    MonitorPredicted_1 = MonitorPredicted;
    MonitorPredicted = reshape(MonitorPredicted,[EnvironPara.N_Day*EnvironPara.N_Site,1]);
    
    i=1;
    try
        load([EnvironPara.OUTPUTPATH,'TempOutput_',EnvironPara.NAME,'_',EnvironPara.SITENAME_MODEL,'_',EnvironPara.CommonStartPoint,'_',EnvironPara.CommonEndPoint,'.mat']);
        disp('previous temporary results found! the code is to continue!');
    catch
        disp('NO previous temporary results found... the code starts from stratch!');
    end
    
    while i<=size(IndexCV,2)
        fprintf('%s\n',repmat('*',[1,50]));
        fprintf('CV in process!....%d\n',i);
       
        Count_temp = Count;
        MonitorPredicted_temp = MonitorPredicted_1;
        
        %%round 1
        EnvironPara.ModelName = 'NeuralNetwork';
        ModelBuilder = Analysis_FitModel_function2(Input_var,IndexCV(:,i),EnvironPara,'General');
        [Input_var,MonitorPredicted_temp,NetPredicted,Count_temp] =...
            Analysis_FitModel_function1(Input_var,MonitorPredicted_temp,NetPredicted,ModelBuilder,EnvironPara,1,Inf,Count_temp);
        
        %%round 2
        EnvironPara.ModelName = 'NeuralNetwork';
        ModelBuilder = Analysis_FitModel_function2(Input_var,IndexCV(:,i),EnvironPara,'General');
        [Input_var,MonitorPredicted_temp,NetPredicted,Count_temp] =...
            Analysis_FitModel_function1(Input_var,MonitorPredicted_temp,NetPredicted,ModelBuilder,EnvironPara,2,Inf,Count_temp);
        
        %%record results in the testing dataset
        MonitorPredicted(IndexCV(:,i)) = MonitorPredicted_temp(IndexCV(:,i));
        
        Analysis_ResultDescriptiveAnalysis(MonitorPredicted(IndexCV(:,i)),MonitorData(IndexCV(:,i)),NetPredicted,EnvironPara,[num2str(i),ModelBuilder.NAME,'_',EnvironPara.OPTION],{true,false,false,false});
        
        i = i+1;
                
        save([EnvironPara.OUTPUTPATH,'TempOutput_',EnvironPara.NAME,'_',EnvironPara.SITENAME_MODEL,'_',EnvironPara.CommonStartPoint,'_',EnvironPara.CommonEndPoint,'.mat'],...
            'i','MonitorPredicted','NetPredicted','MonitorPredicted_1','Count','IndexCV','-v7.3');
    end
    
    MonitorPredicted = reshape(MonitorPredicted,[EnvironPara.N_Day,EnvironPara.N_Site]);
    
   
elseif(strcmp(EnvironPara.OPTION,'CV_round1')) 
    MonitorPredicted_1 = MonitorPredicted;
    MonitorPredicted = reshape(MonitorPredicted,[EnvironPara.N_Day*EnvironPara.N_Site,1]);
    
    i=1;
    try
        load([EnvironPara.OUTPUTPATH,'TempOutput_',EnvironPara.NAME,'_',EnvironPara.SITENAME_MODEL,'_',EnvironPara.CommonStartPoint,'_',EnvironPara.CommonEndPoint,'.mat']);
        disp('previous temporary results found! the code is to continue!');
    catch
        disp('NO previous temporary results found... the code starts from stratch!');
    end
    
    while i<=size(IndexCV,2)
        fprintf('%s\n',repmat('*',[1,50]));
        fprintf('CV in process!....%d\n',i);
       
        Count_temp = Count;
        MonitorPredicted_temp = MonitorPredicted_1;
        %%round 1
        ModelBuilder = Analysis_FitModel_function2(Input_var,IndexCV(:,i),EnvironPara,'General');
        [Input_var,MonitorPredicted_temp,NetPredicted,Count_temp] =...
            Analysis_FitModel_function1(Input_var,MonitorPredicted_temp,NetPredicted,ModelBuilder,EnvironPara,1,Inf,Count_temp);
               
        %%record results in the testing dataset
        MonitorPredicted(IndexCV(:,i)) = MonitorPredicted_temp(IndexCV(:,i));
        
        Analysis_ResultDescriptiveAnalysis(MonitorPredicted(IndexCV(:,i)),MonitorData(IndexCV(:,i)),NetPredicted,EnvironPara,[num2str(i),ModelBuilder.NAME,'_',EnvironPara.OPTION],{true,false,false,false});
 
        i = i+1;
        
        save([EnvironPara.OUTPUTPATH,'TempOutput_',EnvironPara.NAME,'_',EnvironPara.SITENAME_MODEL,'_',EnvironPara.CommonStartPoint,'_',EnvironPara.CommonEndPoint,'.mat'],...
            'i','MonitorPredicted','NetPredicted','MonitorPredicted_1','Count','IndexCV','-v7.3');
    end
    
    MonitorPredicted = reshape(MonitorPredicted,[EnvironPara.N_Day,EnvironPara.N_Site]);
    
    
    %%%%%%%%%%%%% the prediction we are using!!!
elseif(strcmp(EnvironPara.OPTION,'Prediction'))
   %% just for AOD prediction
    fprintf('%s\n',NetPredicted{5});
    ModelBuilder = Analysis_FitModel_function2(Input_var,[],EnvironPara,NetPredicted{5});
    Net = NetPredicted(1:4);
    MonitorPredicted = reshape(MonitorPredicted,[N_Day*N_Site,1]);
    fprintf('prediction...model type:%s\n',ModelBuilder.NAME);
    [MonitorPredicted(ModelBuilder.Index_Day_1),~] = Analysis_NeuralNetwork([],[],Input_var(ModelBuilder.Index_Day_1,ModelBuilder.Index_Var),EnvironPara.N_Layer,Net,EnvironPara);
    MonitorPredicted = reshape(MonitorPredicted,[N_Day,N_Site]);
   
    
elseif(strcmp(EnvironPara.OPTION,'SimplePredictionAOD'))
   %% just for AOD prediction
    fprintf('%s\n',NetPredicted{5});
    ModelBuilder = Analysis_FitModel_function2(Input_var,[],EnvironPara,'NonRandGC');
    Net = NetPredicted(1:4);
    MonitorPredicted = reshape(MonitorPredicted,[N_Day*N_Site,1]);
    fprintf('prediction...model type:%s\n',ModelBuilder.NAME);
    [MonitorPredicted(ModelBuilder.Index_Day_1),~] = Analysis_NeuralNetwork([],[],Input_var(ModelBuilder.Index_Day_1,ModelBuilder.Index_Var),EnvironPara.N_Layer,Net,EnvironPara);
    MonitorPredicted = reshape(MonitorPredicted,[N_Day,N_Site]);
    
elseif(strcmp(EnvironPara.OPTION,'All_roundInf'))
%     while true
%         ModelBuilder = CalibrateGC_FitModel_function2_1(Input_var,EnvironPara);
%         [Input_var,MonitorPredicted,NetPredicted,Count] =...
%             Analysis_FitModel_function1(Input_var,MonitorPredicted,NetPredicted,ModelBuilder{6},EnvironPara,1,Inf,Count);
%     end
elseif(strcmp(EnvironPara.OPTION,'All_round2'))
    % for modeling PM2.5
    
    %%round 1
    EnvironPara.ModelName = 'Hybrid';
    ModelBuilder = Analysis_FitModel_function2(Input_var,[],EnvironPara,'General');
    [Input_var,MonitorPredicted,NetPredicted,Count] =...
        Analysis_FitModel_function1(Input_var,MonitorPredicted,NetPredicted,ModelBuilder,EnvironPara,1,Inf,Count);
    
    %round 2
    EnvironPara.ModelName = 'NeuralNetwork';
    ModelBuilder = Analysis_FitModel_function2(Input_var,[],EnvironPara,'General');
    [Input_var,MonitorPredicted,NetPredicted,Count] =...
        Analysis_FitModel_function1(Input_var,MonitorPredicted,NetPredicted,ModelBuilder,EnvironPara,2,Inf,Count);

elseif(strcmp(EnvironPara.OPTION,'All_round3'))
    %%round 1
    EnvironPara.ModelName = 'NeuralNetwork';
    ModelBuilder = Analysis_FitModel_function2(Input_var,[],EnvironPara,'General');
    [Input_var,MonitorPredicted,NetPredicted,Count] =...
        Analysis_FitModel_function1(Input_var,MonitorPredicted,NetPredicted,ModelBuilder,EnvironPara,1,Inf,Count);
    
    %round 2
    EnvironPara.ModelName = 'NeuralNetwork';
    ModelBuilder = Analysis_FitModel_function2(Input_var,[],EnvironPara,'General');
    [Input_var,MonitorPredicted,NetPredicted,Count] =...
        Analysis_FitModel_function1(Input_var,MonitorPredicted,NetPredicted,ModelBuilder,EnvironPara,2,Inf,Count);

elseif(strcmp(EnvironPara.OPTION,'All_round1'))
    %%round 1
    ModelBuilder = Analysis_FitModel_function2(Input_var,[],EnvironPara,'General');
    [Input_var,MonitorPredicted,NetPredicted,Count] =...
        Analysis_FitModel_function1(Input_var,MonitorPredicted,NetPredicted,ModelBuilder,EnvironPara,1,Inf,Count);

% elseif(strcmp(EnvironPara.OPTION,'All_hybrid'))
%     disp('warning: EnvironPara.ModelName option does not work now');
%     
%     EnvironPara.ModelName = 'Mixed_Effect';
%     ModelBuilder = Analysis_FitModel_function2(Input_var,[],EnvironPara,'General');
%     [Input_var,MonitorPredicted,NetPredicted,Count] =...
%         Analysis_FitModel_function1(Input_var,MonitorPredicted,NetPredicted,ModelBuilder,EnvironPara,1,Inf,Count);
%     
%     Input 
    
elseif(strcmp(EnvironPara.OPTION,'All_2'))
    %% for calibration AOD
%     %%round 0
%     ModelBuilder = Analysis_FitModel_function2(Input_var,[],EnvironPara,'RandGCOnly');
%     [Input_var,MonitorPredicted,NetPredicted,Count] =...
%         Analysis_FitModel_function1(Input_var,MonitorPredicted,NetPredicted,ModelBuilder,EnvironPara,1,Inf,Count);
%     Analysis_ResultDescriptiveAnalysis(MonitorPredicted',Input_var(:,1),NetPredicted,EnvironPara,[EnvironPara.CommonStartPoint,' ',num2str(Count)],{true,false,false,true});
    %%round 1
    ModelBuilder = Analysis_FitModel_function2(Input_var,[],EnvironPara,'RandOMI');
    [Input_var,MonitorPredicted,NetPredicted,Count] =...
        Analysis_FitModel_function1(Input_var,MonitorPredicted,NetPredicted,ModelBuilder,EnvironPara,1,Inf,Count);
    [R1,~,~]=Analysis_ResultDescriptiveAnalysis(MonitorPredicted',Input_var(:,1),NetPredicted,EnvironPara,[EnvironPara.CommonStartPoint,' ',num2str(Count)],{true,false,false,true});
    %%round 2
    if(strcmp(EnvironPara.IsAqua,'Terra'))
        ModelBuilder = Analysis_FitModel_function2(Input_var,[],EnvironPara,'RandAquaNearby');
    elseif(strcmp(EnvironPara.IsAqua,'Aqua'))
        ModelBuilder = Analysis_FitModel_function2(Input_var,[],EnvironPara,'RandTerraNearby');
    end
    [Input_var,MonitorPredicted,NetPredicted,Count] =...
        Analysis_FitModel_function1(Input_var,MonitorPredicted,NetPredicted,ModelBuilder,EnvironPara,1,Inf,Count);
    [R2,~,~]=Analysis_ResultDescriptiveAnalysis(MonitorPredicted',Input_var(:,1),NetPredicted,EnvironPara,[EnvironPara.CommonStartPoint,' ',num2str(Count)],{true,false,false,true});
    %%round 3
    ModelBuilder = Analysis_FitModel_function2(Input_var,[],EnvironPara,'RandGeneral');
    [Input_var,MonitorPredicted,NetPredicted,Count] =...
        Analysis_FitModel_function1(Input_var,MonitorPredicted,NetPredicted,ModelBuilder,EnvironPara,1,Inf,Count);
    [R3,~,~]=Analysis_ResultDescriptiveAnalysis(MonitorPredicted',Input_var(:,1),NetPredicted,EnvironPara,[EnvironPara.CommonStartPoint,' ',num2str(Count)],{true,false,false,true});
    
    NetPredicted{1,1}{6,1} = R1;
    NetPredicted{2,1}{6,1} = R2;
    NetPredicted{3,1}{6,1} = R3;
    
 elseif(strcmp(EnvironPara.OPTION,'Prediction_2'))
    %% prediction code for calibration AOD
    MonitorPredicted = reshape(MonitorPredicted,[N_Day*N_Site,1]);
    
    %%round 1
    ModelBuilder = Analysis_FitModel_function2(Input_var,[],EnvironPara,'RandGC');
    Net_temp = NetPredicted{1};
    Net = Net_temp(1:4);
    if(~isempty(Net{1}))
        fprintf('prediction...model type:%s\n',ModelBuilder.NAME);
        [MonitorPredicted(ModelBuilder.Index_Day_1),~] = Analysis_NeuralNetwork([],[],Input_var(ModelBuilder.Index_Day_1,ModelBuilder.Index_Var),EnvironPara.N_Layer,Net,EnvironPara);
    end
    
    %%round 2
    if(strcmp(EnvironPara.IsAqua,'Terra'))
        ModelBuilder = Analysis_FitModel_function2(Input_var,[],EnvironPara,'RandAquaNearby');
    elseif(strcmp(EnvironPara.IsAqua,'Aqua'))
        ModelBuilder = Analysis_FitModel_function2(Input_var,[],EnvironPara,'RandTerraNearby');
    end
    Net_temp = NetPredicted{2};
    Net = Net_temp(1:4);
    if(~isempty(Net{1}))
        fprintf('prediction...model type:%s\n',ModelBuilder.NAME);
        [MonitorPredicted(ModelBuilder.Index_Day_1),~] = Analysis_NeuralNetwork([],[],Input_var(ModelBuilder.Index_Day_1,ModelBuilder.Index_Var),EnvironPara.N_Layer,Net,EnvironPara);
    end

    %%round 3
    ModelBuilder = Analysis_FitModel_function2(Input_var,[],EnvironPara,'RandGeneral');
    Net_temp = NetPredicted{3};
    Net = Net_temp(1:4);
    if(~isempty(Net{1}))
        fprintf('prediction...model type:%s\n',ModelBuilder.NAME);
        [MonitorPredicted(ModelBuilder.Index_Day_1),~] = Analysis_NeuralNetwork([],[],Input_var(ModelBuilder.Index_Day_1,ModelBuilder.Index_Var),EnvironPara.N_Layer,Net,EnvironPara);
    end
    
    MonitorPredicted = reshape(MonitorPredicted,[N_Day,N_Site]);

elseif(strcmp(EnvironPara.OPTION,'All_3'))

end
    
    
