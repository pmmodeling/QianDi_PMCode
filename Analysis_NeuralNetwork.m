%% the core code to do neural network
% Y_train and other input: column standing for sites, so do output
% version: with some improvement: put train and test data together to do
% normalization; cancel normalization inside neural network
% source: Part3_Algorithm_Network --- almost the same, but enable
% multicore calculation
% version 2016-08-29 change name from CalibrateGC_NeuralNetwork to Analysis_NeuralNetwork
% change how to initiate the neural network

%% Parameters:
% Y_Train: target value
% X_Train: training variables;
% X_Test: testing data
% N_Layer: neural network structure;
% net_1: fitted neural network; only used for prediction;
% EnvironPara some marco definition

%% code
function [PredictedTest,net_1] = Analysis_NeuralNetwork(Y_Train,X_Train,X_Test,N_Layer,net_1,EnvironPara)

global TESTING MULTICORE
% TESTING = false;

% for testing purpose or not
if(true==TESTING)
    TrainEpoch = 1;
 else
    TrainEpoch = 500;
end


%% used fitted model to do prediction
if(strcmp(EnvironPara.OPTION,'Prediction'))% a pre-trained neural network is set as input, so we do not need to train the network. Just use it in simulation
    net = net_1{1};
    PS_Y = net_1{2};
    PS_X_1 = net_1{3};
    PS_X = net_1{4};
    %%% normalization and remove constant values
    %put data together
    %     [Y_Train,PS_Y]=mapminmax(Y_Train');% every row is a varible, each column stands for site
    %     [X_Train,PS_X_1] = removeconstantrows(X_Train');%remove constant rows;
    %     [X_Train,PS_X]=mapminmax(X_Train);% every row is a varible, each column stands for site
    X_Test = removeconstantrows('apply',X_Test',PS_X_1);
    X_Test = mapminmax('apply',X_Test,PS_X) ;
    
    %%% run the model!!!
    if(strcmp(EnvironPara.ModelName,'NeuralNetwork'))
        if(1 == MULTICORE)
            Temp2 = sim(net,X_Test,'useParallel','no');
        elseif(1 < MULTICORE)
            Temp2 = sim(net,X_Test,'useParallel','no');
%             Temp2 = sim(net,X_Test,'useParallel','yes','showResources','yes');
        end
    elseif(strcmp(EnvironPara.ModelName,'GradientBoosting'))
       
    %% hybrid model here means fitting mixed effect model first;
    % the predicted value from the mixed-effect model will be one input
    % variable of the neural network. This model strategy has been depreciated
    elseif(strcmp(EnvironPara.ModelName,'Hybrid'))
        % mixed-effect model
        lme = net{1};
        Neuralnet = net{2};
        % create data table
        TestData = array2table(X_Test');
        TestData.Properties.VariableNames = EnvironPara.VariableList(2:end)';       
        if(strcmp(EnvironPara.NAME,'MeanTemperature')||strcmp(EnvironPara.NAME,'MaxTemperature')||strcmp(EnvironPara.NAME,'MinTemperature'))
            TestData.AOD_PBL = TestData.MOD11A1_LST_Day_1km_Nearest4./TestData.REANALYSIS_hpbl_DailyMean;
            TestData.AOD_HM = TestData.MOD11A1_LST_Day_1km_Nearest4.*TestData.REANALYSIS_shum_2m_DailyMean;
        elseif(strcmp(EnvironPara.NAME,'PM25'))
            TestData.AOD_PBL = TestData.MOD04L2_550./TestData.REANALYSIS_hpbl_DailyMean;
            TestData.AOD_HM = TestData.MOD04L2_550.*TestData.REANALYSIS_shum_2m_DailyMean;        
        elseif(strcmp(EnvironPara.NAME,'NO2'))
            TestData.AOD_PBL = TestData.OMNO2d_ColumnAmountNO2StratoCloudScreened_Mean./TestData.REANALYSIS_hpbl_DailyMean;
            TestData.AOD_HM = TestData.OMNO2d_ColumnAmountNO2StratoCloudScreened_Mean.*TestData.REANALYSIS_shum_2m_DailyMean;            
        elseif(strcmp(EnvironPara.NAME,'Ozone'))
            TestData.AOD_PBL = TestData.OMTO3e_ColumnAmountO3./TestData.REANALYSIS_hpbl_DailyMean;
            TestData.AOD_HM = TestData.OMTO3e_ColumnAmountO3.*TestData.REANALYSIS_shum_2m_DailyMean;
        end
        
        % mixed-effect model
        disp('mixed-effect model');
        Test_Pred = predict(lme,TestData)';
        % neural network
        if(1 == MULTICORE)
            Temp2 = sim(Neuralnet,[Test_Pred;X_Test],'useParallel','no');
        elseif(1 < MULTICORE)
            Temp2 = sim(Neuralnet,[Test_Pred;X_Test],'useParallel','yes','showResources','yes');
        end        
    elseif(strcmp(EnvironPara.ModelName,'Lasso'))
        B = net{1};
        FitInfo = net{2};
        X_Test = X_Test';
        Temp2 = X_Test*B(:,FitInfo.Index1SE)+FitInfo.Intercept(FitInfo.Index1SE);
    end
    
    %%% reverse normalization
    PredictedTest = mapminmax('reverse',Temp2,PS_Y)';
    
%% doing model training
else
    
    if(1 < MULTICORE)
        myCluster = parcluster('local');
        myCluster.NumWorkers = MULTICORE;
        saveProfile(myCluster);
    end
    
    %%% normalization and remove constant variables
    %put data together
    [Y_Train,PS_Y]=mapminmax(Y_Train');% every row is a varible, each column stands for site
    %%% here, I assume X_train and X_test are of the same data source. We
    %%% do not need to put them together to remove the constant rows.
    [X_Train,PS_X_1] = removeconstantrows(X_Train');%remove constant rows;
    [X_Train,PS_X]=mapminmax(X_Train);% every row is a varible, each column stands for site
    X_Test = removeconstantrows('apply',X_Test',PS_X_1);
    X_Test = mapminmax('apply',X_Test,PS_X);
    % Y_Train: 1*114247; X_Train:105*114247; X_Test: 105*483950
    
    %train the network --- X_Train:22*226792; Y_Train:1*226792;X_Test:22*544580
    if(strcmp(EnvironPara.ModelName,'NeuralNetwork'))
        
        % fit neural network using training data and predict testing data
        net = feedforwardnet(N_Layer{1});
        for i = 1:(length(N_Layer{1,1})+1)
            net.layers{1}.transferFcn =  N_Layer{1,2}{i};
        end
        
        net.trainFcn = 'trainlm';
        net = configure(net,X_Train,Y_Train);
        net.trainParam.showCommandLine = true;
        net.trainParam.showWindow = false;
        net.trainParam.goal = 0;
        net.trainParam.min_grad = 1e-6;
        net.trainParam.epochs=TrainEpoch;
        net.trainParam.time = 10000*60;
        
        if(1 == MULTICORE)
            net = train(net,X_Train,Y_Train,'useParallel','no');
            Temp2 = sim(net,X_Test,'useParallel','no');
        elseif(1 < MULTICORE)
            poolobj = parpool('local',MULTICORE);
            net = train(net,X_Train,Y_Train,'useParallel','yes','showResources','yes');
            Temp2 = sim(net,X_Test,'useParallel','yes','showResources','yes');
            delete(poolobj);
        end
    %% hybrid model here means fitting mixed effect model first;
    % the predicted value from the mixed-effect model will be one input
    % variable of the neural network. This model strategy has been depreciated
    elseif(strcmp(EnvironPara.ModelName,'Hybrid'))
        % creat data table       
        TempVariableList = EnvironPara.VariableList(2:end);
        TempVariableList(PS_X_1.remove,:) = [];
        TrainData = array2table([Y_Train',X_Train']);
        TrainData.Properties.VariableNames = [EnvironPara.VariableList{1};TempVariableList]';
        TestData = array2table(X_Test');
        TestData.Properties.VariableNames = TempVariableList';
        
        %mixed-effect model
        disp('mixed-effect model');
        %mixed-effect model for temperature
        if(strcmp(EnvironPara.NAME,'MeanTemperature')||strcmp(EnvironPara.NAME,'MaxTemperature')||strcmp(EnvironPara.NAME,'MinTemperature'))
            TrainData.AOD_PBL = TrainData.MOD11A1_LST_Day_1km_Nearest4./TrainData.REANALYSIS_hpbl_DailyMean;
            TrainData.AOD_HM = TrainData.MOD11A1_LST_Day_1km_Nearest4.*TrainData.REANALYSIS_shum_2m_DailyMean;
            TestData.AOD_PBL = TestData.MOD11A1_LST_Day_1km_Nearest4./TestData.REANALYSIS_hpbl_DailyMean;
            TestData.AOD_HM = TestData.MOD11A1_LST_Day_1km_Nearest4.*TestData.REANALYSIS_shum_2m_DailyMean;
            lme = fitlme(TrainData,'MonitorData ~ MOD11A1_LST_Day_1km_Nearest4 + AOD_PBL +REANALYSIS_air_sfc_DailyMean +REANALYSIS_shum_2m_DailyMean + AOD_HM +(MOD11A1_LST_Day_1km_Nearest4|CalendarDay:PM25_Region)');
            Train_Pred = predict(lme,TrainData)';
            Test_Pred = predict(lme,TestData)';
            fprintf('Mixed-effect model Temperature DONE with R2:%f\n',corr(Train_Pred',TrainData.MonitorData,'rows','complete')^2);
            fprintf(EnvironPara.FID_SUMMARY,'Mixed-effect model Temperature DONE with R2:%f\n',corr(Train_Pred',TrainData.MonitorData,'rows','complete')^2);
        %mixed-effect model for PM2.5, there are five mixed effect model,
        %one for every AOD variable
        elseif(strcmp(EnvironPara.NAME,'PM25'))
            TrainData.AODPBL_MAIACUS_Optical_Depth_055_Aqua_Nearest4 = TrainData.MAIACUS_Optical_Depth_055_Aqua_Nearest4./TrainData.REANALYSIS_hpbl_DailyMean;
            TrainData.AODHM_MAIACUS_Optical_Depth_055_Aqua_Nearest4 = TrainData.MAIACUS_Optical_Depth_055_Aqua_Nearest4.*TrainData.REANALYSIS_shum_2m_DailyMean;
            TestData.AODPBL_MAIACUS_Optical_Depth_055_Aqua_Nearest4 = TestData.MAIACUS_Optical_Depth_055_Aqua_Nearest4./TestData.REANALYSIS_hpbl_DailyMean;
            TestData.AODHM_MAIACUS_Optical_Depth_055_Aqua_Nearest4 = TestData.MAIACUS_Optical_Depth_055_Aqua_Nearest4.*TestData.REANALYSIS_shum_2m_DailyMean;            
            
            TrainData.AODPBL_MAIACUS_Optical_Depth_047_Aqua_Nearest4 = TrainData.MAIACUS_Optical_Depth_047_Aqua_Nearest4./TrainData.REANALYSIS_hpbl_DailyMean;
            TrainData.AODHM_MAIACUS_Optical_Depth_047_Aqua_Nearest4 = TrainData.MAIACUS_Optical_Depth_047_Aqua_Nearest4.*TrainData.REANALYSIS_shum_2m_DailyMean;
            TestData.AODPBL_MAIACUS_Optical_Depth_047_Aqua_Nearest4 = TestData.MAIACUS_Optical_Depth_047_Aqua_Nearest4./TestData.REANALYSIS_hpbl_DailyMean;
            TestData.AODHM_MAIACUS_Optical_Depth_047_Aqua_Nearest4 = TestData.MAIACUS_Optical_Depth_047_Aqua_Nearest4.*TestData.REANALYSIS_shum_2m_DailyMean;            

            TrainData.AODPBL_MAIACUS_Optical_Depth_055_Terra_Nearest4 = TrainData.MAIACUS_Optical_Depth_055_Terra_Nearest4./TrainData.REANALYSIS_hpbl_DailyMean;
            TrainData.AODHM_MAIACUS_Optical_Depth_055_Terra_Nearest4 = TrainData.MAIACUS_Optical_Depth_055_Terra_Nearest4.*TrainData.REANALYSIS_shum_2m_DailyMean;
            TestData.AODPBL_MAIACUS_Optical_Depth_055_Terra_Nearest4 = TestData.MAIACUS_Optical_Depth_055_Terra_Nearest4./TestData.REANALYSIS_hpbl_DailyMean;
            TestData.AODHM_MAIACUS_Optical_Depth_055_Terra_Nearest4 = TestData.MAIACUS_Optical_Depth_055_Terra_Nearest4.*TestData.REANALYSIS_shum_2m_DailyMean;            
            
            TrainData.AODPBL_MAIACUS_Optical_Depth_047_Terra_Nearest4 = TrainData.MAIACUS_Optical_Depth_047_Terra_Nearest4./TrainData.REANALYSIS_hpbl_DailyMean;
            TrainData.AODHM_MAIACUS_Optical_Depth_047_Terra_Nearest4 = TrainData.MAIACUS_Optical_Depth_047_Terra_Nearest4.*TrainData.REANALYSIS_shum_2m_DailyMean;
            TestData.AODPBL_MAIACUS_Optical_Depth_047_Terra_Nearest4 = TestData.MAIACUS_Optical_Depth_047_Terra_Nearest4./TestData.REANALYSIS_hpbl_DailyMean;
            TestData.AODHM_MAIACUS_Optical_Depth_047_Terra_Nearest4 = TestData.MAIACUS_Optical_Depth_047_Terra_Nearest4.*TestData.REANALYSIS_shum_2m_DailyMean;            
            
            TrainData.AODPBL_MOD04L2_550 = TrainData.MOD04L2_550./TrainData.REANALYSIS_hpbl_DailyMean;
            TrainData.AODHM_MOD04L2_550 = TrainData.MOD04L2_550.*TrainData.REANALYSIS_shum_2m_DailyMean;
            TestData.AODPBL_MOD04L2_550 = TestData.MOD04L2_550./TestData.REANALYSIS_hpbl_DailyMean;
            TestData.AODHM_MOD04L2_550 = TestData.MOD04L2_550.*TestData.REANALYSIS_shum_2m_DailyMean;            
            
            
            Train_Pred = nan(5,size(TrainData,1));
            Test_Pred = nan(5,size(TestData,1));
            lme = cell(5,1);
            parfor k = 1:5 % mixed-effect model, on all five AOD variables (not used any more)
                if(1 == k)
                    % MAIACUS_Optical_Depth_055_Aqua_Nearest4
                    tic;
                    lme_temp = fitlme(TrainData,'MonitorData ~ MAIACUS_Optical_Depth_055_Aqua_Nearest4 +AODPBL_MAIACUS_Optical_Depth_055_Aqua_Nearest4+REANALYSIS_air_sfc_DailyMean +REANALYSIS_shum_2m_DailyMean +REANALYSIS_windspeed_10m_DailyMean +AODHM_MAIACUS_Optical_Depth_055_Aqua_Nearest4+(REANALYSIS_air_sfc_DailyMean+MAIACUS_Optical_Depth_055_Aqua_Nearest4|CalendarDay:PM25_Region)');            
                    Train_Pred(k,:) = predict(lme_temp,TrainData)';
                    Test_Pred(k,:) = predict(lme_temp,TestData)';
                    lme{k} = lme_temp;
                    toc;       
                    fprintf('Mixed-effect model MAIACUS_Optical_Depth_055_Aqua_Nearest4 DONE with R2:%f\n',corr(Train_Pred(k,:)',TrainData.MonitorData,'rows','complete')^2);
                    
                elseif(2 == k)
                    % MAIACUS_Optical_Depth_047_Aqua_Nearest4
                    tic;
                    lme_temp = fitlme(TrainData,'MonitorData ~ MAIACUS_Optical_Depth_047_Aqua_Nearest4 +AODPBL_MAIACUS_Optical_Depth_047_Aqua_Nearest4+REANALYSIS_air_sfc_DailyMean +REANALYSIS_shum_2m_DailyMean +REANALYSIS_windspeed_10m_DailyMean +AODHM_MAIACUS_Optical_Depth_047_Aqua_Nearest4+(REANALYSIS_air_sfc_DailyMean+MAIACUS_Optical_Depth_047_Aqua_Nearest4|CalendarDay:PM25_Region)');            
                    Train_Pred(k,:) = predict(lme_temp,TrainData)';
                    Test_Pred(k,:) = predict(lme_temp,TestData)';
                    lme{k} = lme_temp;
                    toc;
                    fprintf('Mixed-effect model MAIACUS_Optical_Depth_047_Aqua_Nearest4 DONE with R2:%f\n',corr(Train_Pred(k,:)',TrainData.MonitorData,'rows','complete')^2);
                    
                elseif(3 == k)
                    % MAIACUS_Optical_Depth_055_Terra_Nearest4
                    tic;
                    lme_temp = fitlme(TrainData,'MonitorData ~ MAIACUS_Optical_Depth_055_Terra_Nearest4 +AODPBL_MAIACUS_Optical_Depth_055_Terra_Nearest4+REANALYSIS_air_sfc_DailyMean +REANALYSIS_shum_2m_DailyMean +REANALYSIS_windspeed_10m_DailyMean +AODHM_MAIACUS_Optical_Depth_055_Terra_Nearest4+(REANALYSIS_air_sfc_DailyMean+MAIACUS_Optical_Depth_055_Terra_Nearest4|CalendarDay:PM25_Region)');            
                    Train_Pred(k,:) = predict(lme_temp,TrainData)';
                    Test_Pred(k,:) = predict(lme_temp,TestData)';
                    lme{k} = lme_temp;
                    toc;
                    fprintf('Mixed-effect model MAIACUS_Optical_Depth_055_Terra_Nearest4 DONE with R2:%f\n',corr(Train_Pred(k,:)',TrainData.MonitorData,'rows','complete')^2);
                    
                elseif(4 == k)
                    % MAIACUS_Optical_Depth_047_Terra_Nearest4
                    tic;
                    lme_temp = fitlme(TrainData,'MonitorData ~ MAIACUS_Optical_Depth_047_Terra_Nearest4 +AODPBL_MAIACUS_Optical_Depth_047_Terra_Nearest4+REANALYSIS_air_sfc_DailyMean +REANALYSIS_shum_2m_DailyMean +REANALYSIS_windspeed_10m_DailyMean +AODHM_MAIACUS_Optical_Depth_047_Terra_Nearest4+(REANALYSIS_air_sfc_DailyMean+MAIACUS_Optical_Depth_047_Terra_Nearest4|CalendarDay:PM25_Region)');            
                    Train_Pred(k,:) = predict(lme_temp,TrainData)';
                    Test_Pred(k,:) = predict(lme_temp,TestData)';
                    lme{k} = lme_temp;
                    toc;
                    fprintf('Mixed-effect model MAIACUS_Optical_Depth_047_Terra_Nearest4 DONE with R2:%f\n',corr(Train_Pred(k,:)',TrainData.MonitorData,'rows','complete')^2);
                elseif(5 == k)
                    % MOD04L2_550
                    tic;
                    lme_temp = fitlme(TrainData,'MonitorData ~ MOD04L2_550 +AODPBL_MOD04L2_550+REANALYSIS_air_sfc_DailyMean +REANALYSIS_shum_2m_DailyMean +REANALYSIS_windspeed_10m_DailyMean +AODHM_MOD04L2_550+(REANALYSIS_air_sfc_DailyMean+MOD04L2_550|CalendarDay:PM25_Region)');            
                    Train_Pred(k,:) = predict(lme_temp,TrainData)';
                    Test_Pred(k,:) = predict(lme_temp,TestData)';
                    lme{k} = lme_temp;
                    toc;
                    fprintf('Mixed-effect model MOD04L2_550 DONE with R2:%f\n',corr(Train_Pred(k,:)',TrainData.MonitorData,'rows','complete')^2);
                end
            end
            poolobj = gcp('nocreate');
            delete(poolobj);
            
            %% report results
            fprintf(EnvironPara.FID_SUMMARY,'Mixed-effect model MAIACUS_Optical_Depth_055_Aqua_Nearest4 DONE with R2:%f\n',corr(Train_Pred(1,:)',TrainData.MonitorData,'rows','complete')^2);
            fprintf(EnvironPara.FID_SUMMARY,'Mixed-effect model MAIACUS_Optical_Depth_047_Aqua_Nearest4 DONE with R2:%f\n',corr(Train_Pred(2,:)',TrainData.MonitorData,'rows','complete')^2);
            fprintf(EnvironPara.FID_SUMMARY,'Mixed-effect model MAIACUS_Optical_Depth_055_Terra_Nearest4 DONE with R2:%f\n',corr(Train_Pred(3,:)',TrainData.MonitorData,'rows','complete')^2);
            fprintf(EnvironPara.FID_SUMMARY,'Mixed-effect model MAIACUS_Optical_Depth_047_Terra_Nearest4 DONE with R2:%f\n',corr(Train_Pred(4,:)',TrainData.MonitorData,'rows','complete')^2);
            fprintf(EnvironPara.FID_SUMMARY,'Mixed-effect model MOD04L2_550 DONE with R2:%f\n',corr(Train_Pred(5,:)',TrainData.MonitorData,'rows','complete')^2);          
            
        %mixed-effect model for NO2
        elseif(strcmp(EnvironPara.NAME,'NO2'))
            TrainData.AOD_PBL = TrainData.OMNO2d_ColumnAmountNO2StratoCloudScreened_Mean./TrainData.REANALYSIS_hpbl_DailyMean;
            TrainData.AOD_HM = TrainData.OMNO2d_ColumnAmountNO2StratoCloudScreened_Mean.*TrainData.REANALYSIS_shum_2m_DailyMean;
            TestData.AOD_PBL = TestData.OMNO2d_ColumnAmountNO2StratoCloudScreened_Mean./TestData.REANALYSIS_hpbl_DailyMean;
            TestData.AOD_HM = TestData.OMNO2d_ColumnAmountNO2StratoCloudScreened_Mean.*TestData.REANALYSIS_shum_2m_DailyMean;
            lme = fitlme(TrainData,'MonitorData ~ OMNO2d_ColumnAmountNO2StratoCloudScreened_Mean + AOD_PBL +REANALYSIS_air_sfc_DailyMean +REANALYSIS_shum_2m_DailyMean +REANALYSIS_windspeed_10m_DailyMean + AOD_HM +(REANALYSIS_air_sfc_DailyMean+OMNO2d_ColumnAmountNO2StratoCloudScreened_Mean|CalendarDay:NO2_Region)');
            Train_Pred = predict(lme,TrainData)';
            Test_Pred = predict(lme,TestData)';
            fprintf('Mixed-effect model NO2 DONE with R2:%f\n',corr(Train_Pred',TrainData.MonitorData,'rows','complete')^2);
            fprintf(EnvironPara.FID_SUMMARY,'Mixed-effect model NO2 DONE with R2:%f\n',corr(Train_Pred',TrainData.MonitorData,'rows','complete')^2);
            
        %mixed-effect model for ozone
        elseif(strcmp(EnvironPara.NAME,'Ozone'))
            TrainData.AOD_PBL = TrainData.OMTO3e_ColumnAmountO3./TrainData.REANALYSIS_hpbl_DailyMean;
            TrainData.AOD_HM = TrainData.OMTO3e_ColumnAmountO3.*TrainData.REANALYSIS_shum_2m_DailyMean;
            TestData.AOD_PBL = TestData.OMTO3e_ColumnAmountO3./TestData.REANALYSIS_hpbl_DailyMean;
            TestData.AOD_HM = TestData.OMTO3e_ColumnAmountO3.*TestData.REANALYSIS_shum_2m_DailyMean;
            lme = fitlme(TrainData,'MonitorData ~ OMTO3e_ColumnAmountO3 + AOD_PBL +REANALYSIS_air_sfc_DailyMean +REANALYSIS_shum_2m_DailyMean +REANALYSIS_windspeed_10m_DailyMean + AOD_HM +(REANALYSIS_air_sfc_DailyMean+OMTO3e_ColumnAmountO3|CalendarDay:Ozone_Region)');
            Train_Pred = predict(lme,TrainData)';
            Test_Pred = predict(lme,TestData)';
            fprintf('Mixed-effect model Ozone DONE with R2:%f\n',corr(Train_Pred',TrainData.MonitorData,'rows','complete')^2);
            fprintf(EnvironPara.FID_SUMMARY,'Mixed-effect model Ozone DONE with R2:%f\n',corr(Train_Pred',TrainData.MonitorData,'rows','complete')^2);
        
        end

        
        %neural network
        Neuralnet = feedforwardnet(N_Layer{1});
        for i = 1:(length(N_Layer{1,1})+1)
            Neuralnet.layers{1}.transferFcn =  N_Layer{1,2}{i};
        end
        
        %%%%%%%%%%%%%%%%%
        %%% load csv results from here
        %%%%%%%%%%%%%%%%%%%
        
        Neuralnet.trainFcn = 'trainlm';
        Neuralnet = configure(Neuralnet,[Train_Pred;X_Train],Y_Train);
        Neuralnet.trainParam.showCommandLine = true;
        Neuralnet.trainParam.showWindow = false;
        Neuralnet.trainParam.goal = 0;
        Neuralnet.trainParam.min_grad = 1e-6;
        Neuralnet.trainParam.epochs=TrainEpoch;
        Neuralnet.trainParam.time = 10000*60;
        
        %% train neural network
        if(1 == MULTICORE)
            Neuralnet = train(Neuralnet,[Train_Pred;X_Train],Y_Train,'useParallel','no');
            Temp2 = sim(Neuralnet,[Test_Pred;X_Test],'useParallel','no');
        elseif(1 < MULTICORE)
            poolobj = parpool('local',MULTICORE);
            Neuralnet = train(Neuralnet,[Train_Pred;X_Train],Y_Train,'useParallel','yes','showResources','yes');
            Temp2 = sim(Neuralnet,[Test_Pred;X_Test],'useParallel','yes','showResources','yes');
            delete(poolobj);
        end
        net = {lme;Neuralnet};
    elseif(strcmp(EnvironPara.ModelName,'GradientBoosting'))
        X_Train = X_Train';
        Y_Train = Y_Train';
        X_Test = X_Test';
        
        t = RegressionTree.template('MinLeaf',5,'Surrogate','all');
        net = fitensemble(X_Train,Y_Train,'LSBoost',500,t,'LearnRate',0.01,'KFold',10);
        A = kfoldLoss(net,'Mode','cumulative','lossfun','MSE');
        
    elseif(strcmp(EnvironPara.ModelName,'Lasso'))
        X_Train = X_Train';
        Y_Train = Y_Train';
        X_Test = X_Test';
        if(1 == MULTICORE)
            [B,FitInfo] = lasso(X_Train,Y_Train,'CV',10);
        elseif(1 < MULTICORE)
            poolobj = parpool('local',MULTICORE);
            opts = statset('UseParallel',true);
            [B,FitInfo] = lasso(X_Train,Y_Train,'CV',10,'Options',opts);
            delete(poolobj);
        end
        Temp2 = X_Test*B(:,FitInfo.Index1SE)+FitInfo.Intercept(FitInfo.Index1SE);
        net = {B,FitInfo};
    end
    % Temp2: 1*483950
    PredictedTest = mapminmax('reverse',Temp2,PS_Y)';
    % PredictedTest:483950*1
    net_1 = {net;PS_Y;PS_X_1;PS_X};
end

end

