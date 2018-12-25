%% data imputation
% EnvironPara.2016-11-26: at prediction stage, "compete" does not do anything
% 2017-05-31: use more dedicate models to fill in missingness;
function Input_var_AOD = Analysis_ImputationInputData(Input_var_AOD,EnvironPara,ImputationOption)
% for testing only
% load('Analysis_ImputationInputData_debug.mat')

Data = array2table(Input_var_AOD);
Data.Properties.VariableNames = EnvironPara.VariableList(1:end)';       
Data.REANALYSIS_windspeed_10m_DailyMean_1 = 1./Data.REANALYSIS_windspeed_10m_DailyMean;
% with missing value here, the model fails
if(ismember('PM25_Region',Data.Properties.VariableNames))
    Data.PM25_Region(isnan(Data.PM25_Region)) = 1;
end

if(ismember('NO2_Region',Data.Properties.VariableNames))
    Data.NO2_Region(isnan(Data.NO2_Region)) = 1;
end

if(ismember('Ozone_Region',Data.Properties.VariableNames))
    Data.Ozone_Region(isnan(Data.Ozone_Region)) = 1;
end

%% get the fillin model
mkdir([EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep]);
% TempFileName = [EnvironPara.DIRPATH_TRAIN,'Temp',EnvironPara.Sep,num2str(EnvironPara.IDNUM),EnvironPara.Sep,'',datestr(datenum(EnvironPara.StartYear,1,1),'yyyymmdd'),'_',datestr(datenum(EnvironPara.EndYear,12,31),'yyyymmdd'),'.mat'];
% if(exist(TempFileName,'file'))
%     LoadData_function(TempFileName);
% else


CommonFormula_weight = 'REANALYSIS_air_sfc_DailyMean+REANALYSIS_shum_2m_DailyMean+REANALYSIS_snowc_DailyMean+REANALYSIS_tcdc_DailyMean+USElevation_mea100';
% CommonFormula = 'Nearby_Peak2_NO2+Nearby_Peak2_PM25+Nearby_Peak2_Ozone+REANALYSIS_air_sfc_DailyMean+REANALYSIS_apcp_DailyMean+REANALYSIS_dlwrf_DailyMean+REANALYSIS_evap_DailyMean+REANALYSIS_hpbl_DailyMean+REANALYSIS_lhtfl_DailyMean+REANALYSIS_shtfl_DailyMean+REANALYSIS_shum_2m_DailyMean+REANALYSIS_snowc_DailyMean+REANALYSIS_tcdc_DailyMean+REANALYSIS_omega_DailyMean+REANALYSIS_windspeed_10m_DailyMean_1+REANALYSIS_weasd_DailyMean+REANALYSIS_prate_DailyMean+REANALYSIS_vis_DailyMean+NLCD_Barren100+NLCD_Developed100+NLCD_Herbaceous100+NLCD_Planted100+NLCD_Shrubland100+NLCD_Water100+NLCD_Wetlands100+RoadDensity_primaryroads1000+RoadDensity_prisecroads1000+RoadDensity_roads1000+USElevation_mea100+(1|CalendarDay)+(1|PM25_Region)';
if(strcmp(EnvironPara.ImputationOption_Model,'Model2'))
    CommonFormula = 'Nearby_Peak2_PM25+Nearby_Peak2_NO2+Nearby_Peak2_Ozone+REANALYSIS_air_sfc_DailyMean+REANALYSIS_evap_DailyMean+REANALYSIS_hpbl_DailyMean+REANALYSIS_lhtfl_DailyMean+REANALYSIS_shtfl_DailyMean+REANALYSIS_shum_2m_DailyMean+REANALYSIS_snowc_DailyMean+REANALYSIS_tcdc_DailyMean+REANALYSIS_omega_DailyMean+REANALYSIS_windspeed_10m_DailyMean_1+REANALYSIS_weasd_DailyMean+REANALYSIS_prate_DailyMean+REANALYSIS_vis_DailyMean+NLCD_Barren100+NLCD_Developed100+NLCD_Herbaceous100+NLCD_Planted100+NLCD_Shrubland100+NLCD_Water100+NLCD_Wetlands100+RoadDensity_primaryroads1000+RoadDensity_prisecroads1000+RoadDensity_roads1000+USElevation_mea100+GFEDFireCarbon+CMAQ_NO2+CMAQ_NO2_Vertical+CMAQ_Ozone+CMAQ_Ozone_Vertical+CMAQ_PM25_TOT+CMAQ_PM25_Vertical+CMAQ_PM25_NO3+CMAQ_PM25_SO4+(1|CalendarDay)+(1|PM25_Region)';
    CurrentVariableList = {'Nearby_Peak2_PM25','Nearby_Peak2_NO2','Nearby_Peak2_Ozone','REANALYSIS_air_sfc_DailyMean','REANALYSIS_evap_DailyMean','REANALYSIS_hpbl_DailyMean','REANALYSIS_lhtfl_DailyMean','REANALYSIS_shtfl_DailyMean','REANALYSIS_shum_2m_DailyMean','REANALYSIS_snowc_DailyMean','REANALYSIS_tcdc_DailyMean','REANALYSIS_omega_DailyMean','REANALYSIS_windspeed_10m_DailyMean_1','REANALYSIS_weasd_DailyMean','REANALYSIS_prate_DailyMean','REANALYSIS_vis_DailyMean','NLCD_Barren100','NLCD_Developed100','NLCD_Herbaceous100','NLCD_Planted100','NLCD_Shrubland100','NLCD_Water100','NLCD_Wetlands100','RoadDensity_primaryroads1000','RoadDensity_prisecroads1000','RoadDensity_roads1000','USElevation_mea100','GFEDFireCarbon','CMAQ_NO2','CMAQ_NO2_Vertical','CMAQ_Ozone','CMAQ_Ozone_Vertical','CMAQ_PM25_TOT','CMAQ_PM25_Vertical','CMAQ_PM25_NO3','CMAQ_PM25_SO4','CalendarDay','PM25_Region'};
    StartYear = 2000;
    EndYear = 2016;
elseif(strcmp(EnvironPara.ImputationOption_Model,'Model1'))
    CommonFormula = 'Nearby_Peak2_PM25+Nearby_Peak2_NO2+Nearby_Peak2_Ozone+REANALYSIS_air_sfc_DailyMean+REANALYSIS_evap_DailyMean+REANALYSIS_hpbl_DailyMean+REANALYSIS_lhtfl_DailyMean+REANALYSIS_shtfl_DailyMean+REANALYSIS_shum_2m_DailyMean+REANALYSIS_snowc_DailyMean+REANALYSIS_tcdc_DailyMean+REANALYSIS_omega_DailyMean+REANALYSIS_windspeed_10m_DailyMean_1+REANALYSIS_weasd_DailyMean+REANALYSIS_prate_DailyMean+REANALYSIS_vis_DailyMean+NLCD_Barren100+NLCD_Developed100+NLCD_Herbaceous100+NLCD_Planted100+NLCD_Shrubland100+NLCD_Water100+NLCD_Wetlands100+RoadDensity_primaryroads1000+RoadDensity_prisecroads1000+RoadDensity_roads1000+USElevation_mea100+(1|CalendarDay)+(1|PM25_Region)';
    CurrentVariableList = {'Nearby_Peak2_PM25','Nearby_Peak2_NO2','Nearby_Peak2_Ozone','REANALYSIS_air_sfc_DailyMean','REANALYSIS_evap_DailyMean','REANALYSIS_hpbl_DailyMean','REANALYSIS_lhtfl_DailyMean','REANALYSIS_shtfl_DailyMean','REANALYSIS_shum_2m_DailyMean','REANALYSIS_snowc_DailyMean','REANALYSIS_tcdc_DailyMean','REANALYSIS_omega_DailyMean','REANALYSIS_windspeed_10m_DailyMean_1','REANALYSIS_weasd_DailyMean','REANALYSIS_prate_DailyMean','REANALYSIS_vis_DailyMean','NLCD_Barren100','NLCD_Developed100','NLCD_Herbaceous100','NLCD_Planted100','NLCD_Shrubland100','NLCD_Water100','NLCD_Wetlands100','RoadDensity_primaryroads1000','RoadDensity_prisecroads1000','RoadDensity_roads1000','USElevation_mea100','CalendarDay','PM25_Region'};
    StartYear = 2000;
    EndYear = 2016;
elseif(strcmp(EnvironPara.ImputationOption_Model,'Model3'))
    CommonFormula = 'Nearby_Peak2_MaxTemperature+Nearby_Peak2_MinTemperature+Nearby_Peak2_MeanTemperature+REANALYSIS_air_sfc_DailyMean+REANALYSIS_evap_DailyMean+REANALYSIS_hpbl_DailyMean+REANALYSIS_lhtfl_DailyMean+REANALYSIS_shtfl_DailyMean+REANALYSIS_shum_2m_DailyMean+REANALYSIS_snowc_DailyMean+REANALYSIS_tcdc_DailyMean+REANALYSIS_omega_DailyMean+REANALYSIS_windspeed_10m_DailyMean_1+REANALYSIS_weasd_DailyMean+REANALYSIS_prate_DailyMean+REANALYSIS_vis_DailyMean+NLCD_Barren100+NLCD_Developed100+NLCD_Herbaceous100+NLCD_Planted100+NLCD_Shrubland100+NLCD_Water100+NLCD_Wetlands100+RoadDensity_primaryroads1000+RoadDensity_prisecroads1000+RoadDensity_roads1000+USElevation_mea100+(1|CalendarDay)+(1|PM25_Region)';
    CurrentVariableList = {'Nearby_Peak2_MaxTemperature','Nearby_Peak2_MinTemperature','Nearby_Peak2_MeanTemperature','REANALYSIS_air_sfc_DailyMean','REANALYSIS_evap_DailyMean','REANALYSIS_hpbl_DailyMean','REANALYSIS_lhtfl_DailyMean','REANALYSIS_shtfl_DailyMean','REANALYSIS_shum_2m_DailyMean','REANALYSIS_snowc_DailyMean','REANALYSIS_tcdc_DailyMean','REANALYSIS_omega_DailyMean','REANALYSIS_windspeed_10m_DailyMean_1','REANALYSIS_weasd_DailyMean','REANALYSIS_prate_DailyMean','REANALYSIS_vis_DailyMean','NLCD_Barren100','NLCD_Developed100','NLCD_Herbaceous100','NLCD_Planted100','NLCD_Shrubland100','NLCD_Water100','NLCD_Wetlands100','RoadDensity_primaryroads1000','RoadDensity_prisecroads1000','RoadDensity_roads1000','USElevation_mea100','CalendarDay','PM25_Region'};
    StartYear = EnvironPara.StartYear;
    EndYear = EnvironPara.EndYear;
elseif(strcmp(EnvironPara.ImputationOption_Model,'Model4'))
    CommonFormula = 'Nearby_Peak2_PM25+Nearby_Peak2_NO2+Nearby_Peak2_Ozone+REANALYSIS_air_sfc_DailyMean+REANALYSIS_evap_DailyMean+REANALYSIS_hpbl_DailyMean+REANALYSIS_lhtfl_DailyMean+REANALYSIS_shtfl_DailyMean+REANALYSIS_shum_2m_DailyMean+REANALYSIS_snowc_DailyMean+REANALYSIS_tcdc_DailyMean+REANALYSIS_omega_DailyMean+REANALYSIS_windspeed_10m_DailyMean_1+REANALYSIS_weasd_DailyMean+REANALYSIS_prate_DailyMean+REANALYSIS_vis_DailyMean+NLCD_Barren100+NLCD_Developed100+NLCD_Herbaceous100+NLCD_Planted100+NLCD_Shrubland100+NLCD_Water100+NLCD_Wetlands100+RoadDensity_primaryroads1000+RoadDensity_prisecroads1000+RoadDensity_roads1000+USElevation_mea100+GFEDFireCarbon+GEOSChem_PM25+GEOSChem_SO4+GEOSChem_NH4+GEOSChem_DUST+GEOSChem_BC+(1|CalendarDay)+(1|PM25_Region)';
    CurrentVariableList = {'Nearby_Peak2_PM25','Nearby_Peak2_NO2','Nearby_Peak2_Ozone','REANALYSIS_air_sfc_DailyMean','REANALYSIS_evap_DailyMean','REANALYSIS_hpbl_DailyMean','REANALYSIS_lhtfl_DailyMean','REANALYSIS_shtfl_DailyMean','REANALYSIS_shum_2m_DailyMean','REANALYSIS_snowc_DailyMean','REANALYSIS_tcdc_DailyMean','REANALYSIS_omega_DailyMean','REANALYSIS_windspeed_10m_DailyMean_1','REANALYSIS_weasd_DailyMean','REANALYSIS_prate_DailyMean','REANALYSIS_vis_DailyMean','NLCD_Barren100','NLCD_Developed100','NLCD_Herbaceous100','NLCD_Planted100','NLCD_Shrubland100','NLCD_Water100','NLCD_Wetlands100','RoadDensity_primaryroads1000','RoadDensity_prisecroads1000','RoadDensity_roads1000','USElevation_mea100','GFEDFireCarbon','CMAQ_NO2','CMAQ_NO2_Vertical','CMAQ_Ozone','CMAQ_Ozone_Vertical','CMAQ_PM25_TOT','CMAQ_PM25_Vertical','CMAQ_PM25_NO3','CMAQ_PM25_SO4','CalendarDay','PM25_Region'};
    StartYear = EnvironPara.StartYear;
    EndYear = EnvironPara.EndYear;
end



% lme = fitlme(Data,'Nearby_Peak2_PM25~Nearby_Peak2Lag1_PM25+Nearby_Peak2Lag3_PM25');
% Pred = predict(lme,Data);
% Index = isnan(Data.Nearby_Peak2_PM25);
% Data.Nearby_Peak2_PM25(Index) = Pred(Index);



if(~isempty(strfind(EnvironPara.OPTION,'CV'))||~isempty(strfind(EnvironPara.OPTION,'All')))

    %REANALYSIS_soilm_DailyMean
    if(ismember('REANALYSIS_soilm_DailyMean',Data.Properties.VariableNames))
        TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'REANALYSIS_soilm_DailyMean_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
        if(exist(TempFileName,'file'))
            %LoadData_function(TempFileName);
        else
            fprintf('imputation fitting: %s\n',['REANALYSIS_soilm_DailyMean ~',CommonFormula]);
            Data.Is_nan = isnan(Data.REANALYSIS_soilm_DailyMean);
            lme_weight = fitglm(Data,['Is_nan~',CommonFormula_weight],'Distribution','Binomial','Link','logit'); 
            weights = predict(lme_weight,Data);
            weights(isnan(weights)) = 0;
            lme = fitlme(Data(:,['REANALYSIS_soilm_DailyMean',CurrentVariableList]),['REANALYSIS_soilm_DailyMean~',CommonFormula],'Weights',1./(1-weights));
            [~,~,lme_rand] = randomEffects(lme);
            lme_coef = lme.Coefficients;
            save(TempFileName,'lme_rand','lme_coef','-v7.3');
            clear lme_weight lme 
        end
    end
    
    %REANALYSIS_gflux_DailyMean
    if(ismember('REANALYSIS_gflux_DailyMean',Data.Properties.VariableNames))
        TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'REANALYSIS_gflux_DailyMean_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
        if(exist(TempFileName,'file'))
            %LoadData_function(TempFileName);
        else
            fprintf('imputation fitting: %s\n',['REANALYSIS_gflux_DailyMean ~',CommonFormula]);
            Data.Is_nan = isnan(Data.REANALYSIS_gflux_DailyMean);
            lme_weight = fitglm(Data,['Is_nan~',CommonFormula_weight],'Distribution','Binomial','Link','logit'); 
            weights = predict(lme_weight,Data);
            weights(isnan(weights)) = 0;
            lme = fitlme(Data(:,['REANALYSIS_gflux_DailyMean',CurrentVariableList]),['REANALYSIS_gflux_DailyMean~',CommonFormula],'Weights',1./(1-weights));
            [~,~,lme_rand] = randomEffects(lme);
            lme_coef = lme.Coefficients;
            save(TempFileName,'lme_rand','lme_coef','-v7.3');
            clear lme_weight lme 
        end
    end
    
    %MAIACUS_Optical_Depth_047_Aqua_Nearest4
    if(ismember('MAIACUS_Optical_Depth_047_Aqua_Nearest4',Data.Properties.VariableNames))
        TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'MAIACUS_Optical_Depth_047_Aqua_Nearest4_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
        if(exist(TempFileName,'file'))
            %LoadData_function(TempFileName);
        else
            fprintf('imputation fitting: %s\n',['MAIACUS_Optical_Depth_047_Aqua_Nearest4 ~',CommonFormula]);
            Data.Is_nan = isnan(Data.MAIACUS_Optical_Depth_047_Aqua_Nearest4);
            lme_weight = fitglm(Data,['Is_nan~',CommonFormula_weight],'Distribution','Binomial','Link','logit'); 
            weights = predict(lme_weight,Data);
            weights(isnan(weights)) = 0;
            lme = fitlme(Data(:,['MAIACUS_Optical_Depth_047_Aqua_Nearest4',CurrentVariableList]),['MAIACUS_Optical_Depth_047_Aqua_Nearest4~',CommonFormula],'Weights',1./(1-weights));
            [~,~,lme_rand] = randomEffects(lme);
            lme_coef = lme.Coefficients;
            save(TempFileName,'lme_rand','lme_coef','-v7.3');
            clear lme_weight lme 
        end
    end

    %MAIACUS_Optical_Depth_055_Aqua_Nearest4
    if(ismember('MAIACUS_Optical_Depth_055_Aqua_Nearest4',Data.Properties.VariableNames))
        TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'MAIACUS_Optical_Depth_055_Aqua_Nearest4_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
        if(exist(TempFileName,'file'))
            %LoadData_function(TempFileName);
        else
            fprintf('imputation fitting: %s\n',['MAIACUS_Optical_Depth_055_Aqua_Nearest4 ~',CommonFormula]);
            Data.Is_nan = isnan(Data.MAIACUS_Optical_Depth_055_Aqua_Nearest4);
            lme_weight = fitglm(Data,['Is_nan~',CommonFormula_weight],'Distribution','Binomial','Link','logit'); 
            weights = predict(lme_weight,Data);
            weights(isnan(weights)) = 0;
            lme = fitlme(Data(:,['MAIACUS_Optical_Depth_055_Aqua_Nearest4',CurrentVariableList]),['MAIACUS_Optical_Depth_055_Aqua_Nearest4~',CommonFormula],'Weights',1./(1-weights));
            [~,~,lme_rand] = randomEffects(lme);
            lme_coef = lme.Coefficients;
            save(TempFileName,'lme_rand','lme_coef','-v7.3');
            clear lme_weight lme 
        end
    end

    %MAIACUS_Optical_Depth_047_Terra_Nearest4
    if(ismember('MAIACUS_Optical_Depth_047_Terra_Nearest4',Data.Properties.VariableNames))
        TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'MAIACUS_Optical_Depth_047_Terra_Nearest4_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
        if(exist(TempFileName,'file'))
            %LoadData_function(TempFileName);
        else
            fprintf('imputation fitting: %s\n',['MAIACUS_Optical_Depth_047_Terra_Nearest4 ~',CommonFormula]);
            Data.Is_nan = isnan(Data.MAIACUS_Optical_Depth_047_Terra_Nearest4);
            lme_weight = fitglm(Data,['Is_nan~',CommonFormula_weight],'Distribution','Binomial','Link','logit'); 
            weights = predict(lme_weight,Data);
            weights(isnan(weights)) = 0;
            lme = fitlme(Data(:,['MAIACUS_Optical_Depth_047_Terra_Nearest4',CurrentVariableList]),['MAIACUS_Optical_Depth_047_Terra_Nearest4~',CommonFormula],'Weights',1./(1-weights));
            [~,~,lme_rand] = randomEffects(lme);
            lme_coef = lme.Coefficients;
            save(TempFileName,'lme_rand','lme_coef','-v7.3');
            clear lme_weight lme 
        end
    end

    %MAIACUS_Optical_Depth_055_Terra_Nearest4
    if(ismember('MAIACUS_Optical_Depth_055_Terra_Nearest4',Data.Properties.VariableNames))
        TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'MAIACUS_Optical_Depth_055_Terra_Nearest4_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
        if(exist(TempFileName,'file'))
            %LoadData_function(TempFileName);
        else
            fprintf('imputation fitting: %s\n',['MAIACUS_Optical_Depth_055_Terra_Nearest4 ~',CommonFormula]);
            Data.Is_nan = isnan(Data.MAIACUS_Optical_Depth_055_Terra_Nearest4);
            lme_weight = fitglm(Data,['Is_nan~REANALYSIS_air_sfc_DailyMean+REANALYSIS_shum_2m_DailyMean+REANALYSIS_snowc_DailyMean+REANALYSIS_tcdc_DailyMean+USElevation_mea100'],'Distribution','Binomial','Link','logit'); 
            weights = predict(lme_weight,Data);
            weights(isnan(weights)) = 0;
            lme = fitlme(Data(:,['MAIACUS_Optical_Depth_055_Terra_Nearest4',CurrentVariableList]),['MAIACUS_Optical_Depth_055_Terra_Nearest4~',CommonFormula],'Weights',1./(1-weights));
            [~,~,lme_rand] = randomEffects(lme);
            lme_coef = lme.Coefficients;
            save(TempFileName,'lme_rand','lme_coef','-v7.3');
            clear lme_weight lme 
        end
    end
    
    %MOD09A1
    if(ismember('MOD09A1',Data.Properties.VariableNames))
        TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'MOD09A1_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
        if(exist(TempFileName,'file'))
            %LoadData_function(TempFileName);
        else
            fprintf('imputation fitting: %s\n',['MOD09A1 ~',CommonFormula]);
            Data.Is_nan = isnan(Data.MOD09A1);
            lme_weight = fitglm(Data,['Is_nan~',CommonFormula_weight],'Distribution','Binomial','Link','logit'); 
            weights = predict(lme_weight,Data);
            weights(isnan(weights)) = 0;
            lme = fitlme(Data(:,['MOD09A1',CurrentVariableList]),['MOD09A1~',CommonFormula],'Weights',1./(1-weights));
            [~,~,lme_rand] = randomEffects(lme);
            lme_coef = lme.Coefficients;
            save(TempFileName,'lme_rand','lme_coef','-v7.3');
            clear lme_weight lme 
        end
    end

    %MOD13A2_Nearest4
    if(ismember('MOD13A2_Nearest4',Data.Properties.VariableNames))
        TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'MOD13A2_Nearest4_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
        if(exist(TempFileName,'file'))
            %LoadData_function(TempFileName);
        else
            fprintf('imputation fitting: %s\n',['MOD13A2_Nearest4 ~',CommonFormula]);
            Data.Is_nan = isnan(Data.MOD13A2_Nearest4);
            lme_weight = fitglm(Data,['Is_nan~',CommonFormula_weight],'Distribution','Binomial','Link','logit'); 
            weights = predict(lme_weight,Data);
            weights(isnan(weights)) = 0;
            lme = fitlme(Data(:,['MOD13A2_Nearest4',CurrentVariableList]),['MOD13A2_Nearest4~',CommonFormula],'Weights',1./(1-weights));
            [~,~,lme_rand] = randomEffects(lme);
            lme_coef = lme.Coefficients;
            save(TempFileName,'lme_rand','lme_coef','-v7.3');
            clear lme_weight lme 
        end
    end

    %MOD11A1_LST_Day_1km_Nearest4
    if(ismember('MOD11A1_LST_Day_1km_Nearest4',Data.Properties.VariableNames))
        TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'MOD11A1_LST_Day_1km_Nearest4_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
        if(exist(TempFileName,'file'))
            %LoadData_function(TempFileName);
        else
            fprintf('imputation fitting: %s\n',['MOD11A1_LST_Day_1km_Nearest4 ~',CommonFormula]);
            Data.Is_nan = isnan(Data.MOD11A1_LST_Day_1km_Nearest4);
            lme_weight = fitglm(Data,['Is_nan~',CommonFormula_weight],'Distribution','Binomial','Link','logit'); 
            weights = predict(lme_weight,Data);
            weights(isnan(weights)) = 0;
            lme = fitlme(Data(:,['MOD11A1_LST_Day_1km_Nearest4',CurrentVariableList]),['MOD11A1_LST_Day_1km_Nearest4~',CommonFormula],'Weights',1./(1-weights));
            [~,~,lme_rand] = randomEffects(lme);
            lme_coef = lme.Coefficients;
            save(TempFileName,'lme_rand','lme_coef','-v7.3');
            clear lme_weight lme 
        end
    end

    %MOD11A1_LST_Night_1km_Nearest4
    if(ismember('MOD11A1_LST_Night_1km_Nearest4',Data.Properties.VariableNames))
        TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'MOD11A1_LST_Night_1km_Nearest4_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
        if(exist(TempFileName,'file'))
            %LoadData_function(TempFileName);
        else
            fprintf('imputation fitting: %s\n',['MOD11A1_LST_Night_1km_Nearest4 ~',CommonFormula]);
            Data.Is_nan = isnan(Data.MOD11A1_LST_Night_1km_Nearest4);
            lme_weight = fitglm(Data,['Is_nan~',CommonFormula_weight],'Distribution','Binomial','Link','logit'); 
            weights = predict(lme_weight,Data);
            weights(isnan(weights)) = 0;
            lme = fitlme(Data(:,['MOD11A1_LST_Night_1km_Nearest4',CurrentVariableList]),['MOD11A1_LST_Night_1km_Nearest4~',CommonFormula],'Weights',1./(1-weights));
            [~,~,lme_rand] = randomEffects(lme);
            lme_coef = lme.Coefficients;
            save(TempFileName,'lme_rand','lme_coef','-v7.3');
            clear lme_weight lme 
        end
    end

    %MOD11A1_Clear_day_cov_Nearest4
    if(ismember('MOD11A1_Clear_day_cov_Nearest4',Data.Properties.VariableNames))
        TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'MOD11A1_Clear_day_cov_Nearest4_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
        if(exist(TempFileName,'file'))
            %LoadData_function(TempFileName);
        else
            fprintf('imputation fitting: %s\n',['MOD11A1_Clear_day_cov_Nearest4 ~',CommonFormula]);
            Data.Is_nan = isnan(Data.MOD11A1_Clear_day_cov_Nearest4);
            lme_weight = fitglm(Data,['Is_nan~',CommonFormula_weight],'Distribution','Binomial','Link','logit'); 
            weights = predict(lme_weight,Data);
            weights(isnan(weights)) = 0;
            lme = fitlme(Data(:,['MOD11A1_Clear_day_cov_Nearest4',CurrentVariableList]),['MOD11A1_Clear_day_cov_Nearest4~',CommonFormula],'Weights',1./(1-weights));
            [~,~,lme_rand] = randomEffects(lme);
            lme_coef = lme.Coefficients;
            save(TempFileName,'lme_rand','lme_coef','-v7.3');
            clear lme_weight lme 
        end
    end

    %MOD11A1_Clear_night_cov_Nearest4
    if(ismember('MOD11A1_Clear_night_cov_Nearest4',Data.Properties.VariableNames))
        TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'MOD11A1_Clear_night_cov_Nearest4_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
        if(exist(TempFileName,'file'))
            %LoadData_function(TempFileName);
        else
            fprintf('imputation fitting: %s\n',['MOD11A1_Clear_night_cov_Nearest4 ~',CommonFormula]);
            Data.Is_nan = isnan(Data.MOD11A1_Clear_night_cov_Nearest4);
            lme_weight = fitglm(Data,['Is_nan~',CommonFormula_weight],'Distribution','Binomial','Link','logit'); 
            weights = predict(lme_weight,Data);
            weights(isnan(weights)) = 0;
            lme = fitlme(Data(:,['MOD11A1_Clear_night_cov_Nearest4',CurrentVariableList]),['MOD11A1_Clear_night_cov_Nearest4~',CommonFormula],'Weights',1./(1-weights));
            [~,~,lme_rand] = randomEffects(lme);
            lme_coef = lme.Coefficients;
            save(TempFileName,'lme_rand','lme_coef','-v7.3');
            clear lme_weight lme 
        end
    end

    %OMTO3e_ColumnAmountO3
    if(ismember('OMTO3e_ColumnAmountO3',Data.Properties.VariableNames))
        TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'OMTO3e_ColumnAmountO3_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
        if(exist(TempFileName,'file'))
            %LoadData_function(TempFileName);
        else
            fprintf('imputation fitting: %s\n',['OMTO3e_ColumnAmountO3 ~',CommonFormula]);
            Data.Is_nan = isnan(Data.OMTO3e_ColumnAmountO3);
            lme_weight = fitglm(Data,['Is_nan~',CommonFormula_weight],'Distribution','Binomial','Link','logit'); 
            weights = predict(lme_weight,Data);
            weights(isnan(weights)) = 0;
            lme = fitlme(Data(:,['OMTO3e_ColumnAmountO3',CurrentVariableList]),['OMTO3e_ColumnAmountO3~',CommonFormula],'Weights',1./(1-weights));
            [~,~,lme_rand] = randomEffects(lme);
            lme_coef = lme.Coefficients;
            save(TempFileName,'lme_rand','lme_coef','-v7.3');
            clear lme_weight lme 
        end
    end

    %OMUVBd_UVindex_Mean
    if(ismember('OMUVBd_UVindex_Mean',Data.Properties.VariableNames))
        TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'OMUVBd_UVindex_Mean_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
        if(exist(TempFileName,'file'))
            %LoadData_function(TempFileName);
        else
            fprintf('imputation fitting: %s\n',['OMUVBd_UVindex_Mean ~',CommonFormula]);
            Data.Is_nan = isnan(Data.OMUVBd_UVindex_Mean);
            lme_weight = fitglm(Data,['Is_nan~',CommonFormula_weight],'Distribution','Binomial','Link','logit'); 
            weights = predict(lme_weight,Data);
            weights(isnan(weights)) = 0;
            lme = fitlme(Data(:,['OMUVBd_UVindex_Mean',CurrentVariableList]),['OMUVBd_UVindex_Mean~',CommonFormula],'Weights',1./(1-weights));
            [~,~,lme_rand] = randomEffects(lme);
            lme_coef = lme.Coefficients;
            save(TempFileName,'lme_rand','lme_coef','-v7.3');
            clear lme_weight lme 
        end
    end

    %OMSO2e_ColumnAmountSO2_PBL_Mean
    if(ismember('OMSO2e_ColumnAmountSO2_PBL_Mean',Data.Properties.VariableNames))
        TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'OMSO2e_ColumnAmountSO2_PBL_Mean_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
        if(exist(TempFileName,'file'))
            %LoadData_function(TempFileName);
        else
            fprintf('imputation fitting: %s\n',['OMSO2e_ColumnAmountSO2_PBL_Mean ~',CommonFormula]);
            Data.Is_nan = isnan(Data.OMSO2e_ColumnAmountSO2_PBL_Mean);
            lme_weight = fitglm(Data,['Is_nan~',CommonFormula_weight],'Distribution','Binomial','Link','logit'); 
            weights = predict(lme_weight,Data);
            weights(isnan(weights)) = 0;
            lme = fitlme(Data(:,['OMSO2e_ColumnAmountSO2_PBL_Mean',CurrentVariableList]),['OMSO2e_ColumnAmountSO2_PBL_Mean~',CommonFormula],'Weights',1./(1-weights));
            [~,~,lme_rand] = randomEffects(lme);
            lme_coef = lme.Coefficients;
            save(TempFileName,'lme_rand','lme_coef','-v7.3');
            clear lme_weight lme 
        end
    end

    %OMNO2d_ColumnAmountNO2StratoCloudScreened_Mean
    if(ismember('OMNO2d_ColumnAmountNO2StratoCloudScreened_Mean',Data.Properties.VariableNames))
        TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'OMNO2d_ColumnAmountNO2StratoCloudScreened_Mean_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
        if(exist(TempFileName,'file'))
            %LoadData_function(TempFileName);
        else
            fprintf('imputation fitting: %s\n',['OMNO2d_ColumnAmountNO2StratoCloudScreened_Mean ~',CommonFormula]);
            Data.Is_nan = isnan(Data.OMNO2d_ColumnAmountNO2StratoCloudScreened_Mean);
            lme_weight = fitglm(Data,['Is_nan~',CommonFormula_weight],'Distribution','Binomial','Link','logit'); 
            weights = predict(lme_weight,Data);
            weights(isnan(weights)) = 0;
            lme = fitlme(Data(:,['OMNO2d_ColumnAmountNO2StratoCloudScreened_Mean',CurrentVariableList]),['OMNO2d_ColumnAmountNO2StratoCloudScreened_Mean~',CommonFormula],'Weights',1./(1-weights));
            [~,~,lme_rand] = randomEffects(lme);
            lme_coef = lme.Coefficients;
            save(TempFileName,'lme_rand','lme_coef','-v7.3');
            clear lme_weight lme 
        end
    end

    %OMAERUVd_UVAerosolIndex_Mean
    if(ismember('OMAERUVd_UVAerosolIndex_Mean',Data.Properties.VariableNames))
        TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'OMAERUVd_UVAerosolIndex_Mean_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
        if(exist(TempFileName,'file'))
            %LoadData_function(TempFileName);
        else
            fprintf('imputation fitting: %s\n',['OMAERUVd_UVAerosolIndex_Mean ~',CommonFormula]);
            Data.Is_nan = isnan(Data.OMAERUVd_UVAerosolIndex_Mean);
            lme_weight = fitglm(Data,['Is_nan~',CommonFormula_weight],'Distribution','Binomial','Link','logit'); 
            weights = predict(lme_weight,Data);
            weights(isnan(weights)) = 0;
            lme = fitlme(Data(:,['OMAERUVd_UVAerosolIndex_Mean',CurrentVariableList]),['OMAERUVd_UVAerosolIndex_Mean~',CommonFormula],'Weights',1./(1-weights));
            [~,~,lme_rand] = randomEffects(lme);
            lme_coef = lme.Coefficients;
            save(TempFileName,'lme_rand','lme_coef','-v7.3');
            clear lme_weight lme 
        end
    end
    
    %OMAEROe_VISAerosolIndex_Mean
    if(ismember('OMAEROe_VISAerosolIndex_Mean',Data.Properties.VariableNames))
        TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'OMAEROe_VISAerosolIndex_Mean_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
        if(exist(TempFileName,'file'))
            %LoadData_function(TempFileName);
        else
            fprintf('imputation fitting: %s\n',['OMAEROe_VISAerosolIndex_Mean ~',CommonFormula]);
            Data.Is_nan = isnan(Data.OMAEROe_VISAerosolIndex_Mean);
            lme_weight = fitglm(Data,['Is_nan~',CommonFormula_weight],'Distribution','Binomial','Link','logit'); 
            weights = predict(lme_weight,Data);
            weights(isnan(weights)) = 0;
            lme = fitlme(Data(:,['OMAEROe_VISAerosolIndex_Mean',CurrentVariableList]),['OMAEROe_VISAerosolIndex_Mean~',CommonFormula],'Weights',1./(1-weights));
            [~,~,lme_rand] = randomEffects(lme);
            lme_coef = lme.Coefficients;
            save(TempFileName,'lme_rand','lme_coef','-v7.3');
            clear lme_weight lme 
        end
    end

    %OMAEROe_UVAerosolIndex_Mean
    if(ismember('OMAEROe_UVAerosolIndex_Mean',Data.Properties.VariableNames))
        TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'OMAEROe_UVAerosolIndex_Mean_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
        if(exist(TempFileName,'file'))
            %LoadData_function(TempFileName);
        else
            fprintf('imputation fitting: %s\n',['OMAEROe_UVAerosolIndex_Mean ~',CommonFormula]);
            Data.Is_nan = isnan(Data.OMAEROe_UVAerosolIndex_Mean);
            lme_weight = fitglm(Data,['Is_nan~',CommonFormula_weight],'Distribution','Binomial','Link','logit'); 
            weights = predict(lme_weight,Data);
            weights(isnan(weights)) = 0;
            lme = fitlme(Data(:,['OMAEROe_UVAerosolIndex_Mean',CurrentVariableList]),['OMAEROe_UVAerosolIndex_Mean~',CommonFormula],'Weights',1./(1-weights));
            [~,~,lme_rand] = randomEffects(lme);
            lme_coef = lme.Coefficients;
            save(TempFileName,'lme_rand','lme_coef','-v7.3');
            clear lme_weight lme 
        end
    end
    
    %OMO3PR
    if(ismember('OMO3PR',Data.Properties.VariableNames))
        TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'OMO3PR_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
        if(exist(TempFileName,'file'))
            %LoadData_function(TempFileName);
        else
            fprintf('imputation fitting: %s\n',['OMO3PR ~',CommonFormula]);
            Data.Is_nan = isnan(Data.OMO3PR);
            lme_weight = fitglm(Data,['Is_nan~',CommonFormula_weight],'Distribution','Binomial','Link','logit'); 
            weights = predict(lme_weight,Data);
            weights(isnan(weights)) = 0;
            lme = fitlme(Data(:,['OMO3PR',CurrentVariableList]),['OMO3PR~',CommonFormula],'Weights',1./(1-weights));
            [~,~,lme_rand] = randomEffects(lme);
            lme_coef = lme.Coefficients;
            save(TempFileName,'lme_rand','lme_coef','-v7.3');
            clear lme_weight lme 
        end
    end
    
    %MOD04L2_550
    if(ismember('MOD04L2_550',Data.Properties.VariableNames))
        TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'MOD04L2_550_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
        if(exist(TempFileName,'file'))
            %LoadData_function(TempFileName);
        else
            fprintf('imputation fitting: %s\n',['MOD04L2_550 ~',CommonFormula]);
            Data.Is_nan = isnan(Data.MOD04L2_550);
            lme_weight = fitglm(Data,['Is_nan~',CommonFormula_weight],'Distribution','Binomial','Link','logit'); 
            weights = predict(lme_weight,Data);
            weights(isnan(weights)) = 0;
            lme = fitlme(Data(:,['MOD04L2_550',CurrentVariableList]),['MOD04L2_550~',CommonFormula],'Weights',1./(1-weights));
            [~,~,lme_rand] = randomEffects(lme);
            lme_coef = lme.Coefficients;
            save(TempFileName,'lme_rand','lme_coef','-v7.3');
            clear lme_weight lme 
        end
    end
    
    %MAIACUS_cosVZA_Terra_Nearest
    if(ismember('MAIACUS_cosVZA_Terra_Nearest',Data.Properties.VariableNames))
        TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'MAIACUS_cosVZA_Terra_Nearest_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
        if(exist(TempFileName,'file'))
            %LoadData_function(TempFileName);
        else
            fprintf('imputation fitting: %s\n',['MAIACUS_cosVZA_Terra_Nearest ~',CommonFormula]);
            Data.Is_nan = isnan(Data.MAIACUS_cosVZA_Terra_Nearest);
            lme_weight = fitglm(Data,['Is_nan~',CommonFormula_weight],'Distribution','Binomial','Link','logit'); 
            weights = predict(lme_weight,Data);
            weights(isnan(weights)) = 0;
            lme = fitlme(Data(:,['MAIACUS_cosVZA_Terra_Nearest',CurrentVariableList]),['MAIACUS_cosVZA_Terra_Nearest~',CommonFormula],'Weights',1./(1-weights));
            [~,~,lme_rand] = randomEffects(lme);
            lme_coef = lme.Coefficients;
            save(TempFileName,'lme_rand','lme_coef','-v7.3');
            clear lme_weight lme 
        end
    end
    
    %MAIACUS_cosVZA_Aqua_Nearest
    if(ismember('MAIACUS_cosVZA_Aqua_Nearest',Data.Properties.VariableNames))
        TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'MAIACUS_cosVZA_Aqua_Nearest_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
        if(exist(TempFileName,'file'))
            %LoadData_function(TempFileName);
        else
            fprintf('imputation fitting: %s\n',['MAIACUS_cosVZA_Aqua_Nearest ~',CommonFormula]);
            Data.Is_nan = isnan(Data.MAIACUS_cosVZA_Aqua_Nearest);
            lme_weight = fitglm(Data,['Is_nan~',CommonFormula_weight],'Distribution','Binomial','Link','logit'); 
            weights = predict(lme_weight,Data);
            weights(isnan(weights)) = 0;
            lme = fitlme(Data(:,['MAIACUS_cosVZA_Aqua_Nearest',CurrentVariableList]),['MAIACUS_cosVZA_Aqua_Nearest~',CommonFormula],'Weights',1./(1-weights));
            [~,~,lme_rand] = randomEffects(lme);
            lme_coef = lme.Coefficients;
            save(TempFileName,'lme_rand','lme_coef','-v7.3');
            clear lme_weight lme 
        end
    end
    
end
	

%% predict!!!  
%REANALYSIS_soilm_DailyMean
if(ismember('REANALYSIS_soilm_DailyMean',Data.Properties.VariableNames))
    TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'REANALYSIS_soilm_DailyMean_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
    Temp = LoadData_function(TempFileName);
    lme_rand = Temp.lme_rand;
    lme_coef = Temp.lme_coef;
    Index = isnan(Data.REANALYSIS_soilm_DailyMean);
    Pred = LocalPredict_function(lme_coef,lme_rand,Data(:,CurrentVariableList));
    Data.REANALYSIS_soilm_DailyMean(Index) = Pred(Index);
    fprintf('imputation prediction!: %s\n','REANALYSIS_soilm_DailyMean',corr(Pred,Data.REANALYSIS_soilm_DailyMean,'rows','complete')^2);
end

%REANALYSIS_gflux_DailyMean
if(ismember('REANALYSIS_gflux_DailyMean',Data.Properties.VariableNames))
    TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'REANALYSIS_gflux_DailyMean_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
    Temp = LoadData_function(TempFileName);
    lme_rand = Temp.lme_rand;
    lme_coef = Temp.lme_coef;
    Index = isnan(Data.REANALYSIS_gflux_DailyMean);
    Pred = LocalPredict_function(lme_coef,lme_rand,Data(:,CurrentVariableList));
    Data.REANALYSIS_gflux_DailyMean(Index) = Pred(Index);
    fprintf('imputation prediction!: %s\n','REANALYSIS_gflux_DailyMean',corr(Pred,Data.REANALYSIS_gflux_DailyMean,'rows','complete')^2);
end

%MAIACUS_Optical_Depth_047_Aqua_Nearest4
if(ismember('MAIACUS_Optical_Depth_047_Aqua_Nearest4',Data.Properties.VariableNames))
    TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'MAIACUS_Optical_Depth_047_Aqua_Nearest4_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
    Temp = LoadData_function(TempFileName);
    lme_rand = Temp.lme_rand;
    lme_coef = Temp.lme_coef;
    Index = isnan(Data.MAIACUS_Optical_Depth_047_Aqua_Nearest4);
    Pred = LocalPredict_function(lme_coef,lme_rand,Data(:,CurrentVariableList));
    Data.MAIACUS_Optical_Depth_047_Aqua_Nearest4(Index) = Pred(Index);
    fprintf('imputation prediction!: %s\n','MAIACUS_Optical_Depth_047_Aqua_Nearest4',corr(Pred,Data.MAIACUS_Optical_Depth_047_Aqua_Nearest4,'rows','complete')^2);
end

%MAIACUS_Optical_Depth_055_Aqua_Nearest4
if(ismember('MAIACUS_Optical_Depth_055_Aqua_Nearest4',Data.Properties.VariableNames))
    TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'MAIACUS_Optical_Depth_055_Aqua_Nearest4_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
    Temp = LoadData_function(TempFileName);
    lme_rand = Temp.lme_rand;
    lme_coef = Temp.lme_coef;
    Index = isnan(Data.MAIACUS_Optical_Depth_055_Aqua_Nearest4);
    Pred = LocalPredict_function(lme_coef,lme_rand,Data(:,CurrentVariableList));
    Data.MAIACUS_Optical_Depth_055_Aqua_Nearest4(Index) = Pred(Index);
    fprintf('imputation prediction!: %s\n','MAIACUS_Optical_Depth_055_Aqua_Nearest4',corr(Pred,Data.MAIACUS_Optical_Depth_055_Aqua_Nearest4,'rows','complete')^2);
end

%MAIACUS_Optical_Depth_047_Terra_Nearest4
if(ismember('MAIACUS_Optical_Depth_047_Terra_Nearest4',Data.Properties.VariableNames))
    TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'MAIACUS_Optical_Depth_047_Terra_Nearest4_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
    Temp = LoadData_function(TempFileName);
    lme_rand = Temp.lme_rand;
    lme_coef = Temp.lme_coef;
    Index = isnan(Data.MAIACUS_Optical_Depth_047_Terra_Nearest4);
    Pred = LocalPredict_function(lme_coef,lme_rand,Data(:,CurrentVariableList));
    Data.MAIACUS_Optical_Depth_047_Terra_Nearest4(Index) = Pred(Index);
    fprintf('imputation prediction!: %s\n','MAIACUS_Optical_Depth_047_Terra_Nearest4',corr(Pred,Data.MAIACUS_Optical_Depth_047_Terra_Nearest4,'rows','complete')^2);
end

%MAIACUS_Optical_Depth_055_Terra_Nearest4
if(ismember('MAIACUS_Optical_Depth_055_Terra_Nearest4',Data.Properties.VariableNames))
    TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'MAIACUS_Optical_Depth_055_Terra_Nearest4_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
    Temp = LoadData_function(TempFileName);
    lme_rand = Temp.lme_rand;
    lme_coef = Temp.lme_coef;
    Index = isnan(Data.MAIACUS_Optical_Depth_055_Terra_Nearest4);
    Pred = LocalPredict_function(lme_coef,lme_rand,Data(:,CurrentVariableList));
    Data.MAIACUS_Optical_Depth_055_Terra_Nearest4(Index) = Pred(Index);
    fprintf('imputation prediction!: %s\n','MAIACUS_Optical_Depth_055_Terra_Nearest4',corr(Pred,Data.MAIACUS_Optical_Depth_055_Terra_Nearest4,'rows','complete')^2);
end

%MOD09A1
if(ismember('MOD09A1',Data.Properties.VariableNames))
    TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'MOD09A1_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
    Temp = LoadData_function(TempFileName);
    lme_rand = Temp.lme_rand;
    lme_coef = Temp.lme_coef;
    Index = isnan(Data.MOD09A1);
    Pred = LocalPredict_function(lme_coef,lme_rand,Data(:,CurrentVariableList));
    Data.MOD09A1(Index) = Pred(Index);
    fprintf('imputation prediction!: %s\n','MOD09A1',corr(Pred,Data.MOD09A1,'rows','complete')^2);
end

%MOD13A2_Nearest4
if(ismember('MOD13A2_Nearest4',Data.Properties.VariableNames))
    TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'MOD13A2_Nearest4_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
    Temp = LoadData_function(TempFileName);
    lme_rand = Temp.lme_rand;
    lme_coef = Temp.lme_coef;
    Index = isnan(Data.MOD13A2_Nearest4);
    Pred = LocalPredict_function(lme_coef,lme_rand,Data(:,CurrentVariableList));
    Data.MOD13A2_Nearest4(Index) = Pred(Index);
    fprintf('imputation prediction!: %s\n','MOD13A2_Nearest4',corr(Pred,Data.MOD13A2_Nearest4,'rows','complete')^2);
end

%MOD11A1_LST_Day_1km_Nearest4
if(ismember('MOD11A1_LST_Day_1km_Nearest4',Data.Properties.VariableNames))
    TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'MOD11A1_LST_Day_1km_Nearest4_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
    Temp = LoadData_function(TempFileName);
    lme_rand = Temp.lme_rand;
    lme_coef = Temp.lme_coef;
    Index = isnan(Data.MOD11A1_LST_Day_1km_Nearest4);
    Pred = LocalPredict_function(lme_coef,lme_rand,Data(:,CurrentVariableList));
    Data.MOD11A1_LST_Day_1km_Nearest4(Index) = Pred(Index);
    fprintf('imputation prediction!: %s\n','MOD11A1_LST_Day_1km_Nearest4',corr(Pred,Data.MOD11A1_LST_Day_1km_Nearest4,'rows','complete')^2);
end

%MOD11A1_LST_Night_1km_Nearest4
if(ismember('MOD11A1_LST_Night_1km_Nearest4',Data.Properties.VariableNames))
    TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'MOD11A1_LST_Night_1km_Nearest4_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
    Temp = LoadData_function(TempFileName);
    lme_rand = Temp.lme_rand;
    lme_coef = Temp.lme_coef;
    Index = isnan(Data.MOD11A1_LST_Night_1km_Nearest4);
    Pred = LocalPredict_function(lme_coef,lme_rand,Data(:,CurrentVariableList));
    Data.MOD11A1_LST_Night_1km_Nearest4(Index) = Pred(Index);
    fprintf('imputation prediction!: %s\n','MOD11A1_LST_Night_1km_Nearest4',corr(Pred,Data.MOD11A1_LST_Night_1km_Nearest4,'rows','complete')^2);
end

%MOD11A1_Clear_day_cov_Nearest4
if(ismember('MOD11A1_Clear_day_cov_Nearest4',Data.Properties.VariableNames))
    TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'MOD11A1_Clear_day_cov_Nearest4_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
    Temp = LoadData_function(TempFileName);
    lme_rand = Temp.lme_rand;
    lme_coef = Temp.lme_coef;
    Index = isnan(Data.MOD11A1_Clear_day_cov_Nearest4);
    Pred = LocalPredict_function(lme_coef,lme_rand,Data(:,CurrentVariableList));
    Data.MOD11A1_Clear_day_cov_Nearest4(Index) = Pred(Index);
    fprintf('imputation prediction!: %s\n','MOD11A1_Clear_day_cov_Nearest4',corr(Pred,Data.MOD11A1_Clear_day_cov_Nearest4,'rows','complete')^2);
end

%MOD11A1_Clear_night_cov_Nearest4
if(ismember('MOD11A1_Clear_night_cov_Nearest4',Data.Properties.VariableNames))
    TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'MOD11A1_Clear_night_cov_Nearest4_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
    Temp = LoadData_function(TempFileName);
    lme_rand = Temp.lme_rand;
    lme_coef = Temp.lme_coef;
    Index = isnan(Data.MOD11A1_Clear_night_cov_Nearest4);
    Pred = LocalPredict_function(lme_coef,lme_rand,Data(:,CurrentVariableList));
    Data.MOD11A1_Clear_night_cov_Nearest4(Index) = Pred(Index);
    fprintf('imputation prediction!: %s\n','MOD11A1_Clear_night_cov_Nearest4',corr(Pred,Data.MOD11A1_Clear_night_cov_Nearest4,'rows','complete')^2);
end

%OMTO3e_ColumnAmountO3
if(ismember('OMTO3e_ColumnAmountO3',Data.Properties.VariableNames))
    TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'OMTO3e_ColumnAmountO3_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
    Temp = LoadData_function(TempFileName);
    lme_rand = Temp.lme_rand;
    lme_coef = Temp.lme_coef;
    Index = isnan(Data.OMTO3e_ColumnAmountO3);
    Pred = LocalPredict_function(lme_coef,lme_rand,Data(:,CurrentVariableList));
    Data.OMTO3e_ColumnAmountO3(Index) = Pred(Index);
    fprintf('imputation prediction!: %s\n','OMTO3e_ColumnAmountO3',corr(Pred,Data.OMTO3e_ColumnAmountO3,'rows','complete')^2);
end

%OMUVBd_UVindex_Mean
if(ismember('OMUVBd_UVindex_Mean',Data.Properties.VariableNames))
    TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'OMUVBd_UVindex_Mean_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
    Temp = LoadData_function(TempFileName);
    lme_rand = Temp.lme_rand;
    lme_coef = Temp.lme_coef;
    Index = isnan(Data.OMUVBd_UVindex_Mean);
    Pred = LocalPredict_function(lme_coef,lme_rand,Data(:,CurrentVariableList));
    Data.OMUVBd_UVindex_Mean(Index) = Pred(Index);
    fprintf('imputation prediction!: %s\n','OMUVBd_UVindex_Mean',corr(Pred,Data.OMUVBd_UVindex_Mean,'rows','complete')^2);
end

%OMSO2e_ColumnAmountSO2_PBL_Mean
if(ismember('OMSO2e_ColumnAmountSO2_PBL_Mean',Data.Properties.VariableNames))
    TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'OMSO2e_ColumnAmountSO2_PBL_Mean_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
    Temp = LoadData_function(TempFileName);
    lme_rand = Temp.lme_rand;
    lme_coef = Temp.lme_coef;
    Index = isnan(Data.OMSO2e_ColumnAmountSO2_PBL_Mean);
    Pred = LocalPredict_function(lme_coef,lme_rand,Data(:,CurrentVariableList));
    Data.OMSO2e_ColumnAmountSO2_PBL_Mean(Index) = Pred(Index);
    fprintf('imputation prediction!: %s\n','OMSO2e_ColumnAmountSO2_PBL_Mean',corr(Pred,Data.OMSO2e_ColumnAmountSO2_PBL_Mean,'rows','complete')^2);
end

%OMNO2d_ColumnAmountNO2StratoCloudScreened_Mean
if(ismember('OMNO2d_ColumnAmountNO2StratoCloudScreened_Mean',Data.Properties.VariableNames))
    TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'OMNO2d_ColumnAmountNO2StratoCloudScreened_Mean_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
    Temp = LoadData_function(TempFileName);
    lme_rand = Temp.lme_rand;
    lme_coef = Temp.lme_coef;
    Index = isnan(Data.OMNO2d_ColumnAmountNO2StratoCloudScreened_Mean);
    Pred = LocalPredict_function(lme_coef,lme_rand,Data(:,CurrentVariableList));
    Data.OMNO2d_ColumnAmountNO2StratoCloudScreened_Mean(Index) = Pred(Index);
    fprintf('imputation prediction!: %s\n','OMNO2d_ColumnAmountNO2StratoCloudScreened_Mean',corr(Pred,Data.OMNO2d_ColumnAmountNO2StratoCloudScreened_Mean,'rows','complete')^2);
end

%OMAERUVd_UVAerosolIndex_Mean
if(ismember('OMAERUVd_UVAerosolIndex_Mean',Data.Properties.VariableNames))
    TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'OMAERUVd_UVAerosolIndex_Mean_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
    Temp = LoadData_function(TempFileName);
    lme_rand = Temp.lme_rand;
    lme_coef = Temp.lme_coef;
    Index = isnan(Data.OMAERUVd_UVAerosolIndex_Mean);
    Pred = LocalPredict_function(lme_coef,lme_rand,Data(:,CurrentVariableList));
    Data.OMAERUVd_UVAerosolIndex_Mean(Index) = Pred(Index);
    fprintf('imputation prediction!: %s\n','OMAERUVd_UVAerosolIndex_Mean',corr(Pred,Data.OMAERUVd_UVAerosolIndex_Mean,'rows','complete')^2);
end

%OMAEROe_VISAerosolIndex_Mean
if(ismember('OMAEROe_VISAerosolIndex_Mean',Data.Properties.VariableNames))
    TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'OMAEROe_VISAerosolIndex_Mean_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
    Temp = LoadData_function(TempFileName);
    lme_rand = Temp.lme_rand;
    lme_coef = Temp.lme_coef;
    Index = isnan(Data.OMAEROe_VISAerosolIndex_Mean);
    Pred = LocalPredict_function(lme_coef,lme_rand,Data(:,CurrentVariableList));
    Data.OMAEROe_VISAerosolIndex_Mean(Index) = Pred(Index);
    fprintf('imputation prediction!: %s\n','OMAEROe_VISAerosolIndex_Mean',corr(Pred,Data.OMAEROe_VISAerosolIndex_Mean,'rows','complete')^2);
end

%OMAEROe_UVAerosolIndex_Mean
if(ismember('OMAEROe_UVAerosolIndex_Mean',Data.Properties.VariableNames))
    TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'OMAEROe_UVAerosolIndex_Mean_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
    Temp = LoadData_function(TempFileName);
    lme_rand = Temp.lme_rand;
    lme_coef = Temp.lme_coef;
    Index = isnan(Data.OMAEROe_UVAerosolIndex_Mean);
    Pred = LocalPredict_function(lme_coef,lme_rand,Data(:,CurrentVariableList));
    Data.OMAEROe_UVAerosolIndex_Mean(Index) = Pred(Index);
    fprintf('imputation prediction!: %s\n','OMAEROe_UVAerosolIndex_Mean',corr(Pred,Data.OMAEROe_UVAerosolIndex_Mean,'rows','complete')^2);
end

%OMO3PR
if(ismember('OMO3PR',Data.Properties.VariableNames))
    TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'OMO3PR_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
    Temp = LoadData_function(TempFileName);
    lme_rand = Temp.lme_rand;
    lme_coef = Temp.lme_coef;
    Index = isnan(Data.OMO3PR);
    Pred = LocalPredict_function(lme_coef,lme_rand,Data(:,CurrentVariableList));
    Data.OMO3PR(Index) = Pred(Index);
    fprintf('imputation prediction!: %s\n','OMO3PR',corr(Pred,Data.OMO3PR,'rows','complete')^2);
end

%MOD04L2_550
if(ismember('MOD04L2_550',Data.Properties.VariableNames))
    TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'MOD04L2_550_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
    Temp = LoadData_function(TempFileName);
    lme_rand = Temp.lme_rand;
    lme_coef = Temp.lme_coef;
    Index = isnan(Data.MOD04L2_550);
    Pred = LocalPredict_function(lme_coef,lme_rand,Data(:,CurrentVariableList));
    Data.MOD04L2_550(Index) = Pred(Index);
    fprintf('imputation prediction!: %s\n','MOD04L2_550',corr(Pred,Data.MOD04L2_550,'rows','complete')^2);
end

%MAIACUS_cosVZA_Aqua_Nearest
if(ismember('MAIACUS_cosVZA_Aqua_Nearest',Data.Properties.VariableNames))
    TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'MAIACUS_cosVZA_Aqua_Nearest_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
    Temp = LoadData_function(TempFileName);
    lme_rand = Temp.lme_rand;
    lme_coef = Temp.lme_coef;
    Index = isnan(Data.MAIACUS_cosVZA_Aqua_Nearest);
    Pred = LocalPredict_function(lme_coef,lme_rand,Data(:,CurrentVariableList));
    Data.MAIACUS_cosVZA_Aqua_Nearest(Index) = Pred(Index);
    fprintf('imputation prediction!: %s\n','MAIACUS_cosVZA_Aqua_Nearest',corr(Pred,Data.MAIACUS_cosVZA_Aqua_Nearest,'rows','complete')^2);
end

%MAIACUS_cosVZA_Terra_Nearest
if(ismember('MAIACUS_cosVZA_Terra_Nearest',Data.Properties.VariableNames))
    TempFileName = [EnvironPara.DIRPATH_ROOT,'TempData',EnvironPara.Sep,'MAIACUS_cosVZA_Terra_Nearest_',EnvironPara.ImputationOption_Model,'_',num2str(StartYear),'_',num2str(EndYear),'.mat'];
    Temp = LoadData_function(TempFileName);
    lme_rand = Temp.lme_rand;
    lme_coef = Temp.lme_coef;
    Index = isnan(Data.MAIACUS_cosVZA_Terra_Nearest);
    Pred = LocalPredict_function(lme_coef,lme_rand,Data(:,CurrentVariableList));
    Data.MAIACUS_cosVZA_Terra_Nearest(Index) = Pred(Index);
    fprintf('imputation prediction!: %s\n','MAIACUS_cosVZA_Terra_Nearest',corr(Pred,Data.MAIACUS_cosVZA_Terra_Nearest,'rows','complete')^2);
end

if(any(strcmp('Is_nan',fieldnames(Data))))
    Data.Is_nan = [];
end

if(any(strcmp('REANALYSIS_windspeed_10m_DailyMean_1',fieldnames(Data))))
    Data.REANALYSIS_windspeed_10m_DailyMean_1 = [];
end

Input_var_AOD = table2array(Data);

end

% a simplified version of prediction in a linear mixed-effect model --- no
% longer in use, too slow
% function Pred = LocalPredict_function(Cofficient, RandomEffect, Data)
%     B = Cofficient.Estimate(2:end);
%     Pred = Cofficient.Estimate(1) + table2array(Data(:,Cofficient.Name(2:end)))*B;
%     RandomLevel = nominal(RandomEffect.Level);
%     RandomeEffectGroups = unique(RandomEffect.Group);
%     Pred_random = zeros(size(Data,1),1);
%     for i = 1:length(RandomeEffectGroups)
%         [Lia,Locb] = ismember(num2str(Data{:,RandomeEffectGroups{i}}),RandomLevel);
%         Temp = zeros(size(Data,1),1);
%         Temp(Lia) = RandomEffect.Estimate(Locb(Lia));
%         Pred_random = Pred_random + Temp;
%     end
%     Pred = Pred + Pred_random;
% end

function Pred = LocalPredict_function(Cofficient, RandomEffect, Data)
    B = Cofficient.Estimate(2:end);
    Pred = Cofficient.Estimate(1) + table2array(Data(:,Cofficient.Name(2:end)))*B;
    RandomLevel = str2double(RandomEffect.Level);
    RandomeEffectGroups = unique(RandomEffect.Group);
    Pred_random = zeros(size(Data,1),1);
    for i = 1:length(RandomeEffectGroups)
        [Lia,Locb] = ismember(Data{:,RandomeEffectGroups{i}},RandomLevel);
        Temp = zeros(size(Data,1),1);
        Temp(Lia) = RandomEffect.Estimate(Locb(Lia));
        Pred_random = Pred_random + Temp;
    end
    Pred = Pred + Pred_random;
end


% imputation prediction!: REANALYSIS_soilm_DailyMean
% imputation prediction!: 7.101163e-01
% imputation prediction!: MOD09A1
% imputation prediction!: 4.939229e-01
% imputation prediction!: MOD11A1_LST_Day_1km_Nearest4
% imputation prediction!: 9.561398e-01
% imputation prediction!: MOD11A1_LST_Night_1km_Nearest4
% imputation prediction!: 9.684766e-01
% imputation prediction!: MOD11A1_Clear_day_cov_Nearest4
% imputation prediction!: 1.793968e-01
% imputation prediction!: MOD11A1_Clear_night_cov_Nearest4
% imputation prediction!: 1.618356e-01
% imputation prediction!: OMTO3e_ColumnAmountO3
% imputation prediction!: 7.111570e-01
% imputation prediction!: OMUVBd_UVindex_Mean
% imputation prediction!: 9.248228e-01
% imputation prediction!: OMSO2e_ColumnAmountSO2_PBL_Mean
% imputation prediction!: 5.659497e-02
% imputation prediction!: OMNO2d_ColumnAmountNO2StratoCloudScreened_Mean
% imputation prediction!: 9.557364e-01
% imputation prediction!: OMAERUVd_UVAerosolIndex_Mean
% imputation prediction!: 2.147329e-01
% imputation prediction!: OMAEROe_VISAerosolIndex_Mean
% imputation prediction!: 5.226963e-01
% imputation prediction!: OMAEROe_UVAerosolIndex_Mean
% imputation prediction!: 1.914576e-01









% imputation prediction!: REANALYSIS_soilm_DailyMean
% imputation prediction!: 7.671059e-01
% imputation prediction!: MAIACUS_Optical_Depth_047_Aqua_Nearest4
% imputation prediction!: 8.590144e-01
% imputation prediction!: MAIACUS_Optical_Depth_055_Aqua_Nearest4
% imputation prediction!: 8.549505e-01
% imputation prediction!: MAIACUS_Optical_Depth_047_Terra_Nearest4
% imputation prediction!: 8.436790e-01
% imputation prediction!: MAIACUS_Optical_Depth_055_Terra_Nearest4
% imputation prediction!: 8.388082e-01
% imputation prediction!: MAIACUS_cosVZA_Aqua_Nearest
% imputation prediction!: 1.879419e-02
% imputation prediction!: MAIACUS_cosVZA_Terra_Nearest
% imputation prediction!: 1.443366e-02
% imputation prediction!: MOD09A1
% imputation prediction!: 4.563496e-01
% imputation prediction!: MOD13A2_Nearest4
% imputation prediction!: 6.732194e-01
% imputation prediction!: MOD11A1_LST_Day_1km_Nearest4
% imputation prediction!: 9.624917e-01
% imputation prediction!: MOD11A1_LST_Night_1km_Nearest4
% imputation prediction!: 9.739550e-01
% imputation prediction!: MOD11A1_Clear_day_cov_Nearest4
% imputation prediction!: 2.291519e-01
% imputation prediction!: MOD11A1_Clear_night_cov_Nearest4
% imputation prediction!: 1.856952e-01
% imputation prediction!: OMTO3e_ColumnAmountO3
% imputation prediction!: 7.933756e-01
% imputation prediction!: OMUVBd_UVindex_Mean
% imputation prediction!: 9.293583e-01
% imputation prediction!: OMSO2e_ColumnAmountSO2_PBL_Mean
% imputation prediction!: 2.379153e-01
% imputation prediction!: OMNO2d_ColumnAmountNO2StratoCloudScreened_Mean
% imputation prediction!: 9.463820e-01
% imputation prediction!: OMAERUVd_UVAerosolIndex_Mean
% imputation prediction!: 3.548365e-01
% imputation prediction!: OMAEROe_VISAerosolIndex_Mean
% imputation prediction!: 6.681533e-01
% imputation prediction!: OMAEROe_UVAerosolIndex_Mean
% imputation prediction!: 2.851441e-01
