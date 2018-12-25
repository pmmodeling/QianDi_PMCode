%% this is the controller file doing interpolation from aggregated data 
% to input data automatically, for all variables

%% version history
%% 2016-11- 12: change INPUTYEAR to INPUTDATE

%% parameters
% SITEDATA_NAME: name grid cells to be interpolated at; for example,'AQRVPM25', PM2.5
% monitors; 'EPACastNetOzone', ozone monitors; 'EPANO2', NO2 monitors
% INPUTDATE: time period of interest;
% EXTRAOPTION: 1 --> only create weight files; 0--> interpolate the data; 2 --> delete all files; 3: create maps while doing interpolation
% SystemPara: this option is only used if this file is called as a function
% by Run_ProcessDataANDRunModel; this option is set to be nan if used as a
% stand alone function

%% return values:
% ReturnStatus==0: still in processing --- some other matlab session is
% still processing data; wait for some time;
% ReturnStatus==1: processing done before --- this is the normal outputs;
% we should check aggregate data
% ReturnStatus==2, some raw data not available and needs to do it again;
% ReturnStatus==3, some input data not available has to quit;

%% example:
% Run_InterpolateDataSummary('AQRVPM25_test1',[datenum(2000,1,1),datenum(2000,12,31)],0,nan)

%% code
function ReturnStatus = Run_InterpolateDataSummary(ID,INPUTDATE,EXTRAOPTION,SystemPara)


SiteData_Name_list = {'USGrid1','USGrid2','USGrid3','USGrid4','USGrid5','USGrid6','USGrid7','USGrid8','USGrid9','USGrid10','USGrid11','USGrid12','USGrid13','USGrid14','USGrid15','USGrid16','USGrid17','USGrid18','USGrid19','USGrid20','USGrid21','USGrid22','USGrid23','USGrid24','USGrid25','USGrid26','USGrid27','USGrid28','USGrid29','USGrid30','USGrid31','USGrid32','USGrid33','USGrid34','USGrid35','USGrid36','USGrid37','USGrid38','USGrid39','USGrid40','USGrid41','USGrid42','USGrid43','USGrid44','USGrid45','USGrid46','USGrid47','USGrid48','USGrid49','USGrid50'};
if(ID>length(SiteData_Name_list))
   ID = mod(ID,length(SiteData_Name_list)); 
end
SITEDATA_NAME = SiteData_Name_list{ID};
%SITEDATA_NAME = 'EPANO2';

FORMAT = 'h5';
GCS = 'North_America_Equidistant_Conic';

%% which data set to process?
% MCD12Q1 satellite-based landuse variables -- not used any more
IsMCD12Q1 = false;% land-cover
% DATA_OPTION_LIST_MCD12Q1 = {'Mean','Var','Peak2','Nearest'};
DATA_OPTION_LIST_MCD12Q1 = {'Nearest4'};

% all kinds of meteorolgoical variables
IsMete = false;% reanalysis meteorological variables
VARNAME_LIST_Mete =  {'air.sfc','apcp','dlwrf','dswrf','evap','hpbl','gflux','lhtfl','shtfl','shum.2m','snowc','soilm','tcdc','ulwrf.sfc','uwnd.10m','vwnd.10m','weasd','prate','vis',...
    'air.sfc','hpbl','shum.2m','uwnd.10m','vwnd.10m','prate','vis',...
    'air.sfc','hpbl','shum.2m','uwnd.10m','vwnd.10m','prate','vis',...
    'air.sfc','hpbl','shum.2m','uwnd.10m','vwnd.10m','prate','vis'};
DATA_OPTION_LIST_Mete = {'DailyMean','DailyMean','DailyMean','DailyMean','DailyMean','DailyMean','DailyMean','DailyMean','DailyMean','DailyMean','DailyMean','DailyMean','DailyMean','DailyMean','DailyMean','DailyMean','DailyMean','DailyMean','DailyMean',...
    'DailyMax','DailyMax','DailyMax','DailyMax','DailyMax','DailyMax','DailyMax',...
    'DailyMin','DailyMin','DailyMin','DailyMin','DailyMin','DailyMin','DailyMin',...
    '1Day','1Day','1Day','1Day','1Day','1Day','1Day'};

% omega vertical wind speed
IsMeteOmega = false;
VARNAME_LIST_Omega =  {'omega'};
DATA_OPTION_LIST_Omega = {'DailyMean'};
   
% MOD13A1 NDVI
IsMOD13A2_MOD09A1 = false;%NDVI
VARNAME_LIST_MOD13A2 = {'MOD13A2','MOD09A1'};
DATA_OPTION_LIST_MOD13A2 = {'Nearest4'};

% MOD11A1 surface temperature
IsMOD11A1 = false;% surface temperature and others;
DATA_OPTION_LIST_MOD11A1 = {'Nearest4','Mean'};

% MAIAC AOD
IsMAIACUSAOD = false;% surface temperature and others;
VARNAME_LIST_MAIACUSAOD = {'MAIACUSAqua','MAIACUSTerra'};
DATA_OPTION_LIST_MAIACUSAOD = {'Nearest4'};

% MAIAC viewing angel data
IsMAIACUScosVZA = false;% surface temperature and others;
VARNAME_LIST_MAIACUScosVZA = {'MAIACUSTerra_cosVZA','MAIACUSAqua_cosVZA'};
DATA_OPTION_LIST_MAIACUScosVZA = {'Nearest'};

% OMI satellite data and other data set with similar structure --
% N_Day*N_Site, with only one matrix per file
IsOMI = true;
%VARNAME_LIST_OMI = {'OMAERUVd_UVAerosolIndex','OMAEROe_UVAerosolIndex','OMAEROe_VISAerosolIndex',...
%    'MOD04L2_Deep_Blue_Aerosol_Optical_Depth_550_Land','OMO3PR','OMTO3e_ColumnAmountO3',...
%    'OMNO2d_ColumnAmountNO2StratoCloudScreened','OMSO2e_ColumnAmountSO2_PBL','OMUVBd_UVindex',...
%    'GFEDFireCarbon',...
%    'CMAQ_PM25_Vertical','CMAQ_PM25_TOT','CMAQ_PM25_SO4','CMAQ_Ozone_Vertical','CMAQ_Ozone','CMAQ_NO2_Vertical','CMAQ_NO2','CMAQ_PM25_NO3','CMAQ_PM25_OC','CMAQ_PM25_EC',...
%    'MERRA2aer_SO4','MERRA2aer_OCPHOBIC','MERRA2aer_OCPHILIC','MERRA2aer_BCPHOBIC','MERRA2aer_BCPHILIC',...
%    'GEOSChem_BC','GEOSChem_NH4','GEOSChem_NO2','GEOSChem_OA','GEOSChem_PM25','GEOSChem_SO4','GEOSChem_O3','GEOSChem_PM25Scaling','GEOSChem_NO2Scaling','GEOSChem_O3Scaling',...
%    'CAMS_NO2'};
VARNAME_LIST_OMI = {'GEOSChem_BC','GEOSChem_NH4','GEOSChem_NO2','GEOSChem_OA','GEOSChem_PM25','GEOSChem_SO4','GEOSChem_O3','GEOSChem_PM25Scaling','GEOSChem_NO2Scaling','GEOSChem_O3Scaling','CAMS_NO2'};
DATA_OPTION_LIST_OMI = {'Mean'};

% nearby monitors
IsNearby = false;
DATA_OPTION_LIST_Nearby = {'Peak2_Lag0','Peak2_Lag1','Peak2_Lag3'};
VARNAME_LIST_Nearby = {'NO2','PM25','Ozone'};

% % % land use variables
IsLandUse = true;

VARNAME_LIST_LandUse = {'RoadDensity_primaryroads1000','RoadDensity_primaryroads10000','RoadDensity_prisecroads1000','RoadDensity_prisecroads10000','RoadDensity_roads1000',...
    'NLCD_Water100','NLCD_Wetlands10000','NLCD_Wetlands100','NLCD_Planted10000','NLCD_Planted100','NLCD_Herbaceous10000','NLCD_Herbaceous100','NLCD_Shrubland10000','NLCD_Shrubland100','NLCD_Barren10000','NLCD_Barren100','NLCD_Developed10000','NLCD_Developed100','NLCD_Water10000',...
    'NLCD_Impervious100','NLCD_Impervious10000',...
    'NLCD_canopy100','NLCD_canopy10000',...
    'Business_Restaurant1000',...
    'TruckRoute_Traffic1000','TruckRoute_Traffic100','TruckRoute_ShortDis1000','TruckRoute_ShortDis100','TruckRoute_Density100','TruckRoute_Density1000',...
    'NightLight',...
    'USElevation_dsc10000','USElevation_max100','USElevation_max10000','USElevation_mea100','USElevation_mea10000','USElevation_med100','USElevation_med10000','USElevation_min100','USElevation_min10000','USElevation_std100','USElevation_std10000','USElevation_bln100','USElevation_bln10000','USElevation_dsc100',...
    'AQRVPM25Site_Region','EPACastNetOzoneSite_Region','EPANO2Site_Region'};
if(2 == EXTRAOPTION)
    VARNAME_LIST_LandUse = {'RoadDensity_primaryroads1000','RoadDensity_primaryroads10000','RoadDensity_prisecroads1000','RoadDensity_prisecroads10000','RoadDensity_roads1000',...
        'NLCD_Water100','NLCD_Wetlands10000','NLCD_Wetlands100','NLCD_Planted10000','NLCD_Planted100','NLCD_Herbaceous10000','NLCD_Herbaceous100','NLCD_Shrubland10000','NLCD_Shrubland100','NLCD_Barren10000','NLCD_Barren100','NLCD_Developed10000','NLCD_Developed100','NLCD_Water10000',...
        'NLCD_Impervious100','NLCD_Impervious10000',...
        'NLCD_canopy100','NLCD_canopy10000',...
        'NightLight',...
        'Business_Restaurant1000'};
end
    
% this is some options may be changed if used as stand alone function; these
% options are overriden if called as a function by Run_ProcessDataANDRunModel
if(~isstruct(SystemPara))
    SystemPara = struct(...
    'NAME',{{'PM25','Ozone','NO2'}},...% air pollutants we will work at
    'SITENAME_TRAIN',{{'AQRVPM25','EPACastNetOzone','EPANO2'}},...% grid cells used for training   
    'GCS','North_America_Equidistant_Conic',...% projection system
    'TrainDataOption','Interpolation',...% we just interpolate aggregate data to input data
    'ProcessDataOption','By-Day');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DO NOT CHANGE ANYTHING BELOW

%% marco-definition
OPTION = SystemPara.ProcessDataOption;
INPUTYEAR = year(INPUTDATE(1));
MAXYEAR = 2016;
MINYEAR = 1999;
if(INPUTYEAR>MAXYEAR)
    INPUTYEAR = INPUTYEAR - MAXYEAR - 1 + MINYEAR;
elseif(INPUTYEAR<MINYEAR)
    INPUTYEAR = INPUTYEAR - MINYEAR + 1 + MAXYEAR;
end

%% path
% DIRAPTH = 'D:/Google Drive/Research/USTemperature/processed_data/DataProcessRecorder/';
% DIRPATH_ROOT = 'D:/Google Drive/Research/USTemperature/processed_data/';
DIRAPTH = '../../processed_data/DataProcessRecorder/';
DIRPATH_ROOT = '../../processed_data/';

Sep = '/';
DIRAPTH = [DIRAPTH,num2str(INPUTYEAR),Sep];
mkdir([DIRPATH_ROOT,SITEDATA_NAME,Sep]);
mkdir(DIRAPTH);

%% create weights or delete lagged terms
% this is used only when called by Run_ProcessDataANDRunModel 
if(1 == EXTRAOPTION)
    SiteName_Train_list = SystemPara.SITENAME_TRAIN;
    for i = 1:length(SiteName_Train_list)
        SiteName_Predict = SITEDATA_NAME;
        SiteName_Train = SiteName_Train_list{i};

        SiteData_Train = LoadData_function([DIRPATH_ROOT,SiteName_Train,Sep,'Location',Sep,SiteName_Train,'Site_',GCS,'.mat'],'SiteData');% the monitor sites we will use for modeling;monitors of interest
        SiteData_Train = SiteData_Train.SiteData;

        SiteData_Predict = LoadData_function([DIRPATH_ROOT,SiteName_Predict,Sep,'Location',Sep,SiteName_Predict,'Site_',GCS,'.mat'],'SiteData');
        SiteData_Predict = SiteData_Predict.SiteData;
        
        Analysis_TwoStep1([DIRPATH_ROOT,SiteName_Predict,Sep,'Temp',Sep],SiteData_Predict,SiteData_Train,SiteName_Predict,SiteName_Train);        
    end
end

%% total output
diary([DIRAPTH,'Run_ProcessDataSummary_',num2str(INPUTYEAR),'.txt'])
FID_Error = fopen([DIRAPTH,'Run_ProcessDataSummary_',num2str(INPUTYEAR),'_Error.txt'],'a');

%% begin to process
% ReturnStatus==1: processing done before;  ReturnStatus==0: still in processing
% ReturnStatus==2, some raw data not available and needs to do it again; ReturnStatus==3, some input data not available has to quit;
TotalStatus = [];
while(true)
    TotalStatus_previous = TotalStatus;
    TotalStatus = [];
    
   
    %% land-cover aka MCD12Q1
    DATA_OPTION_LIST = DATA_OPTION_LIST_MCD12Q1;
    if(IsMCD12Q1)
        for j=1:length(DATA_OPTION_LIST)
            fprintf('%s\n',repmat('*',5));
            
            % is it being processed? if yes, skip this and process nex tvariable!
            FileName = [SITEDATA_NAME,num2str(INPUTYEAR),'MCD12Q1','By_Year',DATA_OPTION_LIST{j}];
            IsProcessFlag = false;
            if(exist([DIRAPTH,FileName,'.mat.part'],'file'))
                list = dir([DIRAPTH,FileName,'.mat.part']);
                if(etime(clock,datevec(list.datenum))>3600)
                    delete([DIRAPTH,FileName,'.mat.part']);
                    fprintf('%s process work started long time ago and still has not completed. we do it again!...\n',FileName);
                    IsProcessFlag = true;
                else
                    fprintf('%s is in the middle processing!...\n',FileName);
                    IsProcessFlag = false;
                end
            else
                IsProcessFlag = true;
            end
            
            % if no, we will process this variable
            if(IsProcessFlag)
                save([DIRAPTH,FileName,'.mat.part'],'Sep');
                fprintf('SiteName:%s,Year:%d,Variable Name:%s,Output format:%s,convolutional option:%s\n',SITEDATA_NAME,INPUTYEAR,'MCD12Q1','By-Year',DATA_OPTION_LIST{j});
                
                % beging to process! if something went wrong return error (ReturnStatus = 2)
                try
                    Status = Interpolate_MCD12Q1_Main(SITEDATA_NAME,INPUTDATE,'MCD12Q1','By-Year',DATA_OPTION_LIST{j},EXTRAOPTION,FORMAT);
                catch exception
                    fprintf('%s\n',exception.message);
                    fprintf(FID_Error,'%s,SiteName:%s,Year:%d,Variable Name:%s,Output format:%s,convolutional option:%s, error message:%s\n',datestr(now),SITEDATA_NAME,INPUTYEAR,'MCD12Q1','By-Year',DATA_OPTION_LIST{j},exception.message);
                    Status = 2;
                end
                % get the return message, get the processing status
                % ReturnStatus==1: processing done before;  ReturnStatus==0: still in processing
                % ReturnStatus==2, some raw data not available and needs to do it again; ReturnStatus==3, some input data not available has to quit;
                TotalStatus = cat(1,TotalStatus,Status);
                fprintf('%s\t%d\n',repmat('*',5),Status);
                delete([DIRAPTH,FileName,'.mat.part']);
            else
                Status = 1;
                TotalStatus = cat(1,TotalStatus,Status);
            end
        end
    end
    pause(2);
    
    %% meteorological variables
    VARNAME_LIST =  VARNAME_LIST_Mete;
    DATA_OPTION_LIST = DATA_OPTION_LIST_Mete;
    if(length(DATA_OPTION_LIST_Mete)==length(DATA_OPTION_LIST))
        if(IsMete)
            for i = 1:length(VARNAME_LIST)
                j = i;
                fprintf('%s\n',repmat('*',5));

                % is it being processed? if yes, skip this and process nex tvariable!
                FileName = [SITEDATA_NAME,num2str(INPUTYEAR),VARNAME_LIST{i},OPTION,DATA_OPTION_LIST{j}];
                IsProcessFlag = false;
                if(exist([DIRAPTH,FileName,'.mat.part'],'file'))
                    list = dir([DIRAPTH,FileName,'.mat.part']);
                    if(etime(clock,datevec(list.datenum))>3600)
                        delete([DIRAPTH,FileName,'.mat.part']);
                        fprintf('%s process work started long time ago and still has not completed. we do it again!...\n',FileName);
                        IsProcessFlag = true;
                    else
                        fprintf('%s is in the middle processing!...\n',FileName);
                        IsProcessFlag = false;
                    end
                else
                    IsProcessFlag = true;
                end

                % if no, we will process this variable
                if(IsProcessFlag)
                    save([DIRAPTH,FileName,'.mat.part'],'Sep');
                    fprintf('SiteName:%s,Year:%d,Variable Name:%s,Output format:%s,convolutional option:%s\n',SITEDATA_NAME,INPUTYEAR,VARNAME_LIST{i},OPTION,DATA_OPTION_LIST{j});

                    % beging to process! if something went wrong return error (ReturnStatus = 2)
                    try
                        Status = Interpolate_3hourReanalysis_Main(SITEDATA_NAME,INPUTDATE,VARNAME_LIST{i},OPTION,DATA_OPTION_LIST{j},EXTRAOPTION,FORMAT);
                    catch exception
                        fprintf('%s\n',exception.message);
                        fprintf(FID_Error,'%s,SiteName:%s,Year:%d,Variable Name:%s,Output format:%s,convolutional option:%s, error message:%s\n',datestr(now),SITEDATA_NAME,INPUTYEAR,VARNAME_LIST{i},OPTION,DATA_OPTION_LIST{j},exception.message);
                        Status = 2;
                    end  
                    % get the return message, get the processing status
                    % ReturnStatus==1: processing done before;  ReturnStatus==0: still in processing
                    % ReturnStatus==2, some raw data not available and needs to do it again; ReturnStatus==3, some input data not available has to quit;
                    TotalStatus = cat(1,TotalStatus,Status);
                    fprintf('%s\t%d\n',repmat('*',5),Status);
                    delete([DIRAPTH,FileName,'.mat.part']);
                else
                    Status = 1;
                    TotalStatus = cat(1,TotalStatus,Status);
                end
            end
        end
    else
        disp('length of DATA_OPTION_LIST_Mete and VARNAME_LIST_Mete shoud be identical!');
    end
    pause(2);
    
    
    %% meteorological variable -- pressure level
    VARNAME_LIST =  VARNAME_LIST_Omega;
    DATA_OPTION_LIST = DATA_OPTION_LIST_Omega;
    if(IsMeteOmega)
        for i = 1:length(VARNAME_LIST)
            for j=1:length(DATA_OPTION_LIST)
                fprintf('%s\n',repmat('*',5));
                
                % is it being processed? if yes, skip this and process nex tvariable!
                FileName = [SITEDATA_NAME,num2str(INPUTYEAR),VARNAME_LIST{i},OPTION,DATA_OPTION_LIST{j}];
                IsProcessFlag = false;
                if(exist([DIRAPTH,FileName,'.mat.part'],'file'))
                    list = dir([DIRAPTH,FileName,'.mat.part']);
                    if(etime(clock,datevec(list.datenum))>3600)
                        delete([DIRAPTH,FileName,'.mat.part']);
                        fprintf('%s process work started long time ago and still has not completed. we do it again!...\n',FileName);
                        IsProcessFlag = true;
                    else
                        fprintf('%s is in the middle processing!...\n',FileName);
                        IsProcessFlag = false;
                    end
                else
                    IsProcessFlag = true;
                end
                
                % if no, we will process this variable
                if(IsProcessFlag)
                    save([DIRAPTH,FileName,'.mat.part'],'Sep');
                    fprintf('SiteName:%s,Year:%d,Variable Name:%s,Output format:%s,convolutional option:%s\n',SITEDATA_NAME,INPUTYEAR,VARNAME_LIST{i},OPTION,DATA_OPTION_LIST{j});
                    
                    % beging to process! if something went wrong return error (ReturnStatus = 2)
                    try
                        Status = Interpolate_3hourPressureLevelReanalysis_Main(SITEDATA_NAME,INPUTDATE,VARNAME_LIST{i},OPTION,DATA_OPTION_LIST{j},EXTRAOPTION,FORMAT);
                    catch exception
                        fprintf('%s\n',exception.message);
                        fprintf(FID_Error,'%s,SiteName:%s,Year:%d,Variable Name:%s,Output format:%s,convolutional option:%s, error message:%s\n',datestr(now),SITEDATA_NAME,INPUTYEAR,VARNAME_LIST{i},OPTION,DATA_OPTION_LIST{j},exception.message);
                        Status = 2;
                    end
                    % get the return message, get the processing status
                    % ReturnStatus==1: processing done before;  ReturnStatus==0: still in processing
                    % ReturnStatus==2, some raw data not available and needs to do it again; ReturnStatus==3, some input data not available has to quit;
                    TotalStatus = cat(1,TotalStatus,Status);
                    fprintf('%s\t%d\n',repmat('*',5),Status);
                    delete([DIRAPTH,FileName,'.mat.part']);
                else
                    Status = 1;
                    TotalStatus = cat(1,TotalStatus,Status);
                end
            end
        end
    end
    pause(2);
    
    %% MOD13A1
    VARNAME_LIST =  VARNAME_LIST_MOD13A2;
    DATA_OPTION_LIST = DATA_OPTION_LIST_MOD13A2;
    if(IsMOD13A2_MOD09A1 && INPUTYEAR>=2000)
        for i = 1:length(VARNAME_LIST)
            for j=1:length(DATA_OPTION_LIST)
                fprintf('%s\n',repmat('*',5));

                % is it being processed? if yes, skip this and process nex tvariable!
                FileName = [SITEDATA_NAME,num2str(INPUTYEAR),VARNAME_LIST{i},OPTION,DATA_OPTION_LIST{j}];
                IsProcessFlag = false;
                if(exist([DIRAPTH,FileName,'.mat.part'],'file'))
                    list = dir([DIRAPTH,FileName,'.mat.part']);
                    if(etime(clock,datevec(list.datenum))>3600)
                        delete([DIRAPTH,FileName,'.mat.part']);
                        fprintf('%s process work started long time ago and still has not completed. we do it again!...\n',FileName);
                        IsProcessFlag = true;
                    else
                        fprintf('%s is in the middle processing!...\n',FileName);
                        IsProcessFlag = false;
                    end
                else
                    IsProcessFlag = true;
                end

                % if no, we will process this variable
                if(IsProcessFlag)
                    save([DIRAPTH,FileName,'.mat.part'],'Sep');
                    fprintf('SiteName:%s,Year:%d,Variable Name:%s,Output format:%s,convolutional option:%s\n',SITEDATA_NAME,INPUTYEAR,VARNAME_LIST{i},OPTION,DATA_OPTION_LIST{j});

                    % beging to process! if something went wrong return error (ReturnStatus = 2)
                    try
                        Status = Interpolate_Generic_Main(SITEDATA_NAME,INPUTDATE,VARNAME_LIST{i},OPTION,DATA_OPTION_LIST{j},EXTRAOPTION,FORMAT);
                    catch exception
                        fprintf('%s\n',exception.message);
                        fprintf(FID_Error,'%s,SiteName:%s,Year:%d,Variable Name:%s,Output format:%s,convolutional option:%s, error message:%s\n',datestr(now),SITEDATA_NAME,INPUTYEAR,VARNAME_LIST{i},OPTION,DATA_OPTION_LIST{j},exception.message);
                        Status = 2;
                    end
                    % get the return message, get the processing status
                    % ReturnStatus==1: processing done before;  ReturnStatus==0: still in processing
                    % ReturnStatus==2, some raw data not available and needs to do it again; ReturnStatus==3, some input data not available has to quit;
                    TotalStatus = cat(1,TotalStatus,Status);
                    fprintf('%s\t%d\n',repmat('*',5),Status);
                    delete([DIRAPTH,FileName,'.mat.part']);
                else
                    Status = 1;
                    TotalStatus = cat(1,TotalStatus,Status);
                end

            end
        end
    end
    pause(2);
    
    
    %% MOD11A1
    DATA_OPTION_LIST = DATA_OPTION_LIST_MOD11A1;
    VARNAME_LIST = {'MOD11A1'};
    if(IsMOD11A1)
        for i = 1:length(VARNAME_LIST)
            for j=1:length(DATA_OPTION_LIST)
                fprintf('%s\n',repmat('*',5));
                
                % is it being processed? if yes, skip this and process nex tvariable!
                FileName = [SITEDATA_NAME,num2str(INPUTYEAR),VARNAME_LIST{i},OPTION];
                IsProcessFlag = false;
                if(exist([DIRAPTH,FileName,'.mat.part'],'file'))
                    list = dir([DIRAPTH,FileName,'.mat.part']);
                    if(etime(clock,datevec(list.datenum))>3600)
                        delete([DIRAPTH,FileName,'.mat.part']);
                        fprintf('%s process work started long time ago and still has not completed. we do it again!...\n',FileName);
                        IsProcessFlag = true;
                    else
                        fprintf('%s is in the middle processing!...\n',FileName);
                        IsProcessFlag = false;
                    end
                else
                    IsProcessFlag = true;
                end
                
                % if no, we will process this variable
                if(IsProcessFlag)
                    save([DIRAPTH,FileName,'.mat.part'],'Sep');
                    fprintf('SiteName:%s,Year:%d,Variable Name:%s,Output format:%s,convolutional option:%s\n',SITEDATA_NAME,INPUTYEAR,VARNAME_LIST{i},OPTION,nan);
                    
                    % beging to process! if something went wrong return error (ReturnStatus = 2)
                    try
                        Status = Interpolate_MOD11A1_Main(SITEDATA_NAME,INPUTDATE,VARNAME_LIST{i},OPTION,DATA_OPTION_LIST{j},EXTRAOPTION,FORMAT);
                    catch exception
                        fprintf('%s\n',exception.message);
                        fprintf(FID_Error,'%s,SiteName:%s,Year:%d,Variable Name:%s,Output format:%s,convolutional option:%s, error message:%s\n',datestr(now),SITEDATA_NAME,INPUTYEAR,VARNAME_LIST{i},OPTION,DATA_OPTION_LIST{j},exception.message);
                        Status = 2;
                    end
                    % get the return message, get the processing status
                    % ReturnStatus==1: processing done before;  ReturnStatus==0: still in processing
                    % ReturnStatus==2, some raw data not available and needs to do it again; ReturnStatus==3, some input data not available has to quit;
                    TotalStatus = cat(1,TotalStatus,Status);
                    fprintf('%s\t%d\n',repmat('*',5),Status);
                    delete([DIRAPTH,FileName,'.mat.part']);
                else
                    Status = 1;
                    TotalStatus = cat(1,TotalStatus,Status);
                end
            end
        end
    end
    pause(2);
    
    %% MAIAC US 1km AOD
    DATA_OPTION_LIST = DATA_OPTION_LIST_MAIACUSAOD;
    VARNAME_LIST = VARNAME_LIST_MAIACUSAOD;
    if(IsMAIACUSAOD)
        for i = 1:length(VARNAME_LIST)
            for j=1:length(DATA_OPTION_LIST)
                fprintf('%s\n',repmat('*',5));
                
                % is it being processed? if yes, skip this and process nex tvariable!
                FileName = [SITEDATA_NAME,num2str(INPUTYEAR),VARNAME_LIST{i},OPTION];
                IsProcessFlag = false;
                if(exist([DIRAPTH,FileName,'.mat.part'],'file'))
                    list = dir([DIRAPTH,FileName,'.mat.part']);
                    if(etime(clock,datevec(list.datenum))>3600)
                        delete([DIRAPTH,FileName,'.mat.part']);
                        fprintf('%s process work started long time ago and still has not completed. we do it again!...\n',FileName);
                        IsProcessFlag = true;
                    else
                        fprintf('%s is in the middle processing!...\n',FileName);
                        IsProcessFlag = false;
                    end
                else
                    IsProcessFlag = true;
                end
                
                % if no, we will process this variable
                if(IsProcessFlag)
                    save([DIRAPTH,FileName,'.mat.part'],'Sep');
                    fprintf('SiteName:%s,Year:%d,Variable Name:%s,Output format:%s,convolutional option:%s\n',SITEDATA_NAME,INPUTYEAR,VARNAME_LIST{i},OPTION,nan);
                    
                    % beging to process! if something went wrong return error (ReturnStatus = 2)
                    try
                        Status = Interpolate_MOD11A1_Main(SITEDATA_NAME,INPUTDATE,VARNAME_LIST{i},OPTION,DATA_OPTION_LIST{j},EXTRAOPTION,FORMAT);
                    catch exception
                        fprintf('%s\n',exception.message);
                        fprintf(FID_Error,'%s,SiteName:%s,Year:%d,Variable Name:%s,Output format:%s,convolutional option:%s, error message:%s\n',datestr(now),SITEDATA_NAME,INPUTYEAR,VARNAME_LIST{i},OPTION,DATA_OPTION_LIST{j},exception.message);
                        Status = 2;
                    end
                    % get the return message, get the processing status
                    % ReturnStatus==1: processing done before;  ReturnStatus==0: still in processing
                    % ReturnStatus==2, some raw data not available and needs to do it again; ReturnStatus==3, some input data not available has to quit;
                    TotalStatus = cat(1,TotalStatus,Status);
                    fprintf('%s\t%d\n',repmat('*',5),Status);
                    delete([DIRAPTH,FileName,'.mat.part']);
                else
                    Status = 1;
                    TotalStatus = cat(1,TotalStatus,Status);
                end
            end
        end
    end
    pause(2);
    
    %% MAIAC US 1km AOD
    DATA_OPTION_LIST = DATA_OPTION_LIST_MAIACUScosVZA;
    VARNAME_LIST = VARNAME_LIST_MAIACUScosVZA;
    if(IsMAIACUScosVZA)
        for i = 1:length(VARNAME_LIST)
            for j=1:length(DATA_OPTION_LIST)
                fprintf('%s\n',repmat('*',5));
                
                % is it being processed? if yes, skip this and process nex tvariable!
                FileName = [SITEDATA_NAME,num2str(INPUTYEAR),VARNAME_LIST{i},OPTION];
                IsProcessFlag = false;
                if(exist([DIRAPTH,FileName,'.mat.part'],'file'))
                    list = dir([DIRAPTH,FileName,'.mat.part']);
                    if(etime(clock,datevec(list.datenum))>3600)
                        delete([DIRAPTH,FileName,'.mat.part']);
                        fprintf('%s process work started long time ago and still has not completed. we do it again!...\n',FileName);
                        IsProcessFlag = true;
                    else
                        fprintf('%s is in the middle processing!...\n',FileName);
                        IsProcessFlag = false;
                    end
                else
                    IsProcessFlag = true;
                end
                
                % if no, we will process this variable
                if(IsProcessFlag)
                    save([DIRAPTH,FileName,'.mat.part'],'Sep');
                    fprintf('SiteName:%s,Year:%d,Variable Name:%s,Output format:%s,convolutional option:%s\n',SITEDATA_NAME,INPUTYEAR,VARNAME_LIST{i},OPTION,nan);
                    
                    % beging to process! if something went wrong return error (ReturnStatus = 2)
                    try
                        Status = Interpolate_MOD11A1_Main(SITEDATA_NAME,INPUTDATE,VARNAME_LIST{i},OPTION,DATA_OPTION_LIST{j},EXTRAOPTION,FORMAT);
                    catch exception
                        fprintf('%s\n',exception.message);
                        fprintf(FID_Error,'%s,SiteName:%s,Year:%d,Variable Name:%s,Output format:%s,convolutional option:%s, error message:%s\n',datestr(now),SITEDATA_NAME,INPUTYEAR,VARNAME_LIST{i},OPTION,DATA_OPTION_LIST{j},exception.message);
                        Status = 2;
                    end
                    % get the return message, get the processing status
                    % ReturnStatus==1: processing done before;  ReturnStatus==0: still in processing
                    % ReturnStatus==2, some raw data not available and needs to do it again; ReturnStatus==3, some input data not available has to quit;
                    TotalStatus = cat(1,TotalStatus,Status);
                    fprintf('%s\t%d\n',repmat('*',5),Status);
                    delete([DIRAPTH,FileName,'.mat.part']);
                else
                    Status = 1;
                    TotalStatus = cat(1,TotalStatus,Status);
                end
            end
        end
    end
    pause(2);
    
    %% OMI
    DATA_OPTION_LIST = DATA_OPTION_LIST_OMI;
    VARNAME_LIST = VARNAME_LIST_OMI;
    if(IsOMI)
        for i = 1:length(VARNAME_LIST)
            for j=1:length(DATA_OPTION_LIST)
                fprintf('%s\n',repmat('*',5));
                
                % is it being processed? if yes, skip this and process nex tvariable!
                FileName = [SITEDATA_NAME,num2str(INPUTYEAR),VARNAME_LIST{i},OPTION];
                IsProcessFlag = false;
                if(exist([DIRAPTH,FileName,'.mat.part'],'file'))
                    list = dir([DIRAPTH,FileName,'.mat.part']);
                    if(etime(clock,datevec(list.datenum))>3600)
                        delete([DIRAPTH,FileName,'.mat.part']);
                        fprintf('%s process work started long time ago and still has not completed. we do it again!...\n',FileName);
                        IsProcessFlag = true;
                    else
                        fprintf('%s is in the middle processing!...\n',FileName);
                        IsProcessFlag = false;
                    end
                else
                    IsProcessFlag = true;
                end
                
                % if no, we will process this variable
                if(IsProcessFlag)
                    save([DIRAPTH,FileName,'.mat.part'],'Sep');
                    fprintf('SiteName:%s,Year:%d,Variable Name:%s,Output format:%s,convolutional option:%s\n',SITEDATA_NAME,INPUTYEAR,VARNAME_LIST{i},OPTION,nan);
                    
                    % beging to process! if something went wrong return error (ReturnStatus = 2)
                    try
                        Status = Interpolate_OMI_Main(SITEDATA_NAME,INPUTDATE,VARNAME_LIST{i},OPTION,DATA_OPTION_LIST{j},EXTRAOPTION,FORMAT);
                    catch exception
                        fprintf('%s\n',exception.message);
                        fprintf(FID_Error,'%s,SiteName:%s,Year:%d,Variable Name:%s,Output format:%s,convolutional option:%s, error message:%s\n',datestr(now),SITEDATA_NAME,INPUTYEAR,VARNAME_LIST{i},OPTION,DATA_OPTION_LIST{j},exception.message);
                        Status = 2;
                    end
                    % get the return message, get the processing status
                    % ReturnStatus==1: processing done before;  ReturnStatus==0: still in processing
                    % ReturnStatus==2, some raw data not available and needs to do it again; ReturnStatus==3, some input data not available has to quit;
                    TotalStatus = cat(1,TotalStatus,Status);
                    fprintf('%s\t%d\n',repmat('*',5),Status);
                    delete([DIRAPTH,FileName,'.mat.part']);
                else
                    Status = 1;
                    TotalStatus = cat(1,TotalStatus,Status);
                end
            end
        end
    end
    pause(2);
    
    %% nearby monitors
    DATA_OPTION_LIST = DATA_OPTION_LIST_Nearby;
    VARNAME_LIST = VARNAME_LIST_Nearby;
    if(IsNearby)
        for i = 1:length(VARNAME_LIST)
            for j=1:length(DATA_OPTION_LIST)
                fprintf('%s\n',repmat('*',5));
                
                % is it being processed? if yes, skip this and process nex tvariable!
                FileName = [SITEDATA_NAME,num2str(INPUTYEAR),VARNAME_LIST{i},OPTION];
                IsProcessFlag = false;
                if(exist([DIRAPTH,FileName,'.mat.part'],'file'))
                    list = dir([DIRAPTH,FileName,'.mat.part']);
                    if(etime(clock,datevec(list.datenum))>3600)
                        delete([DIRAPTH,FileName,'.mat.part']);
                        fprintf('%s process work started long time ago and still has not completed. we do it again!...\n',FileName);
                        IsProcessFlag = true;
                    else
                        fprintf('%s is in the middle processing!...\n',FileName);
                        IsProcessFlag = false;
                    end
                else
                    IsProcessFlag = true;
                end
                
                % if no, we will process this variable
                if(IsProcessFlag)
                    save([DIRAPTH,FileName,'.mat.part'],'Sep');
                    fprintf('SiteName:%s,Year:%d,Variable Name:%s,Output format:%s,convolutional option:%s\n',SITEDATA_NAME,INPUTYEAR,VARNAME_LIST{i},DATA_OPTION_LIST{j},nan);
                    
                    % beging to process! if something went wrong return error (ReturnStatus = 2)
                    try
                        Status = Interpolate_NearbyMonitor_Main(SITEDATA_NAME,INPUTDATE,VARNAME_LIST{i},OPTION,DATA_OPTION_LIST{j},EXTRAOPTION,FORMAT);
                    catch exception
                        fprintf('%s\n',exception.message);
                        fprintf(FID_Error,'%s,SiteName:%s,Year:%d,Variable Name:%s,Output format:%s,convolutional option:%s, error message:%s\n',datestr(now),SITEDATA_NAME,INPUTYEAR,VARNAME_LIST{i},OPTION,DATA_OPTION_LIST{j},exception.message);
                        Status = 2;
                    end
                    % get the return message, get the processing status
                    % ReturnStatus==1: processing done before;  ReturnStatus==0: still in processing
                    % ReturnStatus==2, some raw data not available and needs to do it again; ReturnStatus==3, some input data not available has to quit;
                    TotalStatus = cat(1,TotalStatus,Status);
                    fprintf('%s\t%d\n',repmat('*',5),Status);
                    delete([DIRAPTH,FileName,'.mat.part']);
                else
                    Status = 1;
                    TotalStatus = cat(1,TotalStatus,Status);
                end
            end
        end
    end
    pause(2);
    
    %% LandUse variables from road density and NLCD
    DATA_OPTION_LIST = {'test'};
    VARNAME_LIST = VARNAME_LIST_LandUse;
    if(IsLandUse)
        for i = 1:length(VARNAME_LIST)
            for j=1:length(DATA_OPTION_LIST)
                fprintf('%s\n',repmat('*',5));
                
                % is it being processed? if yes, skip this and process nex tvariable!
                FileName = [SITEDATA_NAME,num2str(INPUTYEAR),VARNAME_LIST{i},OPTION];
                IsProcessFlag = false;
                if(exist([DIRAPTH,FileName,'.mat.part'],'file'))
                    list = dir([DIRAPTH,FileName,'.mat.part']);
                    if(etime(clock,datevec(list.datenum))>3600)
                        delete([DIRAPTH,FileName,'.mat.part']);
                        fprintf('%s process work started long time ago and still has not completed. we do it again!...\n',FileName);
                        IsProcessFlag = true;
                    else
                        fprintf('%s is in the middle processing!...\n',FileName);
                        IsProcessFlag = false;
                    end
                else
                    IsProcessFlag = true;
                end
                
                % if no, we will process this variable
                if(IsProcessFlag)
                    save([DIRAPTH,FileName,'.mat.part'],'Sep');
                    fprintf('SiteName:%s,Year:%d,Variable Name:%s,Output format:%s,convolutional option:%s\n',SITEDATA_NAME,INPUTYEAR,VARNAME_LIST{i},OPTION,nan);
                    
                    % beging to process! if something went wrong return error (ReturnStatus = 2)
                    try
                        Status = Interpolate_LandUse_Main(SITEDATA_NAME,INPUTDATE,VARNAME_LIST{i},OPTION,DATA_OPTION_LIST{j},EXTRAOPTION,FORMAT);
                    catch exception
                        fprintf('%s\n',exception.message);
                        fprintf(FID_Error,'%s,SiteName:%s,Year:%d,Variable Name:%s,Output format:%s,convolutional option:%s, error message:%s\n',datestr(now),SITEDATA_NAME,INPUTYEAR,VARNAME_LIST{i},OPTION,DATA_OPTION_LIST{j},exception.message);
                        Status = 2;
                    end
                    % get the return message, get the processing status
                    % ReturnStatus==1: processing done before;  ReturnStatus==0: still in processing
                    % ReturnStatus==2, some raw data not available and needs to do it again; ReturnStatus==3, some input data not available has to quit;
                    TotalStatus = cat(1,TotalStatus,Status);
                    fprintf('%s\t%d\n',repmat('*',5),Status);
                    delete([DIRAPTH,FileName,'.mat.part']);
                else
                    Status = 1;
                    TotalStatus = cat(1,TotalStatus,Status);
                end
            end
        end
    end
    pause(2);
    
    % if all data have been processed, return;
    if(all(TotalStatus == 1))
        ReturnStatus = 1;
        fclose(FID_Error);
        return;
    end
       
    % if there are some issue with aggregate data cannot be fixed, and no data are being procssed return
    if(sum(TotalStatus == 2) == sum(TotalStatus_previous == 2) && sum(TotalStatus == 0) == 0)
        ReturnStatus = 2;
        fclose(FID_Error);
        return;
    end
end

end


%% one important step of the air pollution modeling is what called "two-step" 
% modeling. It means (1) we fitted the neural network (it can be any
% machine learning algorithm or even regression model) and obtained
% predicted air pollution levels at each monitoring site; then (2) we
% calculated the averaged predicted air pollution levels from neigbouring
% days and nearby monitoring site; and (3) we added those average values at
% additional input variables and fitted the neural network again. Since we
% fitted the neural network twice, we call it "two-step" modeling.
% this function is not supposed to run in a stand alone way

%% parameter
% Data: predicted air pollution level at monitoring sites, a huge matrix;
% SiteData_Predict: the grid cell we will calculate average predicted values at
% SiteData_Data: the grid cell we already have fitted model at;
% SiteName_Predict:  the grid cell name we will calculate average predicted values at
% SiteName_Data: the grid cell name we already have fitted model at;
% EnvironPara: some macro definition
% OPTION: CreateWeight: just create weight; UPDATE: create files to do
% two-step modeling

%% version history
% version 2
% SiteData_Data is the Site Data we have data at, data is the variable "Data"[N_Day,N_Site] fomrat;
% SiteData_Predict is the Site Data will make prediction at;
% if OPTION == CreateWeight --> just create and save those weight matrix,
% which takes a lot of memory
%  version 3
% not long need to update initially
%  version 4
% change weight again
% change name to Analysis_TwoStep

%% code
function Result = Analysis_TwoStep1(DIRPATH,SiteData_Predict,SiteData_Train,SiteName_Predict,SiteName_Train)
disp('Analysis_TwoStep');
Sep = '/';
mkdir(DIRPATH);
%% estimation-based spatial lagged terms
%%after getting better estimation, update those terms recursively

%% just create those weights -- no threshod value/with threshold values
% these weight matries are used to calculate averaged predicted values from
% nearby monitors. Here we used matrix multiplication to do calculate such
% averaged values
    if(iscell(SiteData_Predict))
        Lon_Target = cell2mat(SiteData_Predict(:,3));
        Lat_Target = cell2mat(SiteData_Predict(:,2));
    elseif(istable(SiteData_Predict))
        Lon_Target = SiteData_Predict.Lon;
        Lat_Target = SiteData_Predict.Lat;
    else
        Lon_Target = SiteData_Predict(:,3);
        Lat_Target = SiteData_Predict(:,2);
    end
    N_Traget = length(Lon_Target);
    
    %% raw data location list
    if(iscell(SiteData_Train))
        Lon_Query = cell2mat(SiteData_Train(:,3));
        Lat_Query = cell2mat(SiteData_Train(:,2));
    elseif(istable(SiteData_Train))
        Lon_Query = SiteData_Train.Lon;
        Lat_Query = SiteData_Train.Lat;
    else
        Lon_Query = SiteData_Train(:,3);
        Lat_Query = SiteData_Train(:,2);
    end
    N_Query = length(Lon_Query);
    
    
    THRESHOLD = Inf;
    N_Neighbour = 200;
    
    %% without threshold values
    if(exist([DIRPATH,'SpatialLaggedWeightPeak41_',SiteName_Train,'_',SiteName_Predict,'.h5'],'file') &&...
            exist([DIRPATH,'SpatialLaggedWeightPeak42_',SiteName_Train,'_',SiteName_Predict,'.h5'],'file') &&...
            exist([DIRPATH,'SpatialLaggedWeightPeak43_',SiteName_Train,'_',SiteName_Predict,'.h5'],'file'))
        disp('Weight matrixes Peak2 have been processed!'); 
    else
        %% first find neighbours.
        CurrentFile1 = [DIRPATH,Sep,'SpatialLaggedNeighbourAdjacentList_',['Peak','_Thres',num2str(THRESHOLD)],'_',SiteName_Train,'_',SiteName_Predict,'.mat'];
        if(exist(CurrentFile1,'file'))
            Result = LoadData_function(CurrentFile1);
            d = Result.d;
            n = Result.n;
        else
            fprintf('creating neighourhood files...\n');
            [n,d] = knnsearch([Lat_Query,Lon_Query],[Lat_Target,Lon_Target],'k',N_Neighbour); 
            save(CurrentFile1,'n','d','-v7.3');
        end
        
        d(d(:,1)==0,1)=0;%
        Index_d = d<THRESHOLD;
        %% version 4, power 1
        if(~exist([DIRPATH,'SpatialLaggedWeightPeak41_',SiteName_Train,'_',SiteName_Predict,'.h5'],'file'))     
            d(Index_d)=1./d(Index_d);%1/d;
            d(isinf(d)) = 0;
            d(~Index_d)=0;
            Weight = sparse(reshape(n,[N_Traget*N_Neighbour,1]),repmat((1:1:N_Traget)',[N_Neighbour,1]),reshape(d,[N_Traget*N_Neighbour,1]),N_Query,N_Traget);
            Sum_Weight_matrix = sum(Weight,1);
            Weight = full(Weight./repmat(Sum_Weight_matrix,[size(SiteData_Train,1),1]));
            %save([DIRPATH,'SpatialLaggedWeightPeak41_',SiteName_Train,'_',SiteName_Predict,'.mat'],'Weight','Sum_Weight_matrix','-v7.3');
            h5create([DIRPATH,'SpatialLaggedWeightPeak41_',SiteName_Train,'_',SiteName_Predict,'.h5'],'/Weight',size(Weight),'Datatype','double', 'ChunkSize',size(Weight),'Deflate',9);
            h5write([DIRPATH,'SpatialLaggedWeightPeak41_',SiteName_Train,'_',SiteName_Predict,'.h5'],'/Weight',Weight);
        end
        
        %% version 4, power 2
        if(~exist([DIRPATH,'SpatialLaggedWeightPeak42_',SiteName_Train,'_',SiteName_Predict,'.h5'],'file'))     
            d(Index_d)=1./power(d(Index_d),2);%1/d^2;
            d(isinf(d)) = 0;
            d(~Index_d)=0;
            Weight = sparse(reshape(n,[N_Traget*N_Neighbour,1]),repmat((1:1:N_Traget)',[N_Neighbour,1]),reshape(d,[N_Traget*N_Neighbour,1]),N_Query,N_Traget);
            Sum_Weight_matrix = sum(Weight,1);
            Weight = full(Weight./repmat(Sum_Weight_matrix,[size(SiteData_Train,1),1]));
            %save([DIRPATH,'SpatialLaggedWeightPeak42_',SiteName_Train,'_',SiteName_Predict,'.mat'],'Weight','Sum_Weight_matrix','-v7.3');
            h5create([DIRPATH,'SpatialLaggedWeightPeak42_',SiteName_Train,'_',SiteName_Predict,'.h5'],'/Weight',size(Weight),'Datatype','double', 'ChunkSize',size(Weight),'Deflate',9);
            h5write([DIRPATH,'SpatialLaggedWeightPeak42_',SiteName_Train,'_',SiteName_Predict,'.h5'],'/Weight',Weight);
        end
        
                
        %% version 4, power 3
        if(~exist([DIRPATH,'SpatialLaggedWeightPeak43_',SiteName_Train,'_',SiteName_Predict,'.h5'],'file'))     
            d(Index_d)=1./power(d(Index_d),3);%1/d^3;
            d(isinf(d)) = 0;
            d(~Index_d)=0;
            Weight = sparse(reshape(n,[N_Traget*N_Neighbour,1]),repmat((1:1:N_Traget)',[N_Neighbour,1]),reshape(d,[N_Traget*N_Neighbour,1]),N_Query,N_Traget);
            Sum_Weight_matrix = sum(Weight,1);
            Weight = full(Weight./repmat(Sum_Weight_matrix,[size(SiteData_Train,1),1]));
            %save([DIRPATH,'SpatialLaggedWeightPeak43_',SiteName_Train,'_',SiteName_Predict,'.mat'],'Weight','Sum_Weight_matrix','-v7.3');
            h5create([DIRPATH,'SpatialLaggedWeightPeak43_',SiteName_Train,'_',SiteName_Predict,'.h5'],'/Weight',size(Weight),'Datatype','double', 'ChunkSize',size(Weight),'Deflate',9);
            h5write([DIRPATH,'SpatialLaggedWeightPeak43_',SiteName_Train,'_',SiteName_Predict,'.h5'],'/Weight',Weight);
        end

    end  
end
