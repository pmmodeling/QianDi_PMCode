%% this is the main function to do (1) model training; (2) model prediction
% this code is able to train different air pollutants; (3) create assembled
% data set;

% requirement: (1) produce input data at grid cells of interests; 
% (2) get VariableList_[IDNUM].csv ready, IDNUM is the model ID number to
% distinguish different models;

%% version history
% version2: some changes in reading data (using date not string);
% could be used for reading grid data;
% some changes in fitting the model, could be used in prediction
% version3:
% I had to create version 3 to perserve some testing code in version 2
%  version3.2:
% same as version 3.1; can do 10-fold cross-valiation now; read variable list from external file
% function CalibrateGC_RunAOD_3_2(IDNUM,TIME_ID,OPTION)
% 2016-08-29: change name from CalibrateGC_RunAOD_3_2 to Analysis_RunModel

%% parameters
% Region_ID: index of region to be used for prediction; 
% this parameter is not used in model training, simply set it to be 1;

% IDNUM: model ID number; for example, if you set IDNUM = 9991, then you
% should prepare file called VariableList_9991.csv to specify input
% variables;

% TIME_ID: which year are your fitting?

% OPTION: cross-validation or fitting on the 100% data set? one-step
% modeling or two step modeling? CV: cross-validation; All: training with
% 100% data; round1: one-step modeling; round3: two-step modeling; 
% so "All_round3" means fitting 100% data with two-step modeling strategy;

% NAME: air pollutant to be modelled; only put 'PM25', 'Ozone', 'NO2';
% other air pollutants are under testing;
% ImputationOption = 'AllRecords';
% % ImputationOption = 'CompleteRecords';
% % % ImputationOption = 'MissingRecords';

% ImputationOption: imputation strategy:  AllRecords = doing imputation and
% running on the all available data set; CompleteRecords = do not run
% imputation, run on complete data only; MissingRecords = doing imputation
% but only on imputed data;

% SystemPara: other system parameter, used only when called by
% Run_ProcessDataANDRunModel to do prediction; leave it blank when training
% the model;

%% return value
% ReturnStatus==0: still in processing
% ReturnStatus==1: processing done before; 
% ReturnStatus==2, some raw data not available and needs to do it again;
% ReturnStatus==3, some input data not available has to quit;

% examples:
% Analysis_RunModel(1,99941,2008,'All_round3','PM25','AllRecords','');
% Analysis_RunModel(1,99941,2008,'CV_round3','PM25','AllRecords','');
% Analysis_RunModel(1,99941,2008,'All_round3','Ozone','AllRecords','');
% Analysis_RunModel(1,99941,2008,'CV_round3','Ozone','AllRecords','');
% Analysis_RunModel(1,99941,2008,'All_round3','NO2','AllRecords','');
% Analysis_RunModel(1,99941,2008,'CV_round3','NO2','AllRecords','');

%% code
function ReturnStatue = Analysis_RunModel(Region_ID,IDNUM,TIME_ID,OPTION,NAME,ImputationOption,SystemPara)

% change some of the setting
if(strcmp(SystemPara,''))
    SystemPara = struct(...
    'DIRPATH_ROOT','../../processed_data/',...
    'CODEPATH','./',...%'C:\Users\qid335\Documents\GitHub\AODModelWork2\',...
    'Sep','/',...
    'IDNUM',IDNUM,...
    'GCS','North_America_Equidistant_Conic',...
    'TIME_ID',TIME_ID,...
    ...
    'TESTING',0,...%%% testing mode to run faster...
    'MULTICORE',1,...%%% use multicore? MULTICORE is the number of cores
    'ModelName','NeuralNetwork',...%%% Machine Learning Method, 'NeuralNetwork' or 'Lasso'
    'FIRSTYEAR',TIME_ID,...
    'ENDYEAR',TIME_ID);
end

%% set the imputation model for different year
if(TIME_ID>2013)
    ImputationOption_Model = 'Model1';% does NOT include CMAQ
else
    ImputationOption_Model = 'Model2';% include CMAQ
end
if(strcmp(NAME,'MeanTemperature') || strcmp(NAME,'MinTemperature') || strcmp(NAME,'MaxTemperature'))
    ImputationOption_Model = 'Model3';
end


%% some marco-definition
DIRPATH_ROOT = SystemPara.DIRPATH_ROOT;
CODEPATH =  SystemPara.CODEPATH;
Sep = SystemPara.Sep;

FIRSTYEAR = SystemPara.FIRSTYEAR;
ENDYEAR = SystemPara.ENDYEAR;


%% some global definition
global TESTING MULTICORE

%%% save pictures??
IsSaveCheckingFigure = false;

%%% testing mode to run faster...
TESTING = SystemPara.TESTING; % if true, the neural network only runs several times;

%%% use multicore? MULTICORE is the number of cores
MULTICORE = SystemPara.MULTICORE;

%%% Machine Learning Method
ModelName = SystemPara.ModelName;

%%% how many neighouring points should be considered
N_Neighbour = 100;

%%% record R2 in every loop
IsRecord = true;

%%% projection system
GCS = SystemPara.GCS;


%% specify the grid cells name to do model training
if(~isempty(strfind(OPTION,'CV'))||~isempty(strfind(OPTION,'All')))   
    if(strcmp(NAME,'MeanTemperature'))
        VERSION_MODEL = 'YanWangTemperature';%main site we are working on
        SITENAME_MODEL = 'YanWangTemperature';

        VERSION_DATA = 'YanWangTemperature';%main site the data files have on
        SITENAME_DATA = 'YanWangTemperature';
    elseif(strcmp(NAME,'MinTemperature'))
        VERSION_MODEL = 'YanWangTemperature';%main site we are working on
        SITENAME_MODEL = 'YanWangTemperature';

        VERSION_DATA = 'YanWangTemperature';%main site the data files have on
        SITENAME_DATA = 'YanWangTemperature';
    elseif(strcmp(NAME,'MaxTemperature'))
        VERSION_MODEL = 'YanWangTemperature';%main site we are working on
        SITENAME_MODEL = 'YanWangTemperature';

        VERSION_DATA = 'YanWangTemperature';%main site the data files have on
        SITENAME_DATA = 'YanWangTemperature';
    elseif(strcmp(NAME,'PM25'))
        VERSION_MODEL = 'AQRVPM25';
        SITENAME_MODEL = 'AQRVPM25';

        VERSION_DATA = 'AQRVPM25';%main site the data files have on
        SITENAME_DATA = 'AQRVPM25';
    elseif(strcmp(NAME,'NO2'))
        VERSION_MODEL = 'EPANO2';%main site we are working on
        SITENAME_MODEL = 'EPANO2';

        VERSION_DATA = 'EPANO2';%main site the data files have on
        SITENAME_DATA = 'EPANO2';
    elseif(strcmp(NAME,'Ozone'))
        VERSION_MODEL = 'EPACastNetOzone';%main site we are working on
        SITENAME_MODEL = 'EPACastNetOzone';

        VERSION_DATA = 'EPACastNetOzone';%main site the data files have on
        SITENAME_DATA = 'EPACastNetOzone';
    end

    
    VERSION_PREDICT = VERSION_MODEL;% at what grid cell we are making prediction; for all/cv, this is equvalient to VERSION_MODEL
    SITENAME_PREDICT = SITENAME_MODEL;
    
    VERSION_TRAIN = VERSION_MODEL;
    SITENAME_TRAIN = SITENAME_MODEL;
    
    %% determine date of global starting date and global ending date
    FIRSTDAY = datenum([FIRSTYEAR,1,1]);
    LASTDAY = datenum([ENDYEAR,12,31]);
    
elseif(~isempty(strfind(OPTION,'Prediction')))
    
    % for testing only
%     VERSION_MODEL = 'YanWangTemperature';%main site we are working on
%     SITENAME_MODEL = 'YanWangTemperature';
%     
%     VERSION_DATA = 'YanWangTemperature';%main site the data files have on
%     SITENAME_DATA = 'YanWangTemperature';
    REGIONLIST = SystemPara.REGIONLIST;
    
    % for real prediction only
    VERSION_MODEL = REGIONLIST{Region_ID};%main site we are working on
    SITENAME_MODEL = REGIONLIST{Region_ID};
    
    VERSION_DATA = REGIONLIST{Region_ID};%main site the data files have on
    SITENAME_DATA = REGIONLIST{Region_ID};
    
    VERSION_PREDICT = VERSION_MODEL;% at what grid cell we are making prediction; for all/cv, this is equvalient to VERSION_MODEL
    SITENAME_PREDICT = SITENAME_MODEL;
    
    % for both real prediction and testing
%     if(strcmp(NAME,'MeanTemperature')||strcmp(NAME,'MaxTemperature')||strcmp(NAME,'MinTemperature'))
%         VERSION_TRAIN = 'YanWangTemperature';
%         SITENAME_TRAIN = 'YanWangTemperature';
%     end
    
    if(strcmp(NAME,'MeanTemperature'))
        VERSION_TRAIN = 'YanWangTemperature';%main site we are working on
        SITENAME_TRAIN = 'YanWangTemperature';
    elseif(strcmp(NAME,'MinTemperature'))
        VERSION_TRAIN = 'YanWangTemperature';%main site we are working on
        SITENAME_TRAIN = 'YanWangTemperature';
    elseif(strcmp(NAME,'MaxTemperature'))
        VERSION_TRAIN = 'YanWangTemperature';%main site we are working on
        SITENAME_TRAIN = 'YanWangTemperature';
    elseif(strcmp(NAME,'PM25'))
        VERSION_TRAIN = 'AQRVPM25';
        SITENAME_TRAIN = 'AQRVPM25';
    elseif(strcmp(NAME,'NO2'))
        VERSION_TRAIN = 'EPANO2';%main site we are working on
        SITENAME_TRAIN = 'EPANO2';
    elseif(strcmp(NAME,'Ozone'))
        VERSION_TRAIN = 'EPACastNetOzone';%main site we are working on
        SITENAME_TRAIN = 'EPACastNetOzone';
    end
    
    % determine date of global starting date and global ending date
    FIRSTDAY = datenum([TIME_ID,1,1]);
    LASTDAY = datenum([TIME_ID,12,31]);
end

%% specify the neural network structure
%%% structure of neural network
N_Layer = cell(1,3);
N_Layer{1,1} = [20,20];
N_Layer{1,2} = {'tansig','tansig','purelin'};

%%% more path
DIRPATH_MODEL = [DIRPATH_ROOT,VERSION_MODEL,Sep];%% we are working on this directory
DIRPATH_DATA = [DIRPATH_ROOT,VERSION_DATA,Sep]; %% we go to this directory to find out 1*1 km grid cell data for the whole US
DIRPATH_PREDICT = [DIRPATH_ROOT,VERSION_PREDICT,Sep];%% we make prediction based on data in this directory
DIRPATH_TRAIN = [DIRPATH_ROOT,VERSION_TRAIN,Sep];
% output path
OUTPUTPATH = ['../../assembled_data/',Sep,num2str(year(FIRSTDAY)),'_',num2str(year(LASTDAY)),'_',VERSION_MODEL,'_', ImputationOption,'_',OPTION,'_',NAME,'_',num2str(IDNUM),Sep];%% we store data at this directory
mkdir(OUTPUTPATH);
mkdir([OUTPUTPATH,'BackupCodes',Sep]);
mkdir([OUTPUTPATH,'Pictures',Sep]);
mkdir([DIRPATH_ROOT,'DataProcessRecorder',Sep]);
mkdir([DIRPATH_ROOT,'DataProcessRecorder',Sep,num2str(year(FIRSTDAY)),'_',num2str(year(LASTDAY)),Sep]);

%% export file
FID = fopen([OUTPUTPATH,'Summary_',VERSION_MODEL,'_',NAME,'_',num2str(IDNUM),'.txt'],'a');
FID_check = fopen([OUTPUTPATH,'VariableCheck_',VERSION_MODEL,'_',NAME,'_',num2str(IDNUM),'.txt'],'a');
FID_summary = fopen([DIRPATH_ROOT,'DataProcessRecorder',Sep,num2str(year(FIRSTDAY)),'_',num2str(year(LASTDAY)),Sep,'Summary_',VERSION_MODEL,'_',NAME,'_',num2str(IDNUM),'.txt'],'a');
FID_error_files = fopen([DIRPATH_ROOT,'DataProcessRecorder',Sep,num2str(year(FIRSTDAY)),'_',num2str(year(LASTDAY)),Sep,'ErrorFilesSummary_',VERSION_MODEL,'_',NAME,'_',num2str(IDNUM),'.txt'],'a');
diary([OUTPUTPATH,'DiaryOutput_',VERSION_MODEL,'_',NAME,'_',num2str(IDNUM),'.txt']);


%% specify variables to be used
CSVFile = ['./VariableList_',num2str(IDNUM),'.csv'];
try
    copyfile(['./VariableList_',num2str(IDNUM),'.csv'],OUTPUTPATH);
    
    fprintf('reading %s\n',CSVFile);
    fid_variable = fopen(CSVFile,'r');
    textscan(fid_variable,'%s %s %s %s %s %s %s %s',1);
    C = textscan(fid_variable,'%s %s %s %s %s %s %f %f','Delimiter',',');
    VariableListFull = C{1,1};
    VariableRead = C{1,2};
    VariableUse = C{1,3};
    VariableFolder = C{1,4};
    VariableNameTemplate = C{1,5};
    VariableFormat = C{1,6};
    VariableStartYear = C{1,7};
    VariableEndYear = C{1,8};
    
    % skip line starting with #
    Index = ~arrayfun(@(i) (VariableListFull{i}(1)=='#'),1:1:length(VariableListFull));
    VariableListFull = VariableListFull(Index,:);
    VariableRead = VariableRead(Index,:);
    VariableUse = VariableUse(Index,:);
    VariableFolder = VariableFolder(Index,:);
    VariableNameTemplate = VariableNameTemplate(Index,:);
    VariableFormat = VariableFormat(Index,:);
    VariableStartYear = VariableStartYear(Index,:);
    VariableEndYear = VariableEndYear(Index,:);
    
    VariableList = VariableListFull(ismember(VariableUse,'T'));
    fclose(fid_variable);
catch exception
    fprintf('reading csv error...%s\n',exception.message);
    disp('error! NO VariableList.csv found!');
    
    fclose(FID);
    fclose(FID_check);
    fclose(FID_summary);
    fclose(FID_error_files);
    return;
end
fprintf('variables to be read...\n');
fprintf('%s\n',strjoin(VariableList,'\n'));
fprintf(FID,'variables to be read...\n');
fprintf(FID,'%s\n',strjoin(VariableList,'\n'));

% back up some files...
list = dir(CODEPATH);
for i=3:length(list)
    if(isempty(strfind(list(i).name,'.m')))
        continue;
    end
    try
        copyfile(list(i).name,[OUTPUTPATH,'BackupCodes',Sep]);
    catch exception
        fprintf('%s\n',exception.message);
    end
end

%% Load site data
SiteData = LoadData_function([DIRPATH_MODEL,'Location',Sep,SITENAME_MODEL,'Site_',GCS,'.mat'],'SiteData');% the monitor sites we will use for modeling;monitors of interest
SiteData_Model = SiteData.SiteData;
N_ModelSite = size(SiteData_Model,1);

%% Load shapefile Data
ShapeFile1 = ['../data/shapefile',Sep,'US_WGS_clipcoast_2.shp'];

%% wait time
%%sometime one prediction thread needs to wait for other thread to finish.
%%The waiting time is specified by seconds;
if(N_ModelSite<10000)
    WAITTIME = 10;
else
    WAITTIME = 100;
end

%% define global parameters
EnvironPara = struct(...
    'Sep',Sep,...
    'OPTION',OPTION,...
    ...
    'SOLUTION',ImputationOption,...
    'ImputationOption_Model',ImputationOption_Model,...
    'NAME',NAME,...
    ...
    'DIRPATH_ROOT',DIRPATH_ROOT,...
    'SITENAME_MODEL',SITENAME_MODEL,...%the site name corresponds to DIRPATH_MODEL
    'DIRPATH_MODEL',DIRPATH_MODEL,...
    'SITENAME_PREDICT',SITENAME_PREDICT,...
    'DIRPATH_PREDICT',DIRPATH_PREDICT,...
    'SITENAME_TRAIN',SITENAME_TRAIN,...
    'DIRPATH_TRAIN',DIRPATH_TRAIN,...
    'SITENAME_DATA',SITENAME_DATA,...
    'DIRPATH_DATA',DIRPATH_DATA,...
    'OUTPUTPATH',OUTPUTPATH,...
    'OUTPUTPATH_PIC',[OUTPUTPATH,'Pictures',Sep],...
    ...
    'IDNUM',IDNUM,...
    ...
    'FID',FID,...
    'FID_SUMMARY',[],...
    'FID_CHECK',FID_check,...
    'FID_error_files',FID_error_files,...
    ...
    'N_Neighbour',N_Neighbour,...
    'IsRecord',IsRecord,...
    'IsSaveCheckingFigure',IsSaveCheckingFigure,...
    'IsTest',TESTING,...
    'N_Layer',{N_Layer},...
    'ModelName',ModelName,...
    ...
    'SiteData_Model',{SiteData_Model},...
    'N_Site',N_ModelSite,...
    'GCS',GCS,...
    ...
    'TIME_ID',TIME_ID,...
    'FIRSTDAY',FIRSTDAY,...
    'LASTDAY',LASTDAY,...
    'N_Day',{-1},...
    'StartDate',{-1},...
    'EndDate',{-1},...
    'StartYear',FIRSTYEAR,...
    'EndYear',ENDYEAR,...
    ...
    'VariableListFull',{VariableListFull},...
    'VariableList',{VariableList},...
    'VariableRead',{VariableRead},...
    'VariableUse',{VariableUse},...
    'VariableFolder',{VariableFolder},...
    'VariableNameTemplate',{VariableNameTemplate},...
    'VariableFormat',{VariableFormat},...
    'VariableStartYear',{VariableStartYear},...
    'VariableEndYear',{VariableEndYear},...
    ...
    'N_Var',length(VariableList),...
    ...
    'WAIT_TIME',WAITTIME,...
    ...
    'S_back',geoshape(shaperead('../../raw_data/data/shapefile/US_WGS_clipcoast_2.shp', 'UseGeoCoords', true)),...
    'states',geoshape(shaperead('usastatehi', 'UseGeoCoords', true)),...
    'US_North_America_Equidistant_Conic', '../../raw_data/data/shapefile/US_North_America_Equidistant_Conic.shp'...
    );

EnvironPara.SITENAME_SUBSET = EnvironPara.SITENAME_MODEL;

%% model training
if(~isempty(strfind(OPTION,'CV'))||~isempty(strfind(OPTION,'All')))
    EnvironPara.StartDate = FIRSTDAY;
    EnvironPara.EndDate = LASTDAY;
    EnvironPara.N_Day = datenum(EnvironPara.EndYear,12,31) - datenum(EnvironPara.StartYear,1,1) + 1;
    EnvironPara.CommonStartPoint = datestr(datenum(EnvironPara.StartYear,1,1),'yyyymmdd');%%start date
    EnvironPara.CommonEndPoint = datestr(datenum(EnvironPara.EndYear,12,31),'yyyymmdd');%%
    
    
    
    OutputFileName = ['Output_',EnvironPara.NAME,'_',EnvironPara.SITENAME_MODEL,'_',EnvironPara.CommonStartPoint,'_',EnvironPara.CommonEndPoint];
    
    IsProcessFlag = false;
    try
        MonitorPredicted = LoadData_function([EnvironPara.OUTPUTPATH,OutputFileName,'.mat'],'MonitorPredicted');
        MonitorPredicted = MonitorPredicted.MonitorPredicted;
        IsProcessFlag = false;
    catch
        if(exist([EnvironPara.OUTPUTPATH,OutputFileName,'.mat.part'],'file'))
            list = dir([EnvironPara.OUTPUTPATH,OutputFileName,'.mat.part']);
            if(etime(clock,datevec(list.datenum))>60*120)
                delete([EnvironPara.OUTPUTPATH,OutputFileName,'.mat.part']);
                fprintf('%s process work started long time ago and still has not completed. we do it again!...\n',OutputFileName);
                IsProcessFlag = true;
            else
                fprintf('%s is in the middle processing!...\n',OutputFileName);
                IsProcessFlag = false;
            end
        else
            IsProcessFlag = true;
        end
    end
    
    if(~IsProcessFlag)
        fprintf('%s has been processed!\n',OutputFileName);
    else
         %save([EnvironPara.OUTPUTPATH,OutputFileName,'.mat.part'],'Sep');
        
        %% intially create lagged terms
        if(strcmp(OPTION,'CV_round1')||strcmp(OPTION,'All_round1'))
        
        else
            Analysis_TwoStep([],EnvironPara.SiteData_Model,EnvironPara.SiteData_Model,EnvironPara.SITENAME_MODEL,EnvironPara.SITENAME_MODEL,EnvironPara,'CreateWeight');
        end 
        
        %% reading cross-valiation index
        if(~isempty(strfind(OPTION,'CV')))
            FileName = [DIRPATH_MODEL,Sep,'Temp',Sep,'CrossValidationIndex',datestr(datenum(EnvironPara.StartYear,1,1),'yyyymmdd'),'_',datestr(datenum(EnvironPara.EndYear,12,31),'yyyymmdd'),'.mat'];
            try
%                 need to change!!!
                IndexCV = LoadData_function(FileName);
                IndexCV = IndexCV.IndexCV;
            catch
                IndexCV = Analysis_CreateTenfold(EnvironPara.N_Site,EnvironPara.N_Day);
                save(FileName,'IndexCV');
            end
        else
            IndexCV = [];
        end
        
        %% reading input data and doing imputation
        if(exist([OUTPUTPATH,'InputData_',NAME,'_',SITENAME_MODEL,'_',EnvironPara.CommonStartPoint,'_',EnvironPara.CommonEndPoint,'.mat'],'file'))
            load([OUTPUTPATH,'InputData_',NAME,'_',SITENAME_MODEL,'_',EnvironPara.CommonStartPoint,'_',EnvironPara.CommonEndPoint,'.mat']);
        else
            Notice = [repmat('*',[1,50]),'\nread data set for temperature modeling\n',repmat('*',[1,50]),'\n'];
            fprintf(Notice);
            [Input_var_AOD,VariableList_AOD] = Analysis_ReadInputData(EnvironPara,'first');
            % reshape and descriptive results
            fprintf([repmat('*',[1,50]),'\n']);
            Analysis_InputDataDescriptiveAnalysis(Input_var_AOD,EnvironPara,true,'total dataset',{'WHOLE',false,false});
            fprintf([repmat('*',[1,50]),'\n']);
            Analysis_InputDataDescriptiveAnalysis(Input_var_AOD(~isnan(Input_var_AOD(:,1)),:),EnvironPara,true,'for only monitoring data',{'WHOLE',false,false});
            % imputation
            fprintf([repmat('*',[1,50]),'\n']);
            Input_var_AOD = Analysis_ImputationInputData(Input_var_AOD,EnvironPara,ImputationOption);
            Analysis_InputDataDescriptiveAnalysis(Input_var_AOD,EnvironPara,true,'after imputation',{'WHOLE',false,false});
            %save
            save([OUTPUTPATH,'InputData_',NAME,'_',SITENAME_MODEL,'_',EnvironPara.CommonStartPoint,'_',EnvironPara.CommonEndPoint,'.mat'],'Input_var_AOD','VariableList_AOD');
        end
        
        %% fit main model
        EnvironPara.FID_SUMMARY = FID_summary;
        MonitorData = reshape(Input_var_AOD(:,1),[EnvironPara.N_Day,EnvironPara.N_Site]);
        Notice = [repmat('*',[1,50]),'\nrunning Nain Neural Network models...\n',repmat('*',[1,50]),'\n'];
        fprintf(Notice);
        MonitorPredicted = nan(EnvironPara.N_Day,EnvironPara.N_Site);
        NetPredicted = {};
        % fit mode, save output and visualization
        Input_var_AOD(:,1) = Analysis_Transform_fucntion1(Input_var_AOD(:,1),NAME,'transform');
        [Input_var_AOD,MonitorPredicted,NetPredicted] = Analysis_FitModel(Input_var_AOD,MonitorPredicted,NetPredicted,0,IndexCV,EnvironPara.N_Site,EnvironPara.N_Day,EnvironPara);
        MonitorPredicted = Analysis_Transform_fucntion1(MonitorPredicted,NAME,'inverse');
        %save([OUTPUTPATH,'InputDataPostAOD_',NAME,'_',SITENAME_MODEL,'_',EnvironPara.CommonStartPoint,'_',EnvironPara.CommonEndPoint,'.mat'],'Input_var_AOD','VariableList_AOD');
        save([EnvironPara.OUTPUTPATH,OutputFileName,'.mat'],'MonitorPredicted','MonitorData','NetPredicted','SiteData');
        [R,~,~] = Analysis_ResultDescriptiveAnalysis(MonitorPredicted,MonitorData,NetPredicted,EnvironPara,OutputFileName,{true,true,false,true});    
    end
        
elseif(strcmp(OPTION,'Prediction'))
    
    %% Load stage 1 outputs
    SiteData = LoadData_function([EnvironPara.DIRPATH_TRAIN,'Temp',Sep,num2str(IDNUM),Sep,'Output_',EnvironPara.NAME,'_',EnvironPara.SITENAME_TRAIN,'_',datestr(datenum(TIME_ID,1,1),'yyyymmdd'),'_',datestr(datenum(TIME_ID,12,31),'yyyymmdd'),'.mat']);
    MonitorPredicted_ModelFit = SiteData.MonitorPredicted;
    MonitorData = SiteData.MonitorData;
    NetPredicted = SiteData.NetPredicted;
        
    %% Load Site Data for training data set
    SiteData = LoadData_function([EnvironPara.DIRPATH_TRAIN,'Location',Sep,EnvironPara.SITENAME_TRAIN,'Site_',GCS,'.mat'],'SiteData');
    SiteData_Train = SiteData.SiteData;
    EnvironPara.N_Day = 1;
    
    
    %% intially create lagged terms
    if(length(NetPredicted)>1)
        Analysis_TwoStep([],EnvironPara.SiteData_Model,SiteData_Train,EnvironPara.SITENAME_MODEL,EnvironPara.SITENAME_TRAIN,EnvironPara,'CreateWeight');
    end
    %% use for result check...
    MonitorPredicted_interp = nan([yeardays(TIME_ID),size(SiteData_Train,1)]);
    %% predict recursively -- PM2.5 used to have 2 rounds of modeling, for temperature modeling, there is only one round of modeling
    for l=1:length(NetPredicted)
        if(1 == l)
            %% load input data set
            Input_var_AOD = LoadData_function([EnvironPara.DIRPATH_TRAIN,'Temp',Sep,num2str(IDNUM),Sep,'InputData_',EnvironPara.NAME,'_',EnvironPara.SITENAME_TRAIN,'_',datestr(datenum(TIME_ID,1,1),'yyyymmdd'),'_',datestr(datenum(TIME_ID,12,31),'yyyymmdd'),'.mat'],'Input_var_AOD');
            Input_var_AOD_Train = Input_var_AOD.Input_var_AOD;
            Input_var_AOD_Train = reshape(Input_var_AOD_Train,[datenum(TIME_ID,12,31) - datenum(TIME_ID,1,1) + 1,size(Input_var_AOD_Train,1)/(datenum(TIME_ID,12,31) - datenum(TIME_ID,1,1) + 1),size(Input_var_AOD_Train,2)]);
        elseif(2 == l)
            %% load input data set
            Input_var_AOD = LoadData_function([EnvironPara.DIRPATH_TRAIN,'Temp',Sep,num2str(IDNUM),Sep,'InputDataPostAOD_',EnvironPara.NAME,'_',EnvironPara.SITENAME_TRAIN,'_',datestr(datenum(TIME_ID,1,1),'yyyymmdd'),'_',datestr(datenum(TIME_ID,12,31),'yyyymmdd'),'.mat'],'Input_var_AOD');
            Input_var_AOD_Train = Input_var_AOD.Input_var_AOD;
            Input_var_AOD_Train = reshape(Input_var_AOD_Train,[datenum(TIME_ID,12,31) - datenum(TIME_ID,1,1) + 1,size(Input_var_AOD_Train,1)/(datenum(TIME_ID,12,31) - datenum(TIME_ID,1,1) + 1),size(Input_var_AOD_Train,2)]);

        end
        
        
        %% Load stage 1 result
        NetPredicted_temp = NetPredicted{l};
        MonitorPredicted_temp = NetPredicted_temp{6};  
        if(length(NetPredicted)==l)
            MonitorPredicted_temp = Analysis_Transform_fucntion1(MonitorPredicted_temp,NAME,'inverse');
        end
        fprintf('fitting model %d\t...%s model\n',l,NetPredicted_temp{5});
        
        %% all have been predicted
        TempFileName = arrayfun(@(i)  ([OUTPUTPATH, 'PREDICTION',NetPredicted_temp{5},num2str(l),'_',NAME,'_',EnvironPara.SITENAME_MODEL,'_',datestr(i,'yyyymmdd'),'_',datestr(i,'yyyymmdd'),'.mat']),FIRSTDAY:LASTDAY,'UniformOutput',false); 
        if(cellfun(@(x) exist(x,'file'),TempFileName))%% all have been processed
            fprintf('SKIP prediction! ...%s..%s\n',TempFileName{1},TempFileName{length(TempFileName)});
        else
            %% predict on every day
            i = FIRSTDAY;
            while i>=FIRSTDAY && i<=LASTDAY
                EnvironPara.StartDate = i;
                EnvironPara.EndDate = i;
                EnvironPara.CommonStartPoint = datestr(i,'yyyymmdd');%%start date
                EnvironPara.CommonEndPoint = datestr(i,'yyyymmdd');%%

                %% if pre-processed skip
                %%file name example:PREDICTIONFinal1_PM25_USPM25_20080101_20080101.mat           
                FileName = ['PREDICTION',NetPredicted_temp{5},num2str(l),'_',NAME,'_',EnvironPara.SITENAME_MODEL,'_',EnvironPara.CommonStartPoint,'_',EnvironPara.CommonEndPoint];

                fprintf('%s\n',repmat('.',[1,50]));
                fprintf('predicting...%s\n',FileName);

                IsProcessFlag = false;
                try
                    MonitorPredicted = LoadData_function([OUTPUTPATH,FileName,'.mat'],'MonitorPredicted');
                    MonitorPredicted = MonitorPredicted.MonitorPredicted;
                    IsProcessFlag = false;
                catch
                    if(exist([OUTPUTPATH,FileName,'.mat.part'],'file'))
                        list = dir([OUTPUTPATH,FileName,'.mat.part']);
                        if(etime(clock,datevec(list.datenum))>60)
                            delete([OUTPUTPATH,FileName,'.mat.part']);
                            fprintf('%s process work started long time ago and still has not completed. we do it again!...\n',FileName);
                            IsProcessFlag = true;
                        else
                            fprintf('%s is in the middle processing!...\n',FileName);
                            IsProcessFlag = false;
                        end
                    else
                        IsProcessFlag = true;
                    end
                end

                if(~IsProcessFlag)
                    fprintf('%s...%s has been processed!\n',EnvironPara.CommonStartPoint,EnvironPara.CommonEndPoint);
                else
                    save([OUTPUTPATH,FileName,'.mat.part'],'i');
                    %% begin to predict!!!
                    fprintf('%s...%s is being processed!\n',EnvironPara.CommonStartPoint,EnvironPara.CommonEndPoint);

                    %% reading data
                    Notice = [repmat('*',[1,50]),'\nread data set for temperature modeling\n',repmat('*',[1,50]),'\n'];
                    fprintf(Notice);

                    if(1 == l)
                        [Input_var_AOD,VariableList_AOD] = Analysis_ReadInputData(EnvironPara,'first');
                    elseif(2 == l)
                        [Input_var_AOD,VariableList_AOD] = Analysis_ReadInputData(EnvironPara,'second');
                    end

                    %% check data
                    if(true)
                        if(1==l)
                            if(i==FIRSTDAY)
                                Analysis_InputDataDescriptiveAnalysis(Input_var_AOD,EnvironPara,true,'total dataset',{'WHOLE',true&IsSaveCheckingFigure,true&IsSaveCheckingFigure});
                            else
                                Analysis_InputDataDescriptiveAnalysis(Input_var_AOD,EnvironPara,true,'total dataset',{'WHOLE',true&IsSaveCheckingFigure,false});
                            end
                            Analysis_InputDataDescriptiveAnalysis(Input_var_AOD,EnvironPara,true,'after transform',{'WHOLE',false,false});
                        else
                            Analysis_InputDataDescriptiveAnalysis(Input_var_AOD,EnvironPara,true,'total dataset',{'SUBSET1',true&IsSaveCheckingFigure,false});
                            Analysis_InputDataDescriptiveAnalysis(Input_var_AOD,EnvironPara,true,'after transform',{'SUBSET1',false,false});
                        end
                    end   

                    %% imputation
                    fprintf([repmat('*',[1,50]),'\n']);
                    Input_var_AOD = Analysis_ImputationInputData(Input_var_AOD,EnvironPara,ImputationOption);


                    %% check with input data from stage1 model
                    if(false)
                        if(true)
                            %                     if(1==l)
                            for m = 1:EnvironPara.N_Var
    %                             %%we do not check lagged term the first time
    %                             if(~isempty(strfind(EnvironPara.VariableList{m},'Lagged')) && 1 ==l)
    %                                 continue;
    %                             end
    %                             
    %                             %%we noly check lagged term after the first term
    %                             if(isempty(strfind(EnvironPara.VariableList{m},'Lagged')) && 1<=l)
    %                                 continue;
    %                             end

                                if(size(SiteData_Model,1) == size(SiteData_Train,1))
                                    Analysis_ResultDescriptiveAnalysis(Input_var_AOD_Train(i-FIRSTDAY+1,:,m)',Input_var_AOD(:,m),[],EnvironPara,[EnvironPara.VariableList{m},' ',EnvironPara.CommonStartPoint],{true,false,false,false});
                                    fprintf('compare:%d %d\n',min(Input_var_AOD_Train(i-FIRSTDAY+1,:,m)'-Input_var_AOD(:,m)),max(Input_var_AOD_Train(i-FIRSTDAY+1,:,m)'-Input_var_AOD(:,m)));
                                else
                                    if(iscell(SiteData_Model))
                                        Input_var_AOD_interp = InterpMyData_2(Input_var_AOD(:,m)',cell2mat(SiteData_Model(:,3)),cell2mat(SiteData_Model(:,2)),cell2mat(SiteData_Train(:,3)),cell2mat(SiteData_Train(:,2)),'');
                                    elseif(istable(SiteData_Model))
                                        Input_var_AOD_interp = InterpMyData_2(Input_var_AOD(:,m)',SiteData_Model.Lon,SiteData_Model.Lat,SiteData_Train.Lon,SiteData_Train.Lat,'');
                                    else
                                        Input_var_AOD_interp = InterpMyData_2(Input_var_AOD(:,m)',SiteData_Model(:,3),SiteData_Model(:,2),cell2mat(SiteData_Train(:,3)),cell2mat(SiteData_Train(:,2)),'');
                                    end
                                    Analysis_ResultDescriptiveAnalysis(Input_var_AOD_Train(i-FIRSTDAY+1,:,m)',Input_var_AOD_interp',[],EnvironPara,[EnvironPara.VariableList{m},' ',EnvironPara.CommonStartPoint],{true,false,false,false});
                                    fprintf('compare:%d %d\n',min(Input_var_AOD_Train(i-FIRSTDAY+1,:,m)'-Input_var_AOD_interp'),max(Input_var_AOD_Train(i-FIRSTDAY+1,:,m)'-Input_var_AOD_interp'));
                                end
                            end
                        end
                    end


                    %% Load prior prediction results
                    if(1==l)
                        MonitorPredicted = nan(EnvironPara.N_Day,EnvironPara.N_Site);
                    else
                        MonitorPredicted = LoadData_function([OUTPUTPATH,strrep(EnvironPara.PredictionResultName,'@@@@',EnvironPara.CommonStartPoint),'.mat'],'MonitorPredicted');
                        MonitorPredicted = MonitorPredicted.MonitorPredicted;
                    end

                    %% predict!!!
                    %%pay attention to whether MonitorPredicted is at log scale or not
                    Input_var_AOD(:,1) = Analysis_Transform_fucntion1(Input_var_AOD(:,1),NAME,'transform');
                    [~,MonitorPredicted,~] = Analysis_FitModel(Input_var_AOD,MonitorPredicted,NetPredicted_temp,l,[],EnvironPara.N_Site,EnvironPara.N_Day,EnvironPara);
                    %%the last time is at normal scale
                    if(length(NetPredicted)==l)
                        MonitorPredicted = Analysis_Transform_fucntion1(MonitorPredicted,NAME,'inverse');
                    end
                    %% save result
                    save([OUTPUTPATH,FileName,'.mat'],'MonitorPredicted');
                    delete([OUTPUTPATH,FileName,'.mat.part']);

                    %% plot result
                     if(false)
                         if(EnvironPara.N_Site>10000)
                             try
                                  Visualization_USResult_1(FileName,MonitorPredicted,SiteData_Model,EnvironPara);
                             catch exception
                                 fprintf('%s\n',exception.message);
                             end
                         end
                     end

                    %% check results and save R2
                    if(iscell(SiteData_Model))
                        MonitorPredicted_1 = InterpMyData_2(MonitorPredicted,cell2mat(SiteData_Model(:,3)),cell2mat(SiteData_Model(:,2)),cell2mat(SiteData_Train(:,3)),cell2mat(SiteData_Train(:,2)),'');
                    elseif(istable(SiteData_Model))
                        MonitorPredicted_1 = InterpMyData_2(MonitorPredicted,SiteData_Model.Lon,SiteData_Model.Lat,SiteData_Train.Lon,SiteData_Train.Lat,'');
                    else
                        MonitorPredicted_1 = InterpMyData_2(MonitorPredicted,SiteData_Model(:,3),SiteData_Model(:,2),cell2mat(SiteData_Train(:,3)),cell2mat(SiteData_Train(:,2)),'');
                    end
                    Analysis_ResultDescriptiveAnalysis(MonitorPredicted_1',MonitorPredicted_temp(i-datenum(EnvironPara.StartYear,1,1)+1,:)',[],EnvironPara,[num2str(l),'Comparison with training result!MonitorPredicted ',EnvironPara.CommonStartPoint],{true,false,false,true});
                    fprintf('compare:%d %d\n',min(MonitorPredicted_1'-MonitorPredicted_temp(i-datenum(EnvironPara.StartYear,1,1)+1,:)'),max(MonitorPredicted_1'-MonitorPredicted_temp(i-datenum(EnvironPara.StartYear,1,1)+1,:)'));

                    MonitorPredicted_interp(i-FIRSTDAY+1,:) = MonitorPredicted_1;

                    %% descriptive stat of result
                    Analysis_ResultDescriptiveAnalysis(MonitorPredicted,[],[],EnvironPara,[],{false,false,false,true});
                end

                %% next day of prediction
                i = i + 1;
            end
        end
        
        %% update
        if(l<length(NetPredicted))
            EnvironPara.PredictionResultName = ['PREDICTION',NetPredicted_temp{5},num2str(l),'_',NAME,'_',EnvironPara.SITENAME_MODEL,'_','@@@@','_','@@@@'];
            Analysis_TwoStep(MonitorPredicted_temp,SiteData_Model,SiteData_Train,EnvironPara.SITENAME_MODEL,EnvironPara.SITENAME_TRAIN,EnvironPara,'UPDATE');        
        end
    end
    
    %% prediction results vs. stage 1 model fiting results
    Analysis_ResultDescriptiveAnalysis(MonitorPredicted_interp,MonitorPredicted_ModelFit,[],EnvironPara,'MonitorPredicted_interp vs. MonitorPredicted',{true,false,false,false});
    Analysis_ResultDescriptiveAnalysis(MonitorPredicted_interp,MonitorData,[],EnvironPara,'MonitorData vs. MonitorPredicted',{true,false,false,false});
    MonitorPredicted = MonitorPredicted_interp;
    save([OUTPUTPATH,'PREDICTIONInterp_',NAME,'_',EnvironPara.SITENAME_PREDICT,'_',num2str(EnvironPara.StartYear),'_',num2str(EnvironPara.EndYear),'.mat'],'MonitorPredicted');
    
    %% delete intermediate lagged terms
    %Analysis_TwoStep([],EnvironPara.SiteData_Model,EnvironPara.SiteData_Model,EnvironPara.SITENAME_MODEL,EnvironPara.SITENAME_MODEL,EnvironPara,'DeleteTerms');
    
end
diary off;
fclose all;
ReturnStatue = 1;

end

%% create simple 10-fold cross-valiation on sites
% this is random split on monitoring site, not on days
% this fold is not supposed to run as a standalone code

%% parameter: 
% N_Site: the number of monitoring site;
% N_Day: the number of days
function Index = Analysis_CreateTenfold(N_Site,N_Day)
% N_Site = 1928;
% N_Day = 366;

N = 10;
Index = zeros(N_Site,N);

k = rand(1,N_Site);
[~,n] = sort(k);%n: random number from 1 to N;
Bins = round(linspace(0,N_Site,N+1));

for i=1:N
    Temp = n(Bins(i)+1:Bins(i+1));
    Index(Temp,i) = 1;%%% 1 means testing data;
end

Index = logical(Index);
Index = repmat(Index,[N_Day,1]);

end


