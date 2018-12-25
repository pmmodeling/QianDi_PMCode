%% this is the main function to do assemble data together

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
% 2018-03-26:create Run_AssembleDataSummary based on Analysis_RunModel;
% keep only the data reading part

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

% examples 2000-2014:
% Run_AssembleDataSummary(1,99941,2008,'All_round3','PM25','AllRecords','');
% Run_AssembleDataSummary(1,99941,2008,'All_round3','Ozone','AllRecords','');
% Run_AssembleDataSummary(1,99941,2008,'All_round3','NO2','AllRecords','');

% example 2015
% Run_AssembleDataSummary(1,999411,2015,'All_round3','PM25','AllRecords','');
% Run_AssembleDataSummary(1,999411,2015,'All_round3','Ozone','AllRecords','');
% Run_AssembleDataSummary(1,999411,2015,'All_round3','NO2','AllRecords','');


%% code
function ReturnStatue = Run_AssembleDataSummary(Region_ID,IDNUM,TIME_ID,OPTION,NAME,ImputationOption,SystemPara)

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
    'TESTING',1,...%%% testing mode to run faster...
    'MULTICORE',1,...%%% use multicore? MULTICORE is the number of cores
    'ModelName','NeuralNetwork',...%%% Machine Learning Method, 'NeuralNetwork' or 'Lasso'
    'FIRSTYEAR',TIME_ID,...
    'ENDYEAR',TIME_ID);
end

%% set the imputation model for different year
if(TIME_ID>2014)
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
ShapeFile1 = ['../../raw_data/data/shapefile',Sep,'US_WGS_clipcoast_2.shp'];

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
    'OPTION','All_round3',...
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
    'S_back',geoshape(shaperead(ShapeFile1, 'UseGeoCoords', true)),...
    'states',geoshape(shaperead('usastatehi', 'UseGeoCoords', true)),...
    'US_North_America_Equidistant_Conic',[DIRPATH_ROOT,'Shapefile',Sep,'US_North_America_Equidistant_Conic.shp']...
    );

EnvironPara.SITENAME_SUBSET = EnvironPara.SITENAME_MODEL;


%% model training

EnvironPara.StartDate = FIRSTDAY;
EnvironPara.EndDate = LASTDAY;
EnvironPara.N_Day = datenum(EnvironPara.EndYear,12,31) - datenum(EnvironPara.StartYear,1,1) + 1;
EnvironPara.CommonStartPoint = datestr(datenum(EnvironPara.StartYear,1,1),'yyyymmdd');%%start date
EnvironPara.CommonEndPoint = datestr(datenum(EnvironPara.EndYear,12,31),'yyyymmdd');%%



OutputFileName = ['Output_',EnvironPara.NAME,'_',EnvironPara.SITENAME_MODEL,'_',EnvironPara.CommonStartPoint,'_',EnvironPara.CommonEndPoint];

%% check whether another matlab session is in processing;
% if yes, skip
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
        %save
        save([OUTPUTPATH,'InputData_',NAME,'_',SITENAME_MODEL,'_',EnvironPara.CommonStartPoint,'_',EnvironPara.CommonEndPoint,'.mat'],'Input_var_AOD','VariableList_AOD','SiteData_Model');
    end
end
        
    

diary off;
fclose all;
ReturnStatue = 1;

end