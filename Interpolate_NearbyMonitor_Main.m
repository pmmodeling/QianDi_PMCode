%% calculate averaged monitoring data

%% version history
% 2016-11- 12: change INPUTYEAR to INPUTDATE
% adapted from Interpolate_OMI_Main

%% parameters
% SITENAME: name of the grid cells to be interpolated at;
% INPUTDATE: time period of interest;
% VARNAME: variable name of monitoring data name;
% OPTION: output format: by day or by year;
% DATA_OPTION: what kind of output? daily mean, daily min, daily max, or lagged value?
% EXTRAOPTION: 1 --> only create weight files; 0--> interpolate the data; 2 --> delete all files; 3: create maps while doing interpolation
% FORMAT: 'mat': output file is matlab file; 'h5': output file is hdf file

%% return values:
% ReturnStatus==0: still in processing
% ReturnStatus==1: processing done before; 
% ReturnStatus==2, some raw data not available and needs to do it again;
% ReturnStatus==3, some input data not available has to quit;

%% example:
% Interpolate_NearbyMonitor_Main('AQRVPM25',[datenum(2012,1,1),datenum(2012,12,31)],'PM25','By-Year','Peak2_Lag1',0)

%% code
function ReturnStatus = Interpolate_NearbyMonitor_Main(SITENAME,INPUTDATE,VARNAME,OPTION,DATA_OPTION_All,EXTRAOPTION,FORMAT)


%% marco-definition
% THRESHOLD = 999*1000;
Temp = strsplit(DATA_OPTION_All,'_');
DATA_OPTION = Temp{1};
DATA_LAG = str2num(strrep(Temp{2},'Lag',''));

% output format, output file saved by day or year?
INPUTYEAR = year(INPUTDATE(1));
if(strcmp(OPTION,'By-Date'))
    DATE_TO_PROCESS = (INPUTDATE(1):1:INPUTDATE(2)) - datenum(INPUTYEAR,1,1)+1;%index: day of year
elseif(strcmp(OPTION,'By-Day') || strcmp(OPTION,'By-Year'))
    DATE_TO_PROCESS = 1:1:yeardays(INPUTYEAR);
end

% the distance threshold for each variable --- avoid take too many nearby monitoring site into calcualtion, reducing computation load
if(any(strfind(VARNAME,'MaxTemperature')))
    SITENAME_DATA = 'YanWangTemperature';
    THRESHOLD = Inf;
    N_Neighbour = 2000;
elseif(any(strfind(VARNAME,'MinTemperature')))
    SITENAME_DATA = 'YanWangTemperature';
    THRESHOLD = Inf;
    N_Neighbour = 2000;
elseif(any(strfind(VARNAME,'MeanTemperature')))
    SITENAME_DATA = 'YanWangTemperature';
    THRESHOLD = Inf;
    N_Neighbour = 2000;
elseif(any(strfind(VARNAME,'PM25')))
    SITENAME_DATA = 'AQRVPM25';
    THRESHOLD = Inf;
    N_Neighbour = 3000;
elseif(any(strfind(VARNAME,'NO2')))
    SITENAME_DATA = 'EPANO2';
    THRESHOLD = Inf;
    N_Neighbour = 3000;
elseif(any(strfind(VARNAME,'Ozone')))
    SITENAME_DATA = 'EPACastNetOzone';
    THRESHOLD = Inf;
    N_Neighbour = 3000;
end

%% path
DIRPATH = ['../../processed_data/',SITENAME_DATA,'/Monitor/'];
OUTPUTPATH_ROOT = '../../processed_data/';
Sep = '/';

% file path to the shape file used for making maps

EnvironPara = struct();
EnvironPara.states = geoshape(shaperead('usastatehi', 'UseGeoCoords', true));
EnvironPara.S_back = geoshape(shaperead('../data/shapefile/US_WGS_clipcoast_2.shp', 'UseGeoCoords', true));
EnvironPara.US_North_America_Equidistant_Conic = '../data/shapefile/US_North_America_Equidistant_Conic.shp';
EnvironPara.US = '../data/shapefile/US.shp';
GCS = 'North_America_Equidistant_Conic';%% projection system to be used

FileNameTemplate = 'MONITOR_$VARNAME$_$SITENAME_DATA$_$STARTDATE$_$ENDDATE$.mat';
OutputTemplate = ['MONITOR_$OPTION$_$VARNAME$_$SITENAME$_$STARTDATE$_$ENDDATE$.',FORMAT];

FileNameTemplate = strrep(strrep(FileNameTemplate,'$SITENAME_DATA$',SITENAME_DATA),'$VARNAME$',VARNAME);
OutputTemplate = strrep(strrep(strrep(OutputTemplate,'$OPTION$',[DATA_OPTION,'Lag',num2str(DATA_LAG),'_Thres',num2str(THRESHOLD)]),'$SITENAME$',SITENAME),'$VARNAME$',VARNAME);

%% skip if output files exist
OUTPUTPATH = [OUTPUTPATH_ROOT,SITENAME,Sep,'Nearby',Sep];
EnvironPara.OUTPUTPATH_PIC = OUTPUTPATH;
mkdir(OUTPUTPATH);
% check whether has been processed before
if(strcmp(OPTION,'By-Year'))
    TempFileName = [OUTPUTPATH, strrep(strrep(OutputTemplate,'$STARTDATE$',datestr(datenum(INPUTYEAR,1,1),'yyyymmdd')),'$ENDDATE$',datestr(datenum(INPUTYEAR,12,31),'yyyymmdd'))];
    if(2==EXTRAOPTION)
        delete(TempFileName);
        fprintf('delete ...%s\n',TempFileName);
        ReturnStatus = 1;
        return;
    else
        if(exist(TempFileName, 'file'))
            fprintf('SKIP! ...%s\n',TempFileName);
            ReturnStatus = 1;
            return;
        end
    end
elseif(strcmp(OPTION,'By-Day')|| strcmp(OPTION,'By-Date'))
    StartDate = datenum(INPUTYEAR,1,1);
    TempFileName = arrayfun(@(i)  ([OUTPUTPATH, strrep(strrep(OutputTemplate,'$STARTDATE$',datestr(StartDate+i-1,'yyyymmdd')),'$ENDDATE$',datestr(StartDate+i-1,'yyyymmdd'))]),DATE_TO_PROCESS,'UniformOutput',false);
    if(2==EXTRAOPTION)
        for i=1:length(TempFileName)
            delete(TempFileName{i});
            fprintf('delete ...%s\n',TempFileName{i});
        end
        ReturnStatus = 1;
        return;
    else
        if(cellfun(@(x) exist(x,'file'),TempFileName))%% all have been processed
            fprintf('SKIP! ...%s..%s\n',TempFileName{1},TempFileName{length(TempFileName)});
            ReturnStatus = 1;
            return;
        end
    end
end

%% read points of interests
% these are the grid cells we will interpolate at
SiteData = LoadData_function([OUTPUTPATH_ROOT,SITENAME,Sep,'Location',Sep,SITENAME,'Site_',GCS,'.mat'],'SiteData');
SiteData_Target = SiteData.SiteData;
if(iscell(SiteData_Target))
    Lon_Target = cell2mat(SiteData_Target(:,3));
    Lat_Target = cell2mat(SiteData_Target(:,2));
elseif(istable(SiteData_Target))
    Lon_Target = SiteData_Target.Lon;
    Lat_Target = SiteData_Target.Lat;
else
    Lon_Target = SiteData_Target(:,3);
    Lat_Target = SiteData_Target(:,2);
end
N_Traget = length(Lon_Target);

%% calculate weight --- used for interpolation (weighted average of nearby grid cells)
mkdir([OUTPUTPATH_ROOT,SITENAME,Sep,'Temp',Sep]);
WeightFile = [OUTPUTPATH_ROOT,SITENAME,Sep,'Temp',Sep,'NeighbourWeightMatrix_',[DATA_OPTION,'_Thres',num2str(THRESHOLD)],'_',SITENAME_DATA,'_',SITENAME,'.mat'];
if(exist(WeightFile,'file'))
    Result = LoadData_function(WeightFile);
    if(strcmp(DATA_OPTION,'Peak1')||strcmp(DATA_OPTION,'Peak2')||strcmp(DATA_OPTION,'Nearest')||strcmp(DATA_OPTION,'Nearest4'))
        Sum_Weight_matrix = Result.Sum_Weight_matrix;
        Weight_matrix = Result.Weight_matrix;
    elseif(strcmp(DATA_OPTION,'Mean')||strcmp(DATA_OPTION,'Max')||strcmp(DATA_OPTION,'Min')||strcmp(DATA_OPTION,'Diff')||strcmp(DATA_OPTION,'Var'))
        Index_d = Result.Index_d;
        n = Result.n;
    end
else
    % raw data location list
    SiteData= LoadData_function([OUTPUTPATH_ROOT,SITENAME_DATA,Sep,'Location',Sep,SITENAME_DATA,'Site_',GCS,'.mat'],'SiteData');
    SiteData_Query = SiteData.SiteData;
    if(iscell(SiteData_Query))
        Lon_Query = cell2mat(SiteData_Query(:,3));
        Lat_Query = cell2mat(SiteData_Query(:,2));
    elseif(istable(SiteData_Query))
        Lon_Query = SiteData_Query.Lon;
        Lat_Query = SiteData_Query.Lat;
    else
        Lon_Query = SiteData_Query(:,3);
        Lat_Query = SiteData_Query(:,2);
    end
    N_Query = length(Lon_Query);
    N_Neighbour = min(N_Neighbour,N_Query);
    
    % first find neighbours.
    CurrentFile1 = [OUTPUTPATH_ROOT,SITENAME,Sep,'Temp',Sep,'NeighbourAdjacentList_',[DATA_OPTION,'_Thres',num2str(THRESHOLD)],'_',SITENAME_DATA,'_',SITENAME,'.mat'];
    if(exist(CurrentFile1,'file'))
        Result = LoadData_function(CurrentFile1);
        d = Result.d;
        n = Result.n;
    else
        fprintf('creating neighourhood files...\n');
        [n,d] = knnsearch([Lat_Query,Lon_Query],[Lat_Target,Lon_Target],'k',N_Neighbour);
        save(CurrentFile1,'n','d','-v7.3');
    end
    % create weight matrix -- different options, inverse-distnace^2, mean, etc. 
    d(d(:,1)==0,1)=0;%
    Index_d = d<THRESHOLD;
    if(strcmp(DATA_OPTION,'Peak1'))
        d(Index_d)=1./d(Index_d);%1/d;
        d(isinf(d)) = 0;
        d(~Index_d)=0;
        Weight_matrix = sparse(reshape(n,[N_Traget*N_Neighbour,1]),repmat((1:1:N_Traget)',[N_Neighbour,1]),reshape(d,[N_Traget*N_Neighbour,1]),N_Query,N_Traget);
        if(isinf(THRESHOLD))
            Weight_matrix = full(Weight_matrix);
        end
        Sum_Weight_matrix = sum(Weight_matrix,1);
        save(WeightFile,'Weight_matrix','Sum_Weight_matrix','-v7.3');
        clear d
    elseif(strcmp(DATA_OPTION,'Peak2'))
        d(Index_d)=1./(d(Index_d).^3);%1/d^2
        d(isinf(d)) = 0;
        d(~Index_d)=0;
        Weight_matrix = sparse(reshape(n,[N_Traget*N_Neighbour,1]),repmat((1:1:N_Traget)',[N_Neighbour,1]),reshape(d,[N_Traget*N_Neighbour,1]),N_Query,N_Traget);
        if(isinf(THRESHOLD))
            Weight_matrix = full(Weight_matrix);
        end
        Sum_Weight_matrix = sum(Weight_matrix,1);
        save(WeightFile,'Weight_matrix','Sum_Weight_matrix','-v7.3');
        clear d
    elseif(strcmp(DATA_OPTION,'Nearest4'))
        d(:,5:N_Neighbour) = 0;
        d(:,1)=1;
        d(~Index_d)=0;
        Weight_matrix = sparse(reshape(n,[N_Traget*N_Neighbour,1]),repmat((1:1:N_Traget)',[N_Neighbour,1]),reshape(d,[N_Traget*N_Neighbour,1]),N_Query,N_Traget);
        if(isinf(THRESHOLD))
            Weight_matrix = full(Weight_matrix);
        end
        Sum_Weight_matrix = sum(Weight_matrix,1);
        save(WeightFile,'Weight_matrix','Sum_Weight_matrix','-v7.3');
        clear d
    elseif(strcmp(DATA_OPTION,'Nearest'))
        d(:,2:N_Neighbour) = 0;
        d(:,1)=1;
        d(~Index_d)=0;
        Weight_matrix = sparse(reshape(n,[N_Traget*N_Neighbour,1]),repmat((1:1:N_Traget)',[N_Neighbour,1]),reshape(d,[N_Traget*N_Neighbour,1]),N_Query,N_Traget);
        if(isinf(THRESHOLD))
            Weight_matrix = full(Weight_matrix);
        end
        Sum_Weight_matrix = sum(Weight_matrix,1);
        save(WeightFile,'Weight_matrix','Sum_Weight_matrix','-v7.3');
        clear d
    elseif(strcmp(DATA_OPTION,'Mean')||strcmp(DATA_OPTION,'Max')||strcmp(DATA_OPTION,'Min')||strcmp(DATA_OPTION,'Diff')||strcmp(DATA_OPTION,'Var'))
        save(WeightFile,'Index_d','n','-v7.3');
    end
end

%  1 == EXTRAOPTION, we only process weight files, exit
if(1 == EXTRAOPTION)
    ReturnStatus = 1;
    return;
end

%% Load data, interpolation, and save outputs
TempInputFile = [DIRPATH,strrep(strrep(FileNameTemplate,'$STARTDATE$',num2str(INPUTYEAR)),'$ENDDATE$',num2str(INPUTYEAR))];

if(strcmp(OPTION,'By-Year'))
    TempOutput = TempFileName;
    if(exist(TempOutput,'file'))
        fprintf('skipping...%s\n',TempOutput);
    else
        % read data
        try
            fprintf('reading...%s\n',TempInputFile);
            Data = LoadData_function(TempInputFile,'Result');
            Data = Data.Result;
        catch exception
            fprintf('error reading...%s %s\n',TempInputFile,getReport(exception));
            ReturnStatus = 2;
            return;
        end
        % take time lag
        Data = LagData(Data,DATA_LAG);
        % process data
        Result = nan(yeardays(INPUTYEAR),N_Traget);
        if(strcmp(DATA_OPTION,'Peak1') || strcmp(DATA_OPTION,'Peak2') ||strcmp(DATA_OPTION,'Nearest') ||strcmp(DATA_OPTION,'Nearest4'))
            for i = 1:yeardays(INPUTYEAR)
%                 Result(i,:) = (Data(i,:)*Weight_matrix)./Sum_Weight_matrix;
                Result(i,:) = MultipleWeightMatrix_1(Weight_matrix',Data(i,:)')';
            end
        elseif(strcmp(DATA_OPTION,'Mean'))
            for i = 1:yeardays(INPUTYEAR)
%                 Result(i,:) = (Data(i,:)*Weight_matrix)./Sum_Weight_matrix;
                Temp = Data(n);
                Temp(~Index_d) = nan;
                Result(i,:) = nanmean(Temp,2)';
            end
%         elseif(strcmp(DATA_OPTION,'Max'))
%             Temp = Data(n);
%             Temp(~Index_d) = nan;
%             Result(i,:) = nanmax(Temp,[],2)';
%         elseif(strcmp(DATA_OPTION,'Min'))
%             Temp = Data(n);
%             Temp(~Index_d) = nan;
%             Result(i,:) = nanmin(Temp,[],2)';
%         elseif(strcmp(DATA_OPTION,'Diff'))
%             Temp = Data(n);
%             Temp(~Index_d) = nan;
%             Result(i,:) = nanmax(Temp,[],2)'-nanmin(Temp,[],2)';
%         elseif(strcmp(DATA_OPTION,'Var'))
%             Temp = Data(n);
%             Temp(~Index_d) = nan;
%             Result(i,:) = var(Temp,[],2,'omitnan')';
        end
        % save data
        fprintf('saving...%s\n',TempOutput);
        % visualize
        if(3 == EXTRAOPTION)
            Visualization_USResult_1(strrep(TempOutput,OUTPUTPATH,''),nanmean(Result),SiteData_Target,EnvironPara);
        end

        % save as matlab format or hdf5 format
        if(strcmp(FORMAT,'mat'))
            save(TempOutput,'Result');
        elseif(strcmp(FORMAT,'h5'))
            hdf5write(TempOutput, 'Result', Result);
        end
    end
    
elseif(strcmp(OPTION,'By-Day') || strcmp(OPTION,'By-Date'))
    % read data
    try
        fprintf('reading...%s\n',TempInputFile);
        DataAll = LoadData_function(TempInputFile,'Result');
        DataAll = DataAll.Result;
    catch exception
        fprintf('error reading...%s %s\n',TempInputFile,getReport(exception));
        ReturnStatus = 2;
        return;
    end
    
    %take time lag
    DataAll = LagData(DataAll,DATA_LAG);
    
    % process data
    for i = DATE_TO_PROCESS-DATE_TO_PROCESS(1)+1% make sure the i stands for day of year or day of TempFileName
        TempOutput = TempFileName{i};
        Data = DataAll(DATE_TO_PROCESS(i),:);
        if(exist(TempOutput,'file'))
            fprintf('skipping...%s\n',TempOutput);
        else
            if(strcmp(DATA_OPTION,'Mean') || strcmp(DATA_OPTION,'Peak1') || strcmp(DATA_OPTION,'Peak2') ||strcmp(DATA_OPTION,'Nearest')||strcmp(DATA_OPTION,'Nearest4'))
%                 Result = (Data*Weight_matrix)./Sum_Weight_matrix;
                Result = MultipleWeightMatrix_1(Weight_matrix',Data')';
            elseif(strcmp(DATA_OPTION,'Max'))
                Temp = Data(n);
                Temp(~Index_d) = nan;
                Result = nanmax(Temp,[],2)';
            elseif(strcmp(DATA_OPTION,'Min'))
                Temp = Data(n);
                Temp(~Index_d) = nan;
                Result = nanmin(Temp,[],2)';
            elseif(strcmp(DATA_OPTION,'Diff'))
                Temp = Data(n);
                Temp(~Index_d) = nan;
                Result = nanmax(Temp,[],2)'-nanmin(Temp,[],2)';
            elseif(strcmp(DATA_OPTION,'Var'))
                Temp = Data(n);
                Temp(~Index_d) = nan;
                Result = var(Temp,[],2,'omitnan')';
            end
            % save data
            fprintf('saving...%s\n',TempOutput);
            % visualize
            if(3 == EXTRAOPTION)
                Visualization_USResult_1(strrep(TempOutput,OUTPUTPATH,''),Result,SiteData_Target,EnvironPara);
            end
            
            % save as matlab format or hdf5 format
            if(strcmp(FORMAT,'mat'))
                save(TempOutput,'Result');
            elseif(strcmp(FORMAT,'h5'))
                hdf5write(TempOutput, 'Result', Result);
            end
        end
    end
end


ReturnStatus = 0;


