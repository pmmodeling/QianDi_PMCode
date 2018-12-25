%%%% main function processing reanalysis data
%%% read 3-hour reanalysis data at vertical level and interpolate to locations
% require the grid cells files to be placed under ./processed_data/your_grid_cell_name/Location

%% version history
% 2016-08-14 skip all file that have been processed before
% 2016-10-02 this code can compute 1-day lag, 5-day lag now
% 2016-10-04: read pressure level data (multiple levels)
% 2016-10-14: remove filling value and missing values

%% parameters
% SITENAME: name of the grid cells to be interpolated at;
% INPUTDATE: time period of interest;
% VARNAME: meteorological variable name (here is only one option: omega);
% OPTION: output format: by day or by year;
% DATA_OPTION: what kind of output? daily mean, daily min, daily max, or lagged value?
% EXTRAOPTION: 1 --> only create weight files; 0--> interpolate the data; 2 --> delete all files; 3: create maps while doing interpolation
% FORMAT: 'mat': output file is matlab file; 'h5': output file is hdf file

%% return values:
% ReturnStatus==0: still in processing
% ReturnStatus==1: processing done before; 
% ReturnStatus==2, some raw data not available and needs to do it again;
% ReturnStatus==3, some input data not available has to quit;

%% example
% Interpolate_3hourPressureLevelReanalysis_Main('AQRVPM25',[datenum(2012,1,1),datenum(2012,12,31)],'omega','By-Year','DailyMax',0)

%% code
function ReturnStatus = Interpolate_3hourPressureLevelReanalysis_Main(SITENAME,INPUTDATE,VARNAME,OPTION,DATA_OPTION,EXTRAOPTION,FORMAT)

%% macro definition
GCS = 'North_America_Equidistant_Conic';%% projection system to be used
VARNAME_LIST =  {'omega'};
VARNAME_FIELD_LIST =  {'omega'};
[~,locb] = ismember(VARNAME,VARNAME_LIST);
VARNAME_FIELD = VARNAME_FIELD_LIST{locb};

% output format, output file saved by day or year?
INPUTYEAR = year(INPUTDATE(1));
if(strcmp(OPTION,'By-Date'))   
    DATE_TO_PROCESS = (INPUTDATE(1):1:INPUTDATE(2)) - datenum(INPUTYEAR,1,1)+1;
elseif(strcmp(OPTION,'By-Day')||strcmp(OPTION,'By-Year'))   
    DATE_TO_PROCESS = 1:1:yeardays(INPUTYEAR);
end

% output path
DIRPATHNC_ROOT = '../data/unprocessed/Data_3HourPressureLevelNCEP/';
OUTPUTPATH_ROOT = '../../processed_data/';
Sep = '/';

% file path to the shape file used for making maps
EnvironPara = struct();
EnvironPara.states = geoshape(shaperead('usastatehi', 'UseGeoCoords', true));
EnvironPara.S_back = geoshape(shaperead('../data/shapefile/US_WGS_clipcoast_2.shp', 'UseGeoCoords', true));
EnvironPara.US_North_America_Equidistant_Conic = '../data/shapefile/US_North_America_Equidistant_Conic.shp';
EnvironPara.US = '../data/shapefile/US.shp';
GCS = 'North_America_Equidistant_Conic';

OutputTemplate = ['REANALYSIS_$VARNAME$_$OPTION$_$SITENAME$_$STARTDATE$_$ENDDATE$.',FORMAT];
OutputTemplate = strrep(strrep(strrep(OutputTemplate,'$VARNAME$',VARNAME),'$OPTION$',DATA_OPTION),'$SITENAME$',SITENAME);


%% skip if output files exist
OUTPUTPATH = [OUTPUTPATH_ROOT,SITENAME,Sep,'MeteFields',Sep];
EnvironPara.OUTPUTPATH_PIC = OUTPUTPATH;
mkdir(OUTPUTPATH);
% check whether has been processed before or remove files
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
elseif(strcmp(OPTION,'By-Day')||strcmp(OPTION,'By-Date'))
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


%% current file
DIRPATHNC = [DIRPATHNC_ROOT,VARNAME,Sep];

fprintf('%s\n',SITENAME);
fprintf('%d\n',INPUTYEAR);
fprintf('%s\n',VARNAME);

%% read points of interests
% these are the grid cells we will interpolate at
SiteData = LoadData_function([OUTPUTPATH_ROOT,SITENAME,Sep,'Location',Sep,SITENAME,'Site_',GCS,'.mat'],'SiteData');
SiteData_Target = SiteData.SiteData;
if(iscell(SiteData))
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

%% read raw data and grid location of raw data
% read raw data -- monthly data
TempTime = [];
TempData = [];
for i=1:12
    try
        FileName = [DIRPATHNC_ROOT,VARNAME,Sep,VARNAME,'.',num2str(INPUTYEAR),sprintf('%02d',i),'.nc'];
        fprintf('%s\n',FileName);
        FileInfo = ncinfo(FileName);
        TempTime_temp = double(ncread(FileName,'time'));%size:248*1
        TempTime = cat(1,TempTime,TempTime_temp);
        
        TempData_temp = double(ncread(FileName,VARNAME_FIELD));%size: 349   277    29   248
        TempData_temp = TempData_temp(:,:,1,:);% take the first layer near the ground
        pack;
        TempData_temp = reshape(TempData_temp,[size(TempData_temp,1),size(TempData_temp,2),size(TempData_temp,4)]);%size:349   277   248
        TempData_temp =  Interpolate_3hourReanalysis_function(FileInfo,TempData_temp,VARNAME_FIELD);
        TempData = cat(3,TempData,TempData_temp);
    catch exception
        fprintf('error reading...%s %s\n',FileName,getReport(exception));
        ReturnStatus = 2;
        return;  
    end
end

% read lat and lon of raw data
SiteData= LoadData_function([DIRPATHNC_ROOT,'Data_3HourNCEPSite_',GCS,'.mat'],'SiteData');
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

% we need also to read previous-year data to calculate n-day lagged values;
if(~isempty(strfind(DATA_OPTION,'Day')))
    for i=flip(12:12)% only read previous December
        try
            FileName = [DIRPATHNC_ROOT,VARNAME,Sep,VARNAME,'.',num2str(INPUTYEAR-1),sprintf('%02d',i),'.nc'];
            fprintf('%s\n',FileName);
            FileInfo = ncinfo(FileName);
            TempTime_temp = double(ncread(FileName,'time'));%size:248*1
            TempTime = cat(1,TempTime_temp,TempTime);
            
            TempData_temp = double(ncread(FileName,VARNAME_FIELD));%size: 349   277    29   248
            TempData_temp = TempData_temp(:,:,1,:);% take the first layer near the ground
            TempData_temp = reshape(TempData_temp,[size(TempData_temp,1),size(TempData_temp,2),size(TempData_temp,4)]);%size:349   277   248
            % remove invalid values
            TempData_temp =  Interpolate_3hourReanalysis_function(FileInfo,TempData_temp,VARNAME_FIELD);
            TempData = cat(3,TempData_temp,TempData);
        catch exception
            fprintf('error reading...%s %s\n',FileName,getReport(exception));
            ReturnStatus = 2;
            return;
        end
    end
end

%% calculate weight --- used for interpolation (weighted average of nearby grid cells)
mkdir([OUTPUTPATH_ROOT,SITENAME,Sep,'Temp',Sep]);
THRESHOLD = 50*1000;
N_Neighbour = 4;
WeightFile = [OUTPUTPATH_ROOT,SITENAME,Sep,'Temp',Sep,'NeighbourWeightMatrix_',['Nearest4','_Thres',num2str(THRESHOLD)],'_','Data_3HourNCEPSite','_',SITENAME,'.mat'];
if(exist(WeightFile,'file'))
    Result = LoadData_function(WeightFile);
    Sum_Weight_matrix = Result.Sum_Weight_matrix;
    Weight_matrix = Result.Weight_matrix;
else
    %% first find neighbours.
    CurrentFile1 = [OUTPUTPATH_ROOT,SITENAME,Sep,'Temp',Sep,'NeighbourAdjacentList_',['Nearest4','_Thres',num2str(THRESHOLD)],'_','Data_3HourNCEPSite','_',SITENAME,'.mat'];
    if(exist(CurrentFile1,'file'))
        Result = LoadData_function(CurrentFile1);
        d = Result.d;
        n = Result.n;
    else
        fprintf('creating neighourhood files...\n');
        [n,d] = knnsearch([Lat_Query,Lon_Query],[Lat_Target,Lon_Target],'k',N_Neighbour); 
        save(CurrentFile1,'n','d','-v7.3');
    end
    %% create weight matrix
    fprintf('creating weight matrix, Option:%s...\n','Nearest4');
    d(d(:,1)==0,1)=0;%    
    Index_d = d<THRESHOLD;
    d(Index_d)=1;%equal weight...;
    d(~Index_d)=0;
    Weight_matrix = sparse(reshape(n,[N_Traget*N_Neighbour,1]),repmat((1:1:N_Traget)',[N_Neighbour,1]),reshape(d,[N_Traget*N_Neighbour,1]),N_Query,N_Traget);
    Sum_Weight_matrix = sum(Weight_matrix,1);
    save(WeightFile,'Weight_matrix','Sum_Weight_matrix','-v7.3');
end

%  1 == EXTRAOPTION, we only process weight files, exit
if(1 == EXTRAOPTION)
    ReturnStatus = 1;
    return;
end

% check for possible erros if the raw reanalysis file has a different format
Date = floor(TempTime/24);
if(Date(1)~=Date(8))
    disk('error!!! only process 3-hour reanalysis data');
    ReturnStatus = 2;
    return;
end

%% calculate daily average, daily max, daily min, lagged values
% only process 3-hour data!!!! -- find daily max, daily min or daily mean
TempResult = nan(size(TempData,1),size(TempData,2),size(TempData,3)/8);
fprintf('processing...taking...%s\n',DATA_OPTION);
for i = 1:size(TempResult,3)
    Temp = TempData(:,:,((i-1)*8+1):(i*8));
    if(strcmp(DATA_OPTION,'DailyMean')  || strcmp(DATA_OPTION,'AnnualAverage'))
        TempResult(:,:,i) = nanmean(Temp,3); 
    elseif(strcmp(DATA_OPTION,'DailyMax'))
        TempResult(:,:,i) = nanmax(Temp,[],3); 
    elseif( strcmp(DATA_OPTION,'DailyMin'))
        TempResult(:,:,i) = nanmin(Temp,[],3); 
    elseif(~isempty(strfind(DATA_OPTION,'Day')))% for temporal convolutional layer, it is based on mean values from previous days
        TempResult(:,:,i) = nanmean(Temp,3); 
    end
end


% reshape and interpolation and save
TempResult = reshape(TempResult,[size(TempResult,1)*size(TempResult,2),size(TempResult,3)])';
% TempData = reshape(TempData,[size(TempData,1), size(TempData,2)*size(TempData,3)]);

%% calculate lagged values
if(~isempty(strfind(DATA_OPTION,'Day')))
    fprintf('creating temporal convolutional layer...%s\n',DATA_OPTION);
    DayOffset = str2double(strrep(DATA_OPTION,'Day',''));
    if(~isnumeric(DayOffset))
        disk('error!!! DATA_OPTION must be format like "5Day"!!!');
        ReturnStatus = 2;
        return;
    end
    TemporalSize = yeardays(INPUTYEAR)+31;% only takes December
    TemporalWeight = tril(ones(TemporalSize,TemporalSize),-1)-tril(ones(TemporalSize,TemporalSize),-(1+DayOffset));
    TemporalWeight = TemporalWeight/DayOffset;
    TempResult = TemporalWeight*TempResult;
    % take only this year...
    TempResult = TempResult((31+1):(yeardays(INPUTYEAR)+31),:);% only takes to Previous December
end

%% interpolation and save output
if(strcmp(OPTION,'By-Year'))  
    Result = nan(yeardays(INPUTYEAR),N_Traget);
%     Result =
%     InterpMyData_2(TempResult,TempLon,TempLat,Lon_Target,Lat_Target,'');%too slow...
    for i = 1:yeardays(INPUTYEAR)
        Result(i,:) =  (TempResult(i,:)*Weight_matrix)./Sum_Weight_matrix;
    end
    
    % visualize
    if(3 == EXTRAOPTION)
        Visualization_USResult_1(strrep(TempFileName,OUTPUTPATH,''),nanmean(Result),SiteData_Target,EnvironPara);
    end
    
    if(strcmp(DATA_OPTION,'AnnualAverage'))
        Result = nanmean(Result,1);
    end
    
    fprintf('saving...%d\n',INPUTYEAR);

    % save as matlab format or hdf5 format
    if(strcmp(FORMAT,'mat'))
        save(TempFileName,'Result');
    elseif(strcmp(FORMAT,'h5'))
        hdf5write(TempFileName, 'Result', Result);
    end
elseif(strcmp(OPTION,'By-Day')||strcmp(OPTION,'By-Date'))
    for i = DATE_TO_PROCESS-DATE_TO_PROCESS(1)+1% make sure the i stands for day of year or day of TempFileName
        CurrentDay = INPUTDATE(1)+i-1;
        if(exist(TempFileName{i},'file'))
            fprintf('skipping...%s\n',datestr(CurrentDay));
            continue;
        end
%         Result = InterpMyData_2(TempResult(i,:),TempLon,TempLat,Lon_Target,Lat_Target,'');%too slow...
        Result =  (TempResult(DATE_TO_PROCESS(i),:)*Weight_matrix)./Sum_Weight_matrix;
        fprintf('saving...%s\n',datestr(CurrentDay));
        % visualize
        if(3 == EXTRAOPTION)
            Visualization_USResult_1(strrep(TempFileName{i},OUTPUTPATH,''),Result,SiteData_Target,EnvironPara);
        end
        
        % save as matlab format or hdf5 format
        if(strcmp(FORMAT,'mat'))
            save(TempFileName{i},'Result');
        elseif(strcmp(FORMAT,'h5'))
            hdf5write(TempFileName{i}, 'Result', Result);
        end
    end
end

ReturnStatus = 0;

