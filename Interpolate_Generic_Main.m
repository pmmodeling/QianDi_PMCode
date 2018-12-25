%% extract data from MOD13A2 and MOD09A1 MODIS data
% although this code can process GMTED2010 files, but it is no longer used
% require the grid cells files to be placed under ./processed_data/your_grid_cell_name/Location

%% parameters
% SITENAME: name of the grid cells to be interpolated at;
% INPUTDATE: time period of interest;
% VARNAME: data set to be used;
% OPTION: output format: by day or by year;
% DATA_OPTION: what kind of output? min, max, diff, etc.
% EXTRAOPTION: 1 --> only create weight files; 0--> interpolate the data; 2 --> delete all files; 3: create maps while doing interpolation
% FORMAT: 'mat': output file is matlab file; 'h5': output file is hdf file

%% return values:
% ReturnStatus==0: still in processing
% ReturnStatus==1: processing done before; 
% ReturnStatus==2, some raw data not available and needs to do it again;
% ReturnStatus==3, some input data not available has to quit;

%% example:
% Interpolate_Generic_Main('AQRVPM25',[datenum(2012,1,1),datenum(2012,12,31)],'air.sfc','By-Year','DailyMax',0)

%% code
function ReturnStatus = Interpolate_Generic_Main(SITENAME,INPUTDATE,VARNAME,OPTION,DATA_OPTION,EXTRAOPTION,FORMAT)

%% marco
START_DATE = datenum(2000,1,1);
END_DATE= datenum(2017,1,1);

INPUTYEAR = year(INPUTDATE(1));
if(strcmp(OPTION,'By-Date'))   
    DATE_TO_PROCESS = (INPUTDATE(1):1:INPUTDATE(2)) - datenum(INPUTYEAR,1,1)+1;
elseif(strcmp(OPTION,'By-Day')||strcmp(OPTION,'By-Year'))   
    DATE_TO_PROCESS = 1:1:yeardays(INPUTYEAR);
end

%% specify distance threshold and the nubmer of neighboring cells
if(strcmp(VARNAME,'MOD13A2'))
    THRESHOLD = 7.5*1000;
    FileNameTemplate = 'MOD13A2_$STARTDATE$_$ENDDATE$.mat';
    OutputTemplate = ['MOD13A2_$OPTION$_$SITENAME$_$STARTDATE$_$ENDDATE$.',FORMAT];
    N_Neighbour = floor((THRESHOLD)^2*pi/809081);
    VARNAME_Simple = 'MOD13A2';
elseif(strcmp(VARNAME,'MOD09A1'))
    THRESHOLD = 7.5*1000;
    FileNameTemplate = 'MOD09A1_$STARTDATE$_$ENDDATE$.mat';
    OutputTemplate = ['MOD09A1_$OPTION$_$SITENAME$_$STARTDATE$_$ENDDATE$.',FORMAT];
    N_Neighbour = floor((THRESHOLD)^2*pi/809081);
    VARNAME_Simple = 'MOD09A1';
elseif(strcmp(VARNAME,'GMTED2010med150'))
    % the "area" of a grid cell in GMTED2010Arc300 is 152920 meters 
    THRESHOLD = 50*1000;
    FileNameTemplate = 'GMTED2010med_Evelation_USGMTED2010Arc150.mat';
    OutputTemplate = ['GMTED2010med_$OPTION$_$SITENAME$.',FORMAT];
    N_Neighbour = floor((THRESHOLD)^2*pi/152920);
    VARNAME_Simple = 'GMTED2010Arc150';
elseif(strcmp(VARNAME,'GMTED2010max150'))
    THRESHOLD = 50*1000;
    FileNameTemplate = 'GMTED2010max_Evelation_USGMTED2010Arc150.mat';
    OutputTemplate = ['GMTED2010max_$OPTION$_$SITENAME$.',FORMAT];
    N_Neighbour = floor((THRESHOLD)^2*pi/152920);
    VARNAME_Simple = 'GMTED2010Arc150';
elseif(strcmp(VARNAME,'GMTED2010min150'))
    THRESHOLD = 50*1000;
    FileNameTemplate = 'GMTED2010min_Evelation_USGMTED2010Arc150.mat';
    OutputTemplate = ['GMTED2010min_$OPTION$_$SITENAME$.',FORMAT];
    N_Neighbour = floor((THRESHOLD)^2*pi/152920);
    VARNAME_Simple = 'GMTED2010Arc150';
elseif(strcmp(VARNAME,'GMTED2010std150'))
    THRESHOLD = 50*1000;
    FileNameTemplate = 'GMTED2010std_Evelation_USGMTED2010Arc150.mat';
    OutputTemplate = ['GMTED2010std_$OPTION$_$SITENAME$.',FORMAT];
    N_Neighbour = floor((THRESHOLD)^2*pi/152920);
    VARNAME_Simple = 'GMTED2010Arc150';
elseif(strcmp(VARNAME,'GMTED2010mea150'))
    THRESHOLD = 50*1000;
    FileNameTemplate = 'GMTED2010mea_Evelation_USGMTED2010Arc150.mat';
    OutputTemplate = ['GMTED2010mea_$OPTION$_$SITENAME$.',FORMAT];
    N_Neighbour = floor((THRESHOLD)^2*pi/152920);
    VARNAME_Simple = 'GMTED2010Arc150';
elseif(strcmp(VARNAME,'GMTED2010bln150'))
    THRESHOLD = 50*1000;
    FileNameTemplate = 'GMTED2010bln_Evelation_USGMTED2010Arc150.mat';
    OutputTemplate = ['GMTED2010bln_$OPTION$_$SITENAME$.',FORMAT];
    N_Neighbour = floor((THRESHOLD)^2*pi/152920);
    VARNAME_Simple = 'GMTED2010Arc150';
    
elseif(strcmp(VARNAME,'GMTED2010med300'))
    % the "area" of a grid cell in GMTED2010Arc300 is 611681 meters 
    THRESHOLD = 7.5*1000;
    FileNameTemplate = 'GMTED2010med_Evelation_USGMTED2010Arc300.mat';
    OutputTemplate = ['GMTED2010med_$OPTION$_$SITENAME$.',FORMAT];
    N_Neighbour = floor((THRESHOLD)^2*pi/611681);
    VARNAME_Simple = 'GMTED2010Arc300';
elseif(strcmp(VARNAME,'GMTED2010max300'))
    THRESHOLD = 7.5*1000;
    FileNameTemplate = 'GMTED2010max_Evelation_USGMTED2010Arc300.mat';
    OutputTemplate = ['GMTED2010max_$OPTION$_$SITENAME$.',FORMAT];
    N_Neighbour = floor((THRESHOLD)^2*pi/611681);
    VARNAME_Simple = 'GMTED2010Arc300';
elseif(strcmp(VARNAME,'GMTED2010min300'))
    THRESHOLD = 7.5*1000;
    FileNameTemplate = 'GMTED2010min_Evelation_USGMTED2010Arc300.mat';
    OutputTemplate = ['GMTED2010min_$OPTION$_$SITENAME$.',FORMAT];
    N_Neighbour = floor((THRESHOLD)^2*pi/611681);
    VARNAME_Simple = 'GMTED2010Arc300';
elseif(strcmp(VARNAME,'GMTED2010std300'))
    THRESHOLD = 7.5*1000;
    FileNameTemplate = 'GMTED2010std_Evelation_USGMTED2010Arc300.mat';
    OutputTemplate = ['GMTED2010std_$OPTION$_$SITENAME$.',FORMAT];
    N_Neighbour = floor((THRESHOLD)^2*pi/611681);
    VARNAME_Simple = 'GMTED2010Arc300';
elseif(strcmp(VARNAME,'GMTED2010mea300'))
    THRESHOLD = 7.5*1000;
    FileNameTemplate = 'GMTED2010mea_Evelation_USGMTED2010Arc300.mat';
    OutputTemplate = ['GMTED2010mea_$OPTION$_$SITENAME$.',FORMAT];
    N_Neighbour = floor((THRESHOLD)^2*pi/611681);
    VARNAME_Simple = 'GMTED2010Arc300';
elseif(strcmp(VARNAME,'GMTED2010bln300'))
    THRESHOLD = 7.5*1000;
    FileNameTemplate = 'GMTED2010bln_Evelation_USGMTED2010Arc300.mat';
    OutputTemplate = ['GMTED2010bln_$OPTION$_$SITENAME$.',FORMAT];
    N_Neighbour = floor((THRESHOLD)^2*pi/611681);
    VARNAME_Simple = 'GMTED2010Arc300';
end

DIRPATH = ['../data/aggregate/',VARNAME_Simple,'/'];
OUTPUTPATH_ROOT = '../../processed_data/';
Sep = '/';

%% shapefile path, used to make maps
EnvironPara = struct();
EnvironPara.states = geoshape(shaperead('usastatehi', 'UseGeoCoords', true));
EnvironPara.S_back = geoshape(shaperead('../data/shapefile/US_WGS_clipcoast_2.shp', 'UseGeoCoords', true));
EnvironPara.US_North_America_Equidistant_Conic = '../data/shapefile/US_North_America_Equidistant_Conic.shp';
EnvironPara.US = '../data/shapefile/US.shp';
GCS = 'North_America_Equidistant_Conic';%% projection system to be used

OutputTemplate = strrep(strrep(OutputTemplate,'$OPTION$',[DATA_OPTION,'_Thres',num2str(THRESHOLD)]),'$SITENAME$',SITENAME);

%% output
OUTPUTPATH = [OUTPUTPATH_ROOT,SITENAME,Sep,VARNAME_Simple,Sep];
EnvironPara.OUTPUTPATH_PIC = OUTPUTPATH;
mkdir(OUTPUTPATH);
%% check whether has been processed before
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

%% interpolation location list
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

%% calculate weight for interpolation and save it for future use
mkdir([OUTPUTPATH_ROOT,SITENAME,Sep,'Temp',Sep]);
WeightFile = [OUTPUTPATH_ROOT,SITENAME,Sep,'Temp',Sep,'NeighbourWeightMatrix_',[DATA_OPTION,'_Thres',num2str(THRESHOLD)],'_','US',VARNAME_Simple,'_',SITENAME,'.mat'];
if(exist(WeightFile,'file'))
    Result = LoadData_function(WeightFile);
    if(strcmp(DATA_OPTION,'Mean')||strcmp(DATA_OPTION,'Peak1')||strcmp(DATA_OPTION,'Peak2')||strcmp(DATA_OPTION,'Nearest')||strcmp(DATA_OPTION,'Nearest4')||strcmp(DATA_OPTION,'AnnualAverage'))
        Sum_Weight_matrix = Result.Sum_Weight_matrix;
        Weight_matrix = Result.Weight_matrix;
    elseif(strcmp(DATA_OPTION,'Max')||strcmp(DATA_OPTION,'Min')||strcmp(DATA_OPTION,'Diff')||strcmp(DATA_OPTION,'Var'))
        Index_d = Result.Index_d;
        n = Result.n;
    end
else
    %% raw data location list
    SiteData= LoadData_function([DIRPATH,'US',VARNAME_Simple,'Site_',GCS,'.mat'],'SiteData');
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
    
    %% first find neighbours.
    CurrentFile1 = [OUTPUTPATH_ROOT,SITENAME,Sep,'Temp',Sep,'NeighbourAdjacentList_',[DATA_OPTION,'_Thres',num2str(THRESHOLD)],'_','US',VARNAME_Simple,'_',SITENAME,'.mat'];
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
    % this weight can be used for calculating weight average, maximal
    % values within certain threshold and so on
    d(d(:,1)==0,1)=0;%
    Index_d = d<THRESHOLD;
    if(strcmp(DATA_OPTION,'Mean'))
        d(Index_d)=1;%equal weight...;
        d(~Index_d)=0;
        Weight_matrix = sparse(reshape(n,[N_Traget*N_Neighbour,1]),repmat((1:1:N_Traget)',[N_Neighbour,1]),reshape(d,[N_Traget*N_Neighbour,1]),N_Query,N_Traget);
        Sum_Weight_matrix = sum(Weight_matrix,1);
        save(WeightFile,'Weight_matrix','Sum_Weight_matrix','-v7.3');
    elseif(strcmp(DATA_OPTION,'Peak1'))
        d(Index_d)=1./d(Index_d);%1/d;
        d(isinf(d)) = 0;
        d(~Index_d)=0;
        Weight_matrix = sparse(reshape(n,[N_Traget*N_Neighbour,1]),repmat((1:1:N_Traget)',[N_Neighbour,1]),reshape(d,[N_Traget*N_Neighbour,1]),N_Query,N_Traget);
        Sum_Weight_matrix = sum(Weight_matrix,1);
        save(WeightFile,'Weight_matrix','Sum_Weight_matrix','-v7.3');
    elseif(strcmp(DATA_OPTION,'Peak2'))
        d(Index_d)=1./(d(Index_d).*d(Index_d));%1/d^2
        d(isinf(d)) = 0;
        d(~Index_d)=0;
        Weight_matrix = sparse(reshape(n,[N_Traget*N_Neighbour,1]),repmat((1:1:N_Traget)',[N_Neighbour,1]),reshape(d,[N_Traget*N_Neighbour,1]),N_Query,N_Traget);
        Sum_Weight_matrix = sum(Weight_matrix,1);
        save(WeightFile,'Weight_matrix','Sum_Weight_matrix','-v7.3');
    elseif(strcmp(DATA_OPTION,'Nearest4')||strcmp(DATA_OPTION,'AnnualAverage'))
        d(:,5:N_Neighbour) = 0;
        d(:,1)=1;
        d(~Index_d)=0;
        Weight_matrix = sparse(reshape(n,[N_Traget*N_Neighbour,1]),repmat((1:1:N_Traget)',[N_Neighbour,1]),reshape(d,[N_Traget*N_Neighbour,1]),N_Query,N_Traget);
        Sum_Weight_matrix = sum(Weight_matrix,1);
        save(WeightFile,'Weight_matrix','Sum_Weight_matrix','-v7.3');
    elseif(strcmp(DATA_OPTION,'Nearest'))
        d(:,2:N_Neighbour) = 0;
        d(:,1)=1;
        d(~Index_d)=0;
        Weight_matrix = sparse(reshape(n,[N_Traget*N_Neighbour,1]),repmat((1:1:N_Traget)',[N_Neighbour,1]),reshape(d,[N_Traget*N_Neighbour,1]),N_Query,N_Traget);
        Sum_Weight_matrix = sum(Weight_matrix,1);
        save(WeightFile,'Weight_matrix','Sum_Weight_matrix','-v7.3');
    elseif(strcmp(DATA_OPTION,'Max')||strcmp(DATA_OPTION,'Min')||strcmp(DATA_OPTION,'Diff')||strcmp(DATA_OPTION,'Var'))
        save(WeightFile,'Index_d','n','-v7.3');
    end
    
    
end

%%  1 == EXTRAOPTION, we only process weight files
if(1 == EXTRAOPTION)
    ReturnStatus = 1;
    return;
end

%% read daily files
Result = nan(yeardays(INPUTYEAR),length(Lon_Target));
for i=DATE_TO_PROCESS-DATE_TO_PROCESS(1)+1
    % make sure the i stands for day of year or day of TempFileName
    
    CurrentDay = INPUTDATE(1)+i-1;
    if(strcmp(OPTION,'By-Day')||strcmp(OPTION,'By-Date'))
        if(exist(TempFileName{i},'file'))
            fprintf('skipping...%s\n',datestr(CurrentDay));
            continue;
        end
    end
    fprintf('%s\n',datestr(CurrentDay));
    
    % read daily file a particular day; if it does not exist, read
    % neighouring days and do linear interpolation
    TempFile = strrep(strrep(FileNameTemplate,'$STARTDATE$',datestr(CurrentDay,'yyyymmdd')),'$ENDDATE$',datestr(CurrentDay,'yyyymmdd'));
    try
        Data = LoadData_function([DIRPATH,TempFile],'Result');
        fprintf('reading...%s\n',[DIRPATH,TempFile]);
        Data = Data.Result;
    catch
        %         CurrentDayA = datenum(2000,2,18);% the first day of MOD13A2 is 2000-02-18
        for j=CurrentDay:-1:START_DATE
            TempFileA = strrep(strrep(FileNameTemplate,'$STARTDATE$',datestr(j,'yyyymmdd')),'$ENDDATE$',datestr(j,'yyyymmdd'));
            try
                DataA = LoadData_function([DIRPATH,TempFileA],'Result');
                fprintf('\treading A...%s\n',datestr(j,'yyyymmdd'));
                DataA = DataA.Result;
                CurrentDayA = j;
                break;
            catch
                CurrentDayA = j;
                continue;
            end
        end
        
        %         CurrentDayB = datenum(2015,12,19);% the first day of MOD13A2 is 2000-02-18
        for j=CurrentDay:1:END_DATE
            TempFileB = strrep(strrep(FileNameTemplate,'$STARTDATE$',datestr(j,'yyyymmdd')),'$ENDDATE$',datestr(j,'yyyymmdd'));
            try
                DataB = LoadData_function([DIRPATH,TempFileB],'Result');
                fprintf('\treading B...%s\n',datestr(j,'yyyymmdd'));
                DataB = DataB.Result;
                CurrentDayB = j;
                break;
            catch
                CurrentDayB = j;
                continue;
            end
        end
        
        % do linear interpolation
        fprintf('temporal interpolating...%s\n',datestr(CurrentDay,'yyyymmdd'));
        if(CurrentDayB == END_DATE)%% find not data afterwards
            Data = DataA;
        elseif(CurrentDayA == START_DATE)
            Data = DataB;
        else
            Data = (DataA*abs(CurrentDayB - CurrentDay)+DataB*abs(CurrentDayA - CurrentDay))/abs(CurrentDayA-CurrentDayB);
        end
    end
    
    % match from aggregate grid cells to grid cells of interests; we can
    % choose how to match, based on the mean of four nearest grids? of max
    % of the four nearest grids? or min?
    if(strcmp(DATA_OPTION,'Mean') || strcmp(DATA_OPTION,'Peak1') || strcmp(DATA_OPTION,'Peak2') ||strcmp(DATA_OPTION,'Nearest') ||strcmp(DATA_OPTION,'Nearest4')||strcmp(DATA_OPTION,'AnnualAverage'))
        Temp =  MultipleWeightMatrix_1(Weight_matrix',Data)';%(Data*Weight_matrix)./Sum_Weight_matrix;
    elseif(strcmp(DATA_OPTION,'Max'))
        Temp = Data(n);
        Temp(~Index_d) = nan;
        Temp = nanmax(Temp,[],2)';
    elseif(strcmp(DATA_OPTION,'Min'))
        Temp = Data(n);
        Temp(~Index_d) = nan;
        Temp = nanmin(Temp,[],2)';
    elseif(strcmp(DATA_OPTION,'Diff'))
        Temp = Data(n);
        Temp(~Index_d) = nan;
        Temp = nanmax(Temp,[],2)'-nanmin(Temp,[],2)';
    elseif(strcmp(DATA_OPTION,'Var')) 
        Temp = Data(n);
        Temp(~Index_d) = nan;
        Temp = var(Temp,[],2,'omitnan')';
    end
    
    if(strcmp(OPTION,'By-Year'))
        Result(i,:) = Temp;
    elseif(strcmp(OPTION,'By-Day')||strcmp(OPTION,'By-Date'))
        Result = Temp;
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

%% save
% save the data either by year of by day
if(strcmp(OPTION,'By-Year'))
    % visualize
    if(3 == EXTRAOPTION)
        Visualization_USResult_1(strrep(TempFileName,OUTPUTPATH,''),nanmean(Result),SiteData_Target,EnvironPara);
    end
    
    if(strcmp(DATA_OPTION,'AnnualAverage'))
        Result = nanmean(Result,1);
    end
    
    % save as matlab format or hdf5 format
    if(strcmp(FORMAT,'mat'))
        save(TempFileName,'Result');
    elseif(strcmp(FORMAT,'h5'))
        hdf5write(TempFileName, 'Result', Result);
    end
end

ReturnStatus = 0;





