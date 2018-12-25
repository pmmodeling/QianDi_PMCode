%% main function to extract values from MOD11A1; it can also extract value from MAIAC AOD aggregate files
% require the grid cells files to be placed under ./processed_data/your_grid_cell_name/Location

%% version history
% 2016-07-27: read MOD11A1 processed output and interpolated them to Sites
% 2016-10-14: use projection system; add a mean convolutional layer
% 2016-12-14: generalized this code to process AOD data; and any other daily output; variable form: Data_xxxx; 
% can process multiple variable at once, but they has to have the same SiteData;

%% parameters
% SITENAME: name of the grid cells to be interpolated at;
% INPUTDATE: time period of interest;
% VARNAME: variable name, MOD11A1,MAIACUSTerra, MAIACUSAqua, MAIACUSTerra_cosVZA, MAIACUSAqua_cosVZA ;
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
% Interpolate_MOD11A1_Main('AQRVPM25',[datenum(2012,1,1),datenum(2012,12,31)],'air.sfc','By-Year','DailyMax',0)

%% code
function ReturnStatus = Interpolate_MOD11A1_Main(SITENAME,INPUTDATE,VARNAME,OPTION,DATA_OPTION,EXTRAOPTION,FORMAT)

%% marco-definition for each data set
if(strcmp('MOD11A1',VARNAME))
    VARNAME_LIST = {'LST_Day_1km','LST_Night_1km','Clear_day_cov','Clear_night_cov'};
    THRESHOLD = 7.5*1000;
    FileNameTemplate = 'MOD11A1_USMOD11A1_$STARTDATE$_$ENDDATE$.mat';
    OutputTemplate = ['MOD11A1_$NAME$_$OPTION$_$SITENAME$_$STARTDATE$_$ENDDATE$.',FORMAT];
    VARNAME_Simple = 'MOD11A1';
    SITENAME_DATA = 'USMOD11A1';
elseif(strcmp('MAIACUSAqua',VARNAME))
    VARNAME_LIST = {'Optical_Depth_047','Optical_Depth_055'};
    THRESHOLD = 7.5*1000;
    FileNameTemplate = 'MAIACUSAqua_USMAIACUS5km_$STARTDATE$_$ENDDATE$.mat';
    OutputTemplate = ['MAIACUSAqua_$NAME$_$OPTION$_$SITENAME$_$STARTDATE$_$ENDDATE$.',FORMAT];
    VARNAME_Simple = 'MAIACUS';
    SITENAME_DATA = 'USMAIACUS1km';
elseif(strcmp('MAIACUSTerra',VARNAME))
    VARNAME_LIST = {'Optical_Depth_047','Optical_Depth_055'};
    THRESHOLD = 7.5*1000;
    FileNameTemplate = 'MAIACUSTerra_USMAIACUS5km_$STARTDATE$_$ENDDATE$.mat';
    OutputTemplate = ['MAIACUSTerra_$NAME$_$OPTION$_$SITENAME$_$STARTDATE$_$ENDDATE$.',FORMAT];
    VARNAME_Simple = 'MAIACUS';
    SITENAME_DATA = 'USMAIACUS1km';
elseif(strcmp('MAIACUSTerra_cosVZA',VARNAME))
    VARNAME_LIST = {'cosVZA'};
    THRESHOLD = 7.5*1000;
    FileNameTemplate = 'MAIACUSTerra_USMAIACUS5km_$STARTDATE$_$ENDDATE$.mat';
    OutputTemplate = ['MAIACUSTerra_$NAME$_$OPTION$_$SITENAME$_$STARTDATE$_$ENDDATE$.',FORMAT];
    VARNAME_Simple = 'MAIACUS';
    SITENAME_DATA = 'USMAIACUS5km';
elseif(strcmp('MAIACUSAqua_cosVZA',VARNAME))
    VARNAME_LIST = {'cosVZA'};
    THRESHOLD = 7.5*1000;
    FileNameTemplate = 'MAIACUSAqua_USMAIACUS5km_$STARTDATE$_$ENDDATE$.mat';
    OutputTemplate = ['MAIACUSAqua_$NAME$_$OPTION$_$SITENAME$_$STARTDATE$_$ENDDATE$.',FORMAT];
    VARNAME_Simple = 'MAIACUS';
    SITENAME_DATA = 'USMAIACUS5km';
end

OutputTemplate = strrep(strrep(OutputTemplate,'$SITENAME$',SITENAME),'$OPTION$',[DATA_OPTION,'_Thres',num2str(THRESHOLD)]);

INPUTYEAR = year(INPUTDATE(1));
if(strcmp(OPTION,'By-Date'))
    DATE_TO_PROCESS = (INPUTDATE(1):1:INPUTDATE(2)) - datenum(INPUTYEAR,1,1)+1;
elseif(strcmp(OPTION,'By-Day')||strcmp(OPTION,'By-Year'))
    DATE_TO_PROCESS = 1:1:yeardays(INPUTYEAR);
end

%% path
DIRPATH = ['../data/aggregate/',VARNAME_Simple,'/'];
OUTPUTPATH_ROOT = '../../processed_data/';
Sep = '/';

% file path to the shape file used for making maps
EnvironPara = struct();
EnvironPara.states = geoshape(shaperead('usastatehi', 'UseGeoCoords', true));
EnvironPara.S_back = geoshape(shaperead('../data/shapefile/US_WGS_clipcoast_2.shp', 'UseGeoCoords', true));
EnvironPara.US_North_America_Equidistant_Conic = '../data/shapefile/US_North_America_Equidistant_Conic.shp';
EnvironPara.US = '../data/shapefile/US.shp';
GCS = 'North_America_Equidistant_Conic';%% projection system to be used

% output
OUTPUTPATH = [OUTPUTPATH_ROOT,SITENAME,Sep,VARNAME_Simple,Sep];
EnvironPara.OUTPUTPATH_PIC = OUTPUTPATH;
mkdir(OUTPUTPATH);

%% check whether has been processed before
if(strcmp(OPTION,'By-Year'))
    TempFileName = {};
    for l = 1:length(VARNAME_LIST)
        TempFileName = cat(1,TempFileName,[OUTPUTPATH, strrep(strrep(strrep(OutputTemplate,'$NAME$',VARNAME_LIST{l}),'$STARTDATE$',datestr(datenum(INPUTYEAR,1,1),'yyyymmdd')),'$ENDDATE$',datestr(datenum(INPUTYEAR,12,31),'yyyymmdd'))]);
    end
    
    if(2==EXTRAOPTION)
        cellfun(@(x) delete(x),TempFileName)
        fprintf('delete ...%s...%s\n',TempFileName{1},TempFileName{end});
        ReturnStatus = 1;
        return;
    else
        if(cellfun(@(x) exist(x,'file'),TempFileName))
            fprintf('SKIP ...%s...%s\n',TempFileName{1},TempFileName{end});
            ReturnStatus = 1;
            return;
        end
    end
elseif(strcmp(OPTION,'By-Day') || strcmp(OPTION,'By-Date'))
    StartDate = datenum(INPUTYEAR,1,1);
    
    TempFileName = {};
    for l = 1:length(VARNAME_LIST)
        TempFileName = cat(1,TempFileName,arrayfun(@(i)  ([OUTPUTPATH, strrep(strrep(strrep(OutputTemplate,'$NAME$',VARNAME_LIST{l}),'$STARTDATE$',datestr(StartDate+i-1,'yyyymmdd')),'$ENDDATE$',datestr(StartDate+i-1,'yyyymmdd'))]),DATE_TO_PROCESS,'UniformOutput',false));
    end
    
    if(2==EXTRAOPTION)
        for l = 1:length(VARNAME_LIST)
            for i=1:length(TempFileName)
                delete(TempFileName{l,i});
                fprintf('delete ...%s\n',TempFileName{l,i});
            end
        end
        ReturnStatus = 1;
        return;
    else
        if(cellfun(@(x) exist(x,'file'),TempFileName))%% all have been processed
            fprintf('SKIP! ...%s..%s\n',TempFileName{1,1},TempFileName{length(VARNAME_LIST),length(TempFileName)});
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
N_Neighbour = floor((THRESHOLD)^2*pi/816355);
WeightFile = [OUTPUTPATH_ROOT,SITENAME,Sep,'Temp',Sep,'NeighbourWeightMatrix_',[DATA_OPTION,'_Thres',num2str(THRESHOLD)],'_',SITENAME_DATA,'_',SITENAME,'.mat'];
if(exist(WeightFile,'file'))
    Result = LoadData_function(WeightFile);
    if(strcmp(DATA_OPTION,'Peak1')||strcmp(DATA_OPTION,'Peak2')||strcmp(DATA_OPTION,'Nearest')||strcmp(DATA_OPTION,'Nearest4')||strcmp(DATA_OPTION,'AnnualAverage'))
        Sum_Weight_matrix = Result.Sum_Weight_matrix;
        Weight_matrix = Result.Weight_matrix;
    elseif(strcmp(DATA_OPTION,'Mean')||strcmp(DATA_OPTION,'Max')||strcmp(DATA_OPTION,'Min')||strcmp(DATA_OPTION,'Diff')||strcmp(DATA_OPTION,'Var'))
        Index_d = Result.Index_d;
        n = Result.n;
    end
else
    %% raw data location list
    SiteData= LoadData_function([DIRPATH,SITENAME_DATA,'Site_',GCS,'.mat'],'SiteData');
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
    
    %% create weight matrix
    % this weight can be used for calculating weight average, maximal
    % values within certain threshold and so on
    d(d(:,1)==0,1)=0;%
    Index_d = d<THRESHOLD;
    if(strcmp(DATA_OPTION,'Peak1'))
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
    elseif(strcmp(DATA_OPTION,'Mean')||strcmp(DATA_OPTION,'Max')||strcmp(DATA_OPTION,'Min')||strcmp(DATA_OPTION,'Diff')||strcmp(DATA_OPTION,'Var'))
        save(WeightFile,'Index_d','n','-v7.3');
    end
    
end

%%  1 == EXTRAOPTION, we only process weight files
if(1 == EXTRAOPTION)
    ReturnStatus = 1;
    return;
end

%% output, by year or by day?
if(strcmp(OPTION,'By-Year'))
    Result = nan(yeardays(INPUTYEAR),length(Lon_Target),length(VARNAME_LIST));
elseif(strcmp(OPTION,'By-Day') || strcmp(OPTION,'By-Date'))
    Result = nan(1,length(Lon_Target));
end

%% read daily aggregate file and extract values
for i=DATE_TO_PROCESS-DATE_TO_PROCESS(1)+1% make sure the i stands for day of year or day of TempFileName
    
    CurrentDay = i+INPUTDATE(1)-1;
    % check whether exist; skip if exist
    if(strcmp(OPTION,'By-Day') || strcmp(OPTION,'By-Date'))
        if(cellfun(@(x) exist(x,'file'),TempFileName(:,i)))
            fprintf('skipping...%s\n',datestr(CurrentDay));
            
            continue;
        end
    end
    
    % read data
    TempFile = strrep(strrep(FileNameTemplate,'$STARTDATE$',datestr(CurrentDay,'yyyymmdd')),'$ENDDATE$',datestr(CurrentDay,'yyyymmdd'));
    if(exist([DIRPATH,TempFile],'file'))
        try
            DataAll = LoadData_function([DIRPATH,TempFile]);
            fprintf('reading...%s %s\n',datestr(CurrentDay,'yyyymmdd'),TempFile);
        catch
            fprintf('reading error...%s %s\n',datestr(CurrentDay,'yyyymmdd'),TempFile);
            ReturnStatus = 2;
            return;
        end
        
        % remove data based on QA flag -- no need; done that when
        % processing raw data
%         QAFlag = dec2bin(DataAll.Data_QC_Day,8);
%         TempFlag = bin2dec(QAFlag(:,7:8));% bit:00-01
%         DataAll.Data_LST_Day_1km(TempFlag==3|TempFlag==2)=nan;% remove "10 = LST not produced due to cloud effects" and "11 = LST not produced primarily due to reasons other than clouds" 
%         TempFlag = bin2dec(QAFlag(:,1:2));% bit:06-07
%         DataAll.Data_LST_Day_1km(TempFlag==3)=nan;%remove "11 = Average LST error > 3 K" 
%         
%         QAFlag = dec2bin(DataAll.Data_QC_Night,8);
%         TempFlag = bin2dec(QAFlag(:,7:8));% bit:00-01
%         DataAll.Data_LST_Night_1km(TempFlag==3|TempFlag==2)=nan;% remove "10 = LST not produced due to cloud effects" and "11 = LST not produced primarily due to reasons other than clouds" 
%         TempFlag = bin2dec(QAFlag(:,1:2));% bit:06-07
%         DataAll.Data_LST_Night_1km(TempFlag==3)=nan;%remove "11 = Average LST error > 3 K" 
        
        % process, extract value and match values from aggregate grid cells
        % to grid cells of interests;
        for l = 1:length(VARNAME_LIST)
            Data = eval(['DataAll.Data_',VARNAME_LIST{l}]);
            if(size(Data,1)==1)
                Data = Data';
            end
            if(strcmp(DATA_OPTION,'Mean'))
                Temp = nanmean(Data(n),2)';
            elseif(strcmp(DATA_OPTION,'Nearest')||strcmp(DATA_OPTION,'Nearest4')||strcmp(DATA_OPTION,'AnnualAverage'))
                Temp = MultipleWeightMatrix_1(Weight_matrix',Data)';
            end
            if(strcmp(OPTION,'By-Year'))
                Result(i,:,l) = Temp;
            elseif(strcmp(OPTION,'By-Day') || strcmp(OPTION,'By-Date'))
                Result = Temp;
                % visualize
                if(3 == EXTRAOPTION)
                    Visualization_USResult_1(strrep(TempFileName{l,i},OUTPUTPATH,''),Result,SiteData_Target,EnvironPara);
                end

                % save as matlab format or hdf5 format
                if(strcmp(FORMAT,'mat'))
                    save(TempFileName{l,i},'Result');
                elseif(strcmp(FORMAT,'h5'))
                    hdf5write(TempFileName{l,i}, 'Result', Result);
                end
                
            end
        end
    else
        fprintf('not exist...%s\n',datestr(CurrentDay,'yyyymmdd'));
        if(~strcmp(OPTION,'By-Year'))
            for l = 1:length(VARNAME_LIST)
                Result = nan(1,N_Traget); 

                % save as matlab format or hdf5 format
                if(strcmp(FORMAT,'mat'))
                    save(TempFileName{l,i},'Result');
                elseif(strcmp(FORMAT,'h5'))
                    hdf5write(TempFileName{l,i}, 'Result', Result);
                end
            end
            Temp = nan(1,length(Lon_Target));
        end
    end

    clear DataAll Data Temp
end


%% save results
if(strcmp(OPTION,'By-Year'))
    ResultAll = Result;
    for l = 1:length(VARNAME_LIST)
        Result = ResultAll(:,:,l);
        % visualize
        if(3 == EXTRAOPTION)
            Visualization_USResult_1(strrep(TempFileName{l},OUTPUTPATH,''),nanmean(Result),SiteData_Target,EnvironPara);
        end
        
        if(strcmp(DATA_OPTION,'AnnualAverage'))
            Result = nanmean(Result,1);
        end
        
        % save as matlab format or hdf5 format
        if(strcmp(FORMAT,'mat'))
            save(TempFileName{l},'Result');
        elseif(strcmp(FORMAT,'h5'))
            hdf5write(TempFileName{l}, 'Result', Result);
        end
    end
end

ReturnStatus = 0;
