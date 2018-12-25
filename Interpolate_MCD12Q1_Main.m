%% extract value for grid cells of interests from all types of MCD12Q1

%% parameters
% SITENAME: name of the grid cells to be interpolated at;
% INPUTDATE: time period of interest;
% VARNAME: land use classification number variable name; value 1-16, 254,255 
% OPTION: output format: by day or by year;
% DATA_OPTION: not used here?
% EXTRAOPTION: 1 --> only create weight files; 0--> interpolate the data; 2 --> delete all files; 3: create maps while doing interpolation
% FORMAT: 'mat': output file is matlab file; 'h5': output file is hdf file

%% return values:
% ReturnStatus==0: still in processing
% ReturnStatus==1: processing done before; 
% ReturnStatus==2, some raw data not available and needs to do it again;
% ReturnStatus==3, some input data not available has to quit;

%% example
% Interpolate_MCD12Q1_Main('AQRVPM25',[datenum(2012,1,1),datenum(2012,12,31)],'GMTED2010bln300','By-Year','',0)

%% code
function ReturnStatus = Interpolate_MCD12Q1_Main(SYSTEM,SITENAME,INPUTDATE,VARNAME,OPTION,DATA_OPTION,EXTRAOPTION,FORMAT)

TypeList = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,254];
% distance threshold used to extract values from
THRESHOLD = 5*1000;
% maximal number of grid cells used to extract values from
N_Neighbour = floor((THRESHOLD)^2*pi/203920);

INPUTYEAR = year(INPUTDATE(1));
if(strcmp(OPTION,'By-Date'))   
    DATE_TO_PROCESS = (INPUTDATE(1):1:INPUTDATE(2)) - datenum(INPUTYEAR,1,1)+1;
elseif(strcmp(OPTION,'By-Day')||strcmp(OPTION,'By-Year'))   
    DATE_TO_PROCESS = 1:1:yeardays(INPUTYEAR);
end

%% path
DIRPATH = '../data/processed_data/LANDUSE_aggregate/';
OUTPUTPATH_ROOT = '../../processed_data/';
Sep = '/';

% file path to the shape file used for making maps
EnvironPara = struct();
EnvironPara.states = geoshape(shaperead('usastatehi', 'UseGeoCoords', true));
EnvironPara.S_back = geoshape(shaperead('../data/shapefile/US_WGS_clipcoast_2.shp', 'UseGeoCoords', true));
EnvironPara.US_North_America_Equidistant_Conic = '../data/shapefile/US_North_America_Equidistant_Conic.shp';
EnvironPara.US = '../data/shapefile/US.shp';
GCS = 'North_America_Equidistant_Conic';%% projection system to be used

FileNameTemplate = 'MCD12Q1Post_$TYPE$_USMCD12Q1_$STARTDATE$_$ENDDATE$.mat';
OutputTemplate = ['MCD12Q1_$TYPE$_$OPTION$_$SITENAME$_$STARTDATE$_$ENDDATE$.',FORMAT];

if(INPUTYEAR>2013)
    FileNameTemplate = strrep(strrep(FileNameTemplate,'$STARTDATE$',datestr(datenum(2013,1,1),'yyyymmdd')),'$ENDDATE$',datestr(datenum(2013,12,31),'yyyymmdd'));    
elseif(INPUTYEAR<2001)
    FileNameTemplate = strrep(strrep(FileNameTemplate,'$STARTDATE$',datestr(datenum(2001,1,1),'yyyymmdd')),'$ENDDATE$',datestr(datenum(2001,12,31),'yyyymmdd'));
else
    FileNameTemplate = strrep(strrep(FileNameTemplate,'$STARTDATE$',datestr(datenum(INPUTYEAR,1,1),'yyyymmdd')),'$ENDDATE$',datestr(datenum(INPUTYEAR,12,31),'yyyymmdd'));
end
OutputTemplate = strrep(strrep(strrep(strrep(OutputTemplate,'$STARTDATE$',datestr(datenum(INPUTYEAR,1,1),'yyyymmdd')),'$ENDDATE$',datestr(datenum(INPUTYEAR,12,31),'yyyymmdd')),'$OPTION$',[DATA_OPTION,'_Thres',num2str(THRESHOLD)]),'$SITENAME$',SITENAME);


%% output
OUTPUTPATH = [OUTPUTPATH_ROOT,SITENAME,Sep,VARNAME,Sep];
EnvironPara.OUTPUTPATH_PIC = OUTPUTPATH;
mkdir(OUTPUTPATH);
%% check whether has been processed before
TempFileName = arrayfun(@(i)  ([OUTPUTPATH,strrep(OutputTemplate,'$TYPE$',['Type',num2str(TypeList(i))])]),1:1:length(TypeList),'UniformOutput',false);
if(2==EXTRAOPTION)
    for i=1:length(TempFileName)
        delete(TempFileName{i});
        fprintf('delete ...%s\n',TempFileName{i});
    end
    ReturnStatus = 1;
    return;
else
    if(cellfun(@(x) exist(x,'file'),TempFileName))%% all have been processed
        fprintf('SKIP! ...%s..%s\n',TempFileName{1},TempFileName{length(TypeList)});
        ReturnStatus = 1;
        return;
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

%% calculate weight
mkdir([OUTPUTPATH_ROOT,SITENAME,Sep,'Temp',Sep]);
WeightFile = [OUTPUTPATH_ROOT,SITENAME,Sep,'Temp',Sep,'NeighbourWeightMatrix_',[DATA_OPTION,'_Thres',num2str(THRESHOLD)],'_','US',VARNAME,'_',SITENAME,'.mat'];
if(exist(WeightFile,'file'))
    Result = LoadData_function(WeightFile);
    if(strcmp(DATA_OPTION,'Mean')||strcmp(DATA_OPTION,'Peak1')||strcmp(DATA_OPTION,'Peak2')||strcmp(DATA_OPTION,'Nearest')||strcmp(DATA_OPTION,'Nearest4'))
        Sum_Weight_matrix = Result.Sum_Weight_matrix;
        Weight_matrix = Result.Weight_matrix;
    elseif(strcmp(DATA_OPTION,'Max')||strcmp(DATA_OPTION,'Min')||strcmp(DATA_OPTION,'Diff')||strcmp(DATA_OPTION,'Var'))
        Index_d = Result.Index_d;
        n = Result.n;
    end
else
    %% raw data location list
    SiteData= LoadData_function([DIRPATH,'US',VARNAME,'Site_',GCS,'.mat'],'SiteData');
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
    CurrentFile1 = [OUTPUTPATH_ROOT,SITENAME,Sep,'Temp',Sep,'NeighbourAdjacentList_',[DATA_OPTION,'_Thres',num2str(THRESHOLD)],'_','US',VARNAME,'_',SITENAME,'.mat'];
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
    elseif(strcmp(DATA_OPTION,'Nearest4'))
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

%% read data
for i = 1:length(TypeList)
    CurrentType = ['Type',num2str(TypeList(i))];
    TempInputFile = [DIRPATH,strrep(FileNameTemplate,'$TYPE$',CurrentType)];
    TempOutput = TempFileName{i};
    
    if(exist(TempOutput,'file'))
        fprintf('skipping...%s\n',TempOutput);
        continue;
    end
    
    try
        fprintf('reading...%s\n',TempInputFile);
        Data = LoadData_function(TempInputFile);
        Data = Data.Result;
    catch exception
        fprintf('error reading...%s %s\n',TempInputFile,getReport(exception));
        ReturnStatus = 2;
        return;
    end
    
    if(length(Data)>1)
       Data = Data'; 
    end
    
    if(strcmp(DATA_OPTION,'Mean') || strcmp(DATA_OPTION,'Peak1') || strcmp(DATA_OPTION,'Peak2') ||strcmp(DATA_OPTION,'Nearest') ||strcmp(DATA_OPTION,'Nearest4'))
        Result = (Data*Weight_matrix)./Sum_Weight_matrix;
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

ReturnStatus = 0;


