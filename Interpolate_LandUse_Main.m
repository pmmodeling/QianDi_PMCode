%% extract value for grid cells of interests from land use raster files

%% version history
% 2016-08-14 skip all file that have been processed before
% 2016-10-02 this code can compute 1-day lag, 5-day lag now
% 2016-10-04: read pressure level data (multiple levels)
% 2016-10-14: remove filling value and missing values

%% parameters
% SITENAME: name of the grid cells to be interpolated at;
% INPUTDATE: time period of interest;
% VARNAME: land use variable name;
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
% Interpolate_LandUse_Main('AQRVPM25',[datenum(2012,1,1),datenum(2012,12,31)],'GMTED2010bln300','By-Year','',0)

%% code
function ReturnStatus = Interpolate_LandUse_Main(SITENAME,INPUTDATE,VARNAME,OPTION,DATA_OPTION,EXTRAOPTION,FORMAT)
 
%% marco
START_DATE = datenum(1998,1,1);
END_DATE= datenum(2020,1,1);
OPTION = 'By-Year';

INPUTYEAR = year(INPUTDATE(1));
if(strcmp(OPTION,'By-Date'))   
    DATE_TO_PROCESS = (INPUTDATE(1):1:INPUTDATE(2)) - datenum(INPUTYEAR,1,1)+1;
elseif(strcmp(OPTION,'By-Day')||strcmp(OPTION,'By-Year'))   
    DATE_TO_PROCESS = 1:1:yeardays(INPUTYEAR);
end

if(any(strfind(VARNAME,'Elevation'))||strcmp(VARNAME,'EPANO2Site_Region')||strcmp(VARNAME,'EPACastNetOzoneSite_Region')||strcmp(VARNAME,'AQRVPM25Site_Region')...
        ||strcmp(VARNAME,'TruckRoute_Traffic1000')||strcmp(VARNAME,'TruckRoute_Traffic100')||strcmp(VARNAME,'TruckRoute_ShortDis1000')||strcmp(VARNAME,'TruckRoute_ShortDis100')||strcmp(VARNAME,'TruckRoute_Density100')||strcmp(VARNAME,'TruckRoute_Density1000'))
    FileNameTemplate = '$VARNAME$.tif';
    OutputTemplate = ['$VARNAME$_$SITENAME$.',FORMAT];
else
    FileNameTemplate = '$VARNAME$_$STARTDATE$_$ENDDATE$.tif';
    OutputTemplate = ['$VARNAME$_$SITENAME$_$STARTDATE$_$ENDDATE$.',FORMAT]; 
end

%% path
DIRPATH = '../data/aggregate/LANDUSE/';
OUTPUTPATH_ROOT = '../../processed_data/';
Sep = '/';

% file path to the shape file used for making maps
EnvironPara = struct();
EnvironPara.states = geoshape(shaperead('usastatehi', 'UseGeoCoords', true));
EnvironPara.S_back = geoshape(shaperead('../data/shapefile/US_WGS_clipcoast_2.shp', 'UseGeoCoords', true));
EnvironPara.US_North_America_Equidistant_Conic = '../data/shapefile/US_North_America_Equidistant_Conic.shp';
EnvironPara.US = '../data/shapefile/US.shp';
GCS = 'North_America_Equidistant_Conic';%% projection system to be used

FileNameTemplate = strrep(strrep(FileNameTemplate,'$VARNAME$',VARNAME),'$SITENAME$',SITENAME);
OutputTemplate = strrep(strrep(OutputTemplate,'$VARNAME$',VARNAME),'$SITENAME$',SITENAME);

%% output
OUTPUTPATH = [OUTPUTPATH_ROOT,SITENAME,Sep,'LandUse',Sep];
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

%%  1 == EXTRAOPTION, we only process weight files
if(1 == EXTRAOPTION)
    ReturnStatus = 1;
    return;
end

%% extract value at grid cells of interests from a raster file
SiteData = LoadData_function([OUTPUTPATH_ROOT,SITENAME,Sep,'Location',Sep,SITENAME,'Site_',GCS,'.mat'],'SiteData');
SiteData_Target = SiteData.SiteData;

FileNameTemplateTemp = strrep(strrep(FileNameTemplate,'$STARTDATE$',num2str(year(INPUTDATE(1)))),'$ENDDATE$',num2str(year(INPUTDATE(1))));
if(exist([DIRPATH,FileNameTemplateTemp],'file'))
    % extract value from a raster file
    [Data, R] = geotiffread([DIRPATH,FileNameTemplateTemp]);
    fprintf('reading...%s\n',[DIRPATH,FileNameTemplateTemp]);
    
    % outside the area is coded as nan
    XIndex = ceil(abs(SiteData_Target.Lon - R.XWorldLimits(1))/R.CellExtentInWorldX);
    YIndex = ceil(abs(SiteData_Target.Lat - R.YWorldLimits(2))/R.CellExtentInWorldY);
    XIndex(XIndex<1)=nan;
    XIndex(XIndex>size(Data,2))=nan;
    YIndex(YIndex<1)=nan;
    YIndex(YIndex>size(Data,1))=nan;
    idx = sub2ind(size(Data), YIndex, XIndex);
    Result = zeros(1,size(SiteData_Target,1));
    Result(~isnan(idx)) = Data(idx(~isnan(idx)));
    Result(Result<-1e30)=nan;
    clear 'Data'
else
    for j=year(INPUTDATE(1)):-1:year(START_DATE)
        FileNameTemplateTemp = strrep(strrep(FileNameTemplate,'$STARTDATE$',num2str(j)),'$ENDDATE$',num2str(j));
        DateA = j;
        if(exist([DIRPATH,FileNameTemplateTemp],'file'))
            % extract value from a raster file
            [Data, R] = geotiffread([DIRPATH,FileNameTemplateTemp]);
            % take value
            XIndex = ceil(abs(SiteData_Target.Lon - R.XWorldLimits(1))/R.CellExtentInWorldX);
            YIndex = ceil(abs(SiteData_Target.Lat - R.YWorldLimits(2))/R.CellExtentInWorldY);
            XIndex(XIndex<1)=nan;
            XIndex(XIndex>size(Data,2))=nan;
            YIndex(YIndex<1)=nan;
            YIndex(YIndex>size(Data,1))=nan;
            idx = sub2ind(size(Data), YIndex, XIndex);
            ResultA = zeros(1,size(SiteData_Target,1));
            ResultA(~isnan(idx)) = Data(idx(~isnan(idx)));
            ResultA(ResultA<-1e30)=nan;
              
            fprintf('\treading A...%s\n',[DIRPATH,FileNameTemplateTemp]);
            clear 'Data'
            break;
        end
    end
    
    for j=year(INPUTDATE(1)):1:year(END_DATE)
        FileNameTemplateTemp = strrep(strrep(FileNameTemplate,'$STARTDATE$',num2str(j)),'$ENDDATE$',num2str(j));
        DateB = j;
        if(exist([DIRPATH,FileNameTemplateTemp],'file'))
            % extract value from a raster file
            [Data, R] = geotiffread([DIRPATH,FileNameTemplateTemp]);
            % take value
            XIndex = ceil(abs(SiteData_Target.Lon - R.XWorldLimits(1))/R.CellExtentInWorldX);
            YIndex = ceil(abs(SiteData_Target.Lat - R.YWorldLimits(2))/R.CellExtentInWorldY);
            XIndex(XIndex<1)=nan;
            XIndex(XIndex>size(Data,2))=nan;
            YIndex(YIndex<1)=nan;
            YIndex(YIndex>size(Data,1))=nan;
            idx = sub2ind(size(Data), YIndex, XIndex);
            ResultB = zeros(1,size(SiteData_Target,1));
            ResultB(~isnan(idx)) = Data(idx(~isnan(idx)));
            ResultB(ResultB<-1e30)=nan;
            
            fprintf('\treading B...%s\n',[DIRPATH,FileNameTemplateTemp]);
            clear 'Data'
            break;
        end
    end
    
    fprintf('temporal interpolating...%s\n',num2str(year(INPUTDATE(1))));
    if(DateB == year(END_DATE))%% find not data afterwards
        Result = ResultA;
    elseif(DateA == year(START_DATE))
        Result = ResultB;
    else
        Result = (ResultA*abs(DateB - year(INPUTDATE(1)))+ResultB*abs(DateA - year(INPUTDATE(1))))/abs(DateA-DateB);
    end 
end

% visualize and save result
if(3 == EXTRAOPTION)
    Visualization_USResult_1(strrep(TempFileName,OUTPUTPATH,''),Result,SiteData_Target,EnvironPara);
end

% save as matlab format or hdf5 format
if(strcmp(FORMAT,'mat'))
    save(TempFileName,'Result');
elseif(strcmp(FORMAT,'h5'))
    hdf5write(TempFileName, 'Result', Result);
end

ReturnStatus = 0;





