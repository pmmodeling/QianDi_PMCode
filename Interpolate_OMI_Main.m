%%	Extract values from OMI data set; it can also extract values from MOD04L2, forest fire, CMAQ output
% used for OMI data set, MOD04L2, forest fire, CMAQ output 

%% version history
%% 2016-11- 12: change INPUTYEAR to INPUTDATE

%% parameters
% SYSTEM: on what system this code is to run; some macro definitions;
% SITENAME: name of the grid cells to be interpolated at;
% INPUTDATE: time period of interest;
% VARNAME: variable name of OMI data set, 'OMAERUVd','OMAEROe',...'MOD04L2';
% OPTION: output format: by day or by year;
% DATA_OPTION: what kind of output min, max, mean??
% EXTRAOPTION: 1 --> only create weight files; 0--> interpolate the data; 2 --> delete all files; 3: create maps while doing interpolation
% FORMAT: 'mat': output file is matlab file; 'h5': output file is hdf file

%% return values:
% ReturnStatus==0: still in processing
% ReturnStatus==1: processing done before; 
% ReturnStatus==2, some raw data not available and needs to do it again;
% ReturnStatus==3, some input data not available has to quit;

%% example:
% Interpolate_OMI_Main('EPANO2',[datenum(2004,1,1),datenum(2004,12,31)],'GEOSChem_BC','By-Year','Mean',0)

%% code
function ReturnStatus = Interpolate_OMI_Main(SITENAME,INPUTDATE,VARNAME,OPTION,DATA_OPTION,EXTRAOPTION,FORMAT)

INPUTYEAR = year(INPUTDATE(1));
if(strcmp(OPTION,'By-Date'))
    DATE_TO_PROCESS = (INPUTDATE(1):1:INPUTDATE(2)) - datenum(INPUTYEAR,1,1)+1;%index: day of year
elseif(strcmp(OPTION,'By-Day') || strcmp(OPTION,'By-Date'))
    DATE_TO_PROCESS = 1:1:yeardays(INPUTYEAR);
end

% specify the variable of interest for each data set
if(any(strfind(VARNAME,'OMAERUVd')))
    % no data before 2005, return!
    if(INPUTYEAR<2005)
        ReturnStatus = 1;
        return;
    end
    AOD_OPTION = 'OMAERUVd';
    THRESHOLD = 200*1000;
    N_Neighbour = 4;
elseif(any(strfind(VARNAME,'OMAEROe')))
    % no data before 2005, return!
    if(INPUTYEAR<2005)
        ReturnStatus = 1;
        return;
    end
    AOD_OPTION = 'OMAEROe';
    THRESHOLD = 100*1000;
    N_Neighbour = 4;
elseif(any(strfind(VARNAME,'OMUVBd'))) 
    % no data before 2005, return!
    if(INPUTYEAR<2005)
        ReturnStatus = 1;
        return;
    end
    AOD_OPTION = 'OMUVBd';
    THRESHOLD = 100*1000;
    N_Neighbour = 4;
elseif(any(strfind(VARNAME,'OMO3PR')))
    % no data before 2005, return!
    if(INPUTYEAR<2005)
        ReturnStatus = 1;
        return;
    end
    AOD_OPTION = 'OMO3PR';
    THRESHOLD = 100*1000;
    N_Neighbour = 4;
elseif(any(strfind(VARNAME,'OMNO2d')))
    % no data before 2005, return!
    if(INPUTYEAR<2005)
        ReturnStatus = 1;
        return;
    end
    AOD_OPTION = 'OMNO2d';
    THRESHOLD = 100*1000;
    N_Neighbour = 4;
elseif(any(strfind(VARNAME,'OMSO2e')))
    % no data before 2005, return!
    if(INPUTYEAR<2005)
        ReturnStatus = 1;
        return;
    end
    AOD_OPTION = 'OMSO2e';
    THRESHOLD = 100*1000;
    N_Neighbour = 4;
elseif(any(strfind(VARNAME,'OMTO3e')))
    % no data before 2005, return!
    if(INPUTYEAR<2005)
        ReturnStatus = 1;
        return;
    end
    AOD_OPTION = 'OMTO3e';
    THRESHOLD = 100*1000;
    N_Neighbour = 4;
elseif(any(strfind(VARNAME,'MOD04L2')))
    AOD_OPTION = 'MOD04L2';
    THRESHOLD = 100*1000;
    N_Neighbour = 4;   
elseif(any(strfind(VARNAME,'GFEDFireCarbon')))
    AOD_OPTION = 'GFED';
    THRESHOLD = 100*1000;
    N_Neighbour = 4;   
elseif(any(strfind(VARNAME,'CAMS_NO2')))
    % no data before 2003, return!
    if(INPUTYEAR<2003)
        ReturnStatus = 1;
        return;
    end
    AOD_OPTION = 'CAMS';
    THRESHOLD = 100*1000;
    N_Neighbour = 4;  
elseif(any(strfind(VARNAME,'CMAQ_NO2'))||any(strfind(VARNAME,'CMAQ_NO2_Vertical'))||any(strfind(VARNAME,'CMAQ_NOX'))||...
       any(strfind(VARNAME,'CMAQ_Ozone'))||any(strfind(VARNAME,'CMAQ_Ozone_Vertical'))||any(strfind(VARNAME,'CMAQ_PM25_EC'))||...
       any(strfind(VARNAME,'CMAQ_PM25_NH4'))||any(strfind(VARNAME,'CMAQ_PM25_NO3'))||any(strfind(VARNAME,'CMAQ_PM25_OC'))||...
       any(strfind(VARNAME,'CMAQ_PM25_OM'))||any(strfind(VARNAME,'CMAQ_PM25_SO4'))||any(strfind(VARNAME,'CMAQ_PM25_TOT'))||...
       any(strfind(VARNAME,'CMAQ_RH'))||any(strfind(VARNAME,'CMAQ_TA'))||any(strfind(VARNAME,'CMAQ_AIR_DENS'))||any(strfind(VARNAME,'CMAQ_PM25_Vertical')))
    % no data after 2014, return!
   if(INPUTYEAR>2014)
        ReturnStatus = 1;
        return;
    end
    AOD_OPTION = 'CMAQ';% the folder name
    THRESHOLD = 50*1000;
    N_Neighbour = 4; 
elseif(any(strfind(VARNAME,'MERRA2aer_SO4')) || any(strfind(VARNAME,'MERRA2aer_OCPHOBIC')) || any(strfind(VARNAME,'MERRA2aer_OCPHILIC')) || any(strfind(VARNAME,'MERRA2aer_BCPHOBIC')) || any(strfind(VARNAME,'MERRA2aer_BCPHILIC')))
    AOD_OPTION = 'MERRA2aer';
    THRESHOLD = 100*1000;
    N_Neighbour = 4;
elseif(any(strfind(VARNAME,'GEOSChem_BC')) || any(strfind(VARNAME,'GEOSChem_NH4')) || any(strfind(VARNAME,'GEOSChem_NO2')) || any(strfind(VARNAME,'GEOSChem_OA')) || any(strfind(VARNAME,'GEOSChem_PM25'))|| any(strfind(VARNAME,'GEOSChem_SO4'))|| any(strfind(VARNAME,'GEOSChem_O3')) || any(strfind(VARNAME,'GEOSChem_NO2Scaling')) || any(strfind(VARNAME,'GEOSChem_O3Scaling')) || any(strfind(VARNAME,'GEOSChem_PM25Scaling')))
    AOD_OPTION = 'GEOSChem';
    THRESHOLD = 100*1000;
    N_Neighbour = 4;   
end

DIRPATH = ['../data/aggregate/',AOD_OPTION,'/'];
OUTPUTPATH_ROOT = '../../processed_data/';
Sep = '/';

% file path to shape file used for making maps
EnvironPara = struct();
EnvironPara.states = geoshape(shaperead('usastatehi', 'UseGeoCoords', true));
EnvironPara.S_back = geoshape(shaperead('../data/shapefile/US_WGS_clipcoast_2.shp', 'UseGeoCoords', true));
EnvironPara.US_North_America_Equidistant_Conic = '../data/shapefile/US_North_America_Equidistant_Conic.shp';
EnvironPara.US = '../data/shapefile/US.shp';
GCS = 'North_America_Equidistant_Conic';%% projection system to be used

FileNameTemplate = '$VARNAME$_$SITENAME$_$STARTDATE$_$ENDDATE$.mat';
OutputTemplate = ['$VARNAME$_$OPTION$_$SITENAME$_$STARTDATE$_$ENDDATE$.',FORMAT];

%% specify the file name of aggregate files
if(any(strfind(VARNAME,'OMAERUVd'))||any(strfind(VARNAME,'OMAEROe'))||any(strfind(VARNAME,'OMUVBd'))||...
        any(strfind(VARNAME,'OMO3PR'))||any(strfind(VARNAME,'OMNO2d'))||any(strfind(VARNAME,'OMSO2e'))||...
        any(strfind(VARNAME,'OMTO3e'))||any(strfind(VARNAME,'MOD04L2'))||any(strfind(VARNAME,'GFEDFireCarbon'))||...
        any(strfind(VARNAME,'MERRA2aer_SO4')) || any(strfind(VARNAME,'MERRA2aer_OCPHOBIC')) || any(strfind(VARNAME,'MERRA2aer_OCPHILIC')) || any(strfind(VARNAME,'MERRA2aer_BCPHOBIC')) || any(strfind(VARNAME,'MERRA2aer_BCPHILIC'))||...
        any(strfind(VARNAME,'CAMS_NO2')) || ...
        any(strfind(VARNAME,'GEOSChem_BC')) || any(strfind(VARNAME,'GEOSChem_NH4')) || any(strfind(VARNAME,'GEOSChem_NO2')) || any(strfind(VARNAME,'GEOSChem_OA')) || any(strfind(VARNAME,'GEOSChem_O3')) || any(strfind(VARNAME,'GEOSChem_PM25'))|| any(strfind(VARNAME,'GEOSChem_SO4')) || any(strfind(VARNAME,'GEOSChem_NO2Scaling')) || any(strfind(VARNAME,'GEOSChem_O3Scaling')) || any(strfind(VARNAME,'GEOSChem_PM25Scaling')) )
    
    FileNameTemplate = strrep(strrep(FileNameTemplate,'$SITENAME$',['US',AOD_OPTION]),'$VARNAME$',VARNAME);
    TempInputFile = [DIRPATH,strrep(strrep(FileNameTemplate,'$STARTDATE$',num2str(INPUTYEAR)),'$ENDDATE$',num2str(INPUTYEAR))];
  
% CMAQ outptus
elseif(any(strfind(VARNAME,'CMAQ_NO2'))||any(strfind(VARNAME,'CMAQ_NO2_Vertical'))||any(strfind(VARNAME,'CMAQ_NOX'))||...
       any(strfind(VARNAME,'CMAQ_Ozone'))||any(strfind(VARNAME,'CMAQ_Ozone_Vertical'))||any(strfind(VARNAME,'CMAQ_PM25_EC'))||...
       any(strfind(VARNAME,'CMAQ_PM25_NH4'))||any(strfind(VARNAME,'CMAQ_PM25_NO3'))||any(strfind(VARNAME,'CMAQ_PM25_OC'))||...
       any(strfind(VARNAME,'CMAQ_PM25_OM'))||any(strfind(VARNAME,'CMAQ_PM25_SO4'))||any(strfind(VARNAME,'CMAQ_PM25_TOT'))||...
       any(strfind(VARNAME,'CMAQ_RH'))||any(strfind(VARNAME,'CMAQ_TA'))||any(strfind(VARNAME,'CMAQ_AIR_DENS'))||any(strfind(VARNAME,'CMAQ_PM25_Vertical')))
    FileNameTemplate = strrep(strrep(FileNameTemplate,'$SITENAME$','CMAQ'),'$VARNAME$',VARNAME);
    TempInputFile = [DIRPATH,strrep(strrep(FileNameTemplate,'$STARTDATE$',num2str(INPUTYEAR)),'$ENDDATE$',num2str(INPUTYEAR))];
    
    % for 2002 to 2006, although 12 km CMAQ is available, use the 36 km one,s since it covers the whole country
    if(INPUTYEAR == 2002)
        TempInputFile = [DIRPATH,strrep(strrep(FileNameTemplate,'$STARTDATE$',[num2str(INPUTYEAR),'B']),'$ENDDATE$',[num2str(INPUTYEAR),'B'])];
    elseif(INPUTYEAR == 2003)
        TempInputFile = [DIRPATH,strrep(strrep(FileNameTemplate,'$STARTDATE$',[num2str(INPUTYEAR),'B']),'$ENDDATE$',[num2str(INPUTYEAR),'B'])];
    elseif(INPUTYEAR == 2004)
        TempInputFile = [DIRPATH,strrep(strrep(FileNameTemplate,'$STARTDATE$',[num2str(INPUTYEAR),'A']),'$ENDDATE$',[num2str(INPUTYEAR),'A'])];
    elseif(INPUTYEAR == 2005)
        TempInputFile = [DIRPATH,strrep(strrep(FileNameTemplate,'$STARTDATE$',[num2str(INPUTYEAR),'B']),'$ENDDATE$',[num2str(INPUTYEAR),'B'])];
    elseif(INPUTYEAR == 2006)
        TempInputFile = [DIRPATH,strrep(strrep(FileNameTemplate,'$STARTDATE$',[num2str(INPUTYEAR),'B']),'$ENDDATE$',[num2str(INPUTYEAR),'B'])];
    end     
end
OutputTemplate = strrep(strrep(strrep(OutputTemplate,'$OPTION$',[DATA_OPTION,'_Thres',num2str(THRESHOLD)]),'$SITENAME$',SITENAME),'$VARNAME$',VARNAME);


%% output
OUTPUTPATH = [OUTPUTPATH_ROOT,SITENAME,Sep,AOD_OPTION,Sep];
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
WeightFile = [OUTPUTPATH_ROOT,SITENAME,Sep,'Temp',Sep,'NeighbourWeightMatrix_',[DATA_OPTION,'_Thres',num2str(THRESHOLD)],'_','US',AOD_OPTION,'_',SITENAME,'.mat'];
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
    SiteData= LoadData_function([DIRPATH,'US',AOD_OPTION,'Site_',GCS,'.mat'],'SiteData');
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
    CurrentFile1 = [OUTPUTPATH_ROOT,SITENAME,Sep,'Temp',Sep,'NeighbourAdjacentList_',[DATA_OPTION,'_Thres',num2str(THRESHOLD)],'_','US',AOD_OPTION,'_',SITENAME,'.mat'];
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
    if(strcmp(DATA_OPTION,'Mean')||strcmp(DATA_OPTION,'AnnualAverage'))
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

%% Load data
if(strcmp(OPTION,'By-Year'))
    if(exist(TempFileName,'file'))
        fprintf('skipping...%s\n',TempFileName);
    else
        % read data
        try
            fprintf('reading...%s\n',TempInputFile);            
            Data = LoadData_function(TempInputFile);
            Name = fieldnames(Data);
            Data = getfield(Data,Name{1});
        catch exception
            fprintf('error reading...%s %s\n',TempInputFile,getReport(exception));
            ReturnStatus = 2;
            return;
        end
        
        % process data, and calcualte mean, maximal, minimal, first-order
        % difference and so on.
        Result = nan(yeardays(INPUTYEAR),N_Traget);
        if(strcmp(DATA_OPTION,'Mean') || strcmp(DATA_OPTION,'Peak1') || strcmp(DATA_OPTION,'Peak2') ||strcmp(DATA_OPTION,'Nearest')||strcmp(DATA_OPTION,'Nearest4')||strcmp(DATA_OPTION,'AnnualAverage'))
            for i = 1:yeardays(INPUTYEAR)
                Result(i,:) = MultipleWeightMatrix_1(Weight_matrix',Data(i,:)')';
%                 Result(i,:) = (Data(i,:)*Weight_matrix)./Sum_Weight_matrix;
            end
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
        fprintf('saving...%s\n',TempFileName);
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
    
elseif(strcmp(OPTION,'By-Day') || strcmp(OPTION,'By-Date'))
    
    % read data
    try
        fprintf('reading...%s\n',TempInputFile);
        DataAll = LoadData_function(TempInputFile);
        Name = fieldnames(DataAll);
        DataAll = getfield(DataAll,Name{1});
    catch exception
        fprintf('error reading...%s %s\n',TempInputFile,getReport(exception));
        ReturnStatus = 2;
        return;
    end
    
    % process data
    for i = DATE_TO_PROCESS-DATE_TO_PROCESS(1)+1% make sure the i stands for day of year or day of TempFileName
        Data = DataAll(DATE_TO_PROCESS(i),:);
        if(exist(TempFileName{i},'file'))
            fprintf('skipping...%s\n',TempFileName{i});
        else
            
            % process data, and calcualte mean, maximal, minimal, first-order
            % difference and so on.
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
            fprintf('saving...%s\n',TempFileName{i});
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
end


ReturnStatus = 0;


