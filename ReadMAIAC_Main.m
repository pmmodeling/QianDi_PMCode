%% read MAIAC AOD Data from different tiles and aggregate them into unified grids
% store the measurements from two satellites by date

%% version history
% adapted from ReadMOD11A1_Main.m
% 2016-11-24: read MAIAC data and arrange them by date
% 2017-05-22: OPTION:Aqua or TERRA

%% parameters:
% YEAR = 2000;% year of files to be read
% OPTION = 'Aqua'% Aqua or TERRA

%% example 
%ReadMAIAC_Main(2010,'Aqua')
function ReadMAIAC_Main(YEAR,OPTION)

% path
DIRPATHROOT = '../data/unprocessed/';
OUTPUTPATH = '../data/aggregate/MAIACUS/';
SEP = '/';

% path of the shape file, used for making maps 
EnvironPara = struct();
EnvironPara.states = geoshape(shaperead('usastatehi', 'UseGeoCoords', true));
EnvironPara.S_back = geoshape(shaperead('../data/shapefile/US_WGS_clipcoast_2.shp', 'UseGeoCoords', true));
EnvironPara.US_North_America_Equidistant_Conic = '../data/shapefile/US_North_America_Equidistant_Conic.shp';
EnvironPara.US = '../data/shapefile/US.shp';
GCS = 'North_America_Equidistant_Conic';


EnvironPara.OUTPUTPATH_PIC = OUTPUTPATH;

%% other parameters
DateA = datenum(YEAR,1,1);
DateB = datenum(YEAR,12,31);
TileName = {'h01v03','h01v04','h01v05','h02v03','h02v04','h02v05','h03v03','h03v04','h03v05','h04v03','h04v04','h04v05','h05v03'};

%% obtain site data
VARIABLE_LIST = {'Optical_Depth_047','Optical_Depth_055','cosVZA'};
SITENAME = 'USMAIACUS1km'; 
if(exist([OUTPUTPATH,SITENAME,'Site','.mat'],'file'))
    SiteData_1km = load([OUTPUTPATH,SITENAME,'Site','.mat']);%SiteData, TileName --- tile name is order of tiles while creating the SiteData
    SiteData_1km = SiteData_1km.SiteData;
else
    disp('no site data!!!');
    return;
end

SITENAME = 'USMAIACUS5km'; 
if(exist([OUTPUTPATH,SITENAME,'Site','.mat'],'file'))
    SiteData_5km = load([OUTPUTPATH,SITENAME,'Site','.mat']);%SiteData, TileName --- tile name is order of tiles while creating the SiteData
    SiteData_5km = SiteData_5km.SiteData;
else
    disp('no site data!!!');
    return;
end

%% get all list of all files;
FileList = [];
for l = 1:length(TileName)
    TempFilelist = dir([DIRPATHROOT,'dataportal.nccs.nasa.gov',SEP,'DataRelease',SEP,'NorthAmerica_2000-2016',SEP,TileName{l},SEP,num2str(YEAR),SEP]);
    FileList = cat(1,FileList,TempFilelist);
end
fnames = {FileList.name}.';
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(fnames, str, 'all'));

%% begin to process all files in that year
diary([OUTPUTPATH,'ReadMAIAC_Main',num2str(YEAR),'.txt']);
for i = DateA:DateB
    CurrentDay = i;
    
    fprintf('%s\n',datestr(CurrentDay));
    
    % get julian date and find all files with matching file name pattern
    JulianDay = CurrentDay - datenum(year(CurrentDay),1,1)+1;
    if(strcmp(OPTION,'Aqua'))
        FileName_Template = ['MAIACAAOT.h[\d]{2}v[\d]{2}.',num2str(year(CurrentDay)), sprintf('%03d',JulianDay),'[\d]{4}.hdf'];
    elseif(strcmp(OPTION,'Terra'))
        FileName_Template = ['MAIACTAOT.h[\d]{2}v[\d]{2}.',num2str(year(CurrentDay)), sprintf('%03d',JulianDay),'[\d]{4}.hdf'];
    end
    FileName = fnames(FIND(FileName_Template));
    
    % to store output
    Data_Optical_Depth_047 = nan(size(SiteData_1km,1),1);
    Data_Optical_Depth_055 = nan(size(SiteData_1km,1),1);
    Data_cosVZA = nan(size(SiteData_5km,1),1);
    Data_AOT_QA = nan(size(SiteData_1km,1),1);
    Data_Optical_Depth_047_count = zeros(size(SiteData_1km,1),1);
    Data_Optical_Depth_055_count = zeros(size(SiteData_1km,1),1);
    Data_cosVZA_count = zeros(size(SiteData_5km,1),1);
    Data_AOT_QA_count = zeros(size(SiteData_1km,1),1);
    
    TempOutputFileName = [OUTPUTPATH,'MAIACUS',OPTION,'_',SITENAME,'_',datestr(CurrentDay,'yyyymmdd'),'_',datestr(CurrentDay,'yyyymmdd'),'.mat'];
    
    % to store data
    if(~isempty(FileName))       
        if(exist(TempOutputFileName,'file'))
            fprintf('this day has been processed!!! skip!!!..%s\n',datestr(CurrentDay,'yyyy-mm-dd'));
            continue;
        else
            fprintf('processing..%s\n',datestr(CurrentDay,'yyyy-mm-dd'));
        end
        
        % if this file is for that day
        for j=1:length(FileName)
            TempFileName = FileName{j};
            
            Index = regexp(TempFileName,'h[\d]{2}');
            h = TempFileName(Index(1)+1:Index(1)+2);
            Index = regexp(TempFileName,'v[\d]{2}');
            v = TempFileName(Index(1)+1:Index(1)+2);
            
            % create unique id for each grid cell
            TempSiteCode_1km = int64(repmat(str2double(['2',h,v])*10000000,[1200*1200,1])+(1:1:1200*1200)');
            TempSiteCode_5km = int64(repmat(str2double(['2',h,v])*10000000,[240*240,1])+(1:1:240*240)');
            
            % match the grid cells from each tile to the national
            % consistent grid cells
            [lia_1km,locb_1km]=ismember(SiteData_1km.SiteCode,TempSiteCode_1km);
            [lia_5km,locb_5km]=ismember(SiteData_5km.SiteCode,TempSiteCode_5km);
            locb_1km = locb_1km(lia_1km);
            locb_5km = locb_5km(lia_5km);
            
            fprintf('\t%s\n',TempFileName);
            FILE_NAME = [DIRPATHROOT,'dataportal.nccs.nasa.gov',SEP,'DataRelease',SEP,'NorthAmerica_2000-2016',SEP,['h',h,'v',v],SEP,num2str(YEAR),SEP,TempFileName];
            try
                % read the raw data
                Temp_Data_Optical_Depth_047 = double(hdfread(FILE_NAME,'Optical_Depth_047'));
                Temp_Data_Optical_Depth_055 = double(hdfread(FILE_NAME,'Optical_Depth_055'));
                Temp_Data_cosVZA = double(hdfread(FILE_NAME,'cosVZA'));
                Temp_Data_AOT_QA = double(hdfread(FILE_NAME,'AOT_QA'));
            catch exception
                fprintf('error!%s\n',exception.message);
                continue;
            end
            
            % reshape
            Temp_Data_Optical_Depth_047 = reshape(Temp_Data_Optical_Depth_047,[size(Temp_Data_Optical_Depth_047,1)*size(Temp_Data_Optical_Depth_047,2),1]);
            Temp_Data_Optical_Depth_055 = reshape(Temp_Data_Optical_Depth_055,[size(Temp_Data_Optical_Depth_055,1)*size(Temp_Data_Optical_Depth_055,2),1]);
            Temp_Data_cosVZA = reshape(Temp_Data_cosVZA,[size(Temp_Data_cosVZA,1)*size(Temp_Data_cosVZA,2),1]);
            Temp_Data_AOT_QA = reshape(Temp_Data_AOT_QA,[size(Temp_Data_AOT_QA,1)*size(Temp_Data_AOT_QA,2),1]);
                     
            % remove invalid data, based on QA flag
            Temp_Data_Optical_Depth_047(Temp_Data_Optical_Depth_047==-28672)=nan;
            Temp_Data_Optical_Depth_047(Temp_Data_Optical_Depth_047>5000|Temp_Data_Optical_Depth_047<-100)=nan;
            Temp_Data_Optical_Depth_047 = Temp_Data_Optical_Depth_047*1.000000000000000e-03;
            
            % adjust the value based on scaling factor
            Temp_Data_Optical_Depth_055(Temp_Data_Optical_Depth_055==-28672)=nan;
            Temp_Data_Optical_Depth_055(Temp_Data_Optical_Depth_055>5000|Temp_Data_Optical_Depth_055<-100)=nan;
            Temp_Data_Optical_Depth_055 = Temp_Data_Optical_Depth_055*1.000000000000000e-03;
            
            % remove invalid data and label it as missing
            Temp_Data_cosVZA(Temp_Data_cosVZA==-28672)=nan;
            Temp_Data_cosVZA(Temp_Data_cosVZA>10000|Temp_Data_cosVZA<0)=nan;
            Temp_Data_cosVZA = Temp_Data_cosVZA*1.000000000000000e-04;
            
            Temp_Data_AOT_QA(Temp_Data_AOT_QA==0)=nan;
            Temp_Data_AOT_QA(Temp_Data_AOT_QA>255|Temp_Data_AOT_QA<0)=nan;
            
            % set count
            Temp_Data_Optical_Depth_047_count = ~isnan(Temp_Data_Optical_Depth_047);
            Temp_Data_Optical_Depth_055_count = ~isnan(Temp_Data_Optical_Depth_055);
            Temp_Data_cosVZA_count = ~isnan(Temp_Data_cosVZA);
            Temp_Data_AOT_QA_count = ~isnan(Temp_Data_AOT_QA);

            % put them into the unified SiteData
            Data_Optical_Depth_047(lia_1km) = nansum([Data_Optical_Depth_047(lia_1km),Temp_Data_Optical_Depth_047(locb_1km)],2);
            Data_Optical_Depth_055(lia_1km) = nansum([Data_Optical_Depth_055(lia_1km),Temp_Data_Optical_Depth_055(locb_1km)],2);
            Data_cosVZA(lia_5km) = nansum([Data_cosVZA(lia_5km),Temp_Data_cosVZA(locb_5km)],2);
            Data_AOT_QA(lia_1km) = nansum([Data_AOT_QA(lia_1km),Temp_Data_AOT_QA(locb_1km)],2);
            
            Data_Optical_Depth_047_count(lia_1km) = nansum([Data_Optical_Depth_047_count(lia_1km),Temp_Data_Optical_Depth_047_count(locb_1km)],2);
            Data_Optical_Depth_055_count(lia_1km) = nansum([Data_Optical_Depth_055_count(lia_1km),Temp_Data_Optical_Depth_055_count(locb_1km)],2);
            Data_cosVZA_count(lia_5km) = nansum([Data_cosVZA_count(lia_5km),Temp_Data_cosVZA_count(locb_5km)],2);
            Data_AOT_QA_count(lia_1km) = nansum([Data_AOT_QA_count(lia_1km),Temp_Data_AOT_QA_count(locb_1km)],2);
        end
        
        % for grid cells with more than one daily measurement, take
        % average;
        Data_Optical_Depth_047 = Data_Optical_Depth_047./Data_Optical_Depth_047_count;
        Data_Optical_Depth_055 = Data_Optical_Depth_055./Data_Optical_Depth_055_count;
        Data_cosVZA = Data_cosVZA./Data_cosVZA_count;
        Data_AOT_QA = Data_AOT_QA./Data_AOT_QA_count;
                
        Data_Optical_Depth_047 = Data_Optical_Depth_047';
        Data_Optical_Depth_055 = Data_Optical_Depth_055';
        Data_cosVZA = Data_cosVZA';
        Data_AOT_QA = Data_AOT_QA';

        % take land mask %% big endian; See onenote QA and Big-endian issue
        Data_AOT_QA(isnan(Data_AOT_QA)) = 0;
        Data_AOT_QA = dec2bin(Data_AOT_QA,16);
        Data_AOT_QA_CloudMask = bin2dec(Data_AOT_QA(:,14:16));% only take 001---Clear; Bits:0-2
        Data_AOT_QA_LandMask = bin2dec(Data_AOT_QA(:,12:13));% only 00 --- Land; Bits:3-4
        Data_AOT_QA_AdjacencyMask = bin2dec(Data_AOT_QA(:,9:11));% only 000 ---  Normal condition; Bits:5-7
        Data_AOT_QA_CloudDetection = bin2dec(Data_AOT_QA(:,5:8));% only 0000 ---  clear; Bits:8-11
        % remove bad values
        Data_AOT_QA_Good = Data_AOT_QA_CloudMask==1&Data_AOT_QA_LandMask==0&Data_AOT_QA_AdjacencyMask==0&Data_AOT_QA_CloudDetection==0;
        Data_Optical_Depth_047(~Data_AOT_QA_Good)=nan;
        Data_Optical_Depth_055(~Data_AOT_QA_Good)=nan;
        
        save(TempOutputFileName,'Data_Optical_Depth_047','Data_Optical_Depth_055','Data_cosVZA','Data_AOT_QA');
        
        if(sum(isnan(Data_Optical_Depth_055))<length(Data_Optical_Depth_055)*0.99)
            % make maps
            Visualization_USResult_1(['MAIACUS',OPTION,'_',SITENAME,'_',datestr(CurrentDay,'yyyymmdd'),'_',datestr(CurrentDay,'yyyymmdd')],Data_Optical_Depth_055,SiteData_1km,EnvironPara);
        end
    end
end
diary off
