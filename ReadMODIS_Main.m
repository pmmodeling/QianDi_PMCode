%% main function to read MODIS data;
% this function reads MOD13A2, MOD09A1 extract variable field of interests; 
% aggregate data from different tile together, match the data to consistent grid
% cells, make maps, and save results.

%% version history

%% parameter:
% YEAR: which year of data to be processed?
% FILETYPE = 10; % different MODIS data set: 9: MOD13A2 or 10: MOD09A1 or
% 12: MCD12Q1 (not used anymore)

%% example: 
% ReadMODIS_Main(2015 ,9)

%% code
function ReadMODIS_Main(YEAR,FILETYPE)

% different MODIS data set
if(9 == FILETYPE)
    % NDVI
    OPTION = 'MOD13A2';
    TEMPLATE = 'MOD13A2.A';
elseif(10 == FILETYPE)
    % surface reflectance
    OPTION = 'MOD09A1';
    TEMPLATE = 'MOD09A1.A';
elseif(12 == FILETYPE)
    % satellite-based landuse classification
    OPTION = 'MCD12Q1';
    TEMPLATE = 'MCD12Q1.A';    
end

% path
DIRPATHROOT = '../data/unprocessed/';
SEP = '/';

DIRPATH = [DIRPATHROOT,OPTION,SEP];
OUTPUTPATH = ['../data/aggregate/',OPTION,'/'];

% file path of shape file
EnvironPara = struct();
EnvironPara.states = geoshape(shaperead('usastatehi', 'UseGeoCoords', true));
EnvironPara.S_back = geoshape(shaperead('../data/shapefile/US_WGS_clipcoast_2.shp', 'UseGeoCoords', true));
EnvironPara.US_North_America_Equidistant_Conic = '../data/shapefile/US_North_America_Equidistant_Conic.shp';
EnvironPara.US = '../data/shapefile/US.shp';
GCS = 'North_America_Equidistant_Conic';

mkdir(OUTPUTPATH);
diary([OUTPUTPATH,'ReadMODIS_Main_',OPTION,num2str(YEAR),'.txt']);
EnvironPara.OUTPUTPATH_PIC = OUTPUTPATH;
fprintf('OPTION:%s\n',OPTION);
fprintf('TEMPLATE:%s\n',TEMPLATE);

StartDay = datenum([YEAR,1,1]);

% find out all raw data files
list = dir(DIRPATH);
fnames = {list.name}.';
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(fnames, str, 'all'));

%% load grid cell locations --- this is the unified points
try
    load([OUTPUTPATH,'US',OPTION,'Site_',GCS,'.mat']);% "USMOD13A2Site_OLD is lat/lon"
catch exception
   fprintf('error %s\n',exception.message);
   fprintf('check this file is correct and uploaded: %s\n',[OUTPUTPATH,'US',OPTION,'Site_',GCS,'.mat'])
   return;
end
AllSiteData = SiteData;
N_Site = size(AllSiteData,1);
N_Day = yeardays(YEAR);
clear SiteData

%% read each individual files  and extract variables of interests
for i=1:N_Day
    
    % find out all files with file names matching a certain pattern
    CurrentDay = StartDay+i-1;
    fprintf('%s\n',datestr(CurrentDay,'yyyy-mm-dd')); 
    FileName_Template = [TEMPLATE,num2str(year(CurrentDay),'%02d'),num2str(i,'%03d'),'.h[\d]{2}v[\d]{2}.*'];  
    
    FileName = fnames(FIND(FileName_Template));
    fprintf('%d\n',length(FileName));  
    
    TempOutputFileName = [OUTPUTPATH,OPTION,'_',datestr(CurrentDay,'yyyymmdd'),'_',datestr(CurrentDay,'yyyymmdd'),'.mat'];
    Result = nan(N_Site,1);
    % for each matched file
    if(~isempty(FileName) && ~exist(TempOutputFileName,'file'))
        try
            load(TempOutputFileName);
            fprintf('this day has been processed!!! skip!!!..%s\n',datestr(CurrentDay,'yyyy-mm-dd')); 
            continue;
        catch
            fprintf('processing..%s\n',datestr(CurrentDay,'yyyy-mm-dd')); 
        end
        
        % for all files on that day
        for j=1:length(FileName)
            fprintf('\t%s\n',FileName{j});

            % extract variable of interest from MODIS hdf files; variable
            % of interests are specified inside the function
            [data,CurrentSiteData,date] = ReadMODIS_function1([DIRPATH,FileName{j}],OPTION);

            % match the to the unified aggregated grid points
            [Lia,Locb] = ismember(int64(AllSiteData.SiteCode),int64(CurrentSiteData.SiteCode));
            Locb_1 = Locb(Lia);
            Result(Lia) = data(Locb_1);

            % label abnormal values as missingness           
            if(9 == FILETYPE)
                Result(Result<=-3000)=nan;
            elseif(10 == FILETYPE)
                Result(Result<=-28672)=nan;  
            elseif(12 == FILETYPE)
                Result(Result==255)=nan;
            end

            % descriptive analysis of the extract value
            fprintf('\tcount and percentage of nan data:%d\t%d\n',sum(~isnan(data(Locb_1))),sum(isnan(Result))/length(Result));
            
        end
        
        % save results
        fprintf('\tpercentage of nan:%f\tmax:%f\tmean:%f\tmin:%f\n',sum(isnan(Result))/length(Result),...
        nanmax(Result),nanmean(Result),nanmin(Result));
        save(TempOutputFileName,'Result');  
         Visualization_USResult_1([OPTION,'_',datestr(CurrentDay,'yyyymmdd'),'_',datestr(CurrentDay,'yyyymmdd')],Result,AllSiteData,EnvironPara);
        
    end
end

diary off;
fclose all;
