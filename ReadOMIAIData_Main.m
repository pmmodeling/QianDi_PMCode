%% read all L3 OMI data 
% this function reads OMAERUVd, OMAEROe, OMNO2d, OMSO2e, OMTO3e, OMUVBd and
% extract variable field of interests; match the data to consistent grid
% cells, make maps, and save results.

%% parameter:
% YEAR = 2015;
% OPTION = 10; % different OMI data set: OMAERUVd, OMAEROe, OMNO2d, OMSO2e, OMTO3e, OMUVBd

%% example: 
% ReadOMIAIData_Main(2010 ,'OMAERUVd')

%% code
function ReadOMIAIData_Main(YEAR,OPTION)

% specify the variable of interest for each data set
if(strcmp(OPTION,'OMAERUVd'))
    TEMPLATE = 'OMI-Aura_L3-OMAERUVd';
    VARIABLE_LIST = {'UVAerosolIndex','FinalAerosolSingleScattAlb500'};
elseif(strcmp(OPTION,'OMUVBd'))
    TEMPLATE = 'OMI-Aura_L3-OMUVBd';
    VARIABLE_LIST = {'UVindex','ViewingZenithAngle','SolarZenithAngle'};
elseif(strcmp(OPTION,'OMAEROe'))
    TEMPLATE = 'OMI-Aura_L3-OMAEROe';
    VARIABLE_LIST = {'VISAerosolIndex','UVAerosolIndex','ViewingZenithAngle','SolarZenithAngle'};
elseif(strcmp(OPTION,'OMTO3e'))
    TEMPLATE = 'OMI-Aura_L3-OMTO3e';
    VARIABLE_LIST = {'ColumnAmountO3','ViewingZenithAngle'};
elseif(strcmp(OPTION,'OMSO2e'))
    TEMPLATE = 'OMI-Aura_L3-OMSO2e';
    VARIABLE_LIST = {'ColumnAmountSO2_PBL'};
elseif(strcmp(OPTION,'OMNO2d'))
    TEMPLATE = 'OMI-Aura_L3-OMNO2d';
    VARIABLE_LIST = {'ColumnAmountNO2StratoCloudScreened'};
end

% path
DIRPATH_ROOT = '../data/unprocessed/';
SEP = '/';

% file path to shape file, used to make maps
EnvironPara = struct();
EnvironPara.states = geoshape(shaperead('usastatehi', 'UseGeoCoords', true));
EnvironPara.S_back = geoshape(shaperead('../data/shapefile/US_WGS_clipcoast_2.shp', 'UseGeoCoords', true));
EnvironPara.US_North_America_Equidistant_Conic = '../data/shapefile/US_North_America_Equidistant_Conic.shp';
EnvironPara.US = '../data/shapefile/US.shp';
GCS = 'North_America_Equidistant_Conic';

DIRPATH = [DIRPATH_ROOT,OPTION,SEP];
OUTPUTPATH = ['../data/aggregate/',OPTION,'/'];
EnvironPara.OUTPUTPATH_PIC = OUTPUTPATH;
mkdir(OUTPUTPATH);
SITENAME = ['US',OPTION];
StartDay = datenum([YEAR,1,1]);

list = dir(DIRPATH);
fnames = {list.name}.';
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(fnames, str, 'once'));

%%Regional; the boundary box
LonMax = -60;
LonMin = -135;
LatMax = 60;
LatMin = 15;

% LonMax = Inf;
% LonMin = -Inf;
% LatMax = Inf;
% LatMin = -Inf;

%% load grid cell locations --- this is the unified points
diary([OUTPUTPATH,OPTION,'_',SITENAME,'_',num2str(YEAR),'_',num2str(YEAR),'.txt']);
if(exist([OUTPUTPATH,SITENAME,'Site.mat'],'file'))
    SiteData = load([OUTPUTPATH,SITENAME,'Site.mat']);%SiteData, TileName --- tile name is order of tiles while creating the SiteData
    SiteData = SiteData.SiteData;
else
    disp('creating site data!!!');
    [data,SiteData,date] = ReadOMIAIData_function([DIRPATH,list(4).name],OPTION,VARIABLE_LIST{1});
    Index = (SiteData.Lon<LonMax & SiteData.Lon>LonMin & SiteData.Lat<LatMax & SiteData.Lat>LatMin);
    SiteData = SiteData(Index,:);
    save([OUTPUTPATH,SITENAME,'Site.mat'],'SiteData');
end

N_Site = size(SiteData,1);
N_Day = yeardays(YEAR);
ResultAll = nan(N_Day,N_Site,length(VARIABLE_LIST));

%% read all measurments in that year
for i=1:N_Day
    CurrentDay = StartDay+i-1;
    fprintf('%s\n',datestr(CurrentDay,'yyyy-mm-dd'));
    
    % find all files on that day
    FileName_Template = [TEMPLATE,'_',num2str(year(CurrentDay),'%02d'),'m',num2str(month(CurrentDay),'%02d'),num2str(day(CurrentDay),'%02d'),'_v003-','*'];
    FileName = fnames(FIND(FileName_Template));
    
    if(isempty(FileName))
        continue;
    end
    
    % if multiple files for a single day
    Temp = nan(N_Site,length(VARIABLE_LIST),length(FileName));
    for j=1:length(FileName)
        fprintf('\t%s\n',FileName{j});
        for k = 1:length(VARIABLE_LIST)
            % this function reads OMI L3 data and extract variables of interests
            try
                [data,SiteData_Temp,date] = ReadOMIAIData_function([DIRPATH,FileName{j}],OPTION,VARIABLE_LIST{k});
            catch exception
                fprintf('reading warning:%s\n',exception.message);
                continue;
            end
            
            if(isempty(data))
                disp('option and variable combination error!');
                return;
            end
            
            % match the data to consistent grid cells
            [LiA,LocB] = ismember(SiteData_Temp.SiteCode,SiteData.SiteCode);
            LocB = LocB(LiA);
            Temp(LocB,k,j) = data(LiA);
        end
    end
    
    % conduct descriptive analysis
    Temp = nanmean(Temp,3);
    for k = 1:length(VARIABLE_LIST)
        fprintf('\t\t%s---\tnan percentage:%f\tmax:%f\tmean:%f\tmin:%f\n',VARIABLE_LIST{k},sum(isnan(Temp(:,k)))/length(Temp(:,k)),nanmax(Temp(:,k)),nanmean(Temp(:,k)),nanmin(Temp(:,k)));        
    end
    ResultAll(i,:,:) = Temp;
end

% save result by variable
for k = 1:length(VARIABLE_LIST)
    Result = ResultAll(:,:,k);
    save([OUTPUTPATH,OPTION,'_',VARIABLE_LIST{k},'_',SITENAME,'_',num2str(YEAR),'_',num2str(YEAR),'.mat'],'Result');
    Visualization_USResult_1([OPTION,'_',VARIABLE_LIST{k},'_',SITENAME,'_',num2str(YEAR),'_',num2str(YEAR)],nanmean(Result),SiteData,EnvironPara);
%     Visualization_USResult_1(FileName,Result,SiteData,EnvironPara)
end
diary off
