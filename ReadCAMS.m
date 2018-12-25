%% main function processing CAMS raw data; 
% 

function ReadCAMS(YEAR)
% YEAR = 2004;

%% macro definition
SITENAME = 'USCAMS';

% path
DIRPATH = 'Z:\Heresh\1 Raw data\9 NO2 CAMS reanalysis 2003-2016\Downloaded data\';
OUTPUTPATH = 'D:\Downloads\';

EnvironPara = struct();
EnvironPara.states = geoshape(shaperead('usastatehi', 'UseGeoCoords', true));
EnvironPara.S_back = geoshape(shaperead('D:/GoogleDrive/Research/USTemperature/raw_data/data/shapefile/US_WGS_clipcoast_2.shp', 'UseGeoCoords', true));
EnvironPara.US_North_America_Equidistant_Conic = 'D:/GoogleDrive/Research/USTemperature/raw_data/data/shapefile/US_North_America_Equidistant_Conic.shp';
EnvironPara.US = 'D:/GoogleDrive/Research/USTemperature/raw_data/data/shapefile/US.shp';
GCS = 'North_America_Equidistant_Conic';

mkdir(OUTPUTPATH);
diary([OUTPUTPATH,'ReadCAMS_',num2str(YEAR),'.txt']);
EnvironPara.OUTPUTPATH_PIC = OUTPUTPATH;

% list = dir(DIRPATH);
% fnames = {list.name}.';
% FIND = @(str) cellfun(@(c) ~isempty(c), regexp(fnames, str, 'once'));

%% create site data for current grid
if(exist([OUTPUTPATH,SITENAME,'Site','.mat'],'file'))
    SiteData = LoadData_function([OUTPUTPATH,SITENAME,'Site','.mat']);%SiteData --- uses lat/lon; since raw data are in lat/lon
    SiteData = SiteData.SiteData;
    Lon_Query = SiteData.Lon;
    Lat_Query = SiteData.Lat;
else
    lat = ncread([DIRPATH,'1-4 2003 CAMS RA NO2.nc'],'latitude');
    lat = lat(lat>24 & lat < 52);% 306 - 528: total:223
    lon = ncread([DIRPATH,'1-4 2003 CAMS RA NO2.nc'],'longitude');
    lon = rem((lon+180),360)-180;%1858 - 2368; total: 511
    lon = lon(lon>-128 & lon < -64);
    [Y,X] = meshgrid(lat,lon);
     
    X = double(reshape(X,[size(X,1)*size(X,2),1]));
    Y = double(reshape(Y,[size(X,1)*size(X,2),1]));
    SiteData = table((1:1:(size(X,1)*size(X,2)))',Y,X,'VariableNames',{'SiteCode','Lat','Lon'});
    save([OUTPUTPATH,SITENAME,'Site','.mat'],'SiteData');
    writetable(SiteData,[OUTPUTPATH,SITENAME,'Site','.csv']);
end

% read 1 - 4 month; just read grid cells inside U.S.
Temp = double(ncread([DIRPATH,'1-4 ',num2str(YEAR),' CAMS RA NO2.nc'],'tcno2',[1858,306,1],[511,223,Inf]));
Temp1 = permute(Temp,[3,2,1]);
Temp21 = permute(nanmean(reshape(Temp1,[8,size(Temp1,1)/8,size(Temp1,2),size(Temp1,3)]),1),[4,3,2,1]);

% read 5 - 8 month; just read grid cells inside U.S.
Temp = double(ncread([DIRPATH,'5-8 ',num2str(YEAR),' CAMS RA NO2.nc'],'tcno2',[1858,306,1],[511,223,Inf]));
Temp1 = permute(Temp,[3,2,1]);
Temp22 = permute(nanmean(reshape(Temp1,[8,size(Temp1,1)/8,size(Temp1,2),size(Temp1,3)]),1),[4,3,2,1]);

% read 9 - 12 month; just read grid cells inside U.S.
Temp = double(ncread([DIRPATH,'9-12 ',num2str(YEAR),' CAMS RA NO2.nc'],'tcno2',[1858,306,1],[511,223,Inf]));
Temp1 = permute(Temp,[3,2,1]);
Temp23 = permute(nanmean(reshape(Temp1,[8,size(Temp1,1)/8,size(Temp1,2),size(Temp1,3)]),1),[4,3,2,1]);

Result = cat(3,Temp21,Temp22,Temp23);
Result = reshape(Result,[size(Result,1)*size(Result,2),size(Result,3)])';
save([OUTPUTPATH,'CAMS_NO2_',SITENAME,'_',num2str(YEAR),'_',num2str(YEAR),'.mat'],'Result');
Visualization_USResult_1(['CAMS_NO2_',datestr(datenum(YEAR,1,1),'yyyymmdd'),'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),'.tif'],nanmean(Result),SiteData,EnvironPara);
