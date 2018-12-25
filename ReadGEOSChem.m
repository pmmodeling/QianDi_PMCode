%% main function processing MERRAaer raw data; 
% this code reads MERRAaer data in each; and extract five variables of
% interests: Hydrophilic Black Carbon, Hydrophobic Black Carbon, Hydrophilic Organic Carbon (Particulate Matter), Hydrophobic Organic Carbon (Particulate Matter), and Sulphate aerosol

%% input file: 
% USMEERA2aerSite.mat, and USMEERA2aerSite_North_America_Equidistant_Conic
% put inside output folder

%% parameters:
% YEAR = 2000;

%% example: 
% ReadMERRA2aer(2000)

function ReadGEOSChem(YEAR)
% YEAR = 2004;

%% macro definition
SITENAME = 'USGEOSChem';

% path
DIRPATH = '../data/unprocessed/GEOSChem/';
OUTPUTPATH = '../data/aggregate/GEOSChem/';

EnvironPara = struct();
EnvironPara.states = geoshape(shaperead('usastatehi', 'UseGeoCoords', true));
EnvironPara.S_back = geoshape(shaperead('../data/shapefile/US_WGS_clipcoast_2.shp', 'UseGeoCoords', true));
EnvironPara.US_North_America_Equidistant_Conic = '../data/shapefile/US_North_America_Equidistant_Conic.shp';
EnvironPara.US = '../data/shapefile/US.shp';
GCS = 'North_America_Equidistant_Conic';

mkdir(OUTPUTPATH);
diary([OUTPUTPATH,'ReadGEOSChem_',num2str(YEAR),'.txt']);
EnvironPara.OUTPUTPATH_PIC = OUTPUTPATH;

list = dir(DIRPATH);
fnames = {list.name}.';
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(fnames, str, 'once'));

%% create site data for current grid
if(exist([OUTPUTPATH,SITENAME,'Site','.mat'],'file'))
    SiteData = LoadData_function([OUTPUTPATH,SITENAME,'Site','.mat']);%SiteData --- uses lat/lon; since raw data are in lat/lon
    SiteData = SiteData.SiteData;
    Lon_Query = SiteData.Lon;
    Lat_Query = SiteData.Lat;
else
    lat = ncread([DIRPATH,'MDA8_03.05x0625_NA.2000.nc'],'lat');
    lon = ncread([DIRPATH,'MDA8_03.05x0625_NA.2000.nc'],'lon');
    [Y,X] = meshgrid(lat,lon);
     
    X = double(reshape(X,[size(X,1)*size(X,2),1]));
    Y = double(reshape(Y,[size(X,1)*size(X,2),1]));
    SiteData = table((1:1:(size(X,1)*size(X,2)))',Y,X,'VariableNames',{'SiteCode','Lat','Lon'});
    save([OUTPUTPATH,SITENAME,'Site','.mat'],'SiteData');
    writetable(SiteData,[OUTPUTPATH,SITENAME,'Site','.csv']);
end

% NO2
Temp = double(ncread([DIRPATH,'MDA8_03.05x0625_NA.',num2str(YEAR),'.nc'],'O3'));
Result = reshape(Temp,[size(SiteData,1),size(Temp,3)])';
save([OUTPUTPATH,'GEOSChem_O3_',SITENAME,'_',num2str(YEAR),'_',num2str(YEAR),'.mat'],'Result');
Visualization_USResult_1(['GEOSChem_O3_',datestr(datenum(YEAR,1,1),'yyyymmdd'),'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),'.tif'],nanmean(Result),SiteData,EnvironPara);

Temp = double(ncread([DIRPATH,'NO2.05x0625_NA.',num2str(YEAR),'.nc'],'NO2'));
Result = reshape(Temp,[size(SiteData,1),size(Temp,3)])';
save([OUTPUTPATH,'GEOSChem_NO2_',SITENAME,'_',num2str(YEAR),'_',num2str(YEAR),'.mat'],'Result');
Visualization_USResult_1(['GEOSChem_NO2_',datestr(datenum(YEAR,1,1),'yyyymmdd'),'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),'.tif'],nanmean(Result),SiteData,EnvironPara);

Temp = double(ncread([DIRPATH,'PM25.05x0625_NA.',num2str(YEAR),'.nc'],'PM25'));
Result = reshape(Temp,[size(SiteData,1),size(Temp,3)])';
save([OUTPUTPATH,'GEOSChem_PM25_',SITENAME,'_',num2str(YEAR),'_',num2str(YEAR),'.mat'],'Result');
Visualization_USResult_1(['GEOSChem_PM25_',datestr(datenum(YEAR,1,1),'yyyymmdd'),'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),'.tif'],nanmean(Result),SiteData,EnvironPara);

Temp = double(ncread([DIRPATH,'PM25.05x0625_NA.',num2str(YEAR),'.nc'],'SO4'));
Result = reshape(Temp,[size(SiteData,1),size(Temp,3)])';
save([OUTPUTPATH,'GEOSChem_SO4_',SITENAME,'_',num2str(YEAR),'_',num2str(YEAR),'.mat'],'Result');
Visualization_USResult_1(['GEOSChem_SO4_',datestr(datenum(YEAR,1,1),'yyyymmdd'),'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),'.tif'],nanmean(Result),SiteData,EnvironPara);

Temp = double(ncread([DIRPATH,'PM25.05x0625_NA.',num2str(YEAR),'.nc'],'NH4'));
Result = reshape(Temp,[size(SiteData,1),size(Temp,3)])';
save([OUTPUTPATH,'GEOSChem_NH4_',SITENAME,'_',num2str(YEAR),'_',num2str(YEAR),'.mat'],'Result');
Visualization_USResult_1(['GEOSChem_NH4_',datestr(datenum(YEAR,1,1),'yyyymmdd'),'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),'.tif'],nanmean(Result),SiteData,EnvironPara);

Temp = double(ncread([DIRPATH,'PM25.05x0625_NA.',num2str(YEAR),'.nc'],'BC'));
Result = reshape(Temp,[size(SiteData,1),size(Temp,3)])';
save([OUTPUTPATH,'GEOSChem_BC_',SITENAME,'_',num2str(YEAR),'_',num2str(YEAR),'.mat'],'Result');
Visualization_USResult_1(['GEOSChem_BC_',datestr(datenum(YEAR,1,1),'yyyymmdd'),'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),'.tif'],nanmean(Result),SiteData,EnvironPara);

Temp = double(ncread([DIRPATH,'PM25.05x0625_NA.',num2str(YEAR),'.nc'],'OA'));
Result = reshape(Temp,[size(SiteData,1),size(Temp,3)])';
save([OUTPUTPATH,'GEOSChem_OA_',SITENAME,'_',num2str(YEAR),'_',num2str(YEAR),'.mat'],'Result');
Visualization_USResult_1(['GEOSChem_OA_',datestr(datenum(YEAR,1,1),'yyyymmdd'),'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),'.tif'],nanmean(Result),SiteData,EnvironPara);
