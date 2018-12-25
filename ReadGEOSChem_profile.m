%% main function processing GEOS-Chem vertical profile data

%% input file: 
% USMEERA2aerSite.mat, and USMEERA2aerSite_North_America_Equidistant_Conic
% put inside output folder

%% parameters:
% YEAR = 2000;

%% example: 
% ReadMERRA2aer(2000)

% function ReadGEOSChem_profile(YEAR)
YEAR = 2004;

%% macro definition
SITENAME = 'USGEOSChem';

% path
DIRPATH = 'D:\GoogleDrive\Research\USTemperature\raw_data\data\unprocessed\GEOSChem\';
% DIRPATH = '../data/unprocessed/GEOSChem/';
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
for TempName = {'O3','PM25','NO2'}
    Temp = double(ncread([DIRPATH,'PM25_NO2_O3.PROFILE.05x0625_NA.',num2str(YEAR),'.nc'],TempName{1}));
    if(2004 == YEAR)
        % for year 2004, take Jan of 2004 to replace Dec 2013
        Temp1 = double(ncread([DIRPATH,'PM25_NO2_O3.PROFILE.05x0625_NA.',num2str(YEAR),'.nc'],TempName{1}));
        Temp1 = Temp1(:,:,:,1);
    else
        % take Dec of the previous year
        Temp1 = double(ncread([DIRPATH,'PM25_NO2_O3.PROFILE.05x0625_NA.',num2str(YEAR-1),'.nc'],TempName{1}));
        Temp1 = Temp1(:,:,:,12);
    end

    if(2016 == YEAR)
        % for year 2016, take Dec 2016 to replace Jan 2017
        Temp2 = double(ncread([DIRPATH,'PM25_NO2_O3.PROFILE.05x0625_NA.',num2str(YEAR),'.nc'],TempName{1}));
        Temp2 = Temp2(:,:,:,12);
    else
        % take Jan of the next year
        Temp2 = double(ncread([DIRPATH,'PM25_NO2_O3.PROFILE.05x0625_NA.',num2str(YEAR+1),'.nc'],TempName{1}));
        Temp2 = Temp2(:,:,:,1);
    end

    Temp = cat(4,Temp1,Temp,Temp2);
    Temp_scaling = Temp(:,:,1,:)./sum(Temp,3);
    Temp_scaling = reshape(Temp_scaling,[size(Temp_scaling,1),size(Temp_scaling,2),size(Temp_scaling,4)]);
    Temp_scaling = reshape(Temp_scaling,[size(Temp_scaling,1)*size(Temp_scaling,2),size(Temp_scaling,3)]);
    Temp_scaling = Temp_scaling';

    xi = datenum(YEAR-1,12,1):1:datenum(YEAR+1,1,31);
    x = [15,46,77,105,136,166,197,227,258,289,319,350,380,411]+datenum(YEAR-1,12,1)-1;
    Temp_scaling = interp1(x,Temp_scaling,xi);
    Temp_scaling = Temp_scaling(32:end-31,:);

    Result = Temp_scaling;
    save([OUTPUTPATH,'GEOSChem_',TempName{1},'Scaling_',SITENAME,'_',num2str(YEAR),'_',num2str(YEAR),'.mat'],'Result');
%     Visualization_USResult_1(['GEOSChem_',TempName{1},'Scaling_',datestr(datenum(YEAR,1,1),'yyyymmdd'),'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),'.tif'],nanmean(Result),SiteData,EnvironPara);
end