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

%% code 
function ReadMERRA2aer(YEAR)

%% macro definition
SITENAME = 'USMERRA2aer';

% path
DIRPATH = '../data/unprocessed/MERRA2aer/';
OUTPUTPATH = '../data/aggregate/MERRA2aer/';

EnvironPara = struct();
EnvironPara.states = geoshape(shaperead('usastatehi', 'UseGeoCoords', true));
EnvironPara.S_back = geoshape(shaperead('../data/shapefile/US_WGS_clipcoast_2.shp', 'UseGeoCoords', true));
EnvironPara.US_North_America_Equidistant_Conic = '../data/shapefile/US_North_America_Equidistant_Conic.shp';
EnvironPara.US = '../data/shapefile/US.shp';
GCS = 'North_America_Equidistant_Conic';

mkdir(OUTPUTPATH);
diary([OUTPUTPATH,'ReadMERRA2aer_',num2str(YEAR),'.txt']);
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
    lat = ncread([DIRPATH,strrep('MERRA2_200.inst3_3d_aer_Nv.$Date$.SUB.nc4','$Date$',datestr(datenum(2000,1,1),'yyyymmdd'))],'lat');
    lon = ncread([DIRPATH,strrep('MERRA2_200.inst3_3d_aer_Nv.$Date$.SUB.nc4','$Date$',datestr(datenum(2000,1,1),'yyyymmdd'))],'lon');
    [Y,X] = meshgrid(lat,lon);
     
    X = double(reshape(X,[size(X,1)*size(X,2),1]));
    Y = double(reshape(Y,[size(X,1)*size(X,2),1]));
    SiteData = table((1:1:(size(X,1)*size(X,2)))',Y,X,'VariableNames',{'SiteCode','Lat','Lon'});
    save([OUTPUTPATH,SITENAME,'Site','.mat'],'SiteData');
    writetable(SiteData,[OUTPUTPATH,SITENAME,'Site','.csv']);
    
    Lon_Query = SiteData.Lon;
    Lat_Query = SiteData.Lat;
end

Result_BCPHILIC = nan(yeardays(YEAR),size(SiteData,1));
Result_BCPHOBIC = nan(yeardays(YEAR),size(SiteData,1));
Result_OCPHILIC = nan(yeardays(YEAR),size(SiteData,1));
Result_OCPHOBIC = nan(yeardays(YEAR),size(SiteData,1));
Result_SO4 = nan(yeardays(YEAR),size(SiteData,1));

%% read data and extract variables of interests
StartDay = datenum(YEAR,1,1);
for i = 1:yeardays(YEAR)
  
   CurrentDay = StartDay + i - 1;
   
   try
       FileName = fnames(FIND(['MERRA2_[1234]00.inst3_3d_aer_Nv.',datestr(CurrentDay,'yyyymmdd'),'.SUB.nc4']));   FileName = FileName{1};
       fprintf('processing...%s\n',FileName);
       
       % Hydrophilic Black Carbon 
       Temp = double(ncread([DIRPATH,FileName],'BCPHILIC'));%112*71*72
       Result_BCPHILIC(i,:) = reshape(Temp(:,:,1),[1,size(Temp,1)*size(Temp,2)]);     
       
       % Hydrophobic Black Carbon 
       Temp = double(ncread([DIRPATH,FileName],'BCPHOBIC'));
       Result_BCPHOBIC(i,:) = reshape(Temp(:,:,1),[1,size(Temp,1)*size(Temp,2)]);
       
       % Hydrophilic Organic Carbon (Particulate Matter) 
       Temp = double(ncread([DIRPATH,FileName],'OCPHILIC'));
       Result_OCPHILIC(i,:) = reshape(Temp(:,:,1),[1,size(Temp,1)*size(Temp,2)]);
       
       % Hydrophobic Organic Carbon (Particulate Matter) 
       Temp = double(ncread([DIRPATH,FileName],'OCPHOBIC'));
       Result_OCPHOBIC(i,:) = reshape(Temp(:,:,1),[1,size(Temp,1)*size(Temp,2)]);
       
       % Sulphate aerosol 
       Temp = double(ncread([DIRPATH,FileName],'SO4'));
       Result_SO4(i,:) = reshape(Temp(:,:,1),[1,size(Temp,1)*size(Temp,2)]);

   catch exception
       fprintf('%s\n',exception.message);
   end
   
end

% make map and store results
Result = Result_BCPHILIC;
save([OUTPUTPATH,'MERRA2aer_BCPHILIC_',SITENAME,'_',num2str(YEAR),'_',num2str(YEAR),'.mat'],'Result');
Visualization_USResult_1(['MERRA2aer_BCPHILIC_',datestr(datenum(YEAR,1,1),'yyyymmdd'),'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),'.tif'],nanmean(Result),SiteData,EnvironPara);

Result = Result_BCPHOBIC;
save([OUTPUTPATH,'MERRA2aer_BCPHOBIC_',SITENAME,'_',num2str(YEAR),'_',num2str(YEAR),'.mat'],'Result');
Visualization_USResult_1(['MERRA2aer_BCPHOBIC_',datestr(datenum(YEAR,1,1),'yyyymmdd'),'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),'.tif'],nanmean(Result),SiteData,EnvironPara);

Result = Result_OCPHILIC;
save([OUTPUTPATH,'MERRA2aer_OCPHILIC_',SITENAME,'_',num2str(YEAR),'_',num2str(YEAR),'.mat'],'Result');
Visualization_USResult_1(['MERRA2aer_OCPHILIC_',datestr(datenum(YEAR,1,1),'yyyymmdd'),'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),'.tif'],nanmean(Result),SiteData,EnvironPara);

Result = Result_OCPHOBIC;
save([OUTPUTPATH,'MERRA2aer_OCPHOBIC_',SITENAME,'_',num2str(YEAR),'_',num2str(YEAR),'.mat'],'Result');
Visualization_USResult_1(['MERRA2aer_OCPHOBIC_',datestr(datenum(YEAR,1,1),'yyyymmdd'),'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),'.tif'],nanmean(Result),SiteData,EnvironPara);

Result = Result_SO4;
save([OUTPUTPATH,'MERRA2aer_SO4_',SITENAME,'_',num2str(YEAR),'_',num2str(YEAR),'.mat'],'Result');
Visualization_USResult_1(['MERRA2aer_SO4_',datestr(datenum(YEAR,1,1),'yyyymmdd'),'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),'.tif'],nanmean(Result),SiteData,EnvironPara);
