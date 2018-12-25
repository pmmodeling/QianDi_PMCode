%% read forest fire data; take the monthly forest fire emission data; 
% make sure  make sure USGFEDSite.mat and 
% USGFEDSite_North_America_Equidistant_Conic.mat have been uploaded to the
% output folder

%% parameters:
% YEAR = 2000;% year of files to be read

%% example 
% ReadForestfire(2000)

%% code
function ReadForestfire(YEAR)

SITENAME = 'USGFED';

%% path
DIRPATH = '../data/unprocessed/GFED/';
OUTPUTPATH = '../data/aggregate/GFED/';
SEP = '/';

% file path of the shape file, used for making maps
EnvironPara = struct();
EnvironPara.states = geoshape(shaperead('usastatehi', 'UseGeoCoords', true));
EnvironPara.S_back = geoshape(shaperead('../data/shapefile/US_WGS_clipcoast_2.shp', 'UseGeoCoords', true));
EnvironPara.US_North_America_Equidistant_Conic = '../data/shapefile/US_North_America_Equidistant_Conic.shp';
EnvironPara.US = '../data/shapefile/US.shp';
GCS = 'North_America_Equidistant_Conic';

mkdir(OUTPUTPATH);
EnvironPara.OUTPUTPATH_PIC = OUTPUTPATH;
FileName = [DIRPATH,'GFED4.1s_',num2str(2000),'.hdf5'];

%% create site data for current grid
if(exist([OUTPUTPATH,SITENAME,'Site','.mat'],'file'))
    SiteData = LoadData_function([OUTPUTPATH,SITENAME,'Site','.mat']);%SiteData --- uses lat/lon; since raw data are in lat/lon
    SiteData_Temp = SiteData.SiteData;
else
%     disp('creating site data!!!');
    Lat = double(hdf5read(FileName,'/lat'));
    Lon = double(hdf5read(FileName,'/lon'));
    Lat = reshape(Lat,[size(Lat,1)*size(Lat,2),1]);
    Lon = reshape(Lon,[size(Lon,1)*size(Lon,2),1]);
    SiteData = table((1:1:size(Lon,1))',Lat,Lon,'VariableNames',{'SiteCode','Lat','Lon'});
    save([OUTPUTPATH,SITENAME,'Site','.mat'],'SiteData');
    writetable(SiteData,[OUTPUTPATH,SITENAME,'Site','.csv']);
end

%% read the annual data and take monthly forest fire emission
% interpolate monthly emission data to daily emission
TempAll = [];
x = [];
for i = YEAR-1:YEAR+1
    FileName = [DIRPATH,'GFED4.1s_',num2str(i),'.hdf5'];
    fprintf('reading...%s\n',FileName);
    % read monthly forest fire emission
    for j = 1:12
        try
            Temp = double(hdf5read(FileName,['/emissions/',sprintf('%02d',j),'/C']));
            Temp = reshape(Temp,[size(Temp,1)*size(Temp,2),1])';
        catch
            fprintf('error reading...%s\n',FileName);
        end
        
        TempAll = cat(1,TempAll,Temp);
        x = cat(1,x,datenum(i,j,15));
    end
end

% interpolate monthly data to daily level
xi = datenum(YEAR-1,1,1):1:datenum(YEAR+1,12,31);
ResultAll = interp1(x,TempAll,xi);
ResultAll(end-15:end,:)=repmat(ResultAll(end-16,:),[16,1]);

%% save results
Result = ResultAll(datenum(YEAR,1,1)-datenum(YEAR-1,1,1)+1:datenum(YEAR,12,31)-datenum(YEAR-1,1,1)+1,:);
FileName = [OUTPUTPATH,'GFEDFireCarbon_',SITENAME,'_',num2str(YEAR),'_',num2str(YEAR),'.mat'];
fprintf('saving...%s\n',FileName);
save(FileName,'Result','-v7.3'); 
Visualization_USResult_1(['GFEDFireCarbon_',SITENAME,'_',num2str(YEAR),'_',num2str(YEAR),'.tif'],nanmean(Result),SiteData_Temp,EnvironPara);
