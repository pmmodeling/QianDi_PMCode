%% main function to ozone vertical profile
% this OMI data product is L2 data, which means the grid cells are not
% consistent but change from day to day. I used the grid cells from
% OMAEROe; interpolated the L2 data to the consistent OMAEROe grid cells.

%% parameter:
% YEAR = 2015;

%% example: 
% ReadOMO3PR_Main(2015)

%% code
function ReadOMO3PR_Main(YEAR)


%% O3 profile
OPTION = 'OMO3PR';
TEMPLATE = 'OMI-Aura_L2-OMO3PR';
N_Layer = 18;
%% set to yes and interpolate for L2 data
IsInterp = true;

%% path
DIRPATHROOT = '../data/unprocessed/';
SEP = '/';

%% file path of shape file
EnvironPara = struct();
EnvironPara.states = geoshape(shaperead('usastatehi', 'UseGeoCoords', true));
EnvironPara.S_back = geoshape(shaperead('../data/shapefile/US_WGS_clipcoast_2.shp', 'UseGeoCoords', true));
EnvironPara.US_North_America_Equidistant_Conic = '../data/shapefile/US_North_America_Equidistant_Conic.shp';
EnvironPara.US = '../data/shapefile/US.shp';
GCS = 'North_America_Equidistant_Conic';

DIRPATH = [DIRPATHROOT,OPTION,SEP];
OUTPUTPATH = ['../data/aggregate/',OPTION,'/'];

diary([OUTPUTPATH,'ReadOMIAIData_L3_Main_',OPTION,num2str(YEAR),'.txt']);
fprintf('OPTION:%s\n',OPTION);
fprintf('TEMPLATE:%s\n',TEMPLATE);

EnvironPara.OUTPUTPATH_PIC = OUTPUTPATH;
StartDay = datenum([YEAR,1,1]);
N_Day = yeardays(YEAR);

list = dir(DIRPATH);
fnames = {list.name}.';
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(fnames, str, 'all'));

%% Regional
LonMax = -60;
LonMin = -135;
LatMax = 60;
LatMin = 15;

%% IMPORTANT!!! THIS IS THE GRID CELL USED to aggregate data at
try
    load([OUTPUTPATH,'US',OPTION,'Site','','.mat']);
    disp('loading SiteData...');
    N_Site = size(SiteData,1);
    Result_Lat = SiteData.Lat;
    Result_Lon = SiteData.Lon;
catch 
    disp('creating SiteData!');
    disp('make sure USOMO3PRSite.mat and USOMO3PRSite_North_America_Equidistant_Conic.mat have been uploaded to ./raw_data/data/aggregate/OMO3PR/');
    %%try read the first file to get size of output
    [Result_Lon,Result_Lat] = meshgrid(-134.875:0.25:-60.125,15.125:0.25:59.875);
%     [Result_Lon,Result_Lat] = meshgrid(-180:1:180,-90:1:90);
    N_Site = size(Result_Lon,1)*size(Result_Lon,2);
    N_Day = yeardays(YEAR);
    Result_Lon = reshape(Result_Lon,[size(Result_Lon,1)*size(Result_Lon,2),1]);
    Result_Lat = reshape(Result_Lat,[size(Result_Lat,1)*size(Result_Lat,2),1]);
    SiteData = table((1:1:length(Result_Lon))',Result_Lat,Result_Lon,'VariableNames',{'SiteCode' 'Lat' 'Lon'});
    
    save([OUTPUTPATH,SITEDATA_NAME,'Site.mat'],'SiteData');
    return;
end

Result = nan(N_Day,N_Site);

%% read daily ozone vertical profile measurements
for i=1:N_Day
    CurrentDay = StartDay+i-1;
    
    % fine all files
    FileName_Template = [TEMPLATE,'_',num2str(year(CurrentDay),'%02d'),'m',num2str(month(CurrentDay),'%02d'),num2str(day(CurrentDay),'%02d'),'t[\d]{4}-o[\d]{5}_v003-','*'];        
    FileName = fnames(FIND(FileName_Template));
    if(length(FileName) == 0)
       continue; 
    end
    
    % begin to process!!!
    fprintf('%s\n',datestr(CurrentDay,'yyyy-mm-dd'));
    fprintf('%d\n',length(FileName));
    Temp = nan(N_Site,length(FileName));
    Temp_error = nan(N_Site,length(FileName));
    %  read all files from that day
    for j=1:length(FileName)
        fprintf('\t%s\n',FileName{j});
        [data,lon,lat,date,data_1,data_2] = ReadOMIAIData_L3([DIRPATH,FileName{j}],OPTION);

        Temp(:,j) = InterpMyData_2(data(1,:)./sum(data,1),lon,lat,Result_Lon,Result_Lat,'')';
%         Temp(:,j) = InterpMyData_2(sum(data,1),lon,lat,Result_Lon,Result_Lat,'')';
        Temp_error(:,j) = InterpMyData_2(data_2',lon,lat,Result_Lon,Result_Lat,'')';
        
        fprintf('\tnan:%f\tmax:%f\tmin:%f\n',sum(isnan(Temp(:,j)))/length(Temp(:,j)),nanmax(Temp(:,j)),nanmin(Temp(:,j)));
        fprintf('\t%s\n',datestr(date,'yyyymmdd'));
        if(date~=CurrentDay)
           disp('date mismatch!'); 
        end       
    end

    % weighted sum based on error
%     Temp_error = 1./Temp_error;
%     Index = isnan(Temp)&isnan(Temp_error);
%     Temp(Index) = nan;
%     Temp_error(Index) = 0;
%     Temp_error1 = Temp_error./repmat(sum(Temp_error,2),1,length(FileName));
%     Result(i,:) = sum(Temp.*Temp_error1,2)';
    
    %instead of taking the mean, we take the one with smallest MSE (based on RootMeanSquareErrorOfFit)
%     Temp_error = Temp_error';
%     IndexError = Temp_error == repmat(min(Temp_error),[size(Temp_error,1),1]);
%     Temp = Temp';
%     Temp1 = nan(size(Temp));
%     Temp1(IndexError) = Temp(IndexError);
%     Result(i,:) = nansum(Temp1,1);
    
    % naive one --- calculate the percentage of ground level ozone in the
    % total column ozone
    Result(i,:) = nanmean(Temp,2)';
    fprintf('\tnan:%f\t%f\t%f\t%f\t%f\t%f\t%f\n',sum(isnan(Result(i,:)))/length(Result(i,:)),...
    nanmax(Result(i,:)),nanmin(Result(i,:)),nanmax(Result_Lon),nanmin(Result_Lon),nanmax(Result_Lat),nanmin(Result_Lat));  

%         %% save
%         TempResult = nanmean(double(Temp),3);
%         TempResult1 = nanmean(double(Temp1),3);
%         TempResult2 = nanmean(double(Temp2),3);
%         fprintf('\tdata:%f\t%f\t%f\ttemp:%f\t%f\t%f\talt:%f\t%f\t%f\n',...
%         sum(sum(isnan(TempResult)))/(size(TempResult,2)*size(TempResult,3)),nanmax(nanmax(TempResult)),nanmin(nanmin(TempResult)),...
%         sum(sum(isnan(TempResult1)))/(size(TempResult1,2)*size(TempResult1,3)),nanmax(nanmax(TempResult1)),nanmin(nanmin(TempResult1)),...
%         sum(sum(isnan(TempResult2)))/(size(TempResult2,2)*size(TempResult2,3)),nanmax(nanmax(TempResult2)),nanmin(nanmin(TempResult2))...
%         );
%         save(TempFileName,'TempResult','TempResult1','TempResult2');
%         
%         delete([TempFileName,'.part']);
end

%% save --- ADD SITENAME HERE!!!
save([OUTPUTPATH,OPTION,'_',num2str(YEAR),'_',num2str(YEAR),'.mat'],'Result');
Visualization_USResult_1([OPTION,'_',num2str(YEAR),'_',num2str(YEAR)],nanmean(Result,1),SiteData,EnvironPara);
diary off;
fclose all;

% %% assemble and save
% %ozone profile
% Result = nan(N_Day,N_Site,N_Layer);
% Result1 = nan(N_Day,N_Site,N_Layer);
% Result2 = nan(N_Day,N_Site,N_Layer+1);
% 
% i = 1;
% while i<=N_Day
%     CurrentDay = StartDay+i-1;
%     fprintf('assembling...%s\n',datestr(CurrentDay,'yyyy-mm-dd'));
% 
%     TempFileName = [OUTPUTPATH,OPTION,'_',datestr(CurrentDay,'yyyymmdd'),'_',datestr(CurrentDay,'yyyymmdd'),'.mat'];
%     try
%         load(TempFileName,'TempResult','TempResult1','TempResult2');
%         Result(i,:,:) = TempResult;
%         Result1(i,:,:) = TempResult1;
%         Result2(i,:,:) = TempResult2;
%         i = i+1;
%     catch exception
%        fprintf('%s\n',exception.message);
%        return;
%     end
% 
% end

% save([OUTPUTPATH,OPTION,'O3','_',num2str(YEAR),'_',num2str(YEAR),'.mat'],'Result','Result_Lon','Result_Lat');
% Result = Result1;
% save([OUTPUTPATH,OPTION,'Temperature','_',num2str(YEAR),'_',num2str(YEAR),'.mat'],'Result','Result_Lon','Result_Lat');
% Result = Result2;
% save([OUTPUTPATH,OPTION,'Altitude','_',num2str(YEAR),'_',num2str(YEAR),'.mat'],'Result','Result_Lon','Result_Lat');

diary off;
fclose all;

end


%% extract lon/lon, and other variable of interests, such as O3 (with vertical layer),
% O3 is the most important variable; other variables are not used in current modeling
function [data,lon,lat,date,data_1,data_2] = ReadOMIAIData_L3(FILE_NAME,OPTION)
% % disp('ReadOMIAIData_L3');

% % Open the HDF5 File.
% FILE_NAME = 'D:\Data\Data_AOD\OMI-Aura_L3-OMDOAO3e_2009m0708_v003-2011m1109t111320.he5';
% OPTION = 'OMAEROe';
% h5disp(FILE_NAME);
% info = h5info(FILE_NAME);

% FILE_NAME = 'A:\QianDi\OMTO3e\OMI-Aura_L3-OMTO3e_2004m1001_v003-2012m0409t101417.he5';
% OPTION = 'OMTO3e';
% h5disp(FILE_NAME);

data = nan;
data_1 = nan;
data_2 = nan;

    
Templon = double(h5read(FILE_NAME,'/HDFEOS/SWATHS/O3Profile/Geolocation Fields/Longitude'));
lon = nan([1,size(Templon)]);
lon(1,:,:) = Templon;
lon(lon<-1e30)=nan;

Templat = double(h5read(FILE_NAME,'/HDFEOS/SWATHS/O3Profile/Geolocation Fields/Latitude'));
lat = nan([1,size(Templat)]);
lat(1,:,:) = Templat;
lat(lat<-1e30)=nan;
    
%     data_2 = h5read(FILE_NAME,'/HDFEOS/SWATHS/O3Profile/Geolocation Fields/Altitude');
%     data_2(data_2<-1e30)=nan;

data_2 = h5read(FILE_NAME,'/HDFEOS/SWATHS/O3Profile/Data Fields/RootMeanSquareErrorOfFit');
data_2(data_2<-1e30)=nan;

data_1 = h5read(FILE_NAME,'/HDFEOS/SWATHS/O3Profile/Data Fields/Temperature');
data_1(data_1<-1e30)=nan;

data = h5read(FILE_NAME,'/HDFEOS/SWATHS/O3Profile/Data Fields/O3');
data(data<-1e30)=nan;

info = h5info(FILE_NAME);
Year = double(info.Groups(1).Groups(1).Groups.Attributes(8).Value);
Day = double(info.Groups(1).Groups(1).Groups.Attributes(7).Value);
Month = double(info.Groups(1).Groups(1).Groups.Attributes(6).Value);
date = datenum([Year,Month,Day]);

lon = double(reshape(lon,[size(lon,1),size(lon,2)*size(lon,3)]));
lat = double(reshape(lat,[size(lat,1),size(lat,2)*size(lat,3)]));

data = double(reshape(data,[size(data,1),size(data,2)*size(data,3)]));
data_1 = double(reshape(data_1,[size(data_1,1),size(data_1,2)*size(data_1,3)]));
data_2 = double(reshape(data_2,[size(data_2,1)*size(data_2,2),1]));

%% reshape
lon = reshape(lon,[size(lon,1)*size(lon,2),1]);
lat = reshape(lat,[size(lat,1)*size(lat,2),1]);

end
