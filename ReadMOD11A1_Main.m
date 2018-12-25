%%%% main function processing raw MOD11A1 data, aggregate MOD11A1 data from tiles together
% this file need aggregate files placed at the output folder

%% version history
% adapt from: ReadAODData.m
% 2016-07-27: read MOD11A1 data and arrange them by date

%% Parameters:
% YEAR: the year of data to be processed

%% example:
% ReadMOD11A1_Main(2001)

%% Code
function ReadMOD11A1_Main(YEAR)

% path
DIRPATHROOT = '../data/unprocessed/MOD11A1/';
OUTPUTPATH = '../data/aggregate/MOD11A1/';

% file path of shape file
EnvironPara = struct();
EnvironPara.states = geoshape(shaperead('usastatehi', 'UseGeoCoords', true));
EnvironPara.S_back = geoshape(shaperead('../data/shapefile/US_WGS_clipcoast_2.shp', 'UseGeoCoords', true));
EnvironPara.US_North_America_Equidistant_Conic = '../data/shapefile/US_North_America_Equidistant_Conic.shp';
EnvironPara.US = '../data/shapefile/US.shp';
GCS = 'North_America_Equidistant_Conic';


%% create site data
SITENAME = 'USMOD11A1';
DateA = datenum(YEAR,1,1);
DateB = datenum(YEAR,12,31);
NSite = 1200*1200;

if(exist([OUTPUTPATH,SITENAME,'Site.mat'],'file'))
    SiteData = load([OUTPUTPATH,SITENAME,'Site.mat']);%SiteData, TileName --- tile name is order of tiles while creating the SiteData
    TileName = SiteData.TileName;
    SiteData = SiteData.SiteData;
else
    disp('no site data!!!');
    return;
end

%% read individual files for a year and aggregate them into unified grid cells
list = dir(DIRPATHROOT);
fnames = {list.name}.';
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(fnames, str, 'all'));

diary([OUTPUTPATH,'ReadMOD11A1_Main_Diary.txt']);
for i = DateA:DateB
    CurrentDay = i;
    
    fprintf('%s\n',datestr(CurrentDay));
    
    % get julian date and find all files matching file name pattern
    JulianDay = CurrentDay - datenum(year(CurrentDay),1,1)+1;
    FileName_Template = ['MOD11A1.A',num2str(year(CurrentDay)), sprintf('%03d',JulianDay),'.*'];
    FileName = fnames(FIND(FileName_Template));
    
    % to store output
    Data_LST_Day_1km = nan(size(SiteData,1),1);
    Data_QC_Day = nan(size(SiteData,1),1);
    Data_LST_Night_1km = nan(size(SiteData,1),1);
    Data_QC_Night = nan(size(SiteData,1),1);
    Data_Emis_31 = nan(size(SiteData,1),1);
    Data_Emis_32 = nan(size(SiteData,1),1);
    Data_Clear_day_cov = nan(size(SiteData,1),1);
    Data_Clear_night_cov = nan(size(SiteData,1),1);
    TempOutputFileName = [OUTPUTPATH,'MOD11A1_',SITENAME,'_',datestr(CurrentDay,'yyyymmdd'),'_',datestr(CurrentDay,'yyyymmdd'),'.mat'];
    
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
            
            % get tile
            Index = strfind(TempFileName,'.');
            TileNameTemp =  TempFileName(Index(2)+1:Index(3)-1);
            if(~ismember(TileNameTemp,TileName))
                continue;
            end
            
            % extract variables of interest and match them to unified grid
            % cells
            [~,locb]=ismember(TileNameTemp,TileName);
            fprintf('\t%s\n',TempFileName);
            FILE_NAME = [DIRPATHROOT,TempFileName];
            try
                % read the data and reshape the matrix
                Data_LST_Day_1km_Temp = 0.02*double(hdfread(FILE_NAME,'LST_Day_1km'));
                Data_QC_Day_Temp = double(swapbytes(hdfread(FILE_NAME,'QC_Day')));%dec2bin(swapbytes(uint8(25)),8)
                Data_LST_Night_1km_Temp = 0.02*double(hdfread(FILE_NAME,'LST_Night_1km'));
                Data_QC_Night_Temp = double(swapbytes(hdfread(FILE_NAME,'QC_Night')));
                Data_Emis_31_Temp = 0.002*double(hdfread(FILE_NAME,'Emis_31'));
                Data_Emis_32_Temp = 0.002*double(hdfread(FILE_NAME,'Emis_32'));
                Data_Clear_day_cov_Temp = 0.0001*double(hdfread(FILE_NAME,'Clear_day_cov'));%scaling factor should be 0.0005; but I made a mistake when processing the data
                Data_Clear_night_cov_Temp = 0.0001*double(hdfread(FILE_NAME,'Clear_night_cov'));
            catch exception
                fprintf('%s\n',exception.message);
                break;
            end
            
            % reshape the data matrix
            Data_LST_Day_1km_Temp = reshape(Data_LST_Day_1km_Temp,[NSite,1]);
            Data_QC_Day_Temp = (reshape(Data_QC_Day_Temp,[NSite,1]));
            Data_LST_Night_1km_Temp = reshape(Data_LST_Night_1km_Temp,[NSite,1]);
            Data_QC_Night_Temp = (reshape(Data_QC_Night_Temp,[NSite,1]));
            Data_Emis_31_Temp = reshape(Data_Emis_31_Temp,[NSite,1]);
            Data_Emis_32_Temp = reshape(Data_Emis_32_Temp,[NSite,1]);
            Data_Clear_day_cov_Temp = reshape(Data_Clear_day_cov_Temp,[NSite,1]);
            Data_Clear_night_cov_Temp = reshape(Data_Clear_night_cov_Temp,[NSite,1]);
            
            % remove invalid values
            Data_LST_Day_1km_Temp(Data_LST_Day_1km_Temp==0)=nan;
            Data_LST_Night_1km_Temp(Data_LST_Night_1km_Temp==0)=nan;
            Data_Emis_31_Temp(Data_Emis_31_Temp==0)=nan;
            Data_Emis_32_Temp(Data_Emis_32_Temp==0)=nan;
            Data_Clear_day_cov_Temp(Data_Clear_day_cov_Temp==0)=nan;
            Data_Clear_night_cov_Temp(Data_Clear_night_cov_Temp==0)=nan;
            
            % match extract data to the consistent grid cells
            Data_LST_Day_1km(((locb-1)*NSite+1):(locb*NSite))=Data_LST_Day_1km_Temp;
            Data_QC_Day(((locb-1)*NSite+1):(locb*NSite))=Data_QC_Day_Temp;
            Data_LST_Night_1km(((locb-1)*NSite+1):(locb*NSite))=Data_LST_Night_1km_Temp;
            Data_QC_Night(((locb-1)*NSite+1):(locb*NSite))=Data_QC_Night_Temp;
            Data_Emis_31(((locb-1)*NSite+1):(locb*NSite))=Data_Emis_31_Temp;
            Data_Emis_32(((locb-1)*NSite+1):(locb*NSite))=Data_Emis_32_Temp;
            Data_Clear_day_cov(((locb-1)*NSite+1):(locb*NSite))=Data_Clear_day_cov_Temp;
            Data_Clear_night_cov(((locb-1)*NSite+1):(locb*NSite))=Data_Clear_night_cov_Temp;
        end
        
        Data_LST_Day_1km = Data_LST_Day_1km';
        Data_QC_Day = Data_QC_Day';
        Data_LST_Night_1km = Data_LST_Night_1km';
        Data_QC_Night = Data_QC_Night';
        Data_Emis_31 = Data_Emis_31';
        Data_Emis_32 = Data_Emis_32';
        Data_Clear_day_cov = Data_Clear_day_cov';
        Data_Clear_night_cov = Data_Clear_night_cov';
        
        save(TempOutputFileName,'Data_LST_Day_1km','Data_QC_Day','Data_LST_Night_1km','Data_QC_Night',...
            'Data_Emis_31','Data_Emis_32','Data_Clear_day_cov','Data_Clear_night_cov');
    end
end
diary off





% %% get SiteData
% DirPath = 'D:\Data\Data_Temperature\MOD11A1\';
% list = dir(DirPath);
% SiteData = [];
% TileName = [];
% for i=3:length(list)
%     % get tile name
%     Index = strfind(list(i).name,'.');
%     TileNameTemp =  list(i).name(Index(2)+1:Index(3)-1);
%     % get siteData
%     FileInfo = hdfinfo([DirPath,list(i).name],'eos');
%     SiteDataTemp = ReadMODIS_function1(FileInfo);
%     SiteData = cat(1,SiteData,SiteDataTemp);
%     TileName = cat(1,TileName,{TileNameTemp});
%
%     fprintf('%s\t%s\n',TileNameTemp,list(i).name);
% end
% save('D:\Data\Data_Temperature\MOD11A1\USTemperatureSite.mat','SiteData','TileName');
