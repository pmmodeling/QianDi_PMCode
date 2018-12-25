%%%% main function processing raw MOD04L2 data, aggregate data from tiles together
% this file need aggregate files placed at the output folder; 
% L2 level data, its lat/lon are changing from day to day, so I used
% interpolation to put them into a unique framework

%% Parameters:
% YEAR: the year of data to be processed
% OPTION: 'MOD04L2': MOD04L2 data; 'OMHCHOG': formaldehyde data (still in
% experimentation)

%% example:
% ReadMOD04L2_Main(2001,'MOD04L2')

%% Code
function ReadMOD04L2_Main(YEAR,OPTION)

%% path 
SEP = '/';
DIRPATH = ['../data/unprocessed/',OPTION,SEP];
OUTPUTPATH = ['../data/aggregate/',OPTION,SEP];

% file path of shape file
EnvironPara = struct();
EnvironPara.states = geoshape(shaperead('usastatehi', 'UseGeoCoords', true));
EnvironPara.S_back = geoshape(shaperead('../data/shapefile/US_WGS_clipcoast_2.shp', 'UseGeoCoords', true));
EnvironPara.US_North_America_Equidistant_Conic = '../data/shapefile/US_North_America_Equidistant_Conic.shp';
EnvironPara.US = '../data/shapefile/US.shp';
GCS = 'North_America_Equidistant_Conic';

% specify the variable of interest for each data set
if(strcmp(OPTION,'OMHCHOG'))
    TEMPLATE = 'OMI-Aura_L2G-OMHCHOG';
    VARIABLE_LIST = {'ColumnAmountDestriped'};
elseif(strcmp(OPTION,'MOD04L2'))
    TEMPLATE = 'MOD04_L2.A';
    VARIABLE_LIST = {'Deep_Blue_Aerosol_Optical_Depth_550_Land'};%% AOD at 550 nm
end


EnvironPara.OUTPUTPATH_PIC = OUTPUTPATH;
SITENAME = ['US',OPTION];
StartDay = datenum([YEAR,1,1]);

list = dir(DIRPATH);
fnames = {list.name}.';
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(fnames, str, 'once'));

% load consistent grid cells
if(exist([OUTPUTPATH,SITENAME,'Site','.mat'],'file'))
    SiteData = LoadData_function([OUTPUTPATH,SITENAME,'Site','.mat']);%SiteData --- uses lat/lon; since raw data are in lat/lon
    SiteData = SiteData.SiteData;
else
    disp('creating site data!!!');
    return;
end



N_Site = size(SiteData,1);
N_Day = yeardays(YEAR);
ResultAll = nan(N_Day,N_Site,length(VARIABLE_LIST));

diary([OUTPUTPATH,'ReadMOD04L2_Main',num2str(YEAR),'.txt']);

%% read data for every day
for i=1:N_Day
    CurrentDay = StartDay+i-1;
    fprintf('%s\n',datestr(CurrentDay,'yyyy-mm-dd'));
    
    % find all files for this day, given a file name template
    if(strcmp(OPTION,'OMHCHOG'))
        FileName_Template = [TEMPLATE,'_',num2str(year(CurrentDay),'%02d'),'m',num2str(month(CurrentDay),'%02d'),num2str(day(CurrentDay),'%02d'),'_v003-','*'];
    elseif(strcmp(OPTION,'MOD04L2'))
        FileName_Template = [TEMPLATE,'',num2str(year(CurrentDay)),num2str(i,'%03d'),'.*'];
    end

    FileName = fnames(FIND(FileName_Template));
    if(isempty(FileName))
        continue;
    end
    
    % read multiple files for a single day
    Temp = nan(N_Site,length(VARIABLE_LIST),length(FileName));
    for j=1:length(FileName)
        fprintf('\t%s\n',FileName{j});
        for k = 1:length(VARIABLE_LIST)
            try
                % read a single file
                [data,SiteData_Temp,date] = ReadMOD04L2_function([DIRPATH,FileName{j}],OPTION,VARIABLE_LIST{k});
                if(isempty(data))
                    disp('option and variable combination error!');
                    return;
                end
                % interploate data from irregulate grid cells to consistent grid cells
                Temp(:,k,j) = InterpMyData_2(data',SiteData_Temp.Lon,SiteData_Temp.Lat,SiteData.Lon,SiteData.Lat,'Threshold0.2');
            catch
                fprintf('\t failed to read:%s\n',FileName{j});
            end
        end
    end
    ResultAll(i,:,:) = nanmean(Temp,3);
%     %visualization
%     if(sum(isnan(ResultAll(i,:,:)))==size(SiteData,1))
%         
%     else
%         Visualization_USResult_1([OPTION,'_',datestr(CurrentDay)],ResultAll(i,:,:),SiteData,EnvironPara);
%     end
end

%% save results
for k = 1:length(VARIABLE_LIST)
    Result = ResultAll(:,:,k);
    Visualization_USResult_1([OPTION,'_',VARIABLE_LIST{k},'_',SITENAME,'_',num2str(YEAR),'_',num2str(YEAR)],nanmean(Result),SiteData,EnvironPara);
    save([OUTPUTPATH,OPTION,'_',VARIABLE_LIST{k},'_',SITENAME,'_',num2str(YEAR),'_',num2str(YEAR),'.mat'],'Result');
end
diary off

end

%% read the messy MOD04L2 data
function [data,SiteData,date] = ReadMOD04L2_function(FILE_NAME,OPTION,VARIABLE)
    
if(strcmp(VARIABLE,'Deep_Blue_Aerosol_Optical_Depth_550_Land'))
%     info = hdfinfo('D:\Google Drive\Research\USTemperature\AODData\MOD04L2\MOD04_L2.A2016149.1940.006.2016152210454.hdf','eos')

    data = double(hdfread(FILE_NAME,'Deep_Blue_Aerosol_Optical_Depth_550_Land'));
%     0       0-25% Cloudy pixels \n",
%     "Quality Flag                      1       25-50% cloudy pixels\n",
%     "                                  2       50-75% cloudy pixels\n",
%     "                                  3       75-100%cloudy pixels\n",
    data_QA = double(hdfread(FILE_NAME,'Deep_Blue_Aerosol_Optical_Depth_550_Land_QA_Flag'));
    %remove cloudy pixel
    data(data_QA>2)=nan;
   
    lon = double(hdfread(FILE_NAME,'Longitude'));
    lat = double(hdfread(FILE_NAME,'Latitude'));
    data(data==-9999)=nan;
    data(data>5000)=nan;
    data(data<0)=nan;
    
    fprintf('nan:%f\n',sum(sum(isnan(data)))/(size(data,1)*size(data,2)));
    %         % testing code visualization
    %         property_state = makesymbolspec('Patch',{'Default', 'LineStyle',':','FaceAlpha',0,'LineWidth',2,'EdgeColor',[0.2,0.2,0.2]});
    %
    %         xmax = -60;
    %         xmin = -130;
    %         ymin = 18;
    %         ymax = 52;
    %         xunit = (xmax - xmin)/250;
    %         yunit = (xmax - xmin)/250;
    %         R = [0,yunit;xunit,0;xmin,ymin];
    %         geoshow(lat,lon,data,'DisplayType', 'texturemap');
    %         colormap jet
    %         geoshow(geoshape(shaperead('usastatehi', 'UseGeoCoords', true)),'SymbolSpec',property_state);
    %         axis([xmin xmax ymin ymax]);
    %         caxis([0 1500]);
    %         colorbar();
    %
    %         PicWidth = 17.35;% 8.3, 12.35, 17.35, cm
    %         PicHeight = PicWidth;
    %         set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 PicWidth PicHeight]);
    %         set(gca,'FontSize',4);
    % %         title(strrep(FILE_NAME,'_',' '),'FontSize',8);
    %         print(gcf,'-dtiff','-r300',[strrep(FILE_NAME,'.hdf',''),'.jpeg']);

else
    data = [];
    SiteData = [];
    date = [];
    return;
end

% reshape the data matrix and create data table
date = nan;
lon = reshape(lon,[size(lon,1)*size(lon,2),1]);
lat = reshape(lat,[size(lat,1)*size(lat,2),1]);
data = reshape(data,[size(data,1)*size(data,2),1]);
SiteCode = int64((1:1:length(lat))');
SiteData = table(SiteCode,lat,lon,'VariableNames',{'SiteCode','Lat','Lon'});

end
