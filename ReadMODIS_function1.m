%% read MODIS data and extract value of interests, 
% lat/lon of grid cells; and measurement date

%% version history
% version1: can read multiple types of files (L3)

%% parameters:
% YEAR = 2000;% year of files to be read
% OPTION = 'Aqua'% Aqua or TERRA

%% example 
% this function is called inside ReadMODIS_Main.m, not supposed to be used separately
% as a stand alone code

%% used in ReadMODIS_Main;
function [data,SiteData,date] = ReadMODIS_function1(FILE_NAME,OPTION)

data = nan;
date = nan;


if(strcmp(OPTION,'MOD09A1'))
    %% surface reflectance
    data = double(hdfread(FILE_NAME,'sur_refl_b03'));
    data = reshape(data,[size(data,1)*size(data,2),1]);
    FileInfo = hdfinfo(FILE_NAME,'eos');
    SiteData = ReadMODIS_function2(FileInfo);
    
elseif(strcmp(OPTION,'MOD13A2'))
    %% surface reflectance
    data = double(hdfread(FILE_NAME,'1 km 16 days NDVI'));
    data = reshape(data,[size(data,1)*size(data,2),1]);
    FileInfo = hdfinfo(FILE_NAME,'eos');
    SiteData = ReadMODIS_function2(FileInfo);   
    
elseif(strcmp(OPTION,'MCD12Q1'))
    %% surface reflectance
    data = double(hdfread(FILE_NAME,'Land_Cover_Type_1'));
    data = reshape(data,[size(data,1)*size(data,2),1]);
    FileInfo = hdfinfo(FILE_NAME,'eos');
    SiteData = ReadMODIS_function2(FileInfo);   
    
end



%% extract the lat/lon of MODIS data; based on the tile number and file header
% some MODIS file does not explicitly store lat/lon of every grid cells; we
% need to calculate by ourselves; the algorithm is specified on this
% website page: https://modis-land.gsfc.nasa.gov/MODLAND_grid.html

%% version history
% 2016-10-11, do not do coordinate transform here, lat and lon are just in meters
% 2016-11-11: change name; based on ReadMODIS_function1

%% parameters:
% FileInfo: file header of a MODIS raw data file

%% example 
% this function is called inside ReadMODIS_Main.m, not supposed to be used separately
% as a stand alone code

%% code
function SiteData = ReadMODIS_function2(FileInfo)

%% version of 2016-10-11 --- SiteData is table and use MODIS projection system
Index = regexp(FileInfo.Filename,'h[\d]{2}');
h = FileInfo.Filename(Index(1)+1:Index(1)+2);
Index = regexp(FileInfo.Filename,'v[\d]{2}');
v = FileInfo.Filename(Index(1)+1:Index(1)+2);

[Lon,Lat] = meshgrid(linspace(FileInfo.Grid.UpperLeft(1),FileInfo.Grid.LowerRight(1),FileInfo.Grid.Rows),...
    linspace(FileInfo.Grid.UpperLeft(2),FileInfo.Grid.LowerRight(2),FileInfo.Grid.Columns));
Lat = reshape(Lat,[size(Lat,1)*size(Lat,2),1]);
Lon = reshape(Lon,[size(Lon,1)*size(Lon,2),1]);
SiteCode = int64(repmat(str2double(['2',h,v])*10000000,[length(Lat),1])+(1:1:length(Lat))');
SiteData = table(SiteCode,Lat,Lon,'VariableNames',{'SiteCode','Lat','Lon'});

% % % % % version proir to 2016-10-11 -- SiteData is matrix and use lat lon (GCS)
% R = FileInfo.Grid.Projection.ProjParam(1);
% N = FileInfo.Grid.Rows;
% Diff = FileInfo.Grid.LowerRight(1)-FileInfo.Grid.UpperLeft(1);
% 
% Lat = nan(N,N);
% Lon = nan(N,N);
% 
% for i=1:N
%    for j=1:N 
% %        Lat(i,j) = (FileInfo.Grid.UpperLeft(2)-(i-1+0.5)*Diff/N)*360/2/pi/R;
% %        Lon(i,j) = (FileInfo.Grid.UpperLeft(1)+(j-1+0.5)*Diff/N)*360/2/pi/(R*cos(Lat(i,j)*pi/180));
%        [Lat(i,j),Lon(i,j)]=GetLatLon(FileInfo.Grid.UpperLeft(2),FileInfo.Grid.UpperLeft(1),R,i,j,Diff,N);
%    end
% end
% 
% Lat = reshape(Lat,[size(Lat,1)*size(Lat,2),1]);
% Lon = reshape(Lon,[size(Lon,1)*size(Lon,2),1]);
% 
% %% get horizontal and vertical number
% Index = regexp(FileInfo.Filename,'h[\d]{2}');
% h = FileInfo.Filename(Index(1)+1:Index(1)+2);
% Index = regexp(FileInfo.Filename,'v[\d]{2}');
% v = FileInfo.Filename(Index(1)+1:Index(1)+2);
% SiteCode = repmat(str2double(['2',h,'3',v])*10000000,[length(Lat),1])+(1:1:length(Lat))';
% SiteData = [SiteCode,Lat,Lon];
% 
% %% i: row; j: column
% function [Lat,Lon] = GetLatLon(Upper,Left,R,i,j,Diff,N)
%     Lat = (Upper-(i-1+0.5)*Diff/N)*360/2/pi/R;
%     Lon = (Left+(j-1+0.5)*Diff/N)*360/2/pi/(R*cos(Lat*pi/180));
% end
% 
% end
