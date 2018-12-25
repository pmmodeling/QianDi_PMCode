%% read OMI aerosol index data
% this function reads different OMI data set: OMAERUVd, OMAEROe, OMNO2d,
% OMSO2e, OMTO3e, OMUVBd and extract variables of interests. The variable
% of interest for each data set is hard coded inside.
% for return values, this function return the data for variable of interests (data), return corresponding
% locations (SiteData) and date of measurements (date)

%% version history
% version1: can read multiple types of files (L3)
% date: datenum format; all outputs are in n*1 format;
% 2016-11-09: add more options to read; coordinate projection system;based on ReadOMIAIData_1

%% example 
% this function is called inside ReadMODIS_Main.m, not supposed to be used separately
% as a stand alone code

%% code
function [data,SiteData,date] = ReadOMIAIData_function(FILE_NAME,OPTION,VARIABLE)
% disp('ReadOMIAIData_function');

% %% for testing only
% % % % Open the HDF5 File.
% FILE_NAME = 'D:\Google Drive\Research\USTemperature\AODData\OMAEROe\OMI-Aura_L3-OMAEROe_2009m0609_v003-2011m1109t114126.he5';
% % FILE_NAME = 'D:\Google Drive\Research\USTemperature\AODData\OMUVBd\OMI-Aura_L3-OMUVBd_2010m0628_v003-2016m0706t221231.he5';
% % FILE_NAME = 'D:\Google Drive\Research\USTemperature\AODData\OMAERUVd\OMI-Aura_L3-OMAERUVd_2005m0410_v003-2012m0325t103235.he5';
% % FILE_NAME = 'D:\Google Drive\Research\USTemperature\AODData\OMSO2e\OMI-Aura_L3-OMSO2e_2016m1011_v003-2016m1013t041205.he5';
% % FILE_NAME = 'D:\Google Drive\Research\USTemperature\AODData\OMNO2d\OMI-Aura_L3-OMNO2d_2016m0819_v003-2016m0827t002655.he5';
% 
% OPTION = 'UVindex';%SolarZenithAngle;ViewingZenithAngle
% h5disp(FILE_NAME);

%% ColumnAmountO3
if(strcmp(OPTION,'OMTO3e'))
    if(strcmp(VARIABLE,'ColumnAmountO3'))
        data = double(h5read(FILE_NAME,'/HDFEOS/GRIDS/OMI Column Amount O3/Data Fields/ColumnAmountO3'));
        data(data<-1e30)=nan;
    elseif(strcmp(VARIABLE,'ViewingZenithAngle'))
        data = double(h5read(FILE_NAME,'/HDFEOS/GRIDS/OMI Column Amount O3/Data Fields/ViewingZenithAngle'));
        data(data<-1e30)=nan;
    end
    %     http://disc.sci.gsfc.nasa.gov/Aura/data-holdings/OMI/documents/v003/OMTO3e_FileSpec_V003.txt
    %     The center of the
    % ! first grid cell is located at longitude -179.875 and latitude -89.875. The
    % ! center of the final grid cell is located at longitude 179.875 and latitude
    % ! 89.875. The center of the grid itself is located at longitude 0.0 and
    % ! latitude 0.0, and corresponds to the corners of four grid cells.
    [lat,lon] = meshgrid(-89.875:0.25:89.875,-179.875:0.25:179.875);%
    
    info = h5info(FILE_NAME);
    Year = double(info.Groups(1).Groups(1).Groups.Attributes(7).Value);
    Day = double(info.Groups(1).Groups(1).Groups.Attributes(6).Value);
    Month = double(info.Groups(1).Groups(1).Groups.Attributes(5).Value);
    date = datenum([Year,Month,Day]);
    
    data = reshape(data,[size(data,1)*size(data,2),1]);
    lon = reshape(lon,[size(lon,1)*size(lon,2),1]);
    lat = reshape(lat,[size(lat,1)*size(lat,2),1]);
    data = reshape(data,[size(data,1)*size(data,2),1]);
    
    SiteCode = int64((1:1:length(lat))');
    SiteData = table(SiteCode,lat,lon,'VariableNames',{'SiteCode','Lat','Lon'});
    
    
elseif(strcmp(OPTION,'OMAERUVd'))
    if(strcmp(VARIABLE,'UVAerosolIndex'))
        data = double(h5read(FILE_NAME,'/HDFEOS/GRIDS/Aerosol NearUV Grid/Data Fields/UVAerosolIndex'));
        data(data<-1e30)=nan;
    elseif(strcmp(VARIABLE,'FinalAerosolSingleScattAlb500'))
        data = double(h5read(FILE_NAME,'/HDFEOS/GRIDS/Aerosol NearUV Grid/Data Fields/FinalAerosolSingleScattAlb500'));
        data(data<-1e30)=nan;
    else
        data = [];
        SiteData = [];
        date = [];
        return;
    end
    
    %% document: http://disc.sci.gsfc.nasa.gov/Aura/data-holdings/OMI/omaeruvd_v003.shtml
    [lat,lon] = meshgrid(-89.5:1:89.5,-179.5:1:179.5);%
    
    info = h5info(FILE_NAME);
    Year = double(info.Groups(1).Groups(1).Groups.Attributes(7).Value);
    Day = double(info.Groups(1).Groups(1).Groups.Attributes(6).Value);
    Month = double(info.Groups(1).Groups(1).Groups.Attributes(5).Value);
    date = datenum([Year,Month,Day]);
    
    lon = reshape(lon,[size(lon,1)*size(lon,2),1]);
    lat = reshape(lat,[size(lat,1)*size(lat,2),1]);
    data = reshape(data,[size(data,1)*size(data,2),1]);
    
    SiteCode = int64((1:1:length(lat))');
    SiteData = table(SiteCode,lat,lon,'VariableNames',{'SiteCode','Lat','Lon'});
    
elseif(strcmp(OPTION,'OMUVBd'))
    if(strcmp(VARIABLE,'UVindex'))
        data = double(h5read(FILE_NAME,'/HDFEOS/GRIDS/OMI UVB Product/Data Fields/UVindex'));
        data(data<-1e30)=nan;
    elseif(strcmp(VARIABLE,'SolarZenithAngle'))
        data = double(h5read(FILE_NAME,'/HDFEOS/GRIDS/OMI UVB Product/Data Fields/SolarZenithAngle'));
        data(data<-1e30)=nan;
    elseif(strcmp(VARIABLE,'ViewingZenithAngle'))
        data = double(h5read(FILE_NAME,'/HDFEOS/GRIDS/OMI UVB Product/Data Fields/ViewingZenithAngle'));
        data(data<-1e30)=nan;
    else
        data = [];
        SiteData = [];
        date = [];
        return;
    end
    
    %% http://disc.sci.gsfc.nasa.gov/Aura/data-holdings/OMI/documents/v003/OMUVBd_README_V003.pdf
    [lat,lon] = meshgrid(-89.5:1:89.5,-179.5:1:179.5);%
    
    info = h5info(FILE_NAME);
    Year = double(info.Groups(1).Groups(1).Groups.Attributes(6).Value);
    Day = double(info.Groups(1).Groups(1).Groups.Attributes(3).Value);
    Month = double(info.Groups(1).Groups(1).Groups.Attributes(5).Value);
    date = datenum([Year,Month,Day]);
    
    lon = reshape(lon,[size(lon,1)*size(lon,2),1]);
    lat = reshape(lat,[size(lat,1)*size(lat,2),1]);
    data = reshape(data,[size(data,1)*size(data,2),1]);
    
    SiteCode = int64((1:1:length(lat))');
    SiteData = table(SiteCode,lat,lon,'VariableNames',{'SiteCode','Lat','Lon'});
    
elseif(strcmp(OPTION,'OMAEROe'))
    if(strcmp(VARIABLE,'VISAerosolIndex'))
        data = double(h5read(FILE_NAME,'/HDFEOS/GRIDS/ColumnAmountAerosol/Data Fields/VISAerosolIndex'));
        data(data<-1e30)=nan;
    elseif(strcmp(VARIABLE,'UVAerosolIndex'))
        data = double(h5read(FILE_NAME,'/HDFEOS/GRIDS/ColumnAmountAerosol/Data Fields/UVAerosolIndex'));
        data(data<-1e30)=nan;
    elseif(strcmp(VARIABLE,'ViewingZenithAngle'))
        data = double(h5read(FILE_NAME,'/HDFEOS/GRIDS/ColumnAmountAerosol/Data Fields/ViewingZenithAngle'));
        data(data<-1e30)=nan;
    elseif(strcmp(VARIABLE,'SolarZenithAngle'))
        data = double(h5read(FILE_NAME,'/HDFEOS/GRIDS/ColumnAmountAerosol/Data Fields/SolarZenithAngle'));
        data(data<-1e30)=nan;
    else
        data = [];
        SiteData = [];
        date = [];
        return;
    end
    
    %% document: http://disc.sci.gsfc.nasa.gov/Aura/data-holdings/OMI/omaeruvd_v003.shtml
    [lat,lon] = meshgrid(-89.875:0.25:89.875,-179.875:0.25:179.875);%
    
    info = h5info(FILE_NAME);
    Year = double(info.Groups(1).Groups(1).Groups.Attributes(7).Value);
    Day = double(info.Groups(1).Groups(1).Groups.Attributes(6).Value);
    Month = double(info.Groups(1).Groups(1).Groups.Attributes(5).Value);
    date = datenum([Year,Month,Day]);
    
    lon = reshape(lon,[size(lon,1)*size(lon,2),1]);
    lat = reshape(lat,[size(lat,1)*size(lat,2),1]);
    data = reshape(data,[size(data,1)*size(data,2),1]);
    
    SiteCode = int64((1:1:length(lat))');
    SiteData = table(SiteCode,lat,lon,'VariableNames',{'SiteCode','Lat','Lon'});
    
elseif(strcmp(OPTION,'OMSO2e'))
    if(strcmp(VARIABLE,'ColumnAmountSO2_PBL'))
        data = double(h5read(FILE_NAME,'/HDFEOS/GRIDS/OMI Total Column Amount SO2/Data Fields/ColumnAmountSO2_PBL'));
        data(data<-1e30)=nan;
    else
        data = [];
        SiteData = [];
        date = [];
        return;
    end
    
    %% document: http://disc.sci.gsfc.nasa.gov/Aura/data-holdings/OMI/omaeruvd_v003.shtml
    [lat,lon] = meshgrid(-89.875:0.25:89.875,-179.875:0.25:179.875);%
    
    info = h5info(FILE_NAME);
    Year = double(info.Groups(1).Groups(1).Groups.Attributes(7).Value);
    Day = double(info.Groups(1).Groups(1).Groups.Attributes(6).Value);
    Month = double(info.Groups(1).Groups(1).Groups.Attributes(5).Value);
    date = datenum([Year,Month,Day]);
    
    lon = reshape(lon,[size(lon,1)*size(lon,2),1]);
    lat = reshape(lat,[size(lat,1)*size(lat,2),1]);
    data = reshape(data,[size(data,1)*size(data,2),1]);
    
    SiteCode = int64((1:1:length(lat))');
    SiteData = table(SiteCode,lat,lon,'VariableNames',{'SiteCode','Lat','Lon'});
    
elseif(strcmp(OPTION,'OMNO2d'))
    %     FILE_NAME = 'D:\Google Drive\Research\USTemperature\AODData\OMNO2d\OMI-Aura_L3-OMNO2d_2016m0828_v003-2016m0829t151955.he5';
    if(strcmp(VARIABLE,'ColumnAmountNO2TropCloudScreened'))
        data = double(h5read(FILE_NAME,'/HDFEOS/GRIDS/ColumnAmountNO2/Data Fields/ColumnAmountNO2TropCloudScreened'));
        data(data<-1e30)=nan;
    elseif(strcmp(VARIABLE,'ColumnAmountNO2CloudScreened'))
        data = double(h5read(FILE_NAME,'/HDFEOS/GRIDS/ColumnAmountNO2/Data Fields/ColumnAmountNO2CloudScreened'));
        data(data<-1e30)=nan;
    elseif(strcmp(VARIABLE,'ColumnAmountNO2StratoCloudScreened'))% NO2 in the stratosphere
        data1 = double(h5read(FILE_NAME,'/HDFEOS/GRIDS/ColumnAmountNO2/Data Fields/ColumnAmountNO2CloudScreened'));
        data1(data1<-1e30)=nan;
        
        data2 = double(h5read(FILE_NAME,'/HDFEOS/GRIDS/ColumnAmountNO2/Data Fields/ColumnAmountNO2TropCloudScreened'));
        data2(data2<-1e30)=nan;
        
        data = data1 - data2;
    else
        data = [];
        SiteData = [];
        date = [];
        return;
    end
    
    %% document: http://disc.sci.gsfc.nasa.gov/Aura/data-holdings/OMI/omaeruvd_v003.shtml
    [lat,lon] = meshgrid(-89.875:0.25:89.875,-179.875:0.25:179.875);%
    
    info = h5info(FILE_NAME);
    Year = double(info.Groups(1).Groups(1).Groups.Attributes(5).Value);
    Day = double(info.Groups(1).Groups(1).Groups.Attributes(7).Value);
    Month = double(info.Groups(1).Groups(1).Groups.Attributes(6).Value);
    date = datenum([Year,Month,Day]);
    
    lon = reshape(lon,[size(lon,1)*size(lon,2),1]);
    lat = reshape(lat,[size(lat,1)*size(lat,2),1]);
    data = reshape(data,[size(data,1)*size(data,2),1]);
    
    SiteCode = int64((1:1:length(lat))');
    SiteData = table(SiteCode,lat,lon,'VariableNames',{'SiteCode','Lat','Lon'});
elseif(strcmp(OPTION,'OMHCHO'))%% not used any more
    if(strcmp(VARIABLE,'ReferenceSectorCorrectedVerticalColumn'))
        data = double(h5read(FILE_NAME,'/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Data Fields/ReferenceSectorCorrectedVerticalColumn'));
        data_QA =  double(h5read(FILE_NAME,'/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Data Fields/MainDataQualityFlag'));
        data(data_QA~=0) = nan;
%         data(data_QA<=-1) = nan;
        data_weight = double(h5read(FILE_NAME,'/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Data Fields/ScatteringWeights'));
        data(data<=-1e30|data>=1e30)=nan;

%         data = double(h5read(FILE_NAME,'/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Geolocation Fields/TerrainHeight'));

    else
        data = [];
        SiteData = [];
        date = [];
        return;
    end
    
    info = h5info(FILE_NAME);
    Year = double(info.Groups(1).Groups(1).Groups.Attributes(13).Value);
    Day = double(info.Groups(1).Groups(1).Groups.Attributes(11).Value);
    Month = double(info.Groups(1).Groups(1).Groups.Attributes(12).Value);
    date = datenum([Year,Month,Day]);
    
%     lon = double(h5read(FILE_NAME,'/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Data Fields/PixelCornerLongitudes'));
%     lon(lon<=-180|lon>=180)=nan;
%     lat = double(h5read(FILE_NAME,'/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Data Fields/PixelCornerLatitudes'));
%     lat(lat<=-90|lat>=90)=nan;

    lon = double(h5read(FILE_NAME,'/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Geolocation Fields/Longitude'));
    lon(lon<=-180|lon>=180)=nan;
    lat = double(h5read(FILE_NAME,'/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Geolocation Fields/Latitude'));
    lat(lat<=-90|lat>=90)=nan;

% %     % TAI 93 unit 1993-1-1 is the reference
% %     %     time = double(h5read(FILE_NAME,'/HDFEOS/GRIDS/OMI Total Column Amount HCHO/Data Fields/Time'));
% %     %     time(time<=-1e30)=nan;
% %     %     time1 = datestr(datenum(1993,1,1)+7.166646144556409e+08/86400);
% %     % visalization
% %     property_state = makesymbolspec('Patch',{'Default', 'LineStyle',':','FaceAlpha',0,'LineWidth',2,'EdgeColor',[0.2,0.2,0.2]});
% %     
%     xmax = 180;
%     xmin = -180;
%     ymin = -90;
%     ymax = 90;
%     
% %     data(~isnan(data))=1;
%     data = data/10^19;
%     geoshow(lat,lon,data,'DisplayType', 'texturemap');
%     colormap jet
% %     geoshow(geoshape(shaperead('usastatehi', 'UseGeoCoords', true)),'SymbolSpec',property_state);
%     axis([xmin xmax ymin ymax]);
%     caxis([min(min(data)) max(max(data))]);
%     colorbar();
% 
%     PicWidth = 17.35;% 8.3, 12.35, 17.35, cm
%     PicHeight = PicWidth;
%     set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 PicWidth PicHeight]);
%     set(gca,'FontSize',4);
% %     title(strrep(FILE_NAME,'_',' '),'FontSize',8);
%     print(gcf,'-dtiff','-r300',[strrep(strrep(FILE_NAME,'.he5',''),'.hdf',''),'_.jpeg']);
%     
    
    
    data = reshape(data,[size(data,1)*size(data,2),1]);
    lon = reshape(lon,[size(lon,1)*size(lon,2),1]);
    lat = reshape(lat,[size(lat,1)*size(lat,2),1]);
    %     time = reshape(time,[size(time,1)*size(time,2),size(time,3)]);
    
    SiteData = table(lat,lon,'VariableNames',{'Lat','Lon'});
    
        % visualization
    data(data<quantile(data,0.05))=nan;
    data(data>quantile(data,0.95))=nan;
    Index = ~isnan(lat)&~isnan(lon)&~isnan(data);
    scatter(lon(Index),lat(Index),0.5,mapminmax(data(Index)',0,1));
    print(gcf,'-dtiff','-r300',[strrep(strrep(FILE_NAME,'.he5',''),'.hdf',''),'_.jpeg']);

elseif(strcmp(OPTION,'MOD04L2'))
    if(strcmp(VARIABLE,'Deep_Blue_Aerosol_Optical_Depth_550_Land'))
        data = double(hdfread(FILE_NAME,'Deep_Blue_Aerosol_Optical_Depth_550_Land'));
        lon = double(hdfread(FILE_NAME,'Longitude'));
        lat = double(hdfread(FILE_NAME,'Latitude'));
        data(data==-9999)=nan;
        data(data>5000)=nan;
        data(data<0)=nan;
        %         % visualization
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
    date = nan;
    lon = reshape(lon,[size(lon,1)*size(lon,2),1]);
    lat = reshape(lat,[size(lat,1)*size(lat,2),1]);
    data = reshape(data,[size(data,1)*size(data,2),1]);
    SiteCode = int64((1:1:length(lat))');
    SiteData = table(SiteCode,lat,lon,'VariableNames',{'SiteCode','Lat','Lon'});
else
    data = [];
    SiteData = [];
    date = [];
    return;
end

% fclose all;
