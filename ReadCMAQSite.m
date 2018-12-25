%% read CMAQ outputs by year; the file name of CMAQ outputs is different by
% year; even one year has two types of outputs.
% the output is annual output at unified points, organized by variables of
% interests; also I put CMAQ outputs in different places, the path is also
% hard coded inside this function

%% parameter;
% YEAR: year of CMAQ simulation we will process;
% OPTION: for some years, there are two versions of CMAQ output; I denoted
% them as A and B

%% example:
% ReadCMAQSite(2002,'A');
% ReadCMAQSite(2002,'B');
% ReadCMAQSite(2003,'A');
% ReadCMAQSite(2003,'B');
% ReadCMAQSite(2004,'A');
% ReadCMAQSite(2004,'B');
% ReadCMAQSite(2005,'A');
% ReadCMAQSite(2005,'B');
% ReadCMAQSite(2006,'A');
% ReadCMAQSite(2006,'B');
% ReadCMAQSite(2007,'');
% ReadCMAQSite(2008,'');
% ReadCMAQSite(2009,'');
% ReadCMAQSite(2010,'');
% ReadCMAQSite(2011,'');
% ReadCMAQSite(2012,'');
% ReadCMAQSite(2013,'');

%% code
function ReadCMAQSite(YEAR,OPTION)


%% marco definition
% specify the file name pattern for each year
SITENAME = ['CMAQ',num2str(YEAR),OPTION];
SITENAME_UNI = 'CMAQ';
DIRPATHROOT = '../data/unprocessed/'; 
if(2000 == YEAR)
    TEMPLATE = 'output_us_doe_sf_rrtmg_20_5_1_v3450_2000_36km.36US.35L.cmaq.conc.$Date$';
elseif(2001 == YEAR)
    TEMPLATE = 'output_us_doe_nf_rrtmg_20_5_1_v3450_2001_36km.36US.35L.cmaq.conc.$Date$';
elseif(2002 == YEAR)
    if(strcmp(OPTION,'A'))
        TEMPLATE = 'cctm_n1a_2002af_med_v32soa_v3.4beta3_newbiog.12eus1.24.cmaq.conc.$Date$';  
    elseif(strcmp(OPTION,'B'))
        TEMPLATE = 'cctm_n1a_2002af_med_v32soa_v3.4beta3_newbiog.36us1.24.cmaq.conc.$Date$'; 
    end
elseif(2003 == YEAR)
    if(strcmp(OPTION,'A'))
        TEMPLATE = 'cctm_m2g_v14soa_v3.4beta3_2003.12eus1.24.cmaq.conc.$Date$';
    elseif(strcmp(OPTION,'B'))
        TEMPLATE = 'cctm_m2f_v14soa_v3.4beta3_2003.36us1.24.cmaq.conc.$Date$';
    end   
elseif(2004 == YEAR)
    if(strcmp(OPTION,'A'))
        TEMPLATE = 'cctm_m2f_v14soa_v3.4beta3_2004.36us1.24.cmaq.conc.$Date$';
    elseif(strcmp(OPTION,'B'))
        TEMPLATE = 'cctm_m2g_v14soa_v3.4beta3_2004.12eus1.24.cmaq.conc.$Date$';
    end     
elseif(2005 == YEAR)
    if(strcmp(OPTION,'A'))
        TEMPLATE = 'cctm_n1a_2005af_05b_v3.4beta3.12eus1.24.cmaq.conc.$Date$';
    elseif(strcmp(OPTION,'B'))
        TEMPLATE = 'cctm_n1a_2005af_05b_v3.4beta3.36us1.24.cmaq.conc.$Date$';
    end  
elseif(2006 == YEAR)
    if(strcmp(OPTION,'A'))
        TEMPLATE = 'cctm_n1a_2006af_06b2_med_v3.4beta3.12eus1.24.cmaq.conc.$Date$';
    elseif(strcmp(OPTION,'B'))
        TEMPLATE = 'cctm_n1a_2006af_06b2_med_v3.4beta3.36us1.24.cmaq.conc.$Date$';
    end
elseif(2007 == YEAR)
    TEMPLATE = '2007aq_07c_n5ao_inline.12us1.24.cmaq.conc.$Date$';
elseif(2008 == YEAR)
    TEMPLATE = '2008aa_08c_n5ao_inline.12us1.24.cmaq.conc.$Date$';
elseif(2009 == YEAR)
    TEMPLATE = '2009ef2_09d_cmaq471_n5ao.12us1.24.cmaq.conc.$Date$';
elseif(2010 == YEAR)
    TEMPLATE = '$Date$_harvard.ncf';
elseif(2011 == YEAR)
    TEMPLATE = 'CCTM_v5.0.2_sol_Linux2_x86_64intel_cb05_2011el.12US2.35L.cmaq.conc.$Date$';
elseif(2012 == YEAR)
    TEMPLATE = '2012eh_cb05v2_v6_12g.12US2.25L.cmaq.conc.$Date$';
elseif(2013 == YEAR)
    TEMPLATE = 'CCTM_CMAQv5.1_CCTM_06Nov2015_v51_06Nov2015_Linux2_x86_64intel_.12US1.35L.cmaq.conc.$Date$';
elseif(2014 == YEAR)
    TEMPLATE = 'WRFV381_CMAQ_DEVv5_2Gamma_14Mar2017_cb6r3_ae6nvPOA_aq_2014fb.12US1.35L.cmaq.conc.$Date$';
end

%% path
DIRPATH = [DIRPATHROOT,'CMAQ/CMAQ_',num2str(YEAR),OPTION,'/'];
OUTPUTPATH = '../data/aggregate/CMAQ/';

% file path of the shape file used for making maps
EnvironPara = struct();
EnvironPara.states = geoshape(shaperead('usastatehi', 'UseGeoCoords', true));
EnvironPara.S_back = geoshape(shaperead('../data/shapefile/US_WGS_clipcoast_2.shp', 'UseGeoCoords', true));
EnvironPara.US_North_America_Equidistant_Conic = '../data/shapefile/US_North_America_Equidistant_Conic.shp';
EnvironPara.US = '../data/shapefile/US.shp';
GCS = 'North_America_Equidistant_Conic';

mkdir(OUTPUTPATH);
diary([OUTPUTPATH,'ReadCMAQSite_',num2str(YEAR),'.txt']);
EnvironPara.OUTPUTPATH_PIC = OUTPUTPATH;
fprintf('TEMPLATE:%s\n',TEMPLATE);

%% extract lat/lon of the simulation outputs
info = ncinfo([DIRPATH,strrep(TEMPLATE,'$Date$',datestr(datenum(YEAR,7,1),'yyyymmdd'))]);
NCOLS = info.Attributes(find(strcmp({info.Attributes.Name},'NCOLS'))).Value;
NROWS = info.Attributes(find(strcmp({info.Attributes.Name},'NROWS'))).Value;
NLAYS = info.Attributes(find(strcmp({info.Attributes.Name},'NLAYS'))).Value;
XORIG = info.Attributes(find(strcmp({info.Attributes.Name},'XORIG'))).Value;
YORIG = info.Attributes(find(strcmp({info.Attributes.Name},'YORIG'))).Value;
XCELL = info.Attributes(find(strcmp({info.Attributes.Name},'XCELL'))).Value;
YCELL = info.Attributes(find(strcmp({info.Attributes.Name},'YCELL'))).Value;

%% create site data for current grid
if(exist([OUTPUTPATH,SITENAME,'Site','.mat'],'file'))
    SiteData = LoadData_function([OUTPUTPATH,SITENAME,'Site','.mat']);%SiteData --- uses lat/lon; since raw data are in lat/lon
    SiteData_Temp = SiteData.SiteData;
    Lon_Query = SiteData_Temp.Lon;
    Lat_Query = SiteData_Temp.Lat;
else
%     disp('creating site data!!!');
    [Y,X] = meshgrid((YORIG+YCELL/2):YCELL:(YORIG+YCELL*NROWS-YCELL/2),(XORIG+XCELL/2):XCELL:(XORIG+XCELL*NCOLS-XCELL/2));
    X = double(reshape(X,[NCOLS*NROWS,1]));
    Y = double(reshape(Y,[NCOLS*NROWS,1]));
    SiteData = table((1:1:NCOLS*NROWS)',Y,X,'VariableNames',{'SiteCode','Lat','Lon'});
    save([OUTPUTPATH,SITENAME,'Site','.mat'],'SiteData');
    writetable(SiteData,[OUTPUTPATH,SITENAME,'Site','.csv']);
    
    SiteData_Temp = SiteData;
    Lon_Query = SiteData_Temp.Lon;
    Lat_Query = SiteData_Temp.Lat;
end

%% create/read universal grid cells; --- also the grid cell to be interpolated to
if(exist([OUTPUTPATH,SITENAME_UNI,'Site','.mat'],'file'))
    SiteData = LoadData_function([OUTPUTPATH,SITENAME_UNI,'Site','.mat']);%SiteData --- uses lat/lon; since raw data are in lat/lon
    SiteData_uni = SiteData.SiteData;
    Lon_Target = SiteData_uni.Lon;
    Lat_Target = SiteData_uni.Lat;
else
    disp('creating universal grid cells data!!!');
    return;
end

Result_PM25_EC = nan(yeardays(YEAR),size(SiteData_uni,1));
Result_PM25_NH4 = nan(yeardays(YEAR),size(SiteData_uni,1));
Result_PM25_NO3 = nan(yeardays(YEAR),size(SiteData_uni,1));
Result_PM25_OC = nan(yeardays(YEAR),size(SiteData_uni,1));
Result_PM25_OM = nan(yeardays(YEAR),size(SiteData_uni,1));
Result_PM25_SO4 = nan(yeardays(YEAR),size(SiteData_uni,1));
Result_PM25_TOT = nan(yeardays(YEAR),size(SiteData_uni,1));
Result_PM25_Vertical = nan(yeardays(YEAR),size(SiteData_uni,1));
Result_NOX = nan(yeardays(YEAR),size(SiteData_uni,1));
Result_NO2 = nan(yeardays(YEAR),size(SiteData_uni,1));
Result_NO2_Vertical = nan(yeardays(YEAR),size(SiteData_uni,1));
Result_Ozone = nan(yeardays(YEAR),size(SiteData_uni,1));
Result_Ozone_Vertical = nan(yeardays(YEAR),size(SiteData_uni,1));
Result_RH = nan(yeardays(YEAR),size(SiteData_uni,1));
Result_TA = nan(yeardays(YEAR),size(SiteData_uni,1));
Result_AIR_DENS = nan(yeardays(YEAR),size(SiteData_uni,1));

%% read data
% extract the variables of interest
StartDay = datenum(YEAR,1,1);
PreviousDayOzone_Sum = nan([NCOLS*NROWS,7]);
PreviousDayOzone_Ground = nan([NCOLS*NROWS,7]);
for i = 1:yeardays(YEAR)
   tic;
   CurrentDay = StartDay + i - 1;
   FileName = strrep(TEMPLATE,'$Date$',datestr(CurrentDay,'yyyymmdd'));
   fprintf('processing...%s\n',FileName);
   
   %% ground level of primary air pollutant
   % PM2.5; Temp size: 396   246    25    24
   try
       Temp = nanmean(double(ncread([DIRPATH,FileName],'PM25_TOT')),4);
       Temp_Ground = reshape(Temp(:,:,1),[1,NCOLS*NROWS]);
       Temp_Sum = reshape(nansum(Temp,3),[1,NCOLS*NROWS]);
       Result_PM25_TOT(i,:) = InterpMyData_2(Temp_Ground,Lon_Query,Lat_Query,Lon_Target,Lat_Target,'');
       Result_PM25_Vertical(i,:) = InterpMyData_2(Temp_Ground./Temp_Sum,Lon_Query,Lat_Query,Lon_Target,Lat_Target,'');

       % NO2
       Temp = double(ncread([DIRPATH,FileName],'NO2'));% size(Temp) = 396   246    25    24
       Temp = reshape(Temp,[NCOLS*NROWS,NLAYS,24]);
       % the hourly maximal at ground level;Temp_Max: daily 1-hour max NO2;
       [Temp_Ground,Index] = nanmax(Temp(:,1,:),[],3); 
       Temp_Sum = reshape(nansum(Temp,2),[NCOLS*NROWS,24]);
       Temp_Sum = Temp_Sum(sub2ind(size(Temp_Sum),[1:length(Index)]',Index));
       % the hourly maximal at ground level;
       Result_NO2(i,:) = InterpMyData_2(Temp_Ground',Lon_Query,Lat_Query,Lon_Target,Lat_Target,'');
       Result_NO2_Vertical(i,:) = InterpMyData_2((Temp_Ground./Temp_Sum)',Lon_Query,Lat_Query,Lon_Target,Lat_Target,'');

       % ozone!!!??!!
       Temp = double(ncread([DIRPATH,FileName],'O3'));% Temp: 396   246    25    24
       Temp = reshape(Temp,[NCOLS*NROWS,NLAYS,24]);%97416          25          24
       Temp_Ground = reshape(Temp(:,1,:),[NCOLS*NROWS,24]);% take ground level ozone, Temp_Ground: 97416*24;
       Temp_Ground = cat(2,PreviousDayOzone_Ground,Temp_Ground);% Temp_Ground: 97416*31;
       PreviousDayOzone_Ground = Temp_Ground(:,25:31);% record the last 7 hour 1-hour ozone level: to calculate 8-hour max ozone for the next day
       Temp_Ground = tsmovavg(Temp_Ground,'s',8,2);
       Temp_Ground = Temp_Ground(:,8:31);% now it is 8-hour average ozone for the day; size(Temp_Ground)=97416,24
       [Temp_Ground,Index] = nanmax(Temp_Ground,[],2);%size(Temp_Ground) 97416           1
       % the 8h max ozoneL 
       Result_Ozone(i,:) = InterpMyData_2(Temp_Ground',Lon_Query,Lat_Query,Lon_Target,Lat_Target,'');

       % now it is vertical profile of ozone
       Temp_Sum = reshape(nansum(Temp,2),[NCOLS*NROWS,24]);%size(Temp_Sum) 97416          24
       Temp_Sum = cat(2,PreviousDayOzone_Sum,Temp_Sum);
       PreviousDayOzone_Sum = Temp_Sum(:,25:31);
       Temp_Sum = tsmovavg(Temp_Sum,'s',8,2);
       Temp_Sum = Temp_Sum(:,8:31);
       Temp_Sum = Temp_Sum(sub2ind(size(Temp_Sum),[1:length(Index)]',Index));
       Result_Ozone_Vertical(i,:) = InterpMyData_2((Temp_Ground./Temp_Sum)',Lon_Query,Lat_Query,Lon_Target,Lat_Target,'');

       %% ground level for other variables
       Temp = nanmean(double(ncread([DIRPATH,FileName],'PM25_EC')),4);
       Temp = reshape(Temp(:,:,1),[1,NCOLS*NROWS]);
       Result_PM25_EC(i,:) = InterpMyData_2(Temp,Lon_Query,Lat_Query,Lon_Target,Lat_Target,'');

       Temp = nanmean(double(ncread([DIRPATH,FileName],'PM25_NH4')),4);
       Temp = reshape(Temp(:,:,1),[1,NCOLS*NROWS]);
       Result_PM25_NH4(i,:) = InterpMyData_2(Temp,Lon_Query,Lat_Query,Lon_Target,Lat_Target,'');

       Temp = nanmean(double(ncread([DIRPATH,FileName],'PM25_NO3')),4);
       Temp = reshape(Temp(:,:,1),[1,NCOLS*NROWS]);
       Result_PM25_NO3(i,:) = InterpMyData_2(Temp,Lon_Query,Lat_Query,Lon_Target,Lat_Target,'');

       Temp = nanmean(double(ncread([DIRPATH,FileName],'PM25_OC')),4);
       Temp = reshape(Temp(:,:,1),[1,NCOLS*NROWS]);
       Result_PM25_OC(i,:) = InterpMyData_2(Temp,Lon_Query,Lat_Query,Lon_Target,Lat_Target,'');

       Temp = nanmean(double(ncread([DIRPATH,FileName],'PM25_OM')),4);
       Temp = reshape(Temp(:,:,1),[1,NCOLS*NROWS]);
       Result_PM25_OM(i,:) = InterpMyData_2(Temp,Lon_Query,Lat_Query,Lon_Target,Lat_Target,'');

       Temp = nanmean(double(ncread([DIRPATH,FileName],'PM25_SO4')),4);
       Temp = reshape(Temp(:,:,1),[1,NCOLS*NROWS]);
       Result_PM25_SO4(i,:) = InterpMyData_2(Temp,Lon_Query,Lat_Query,Lon_Target,Lat_Target,'');

       Temp = nanmean(double(ncread([DIRPATH,FileName],'NOX')),4);
       Temp = reshape(Temp(:,:,1),[1,NCOLS*NROWS]);
       Result_NOX(i,:) = InterpMyData_2(Temp,Lon_Query,Lat_Query,Lon_Target,Lat_Target,'');

       Temp = nanmean(double(ncread([DIRPATH,FileName],'RH')),4);
       Temp = reshape(Temp(:,:,1),[1,NCOLS*NROWS]);
       Result_RH(i,:) = InterpMyData_2(Temp,Lon_Query,Lat_Query,Lon_Target,Lat_Target,'');

       Temp = nanmean(double(ncread([DIRPATH,FileName],'TA')),4);
       Temp = reshape(Temp(:,:,1),[1,NCOLS*NROWS]);
       Result_TA(i,:) = InterpMyData_2(Temp,Lon_Query,Lat_Query,Lon_Target,Lat_Target,'');

       Temp = nanmean(double(ncread([DIRPATH,FileName],'AIR_DENS')),4);
       Temp = reshape(Temp(:,:,1),[1,NCOLS*NROWS]);
       Result_AIR_DENS(i,:) = InterpMyData_2(Temp,Lon_Query,Lat_Query,Lon_Target,Lat_Target,'');
   catch exception
       fprintf('%s\n',exception.message);
   end
   toc;
end

if(exist([OUTPUTPATH,SITENAME_UNI,'Site_',GCS,'.mat'],'file'))
    SiteData = LoadData_function([OUTPUTPATH,SITENAME_UNI,'Site_',GCS,'.mat']);%SiteData --- uses lat/lon; since raw data are in lat/lon
    SiteData_uni_map = SiteData.SiteData;
else
    disp('not creating universal North America Equidistant Conic grid cells data!!!');
    return;
end

% make map and store results
save([OUTPUTPATH,'CMAQ_PM25_EC_',SITENAME_UNI,'_',num2str(YEAR),OPTION,'_',num2str(YEAR),OPTION,'.mat'],'Result_PM25_EC');
Visualization_USResult_1(strrep(TEMPLATE,'$Date$',['Result_PM25_EC',datestr(datenum(YEAR,1,1),'yyyymmdd'),OPTION,'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),OPTION]),nanmean(Result_PM25_EC),SiteData_uni_map,EnvironPara);

save([OUTPUTPATH,'CMAQ_PM25_NH4_',SITENAME_UNI,'_',num2str(YEAR),OPTION,'_',num2str(YEAR),OPTION,'.mat'],'Result_PM25_NH4');
Visualization_USResult_1(strrep(TEMPLATE,'$Date$',['Result_PM25_NH4',datestr(datenum(YEAR,1,1),'yyyymmdd'),OPTION,'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),OPTION]),nanmean(Result_PM25_NH4),SiteData_uni_map,EnvironPara);

save([OUTPUTPATH,'CMAQ_PM25_NO3_',SITENAME_UNI,'_',num2str(YEAR),OPTION,'_',num2str(YEAR),OPTION,'.mat'],'Result_PM25_NO3');
Visualization_USResult_1(strrep(TEMPLATE,'$Date$',['Result_PM25_NO3',datestr(datenum(YEAR,1,1),'yyyymmdd'),OPTION,'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),OPTION]),nanmean(Result_PM25_NO3),SiteData_uni_map,EnvironPara);

save([OUTPUTPATH,'CMAQ_PM25_OC_',SITENAME_UNI,'_',num2str(YEAR),OPTION,'_',num2str(YEAR),OPTION,'.mat'],'Result_PM25_OC');
Visualization_USResult_1(strrep(TEMPLATE,'$Date$',['Result_PM25_OC',datestr(datenum(YEAR,1,1),'yyyymmdd'),OPTION,'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),OPTION]),nanmean(Result_PM25_OC),SiteData_uni_map,EnvironPara);

save([OUTPUTPATH,'CMAQ_PM25_OM_',SITENAME_UNI,'_',num2str(YEAR),OPTION,'_',num2str(YEAR),OPTION,'.mat'],'Result_PM25_OM');
Visualization_USResult_1(strrep(TEMPLATE,'$Date$',['Result_PM25_OM',datestr(datenum(YEAR,1,1),'yyyymmdd'),OPTION,'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),OPTION]),nanmean(Result_PM25_OM),SiteData_uni_map,EnvironPara);

save([OUTPUTPATH,'CMAQ_PM25_SO4_',SITENAME_UNI,'_',num2str(YEAR),OPTION,'_',num2str(YEAR),OPTION,'.mat'],'Result_PM25_SO4');
Visualization_USResult_1(strrep(TEMPLATE,'$Date$',['Result_PM25_SO4',datestr(datenum(YEAR,1,1),'yyyymmdd'),OPTION,'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),OPTION]),nanmean(Result_PM25_SO4),SiteData_uni_map,EnvironPara);

save([OUTPUTPATH,'CMAQ_PM25_TOT_',SITENAME_UNI,'_',num2str(YEAR),OPTION,'_',num2str(YEAR),OPTION,'.mat'],'Result_PM25_TOT');
Visualization_USResult_1(strrep(TEMPLATE,'$Date$',['Result_PM25_TOT',datestr(datenum(YEAR,1,1),'yyyymmdd'),OPTION,'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),OPTION]),nanmean(Result_PM25_TOT),SiteData_uni_map,EnvironPara);

save([OUTPUTPATH,'CMAQ_PM25_Vertical_',SITENAME_UNI,'_',num2str(YEAR),OPTION,'_',num2str(YEAR),OPTION,'.mat'],'Result_PM25_Vertical');
Visualization_USResult_1(strrep(TEMPLATE,'$Date$',['Result_PM25_Vertical',datestr(datenum(YEAR,1,1),'yyyymmdd'),OPTION,'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),OPTION]),nanmean(Result_PM25_Vertical),SiteData_uni_map,EnvironPara);

save([OUTPUTPATH,'CMAQ_NOX_',SITENAME_UNI,'_',num2str(YEAR),OPTION,'_',num2str(YEAR),OPTION,'.mat'],'Result_NOX');
Visualization_USResult_1(strrep(TEMPLATE,'$Date$',['Result_NOX',datestr(datenum(YEAR,1,1),'yyyymmdd'),OPTION,'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),OPTION]),nanmean(Result_NOX),SiteData_uni_map,EnvironPara);

save([OUTPUTPATH,'CMAQ_NO2_',SITENAME_UNI,'_',num2str(YEAR),OPTION,'_',num2str(YEAR),OPTION,'.mat'],'Result_NO2');
Visualization_USResult_1(strrep(TEMPLATE,'$Date$',['Result_NO2',datestr(datenum(YEAR,1,1),'yyyymmdd'),OPTION,'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),OPTION]),nanmean(Result_NO2),SiteData_uni_map,EnvironPara);

save([OUTPUTPATH,'CMAQ_NO2_Vertical_',SITENAME_UNI,'_',num2str(YEAR),OPTION,'_',num2str(YEAR),OPTION,'.mat'],'Result_NO2_Vertical');
Visualization_USResult_1(strrep(TEMPLATE,'$Date$',['Result_NO2_Vertical',datestr(datenum(YEAR,1,1),'yyyymmdd'),OPTION,'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),OPTION]),nanmean(Result_NO2_Vertical),SiteData_uni_map,EnvironPara);

save([OUTPUTPATH,'CMAQ_Ozone_',SITENAME_UNI,'_',num2str(YEAR),OPTION,'_',num2str(YEAR),OPTION,'.mat'],'Result_Ozone');
Visualization_USResult_1(strrep(TEMPLATE,'$Date$',['Result_Ozone',datestr(datenum(YEAR,1,1),'yyyymmdd'),OPTION,'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),OPTION]),nanmean(Result_Ozone),SiteData_uni_map,EnvironPara);

save([OUTPUTPATH,'CMAQ_Ozone_Vertical_',SITENAME_UNI,'_',num2str(YEAR),OPTION,'_',num2str(YEAR),OPTION,'.mat'],'Result_Ozone_Vertical');
Visualization_USResult_1(strrep(TEMPLATE,'$Date$',['Result_Ozone_Vertical',datestr(datenum(YEAR,1,1),'yyyymmdd'),OPTION,'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),OPTION]),nanmean(Result_Ozone_Vertical),SiteData_uni_map,EnvironPara);

save([OUTPUTPATH,'CMAQ_RH_',SITENAME_UNI,'_',num2str(YEAR),OPTION,'_',num2str(YEAR),OPTION,'.mat'],'Result_RH');
Visualization_USResult_1(strrep(TEMPLATE,'$Date$',['Result_RH',datestr(datenum(YEAR,1,1),'yyyymmdd'),OPTION,'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),OPTION]),nanmean(Result_RH),SiteData_uni_map,EnvironPara);

save([OUTPUTPATH,'CMAQ_TA_',SITENAME_UNI,'_',num2str(YEAR),OPTION,'_',num2str(YEAR),OPTION,'.mat'],'Result_TA');
Visualization_USResult_1(strrep(TEMPLATE,'$Date$',['Result_TA',datestr(datenum(YEAR,1,1),'yyyymmdd'),OPTION,'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),OPTION]),nanmean(Result_TA),SiteData_uni_map,EnvironPara);

save([OUTPUTPATH,'CMAQ_AIR_DENS_',SITENAME_UNI,'_',num2str(YEAR),OPTION,'_',num2str(YEAR),OPTION,'.mat'],'Result_AIR_DENS');
Visualization_USResult_1(strrep(TEMPLATE,'$Date$',['Result_AIR_DENS',datestr(datenum(YEAR,1,1),'yyyymmdd'),OPTION,'_',datestr(datenum(YEAR,12,31),'yyyymmdd'),OPTION,]),nanmean(Result_AIR_DENS),SiteData_uni_map,EnvironPara);

