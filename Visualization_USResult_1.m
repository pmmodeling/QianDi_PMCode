%% the main function to visualize results

%% parameter;
% FileName: name for output map
% Result: data to be visualized, a 1*n matrix;
% SiteData: a n*3 table, with ID, lat, lon
% EnvironPara: specify where to store output picture; where to read shape file; 
% an example is :
% EnvironPara = struct(...
%     'OUTPUTPATH_PIC','C:\Users\qid335\Downloads\',...
%     'US_North_America_Equidistant_Conic','D:\Google Drive\Research\USTemperature\DataPaper3\Shapefile\US_North_America_Equidistant_Conic.shp',...
%     'US','D:\Google Drive\Research\USTemperature\DataPaper3\Shapefile\US.shp'...
%     );

%% code
function Visualization_USResult_1(FileName,Result,SiteData,EnvironPara)

%% if SiteData not in the proper format -- use the only first 10 records
if(iscell(SiteData))
    SiteData = table((1:1:size(SiteData,1))',cell2mat(SiteData(:,2)),cell2mat(SiteData(:,3)),'VariableNames',{'SiteCode','Lat','Lon'});
end


%% the code for visualization
[~,D] = knnsearch([SiteData.Lon(2:10),SiteData.Lat(2:10)],[SiteData.Lon(1),SiteData.Lat(1)],'k',1);

if(nanmax(SiteData.Lat)>90||nanmin(SiteData.Lat)<-90||nanmin(SiteData.Lon)<-180||nanmax(SiteData.Lon)>180)
    %% the whole country North_America_Equidistant_Conic
    LonMax =  2224207.124100;
    LonMin = -2330934.200700;
    LatMin = -1713278.696000;
    LatMax = 1377766.963500;
    
    % %% for NE only
%     LonMax = 2.2e+06;
%     LonMin = 1.0e+06;
%     LatMin = -3.1+05;
%     LatMax = 1.2e+05;
    
    S1 = shaperead(EnvironPara.US_North_America_Equidistant_Conic);
    Radius = sqrt(D/1000*1)*1;
else%% lat/lon
    LonMax = -64;
    LonMin = -128;
    LatMin = 24;
    LatMax = 52;
    S1 = shaperead(EnvironPara.US);
    Radius = sqrt(D/0.01*1)*1;
end

% constraint within the study area
Index2 = SiteData.Lon>LonMin&SiteData.Lon<LonMax&SiteData.Lat>LatMin&SiteData.Lat<LatMax;
SiteData = SiteData(Index2,:);
Result = Result(Index2);

% take a sample to visualize
Prob = min(1,1000000/size(SiteData,1));
Index1 = ismember((1:1:size(SiteData,1))',randsample(size(SiteData,1),floor(size(SiteData,1)*Prob)));%1:100:size(SiteData,1);
SiteData = SiteData(Index1,:);
Result = Result(Index1);

%% plot map and add colorbar
N = 64;
ColorScale = jet(N);
A2 = discretize(Result,[min(Result),linspace(quantile(Result,0.01),quantile(Result,0.99),N-2),max(Result)]);
MyColorScale = ones(length(Result),3);
MyColorScale(~isnan(A2),:) = ColorScale(A2(~isnan(A2)),:);
scatter(SiteData.Lon,SiteData.Lat,Radius,MyColorScale,'filled','s');hold on
colormap jet
MyColorScaleLabel = [quantile(Result,0.01),quantile(Result,0.25),nanmean(Result),quantile(Result,0.75),quantile(Result,0.99)];
% MyColorScaleLabel = [0.05,0.1,0.12,0.13,0.18];
h = colorbar('Ticks',[0.01,0.25,0.5,0.75,0.99],'TickLabels',num2cell(MyColorScaleLabel));
% h.Label.String = 'Temperature (Kelvin)';

%% plot us boundary
plot([S1.X],[S1.Y],'black');hold off


%% export the map into tif format
axis([LonMin LonMax LatMin LatMax]);

PicWidth = 50;% 8.3, 12.35, 17.35, cm
PicHeight = PicWidth/2;
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 PicWidth PicHeight]);
% set(gca,'FontSize',4);
title(strrep(FileName,'_',' '),'FontSize',8);
print(gcf,'-dtiff','-r300',[EnvironPara.OUTPUTPATH_PIC,FileName,'.tif']);
