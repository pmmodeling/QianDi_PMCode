%% conduct descriptive analysis on the output and save them
% calcualte R2, make scatter plot
% this function is used in the process of model training

%% paramter
% MonitorPredicted: fitted value: N*M format; N is the number of days and M is the
% number of monitoring site;
% MonitorData: monitoring data: N*M format; N is the number of days and M is the
% number of monitoring site;
%%option: 1: calculate R2 and write to file;2: if 1 is true, scatter plot;
%%3: save result (not used);4: descriptive stat

%% code
function [R_CV,R_spatial,R_temporal]=Analysis_ResultDescriptiveAnalysis(MonitorPredicted,MonitorData,NetPredicted,EnvironPara,NOTICE,IsSaved)
disp('Analysis_ResultDescriptiveAnalysis');

%% calculate write R2 to files
if(IsSaved{1})
    SIZE = size(MonitorPredicted);
    X = reshape(MonitorPredicted,[SIZE(1)*SIZE(2),1]);
    Y = reshape(MonitorData,[SIZE(1)*SIZE(2),1]);

    if(sum(isnan(X))==length(X))
        R_CV = -9999;
        R_spatial = -9999;
        R_temporal = -9999;
       return; 
    end
    [b,bint,r,rint,stats] = regress(Y,[ones(length(X),1),X]);

    [R_CV,P_CV,MSE_CV] = CalculateRsquare(MonitorData,MonitorPredicted,'out-of-sample');
    [R_spatial,P_spatial,MSE_spatial] = CalculateRsquare(MonitorData,MonitorPredicted,'spatial');
    [R_temporal,P_temporal,MSE_temporal] = CalculateRsquare(MonitorData,MonitorPredicted,'temporal');
    %%% print result;
    fprintf('%s\n',repmat('-',[1,50]));
    fprintf('%s\t',NOTICE);
    fprintf('%s\t%s\t',EnvironPara.CommonStartPoint,EnvironPara.CommonEndPoint);
    fprintf('%s\t%s\t',EnvironPara.NAME);
    fprintf('out-of-sample\t%f\t%d\t%d\t',R_CV,P_CV,MSE_CV);
    fprintf('spatial\t%f\t%d\t%d\t',R_spatial,P_spatial,MSE_spatial);
    fprintf('temporal\t%f\t%d\t%d\t',R_temporal,P_temporal,MSE_temporal);
    fprintf('b0\t%d\tb1\t%d\tR2\t%d\n',b(1),b(2),stats(1));
    
    fprintf(EnvironPara.FID,'%s\tout-of-sample\t%f\t%d\t%d\tspatial\t%f\t%d\t%d\ttemporal\t%f\t%d\t%d\tb0\t%d\tb1\t%d\tR2\t%d\n',NOTICE,R_CV,P_CV,MSE_CV,R_spatial,P_spatial,MSE_spatial,R_temporal,P_temporal,MSE_temporal,b(1),b(2),stats(1));
    if(~isempty(EnvironPara.FID_SUMMARY))
        fprintf(EnvironPara.FID_SUMMARY,'%s\tout-of-sample\t%f\t%d\t%d\tspatial\t%f\t%d\t%d\ttemporal\t%f\t%d\t%d\tb0\t%d\tb1\t%d\tR2\t%d\n',NOTICE,R_CV,P_CV,MSE_CV,R_spatial,P_spatial,MSE_spatial,R_temporal,P_temporal,MSE_temporal,b(1),b(2),stats(1));
    end
end

%% write to images
if(IsSaved{1} && IsSaved{2})
    Analysis_ScatterPlot(X,Y,b(1),b(2),EnvironPara.NAME,EnvironPara.NAME,EnvironPara.OUTPUTPATH_PIC,datestr(EnvironPara.FIRSTDAY,'yyyymmdd'),datestr(EnvironPara.LASTDAY,'yyyymmdd'),NOTICE);
end

%% store to disk
if(IsSaved{3})
    SiteData = EnvironPara.SiteData_Model;
    save([EnvironPara.OUTPUTPATH,NOTICE,'Output_',EnvironPara.NAME,'_',EnvironPara.SITENAME_MODEL,'_',EnvironPara.CommonStartPoint,'_',EnvironPara.CommonEndPoint,'.mat'],'MonitorPredicted','MonitorData','NetPredicted','SiteData');  
end

%% descriptive analysis
if(IsSaved{4})    
    SIZE = size(MonitorPredicted);
    X = reshape(MonitorPredicted,[SIZE(1)*SIZE(2),1]);
    
    fprintf('%s\t%s\t%s\tnan count:%d\tnon-nan percentage:%f\tis real:%d\t',NOTICE,EnvironPara.CommonStartPoint,'Result',sum(isnan(X)),sum(~isnan(X))/length(X),isreal(X)); 
    fprintf('dist:%f\t%f\t%f\t%f\t%f\t%f\t%f\n',nanmin(X),quantile(X,0.05),quantile(X,0.25),nanmean(X),quantile(X,0.75),quantile(X,0.95),nanmax(X)); 
    
    fprintf(EnvironPara.FID,'%s\t%s\t%s\tnan count:%d\tnon-nan percentage%f\tis real:%d\t',NOTICE,EnvironPara.CommonStartPoint,'Result',sum(isnan(X)),sum(~isnan(X))/length(X),isreal(X)); 
    fprintf(EnvironPara.FID,'dist:%f\t%f\t%f\t%f\t%f\t%f\t%f\n',nanmin(X),quantile(X,0.05),quantile(X,0.25),nanmean(X),quantile(X,0.75),quantile(X,0.95),nanmax(X)); 
end

end        


%% make scatter plot between monitoring data and fitted data
% add fitted linear to the scatter plot as well
function Analysis_ScatterPlot(X,Y,b0,b1,NAME,TRACERNAME,DIRPATH,FIRSTDAY,LASTDAY,variable)

Notice = variable;

x = linspace(min(X),max(X),10);
y = b1*x+b0;
scatter(X,Y,5,'o');hold on
plot(x,y,'Color','black','LineWidth',1);hold off
% title(['GEOSCHEM:',regexprep(TRACERNAME,'_',' '),' VS ','MONITOR:',regexprep(NAME,'_',' '),' FROM ',datestr(FIRSTDAY),' TO ',datestr(LASTDAY)],'FontSize',10);
% title(['GEOSCHEM:',regexprep(TRACERNAME,'_',' '),' VS ','MONITOR:',regexprep(NAME,'_',' '),' FROM ',num2str(year(FIRSTDAY)),' TO ',num2str(year(LASTDAY))],'FontSize',8);
title([TRACERNAME,' Time ',FIRSTDAY],'FontSize',20);
% title([regexprep(TRACERNAME,'_',' '),' vs. ', regexprep(NAME,'_',' '),'year:',num2str(year(FIRSTDAY))],'FontSize',12);
% axis([-(max(x)-min(x))/10 max(x)+(max(x)-min(x))/10 -(max(x)-min(x))/10 max(x)+(max(x)-min(x))/10])
xmin = -quantile(X,0.10);
xmax = quantile(X,0.99)+quantile(X,0.10);
ymin = -quantile(Y,0.10);
ymax = quantile(Y,0.99)+quantile(Y,0.10);
% axis([xmin xmax ymin ymax]);
ylabel('Calibrated','FontSize',18);
xlabel('Monitor','FontSize',18);
 

%% save figure
PicWidth = 12.35;% 8.3, 12.35, 17.35, 23.35cm
PicHeight = PicWidth/4*3;
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 PicWidth PicHeight]);
set(gca,'FontSize',14);
text(0,ymax,['y=',num2str(b1),'*x+',num2str(b0)],'FontSize',14);
print(gcf,'-dtiff','-r600',[DIRPATH,'Figure',TRACERNAME,NAME,'_',FIRSTDAY,'_',LASTDAY,'.tif']);
% saveas(h,[DIRPATH,'Figure',TRACERNAME,NAME,'_',FIRSTDAY,'_',LASTDAY,'.fig'],'fig');
end
