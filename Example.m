%% how to read and write in hdf5

% # in matlab
% hdf5write('myfile.h5', 'dataset1', uint8(magic(5)));
% 
% h5create('myfile.h5','/Result',[1 2156],'Datatype','double', 'ChunkSize',[1 2156],'Deflate',0)
% h5write('myfile.h5','/Result',Result)

% # in R
% source("https://bioconductor.org/biocLite.R")
% biocLite("rhdf5")
% library(rhdf5)
% f <- h5read("C:\\Users\\qid335\\Documents\\GitHub\\AODModelWork2\\myfile.h5","dataset1")

% for i = 2000:2016
%    Result = load(['D:\Google Drive\Research\USTemperature\processed_data\EPACastNetOzone\Monitor\MONITOR_Ozone_EPACastNetOzone_',num2str(i),'_',num2str(i),'.mat'],'Result');
%    Result = Result.Result;
%    hdf5write(['D:\Google Drive\Research\USTemperature\processed_data\EPACastNetOzone\Monitor\MONITOR_Ozone_EPACastNetOzone_',num2str(i),'_',num2str(i),'.h5'],'Result',Result) 
% end

TempFile = 'D:\Google Drive\Research\USTemperature\processed_data\EPANO2\Location\EPANO2Site_North_America_Equidistant_Conic';
load([TempFile,'.mat']);
hdf5write([TempFile,'.h5'], 'SiteCode',SiteData.SiteCode, 'Latitude',SiteData.Lat,'Longitude',SiteData.Lon);