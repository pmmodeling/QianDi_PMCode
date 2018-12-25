%% read Infogroup US Historical Business Data, extract restaurant data, and export in this format: lat, lon, sales volume

%% parameters:
%YEAR: year of US Historical Business Data to be processed

%% example:
% Read_RestaurantData(2000)
function Read_RestaurantData(YEAR)

%% macro definition
DIRPATH = '../data/unprocessed/LANDUSE/Restaurant/';
fid = fopen([DIRPATH,'Infogroup_',num2str(YEAR),'_Business_QCQ.txt'],'r');

%% read US Historical Business Data
C = textscan(fid,'%q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q',1,'Delimiter',',');

% constrain locations to places within the continental U.S.
Index_Lat = 47;
Index_Lon = 48;
Index_Code = 14;
Index_sales = 29;


%% only extract restaurant, its lat/on and sales, and output into csv
% the restaurant code starts with "722511", so this file goes through each
% line of the raw data, and takes the records belong to restaurant category
Result = cell2table(cell(0,3), 'VariableNames', {'Sales','Lon','Lat'});
k = 0;
while(true)
    try
       k = k + 1;
       C = textscan(fid,'%q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q %q',1,'Delimiter',','); 
       if(length(C{Index_Code}{1})>6)
           if(strcmp(C{Index_Code}{1}(1:6),'722511'))
               Temp.Lon = str2double(C{Index_Lon}{1});
               Temp.Lat = str2double(C{Index_Lat}{1});
               Temp.Sales = str2double(C{Index_sales}{1});
               Result = [Result;struct2table(Temp)];
           else
              continue; 
           end
       end
    catch exception
        fprintf('%d\t%s\n',k,exception.message);
       break; 
    end
    
end
% save the data
save([DIRPATH,num2str(YEAR),'_Business_722511_Academic_QCQ.mat'],'Result');
writetable(Result,[DIRPATH,num2str(YEAR),'_Business_722511_Academic_QCQ.csv']);
