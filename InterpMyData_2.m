%% the main function to do 2D interpolation

%% parameter;
% MonitorData: data to be interpolated n*m matrix; n is number of days; m is the number of grid cells;
% lon1: m*1 matrix, longitude of input data 
% lat1: m*1 matrix, latitude of input data
% lon2: longitude of grid cells on which we will interpolate at;
% lat2: latitude of grid cells on which we will interpolate at;
% option: what kind of interpolation we need? 
% 'nearest': interpolating the value based on the nearest monitoring site with available data; 
% 'the_nearest': interpolating the value based on the nearest monitoring, regardless whether it has available data or missing data
% 'Threshold': normal 2D linear interpolation within a distance;
% '': normal 2D linear interpolation (default)

%% return value
% ResultData: interpolated value at [lon2, lat2];

%% code
function ResultData = InterpMyData_2(MonitorData,lon1,lat1,lon2,lat2,option)

if(size(MonitorData,1)==1)
    ONEDAY = true; 
else
    ONEDAY = false;
end

N_Sites = size(MonitorData,2);
N_Days = size(MonitorData,1);
N_Location = length(lon2);

if(N_Sites~=length(lon1))
   fprintf('the number of monitors should be consistent!\n');
   return;
end

%the output data
ResultData = nan(N_Days,length(lon2));

%% nearest interpolation
% there are some point that is assigned to nan value after natural interpolation,
% we use nearest neighbour to assign value
if(strcmp(option,'nearest'))
    if(true == ONEDAY)%%%MonitorData have only one day data
        
        
    else
        
    end
    
    for i=1:N_Days
       Index = isnan(MonitorData(i,:));
       
       if(sum(Index)==length(Index))
           continue; 
       end
       
%        fprintf('%d\n',i);
       F1 = TriScatteredInterp(lon1(~Index),lat1(~Index),MonitorData(i,~Index)');
       F2 = TriScatteredInterp(lon1(~Index),lat1(~Index),MonitorData(i,~Index)','nearest');
       try
           Temp1 = F1(lon2,lat2);
           Temp2 = F2(lon2,lat2);
           Index1 = isnan(Temp1);
           Temp1(Index1) = Temp2(Index1);
           
           ResultData(i,:) = Temp1;
       catch exception
           fprintf('%s\n',exception.message);
           fprintf('something goes wrong in the interpolation part 1! Error from InterpMyData_2!\n');
       end
    end

%% interpolating the value based on the nearest monitoring, regardless whether it has available data or missing data
elseif(strcmp(option,'the_nearest'))  
    
    for i=1:N_Days
       Index = isnan(MonitorData(i,:));
%        fprintf('%d\n',i);

       try
           F1 = TriScatteredInterp(lon1(~Index),lat1(~Index),MonitorData(i,~Index)','nearest');
           ResultData(i,:) = F1(lon2,lat2);
       catch exception
           fprintf('%s\n',exception.message);
           fprintf('something goes wrong in the interpolation part 1.5! Error from InterpMyData_2!\n');
       end
    end
    
 elseif(strcmp(option,'Distance'))%%%
    
    if(sum(sum(MonitorData<=0))>0)
       disp('negative values detected!\n');
       return;
    end
    
    
    F1 = TriScatteredInterp(lon1,lat1,(1:N_Sites)','nearest');
    ListIndex = F1(lon2,lat2);
    
    Points = [lon2,lat2];
    NearPoints = [lon1(ListIndex),lat1(ListIndex)];
    NearDist = 1-tansig(sqrt(sum((Points-NearPoints).^2,2)')*15);
    
    for i=1:N_Days
        TempData = MonitorData(i,ListIndex);
        ResultData(i,:)=NearDist.*sqrt(TempData); 
    end
    
%% normal 2D linear interpolation within a distance;
elseif(any(strfind(option,'Threshold')))%%% nearest neighourhood interpolate within a certain threshold\
    THRESHOLD = str2double(strrep(option,'Threshold',''));% threshold value is written as part of option
%     fprintf('threshold value is %d\n',THRESHOLD);
%     disp('remember to set threshold value!');
    
    for i=1:N_Days
        % find out which points are within the threshold distance of any monitor
        Index = isnan(MonitorData(i,:));
        if(sum(~Index)<5)
            continue;
        end
        
        [~,D] = knnsearch([lon1(~Index),lat1(~Index)],[lon2,lat2]);
        Index_Outside_Threshold = D>THRESHOLD;
        
        % nearest neighourhood interpolation
        try
           F1 = TriScatteredInterp(lon1(~Index),lat1(~Index),MonitorData(i,~Index)','nearest');
           ResultData(i,:) = F1(lon2,lat2);
       catch exception
           fprintf('%s\n',exception.message);
           fprintf('something goes wrong in the Threshold interpolation! Error from InterpMyData_2!\n');
        end
        ResultData(i,Index_Outside_Threshold) = nan;
        
    end
  
%% normal 2D linear interpolation
else
    for i=1:N_Days
       Index = isnan(MonitorData(i,:));
       
      if(sum(Index)==length(Index))
           continue; 
       end
       
    %    fprintf('%d\n',i);
       F = TriScatteredInterp(lon1(~Index),lat1(~Index),MonitorData(i,~Index)');
       try
           ResultData(i,:) = F(lon2,lat2);    
       catch exception
           fprintf('%s\n',exception.message);
           fprintf('something goes wrong in the interpolation part 2! Error from InterpMyData_2!\n');   
       end
    end
end