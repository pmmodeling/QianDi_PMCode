DirPath = 'C:\Users\qid335\Downloads\Temp\';
zipcode = readtable([DirPath,'zipcode_SUBSETSite.csv']);

YEAR = 2009;
% get monitoring data
T = readtable([DirPath,'annual_conc_by_monitor_',num2str(YEAR),'.csv']);
T = T(strcmp(T.ParameterCode,'88101')&strcmp(T.PollutantStandard,'PM25 Annual 2012'),:);

for i=1:yeardays(YEAR)
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
