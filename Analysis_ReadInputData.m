%% this function read all input variable based on the VariableList.csv
% the information inside the VariableList.csv has been stored in the
% variable "EnvironPara", including folder name, file name, time duration.
% and so on. 

%% parameter:
% EnvironPara: variable that store which varible to read;
% OPTION: 'Prediction': read the data by day; OPTION,'All' or 'CV': read
% the data by year

%% return value:
% Input_var: the merged data; each column stands for different variabe;
% there are N*M rows; N is the number of days; and M is the number of
% monitoring site;
% VariableList: variable name for each column

%% version history
%% note
%2015-07-07 add urbanness
%% version 3
%for convultion Neural N
%change the data format into N_Day*N_Site format (also the final data will be reshaped to [N_Day*N_Site,1])
%% version 3.1
%add covolution template of land-use and AOD variables.
%% version 3.2
%DATE_A, DATE_B as date/year, not string;
%covolution template could be read from disk
%% version 3.3
%organize the global parameters
%% version 3.4
%since I used another script to produce neighbour files, there should be
%consequent changes in reading files
%code is much more simplified
%% version 3.5
%can read MonitorPredicted as input
%% version 3.6
%read simulated AOD data
%% version 3.7
% let EnvironPara to determine variable list
% rename Aqua*GC_scale as Aqua_scale;
%% version 4
% only read variables in the VariableList; change format
%% 2016-08-29: change name from CalibrateGC_ReadAODData_4 to Analysis_ReadInputData
%% 2016-10-02: calculate windspeed is smarter now!
%% 2017-05-25: read data from multiple years
function [Input_var,VariableList] = Analysis_ReadInputData(EnvironPara,OPTION)
disp('Analysis_ReadInputData');

Sep = EnvironPara.Sep;

%% start!!!

Input_var = [];

StartYear = year(EnvironPara.StartDate);
EndYear = year(EnvironPara.EndDate);
for j = StartYear:EndYear
    if(~isempty(strfind(EnvironPara.OPTION,'CV'))||~isempty(strfind(EnvironPara.OPTION,'All')))
        N_Day = yeardays(j);
        StartDate = datenum(j,1,1);
        EndDate = datenum(j,12,31);
        StartDateString = datestr(StartDate,'yyyymmdd');
        EndDateString = datestr(EndDate,'yyyymmdd');
        StartYearString = num2str(j);
        EndYearString = num2str(j);
    elseif(strcmp(EnvironPara.OPTION,'Prediction'))
        N_Day = 1;
        StartDate = EnvironPara.StartDate;
        EndDate = EnvironPara.EndDate;
        StartDateString = datestr(StartDate,'yyyymmdd');
        EndDateString = datestr(EndDate,'yyyymmdd');
        StartYearString = num2str(j);
        EndYearString = num2str(j);
    end

    
    fprintf('reading temperature and related data...%s %s\n',StartDateString,EndDateString);
    %% reading Site data
    SiteData_Model = EnvironPara.SiteData_Model;
    % read location data;
    SiteData = LoadData_function([EnvironPara.DIRPATH_DATA,'Location',Sep,EnvironPara.SITENAME_DATA,'Site_',EnvironPara.GCS,'.mat'],'SiteData');
    SiteData_Data = SiteData.SiteData;

    if(isequal(SiteData_Model,SiteData_Data))
        Lia = true(size(SiteData_Model,1),1);
        Locb1 = (1:1:size(SiteData_Model,1))';
    else
        if(iscell(SiteData_Data))
            [Lia,Locb] =  ismember(cell2mat(SiteData_Data(:,1)),cell2mat(SiteData_Model(:,1)));
            Locb1 = Locb(Lia);
        elseif(istable(SiteData_Data))
            [Lia,Locb] =  ismember(SiteData_Data.SiteCode,SiteData_Model.SiteCode);
            Locb1 = Locb(Lia);
        else
            [Lia,Locb] = (ismember(SiteData_Data(:,1),SiteData_Model(:,1)));
            Locb1 = Locb(Lia);
        end
    end
    N_Site_Model = size(SiteData_Model,1);

    %% for most variables
    for i=1:length(EnvironPara.VariableListFull)
        %% can be directly read from disk
        if(strcmp(EnvironPara.VariableRead{i},'T'))
            %%generate file name based on the general file name template
            TempFileName = strrep(EnvironPara.VariableNameTemplate{i},'$NAME$',EnvironPara.NAME);
            TempFileName = strrep(TempFileName,'$EndYear$',EndYearString);
            TempFileName = strrep(TempFileName,'$StartYear$',StartYearString);
            TempFileName = strrep(TempFileName,'$SITENAME_DATA$',EnvironPara.SITENAME_DATA);
            TempFileName = strrep(TempFileName,'$SITENAME_MODEL$',EnvironPara.SITENAME_MODEL);

            if(strcmp(EnvironPara.VariableFormat{i},'1') || strcmp(EnvironPara.VariableFormat{i},'3'))
                TempFileName = strrep(TempFileName,'$StartDate$',StartDateString);
                TempFileName = strrep(TempFileName,'$EndDate$',EndDateString);
            elseif(strcmp(EnvironPara.VariableFormat{i},'2'))
                TempFileName = strrep(TempFileName,'$StartDate$',datestr(datenum(j,1,1),'yyyymmdd'));
                TempFileName = strrep(TempFileName,'$EndDate$',datestr(datenum(j,12,31),'yyyymmdd'));
            end

            %%generate folder
            TempFolderName = [EnvironPara.DIRPATH_MODEL,EnvironPara.VariableFolder{i},EnvironPara.Sep];
            %%Load file
            try
                fprintf('reading %s\npath %s\n',EnvironPara.VariableListFull{i},[TempFolderName,TempFileName]);

                Result = LoadData_function([TempFolderName,TempFileName],'Result');
                Result = Result.Result;
                if(strcmp(EnvironPara.VariableFormat{i},'2'))
                    Temp1 = repmat(Result,[N_Day,1]);
                elseif(strcmp(EnvironPara.VariableFormat{i},'3'))
                    Temp1 = repmat(Result',[N_Day,1]);
                else
                    Temp1 = Result;
                end

                if(~isempty(strfind(EnvironPara.VariableNameTemplate{i},'$SITENAME_DATA$')))    
                    %%on full dataset, we need clip
                    Temp1 = Temp1(:,Lia);
                    Temp2 = nan(size(Temp1));% Temp2 is introduced to deal with the issues of index mismatch in subset
                    Temp2(:,Locb1) = Temp1 ;
                elseif(~isempty(strfind(EnvironPara.VariableNameTemplate{i},'$SITENAME_MODEL$')))
                    %%on subset, we don't need to
                    Temp2 = Temp1;
                end
                eval([EnvironPara.VariableListFull{i},'=Temp2;']);
                clear Temp1 Temp2 Result

            catch exception% at prediction mode, we do not need MonitorData, replace it with nan
                if(any(strfind(exception.message,'file not existed')) && (j < EnvironPara.VariableStartYear(i) ||  j > EnvironPara.VariableEndYear(i)))        
                    fprintf('oops!%s\n',exception.message);
                    Temp2 = nan([N_Day,N_Site_Model]);
                    eval([EnvironPara.VariableListFull{i},'=Temp2;']);
                    clear Temp1 Temp2 Result
                elseif(any(strfind(exception.message,'file not existed')) && (j >= EnvironPara.VariableStartYear(i) &&  j <= EnvironPara.VariableEndYear(i)))
                    % record which file has corrupted
                    if(any(strfind(EnvironPara.OPTION,'Prediction')) && any(strfind(EnvironPara.VariableListFull{i},'MonitorData')))
                        Temp2 = nan([N_Day,N_Site_Model]);
                        eval([EnvironPara.VariableListFull{i},'=Temp2;']);
                        clear Temp1 Temp2 Result
                    else
                        fprintf(EnvironPara.FID_error_files,'%s,%s\n',datestr(now),[TempFolderName,TempFileName]);
                        fprintf('error!%s\n',exception.message);
                        Input_var = nan;
                        VariableList = nan;
                        return;
                    end
                elseif(any(strfind(exception.message,'file corrupted')))
                    % record which file has corrupted
                    fprintf(EnvironPara.FID_error_files,'%s,%s\n',datestr(now),[TempFolderName,TempFileName]);
                    fprintf('error!%s\n',exception.message);
                    Input_var = nan;
                    VariableList = nan;
                    return;
                end            
            end
        end
    end

    %% Calculate wind speed
    windspeed_list = arrayfun(@(i) strfind(EnvironPara.VariableList{i},'REANALYSIS_windspeed_10m'),1:length(EnvironPara.VariableList),'UniformOutput',false);

    for i = 1:length(windspeed_list)
       if(1 == windspeed_list{i})%% it is wind speed
           TempVar = EnvironPara.VariableList{i};
           fprintf('calculating...%s\n',TempVar);
           TempVar1 = strrep(TempVar,'windspeed','vwnd');
           TempVar2 = strrep(TempVar,'windspeed','uwnd');
           TempValue = sqrt((eval(TempVar1)).^2 + (eval(TempVar2)).^2);
           eval([TempVar,'=TempValue;']);
       end
    end

    %%other variable, dummy variables for month, day, and calendar day
    if(ismember('Other_Month',EnvironPara.VariableList))
        disp('indicator variable for month...');
        Other_Month = repmat(month(StartDate:1:EndDate)',[1,EnvironPara.N_Site]);
    end

    if(ismember('Other_Day',EnvironPara.VariableList))
        disp('indicator variable for day...');
        Other_Day = repmat(day(StartDate:1:EndDate)',[1,EnvironPara.N_Site]);
    end
    
    if(ismember('CalendarDay',EnvironPara.VariableList))
        disp('indicator variable for CalendarDay...');
        CalendarDay = repmat((StartDate:1:EndDate)',[1,EnvironPara.N_Site]);
    end
    
    %% create monthly dummy variable
    if(ismember('DM1',EnvironPara.VariableList)&&ismember('DM2',EnvironPara.VariableList)&&...
            ismember('DM3',EnvironPara.VariableList)&&ismember('DM4',EnvironPara.VariableList)&&...
            ismember('DM5',EnvironPara.VariableList)&&ismember('DM6',EnvironPara.VariableList)&&...
            ismember('DM7',EnvironPara.VariableList)&&ismember('DM8',EnvironPara.VariableList)&&...
            ismember('DM9',EnvironPara.VariableList)&&ismember('DM10',EnvironPara.VariableList)&&...
            ismember('DM11',EnvironPara.VariableList)&&ismember('DM12',EnvironPara.VariableList))
        b = dummyvar(month(StartDate:1:EndDate));
        b = reshape(b,[size(b,1),1,size(b,2)]);
        a = repmat(b,[1,N_Site_Model,1]);
        DM1 = a(:,:,1);
        DM2 = a(:,:,2);
        DM3 = a(:,:,3);
        DM4 = a(:,:,4);
        DM5 = a(:,:,5);
        DM6 = a(:,:,6);
        DM7 = a(:,:,7);
        DM8 = a(:,:,8);
        DM9 = a(:,:,9);
        DM10 = a(:,:,10);
        DM11 = a(:,:,11);
        DM12 = a(:,:,12);
    end

    %% add lat and lon for each monitoring site
    if(ismember('Other_Lat',EnvironPara.VariableList))
        disp('location variable for Lat...');
        Other_Lat = repmat((EnvironPara.SiteData_Model.Lat)',[N_Day,1]);
    end

    if(ismember('Other_Lon',EnvironPara.VariableList))
        disp('location variable for Lon...');
        Other_Lon = repmat((EnvironPara.SiteData_Model.Lon)',[N_Day,1]);
    end
    
    if(ismember('Spatial_Lagged_1',EnvironPara.VariableList) && ...
            ismember('Spatial_Lagged_2',EnvironPara.VariableList) && ...
            ismember('Spatial_Lagged_3',EnvironPara.VariableList) && ...
            ismember('Temporal_Lagged_1',EnvironPara.VariableList) && ...
            ismember('Temporal_Lagged_2',EnvironPara.VariableList) && ...
            ismember('Temporal_Lagged_3',EnvironPara.VariableList))
       
        if(strcmp(OPTION,'first'))
            disp('read estimation-based spatial/temporal lagged terms...default values');
        %     Result = CalibrateGC_CreateFilter_3(log(ModelDataGEOSChem(:,Index_SiteCode)),[],[],EnvironPara.SITENAME_PREDICT,EnvironPara.SITENAME_MODEL,EnvironPara,'UPDATE_IN_TRAINING');
            Result = ones(N_Day,EnvironPara.N_Site,6);
            Spatial_Lagged_1 = Result(:,:,1);
            Spatial_Lagged_2 = Result(:,:,2);
            Spatial_Lagged_3 = Result(:,:,3);
            Temporal_Lagged_1 = Result(:,:,4);
            Temporal_Lagged_2 = Result(:,:,5);
            Temporal_Lagged_3 = Result(:,:,6);
            clear Result

        elseif(strcmp(OPTION,'second')||strcmp(OPTION,'final'))

            ResultTemp = LoadData_function([EnvironPara.DIRPATH_MODEL,'LaggedTerms',Sep,'SpatialLaggedTerms1_',EnvironPara.NAME,'_',EnvironPara.SITENAME_PREDICT,'_',StartDateString,'_',EndDateString,'.mat'],'Result');
            Spatial_Lagged_1 = ResultTemp.Result';
            clear ResultTemp
            ResultTemp = LoadData_function([EnvironPara.DIRPATH_MODEL,'LaggedTerms',Sep,'SpatialLaggedTerms2_',EnvironPara.NAME,'_',EnvironPara.SITENAME_PREDICT,'_',StartDateString,'_',EndDateString,'.mat'],'Result');
            Spatial_Lagged_2 = ResultTemp.Result';
            clear ResultTemp
            ResultTemp = LoadData_function([EnvironPara.DIRPATH_MODEL,'LaggedTerms',Sep,'SpatialLaggedTerms3_',EnvironPara.NAME,'_',EnvironPara.SITENAME_PREDICT,'_',StartDateString,'_',EndDateString,'.mat'],'Result');
            Spatial_Lagged_3 = ResultTemp.Result';
            clear ResultTemp
            ResultTemp = LoadData_function([EnvironPara.DIRPATH_MODEL,'LaggedTerms',Sep,'TemporalLaggedTerms_',EnvironPara.NAME,'_',EnvironPara.SITENAME_PREDICT,'_',StartDateString,'_',EndDateString,'.mat'],'Result');
            Temporal_Lagged_1 = ResultTemp.Result(:,1)';
            Temporal_Lagged_2 = ResultTemp.Result(:,2)';
            Temporal_Lagged_3 = ResultTemp.Result(:,3)';
            clear ResultTemp
        end
        
    end
    
    %% form matrix of input
    Input_var_temp = eval(['cat(3,',strjoin(EnvironPara.VariableList,','),')']);
    Input_var = cat(1,Input_var,Input_var_temp);
end

Input_var = reshape(Input_var,[EnvironPara.N_Day*EnvironPara.N_Site,EnvironPara.N_Var]);
VariableList = EnvironPara.VariableList;

