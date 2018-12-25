%% one important step of the air pollution modeling is what called "two-step" 
% modeling. It means (1) we fitted the neural network (it can be any
% machine learning algorithm or even regression model) and obtained
% predicted air pollution levels at each monitoring site; then (2) we
% calculated the averaged predicted air pollution levels from neigbouring
% days and nearby monitoring site; and (3) we added those average values at
% additional input variables and fitted the neural network again. Since we
% fitted the neural network twice, we call it "two-step" modeling.
% this function is not supposed to run in a stand alone way

%% parameter
% Data: predicted air pollution level at monitoring sites, a huge matrix;
% SiteData_Predict: the grid cell we will calculate average predicted values at
% SiteData_Data: the grid cell we already have fitted model at;
% SiteName_Predict:  the grid cell name we will calculate average predicted values at
% SiteName_Data: the grid cell name we already have fitted model at;
% EnvironPara: some macro definition
% OPTION: CreateWeight: just create weight; UPDATE: create files to do
% two-step modeling

%% version history
% version 2
% SiteData_Data is the Site Data we have data at, data is the variable "Data"[N_Day,N_Site] fomrat;
% SiteData_Predict is the Site Data will make prediction at;
% if OPTION == CreateWeight --> just create and save those weight matrix,
% which takes a lot of memory
%  version 3
% not long need to update initially
%  version 4
% change weight again
% change name to Analysis_TwoStep

%% code
function Result = Analysis_TwoStep(Data,SiteData_Predict,SiteData_Data,SiteName_Predict,SiteName_Data,EnvironPara,OPTION)
disp('Analysis_TwoStep');
%% estimation-based spatial lagged terms
%%after getting better estimation, update those terms recursively
DIRPATH = [EnvironPara.DIRPATH_MODEL,'Temp',EnvironPara.Sep];
Sep = EnvironPara.Sep;
mkdir(DIRPATH);


%% just create those weights -- no threshod value/with threshold values
% these weight matries are used to calculate averaged predicted values from
% nearby monitors. Here we used matrix multiplication to do calculate such
% averaged values
if(strcmp(OPTION,'CreateWeight'))

    if(iscell(SiteData_Predict))
        Lon_Target = cell2mat(SiteData_Predict(:,3));
        Lat_Target = cell2mat(SiteData_Predict(:,2));
    elseif(istable(SiteData_Predict))
        Lon_Target = SiteData_Predict.Lon;
        Lat_Target = SiteData_Predict.Lat;
    else
        Lon_Target = SiteData_Predict(:,3);
        Lat_Target = SiteData_Predict(:,2);
    end
    N_Traget = length(Lon_Target);
    
    %% raw data location list
    if(iscell(SiteData_Data))
        Lon_Query = cell2mat(SiteData_Data(:,3));
        Lat_Query = cell2mat(SiteData_Data(:,2));
    elseif(istable(SiteData_Data))
        Lon_Query = SiteData_Data.Lon;
        Lat_Query = SiteData_Data.Lat;
    else
        Lon_Query = SiteData_Data(:,3);
        Lat_Query = SiteData_Data(:,2);
    end
    N_Query = length(Lon_Query);
    
    
    THRESHOLD = Inf;
    N_Neighbour = 200;
    
    %% without threshold values
    if(exist([DIRPATH,'SpatialLaggedWeightPeak41_',SiteName_Data,'_',SiteName_Predict,'.mat'],'file') &&...
            exist([DIRPATH,'SpatialLaggedWeightPeak42_',SiteName_Data,'_',SiteName_Predict,'.mat'],'file') &&...
            exist([DIRPATH,'SpatialLaggedWeightPeak43_',SiteName_Data,'_',SiteName_Predict,'.mat'],'file'))
        disp('Weight matrixes Peak2 have been processed!'); 
    else
        %% first find neighbours.
        CurrentFile1 = [DIRPATH,Sep,'SpatialLaggedNeighbourAdjacentList_',['Peak','_Thres',num2str(THRESHOLD)],'_',SiteName_Data,'_',SiteName_Predict,'.mat'];
        if(exist(CurrentFile1,'file'))
            Result = LoadData_function(CurrentFile1);
            d = Result.d;
            n = Result.n;
        else
            fprintf('creating neighourhood files...\n');
            [n,d] = knnsearch([Lat_Query,Lon_Query],[Lat_Target,Lon_Target],'k',N_Neighbour); 
            save(CurrentFile1,'n','d','-v7.3');
        end
        
        d(d(:,1)==0,1)=0;%
        Index_d = d<THRESHOLD;
        %% version 4, power 1
        if(~exist([DIRPATH,'SpatialLaggedWeightPeak41_',SiteName_Data,'_',SiteName_Predict,'.mat'],'file'))     
            d(Index_d)=1./d(Index_d);%1/d;
            d(isinf(d)) = 0;
            d(~Index_d)=0;
            Weight = sparse(reshape(n,[N_Traget*N_Neighbour,1]),repmat((1:1:N_Traget)',[N_Neighbour,1]),reshape(d,[N_Traget*N_Neighbour,1]),N_Query,N_Traget);
            Sum_Weight_matrix = sum(Weight,1);
            save([DIRPATH,'SpatialLaggedWeightPeak41_',SiteName_Data,'_',SiteName_Predict,'.mat'],'Weight','Sum_Weight_matrix','-v7.3');
        end
        
        %% version 4, power 2
        if(~exist([DIRPATH,'SpatialLaggedWeightPeak42_',SiteName_Data,'_',SiteName_Predict,'.mat'],'file'))     
            d(Index_d)=1./power(d(Index_d),2);%1/d^2;
            d(isinf(d)) = 0;
            d(~Index_d)=0;
            Weight = sparse(reshape(n,[N_Traget*N_Neighbour,1]),repmat((1:1:N_Traget)',[N_Neighbour,1]),reshape(d,[N_Traget*N_Neighbour,1]),N_Query,N_Traget);
            Sum_Weight_matrix = sum(Weight,1);
            save([DIRPATH,'SpatialLaggedWeightPeak42_',SiteName_Data,'_',SiteName_Predict,'.mat'],'Weight','Sum_Weight_matrix','-v7.3');
        end
        
                
        %% version 4, power 3
        if(~exist([DIRPATH,'SpatialLaggedWeightPeak43_',SiteName_Data,'_',SiteName_Predict,'.mat'],'file'))     
            d(Index_d)=1./power(d(Index_d),3);%1/d^3;
            d(isinf(d)) = 0;
            d(~Index_d)=0;
            Weight = sparse(reshape(n,[N_Traget*N_Neighbour,1]),repmat((1:1:N_Traget)',[N_Neighbour,1]),reshape(d,[N_Traget*N_Neighbour,1]),N_Query,N_Traget);
            Sum_Weight_matrix = sum(Weight,1);
            save([DIRPATH,'SpatialLaggedWeightPeak43_',SiteName_Data,'_',SiteName_Predict,'.mat'],'Weight','Sum_Weight_matrix','-v7.3');
        end


        
%         %% use sparse matrix to calculate distance matrix (convolution filter template)
%         II = cell(N_Data,1);
%         JJ = cell(N_Data,1);
%         Distance = cell(N_Data,1);
%         for i=1:N_Data
%             sqdists= (SiteData_Data(i,1)-SiteData_Predict(:,1)).^2 + (SiteData_Data(i,2)-SiteData_Predict(:,2)).^2;
%             II{i}=find(sqdists>0).';%
%             JJ{i}=repmat(i,[1,length(II{i})]);
%             Distance{i}=sqdists(II{i}).'; 
%         end
%         %% take distance decay weight --- Peak42 (no threshold) just for monitoring data
%         %% version 4, power 1
%         if(~exist([DIRPATH,'SpatialLaggedWeightPeak41_',SiteName_Data,'_',SiteName_Predict,'.mat'],'file'))
%             fprintf('calculating Weight Peak41\n');
%             Weight=sparse([II{:}],[JJ{:}],1./[Distance{:}],N_Predict,N_Data);
%             save([DIRPATH,'SpatialLaggedWeightPeak41_',SiteName_Data,'_',SiteName_Predict,'.mat'],'Weight','-v7.3');
%         end
%         
%         %% version 4, power 2
%         if(~exist([DIRPATH,'SpatialLaggedWeightPeak42_',SiteName_Data,'_',SiteName_Predict,'.mat'],'file'))
%             fprintf('calculating Weight Peak42\n');
%             Weight=sparse([II{:}],[JJ{:}],1./power([Distance{:}],2),N_Predict,N_Data);
%             save([DIRPATH,'SpatialLaggedWeightPeak42_',SiteName_Data,'_',SiteName_Predict,'.mat'],'Weight','-v7.3');
%         end
%         
%         %% version 4, power 3
%         if(~exist([DIRPATH,'SpatialLaggedWeightPeak43_',SiteName_Data,'_',SiteName_Predict,'.mat'],'file'))
%             fprintf('calculating Weight Peak43\n');
%             Weight=sparse([II{:}],[JJ{:}],1./power([Distance{:}],3),N_Predict,N_Data);
%             save([DIRPATH,'SpatialLaggedWeightPeak43_',SiteName_Data,'_',SiteName_Predict,'.mat'],'Weight','-v7.3');
%         end
        

    end  
    return;
end

N_Site = size(EnvironPara.SiteData_Model,1);

%% temporal weights
% calculate the average predicted values from nearby days.
N_Neighbor_Day = 9;%odd number
Weight_matrix_temporal_1 = ones(N_Neighbor_Day,1);
Weight_matrix_temporal_2 = [1:1:(N_Neighbor_Day-1)/2,0,flip(1:1:(N_Neighbor_Day-1)/2)]'; % linear template
Weight_matrix_temporal_3 = Weight_matrix_temporal_2.^2;

%% estimation-based temporal lagged terms
%%after getting better estimation, update those terms recursively
%%update estimation-based terms in disk
%%this is for updating in prediction session...large data size ahead!!!

%% update for temporal ones
% calculate the average predicted values from nearby days; it is used in
% prediction stage; since the nubmer of grid cells is so large, and we need
% to calculate such temporal averages by day
if(strcmp(OPTION,'UPDATE'))
    N_Day = 1;
    %% all have been predicted
    TempFileName = arrayfun(@(i)  ( [EnvironPara.DIRPATH_MODEL,'LaggedTerms',Sep,'TemporalLaggedTerms_',EnvironPara.NAME,'_',EnvironPara.SITENAME_PREDICT,'_',datestr(i,'yyyymmdd'),'_',datestr(i,'yyyymmdd'),'.mat']),EnvironPara.FIRSTDAY:EnvironPara.LASTDAY,'UniformOutput',false); 
    if(cellfun(@(x) exist(x,'file'),TempFileName))%% all have been processed
        fprintf('SKIP updating temporal lagged terms! ...%s..%s\n',TempFileName{1},TempFileName{length(TempFileName)});
    else
        mkdir([EnvironPara.DIRPATH_MODEL,'LaggedTerms',Sep]);

        %% initial set up
        DirPath = EnvironPara.OUTPUTPATH;
        FilePrefix = EnvironPara.PredictionResultName;
        LoadedVar = 'MonitorPredicted';

         %% for each day, update --tested on 2015-07-14
        %% load initial 9 days...
        CurrentDay = EnvironPara.FIRSTDAY;
        Temporal_Lagged = [];
        for i=1:N_Neighbor_Day
            TempDay = max(EnvironPara.FIRSTDAY,CurrentDay + i - (N_Neighbor_Day+1)/2);
            FileName = strrep(FilePrefix,'@@@@',datestr(TempDay,'yyyymmdd'));

            %%read next day data, wait if other prediction thread hasn't done
            for k = 1:10
                try
                    load([DirPath,FileName,'.mat'],LoadedVar);
                    break;
                catch
                    pause(EnvironPara.WAIT_TIME);
                    fprintf('\twaiting for prediction...%s! pause %d seconds...\n',FileName,EnvironPara.WAIT_TIME);                              
                end
            end
            if(10 == k)
                return;
            end

            Temporal_Lagged = cat(3,Temporal_Lagged,eval(LoadedVar));    

        end

        %%for every day -- estimation-based temporal lagged terms
        disp('updating estimation-based temporal lagged terms...');
        for i=EnvironPara.FIRSTDAY:EnvironPara.LASTDAY
            CurrentDay = datestr(i,'yyyymmdd');
            FileName = [EnvironPara.DIRPATH_MODEL,'LaggedTerms',Sep,'TemporalLaggedTerms_',EnvironPara.NAME,'_',EnvironPara.SITENAME_PREDICT,'_',CurrentDay,'_',CurrentDay];

            IsProcessFlag = false;
            try
                load([FileName,'.mat'],'Result');
                IsProcessFlag = false;
            catch
                if(exist([FileName,'.mat.part'],'file'))
                    list = dir([FileName,'.mat.part']);
                    if(etime(clock,datevec(list.datenum))>60)
                        delete([FileName,'.mat.part']);
                        fprintf('%s process work started long time ago and still has not completed. we do it again!...\n',FileName);
                        IsProcessFlag = true;
                    else
                        fprintf('%s is in the middle processing!...\n',FileName);
                        IsProcessFlag = false;
                    end
                else
                    IsProcessFlag = true;
                end
            end


            if(~IsProcessFlag)
                fprintf('\tthis day has been processed!%s\n',FileName);
            else
                save([FileName,'.mat.part'],'i');
    %             if(rem(i-EnvironPara.FIRSTDAY+1-EnvironPara.PredictionThread_ID,EnvironPara.PredictionThread)==0)
                if(true)
                    %%calculate and save
                    fprintf('\tprocessing Temporal Lagged terms...%s\n',FileName);

                    %%temporal lagged terms
                    Temporal_Lagged_1 = reshape(Analysis_MultipleMatrix(Temporal_Lagged,Weight_matrix_temporal_1,1),[N_Day*N_Site,1]);
                    Temporal_Lagged_2 = reshape(Analysis_MultipleMatrix(Temporal_Lagged,Weight_matrix_temporal_2,1),[N_Day*N_Site,1]);
                    Temporal_Lagged_3 = reshape(Analysis_MultipleMatrix(Temporal_Lagged,Weight_matrix_temporal_3,1),[N_Day*N_Site,1]);

                    %%save
                    Result = [Temporal_Lagged_1,Temporal_Lagged_2,Temporal_Lagged_3];
                    save([FileName,'.mat'],'Result');
                    delete([FileName,'.mat.part']);
                end
            end

            %%update temporal-lagged terms---remove old and add new day
            TempDay = min(i+(N_Neighbor_Day+1)/2,EnvironPara.LASTDAY);
            FileName = strrep(FilePrefix,'@@@@',datestr(TempDay,'yyyymmdd'));
            %%read next day data, wait if other prediction thread hasn't done
            for k = 1:10
                try
                    load([DirPath,FileName,'.mat'],LoadedVar);
                    break;
                catch
                    pause(EnvironPara.WAIT_TIME);
                    fprintf('\twaiting for prediction...%s! pause %d seconds...\n',FileName,EnvironPara.WAIT_TIME);                
                end
            end
            if(10==k)
                return;
            end
            Temporal_Lagged(:,:,1)=[];
            Temporal_Lagged = cat(3,Temporal_Lagged,eval(LoadedVar));
        end
    end

end

%% update for spatial ones
% calculate the average predicted values from nearby monitors; it is used in
% prediction stage; since the nubmer of grid cells is so large, and we need
% to calculate such spatial averages by day
if(strcmp(OPTION,'UPDATE'))  
% if(false)
    N_Day = 1;
    %% spatial lagged terms; the first one
    disp('updating estimation-based spatial lagged terms 1 (NO. 41)...');
    TempFileName = arrayfun(@(i)  ([EnvironPara.DIRPATH_MODEL,'LaggedTerms',Sep,'SpatialLaggedTerms1_',EnvironPara.NAME,'_',EnvironPara.SITENAME_PREDICT,'_',datestr(i,'yyyymmdd'),'_',datestr(i,'yyyymmdd'),'.mat']),EnvironPara.FIRSTDAY:EnvironPara.LASTDAY,'UniformOutput',false); 
    if(cellfun(@(x) exist(x,'file'),TempFileName))%% all have been processed
        fprintf('SKIP updating estimation-based spatial lagged terms 1 (NO. 41)! ...%s..%s\n',TempFileName{1},TempFileName{length(TempFileName)});
    else
        try
            ResultTemp = LoadData_function([DIRPATH,'SpatialLaggedWeightPeak41_',SiteName_Data,'_',SiteName_Predict,'.mat'],'Weight');
            Weight = ResultTemp.Weight;
        catch exception
            fprintf('error!..%s\n',exception.message);
            disp('perhaps weight matrix file has NOT been created! stupid mistake! Ah~~~');

            Analysis_TwoStep(Data,SiteData_Predict,SiteData_Data,SiteName_Predict,SiteName_Data,EnvironPara,'CreateWeight');
            ResultTemp = LoadData_function([DIRPATH,'SpatialLaggedWeightPeak41_',SiteName_Data,'_',SiteName_Predict,'.mat'],'Weight');
            Weight = ResultTemp.Weight;
        end

        %%calculate, process and save them for every day.
        i = 1;%EnvironPara.PredictionThread_ID;
        while i>=1 && i<=EnvironPara.LASTDAY-EnvironPara.FIRSTDAY+1
            CurrentDay = datestr(i+EnvironPara.FIRSTDAY-1,'yyyymmdd');
            FileName = [EnvironPara.DIRPATH_MODEL,'LaggedTerms',Sep,'SpatialLaggedTerms1_',EnvironPara.NAME,'_',EnvironPara.SITENAME_PREDICT,'_',CurrentDay,'_',CurrentDay];

            IsProcessFlag = false;
            try
                load([FileName,'.mat'],'Result');
                IsProcessFlag = false;
            catch
                if(exist([FileName,'.mat.part'],'file'))
                    list = dir([FileName,'.mat.part']);
                    if(etime(clock,datevec(list.datenum))>60)
                        delete([FileName,'.mat.part']);
                        fprintf('%s process work started long time ago and still has not completed. we do it again!...\n',FileName);
                        IsProcessFlag = true;
                    else
                        fprintf('%s is in the middle processing!...\n',FileName);
                        IsProcessFlag = false;
                    end
                else
                    IsProcessFlag = true;
                end
            end

            if(~IsProcessFlag)
                fprintf('\tthis day has been processed!%s\n',FileName);
            else
                save([FileName,'.mat.part'],'i');
                fprintf('\tprocessing spatial lagged terms; the first one...%s\n',FileName);
                Result= reshape(MultipleWeightMatrix_1(Weight',Data(i,:)')',[N_Day*N_Site,1]);
                save([FileName,'.mat'],'Result');
                delete([FileName,'.mat.part']);
            end
            i = i + 1;
        end
        clear Weight
    end
    
    %% spatial lagged terms; the second one
    disp('updating estimation-based spatial lagged terms 2 (NO.42)...');
    TempFileName = arrayfun(@(i)  ([EnvironPara.DIRPATH_MODEL,'LaggedTerms',Sep,'SpatialLaggedTerms2_',EnvironPara.NAME,'_',EnvironPara.SITENAME_PREDICT,'_',datestr(i,'yyyymmdd'),'_',datestr(i,'yyyymmdd'),'.mat']),EnvironPara.FIRSTDAY:EnvironPara.LASTDAY,'UniformOutput',false); 
    if(cellfun(@(x) exist(x,'file'),TempFileName))%% all have been processed
        fprintf('SKIP updating estimation-based spatial lagged terms 2 (NO. 42)! ...%s..%s\n',TempFileName{1},TempFileName{length(TempFileName)});
    else

        try
            ResultTemp = LoadData_function([DIRPATH,'SpatialLaggedWeightPeak42_',SiteName_Data,'_',SiteName_Predict,'.mat'],'Weight');
            Weight = ResultTemp.Weight;
        catch exception
            fprintf('error!..%s\n',exception.message);
            disp('perhaps weight matrix file has NOT been created! stupid mistake! Ah~~~');
            Analysis_TwoStep(Data,SiteData_Predict,SiteData_Data,SiteName_Predict,SiteName_Data,EnvironPara,'CreateWeight');
            ResultTemp = LoadData_function([DIRPATH,'SpatialLaggedWeightPeak42_',SiteName_Data,'_',SiteName_Predict,'.mat'],'Weight');
            Weight = ResultTemp.Weight;
        end

        %%calculate, process and save them for every day.
        i = 1;%EnvironPara.PredictionThread_ID;
        while i>=1 && i<=EnvironPara.LASTDAY-EnvironPara.FIRSTDAY+1
            CurrentDay = datestr(i+EnvironPara.FIRSTDAY-1,'yyyymmdd');
            FileName = [EnvironPara.DIRPATH_MODEL,'LaggedTerms',Sep,'SpatialLaggedTerms2_',EnvironPara.NAME,'_',EnvironPara.SITENAME_PREDICT,'_',CurrentDay,'_',CurrentDay];

            IsProcessFlag = false;
            try
                load([FileName,'.mat'],'Result');
                IsProcessFlag = false;
            catch
                if(exist([FileName,'.mat.part'],'file'))
                    list = dir([FileName,'.mat.part']);
                    if(etime(clock,datevec(list.datenum))>60)
                        delete([FileName,'.mat.part']);
                        fprintf('%s process work started long time ago and still has not completed. we do it again!...\n',FileName);
                        IsProcessFlag = true;
                    else
                        fprintf('%s is in the middle processing!...\n',FileName);
                        IsProcessFlag = false;
                    end
                else
                    IsProcessFlag = true;
                end
            end

            if(~IsProcessFlag)
                fprintf('\tthis day has been processed!%s\n',FileName);
            else
                save([FileName,'.mat.part'],'i');
                fprintf('\tprocessing spatial lagged terms; the second one...%s\n',FileName);
                Result= reshape(MultipleWeightMatrix_1(Weight',Data(i,:)')',[N_Day*N_Site,1]);
                save([FileName,'.mat'],'Result');
                delete([FileName,'.mat.part']);
            end
            i = i + 1;
        end
        clear Weight
    end
    
    %% spatial lagged terms; the third one
    disp('updating estimation-based spatial lagged terms 3 (NO. 43)...');
    TempFileName = arrayfun(@(i)  ([EnvironPara.DIRPATH_MODEL,'LaggedTerms',Sep,'SpatialLaggedTerms3_',EnvironPara.NAME,'_',EnvironPara.SITENAME_PREDICT,'_',datestr(i,'yyyymmdd'),'_',datestr(i,'yyyymmdd'),'.mat']),EnvironPara.FIRSTDAY:EnvironPara.LASTDAY,'UniformOutput',false); 
    if(cellfun(@(x) exist(x,'file'),TempFileName))%% all have been processed
        fprintf('SKIP updating estimation-based spatial lagged terms 1 (NO. 43)! ...%s..%s\n',TempFileName{1},TempFileName{length(TempFileName)});
    else

        try
            ResultTemp = LoadData_function([DIRPATH,'SpatialLaggedWeightPeak43_',SiteName_Data,'_',SiteName_Predict,'.mat'],'Weight');
            Weight = ResultTemp.Weight;
        catch exception
            fprintf('error!..%s\n',exception.message);
            disp('perhaps weight matrix file has NOT been created! stupid mistake! Ah~~~');
            Analysis_TwoStep(Data,SiteData_Predict,SiteData_Data,SiteName_Predict,SiteName_Data,EnvironPara,'CreateWeight');
            ResultTemp = LoadData_function([DIRPATH,'SpatialLaggedWeightPeak43_',SiteName_Data,'_',SiteName_Predict,'.mat'],'Weight');
            Weight = ResultTemp.Weight;
        end

        %%calculate, process and save them for every day.
        i = 1;
        while i>=1 && i<=EnvironPara.LASTDAY-EnvironPara.FIRSTDAY+1
            CurrentDay = datestr(i+EnvironPara.FIRSTDAY-1,'yyyymmdd');
            FileName = [EnvironPara.DIRPATH_MODEL,'LaggedTerms',Sep,'SpatialLaggedTerms3_',EnvironPara.NAME,'_',EnvironPara.SITENAME_PREDICT,'_',CurrentDay,'_',CurrentDay];

            IsProcessFlag = false;
            try
                load([FileName,'.mat'],'Result');
                IsProcessFlag = false;
            catch
                if(exist([FileName,'.mat.part'],'file'))
                    list = dir([FileName,'.mat.part']);
                    if(etime(clock,datevec(list.datenum))>60)
                        delete([FileName,'.mat.part']);
                        fprintf('%s process work started long time ago and still has not completed. we do it again!...\n',FileName);
                        IsProcessFlag = true;
                    else
                        fprintf('%s is in the middle processing!...\n',FileName);
                        IsProcessFlag = false;
                    end
                else
                    IsProcessFlag = true;
                end
            end

            if(~IsProcessFlag)
                fprintf('\tthis day has been processed!%s\n',FileName);
            else
                save([FileName,'.mat.part'],'i');
                fprintf('\tprocessing spatial lagged terms; the third one...%s\n',FileName);
                Result= reshape(MultipleWeightMatrix_1(Weight',Data(i,:)')',[N_Day*N_Site,1]);
                save([FileName,'.mat'],'Result');
                delete([FileName,'.mat.part']);
            end
            i = i + 1;
        end
        clear Weight   
    end
end

%% delete lagged terms for that year
% after the prediction is done, we need to delete these temporal and
% spatial averages
if(strcmp(OPTION,'DeleteTerms'))
    OutputTemplate = 'SpatialLaggedTerms1_$NAME$_$SiteName_Predict$_$STARTDATE$_$ENDDATE$.mat';
    OutputTemplate = strrep(strrep(OutputTemplate,'$NAME$',EnvironPara.NAME),'$SiteName_Predict$',EnvironPara.SITENAME_PREDICT);  
    TempFileName = arrayfun(@(i)  ([EnvironPara.DIRPATH_MODEL,'LaggedTerms',Sep, strrep(strrep(OutputTemplate,'$STARTDATE$',datestr(EnvironPara.FIRSTDAY+i-1,'yyyymmdd')),'$ENDDATE$',datestr(EnvironPara.FIRSTDAY+i-1,'yyyymmdd'))]),1:1:EnvironPara.LASTDAY-EnvironPara.FIRSTDAY+1,'UniformOutput',false); 
    for i=1:length(TempFileName)
        delete(TempFileName{i});
        fprintf('delete ...%s\n',TempFileName{i});
    end
    
    OutputTemplate = 'SpatialLaggedTerms2_$NAME$_$SiteName_Predict$_$STARTDATE$_$ENDDATE$.mat';
    OutputTemplate = strrep(strrep(OutputTemplate,'$NAME$',EnvironPara.NAME),'$SiteName_Predict$',EnvironPara.SITENAME_PREDICT);  
    TempFileName = arrayfun(@(i)  ([EnvironPara.DIRPATH_MODEL,'LaggedTerms',Sep, strrep(strrep(OutputTemplate,'$STARTDATE$',datestr(EnvironPara.FIRSTDAY+i-1,'yyyymmdd')),'$ENDDATE$',datestr(EnvironPara.FIRSTDAY+i-1,'yyyymmdd'))]),1:1:EnvironPara.LASTDAY-EnvironPara.FIRSTDAY+1,'UniformOutput',false); 
    for i=1:length(TempFileName)
        delete(TempFileName{i});
        fprintf('delete ...%s\n',TempFileName{i});
    end
    
    OutputTemplate = 'SpatialLaggedTerms3_$NAME$_$SiteName_Predict$_$STARTDATE$_$ENDDATE$.mat';
    OutputTemplate = strrep(strrep(OutputTemplate,'$NAME$',EnvironPara.NAME),'$SiteName_Predict$',EnvironPara.SITENAME_PREDICT);  
    TempFileName = arrayfun(@(i)  ([EnvironPara.DIRPATH_MODEL,'LaggedTerms',Sep, strrep(strrep(OutputTemplate,'$STARTDATE$',datestr(EnvironPara.FIRSTDAY+i-1,'yyyymmdd')),'$ENDDATE$',datestr(EnvironPara.FIRSTDAY+i-1,'yyyymmdd'))]),1:1:EnvironPara.LASTDAY-EnvironPara.FIRSTDAY+1,'UniformOutput',false); 
    for i=1:length(TempFileName)
        delete(TempFileName{i});
        fprintf('delete ...%s\n',TempFileName{i});
    end
    
    OutputTemplate = 'TemporalLaggedTerms_$NAME$_$SiteName_Predict$_$STARTDATE$_$ENDDATE$.mat';
    OutputTemplate = strrep(strrep(OutputTemplate,'$NAME$',EnvironPara.NAME),'$SiteName_Predict$',EnvironPara.SITENAME_PREDICT);  
    TempFileName = arrayfun(@(i)  ([EnvironPara.DIRPATH_MODEL,'LaggedTerms',Sep, strrep(strrep(OutputTemplate,'$STARTDATE$',datestr(EnvironPara.FIRSTDAY+i-1,'yyyymmdd')),'$ENDDATE$',datestr(EnvironPara.FIRSTDAY+i-1,'yyyymmdd'))]),1:1:EnvironPara.LASTDAY-EnvironPara.FIRSTDAY+1,'UniformOutput',false); 
    for i=1:length(TempFileName)
        delete(TempFileName{i});
        fprintf('delete ...%s\n',TempFileName{i});
    end
end


%% this is used in CV/All fitting session, data size is small
% calculate the average predicted values from nearby monitors and nearby days; it is used in
% model training stage; since the nubmer of grid cells is small, and we can
% calculate them altogther
if(~isempty(Data) && (strcmp(OPTION,'UPDATE_IN_TRAINING')))    
    disp('updating estimation-based terms...');
    
    disp('loading weight matrix files...');
    Weight_Cell = {};
    if(exist([DIRPATH,'SpatialLaggedWeightPeak41_',SiteName_Data,'_',SiteName_Predict,'.mat'],'file'))
        load([DIRPATH,'SpatialLaggedWeightPeak41_',SiteName_Data,'_',SiteName_Predict,'.mat'],'Weight');
        Weight_Cell = cat(1,Weight_Cell,{Weight});
    else
        disp('perhaps weight matrix file has NOT been created! stupid mistake! Ah~~~');
    end
    
    if(exist([DIRPATH,'SpatialLaggedWeightPeak42_',SiteName_Data,'_',SiteName_Predict,'.mat'],'file'))
        load([DIRPATH,'SpatialLaggedWeightPeak42_',SiteName_Data,'_',SiteName_Predict,'.mat'],'Weight');
        Weight_Cell = cat(1,Weight_Cell,{Weight});
    else
        disp('perhaps weight matrix file has NOT been created! stupid mistake! Ah~~~');
    end
    
    if(exist([DIRPATH,'SpatialLaggedWeightPeak43_',SiteName_Data,'_',SiteName_Predict,'.mat'],'file'))
        load([DIRPATH,'SpatialLaggedWeightPeak43_',SiteName_Data,'_',SiteName_Predict,'.mat'],'Weight');
        Weight_Cell = cat(1,Weight_Cell,{Weight});
    else
        disp('perhaps weight matrix file has NOT been created! stupid mistake! Ah~~~');
    end
    
    %% spatial estimation-based term
    Spatial_LaggedTerms = Analysis_CalculateSpatialAverage(Data,Weight_Cell);

    %% temporal ones
    N_Day = size(Data,1);
    Temporal_Lagged = nan([N_Day,N_Site,9]);
    
    Temporal_Lagged(:,:,1) = LagData(Data,4);
    Temporal_Lagged(:,:,2) = LagData(Data,3);
    Temporal_Lagged(:,:,3) = LagData(Data,2);
    Temporal_Lagged(:,:,4) = LagData(Data,1);
    Temporal_Lagged(:,:,5) = LagData(Data,0);
    Temporal_Lagged(:,:,6) = LagData(Data,-1);
    Temporal_Lagged(:,:,7) = LagData(Data,-2);
    Temporal_Lagged(:,:,8) = LagData(Data,-3);
    Temporal_Lagged(:,:,9) = LagData(Data,-4);
    
    Temporal_Lagged_1 = reshape(Analysis_MultipleMatrix(Temporal_Lagged,Weight_matrix_temporal_1,1),[N_Day*N_Site,1]);
    Temporal_Lagged_2 = reshape(Analysis_MultipleMatrix(Temporal_Lagged,Weight_matrix_temporal_2,1),[N_Day*N_Site,1]);
    Temporal_Lagged_3 = reshape(Analysis_MultipleMatrix(Temporal_Lagged,Weight_matrix_temporal_3,1),[N_Day*N_Site,1]);
    
    %% return results
    Result = [Spatial_LaggedTerms,Temporal_Lagged_1,Temporal_Lagged_2,Temporal_Lagged_3];
    
end

end

%% matrix multiply, M1*M2, excluding nan, normalization
% M1:Weight:[N_Site_Prediction,N_Site_Training] format; M2: training
% dataset; [N_Site_Training,N_Day] format;
% this file is an improved version of MultipleWeightMatrix_1
% this file
function Result = Analysis_MultipleMatrix(M1,M2,OPTION)

if(1 == OPTION)
% M1: N_Day*N_Site*100, M2: 100*1; Result: N_Day*N_Site
   SIZE = size(M1);
   Result = nan(SIZE(1),SIZE(2));
   for i=1:SIZE(1)
      for j=1:SIZE(2)
          %excluding nan and normalization
          Index = isnan(M1(i,j,:));
          Result(i,j) = nansum(reshape(M1(i,j,:),[SIZE(3),1]).*M2)/nansum(M2(~Index));
      end
   end
    
elseif(2 == OPTION)
% M1: N_Site*100, M2: 100*1; Result: N_Day*N_Site
    SIZE = size(M1);
    Result = nan(SIZE(1),1);
    for i=1:SIZE(1)
       Index = isnan(M1(i,:));
       Result(i,:) = nansum(reshape(M1(i,:),[SIZE(2),1]).*M2)/nansum(M2(~Index));
    end
elseif(3==OPTION)




end
    
end





