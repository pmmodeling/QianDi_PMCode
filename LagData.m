%%  get n-day lagged of the data
% missing days are filled by averages
% Day>0 xx-Day lag; 
% Day<0 future values

%% paramaters
% Data: input data; N_Day*N_Site format;
% Day: number of lagged days

%% code
function Result = LagData(Data,Day)

SIZE = size(Data);
Result = nan(SIZE);

if(Day>0)
    Result(1+Day:SIZE(1),:) = Data(1:SIZE(1)-Day,:);
    Result(1:1+Day-1,:)=repmat(Data(1,:),[Day,1]);
elseif(Day<0)
    Result(1:SIZE(1)+Day,:) = Data(1-Day:SIZE(1),:);
    Result(SIZE(1)+Day+1:SIZE(1),:) = repmat(Data(SIZE(1),:),[abs(Day),1]);
else
   Result = Data; 
end