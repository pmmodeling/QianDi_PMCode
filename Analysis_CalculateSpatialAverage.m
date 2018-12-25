%% calculate estimation-based spatial weight
%%return is in [N_Site*N_Day,1] format
function Result = Analysis_CalculateSpatialAverage(Data,Weight_Cell)
SIZE = size(Data);

N_Weight = length(Weight_Cell);
Result = nan(SIZE(1)*SIZE(2),N_Weight);

for i=1:N_Weight
    Weight = Weight_Cell{i};
    Result(:,i)= reshape(MultipleWeightMatrix_1(Weight,Data')',[SIZE(1)*SIZE(2),1]);
end
