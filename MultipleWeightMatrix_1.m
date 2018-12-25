%% a function to do matrix multiplication while ignoring nan values
% weight is a m*n matrix; Y is a n*1 matrix; multiple them, with each rows
% in Weight been normalized(with sum to 1); but elements correspnding to
% NaN value in Y are not included in normalization
function Result = MultipleWeightMatrix_1(Weight,Y)

Result = nan(size(Weight,1),size(Y,2));
for i=1:size(Result,2)
    Result(:,i) = MultipleWeightMatrix(Weight,Y(:,i));
end

% %% new and fast
% NaNIndex = isnan(Y);
% Y(NaNIndex) = 0;
% Result = (Weight*Y)./sum(Weight(:,~NaNIndex),2);