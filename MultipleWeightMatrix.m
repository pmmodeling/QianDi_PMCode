% weight is a m*n matrix; Y is a n*1 matrix; multiple them, with each rows
% in Weight been normalized(with sum to 1); but elements correspnding to
% NaN value in Y are not included in normalization
function Result = MultipleWeightMatrix(Weight,Y)
SIZE = size(Weight);
Index = isnan(Y);% index of NaN value;
Y(Index) = 0;



% Weight_1 = Weight;
% Weight_1(:,Index) = 0;
% WeightTemp = Weight_1./repmat(sum(Weight_1,2),[1,SIZE(2)]);
% Result = WeightTemp*Y;

Result = Weight*Y;
TempSum = sum(Weight(:,~Index),2);
Result = Result./TempSum;