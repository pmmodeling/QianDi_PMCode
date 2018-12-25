%% return the index of rows with complete data

%% Parameter
% Data: input data matrix

%% return value
% the index that data matrix has complete values

%% code
function Index = isnan_matrix(Data)
SIZE = size(Data);
Index = true(SIZE(1),1);

% go over each column and decide whether it has any missing value
for i=1:SIZE(2)
   TempIndex = ~isnan(Data(:,i));
   Index = Index&TempIndex;
end

% change to logical values
Index = logical(Index);

