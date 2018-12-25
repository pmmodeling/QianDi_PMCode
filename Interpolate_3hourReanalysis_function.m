%% a local function called while reading meteorological data remove missing, filling, max, min values from meteorological files
% the valid value range is taken from the meteorological file headers;
% this file is NOT supposed to be used in a standalone way

%% Parameter
% FileInfo: file header of the nc file
% TempData: raw meteorological value extracted from the nc file
% VarName: the short name of meteorological variable name

%% Return value
% TempData: meteorological value with filled value, missing values replaced
% with "NaN"

%% code
function TempData =  Interpolate_3hourReanalysis_function(FileInfo,TempData,VarName)

Index1 = find(strcmp({FileInfo.Variables.Name},VarName));

%% find out filling and missing values
Index2 = find(strcmp({FileInfo.Variables(Index1).Attributes.Name},'missing_value'));
MissingValue = FileInfo.Variables(Index1).Attributes(Index2).Value;
Index2 = find(strcmp({FileInfo.Variables(Index1).Attributes.Name},'_FillValue'));
FillingValue = FileInfo.Variables(Index1).Attributes(Index2).Value;
Index2 = find(strcmp({FileInfo.Variables(Index1).Attributes.Name},'valid_range'));
ValidRange1 = FileInfo.Variables(Index1).Attributes(Index2).Value(1);
ValidRange2 = FileInfo.Variables(Index1).Attributes(Index2).Value(2);

fprintf('removing Fill Value...%d\n',FillingValue);
fprintf('removing missing Value...%d\n',MissingValue);
fprintf('removing values outside valid range...%d %d\n',ValidRange1,ValidRange2);
%% remove filling and missing values
TempData(TempData == MissingValue) = nan;
TempData(TempData == FillingValue) = nan;
TempData(TempData<ValidRange1) = nan;
TempData(TempData>ValidRange2) = nan;
