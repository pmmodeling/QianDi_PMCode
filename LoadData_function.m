%% load matlab file and return error file if loading fails
% sometime the files on super computer corrputed, and we neeed to remove
% the corrupted file. So, this file load matlab file and remove corrupted
% file if loading fails

%% Parameter
% FilePath: the full file path of the matlab file
% varargin: parameters to be passed to "load"

%% return value
% Data: loaded value

%% code
function Data = LoadData_function(FilePath,varargin)
try
    if(isempty(varargin))
        Data = load(FilePath);
    else
        Data = load(FilePath,varargin{1});
    end
catch exception
    if(exist(FilePath,'file'))%file must be corrupt
        delete(FilePath);
        error(['file corrupted...',FilePath]);
    else
        error(['file not existed...',FilePath]);    
    end
end
