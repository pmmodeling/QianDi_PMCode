%% This is the main function to process raw data one by one to 
% obtain aggregate file. 
% raw data have to be download to ./raw_data/data/unprocessed/
% aggregate data will be produce here: ./raw_data/data/aggregate/
% remember to upload grid cell location files of aggregate data; check
% readme file for "file dependence"

%% Parameter:
% YEAR: which year of data you would like to process?

%% code
function Run_AggregateDataSummary(YEAR)

% process restaurant data
IsRestaurant = true;

% process MERRA 2 data
IsMerra2 = true;

% process MOD11A1 data
IsMOD11A1  = true;

% process MOD04L2 data
IsMOD04L2 = true;

% process MAIACAqua, MAIACTerra
IsMAIACAqua = true;
IsMAIACTerra = true;

% process MODIS data
IsMOD13A2 = true;
IsMOD09A1 = true;

% process OMI data
IsOMAERUVd = true;
IsOMAEROe = true;
IsOMNO2d = true;
IsOMSO2e = true;
IsOMTO3e = true;
IsOMUVBd = true;

% OMI L2 data
IsOMO3PR = true;

% CMAQ data
IsCMAQ = true;

% process CMAQ data
IsForestFire = true;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% beging to process
%% process restaurant data
if(IsRestaurant)
    try
        if(YEAR >= 2000 && YEAR <= 2016)
            Read_RestaurantData(YEAR);
        end
    catch exception
        fprintf('error, %s\n',exception.message); 
    end
end

%% process MERRA 2 data
if(IsMerra2)
    try
        if(YEAR >= 2000 && YEAR <= 2016)
            ReadMERRA2aer(YEAR);
        end
    catch exception
        fprintf('error, %s\n',exception.message); 
    end
end

%% process MOD11A1 data
if(IsMOD11A1)
    try
        if(YEAR >= 2000 && YEAR <= 2016)
            ReadMOD11A1_Main(YEAR);
        end
    catch exception
        fprintf('error, %s\n',exception.message); 
    end
end

%% process MOD04L2 data
if(IsMOD04L2)
    try
        if(YEAR >= 2000 && YEAR <= 2016)
            ReadMOD04L2_Main(YEAR,'MOD04L2');
        end
    catch exception
        fprintf('error, %s\n',exception.message); 
    end
end

%% process MAIAC Aqua data
if(IsMAIACAqua)
    try
        if(YEAR >= 2002 && YEAR <= 2016)
            ReadMAIAC_Main(YEAR,'Aqua');
        end
    catch exception
        fprintf('error, %s\n',exception.message); 
    end
end

%% process MAIAC Terra data
if(IsMAIACTerra)
    try
        if(YEAR >= 2000 && YEAR <= 2016)
            ReadMAIAC_Main(YEAR,'Terra');
        end
    catch exception
        fprintf('error, %s\n',exception.message); 
    end
end

%% process MOD013A2 data
if(IsMOD13A2)
    try
        if(YEAR >= 2000 && YEAR <= 2016)
            ReadMODIS_Main(YEAR,9);
        end
    catch exception
        fprintf('error, %s\n',exception.message); 
    end
end

%% process MOD09A1 data
if(IsMOD09A1)
    try
        if(YEAR >= 2000 && YEAR <= 2016)
            ReadMODIS_Main(YEAR,10);
        end
    catch exception
        fprintf('error, %s\n',exception.message); 
    end
end

%% process OMAERUVd data
if(IsOMAERUVd)
    try
        if(YEAR >= 2005 && YEAR <= 2016)
            ReadOMIAIData_Main(YEAR,'OMAERUVd');
        end
    catch exception
        fprintf('error, %s\n',exception.message); 
    end
end

%% process OMAEROe data
if(IsOMAEROe)
    try
        if(YEAR >= 2005 && YEAR <= 2016)
            ReadOMIAIData_Main(YEAR,'OMAEROe');
        end
    catch exception
        fprintf('error, %s\n',exception.message); 
    end
end

%% process OMNO2d data
if(IsOMNO2d)
    try
        if(YEAR >= 2005 && YEAR <= 2016)
            ReadOMIAIData_Main(YEAR,'OMNO2d');
        end
    catch exception
        fprintf('error, %s\n',exception.message); 
    end
end

%% process OMSO2e data
if(IsOMSO2e)
    try
        if(YEAR >= 2005 && YEAR <= 2016)
            ReadOMIAIData_Main(YEAR,'OMSO2e');
        end
    catch exception
        fprintf('error, %s\n',exception.message); 
    end
end

%% process OMTO3e data
if(IsOMTO3e)
    try
        if(YEAR >= 2005 && YEAR <= 2016)
            ReadOMIAIData_Main(YEAR,'OMTO3e');
        end
    catch exception
        fprintf('error, %s\n',exception.message); 
    end
end

%% process OMUVBd data
if(IsOMUVBd)
    try
        if(YEAR >= 2005 && YEAR <= 2016)
            ReadOMIAIData_Main(YEAR,'OMUVBd');
        end
    catch exception
        fprintf('error, %s\n',exception.message); 
    end
end

%% process OMO3PR data
if(IsOMO3PR)
    try
        if(YEAR >= 2005 && YEAR <= 2016)
            ReadOMO3PR_Main(YEAR);
        end
    catch exception
        fprintf('error, %s\n',exception.message); 
    end
end

%% process CMAQ data
if(IsCMAQ)
    try
        if(YEAR >= 2002 && YEAR <= 2006)
            ReadCMAQSite(YEAR,'A');
            ReadCMAQSite(YEAR,'B');
        elseif(YEAR <2002 && YEAR>= 2000)
            ReadCMAQSite(YEAR,'');
        elseif(YEAR >2006 && YEAR <= 2014)
            ReadCMAQSite(YEAR,'');
        end
    catch exception
        fprintf('error, %s\n',exception.message); 
    end
end

%% process forest fire data
if(IsForestFire)
    try
        if(YEAR >= 2000 && YEAR <= 2016)
            ReadForestfire(YEAR);
        end
    catch exception
        fprintf('error, %s\n',exception.message); 
    end
end
