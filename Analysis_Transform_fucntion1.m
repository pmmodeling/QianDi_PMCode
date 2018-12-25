%% transform the data to avoid extreme values
function Output = Analysis_Transform_fucntion1(Input,DataOption,TransformOption)

% %% for testing only
% Input = load('D:\Google Drive\Research\USTemperature\DataPaper3\EPATemperature\Monitor\MONITOR_MeanTemperature_EPATemperature_1999_2015.mat');
% Input = Input.Result;
% Input = reshape(Input,[size(Input,1)*size(Input,2),1]);
% DataOption = 'MeanTemperature';
% TransformOption = 'transform';

%% code
if(any(strfind(DataOption,'MeanTemperature')))
    if(nanmean(nanmean(Input))<100)
        xmin = 250-273.15;% lower limit
        xmax = 340-273.15;% upper limit
        X_20 = 278.27-273.15;%the 20th percentile of the data
        X_80 = 296.11-273.15;%the 80th percentile of the data
        X_mean = 287.01-273.15;% the mean value of the data
    else
        xmin = 250;% lower limit
        xmax = 340;% upper limit
        X_20 = 278.27;%the 20th percentile of the data
        X_80 = 296.11;%the 80th percentile of the data
        X_mean = 287.01;% the mean value of the data
    end
elseif(any(strfind(DataOption,'MaxTemperature')))
    if(nanmean(nanmean(Input))<100)
        xmin = 250-273.15;
        xmax = 360-273.15;
        X_20 = 282.59-273.15;
        X_80 = 302.04-273.15;
        X_mean = 292.31-273.15;
    else
        xmin = 250;
        xmax = 360;
        X_20 = 282.59;
        X_80 = 302.04;
        X_mean = 292.31;
    end
elseif(any(strfind(DataOption,'MinTemperature')))
    if(nanmean(nanmean(Input))<100)
        xmin = 250-273.15;% lower limit
        xmax = 340-273.15;% upper limit
        X_20 = 273.71-273.15;%the 20th percentile of the data
        X_80 = 290.93-273.15;%the 80th percentile of the data
        X_mean = 282.09-273.15;% the mean value of the data
    else
        xmin = 250;% lower limit
        xmax = 340;% upper limit
        X_20 = 273.71;%the 20th percentile of the data
        X_80 = 290.93;%the 80th percentile of the data
        X_mean = 282.09;% the mean value of the data
    end
elseif(any(strfind(DataOption,'PM25')))
    xmin = 0;% lower limit
    xmax = 500;% upper limit
    X_20 = 5.1271;%the 20th percentile of the data
    X_80 = 15.2000;%the 80th percentile of the data
    X_mean = 9;% the mean value of the data
elseif(any(strfind(DataOption,'NO2')))
    xmin = 0;% lower limit
    xmax = 700;% upper limit
    X_20 = 8;%the 20th percentile of the data
    X_80 = 35;%the 80th percentile of the data
    X_mean = 22.5;% the mean value of the data
elseif(any(strfind(DataOption,'Ozone')))
    xmin = 0;% lower limit
    xmax = 200;% upper limit
    X_20 = 30;%the 20th percentile of the data
    X_80 = 55;%the 80th percentile of the data
    X_mean = 40;% the mean value of the data
elseif(any(strfind(DataOption,'soilm')))
    xmin = 0;% lower limit
    xmax = 1200;% upper limit
    X_20 = 387.7927;%the 20th percentile of the data
    X_80 = 604.6529;%the 80th percentile of the data
    X_mean = 494.1644;% the mean value of the data
elseif(any(strfind(DataOption,'MOD11A1_LST_Day_1km')))
    xmin = 230;% lower limit
    xmax = 360;% upper limit
    X_20 = 288.5600;%the 20th percentile of the data
    X_80 = 311.1400;%the 80th percentile of the data
    X_mean = 299.5988;% the mean value of the data
    
elseif(any(strfind(DataOption,'MOD11A1_LST_Night')))
    xmin = 230;% lower limit
    xmax = 340;% upper limit
    X_20 = 274.6150;%the 20th percentile of the data
    X_80 = 292.5900;%the 80th percentile of the data
    X_mean = 283.4293;% the mean value of the data
    
elseif(any(strfind(DataOption,'MOD11A1_Clear_day')))
    xmin = 0;% lower limit; scaling factor should be 0.0005; but I made a mistake when processing the data
    xmax = 1;% upper limit
    X_20 = 0.1936;%the 20th percentile of the data
    X_80 = 0.3275;%the 80th percentile of the data
    X_mean = 0.2311;% the mean value of the data
    
elseif(any(strfind(DataOption,'MOD11A1_Clear_night')))
    xmin = 0;% lower limit; scaling factor should be 0.0005; but I made a mistake when processing the data
    xmax = 1;% upper limit
    X_20 = 0.1779;%the 20th percentile of the data
    X_80 = 0.3361;%the 80th percentile of the data
    X_mean = 0.2308;% the mean value of the data
end

k = (X_80 - X_20)/2/1.0986;

if(strcmp('transform',TransformOption))%change monitored to modeled values
    %fprintf('transform data...%s;min:%f,20th:%f,mean:%f,80th:%f,max:%f\n',DataOption,xmin,X_20,X_mean,X_80,xmax);
    if(any(Input>xmax) || any(Input<xmin))
       disp(['error!!! the value of ',DataOption,' exceeds the valid range']); 
    end
    Input(Input>xmax) = nan;
    Input(Input<xmin) = nan;
    
    Output = X_mean + k*atanh(2*(Input - (xmin+xmax)/2)/(xmax-xmin));

    
elseif(strcmp('inverse',TransformOption))% change modeled to real monitored
    Output = (xmax-xmin)/2*tanh((Input - X_mean)/k)+(xmin+xmax)/2;
    %fprintf('inverse transform data...%s;min:%f,20th:%f,mean:%f,80th:%f,max:%f\n',DataOption,nanmin(Output),quantile(Output,0.2),nanmean(Output),quantile(Output,0.8),nanmax(Output));
else
    disp('input option error!');
    Output = nan;
    return;
end


% a = 250;
% b = 340;
% x: value used for modeling; y: real monitored data
% x = 250:0.1:340;
% y = (b-a)/2*tanh((x- 282)/7.8383)+(a+b)/2;
% plot(x,y);
%
% % x = 250:0.1:340
% % y = tanh((x- 282)/7.8383);
% % plot(x,y);
