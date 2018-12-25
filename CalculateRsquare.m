%% calculate out-of-sample R square; temporal R square; spatial R square
% Assuming Data1 and Data2 have similar structure with rows standing for
% date and column standing for sites.
% R2 is computed from regression Data2 = Y, Data1 = X;

%% Parameter
% Data1: data set 1: N*M format; N is the number of days and M is the
% number of monitoring site; 
% Data: data set 2: the same format;
% option: what kind of R2 to calculate? 'out-of-sample': total R2;
% 'spatial': spatial R2;'temporal': temporal R2;'each-temporal': temporal 
% R2 for each monitoring site; 'each-spatial': spatial R2 for each
% monitoring site

%% return value:
% R2: R2;
% pvalue: associated p-value;
% MSE: mean square error

function [R2,pvalue,MSE] = CalculateRsquare(Data1,Data2,option)

SIZE1 = size(Data1);
SIZE2 = size(Data2);

if(SIZE1(1)~=SIZE2(1) || SIZE1(2)~=SIZE2(2))
   fprintf('inconsistency matrix input in CalculateRsquare!');
   R2=-9999;
   pvalue= -9999;
   MSE=-9999;
   return;
end

if(strcmp(option,'out-of-sample'))
    
    Temp1 = reshape(Data1,[SIZE1(1)*SIZE1(2),1]);
    Temp2 = reshape(Data2,[SIZE1(1)*SIZE1(2),1]);
    
    %%% check if there are all nan values
    if(sum(isnan(Temp1))==length(Temp1)||sum(isnan(Temp2))==length(Temp2))
       fprintf('all nan values!\n');
       R2 = -9999;
       pvalue = -9999;
       MSE = -9999;
       return;
    end
    
    [b,bint,r,rint,stats] = regress(Temp2,[ones(size(Temp1,1),1),Temp1]);
    
    R2 = stats(1);
    pvalue = stats(3);
    MSE = sqrt(nanmean((Temp1-Temp2).^2));

elseif(strcmp(option,'each-spatial'))
    % calculate R2 in each site
    R2 = nan(1,size(Data1,2));
    pvalue = nan(1,size(Data1,2));
    MSE = nan(1,size(Data1,2));
    
    for i = 1:size(Data1,2)
        [R2_,pvalue_,MSE_] = CalculateRsquare(Data1(:,i),Data2(:,i),'out-of-sample');
        R2(i) = R2_;
        pvalue(i) = pvalue_;
        MSE(i) = MSE_;
    end
    
elseif(strcmp(option,'spatial'))
    Index = isnan(Data1)|isnan(Data2);
    Data1(Index) = nan;
    Data2(Index) = nan;
    
    Temp1 = nanmean(Data1)';
    Temp2 = nanmean(Data2)';
    
    %%% check if there are all nan values
    if(sum(isnan(Temp1))==length(Temp1)||sum(isnan(Temp1))==length(Temp2))
       fprintf('all nan values!\n');
       R2 = NaN;
       pvalue = NaN;
       MSE = NaN;
       return;
    end
    
    [b,bint,r,rint,stats] = regress(Temp2,[ones(size(Temp1,1),1),Temp1]);
    R2 = stats(1);
    pvalue = stats(3);
    MSE = sqrt(nanmean((Temp1-Temp2).^2));

elseif(strcmp(option,'each-temporal'))
    % R2 in each day
    R2 = nan(size(Data1,1),1);
    pvalue = nan(size(Data1,1),1);
    MSE = nan(size(Data1,1),1);
    
    for i = 1:size(Data1,1)
        [R2_,pvalue_,MSE_] = CalculateRsquare(Data1(i,:),Data2(i,:),'out-of-sample');
        R2(i) = R2_;
        pvalue(i) = pvalue_;
        MSE(i) = MSE_;
    end
    
elseif(strcmp(option,'temporal'))  
    
    Temp1 = nanmean(Data1);
    Temp2 = nanmean(Data2);
    
    %%% check if there are all nan values
    if(sum(isnan(Temp1))==length(Temp1)||sum(isnan(Temp1))==length(Temp2))
       fprintf('all nan values!\n');
       R2 = NaN;
       pvalue = NaN;
       MSE = NaN;
       return;
    end
    
    Data1_temp = Data1 - repmat(Temp1,[size(Data1,1),1]);
    Data2_temp = Data2 - repmat(Temp2,[size(Data2,1),1]);
    
    Temp1 = reshape(Data1_temp,[size(Data1_temp,1)*size(Data1_temp,2),1]);
    Temp2 = reshape(Data2_temp,[size(Data2_temp,1)*size(Data2_temp,2),1]);
    
    Index = ~isnan(Temp1)&~isnan(Temp2);
    Temp1 = Temp1(Index);
    Temp2 = Temp2(Index);
    
    [b,bint,r,rint,stats] = regress(Temp2,[ones(size(Temp1,1),1),Temp1]);
    R2 = stats(1);
    pvalue = stats(3);
    MSE = sqrt(nanmean((Temp1-Temp2).^2));

elseif(strcmp(option,'overall'))  
    Temp1 = reshape(Data1,[size(Data1,1)*size(Data1,2),1]);
    Temp2 = reshape(Data2,[size(Data2,1)*size(Data2,2),1]);
    
    [b,bint,r,rint,stats] = regress(Temp2,[ones(size(Temp1,1),1),Temp1]);
    R2 = stats(1);
    pvalue = stats(3);  
else
   fprintf('this option is not supported in CalculateRsquare!');  
   return
end
