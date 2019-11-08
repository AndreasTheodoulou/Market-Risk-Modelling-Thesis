%% Value-at-Risk and Expected Shortfall Estimation and Backtesting 
%
% This code estimates the value-at-risk (VaR) and Expected Shortfall(ES)
%using the following methods, and performs a VaR and ES backtesting 
%analysis. The methods are:

% # Historical simulation 
% # Exponential weighted moving average (EWMA) - Normal Distribution
% # GARCH(1,1) - Normal Distribution
% # GJR - Normal Distribution
% # Exponential weighted moving average (EWMA) - t-Distribution
% # GARCH(1,1) - t-Distribution
% # GJR - t-Distribution
% 
% The estimation methods used, estimate the VaR at 95% and 99% 
% confidence levels.

%% Start Function
function [] = VaRestimates(filename, estimationyear, estimationmonth,PortfolioID)

%filename - Data to be imported from excel, 1 column of dates and 1 column
%of closing prices.
% estimation year - year on which the VaR and ES estimations begin.
%estimation month - month on which the VaR and ES estimations begin.

%% Load the Data and Define the Test Window

% Get daily close and date data
[~,~ , rawData ]= xlsread (filename,'','','basic');
date_cell = rawData (2:length(rawData),1);
close_cell = rawData (2:length(rawData),2);
date = zeros (1,length(date_cell));
close = zeros (1,length(close_cell));

% Convert the data to string format
for i =1: length (date)
date(i) = datenum(cell2mat(date_cell(i)));
close_i = cell2mat(close_cell(i));
if iscellstr(close_i)
close (i)= str2num (close_i);
else
close (i)= close_i;
end
end
date = fliplr(date);
% Convert the date format to a correct one
DateReturns = datestr(x2mdate(date));
DateReturns = flipud(DateReturns);

Returns = zeros(length(DateReturns),1);
for i= 2:length(DateReturns)
Returns(i)= log(close(i)/ close(i-1));
end

%Define samplesize
SampleSize = length(Returns); 

%%
% Define the estimation window as 250 trading days for HS and 1000 for
% parametric methods. For parametric methods the estimation window is used
% to calculate the coefficients of the GARCH family models.

%%
TestWindowStart      = find(year(DateReturns)==estimationyear & month(DateReturns)==estimationmonth,1);
TestWindow           = TestWindowStart : SampleSize;
EstimationWindowSizeHS = 250;
EstimationWindowSizeGarch = 1000;

%%
% For a VaR confidence level of 95% and 99%, set the complement of the VaR level.
q = [0.05 0.01];

%% Compute the VaR Using the Historical Simulation Method

%Create the vector for the Historical Simulation VaR
VaR_Historical95 = zeros(length(TestWindow),1);
VaR_Historical99 = zeros(length(TestWindow),1);
ES_Historical95 = zeros(length(TestWindow),1);
ES_Historical99 = zeros(length(TestWindow),1);

%run iterations and save each of the length(TestWindow) Historical Simulation VaR
for t = TestWindow
    i = t - TestWindowStart + 1;
    EstimationWindow = t-EstimationWindowSizeHS:t;
    X = Returns(EstimationWindow);
    %VaR estimates
    VaR_Historical95(i) = -quantile(X,q(1)); 
    VaR_Historical99(i) = -quantile(X,q(2)); 
    %ES estimates method
    N = length(X);  
    k1 = ceil(N*(1-q(1))); 
    k2 = ceil(N*(1-q(2)));
    z = sort(X);
      
     if k1 < N
       ES95 = ((k1 - N*(1-q(1)))*z(k1) + sum(z(k1 +1:N)))/(N*q(1));
    else
       ES95 = z(k1);
     end
     
    if k2 < N
       ES99 = ((k2 - N*(1-q(2)))*z(k2) + sum(z(k2+1:N)))/(N*q(2));
    else
       ES99 = z(k2);
     end
    %ES estimate
    ES_Historical95(i) = ES95;
    ES_Historical99(i) = ES99;
     
end

%% Compute VaR using the Exponential Weighted Moving Average(EWMA) model.

%Set lambda = 0.94
lambda = 0.94;
Zscore   = norminv(q);

%Create EwmaN and GwmaT models
EwmaN = garch ('GARCHLags',1,'ARCHLags',1,'Distribution','Gaussian','Constant',1e-322 ,'ARCH',(1-lambda)-(1e-14) ,'GARCH',lambda-(1e-14));
EwmaT = garch ('GARCHLags',1,'ARCHLags',1,'Distribution','t','Constant',1e-322,'ARCH',(1-lambda)-(1e-14),'GARCH',lambda -(1e-14));

%  Create vectors for GarchN and GarchT VaR estimates, sigma2, Degrees of
%  Freedom.
VaR_EwmaN95 = zeros(length(TestWindow),1);
VaR_EwmaN99 = zeros(length(TestWindow),1);
VaR_EwmaT95 = zeros(length(TestWindow),1);
VaR_EwmaT99 = zeros(length(TestWindow),1);

ES_EwmaN95 = zeros(length(TestWindow),1);
ES_EwmaN99 = zeros(length(TestWindow),1);
ES_EwmaT95 = zeros(length(TestWindow),1);
ES_EwmaT99 = zeros(length(TestWindow),1);

Sigma2EwmaN = zeros(length(TestWindow),1);
Sigma2EwmaT = zeros(length(TestWindow),1);
DoFs1 = zeros(length(TestWindow),1);

% Estimate coefficients every day and compute GARCHN VaR estimates
for t= TestWindow 
    %create a variable that starts from 1
    i = t - TestWindowStart + 1;
    
% Vector with relevant returns
EstimationWindow = t-EstimationWindowSizeGarch:t;
R_t = Returns(EstimationWindow);
% Estimate model coefficients
estEwmaN = estimate(EwmaN,R_t,'Display','off');
estEwmaT = estimate(EwmaT,R_t, 'Display','off');

% Forecast volatility
Sigma2EwmaN(i)= forecast(estEwmaN,1,'Y0',R_t);
SigmaEwmaN = sqrt(Sigma2EwmaN(i));

Sigma2EwmaT(i) = forecast(estEwmaT,1,'Y0',R_t);
SigmaEwmaT = sqrt(Sigma2EwmaT(i));

%degrees of freedom
DoF = estEwmaT.Distribution.DoF;
DoFs1(i) = DoF;
%VaR estimates
    VaR_EwmaN95(i) = -Zscore(1)*SigmaEwmaN;
    VaR_EwmaN99(i) = -Zscore(2)*SigmaEwmaN;
    VaR_EwmaT95(i) = -Zscore(1)*SigmaEwmaT;
    VaR_EwmaT99(i) = -Zscore(2)*SigmaEwmaT;

%ES estimates
   ES_EwmaN95(i) = SigmaEwmaN*(normpdf(Zscore(1))/q(1));
   ES_EwmaN99(i) = SigmaEwmaN*(normpdf(Zscore(2))/q(2));
   ES_EwmaT95(i) = SigmaEwmaT*(tpdf(tinv(1-q(1),DoF),DoF)/(q(1))*((DoF + tinv(1-q(1),DoF)^2)/(DoF-1))*sqrt((DoF-2)/DoF));
   ES_EwmaT99(i) = SigmaEwmaT*(tpdf(tinv(1-q(2),DoF),DoF)/(q(2))*((DoF + tinv(1-q(2),DoF)^2)/(DoF-1))*sqrt((DoF-2)/DoF));
    
end

%% Compute VaR using the GARCH(1,1).

%Create GarchN and GarchT model
GarchN = garch('GARCHLags', 1, 'ARCHLags', 1, 'Distribution', 'Gaussian' );
GarchT = garch('GARCHLags', 1, 'ARCHLags', 1, 'Distribution', 't' );

%  Create vectors for GarchN and GarchT VaR estimates, sigma2, Degrees of
%  Freedom
VaR_GarchN95 = zeros(length(TestWindow),1);
VaR_GarchN99 = zeros(length(TestWindow),1);
VaR_GarchT95 = zeros(length(TestWindow),1);
VaR_GarchT99 = zeros(length(TestWindow),1);

ES_GarchN95 = zeros(length(TestWindow),1);
ES_GarchN99 = zeros(length(TestWindow),1);
ES_GarchT95 = zeros(length(TestWindow),1);
ES_GarchT99 = zeros(length(TestWindow),1);

Sigma2GarchN = zeros(length(TestWindow),1);
Sigma2GarchT = zeros(length(TestWindow),1);
DoFs2 = zeros(length(TestWindow),1);

% Estimate coefficients every day and compute GARCHN VaR estimates
for t= TestWindow 
    %create a variable that starts from 1
    i = t - TestWindowStart +1;
% Vector with relevant returns
EstimationWindow = t-EstimationWindowSizeGarch:t;
R_t = Returns(EstimationWindow);
% Estimate model coefficients
estGarchN = estimate(GarchN,R_t,'Display','off');
estGarchT = estimate(GarchT,R_t, 'Display','off');

% Forecast volatility and compute standard deviation
Sigma2GarchN(i)= forecast(estGarchN,1,'Y0',R_t);
SigmaGarchN = sqrt(Sigma2GarchN(i));

Sigma2GarchT(i) = forecast(estGarchT,1,'Y0',R_t);
SigmaGarchT = sqrt(Sigma2GarchT(i));

%store degrees of freedom 
DoF = estGarchT.Distribution.DoF;
DoFs2(i) = DoF;

%VaR estimates
    VaR_GarchN95(i) = -SigmaGarchN*Zscore(1);
    VaR_GarchN99(i) = -SigmaGarchN*Zscore(2);
    VaR_GarchT95(i) = -SigmaGarchT*(tinv(q(1),DoF)*sqrt((DoF-2)/DoF));
    VaR_GarchT99(i) = -SigmaGarchT*(tinv(q(2),DoF)*sqrt((DoF-2)/DoF));
    
%ES estimates
   ES_GarchN95(i) = SigmaGarchN*(normpdf(Zscore(1))/q(1));
   ES_GarchN99(i) = SigmaGarchN*(normpdf(Zscore(2))/q(2));
   ES_GarchT95(i) = SigmaGarchT*(tpdf(tinv(1-q(1),DoF),DoF)/(q(1))*((DoF + tinv(1-q(1),DoF)^2)/(DoF-1))*sqrt((DoF-2)/DoF));
   ES_GarchT99(i) = SigmaGarchT*(tpdf(tinv(1-q(2),DoF),DoF)/(q(2))*((DoF + tinv(1-q(2),DoF)^2)/(DoF-1))*sqrt((DoF-2)/DoF));
   
end

%% Compute VaR using the GJR(1,1).

%Create GarchN and GarchT model
GjrN = gjr('GARCHLags',1,'ARCHLags',1,'LeverageLags',1,'Distribution', 'Gaussian' );
GjrT = gjr('GARCHLags',1,'ARCHLags',1,'LeverageLags',1,'Distribution','t' );

%  Create vectors for GarchN and GarchT VaR estimates, ES estimates, sigma2, Degrees of
%  Freedom
VaR_GjrN95 = zeros(length(TestWindow),1);
VaR_GjrN99 = zeros(length(TestWindow),1);
VaR_GjrT95 = zeros(length(TestWindow),1);
VaR_GjrT99 = zeros(length(TestWindow),1);

ES_GjrN95 = zeros(length(TestWindow),1);
ES_GjrN99 = zeros(length(TestWindow),1);
ES_GjrT95 = zeros(length(TestWindow),1);
ES_GjrT99 = zeros(length(TestWindow),1);

Sigma2GjrN = zeros(length(TestWindow),1);
Sigma2GjrT = zeros(length(TestWindow),1);
DoFs3 = zeros(length(TestWindow),1);

% Estimate coefficients every day and compute GARCHN VaR estimates
for t= TestWindow 
    %create a variable that starts from index no.1
    i = t - TestWindowStart + 1;
% Vector with relevant returns
EstimationWindow = t-EstimationWindowSizeGarch:t;
R_t = Returns(EstimationWindow);
% Estimate model coefficients
estGjrN = estimate(GjrN,R_t,'Display','off');
estGjrT = estimate(GjrT,R_t, 'Display','off');
% Forecast volatility
Sigma2GjrN(i)= forecast(estGjrN,1,'Y0',R_t);
SigmaGjrN = sqrt(Sigma2GjrN(i));

Sigma2GjrT(i) = forecast(estGjrT,1,'Y0',R_t);
SigmaGjrT = sqrt(Sigma2GjrT(i));

%degrees of freedom
DoF = estGarchT.Distribution.DoF;
DoFs3(i) = DoF;

%VaR estimates
    VaR_GjrN95(i) = -Zscore(1)*SigmaGjrN;
    VaR_GjrN99(i) = -Zscore(2)*SigmaGjrN;
    VaR_GjrT95(i) = -Zscore(1)*SigmaGjrT;
    VaR_GjrT99(i) = -Zscore(2)*SigmaGjrT;
    
%ES estimates
   ES_GjrN95(i) = SigmaGjrN*(normpdf(Zscore(1))/q(1));
   ES_GjrN99(i) = SigmaGjrN*(normpdf(Zscore(2))/q(2));
   ES_GjrT95(i) = SigmaGjrT*(tpdf(tinv(1-q(1),DoF),DoF)/(q(1))*((DoF + tinv(1-q(1),DoF)^2)/(DoF-1))*sqrt((DoF-2)/DoF));
   ES_GjrT99(i) = SigmaGjrT*(tpdf(tinv(1-q(2),DoF),DoF)/(q(2))*((DoF + tinv(1-q(2),DoF)^2)/(DoF-1))*sqrt((DoF-2)/DoF));
    
end

%% VaR Backtesting
%Set up the returns and dates vector for our the time window to be
%backtested.
ReturnsTest = Returns(TestWindow);
DatesTest   = DateReturns(TestWindow);

%Name VaR data, models
VaRData = [VaR_Historical95 VaR_EwmaN95 VaR_EwmaT95 VaR_GarchN95 VaR_GarchT95 VaR_GjrN95 VaR_GjrT95 VaR_Historical99 VaR_EwmaN99 VaR_EwmaT99 VaR_GarchN99 VaR_GarchT99 VaR_GjrN99 ...
    VaR_GjrT99];
Models = {'HS95','EwmaN95','EwmaT95','GarchN95','GarchT95', 'GjrN95', 'GjrT95' ...
    'HS99','EwmaN99','EwmaT99', 'GarchN99','GarchT99', 'GjrN99', 'GjrT99'};
VaRLevel = [0.95 0.95 0.95 0.95 0.95 0.95 0.95 0.99 0.99 0.99 0.99 0.99 0.99 0.99];

% Run the Unconditional Coverage backtest on our methods.

%Set up the backtesting results
vbt = varbacktest(ReturnsTest,VaRData,'PortfolioID',PortfolioID,'VaRID',Models,'VaRLevel',VaRLevel);
%Report the Unconditional coverage backtesting result
pof(vbt)

%% ES Backtesting

% Run the Unconditional backtest on our methods.

rng('default'); % for reproducibility; the esbacktestbysim constructor runs a simulation

%Name VaR data, ES data, and degrees of freedom. 

ESData = [ES_Historical95 ES_EwmaN95 ES_EwmaT95 ES_GarchN95 ES_GarchT95 ES_GjrN95 ES_GjrT95 ES_Historical99 ES_EwmaN99 ES_EwmaT99 ES_GarchN99 ES_GarchT99 ES_GjrN99 ...
    ES_GjrT99];

DoFs = [round(mean(DoFs1)), round(mean(DoFs2)), round(mean(DoFs3)),round(mean(DoFs1)), round(mean(DoFs2)), round(mean(DoFs3))];


%set up the ES backtesting results
VaRDataT = [VaR_EwmaT95 VaR_GarchT95 VaR_GjrT95 VaR_EwmaT99 VaR_GarchT99 VaR_GjrT99];
ESDataT = [ES_EwmaT95 ES_GarchT95 ES_GjrT95 ES_EwmaT99 ES_GarchT99 ES_GjrT99];
ModelsT = {'EwmaT95', 'GarchT95', 'GjrT95', 'EwmaT99', 'GarchT99', 'GjrT99'};
VaRLevelT = [0.95 0.95 0.95 0.99 0.99 0.99];

VaRDataN = [VaR_EwmaN95 VaR_GarchN95 VaR_GjrN95 VaR_EwmaN99 VaR_GarchN99 VaR_GjrN99];
ESDataN = [ES_EwmaN95 ES_GarchN95 ES_GjrN95 ES_EwmaN99 ES_GarchN99 ES_GjrN99];
ModelsN = {'EwmaN95', 'GarchN95', 'GjrN95', 'EwmaN99', 'GarchN99', 'GjrN99'};
VaRLevelN = [0.95 0.95 0.95 0.99 0.99 0.99];

VaRDataHS = [VaR_Historical95 VaR_Historical99];
ESDataHS = [ES_Historical95 ES_Historical99];
ModelsHS = {'HS95','HS99'};
VaRLevelHS = [0.95 0.99];

ebtsT = esbacktestbysim(ReturnsTest,VaRDataT,ESDataT,"t",...
       'DegreesOfFreedom',10,...
       'Location',0,...
       'Scale',1,...
       'PortfolioID',PortfolioID,...
       'VaRID',ModelsT,...
       'VaRLevel',VaRLevelT);
   
ebtsN = esbacktestbysim(ReturnsTest,VaRDataN,ESDataN,"normal",...
       'Mean',0,...
       'StandardDeviation',1,...
       'PortfolioID',PortfolioID,...
       'VaRID',ModelsN,...
       'VaRLevel',VaRLevelN);  
   
ebtsHS = esbacktest(ReturnsTest,VaRDataHS,ESDataHS,'PortfolioID',PortfolioID,'VaRID',ModelsHS,'VaRLevel',VaRLevelHS);     
   %Report the Unconditional test 
   
unconditional(ebtsT)
unconditional(ebtsN)
unconditionalNormal(ebtsHS)

%% Compute Sample Statistics for the Returns
Mean = mean(Returns);
Standard_Deviation = std(Returns);
Kurtosis = kurtosis(Returns);
Skewness = skewness(Returns);
Min = min(Returns);
Max = max(Returns);

SampleStatistics = table(Mean,Standard_Deviation ,Kurtosis,Skewness,Min,Max)
end
