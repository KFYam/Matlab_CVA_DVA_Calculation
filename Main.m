%% Market Environment Development
addpath(genpath(pwd));
timeLog = {};

% load raw data
tic
[fileName,pathName] = uigetfile('*.*','Select the input data file');
fprintf('Loading Source File\n')
raw = importdata([pathName,fileName]);
t = toc;
timeLog = [timeLog;{'load raw data' num2str(t/60)}];

% read market data and model construction
tic
mainData = raw.textdata.InputFilePath;
val_date = datenum(mainData(2,2));
env = Environment.getInstance;
clearEnvironment(env);
env.valuationDate = val_date;

[simEng,pricingGrid] = CreateModels(raw);
t = toc;
timeLog = [timeLog;{'load market data' num2str(t/60)}];

% read portfolio data and counterparty information
tic
[fwdDeals,swpDeals] = readPortfolio(raw);
t = toc;
timeLog = [timeLog;{'load portfolio data' num2str(t/60)}];

%% Market Scenarios Simulation
tic
display('Performing simulation.....');
simulate(simEng);
t = toc;
timeLog = [timeLog;{'simulation'} {num2str(t/60)}];

%% Portfolio Mark to Market Profiles Valuation
% tic
% display('Performing deal valuation.....');
allDeals = [fwdDeals;swpDeals];
keyFwd = fwdDeals.keys;
keySwp = swpDeals.keys;
keyAll = [keyFwd'; keySwp'];
%mtmProfile = containers.Map();
notPriced = {}; % log down not priced positions

bar = waitbar(0,'Instruments Pricing','Name','Pricing');
nbDeals = length(keyAll);
for i = 1:nbDeals
    waitbar(i/nbDeals,bar,sprintf('%d Instrument priced',i))
    try
        [~,mtm] = pathPricer(allDeals(char(keyAll(i))),val_date,env,[0;pricingGrid]);
        %mtmProfile(char(keyAll(i))) = mtm;
        display(strcat('Deal ',char(keyAll(i)), ' calculated'));
    catch err
        display(strcat('Deal ',char(keyAll(i)), ' goes wrong: ', err.message));
        notPriced = [notPriced;{i} keyAll(i)];
    end    
end
close(bar)
t = toc;
timeLog = [timeLog;{'pricing' num2str(t/60)}];

%% CVA/DVA Calculation and Allocation
tic
display('Performing CVA/DVA calulation.....');
keyNS = env.NettingSetCollection.keys;
NSCVA = [];
NSDVA = [];
NSPFE = [];
NSEPE = [];
NSPNE = [];
NSENE = [];
NSE = [];
for i = 1:length(keyNS)
    try
        nettingSet = getNettingSetHandle(env,char(keyNS(i)));
        calculation(nettingSet, [0;pricingGrid]);
        display(strcat(char(keyNS(i)), ' calculated'));
        NSCVA = [NSCVA;nettingSet.CVA];
        NSDVA = [NSDVA;nettingSet.DVA];
        NSPFE = [NSPFE;nettingSet.PFEToEntity];
        NSPNE = [NSPNE;nettingSet.PFEToCounterparty];
        NSEPE = [NSEPE;nettingSet.EPE];
        NSENE = [NSENE;nettingSet.ENE];
        NSE = [NSE;nettingSet.Exposure];
%         exposure = nettingSet.EPEProfile - nettingSet.ENEProfile;
%         NSE = [NSE;max(exposure)]; 
    catch err
        display(strcat(char(keyNS(i)), ' goes wrong: ', err.message));
        NSCVA = [NSCVA;0];
        NSDVA = [NSDVA;0];
        NSPFE = [NSPFE;0];
        NSPNE = [NSPNE;0];
        NSEPE = [NSEPE;0];
        NSENE = [NSENE;0];
        NSE = [NSE;0];
    end
end  
t = toc;
timeLog = [timeLog;{'CVA calculation' num2str(t/60)}];

%%Printing results in xls file
filename = strcat(pathName,'\OutputFile.xlsx');
%by netting set
xlswrite(filename,{'NettingSet'},1,'A1');
xlswrite(filename,keyNS',1,'A2');
xlswrite(filename,{'EPE'},1,'B1');
xlswrite(filename,NSEPE,1,'B2');
xlswrite(filename,{'PFE'},1,'C1');
xlswrite(filename,NSPFE,1,'C2');
xlswrite(filename,{'CVA'},1,'D1');
xlswrite(filename,NSCVA,1,'D2');
xlswrite(filename,{'ENE'},1,'E1');
xlswrite(filename,NSENE,1,'E2');
xlswrite(filename,{'PNE'},1,'F1');
xlswrite(filename,NSPNE,1,'F2');
xlswrite(filename,{'DVA'},1,'G1');
xlswrite(filename,NSDVA,1,'G2');
xlswrite(filename,{'Peak Exposure'},1,'H1');
xlswrite(filename,NSE,1,'H2');

%by deal
dealMtM = [];
dealCVA = [];
dealDVA = [];
dealPFE = [];
dealEPE = [];
dealPNE = [];
dealENE = [];
dealExposure = [];
for i = 1:length(keyAll)
    deal = allDeals(char(keyAll(i)));
    dealMtM = [dealMtM;deal.npv];
    dealCVA = [dealCVA;deal.standaloneCVA];
    dealDVA = [dealDVA;deal.standaloneDVA];
    dealPFE = [dealPFE;deal.standalonePFE];
    dealEPE = [dealEPE;deal.standaloneEPE];
    dealPNE = [dealPNE;deal.standalonePNE];
    dealENE = [dealENE;deal.standaloneENE];
    dealExposure = [dealExposure;deal.standaloneExposure];
    
end

xlswrite(filename,{'Deal'},2,'A1');
xlswrite(filename,keyAll,2,'A2');
xlswrite(filename,{'MtM'},2,'B1');
xlswrite(filename,dealMtM,2,'B2');
xlswrite(filename,{'EPE'},2,'C1');
xlswrite(filename,dealEPE,2,'C2');
xlswrite(filename,{'PFE'},2,'D1');
xlswrite(filename,dealPFE,2,'D2');
xlswrite(filename,{'CVA'},2,'E1');
xlswrite(filename,dealCVA,2,'E2');
xlswrite(filename,{'ENE'},2,'F1');
xlswrite(filename,dealENE,2,'F2');
xlswrite(filename,{'PNE'},2,'G1');
xlswrite(filename,dealPNE,2,'G2');
xlswrite(filename,{'DVA'},2,'H1');
xlswrite(filename,dealDVA,2,'H2');
xlswrite(filename,{'Peak Exposure'},2,'I1');
xlswrite(filename,dealExposure,2,'I2');

%% CVA/DVA deal level allocation
% tic
% display('Performing CVA/DVA calulation.....');
% for i = 1:length(keyNS)
%     try
%         nettingSet = getNettingSetHandle(env,char(keyNS(i)));
%         [ICVA,IDVA] = allocation(getNettingSetHandle(env,char(keyNS(i))));
%     catch err
%         display(strcat(char(keyNS(i)), ' goes wrong'));
%     end
% end 
% t = toc;

