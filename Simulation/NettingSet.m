classdef NettingSet < handle
    %NETTINGSET Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % attributes
        name 
        attributes
        application
        %trade handles collection
        tradeRepository
        
        %pd curves for the entity and its counterparty
        pdCurveEnt
        pdCurveCtp
        
        %aggregated mtm profile, N * M (N = number of time grid, M = number
        %of path)
        mtmTimeGrid
        mtmProfile %not discounted
        mtmProfile_CSAAdjusted %not discounted
        mtmDisFactor
        CSAHandler
        
        %Aggregated cva/dva profile
        EPEProfile
        ENEProfile
        PFEProfileToEntity
        PFEProfileToCounterparty
        CVA
        DVA
        PFEToEntity
        PFEToCounterparty
        EPE
        ENE
        Exposure
        %Allocated cva/dva
        ICVA
        IDVA
    end
    
    methods
        %constructor
        function  obj = NettingSet(name, application)
            %initilaize trade repository
%             obj.tradeRepository = containers.Map('KeyType','double','ValueType','any');
            obj.tradeRepository = containers.Map();
            
            %parse key
            if iscell(name)
                name = char(name);
            end
            if ~isempty(name)
                obj.name = name;
                obj.attributes = yyStrSplit(name,'-');
                obj.application = application;
            else
                % set default value
                obj.attributes = {};
                obj.application = 'NONE-NONE';
            end %end if
        
        end
        
        %register trade into the netting set collection
        function c = setTrade(obj, tradeKey, tradeHandle)
            if ~isa(tradeHandle, 'Instrument')
                error('NettingSet: the trade is not an instrument');
            end
            if isnumeric(tradeKey)
                tradeKey = num2str(tradeKey);
            end
            
            obj.tradeRepository(tradeKey) = tradeHandle;
            
            if isa(tradeHandle.CSAHandler, 'CSAInfo')
                obj.CSAHandler = tradeHandle.CSAHandler;
            end
            
            attachNettingSet(tradeHandle, obj);
        end
 
        %get trade for a netting set
        function c = getTrade(obj)
                c = keys(obj.tradeRepository)';
        end
        
        %check if the attributes of the deal meet the netting set settings
        function c = checkEligibility(obj, varargin)
            
        end
        
        %aggregate the deal mtm profile to into the overall profile in the
        %netting set
        function [CVA, DVA, PFE, PNE, EPE, ENE, PFEProfile, PNEProfile, EPEProfile, ENEProfile, Exposure] = registerMTMProfile(obj, mtmTimeGrid, mtmProfile, mtmProfile_CSAAdjusted)
            if isempty(obj.mtmProfile)
                obj.mtmProfile = mtmProfile;
                obj.mtmTimeGrid = mtmTimeGrid;
                obj.mtmDisFactor = yyMtMDf(obj.mtmTimeGrid);
            elseif isempty(mtmProfile) || isempty(mtmTimeGrid)
                warning('NettingSet: no mtmProfile is registered, mtmProfile or the time grid is empty')
            else
                try
                    index = obj.mtmTimeGrid==mtmTimeGrid;
                catch ME
                    newMsg = message('NettingSet:time grid dimension does not match');
                    newME = MException(newMsg.Identifier,getString(newMsg));
                    newME = addCause(newME,ME);
                    throw(newME)
                end
                
                if sum(index)~=length(obj.mtmTimeGrid)
                    error('NettingSet: time grid not identical, cannot register mtmProfile')
                else
                    if size(obj.mtmProfile) == size(mtmProfile)
                        obj.mtmProfile = obj.mtmProfile + mtmProfile;                        
                    else
                        error('NettingSet: mtmProfile dimension does not match');
                    end
                end
            end
            
            %PERFORM STANDALON CALCULATION
            disF = obj.mtmDisFactor;
            PD = CumPDByTime(obj.pdCurveCtp,obj.mtmTimeGrid);
            r = obj.pdCurveCtp.recoveryRate;
            [CVA, PFE, EPE, EPEProfile, PFEProfile,Exposure] = yyNSCalculate(mtmProfile_CSAAdjusted, disF, PD, r, 1);
            
            
            PD = CumPDByTime(obj.pdCurveEnt,obj.mtmTimeGrid);
            r = obj.pdCurveEnt.recoveryRate;
            [DVA, PNE,ENE,ENEProfile, PNEProfile,~] = yyNSCalculate(mtmProfile_CSAAdjusted, disF, PD, r, -1);
        end
        
        %register PDcurve
        function c = registerPDCurve(obj, EntName, CtpName, env)
            tmp = getPDCurveHandle(env,EntName);
            if ~isempty(tmp)
                obj.pdCurveEnt = tmp;
            end
            tmp = getPDCurveHandle(env,CtpName);
            if ~isempty(tmp)
                obj.pdCurveCtp = tmp;
            end
        end
        
        %perform EPE, ENE, CVA, DVA, PFE calculation
        function [] = calculation(obj, pricingGrid)
            %current setting: row, time grid;column, path;
            if isempty(obj.mtmProfile)
                error('NettingSet: mtmProfile is empty');
            else
%                 MtM = obj.mtmProfile;
                disF = obj.mtmDisFactor;
            end
            nPath = size(obj.mtmProfile, 2);
            %No of pricing gird point
            nPricingGrid = length(pricingGrid);
            if strcmp(obj.application, 'ALL-ALL')
                %apply CSA
                %deal with non-csa support deal first, just do projection using
                %linear interpolation
                %No of simulated paths
                nPath = size(obj.mtmProfile, 2);
                %No of pricing gird point
                nPricingGrid = length(pricingGrid);

                if ~isa(obj.CSAHandler, 'CSAInfo')
                    obj.mtmProfile_CSAAdjusted = zeros(nPricingGrid, nPath);
                     for iPath = 1 : nPath
                        mtmPath = obj.mtmProfile(:, iPath);
                        obj.mtmProfile_CSAAdjusted(:, iPath) = interp1q(obj.mtmTimeGrid, mtmPath, pricingGrid);
                     end
                else
                     obj.mtmProfile_CSAAdjusted = applyCSA(obj.CSAHandler, obj.mtmProfile, ...
                        obj.mtmTimeGrid, pricingGrid, false);
                end

                %calculate cva
                if isempty(obj.pdCurveCtp)
                    error('NettingSet: PD curve for the counterparty is empty');
                else
                    PD = CumPDByTime(obj.pdCurveCtp,obj.mtmTimeGrid);
                    r = obj.pdCurveCtp.recoveryRate;
                    [obj.CVA, obj.PFEToEntity,obj.EPE, obj.EPEProfile, obj.PFEProfileToEntity,obj.Exposure] = ...
                        yyNSCalculate(obj.mtmProfile_CSAAdjusted, disF, PD, r, 1);
                end

                if isempty(obj.pdCurveEnt)
                    error('NettingSet: PD curve for the entity is empty');
                else
                    PD = CumPDByTime(obj.pdCurveEnt,obj.mtmTimeGrid);
                    r = obj.pdCurveEnt.recoveryRate;
                    [obj.DVA, obj.PFEToCounterparty,obj.ENE, obj.ENEProfile, obj.PFEProfileToCounterparty,~] = ...
                        yyNSCalculate(obj.mtmProfile_CSAAdjusted, disF, PD, r, -1);    
                end
            else
                %dummy netting set, composed of non netted deals
                keyDeals = obj.tradeRepository.keys;
                nbDeals = length(keyDeals);
                obj.CVA = 0;
                obj.DVA = 0;
                obj.EPEProfile = zeros(nPricingGrid,1);
                obj.ENEProfile = zeros(nPricingGrid,1);
                obj.PFEProfileToEntity =zeros(nPricingGrid,1);
                obj.PFEProfileToCounterparty = zeros(nPricingGrid,1);
                obj.Exposure = 0;                
                for i = 1:nbDeals
                    deal = obj.tradeRepository(char(keyDeals(i)));
                    obj.CVA = obj.CVA + deal.standaloneCVA;
                    obj.DVA = obj.DVA + deal.standaloneDVA;                    
                    obj.EPEProfile = obj.EPEProfile + deal.standaloneEPEProfile;
                    obj.ENEProfile = obj.ENEProfile + deal.standaloneENEProfile;
                    obj.PFEProfileToEntity = deal.standalonePFEProfile;
                    obj.PFEProfileToCounterparty = obj.PFEProfileToCounterparty + deal.standalonePNEProfile;
                    obj.Exposure = obj.Exposure + deal.standaloneExposure;                           
                end
                obj.PFEToEntity = max(obj.PFEProfileToEntity);
                obj.PFEToCounterparty = max(obj.PFEProfileToCounterparty);
                obj.EPE = mean(obj.EPEProfile);
                obj.ENE = mean(obj.ENEProfile);
            end
        end   
   
        
        %perform CVA, DVA allocation
        function [ICVA,IDVA] = allocation(obj)
            if isempty(obj.CVA) || isempty(obj.DVA)
                error('NettingSet: The CVA/DVA must be calculated before allocation')
            else
                env = Environment.getInstance;
                tradeKeys = obj.tradeRepository.keys;
                ICVA = cell(length(tradeKeys), 2);
                IDVA = cell(length(tradeKeys), 2);
                
                for i = 1:length(tradeKeys)
                    [~,tradeMtM] = pathPricer(obj.tradeRepository(char(tradeKeys(i))),env.valuationDate,env,obj.mtmTimeGrid);
                    mtmProfile = obj.mtmProfile - tradeMtM;

                    %Incremental CVA
                    pdCtp = CumPDByTime(obj.pdCurveCtp,obj.mtmTimeGrid);
                    ICVAi = obj.CVA - yyNSCalculate(mtmProfile,obj.mtmDisFactor, pdCtp, obj.pdCurveCtp.recoveryRate, 1);
                    ICVA(i,:) = {tradeKeys(i) ICVAi};
                    
                    %Incremental DVA
                    pdEnt = CumPDByTime(obj.pdCurveEnt,obj.mtmTimeGrid);
                    IDVAi = obj.DVA - yyNSCalculate(mtmProfile,obj.mtmDisFactor, pdEnt, obj.pdCurveEnt.recoveryRate, -1);
                    IDVA(i,:) = {tradeKeys(i) IDVAi};
                end
                obj.ICVA = ICVA;
                obj.IDVA = IDVA;
            end
        end     
        
    end  
end

