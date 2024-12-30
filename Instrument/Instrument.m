classdef (Abstract) Instrument < handle
    %INSTRUMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = protected, GetAccess = public)
        npv 
        standaloneCVA
        standaloneDVA
        standalonePFE
        standaloneEPE
        standaloneENE
        standalonePNE
        standaloneExposure
        standaloneEPEProfile
        standaloneENEProfile
        standalonePFEProfile
        standalonePNEProfile
        %CSA Related
%         CSAType % bilectaral unilecteral
%         thresholdI % payment threshhold for us
%         thresholdC % threshold of counterparty
%         minTransAmt % mininum transfer amount
%         indAmt % independent amount
%         marginPeriod  % marginPeriod in number of days
        CSAHandler
        
        %Netting set collection
        nettingSets
    end
    
%     properties (SetAccess = protected, GetAccess = protected)
%         modelAssignment %pricing engine
%     end
    
    methods
        
        function obj = Instrument()
            obj.nettingSets = containers.Map();
            obj.npv = 0;
        end
        
        function c = get.npv(obj)
            c = obj.npv;
        end
        
        % attah netting set to the instrument
        function attachNettingSet(obj, nettingSet)
            if ~isa(nettingSet, 'NettingSet')
                error('addNettingSet:this is not a netting set object');
            end
            obj.nettingSets(nettingSet.name) = nettingSet;
        end
        
        % attach CSA 
        function c = attachCSA(obj, CSAHandler)
            if ~isa(CSAHandler, 'CSAInfo')
                error('attachCSA:this is not a CSA Info object');
            end
            obj.CSAHandler = CSAHandler;
        end
        
        % apply CSA to simulated mtm profile and project onto pricing grid
        function [pricingGrid, results] = applyCSA(obj, mtmGrid, mtmProfile, pricingGrid)
               
            %deal with non-csa support deal first, just do projection using
            %linear interpolation
            if ~isa(obj.CSAHandler, 'CSAInfo')
                %No of simulated paths
                nPath = size(mtmProfile, 2);
                %No of pricing gird point
                nPricingGrid = length(pricingGrid);
                results = zeros(nPricingGrid, nPath);
                for iPath = 1 : nPath
                   mtmPath = mtmProfile(:, iPath);
                   results(:, iPath) = interp1q(mtmGrid, mtmPath, pricingGrid);
                end
            else
                results = applyCSA(obj.CSAHandler, mtmProfile, mtmGrid, pricingGrid, false);
            end
                
            %register profile
            if obj.nettingSets ~= 0 
                ns = values(obj.nettingSets);          
                for i = 1:length(ns)                   
                    if strcmp(ns{i}.application,'NONE-NONE')
                        %we save the epe and pfe profile for PFE and EPE
                        %calculation when no netting applied
                        [obj.standaloneCVA, obj.standaloneDVA, obj.standalonePFE,obj.standalonePNE,obj.standaloneEPE,...
                        obj.standaloneENE,obj.standalonePFEProfile,obj.standalonePNEProfile,...
                        obj.standaloneEPEProfile,obj.standaloneENEProfile,obj.standaloneExposure] = ...
                        registerMTMProfile(ns{i}, pricingGrid, mtmProfile, results);
                    else
                        [obj.standaloneCVA, obj.standaloneDVA, obj.standalonePFE,obj.standalonePNE,obj.standaloneEPE,...
                                    obj.standaloneENE,~,~,~,~,obj.standaloneExposure] = ...
                                registerMTMProfile(ns{i}, pricingGrid, mtmProfile, results);
                    end
                    
                end
            end                    
        end %end of apply CSA
        
    end
    
    methods (Abstract)
        calculate(obj, valueDate, environment)
    end
    
end

