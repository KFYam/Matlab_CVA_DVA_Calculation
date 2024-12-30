classdef CSAInfo < handle
    %CSAINFO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        counterPartyName %name of the counterparty
        CSAType % bilectaral unilecteral
        CSACurrency %currency of the colateral paid
        thresholdI % payment threshhold for us
        thresholdC % threshold of counterparty
        minTransAmtI % mininum transfer amount for us
        minTransAmtC % mininum transfer amount for us
        indAmt % independent amount
        marginPeriod  % marginPeriod in number of days
        existingCollateral % collateral already posted
    end
    
    methods
        function obj = CSAInfo(CPName, CSAType, ccy, thresholdI, thresholdC, minTransAmtI, minTransAmtC, indAmt, marginPeriod, existingCollateral)
            if strcmp(CSAType, 'Bilateral') == false && strcmp(CSAType, 'Unilateral (Reciever)') == false ...
                && strcmp(CSAType, 'Unilateral (Payer)') == false
                error(message('instrument:csa type can only be bilateral or unilateral'));
            else
                obj.CSAType = CSAType;
            end
            
            obj.counterPartyName = CPName;
            obj.CSACurrency = ccy;
            
            if isnumeric(thresholdI) == false || thresholdI < 0
                error(message('instrument:threshold must not be nagative'));
            else
                obj.thresholdI = thresholdI;
            end
            
            if isnumeric(thresholdC) == false || thresholdC < 0
                error(message('instrument:threshold must not be nagative'));
            else
                obj.thresholdC = thresholdC;
            end
            
            if isnumeric(minTransAmtI) == false || minTransAmtI < 0
                error(message('instrument:minTransAmt must not be nagative'));
            else
                obj.minTransAmtI = minTransAmtI;
            end
            
            if isnumeric(minTransAmtC) == false || minTransAmtC < 0
                error(message('instrument:minTransAmt must not be nagative'));
            else
                obj.minTransAmtC = minTransAmtC;
            end
            
            if isnumeric(indAmt) == false 
                error(message('instrument:indAmt must be numerical'));
            else
                obj.indAmt = indAmt;
            end
            
            if isnumeric(marginPeriod) == false 
                error(message('instrument:marginPeriod must be numerical'));
            else
                obj.marginPeriod =  marginPeriod;
            end
            
            if isnumeric(existingCollateral) == false 
                error(message('instrument:existingCollateral must be numerical'));
            else
                obj.existingCollateral =  existingCollateral;
            end
        end
        
        function results = applyCSA(obj, mtmProfile, mtmTimeGrid, pricingGrid, applyExistingCollateral)
           
            %No of simulated paths
            nPath = size(mtmProfile, 2);
            %No of pricing gird point
            nPricingGrid = length(pricingGrid);
            %result matrix
            results = zeros(nPricingGrid, nPath);                
            %margin period in year fraction
            mPeriod = obj.marginPeriod / 365;
            
            %check the boundary of simulated profile
            if mtmTimeGrid(end) < pricingGrid(end)
                error(message('Instrument:Pricing grid is beyond simulated grid'));
            end
            
            %check the size of mtm profile and simulation grid size
            if size(mtmProfile, 1) ~= length(mtmTimeGrid)
                error(message('Instrument:simulation grid and mtm profile size doesnt match'));
            end
            
            %margin Grid, if the time if less than value date, adjust it to
            %zero
            marginGrid = max(pricingGrid - mPeriod, 0);
            
            %loop through each path
            for iPath = 1 : nPath
                mtmPath = mtmProfile(:, iPath);
                %substract independent amount
                mtmPath = mtmPath - obj.indAmt;
                %get MTM at beginning of margin period
                mtmA = interp1q(mtmTimeGrid, mtmPath, marginGrid);
                %get MTM at pricing grid
                mtmB = interp1q(mtmTimeGrid, mtmPath, pricingGrid);
                
                %calculate colateral posted at the beginning of the margin
                %risk period
                if strcmp(obj.CSAType, 'Bilateral')
                    %bilateral CSA
                    collateralProfile = max(mtmA - (obj.thresholdI + obj.minTransAmtI), 0) - ...
                        max(-mtmA - (obj.thresholdC + obj.minTransAmtC), 0);
                elseif strcmp(obj.CSAType, 'Unilateral (Reciever)')
                    %Unilateral CSA reciever
                    collateralProfile = max(mtmA - (obj.thresholdI + obj.minTransAmtI), 0);
                else
                    %Unilateral CSA payer
                    collateralProfile = - max(-mtmA - (obj.thresholdC + obj.minTransAmtC), 0);
                end
                
                %apply initial collateral amount
                if marginGrid(1) == 0 && applyExistingCollateral == true
                    collateralProfile(1) = obj.existingCollateral;
                end
                
                % Calculate CSA adjusted exposure at pricing grid
                CSAAdjustedMTM = mtmB - collateralProfile;
                
                results(:, iPath) = CSAAdjustedMTM;
            end %end of path loop
           
        end
    end
    
end

