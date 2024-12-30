classdef HazardRateTermStructure < TermStructure
    %HAZARDRATETERMSTRUCTURE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dates %dates array
        hazardRates
        defaultProbs
        times %times array
        recoveryRate    
        %currency only support piecewise constant interpolation
    end
    
    methods
        function obj = HazardRateTermStructure(refDate,dayCounter, dates, hazardRates, recoveryRate)
            if nargin < 4 
                error(message('HazardRateTermStructure:Too few parameters'));
            end
            
            obj = obj@TermStructure(dayCounter, refDate);            
            obj.dates = dates;
            obj.rates = hazardRates;
            obj.recoveryRate = recoveryRate;
            obj.times = timeFromReference(obj, dates);
        end
        
        function c = defautProbability(obj, d)
            c = defaultProbabilityByTime(timeFromReference(d));
        end
        
        function c = defautProbabilityByTime(obj, t)
        end
    end
    
end

