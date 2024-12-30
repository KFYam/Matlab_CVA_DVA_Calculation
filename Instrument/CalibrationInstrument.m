classdef (Abstract) CalibrationInstrument < handle
    %INSTRUMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = protected, GetAccess = public)
        marketPrice %        
        maturityYF
        modelPrice
    end
    properties (SetAccess = public, GetAccess = public)
        weight %determines the importance of the instrument in the calibration, btw 0 and 1        
    end
    
    properties (SetAccess = protected, GetAccess = protected)
        modelAssignment %pricing engine
    end
    
    methods
        function c = getMarketPrice(obj)
            c = obj.marketPrice;
        end
        function obj = setWeight(obj, w)
            obj.weight = w;
        end
    end
    
    methods (Abstract)
        MarketPrice(obj)
        ModelPrice(obj,inCurve,a,vSigma,vYFVol)
        
    end
    
end

