classdef YieldTermStructure < TermStructure
    %YIELDTERMSTRUCTURE Summary of this class goes here
    %   Detailed explanation goes here
    properties (SetAccess = protected, GetAccess = public)
        compounding = -1
    end
    
    %Public methods
    methods 
        % Constructor
        function obj = YieldTermStructure(dayCount, refDate)
            
            % call super class constructor
            obj = obj@TermStructure(dayCount, refDate);
            
        end % end constructor

        % discount factor by dates
        function c = discountFactor(obj, d)
                c = discountFactorByTime(obj, timeFromReference(obj,d));
        end
        
        function c = forwardDiscountFactor(obj, startDate, endDate)
                c = discountFactor(obj, endDate)./discountFactor(obj, startDate);
        end
        
        function c = discountFactorByTime(obj, T)
            c = discF(obj, T);
        end
    end
     
    % Abstract methods
    methods (Abstract, Access = protected)
        discF(obj, T)
    end
    
end

