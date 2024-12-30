classdef (Abstract) TermStructure < handle
    %TERMSTRUCTURE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = protected, GetAccess = public)
        dayCountConvention  % dayCountConvention
        referenceDate   % reference date
    end
    
    methods
        % Constructor
        function obj = TermStructure(dc, refDate)
            obj.dayCountConvention = dc;
            obj.referenceDate = refDate;
        end
        
        % Property set methods
        function set.dayCountConvention(obj, dc)
            if dc < 0 || dc > 13
                throw(MException('ERROR:InputParameters',...
                    'Matlab day count convention code need to be between 0 - 13'));
            else
                obj.dayCountConvention = dc;
            end
        end
        
        function set.referenceDate(obj, refDate)
            obj.referenceDate = refDate;
        end
        
        % Property get methods
        function c = get.dayCountConvention(obj)
            c = obj.dayCountConvention;
        end
        function c = get.referenceDate(obj)
            c = obj.referenceDate;
        end
        
        % member functions
        function c = timeFromReference(obj, d)
            %c = yearfrac(obj.referenceDate, d, obj.dayCountConvention);
            c = (d - obj.referenceDate)/365;
        end
    end
    
end

