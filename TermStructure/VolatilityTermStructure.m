classdef VolatilityTermStructure < TermStructure
    %VOLATILITYTERMSTRUCTURE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        businessDayConvention
        finCalendar
    end
    
    methods
        function obj = VolatilityTermStructure(refDate, dayCountConvention, businessDayConvention, finCalendar)
            obj = obj@TermStructure(dayCountConvention, refDate);
            
            if ~ismember(businessDayConvention, ['follow','actual','modifiedfollow','previous','modifiedprevious'])
                error(message('VolatilityTS: wrong business convention type'));
            else
                obj.businessDayConvention = businessDayConvention;
            end
            
            if ~isa(finCalendar, 'FinCalendar')
                error(message('VolatilityTS: Need a calendar handle obj'));
            else
                obj.finCalendar = finCalendar;
            end
        end
    end
    
end

