classdef BlackATMVolCurve < VolatilityTermStructure
    %BLACKATMFXVOLCURVE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dates
        volatilities
        variances
        times
        interpolation % linearOnVariance, linearOnVolatility, flatOnVolatility
    end
    
    methods
        function obj = BlackATMVolCurve(refDate, dayCountConvention, ...
                dates, vols, varargin)
            
            if nargin < 4
                error(message('BlackATMVolCurve:Too few parameters'));
            end
            
            %obj = obj@YieldTermStructure(dc, refDate);
            businessDayConvention = 'modifiedfollow';
            finCalendar = FinCalendar('Weekends',[]);
            interp = 'flatOnVolatility';
            
            if ~isempty(varargin) && ischar(varargin{1})
                p = inputParser;
                
                p.addParamValue('interpolation',{'flatOnVolatility'},...
                @(x) all(ismember(x,{'flatOnVolatility','linearOnVolatility','linearOnVariance'})))
                p.addParamValue('finCalendar',FinCalendar('Weekends',[]));
                p.addParamValue('businessDayConvention','modifiedfollow');
             
                try
                    p.parse(varargin{:});
                catch ME
                    newMsg = message('BlackATMVolCurve:variable input error');
                    newME = MException(newMsg.Identifier,getString(newMsg));
                    newME = addCause(newME,ME);
                    throw(newME)
                end
               
                interp = p.Results.interpolation;
                finCalendar = p.Results.finCalendar;
                businessDayConvention =  p.Results.businessDayConvention;
            end % end if
            
            obj = obj@VolatilityTermStructure(refDate, dayCountConvention,businessDayConvention, finCalendar);
            obj.interpolation =  interp;
            
            if ~isvector(dates)||~isvector(vols)
                 error(message('InterpolatedYieldCurve:need vector inputs'));
            end
            
            if size(dates) ~= size(vols)
                error(message('InterpolatedYieldCurve:vector size need to be the same'));
            end
            
            obj.dates = dates;
            obj.volatilities = vols;
            obj.times = timeFromReference(obj, dates);
            obj.variances = obj.volatilities .* obj.volatilities .* obj.times;
        end % end of constructor
        
        function c = blackVol(obj, d)
            switch obj.interpolation
                case 'flatOnVolatility'
                    %flat extrapolation on both sides
                    c = d;
                    c(:) = obj.volatilities(1);
                    for  i = 1 : length(obj.dates)
                        temp = d > obj.dates(i);
                        c(temp) = obj.volatilities(i);
                    end
                case 'linearOnVolatility'
                    %flat extrapolation
                    temp = d < obj.dates(1);
                    d(temp) = obj.dates(1);
            
                    temp = d > obj.dates(length(obj.dates));
                    d(temp) = obj.dates(length(obj.dates));
                    
                    c = interp1(obj.times, obj.volatilities, timeFromReference(obj, d));
                case 'linearOnVariance'
                    c = (blackVariance(obj, d) / timeFromReference(obj, d)).^0.5;
            end
        end
        
        function c = forwardBlackVol(obj, startDate, endDate)
            fwdVar = forwardBlackVariance(obj, startDate, endDate);
            c = (fwdVar ./ (timeFromReference(obj, endDate) - timeFromReference(obj, startDate))) .^0.5;
        end
        
        function c = blackVariance(obj, d)
            switch obj.interpolation
                case {'flatOnVolatility', 'linearOnVolatility'}
                    vol = blackVol(obj, d);
                    c = vol .*  vol .* timeFromReference(obj, d);
                case 'linearOnVariance'
                    % linear extrapolation on variance
                    c = interp1(obj.times,obj.variances,timeFromReference(obj, d),'linear','extrap');
            end
        end
        
        function c = forwardBlackVariance(obj, startDate, endDate)
            c = blackVariance(obj, endDate) - blackVariance(obj, startDate);
        end
        
    end
    
end

