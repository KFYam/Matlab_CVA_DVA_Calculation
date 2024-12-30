classdef InterpolatedZeroCurve < YieldTermStructure
    %INTERPOLATEDZEROCURVE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = protected, GetAccess = public)
        dates %dates array
        rates %zero rate array
        times %times array
        interpolation %interpolation method
    end
    
    methods
        function obj = InterpolatedZeroCurve(dc, refDate, dates, rates, varargin)
            
            if nargin < 4 
                error(message('InterpolatedYieldCurve:Too few parameters'));
            end
            
            obj = obj@YieldTermStructure(dc, refDate);
                            
            if ~isempty(varargin) && ischar(varargin{1})
                p = inputParser;
                
                p.addParamValue('interpolation',{'linear'},...
                @(x) all(ismember(x,{'linear','cubic','spline'})))
                p.addParamValue('compounding',-1,@(x) ismember(x,[-1 0 1 2 3 4 6 12]));
             
                try
                    p.parse(varargin{:});
                catch ME
                    newMsg = message('InterpolatedYieldTermStructure:variable input error');
                    newME = MException(newMsg.Identifier,getString(newMsg));
                    newME = addCause(newME,ME);
                    throw(newME)
                end
               
                obj.interpolation = p.Results.interpolation;
                obj.compounding = p.Results.compounding;
                
            else
                % set default value
                obj.interpolation = 'linear';
                obj.compounding = -1;
            end % end if
            
            if ~isvector(dates)||~isvector(rates)
                 error(message('InterpolatedYieldCurve:need vector inputs'));
            end
            
            if size(dates) ~= size(rates)
                error(message('InterpolatedYieldCurve:vector size need to be the same'));
            end
            
            obj.dates = dates;
            obj.rates = rates;
            obj.times = timeFromReference(obj, dates);
            
        end % end constructor
        
        function c = getRateSpec(obj)
            c = intenvset('Rates', obj.rates', 'StartDates',obj.referenceDate,'EndDates', ...
                obj.dates', 'Compounding', obj.compounding, 'Basis', obj.dayCountConvention);
        end
    end
    
    methods (Access = protected)
        %Interpolate zero rates
        function  df = discF(obj, T)
            
            expired = T <= 0;
            unexpired = T > 0;
            unexpiredT = T(unexpired);
            
            if isempty(unexpiredT)
                df = ones(1, length(T));
                return
            end
            
            %flat extrapolation
            temp = unexpiredT < obj.times(1);
            unexpiredT(temp) = obj.times(1);
            temp = unexpiredT > obj.times(end);
            unexpiredT(temp) = obj.times(end);
            
            %interpolation, all in vertical array
            temp = interp1q((obj.times)', (obj.rates)', unexpiredT')';
            %df = transpose(rate2disc(obj.compounding, transpose(df), transpose(T)));
            %df = zero2disc(df, d, obj.referenceDate, obj.compounding, obj.dayCountConvention);
            df = ones(1, length(T));
            df(unexpired) = exp(- temp .* T(unexpired)); % For CVA, all assming continuous compounding
            
            if ~isempty(expired)
                df(expired) = 1;
            end
        end
        
    end
    
end

