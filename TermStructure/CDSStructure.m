classdef CDSStructure < PDTermStructure
    %CDSSTRUCTURE Summary of this class goes here
    %   Detailed explanation goes here
    properties
        CDSRates
        interpolation
    end
    
    methods
        function obj = CDSStructure(dayCounter, name, refDate, dates, CDSRates, varargin)
            if nargin < 5 
                error('CDSStructure:Too few parameters');
            end

            obj = obj@PDTermStructure(dayCounter, refDate);
            obj.name = name;
            obj.dates = dates;
            obj.times = timeFromReference(obj, dates);
            obj.CDSRates = CDSRates/10000;

            if ~isempty(varargin) && ischar(varargin{1})
                p = inputParser;
                
                p.addParamValue('interpolation',{'linear'},...
                @(x) all(ismember(x,{'linear','cubic','spline'})))
                p.addParamValue('recoveryRate',0.4,@(x) (isnumeric(x) &&...
                all(x>=0) && all(x<=1) ));

                try
                    p.parse(varargin{:});
                catch ME
                    newMsg = message('CDSStructure:variable input error');
                    newME = MException(newMsg.Identifier,getString(newMsg));
                    newME = addCause(newME,ME);
                    throw(newME)
                end
            
                obj.interpolation = p.Results.interpolation;
                obj.recoveryRate = p.Results.recoveryRate;
            else
                obj.interpolation = 'linear';
                obj.recoveryRate = 0.4;
            end
            
            obj.survivalProbs = Survival(obj, obj.times);
            obj.defaultProbs = 1 - obj.survivalProbs;

        end
    end
    
    methods(Access = protected)
        function  surv = Survival(obj, T)

            if size(T,1) < size(T,2)
                T = T';
            end
            expired = T <= 0;
            unexpired = T > 0;
            unexpiredT = T(unexpired);
            times = obj.times;
            cds = obj.CDSRates;
            r = obj.recoveryRate;

            if isempty(unexpiredT)
                return
            end
            
            %flat extrapolation
            temp = unexpiredT < times(1);
            unexpiredT(temp) = times(1);
            temp = unexpiredT > times(end);
            unexpiredT(temp) = times(end);

            %interpolation, all in vertical array
            temp = interp1q(times, cds, unexpiredT);
            surv = ones(length(T),1);
            surv(unexpired) = exp(- temp .* T(unexpired)./(1-r)); 
            
            if ~isempty(expired)
                surv(expired) = 1;
            end
        end
        
        function  pd = PD(obj, T)
            pd = 1-Survival(obj,T);
        end
        
    end
end



