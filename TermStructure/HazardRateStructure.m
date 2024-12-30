classdef HazardRateStructure < PDTermStructure
    %HAZARDRATESTRUCTURE Summary of this class goes here
    %   Detailed explanation goes here
    properties
        hazardRates
    end
    
    methods
        function obj = HazardRateStructure(dayCounter, name, refDate, dates, hazardRates, varargin)
            if nargin < 5 
                error('HazardRateTermStructure:Too few parameters');
            end
            
            obj = obj@PDTermStructure(dayCounter, refDate);
            obj.name = name;
            obj.dates = dates;
            obj.times = timeFromReference(obj, dates);
            obj.hazardRates = hazardRates;

            if ~isempty(varargin) && ischar(varargin{1})
                p = inputParser;
                
                p.addParamValue('defaultProbs',[],@(x) (isnumeric(x) &&...
                all(x>=0) && all(x<=1) ));
                p.addParamValue('recoveryRate',0.4,@(x) (isnumeric(x) &&...
                all(x>=0) && all(x<=1) ));

                try
                    p.parse(varargin{:});
                catch ME
                    newMsg = message('HazardRateTermStructure:Variable input error');
                    newME = MException(newMsg.Identifier,getString(newMsg));
                    newME = addCause(newME,ME);
                    throw(newME)
                end
                
                obj.recoveryRate = p.Results.recoveryRate;
                obj.defaultProbs = p.Results.defaultProbs;
            else
                obj.recoveryRate = 0.4;
                obj.defaultProbs = [];
            end
            
            if ~isempty(hazardRates)
                obj.survivalProbs = Survival(obj, obj.times);
                obj.defaultProbs = 1 - obj.survivalProbs;
            else
                [obj.survivalProbs,obj.hazardRates] = internal.fininst.cdsstdsurvprob(...
                    obj.times,obj.defaultProbs,...
                    obj.times);
            end

        end
        
%         function c = getPDCurve(obj, 
    end
    
    methods(Access = protected)
        function  surv = Survival(obj, T)

            if size(T,1) < size(T,2)
                T = T';
            end
            expired = T <= 0;
            unexpired = T > 0;
            unexpiredT = [zeros(size(T(expired)));T(unexpired)];
            times = obj.times;
            h = obj.hazardRates;

            surv = ones(length(T),1);
            if isempty(unexpiredT)
                return
            end

            time0 = [0;times];
            dtime = diff(time0);

            %piecewise constant interpolation/extrapolation
            for i = 1:length(times)
                if i < length(times)
                    itmp = unexpiredT > time0(i) & unexpiredT <=time0(i+1);
                else
                    itmp = unexpiredT > time0(i);
                end
                H = 0;
                if i >1
                    H = dtime(1:i-1)'*h(1:i-1);
                end
                H = H + (unexpiredT(itmp)-time0(i))*h(i);
                surv(itmp) = exp(-H);
            end
        end
        
        function  pd = PD(obj, T)
            pd = 1-Survival(obj,T);
        end
        
    end
end

