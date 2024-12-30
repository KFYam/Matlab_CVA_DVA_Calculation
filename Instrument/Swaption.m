classdef Swaption < CalibrationInstrument
    %SWAPTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = protected, GetAccess = public)
        TenorYF
        currency        
        volatility
        strike
        inCurve
        type
        notional
        
    end
        
    methods

        function obj = Swaption(maturityYF,TenorYF,currency,volatility,inCurve,type,notional, varargin)
                                
            if nargin < 6
                error(message('Swaption:too few input parameters'));
            end
                    
            obj.maturityYF = maturityYF;
            obj.TenorYF = TenorYF;
            obj.currency = currency;        
            obj.volatility = volatility;            
            obj.notional = notional;
            obj.type = type;
            obj.inCurve = inCurve;
            if ~isempty(varargin) && ischar(varargin{1})
                p = inputParser;
                p.addParamValue('strike',[]);                
                try
                    p.parse(varargin{:});
                    catch ME
                        newMsg = message('Swaption:Can not parse');
                        newME = MException(newMsg.Identifier,getString(newMsg));
                        newME = addCause(newME,ME);
                        throw(newME)
                        throw(newME)
                    end      
                obj.strike = p.Results.strike;                    
            else
                obj.strike =GetATMStrike(obj);                    
            end
            obj.marketPrice = MarketPrice(obj);
            
        end
        function c = ModelPrice(obj,HW2F_obj,x)
            HW2F_obj.a = x(1);
            HW2F_obj.b = x(2);
            HW2F_obj.sigma = x(3);
            HW2F_obj.eta = x(4);
            c = PriceSwaption2F([obj.inCurve],HW2F_obj,[obj.strike],[obj.maturityYF],[obj.TenorYF]);       
%              curve = getRateSpec(obj(1).inCurve); 
%              c = swaptionbylg2f(curve,x(1),x(2),x(3),x(4),HW2F_obj.rho,[obj.strike],...
%                             obj(1).inCurve.referenceDate+365*[obj.maturityYF],obj(1).inCurve.referenceDate+365*[obj.maturityYF]+365*[obj.TenorYF])/100;
        end
        function c = ModelPrice2(obj,a,b,sigma,eta,rho,vYF)
            c = PriceSwaption([obj.inCurve],[obj.strike],[obj.maturityYF],[obj.TenorYF],a,b,sigma,eta,rho,vYF);      
            for i = 1 : length(c)
                obj(i).modelPrice = c(i);
            end
        end

        function c = MarketPrice(obj)
            Settle = obj.inCurve.referenceDate;            
            dExpiryDate = datemnth(Settle,12*obj.maturityYF);
            dMatDate = datemnth(dExpiryDate,12*obj.TenorYF);
            curve = getRateSpec(obj.inCurve); 
            c = swaptionbyblk(curve, {obj.type}, obj.strike,Settle,...
                                dExpiryDate, dMatDate, obj.volatility,'Principal',1);
        end
        function K = GetATMStrike(obj)
            Settle = obj.inCurve.referenceDate; 
            dExpiryDate = datemnth(Settle,12*obj.maturityYF);
            dMatDate = datemnth(dExpiryDate,12*obj.TenorYF);
            curve = getRateSpec(obj.inCurve); 
            [~,K] = swapbyzero(curve,[NaN 0], Settle, dMatDate,'StartDate',dExpiryDate);
        end
        function T = GetMaxMat(obj)
            T = max(max([obj.maturityYF]),max([obj.TenorYF]));
        end
    end       
end

        %calculate swaption price under G2 model
%         function SwaptionPrice = ModelPrice(obj, X, ExerciseDate, Maturity,varargin)
%      
%             if any(ExerciseDate > Maturity)
%                 error(message('G2Model:MaturityBeforeExercise'));
%             end
% 
%             p = inputParser;
% 
%             p.addParamValue('reset',1);
%             p.addParamValue('notional',100);
%             p.addParamValue('optspec',{'call'});
% 
%             try
%                 p.parse(varargin{:});
%             catch ME
%                 newMsg = message('G2Model:optionalInputError');
%                 newME = MException(newMsg.Identifier, getString(newMsg));
%                 newME = addCause(newME, ME);
%                 throw(newME)
%             end
% 
%             Reset = p.Results.reset;
%             Notional = p.Results.notional;
% 
%             if ischar(p.Results.optspec)
%                 OptSpec = cellstr(p.Results.optspec);
%             elseif iscell(p.Results.optspec)
%                 OptSpec = p.Results.optspec;
%             else
%                 error(message('G2Model:invalidOptSpec'));
%             end
% 
%             if ~all(ismember(lower(OptSpec),{'call','put'}))
%                 error(message('G2Model:invalidOptSpec'));
%             end
% 
%             try
%                 [X, ExerciseDate, Maturity, Reset, Notional, OptSpec] = finargsz(1, X(:), ExerciseDate(:),Maturity(:),...
%                     Reset(:), Notional(:), OptSpec(:));
%             catch ME
%                 throwAsCaller(ME)
%             end
% 
%             w = double(strcmpi(OptSpec,'call'));
%             w(w == 0) = -1;
% 
%             PM = @(t)discountFactorByTime(obj.irTS, t);
% 
%             V = @(t,T) obj.sigma^2/obj.a^2*(T - t + 2/obj.a*exp(-obj.a*(T-t)) - 1/(2*obj.a)*exp(-2*obj.a*(T-t)) - 3/2/obj.a) + ...
%                 obj.eta^2/obj.b^2*(T - t + 2/obj.b*exp(-obj.b*(T-t)) - 1/(2*obj.b)*exp(-2*obj.b*(T-t)) - 3/2/obj.b) + ...
%                 2*obj.rho*obj.sigma*obj.eta/(obj.a*obj.b)*(T - t + (exp(-obj.a*(T-t)) - 1)/obj.a + (exp(-obj.b*(T-t)) - 1)/obj.b - ...
%                 (exp(-(obj.a + obj.b)*(T-t)) - 1)/(obj.a + obj.b));
% 
%             A = @(t,T) PM(T)./PM(t) .*exp(1/2*(V(t,T) - V(0,T) + V(0,t)));
% 
%             B = @(z,t,T) (1 - exp(-bsxfun(@times,z,(T-t))))/z;
% 
%             nSwaptions = length(Maturity);
%             SwaptionPrice = zeros(nSwaptions,1);
% 
%             optOptions = optimset('Jacobian','on','display','off');
% 
%             for swapidx=1:nSwaptions
% 
%                 T = round(timeFromReference(obj.irTS, ExerciseDate(swapidx)));
%                 Tenor = round(-timeFromReference(obj.irTS, ExerciseDate(swapidx)) + timeFromReference(obj.irTS, Maturity(swapidx)));
% 
%                 ti = T:1/Reset(swapidx):(Tenor + T);
%                 tau = diff(ti);
%               
%                 c = X(swapidx)*tau;  
%                 c(end) = c(end) + 1;
%                 ti(1) = [];
% 
%                 ux = -(obj.sigma^2/obj.a^2 + obj.rho*obj.sigma*obj.eta/obj.a/obj.b)*(1 - exp(-obj.a*T)) + ...
%                     obj.sigma^2/(2*obj.a^2)*(1 - exp(-2*obj.a*T)) + ...
%                     obj.rho*obj.sigma*obj.eta/(obj.b*(obj.a+obj.b))*(1 - exp(-obj.b*T - obj.a*T));
% 
%                 uy = -(obj.eta^2/obj.b^2 + obj.rho*obj.sigma*obj.eta/obj.a/obj.b)*(1 - exp(-obj.b*T)) + ...
%                     obj.eta^2/(2*obj.b^2)*(1 - exp(-2*obj.b*T)) + ...
%                     obj.rho*obj.sigma*obj.eta/(obj.a*(obj.a+obj.b))*(1 - exp(-obj.b*T - obj.a*T));
% 
%                 sigx = obj.sigma*sqrt((1-exp(-2*obj.a*T))/2/obj.a);
%                 sigy = obj.eta*sqrt((1-exp(-2*obj.b*T))/2/obj.b);
%                 rhoxy = obj.rho*obj.sigma*obj.eta/((obj.a+obj.b)*sigx*sigy)*(1-exp(-(obj.a+obj.b)*T));
% 
%                 x = linspace(ux - 10*sigx,ux + 10*sigx,1001)';
% 
%                 cA = c.*A(T,ti);
%                 [ybar,~,exitflag] = fsolve(@(ybar) localObjFcn(ybar,x,cA,...
%                     B(obj.a,T,ti),B(obj.b,T,ti)),-x,optOptions);
% 
%                 if exitflag <= 0
%                     error(message('fininst:swaptionbylg2f:rootFailure'));
%                 end
% 
%                 h1 = (ybar - uy)./(sigy*sqrt(1 - rhoxy^2)) - rhoxy*(x - ux)./(sigx*sqrt(1 - rhoxy^2));
%                 h2 = bsxfun(@plus,B(obj.b,T,ti).*sigy*sqrt(1 - rhoxy^2),h1);
% 
%                 lambda = bsxfun(@times,A(T,ti).*c,exp(-bsxfun(@times,B(obj.a,T,ti),x)));
% 
%                 k = bsxfun(@times,-B(obj.b,T,ti),bsxfun(@plus,uy - .5*(1 - rhoxy.^2)*sigy^2*B(obj.b,T,ti),...
%                     rhoxy*sigy*(x-ux)/sigx));
% 
%                 Y = exp(-1/2*((x - ux)./sigx).^2)./(sigx*sqrt(2*pi)) .* ...
%                     (normcdf(-w(swapidx)*h1) - sum(lambda.*exp(k).*normcdf(-w(swapidx)*h2),2));
% 
%                 TempVal = trapz(x,Y);
% 
%                 SwaptionPrice(swapidx) = w(swapidx)*Notional(swapidx)*TempVal*PM(T);
%             end
%         end
        