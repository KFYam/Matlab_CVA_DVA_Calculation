classdef Caplet < CalibrationInstrument
    %CAPLET Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = protected, GetAccess = public)
        fixingYF
        currency        
        volatility
        strike  
        inCurve
        notional
    end
        
    methods

        function obj = Caplet(fixingYF,maturityYF,currency,volatility,strike,inCurve,notional, varargin)
                                
            if nargin < 7
                error(message('Caplet:too few input parameters'));
            end
            
            if maturityYF < fixingYF
                error(message('Caplet:maturity date cant be larger than fixing date'));
            end
            
            obj.fixingYF = fixingYF;
            obj.maturityYF = maturityYF;
            obj.currency = currency;        
            obj.volatility = volatility;
            obj.strike = strike ;   
            obj.notional = notional;
            obj.inCurve = inCurve;
            obj.marketPrice = MarketPrice(obj);
            
        end
        function c = ModelPrice(obj,HW2F_obj,x)
            HW2F_obj.a = x(1);
            HW2F_obj.b = x(2);
            HW2F_obj.sigma = x(3);
            HW2F_obj.eta = x(4);
            c = PriceCaplet2F([obj.inCurve],HW2F_obj,[obj.strike],[obj.fixingYF],[obj.maturityYF]);            
        end
        function c = ModelPrice2(obj,a,b,sigma,eta,rho)
            c = PriceCaplet([obj.inCurve],[obj.strike],[obj.fixingYF],[obj.maturityYF],a,b,sigma,eta,rho);
        end

        function c = MarketPrice(obj)
            PM = @(t)discountFactorByTime(obj.inCurve, t);            
            FwdRate = PM(obj.maturityYF)/PM(obj.fixingYF);
            dt = obj.maturityYF - obj.fixingYF;
            d1 = (log(FwdRate / obj.strike) + obj.volatility^2 * (obj.fixingYF)/2) ./ ...
                (obj.volatility * sqrt(obj.fixingYF));
            d2 = d1 - obj.volatility * sqrt(obj.fixingYF);
            c = dt * PM(obj.fixingYF)*(FwdRate * normcdf(d1) - obj.strike * normcdf(d2));
        end
        function T = GetMaxMat(obj)
            T = max([obj.maturityYF]);
        end
    end       
end

