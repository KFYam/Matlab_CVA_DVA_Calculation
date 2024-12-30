classdef FXCall < CalibrationInstrument
    %CAPLET Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = protected, GetAccess = public)        
        spot
        currency        
        volatility
        strike  
        inCurve_d
        inCurve_f
    end
        
    methods

        function obj = FXCall(spot,maturityYF,currency,volatility,strike,inCurve_d,inCurve_f,varargin)
                                
            if nargin < 7
                error(message('FXCall:too few input parameters'));
            end            
            
            obj.spot= spot;
            obj.maturityYF = maturityYF;
            obj.currency = currency;        
            obj.volatility = volatility;
            obj.strike = strike ;       
            obj.inCurve_d = inCurve_d;
            obj.inCurve_f = inCurve_f;
            obj.marketPrice = MarketPrice(obj);
            
        end
        function price = ModelPrice(obj, HW2Model_for,x)
            HW2Model_for.volTS = x;            
            price = yyPriceFXCall([obj.spot]',[obj.inCurve_d]',[obj.inCurve_f]',[obj.maturityYF]',HW2Model_for);
        end

        function c = MarketPrice(obj)            
            %r_d = intenvget(intenvset(obj.inCurve_d,'EndDates',datemnth(obj.inCurve_d.ValuationDate, 12*obj.maturityYF)),'Rates');  
            %r_f = intenvget(intenvset(obj.inCurve_f,'EndDates',datemnth(obj.inCurve_f.ValuationDate, 12*obj.maturityYF)),'Rates');              
            eps = 0.00001;
            r_d = max(intenvget(intenvset(getRateSpec(obj.inCurve_d),'EndDates',obj.inCurve_d.referenceDate + 365*obj.maturityYF),'Rates'),eps);  
            r_f = max(intenvget(intenvset(getRateSpec(obj.inCurve_f),'EndDates',obj.inCurve_f.referenceDate + 365*obj.maturityYF),'Rates'),eps);  
            [c,~] = blsprice(obj.spot,obj.spot,r_d,obj.maturityYF,obj.volatility,r_f);
        end
        function T = GetMaxMat(obj)
            T = max([obj.maturityYF]);
        end
    end       
end

