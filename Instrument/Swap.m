classdef Swap < Instrument
    %SWAP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        leg1Handle
        leg2Handle
        payrec % 2x1 vector [-1; 1] means leg1 is payleg, leg 2 is recieve leg
    end
    
    methods
        function obj = Swap(leg1, leg2, payrec)
            obj = obj@Instrument();
            obj.leg1Handle = leg1;
            obj.leg2Handle = leg2;
            obj.payrec = payrec;
        end
            
        
        function c = calculate(obj, valueDate, env)
            v1 = calculate(obj.leg1Handle, valueDate, env);
            v2 = calculate(obj.leg2Handle, valueDate, env);
            
            obj.npv = [getFXSpot(env, obj.leg1Handle.currency), getFXSpot(env, obj.leg2Handle.currency)] .* [v1,v2] * obj.payrec;
            c = obj.npv;
        end
        
        function [mtmTimeGrid, results] = pathPricer(obj, valueDate, env, pricingGrid)
            %calculate MTM profile at simulation grid
            [mtmTimeGrid, mtmProfile1] = pathPricer(obj.leg1Handle, valueDate, env, []);
            [~, mtmProfile2] = pathPricer(obj.leg2Handle, valueDate, env, []);
            
            %aggregate mtm profile of two legs, expressed all in base ccy
            mtmProfile = mtmProfile1  * obj.payrec(1) + mtmProfile2 * obj.payrec(2);
            
            %aggregate npv of two legs, expressed in leg ccy
            obj.npv = obj.leg1Handle.npv * obj.payrec(1) * getFXSpot(env, obj.leg1Handle.currency) ...
                + obj.leg2Handle.npv * obj.payrec(2) * getFXSpot(env, obj.leg2Handle.currency); 
            
            %apply CSA and project MTM profile onto pricing grid
            [mtmTimeGrid, results] = applyCSA(obj, mtmTimeGrid, mtmProfile, pricingGrid);
        end        
        
    end
    
end

