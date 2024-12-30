classdef FXForward < Instrument
    
    properties(SetAccess = private)
        dMaturityDate_
        ccy1_
        ccy2_
        notional1_
        notional2_
    end
    
    properties(Dependent = true, SetAccess = private)
        strike_
    end
    
    methods
        function c = FXForward(ccy1, ccy2, notional1, notional2, dMaturityDate)
            
            c = c@Instrument();
            
            if nargin > 0
                if ~isa(ccy1, 'char') || ~isa(ccy2, 'char') ...
                || ~isa(notional1, 'numeric')|| ~isa(notional2, 'numeric')
                    err = MException('objCreation:missingArguments', 'Empty argument list');
                    throw(err);
                else 
                    c.ccy1_=ccy1;
                    c.ccy2_ =ccy2;
                    c.notional1_ = notional1;
                    c.notional2_ = notional2;
                    c.dMaturityDate_ = dMaturityDate;
                end
            else
                err = MException('objCreation:missingArguments', 'Empty argument list');
                throw(err);
            end
        end
        
        function c = calculate(obj, valueDate, env)
            
            if ~isa(env, 'Environment')
                error(message('FXSwap:Need environment obj'));
            end
            
            modelAssignment = GlobalSetting.getModelAssignment('FX_FXFORWARD_MTM');
            
            switch modelAssignment
                case 'ANALYTICAL'
                        c = fxforwardAnalyticalEngine(obj, valueDate, env);                  
                otherwise
                    error(message('FXForward:pricing model not supported'));
            end %End switch
            
            obj.npv = c;
        end
        
        function c = fxforwardAnalyticalEngine(obj, valueDate, env)
            
            ccy1CurveName = GlobalSetting.getCurveAssignment('FXFORWARD', obj.ccy1_, 'DISCOUNT', '');
            ccy2CurveName = GlobalSetting.getCurveAssignment('FXFORWARD', obj.ccy2_, 'DISCOUNT', '');
            ccy1DisHandle = getIRCurveHandle(env, ccy1CurveName);
            ccy2DisHandle = getIRCurveHandle(env, ccy2CurveName);
            ccy1FxSpot = getFXSpot(env, obj.ccy1_);
            ccy2FxSpot = getFXSpot(env, obj.ccy2_);
                    
            if obj.dMaturityDate_ <= valueDate
                c = 0;
            else
                c = obj.notional1_ * forwardDiscountFactor(ccy1DisHandle, valueDate, obj.dMaturityDate_) * ccy1FxSpot + ...
                    obj.notional2_ * forwardDiscountFactor(ccy2DisHandle, valueDate, obj.dMaturityDate_) * ccy2FxSpot;
            end
        end
        
        function [mtmTimeGrid, results] = pathPricer(obj, valueDate, env, pricingGrid)
%              modelAssignment = GlobalSetting.getModelAssignment('FX_FXFORWARD_CVA');
%             
%             switch modelAssignment
%                 case 'G2'
%                     model1Handle = getIRModelHandle(env, strcat(obj.ccy1_, '_', modelAssignment));
%                     model2Handle = getIRModelHandle(env, strcat(obj.ccy2_, '_', modelAssignment));
%                 otherwise
%                     error(message('FXForward:model not supported'));
%             end
%              
            if isKey(env.IRModelCollection,  strcat(obj.ccy1_, '_G2'))
                model1Handle = getIRModelHandle(env, strcat(obj.ccy1_, '_G2'));
            else 
                model1Handle = getIRModelHandle(env, strcat(obj.ccy1_, '_CRM'));
            end
            
            if isKey(env.IRModelCollection,  strcat(obj.ccy2_, '_G2'))
                model2Handle = getIRModelHandle(env, strcat(obj.ccy2_, '_G2'));
            else 
                model2Handle = getIRModelHandle(env, strcat(obj.ccy2_, '_CRM'));
            end
            
            timeGrid = model1Handle.timeGrid;
            
            if isempty(timeGrid)
                error(message('FXForward:Model need a simulated path'));
            end
            
            nTimeGrid = length(timeGrid);
            nPath = size(model1Handle.xt, 2);
            results = zeros(nTimeGrid, nPath);
            
            %for iPath = 1:nPath
            for iTimeGrid = 1:nTimeGrid
                %if the dela already expired
                if timeGrid(iTimeGrid) > timeFromReference(model1Handle.irTS, obj.dMaturityDate_)
                    results(iTimeGrid:end, :) = 0;
                    break
                end
                %if not price ccy1 leg
                xt = model1Handle.xt(iTimeGrid, :);
                yt = model1Handle.yt(iTimeGrid, :);

                if length(model1Handle.st) == 0
                    st1 = ones(1, nPath);
                else
                    st1 = model1Handle.st(iTimeGrid, :);
                end
                
                dfs1 = P(model1Handle, timeGrid(iTimeGrid),  ...
                    timeFromReference(model1Handle.irTS, obj.dMaturityDate_), xt, yt);
                
                %price ccy2 leg
                xt = model2Handle.xt(iTimeGrid, :);
                yt = model2Handle.yt(iTimeGrid, :);
        
                if length(model2Handle.st) == 0
                    st2 = ones(1, nPath);
                else
                    st2 = model2Handle.st(iTimeGrid, :);
                end
                
                dfs2 = P(model2Handle, timeGrid(iTimeGrid), ...
                    timeFromReference(model2Handle.irTS, obj.dMaturityDate_), xt, yt);
                    
                results(iTimeGrid, :) = obj.notional1_ * (dfs1.* st1) + obj.notional2_ * (dfs2.* st2);                
            end
            
            %add exposure at time 0
            currentMTM  = fxforwardAnalyticalEngine(obj, valueDate, env);
            obj.npv = currentMTM; 
            rowOne = zeros(1, nPath);
            rowOne(:) = currentMTM;
            results = [rowOne; results];
            mtmTimeGrid = [0; timeGrid];
            
            %apply csa and project exposure
            [mtmTimeGrid, results] = applyCSA(obj, mtmTimeGrid, results, pricingGrid);
        end
        
        function strike_ = get.strike_(obj)
            strike_ = abs(obj.notional1_ / obj.notional2_);
        end
    end
    
end

