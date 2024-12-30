classdef SwapFloatingLeg < Instrument
    %SWAPFLOATINGLEG Summary of this class goes here
    %   Detailed explanation goes here
    
    
    properties (SetAccess = protected, GetAccess = public)
        startDate
        maturityDate
        dayCountConvention
        notional
        indexTenor
        spread
        paymentFrequency
        currency
        businessDayConvention
        calendar
        endOfMonthRule
        lastFixing
        notionalAtStart
        notionalAtMaturity
        swapType
        payAtMaturity  %if coupon payment delayed to maturity date
        exAmtNotional %whether exchange armotized notional
        
        % leg schedule
        lastFixingDate
        lastFixingIndex
        CFAmounts
        CFDates
        CFTimes
        
    end
        
    methods
        
        function obj = SwapFloatingLeg(swapType, startDate, maturityDate, dayCountConvention, notional, ...
                                       indexTenor, spread, paymentFrequency, currency, varargin)
            obj = obj@Instrument();
            
            if nargin < 9
                error(message('SwapFloatLeg:too few input parameters'));
            end
            
            if maturityDate < startDate
                error(message('SwapFixedLeg:maturity date cant be larger than start date'));
            end
            
            if ~ismember(paymentFrequency, [0, 1, 2, 3, 4, 6, 12])
                error(message('SwapFixedLeg:wrong payment frequency'));
            end
            
            if ~ismember(swapType, ['IRS', 'CCS'])
                error(message('SwapFixedLeg:wrong swap type input, only support IRS or CCS'));
            end
            
            obj.swapType = swapType;
            obj.startDate = startDate;
            obj.maturityDate = maturityDate;
            obj.dayCountConvention = dayCountConvention;
            obj.indexTenor = indexTenor;
            obj.spread = spread;
            obj.paymentFrequency = paymentFrequency;
            obj.currency = currency;            
            
            if ~isempty(varargin) && ischar(varargin{1})
                p = inputParser;
   
                p.addParamValue('Calendar',[]);
                p.addParamValue('BusinessDayConvention',{'actual'},...
                    @(x) all(ismember(x,{'actual','follow','previous','modifiedfollow',...
                    'modifiedprevious'})))
                p.addParamValue('EndOfMonthRule', false);
                p.addParamValue('LastFixing',[]);
                p.addParamValue('notionalAtStart', false);
                p.addParamValue('notionalAtMaturity', false);
                p.addParamValue('payAtMaturity', false);
                p.addParamValue('exAmtNotional', false);
                
                try
                    p.parse(varargin{:});
                catch ME
                    newMsg = message('SwapFloatingleg:Can not parse');
                    newME = MException(newMsg.Identifier,getString(newMsg));
                    newME = addCause(newME,ME);
                    throw(newME)
                end      
                obj.calendar = p.Results.Calendar;
                obj.businessDayConvention = p.Results.BusinessDayConvention;
                obj.endOfMonthRule = p.Results.EndOfMonthRule;  
                if isnan(p.Results.LastFixing)
                    obj.lastFixing = [];
                else
                    obj.lastFixing = p.Results.LastFixing;
                end
                obj.notionalAtStart = p.Results.notionalAtStart;
                obj.notionalAtMaturity = p.Results.notionalAtMaturity;
                obj.payAtMaturity = p.Results.payAtMaturity;
                obj.exAmtNotional = p.Results.exAmtNotional;
            else
                % all default values
                obj.calendar = [];
                obj.businessDayConvention = 'actual';
                obj.endOfMonthRule = false;
                obj.lastFixing = [];
                obj.notionalAtStart = false;
                obj.notionalAtMaturity = false;
                obj.payAtMaturity = false;
                obj.exAmtNotional = false;
            end %end if
                                
            % generate leg schedule using fixed rate 1
            if length(notional) > 1
                %for armotizing swap            
                if ~isvector(startDate) || ~isvector(startDate) || size(startDate, 2) ~=1 || size(maturityDate, 2) ~=1
                    error(message('FloatingLeg:Need vertical vector of dates for armotizing swap'));
                end
                
                if size(notional, 2)~=1 || size(notional, 1) ~= size(startDate, 1) || size(maturityDate, 1) ~= size(startDate, 1)
                    error(message('FloatingLeg:start date, end date and notional profile size need to be the same'));
                end
                
                obj.CFDates = startDate;
                %add last cash flow date
                obj.CFDates(length(startDate)+1) = maturityDate(length(maturityDate));
                %calculate year fraction
                obj.CFTimes = yearfrac(startDate(1),obj.CFDates, dayCountConvention);
                
                %ser the notional profile
                obj.notional(1) = notional(1);
                obj.notional(2:length(notional)+1,1) = notional;
            else
                %generate swap schedule using bbg buildin function
                 [~, obj.CFDates]= cfamounts(1, startDate, maturityDate, ...
                 'Period',paymentFrequency, 'Basis', dayCountConvention, 'EndMonthRule', obj.endOfMonthRule, ...
                 'Face', notional, 'StartDate', startDate, 'BusinessDayConvention', obj.businessDayConvention, 'holiday', holidays);
                obj.CFDates = obj.CFDates';
                obj.CFTimes = yearfrac(obj.CFDates(1),obj.CFDates, dayCountConvention);
                
                obj.notional = notional*ones(length(obj.CFTimes),1);
            end
            
            %varing spread swap
            if size(spread, 1)~=1
                if size(spread, 1) ~= length(obj.CFTimes)-1
                     error(message('FloatingLeg:variable spread must have same size as the coupon cash flow'));
                end
            else
                obj.spread = repmat(obj.spread, length(obj.CFTimes)-1 , 1);
            end
            
            %set last fixing date if there is fixing input
            env = Environment.getInstance;
            today = env.valuationDate;
            
            if obj.CFDates(1) <= today && today < obj.CFDates(end)              
                %estimate the fixing date
                for i = 2 : length(obj.CFDates)
                    if obj.CFDates(i) >= today
                        obj.lastFixingDate = obj.CFDates(i-1);
                        obj.lastFixingIndex = i;
                        break;
                    end
                end
            end % end of fixing setting
            
        end % end constructor
        
        
        %calculate MTM for floating leg
        function c = calculate(obj, valueDate, env)
            
             if ~isa(env, 'Environment')
                error(message('SwapFloatingLeg:Need environment obj'));
            end
                
            modelAssignment = GlobalSetting.getModelAssignment('IR_IRS_MTM');
            
            switch modelAssignment
                case 'ANALYTICAL'
                    c = swapAnalyticalEngine(obj, valueDate, env);                    
                otherwise
                    error(message('Fixedleg:pricing model not supported'));
            end %End switch
                               
            obj.npv = c;
        end
        
        %--------------------------------------------%
        %----    analytical IRS pricing engine   ----%
        %--------------------------------------------%% analytical IRS pricing engine
        function c = swapAnalyticalEngine(obj, valueDate, env)
                    
            refCurve = GlobalSetting.getCurveAssignment(obj.swapType, obj.currency, 'FORECAST', obj.indexTenor);
            disCurve = GlobalSetting.getCurveAssignment(obj.swapType, obj.currency, 'DISCOUNT', '');
            refHandle = getIRCurveHandle(env, refCurve);
            discountHandle = getIRCurveHandle(env, disCurve);
                    
            % Set the cashflow date to value date if expired
            expired = obj.CFDates <= valueDate;
            unexpired = obj.CFDates > valueDate;
            unexpiredCFDates = obj.CFDates;
            unexpiredCFDates(expired) = refHandle.referenceDate;
            
            % calculate the floating coupon payment
            % this is an approximation for CVA purpose only that use the cash flow payment
            % dates instead of the actual accrual period          
            obj.CFAmounts(2:length(unexpiredCFDates), 1) = forwardDiscountFactor(refHandle, unexpiredCFDates(1:end-1)', ...
                unexpiredCFDates(2:end)');
            obj.CFAmounts(2:end) = (1./obj.CFAmounts(2:end)) - 1;
                       
            %adjust for last fixing
            if length(obj.lastFixing) ~= 0
                %with given fixing, override
                obj.CFAmounts(obj.lastFixingIndex) = obj.lastFixing * (obj.CFTimes(obj.lastFixingIndex) -  ...
                 obj.CFTimes(obj.lastFixingIndex - 1));
%             else
%                 %without given fixing, approximated using yield curve
%                 obj.CFAmounts(obj.lastFixingIndex) = obj.CFAmounts(obj.lastFixingIndex) * ...
%                     (obj.CFDates(obj.lastFixingIndex) - obj.CFDates(obj.lastFixingIndex - 1)) / ...
%                     (obj.CFDates(obj.lastFixingIndex) - refHandle.referenceDate);
            end
            
            % adjust the notioanl payment at the start or at the end
            if obj.notionalAtStart == true && obj.CFDates(1) >= valueDate
                obj.CFAmounts(1) = obj.CFAmounts(1) - 1;
            end
             
            if obj.notionalAtMaturity == true && obj.CFDates(end) >= valueDate
                obj.CFAmounts(end) = obj.CFAmounts(end) + 1;
            end
            
            %Add the spread
            yf = diff(obj.CFTimes);
            if obj.CFDates(1) <= valueDate
                yf = diff(obj.CFTimes);
                obj.CFAmounts(unexpired) = obj.CFAmounts(unexpired) + ...
                    obj.spread(end - length(obj.CFAmounts(unexpired))+1 : end ).* ...
                     yf(end - length(obj.CFAmounts(unexpired))+1 : end );
            else %forward starting
                obj.CFAmounts(unexpired) = obj.CFAmounts(unexpired) + [0;obj.spread].*[0; yf];
            end
            
            %times the notional profile
            obj.CFAmounts = obj.CFAmounts .* obj.notional;
              
            %handle periodical exchange of armotized notional
            if obj.exAmtNotional == true
                notionalAmt = diff(obj.notional(unexpired));
                obj.CFAmounts(unexpired) = obj.CFAmounts(unexpired) - [notionalAmt; 0];
%                 notionalAmt = diff(obj.notional(2:end));
%                 if size(obj.CFAmounts, 1) == length(obj.CFDates)
%                     %all cash flow not expired, apply the notional
%                     %exchange from second cashflow
%                     obj.CFAmounts(2:end - 1) = obj.CFAmounts(2:end - 1, :) - notionalAmt;
%                 else
%                     obj.CFAmounts(1:end - 1, :) = obj.CFAmounts(1:end - 1, :) - ...
%                         notionalAmt(end - length(obj.CFAmounts, 1) + 1:end);
%                 end
            end
            
            % discount all the cash flows
            if obj.payAtMaturity == false
                dfs = forwardDiscountFactor(discountHandle, valueDate, unexpiredCFDates');  
            else
                df = forwardDiscountFactor(discountHandle, valueDate, unexpiredCFDates(end));
                dfs = unexpiredCFDates';
                dfs(:) = df;
            end
            
            obj.npv = obj.CFAmounts' * dfs';      
            c = obj.npv;        
        end
        
        %--------------------------------------------%
        %----      calculate cva profile         ----%
        %--------------------------------------------%
        function [mtmTimeGrid, results] = pathPricer(obj, valueDate, env, pricingGrid)
%             modelAssignment = GlobalSetting.getModelAssignment('IR_IRS_CVA');
%             
%             switch modelAssignment
%                 case 'G2'
%                     modelHandle = getIRModelHandle(env, strcat(obj.currency, '_', modelAssignment));
%                 otherwise
%                     error(message('SwapFixedLeg:model not supported'));
%             end
            if isKey(env.IRModelCollection,  strcat(obj.currency, '_G2'))
                modelHandle = getIRModelHandle(env, strcat(obj.currency, '_G2'));
            else 
                modelHandle = getIRModelHandle(env, strcat(obj.currency, '_CRM'));
            end
             
            timeGrid = modelHandle.timeGrid;
            
            if isempty(timeGrid)
                error(message('SwapFloatingLeg:Model need a simulated path'));
            end
            
            nTimeGrid = length(timeGrid);
            nPath = size(modelHandle.xt, 2);
            results = zeros(nTimeGrid, nPath);
            
            %for iPath = 1:nPath
            for iTimeGrid = 1  :nTimeGrid
                %if the dela already expired
                                
                if timeGrid(iTimeGrid) >= timeFromReference(modelHandle.irTS, obj.maturityDate)
                    results(iTimeGrid:end, :) = 0;
                    break
                end
                
                %if not
                xt = modelHandle.xt(iTimeGrid, :);
                yt = modelHandle.yt(iTimeGrid, :);   
                
                if length(modelHandle)~= 0
                    st = modelHandle.st(iTimeGrid, :);  
                else
                    st = ones(1, nPath);
                end
                
                unexpired = timeFromReference(modelHandle.irTS, obj.CFDates) > timeGrid(iTimeGrid);                           
                nUnexpired = length(obj.CFTimes(unexpired));
                
                %calculate DF
%                 dfs = P(modelHandle, timeGrid(iTimeGrid),  obj.CFTimes(unexpired)', xt, yt);    
                dfs = P(modelHandle, timeGrid(iTimeGrid),  timeFromReference(modelHandle.irTS, obj.CFDates(unexpired))', xt, yt); 
                cfAmounts = ones(size(dfs,1),size(dfs, 2));
                cfAmounts(2 : end,:) = dfs(1 : end - 1,:);
                cfAmounts = cfAmounts ./ dfs;            
                cfAmounts = (cfAmounts - 1); 
                
                %handl first fixing
                if unexpired(obj.lastFixingIndex) == 1
                    if length(obj.lastFixing) ~= 0
                        if unexpired(1) == 1
                            %if the first notional exchange hasnt expired
                            %(forward start IRS, then adjust second cash
                            %flow
                            cfAmounts(2, :) = obj.lastFixing * (obj.CFTimes(obj.lastFixingIndex) - obj.CFTimes(obj.lastFixingIndex - 1));
                        else
                            %otherwise adjust the first one
                            cfAmounts(1, :) = obj.lastFixing * (obj.CFTimes(obj.lastFixingIndex) - obj.CFTimes(obj.lastFixingIndex - 1));
                        end
                    end
                end
                                   
                %adjust start and end cash flow                
                if obj.notionalAtStart == true && unexpired(1) == 1
                    cfAmounts(1, :) = cfAmounts(1, :) - 1;
                end
                
                %adjust the last cashflow for armotizing swap
                if obj.notionalAtMaturity == true
                    cfAmounts(end, : ) =  cfAmounts(end, : ) + 1;
                end
                
                %Add the spread
                if unexpired(1) == 1     
                    yf = diff(obj.CFTimes);
                    try
                        cfAmounts = cfAmounts + repmat([0, obj.spread.*yf], 1, nPath);                    
                    catch me
                        c = 0;
                    end
                else
                    yf = diff([obj.CFTimes(end - nUnexpired);obj.CFTimes(unexpired)]);
                    try
                        cfAmounts = cfAmounts + repmat(obj.spread(end - size(cfAmounts, 1)+1: end).*yf, 1, nPath);
                    catch me
                        c=0;
                    end
                end
                
                %calculate the cf amount
                cfAmounts = cfAmounts.* repmat(obj.notional(unexpired), 1, nPath);       
                
                %handle periodical exchange of armotized notional
                if obj.exAmtNotional == true
                    notionalAmt = diff(obj.notional(unexpired));
                    cfAmounts = cfAmounts - repmat([notionalAmt; 0],1,nPath);
%                     notionalAmt = diff(obj.notional(2:end));
%                     if size(cfAmounts, 1) == length(obj.CFDates)
%                         %all cash flow not expired, apply the notional
%                         %exchange from second cashflow
%                         cfAmounts(2:end - 1, :) = cfAmounts(2:end - 1, :) - repmat(notionalAmt,1,nPath);
%                     else
%                         %first
%                         cfAmounts(1:end - 1, :) = cfAmounts(1:end - 1, :) - ...
%                             repmat(notionalAmt(end - size(cfAmounts,1) + 2:end),1,nPath);
%                     end
                end
                
                %handle payment at maturity case
                if obj.payAtMaturity == true    
                    df = P(modelHandle, timeGrid(iTimeGrid),  obj.CFTimes(end), xt, yt);
                    dfs = repmat(df,size( obj.CFTimes(unexpired),1),1);
                end
                
                %discount the cash flow
                if size(dfs,1) ~= 1
                    results(iTimeGrid, :) = (sum(dfs .* cfAmounts)')' .* st;     
                else
                    results(iTimeGrid, :) = ((dfs .* cfAmounts)')' .* st;
                end
            end
            
            %add exposure at time 0
            currentMTM  = swapAnalyticalEngine(obj, valueDate, env);
            obj.npv = currentMTM; 
            rowOne = zeros(1, nPath);
            rowOne(:) = currentMTM * getFXSpot(env, obj.currency);
            results = [rowOne; results];
            mtmTimeGrid = [0; timeGrid];
        end
        
        function c = yyPrice(obj, vYF, mDiscMatrix)
            %we assume that the floating rate for the first cashflow is
            %fixed at the pricing date and will be equal to (1/P(t,t+YF(1)) -1)/YF(1)
            [iNbRow,iNbIteration] = size(mDiscMatrix);            
            iNbCashFlow = length(obj.notional);
            
            c = (obj.notional(iNbCashFlow - iNbRow + 1) / vYF(1)) * ones(1,iNbIteration) .*(1 - mDiscMatrix(1,:)) ;
            if iNbRow > 1
                c = [c ; -((obj.notional(iNbCashFlow - iNbRow + 2:end) ./ vYF(2:end)) * ones(1,iNbIteration)) .*diff(mDiscMatrix)] ;
            end
            if obj.spread ~= 0
                c = c + obj.spread *mDiscMatrix;
            end
            c = sum(c,2)/ iNbIteration;
            
        end
        
        function mC = yyPriceForProfile(obj, vYF, mDiscMatrix)
            %we assume that the floating rate for the first cashflow is
            %fixed at the pricing date and will be equal to (1/P(t,t+YF(1)) -1)/YF(1)
            [iNbRow,iNbIteration] = size(mDiscMatrix);            
            iNbCashFlow = length(obj.notional);
            
            mC = (obj.notional(iNbCashFlow - iNbRow + 1) / vYF(1)) * ones(1,iNbIteration) .*(1 - mDiscMatrix(1,:)) ;
            if iNbRow > 1
                mC = [mC ; -((obj.notional(iNbCashFlow - iNbRow + 2:end) ./ vYF(2:end)) * ones(1,iNbIteration)) .*diff(mDiscMatrix)] ;
            end
            if obj.spread ~= 0
                mC = mC + obj.spread *mDiscMatrix;
            end                        
        end                   
    end
    
end

