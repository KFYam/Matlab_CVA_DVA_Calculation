function SwaptionPrice = PriceSwaption2F(vInCurve,HW2F_obj,vStrike,vYFExerciseDate,vYFMaturity)

    nSwaptions = length(vYFMaturity);
    SwaptionPrice = zeros(nSwaptions,1);
    
    for swapidx = 1:nSwaptions
        curve = vInCurve(swapidx);        
        
        T = vYFExerciseDate(swapidx);
        Tenor = vYFMaturity(swapidx);

        ti = [T:(Tenor + T)]';
        YF = diff(ti);
        
        P0 = discountFactorByTime(curve, ti')';

        LVL = YF'*P0(2:end);
        dSwpRate = (P0(1) - P0(end))/LVL;
        
        dV = sqrt(VarSwap(HW2F_obj,ti,P0,LVL, dSwpRate));
        dD = (dSwpRate - vStrike(swapidx))/dV;
        
        SwaptionPrice(swapidx) = (((dSwpRate - vStrike(swapidx))*normcdf(dD) + dV*normpdf(dD))*LVL);
    end
end


function V = VarSwap(HW2F_obj,vTSwp,P0,LVL, dSwpRate)
    dExpiry = vTSwp(1);
    vDiff = diff(vTSwp)'; %line
    
    vExp1 = exp(-HW2F_obj.a*vTSwp);
    dA1 = vExp1(1)*P0(1) - vExp1(end)*P0(end) - dSwpRate*sum(vDiff*(vExp1(2:end).*P0(2:end)));
    dA1 = dA1 /(HW2F_obj.a*LVL);
    
    vExp2 = exp(-HW2F_obj.b*vTSwp);
    dA2 = vExp2(1)*P0(1) - vExp2(end)*P0(end) - dSwpRate*sum(vDiff*(vExp2(2:end).*P0(2:end)));
    dA2 = dA2 /(HW2F_obj.b*LVL);
        
    V  = dA1*dA1* HW2F_obj.sigma * HW2F_obj.sigma *( exp(2*HW2F_obj.a*dExpiry) - 1) / (2*HW2F_obj.a);        
    V  = V + dA2*dA2* HW2F_obj.eta * HW2F_obj.eta *( exp(2*HW2F_obj.b*dExpiry) - 1) / (2*HW2F_obj.b);        
    V  = V + 2*HW2F_obj.rho*dA1*dA2*HW2F_obj.sigma*HW2F_obj.eta*(exp((HW2F_obj.a+HW2F_obj.b)*dExpiry) - 1)/(HW2F_obj.a+HW2F_obj.b);
end