function SwaptionPrice = PriceSwaption(vInCurve,vStrike,vYFExerciseDate,vYFMaturity,a,b,sigma,eta,rho,vYF)
   
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
        
        dV = sqrt(VarSwap2(a,b,sigma,eta,rho,vYF,ti,P0,LVL, dSwpRate));
        dD = (dSwpRate - vStrike(swapidx))/dV;
        
        SwaptionPrice(swapidx) = (((dSwpRate - vStrike(swapidx))*normcdf(dD) + dV*normpdf(dD))*LVL);
    end
end
%assuming piecewise constant volatility
function V = VarSwap2(a,b,sigma,eta,rho,vYF,vTSwp,P0,LVL, dSwpRate)
    iNumVol =1;    
    V1 = 0;
    V2 = 0;
    V12 = 0;
    iNbVolDates = length(vYF);    
    dExpiry = vTSwp(1);
    vDiff = diff(vTSwp)'; %line
    
    vExp1 = exp(-a*vTSwp);
    dA1 = vExp1(1)*P0(1) - vExp1(end)*P0(end) - dSwpRate*vDiff*(vExp1(2:end).*P0(2:end));
    dA1 = dA1 /(a*LVL);
    vExp2 = exp(-b*vTSwp);
    dA2 = vExp2(1)*P0(1) - vExp2(end)*P0(end) - dSwpRate*vDiff*(vExp2(2:end).*P0(2:end));
    dA2 = dA2 /(b*LVL);
    dTprev = 0;
    while vYF(iNumVol) <= dExpiry && iNumVol <= iNbVolDates
        V1  = V1 + sigma(iNumVol)*sigma(iNumVol)*(exp(2*a*vYF(iNumVol)) - exp(2*a*dTprev));        
        V2  = V2 + eta(iNumVol)*eta(iNumVol)*(exp(2*b*vYF(iNumVol)) - exp(2*b*dTprev));        
        V12  = V12 + sigma(iNumVol)*eta(iNumVol)*(exp((a+b)*vYF(iNumVol)) - exp((a+b)*dTprev));
        dTprev = vYF(iNumVol);
        if iNumVol== iNbVolDates
            break;
        else
            iNumVol= iNumVol + 1;
        end        
    end
    if dTprev < dExpiry
        V1  = V1 + sigma(iNumVol)*sigma(iNumVol)*(exp(2*a*dExpiry) - exp(2*a*dTprev));        
        V2  = V2 + eta(iNumVol)*eta(iNumVol)*(exp(2*b*dExpiry) - exp(2*b*dTprev));        
        V12  = V12 + sigma(iNumVol)*eta(iNumVol)*(exp((a+b)*dExpiry) - exp((a+b)*dTprev));        
    end    
    V = dA1*dA1*V1/2/a + dA2*dA2*V2/2/b + 2*rho*dA1*dA2*V12/(a+b);
end
%assuming constant volatility
function V = VarSwap(a,b,sigma,eta,rho,vTSwp,P0,LVL, dSwpRate)
    dExpiry = vTSwp(1);
    vDiff = diff(vTSwp)'; %line
    
    vExp1 = exp(-a*vTSwp);
    dA1 = vExp1(1)*P0(1) - vExp1(end)*P0(end) - dSwpRate*vDiff*(vExp1(2:end).*P0(2:end));
    dA1 = dA1 /(a*LVL);
    vExp2 = exp(-b*vTSwp);
    dA2 = vExp2(1)*P0(1) - vExp2(end)*P0(end) - dSwpRate*vDiff*(vExp2(2:end).*P0(2:end));
    dA2 = dA2 /(b*LVL);
        
    V  = dA1*dA1*sigma*sigma*(exp(2*a*dExpiry) - 1)/(2*a);        
    V  = V + dA2*dA2*eta*eta*(exp(2*b*dExpiry) - 1)/(2*b);        
    V  = V + 2*rho*dA1*dA2*sigma*eta*(exp((a+b)*dExpiry) - 1)/(a+b);
end