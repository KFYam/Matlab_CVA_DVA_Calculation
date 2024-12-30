function price = yyPriceFXCall(vVectSpot,vInCurve_d,vInCurve_f,vMat,HW2Model_for)
    iNbCall = length(vVectSpot);   
    vVar = zeros(iNbCall,1);
    
    for iNumCall = 1:iNbCall
        %pricing
        vVar(iNumCall) = yyVarFX(HW2Model_for,0, vMat(iNumCall));
        if vVar(iNumCall)<0
            fprintf('Negative vol for FX call pricing\n');
            vVar(iNumCall) = 0.0000001;
        end 
    end 
    vP0T_d  = discountFactorByTime(vInCurve_d(1), vMat')';
    vP0T_f  = discountFactorByTime(vInCurve_f(1), vMat')';
%    r_d = intenvset(getRateSpec(vInCurve_d(1)),'EndDates',vInCurve_d(1).referenceDate + 365*vMat);  
%    r_f = intenvset(getRateSpec(vInCurve_f(1)),'EndDates',vInCurve_f(1).referenceDate + 365*vMat);  
%     vP0T_d = intenvget(r_d,'Disc');  
%     vP0T_f = intenvget(r_f,'Disc');  
    dD1 = (log(vP0T_f ./vP0T_d) + 0.5*vVar) ./ sqrt(vVar);
    dD2 = (log(vP0T_f ./vP0T_d) - 0.5*vVar) ./ sqrt(vVar);
    price = vVectSpot .* (vP0T_f .*normcdf(dD1) - vP0T_d .*normcdf(dD2));
end