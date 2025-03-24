function  [Fvalue,pvalue] = Fstat(Xresid, Xorig, p, n)
    SStotal = sum(Xorig - mean(Xorig).^2);
    SS_residual = sum((Xresid).^2);
    Fvalue = ((SStotal-SS_residual)/p) / (SS_residual / (n - p - 1));
    pvalue = 1 - fcdf(Fvalue, p, n - p - 1);
    
end