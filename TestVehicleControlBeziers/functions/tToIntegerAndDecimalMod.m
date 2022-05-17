function [tInt, tDec] = tToIntegerAndDecimalMod(t, maxInteger,t_c0)
    % Perform truncation and cyclic function   
    tTr = triangularWave(t-t_c0,2*(maxInteger - t_c0),maxInteger - t_c0) + t_c0;
    % Extracs integer and decimal part
    tInt = max(maxInteger - fix((maxInteger + 1) - tTr), 0);
    tDec = tTr - tInt;
    
    % Zero-base indexing to One-base indexing
    %tInt = tInt + 1;
    
    % Data Type Convert
    %tInt = int16(tInt);
end

