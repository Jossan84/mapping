function [tInt, tDec] = tToIntegerAndDecimal(t, maxInteger)
    % Perform truncation and cyclic function
    tTr = min(mod(t, 2 * maxInteger),...
        mod(2 * maxInteger - t, 2 * maxInteger));
    % Extracs integer and decimal part
    tInt = max(maxInteger - fix((maxInteger + 1) - tTr), 0);
    tDec = tTr - tInt;
       
    % Zero-base indexing to One-base indexing
    tInt = tInt + 1;
    
    % Data Type Convert
    %tInt = int16(tInt);
end

