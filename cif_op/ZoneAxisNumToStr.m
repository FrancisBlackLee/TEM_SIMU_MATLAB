function [str] = ZoneAxisNumToStr(zoneAxis)
%ZoneAxisNumToStr.m convert the input array for zone axis to a string that
%can be interpreted with latex.

str = ['$', SingleNumToStr(zoneAxis(1)), SingleNumToStr(zoneAxis(2)),...
    SingleNumToStr(zoneAxis(3)), '$'];

    function [singleStr] = SingleNumToStr(singleNum)
        if singleNum < 0
            singleStr = ['\bar{', num2str(-singleNum), '}'];
        else
            singleStr = num2str(singleNum);
        end
    end

end