function [tightTextLine] = TightenCifTextLine(looseTextLine)
%TightenCifTextLine() tightens the scanned text line from CIF file.
% NOTE:
%   This function serves to remove the leading and trailing spaces and
%   spaces within single quote marks.

tightTextLine = looseTextLine;
if ischar(tightTextLine)
    % delete spaces in the front of the text line
    tightTextLine = strtrim(tightTextLine);
    
    % delete parentheses and the content inside
    if contains(tightTextLine, '(') && contains(tightTextLine, ')')
        leftParenIdx = strfind(tightTextLine, '(');
        rightParenIdx = strfind(tightTextLine, ')');
        tightTextLine(leftParenIdx : rightParenIdx) = [];
    end
    
    % find single quote marks
    singleQuoteIndices = strfind(tightTextLine, '''');
    if length(singleQuoteIndices) == 2
        if singleQuoteIndices(2) - singleQuoteIndices(1) >= 2
            quotedText = tightTextLine(singleQuoteIndices(1) + 1 : singleQuoteIndices(2) - 1);
            quotedText = quotedText(~isspace(quotedText));
            tightTextLine = [tightTextLine(1 : singleQuoteIndices(1) - 1),...
                quotedText];
        end
    end
end

end

