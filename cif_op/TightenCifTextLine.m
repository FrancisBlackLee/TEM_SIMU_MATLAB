function [tightTextLine] = TightenCifTextLine(looseTextLine)
%TightenCifTextLine() tightens the scanned text line from CIF file.
% NOTE:
%   This function serves to remove the leading and trailing spaces and
%   spaces within single quote marks.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2021  Francis Black Lee and Li Xian

%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   any later version.

%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.

%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <https://www.gnu.org/licenses/>.

%   Email: warner323@outlook.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tightTextLine = looseTextLine;
if ischar(tightTextLine)
    % delete spaces in the front of the text line
    tightTextLine = strtrim(tightTextLine);
    
    % delete parentheses and the content inside 
    while contains(tightTextLine, '(') && contains(tightTextLine, ')')
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

