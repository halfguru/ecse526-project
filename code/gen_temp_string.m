function str = gen_temp_string(nChar)

possibleChar = 'ABCDEFGHIJKLMNOP01234569789';

str = possibleChar(randi([1,27],[1,nChar]));

% str = repmat(' ',[1,nChar]);
% for iChar = 1:nChar
%     str(iChar) = possibleChar(randi([1,27]));
% end

end
