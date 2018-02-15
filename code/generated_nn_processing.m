function generated_nn_processing(funPath, funName, nPreviousStates)
% This function add lines of code at the beginning and end of the generated
% neural network functions in order to store into matrices the previous
% recurrent layer states in place of havin a seperate input for every state.

% Read in the file as binary and convert to chars.
if nPreviousStates >= 1
    fid = fopen(funPath);
    text = fread(fid, inf, '*char')';
    fclose(fid);
end

% Find recurrent layer ID
if nPreviousStates >= 1
    if nPreviousStates > 1
        out1 = regexp(text,'ai(\d*),','tokens');
        out2 = regexp(text,',ai(\d*)\)','tokens');
    else
        out1 = {{}};
        out2 = regexp(text,',ai(\d*)\)','tokens');
    end
    out1 = [out1{:}];
    out2 = [out2{:}];
    out = [out1,out2];
    out = unique(out);
    recurrLayerId = sort(str2double(out));
end

if nPreviousStates >= 1
    newCodeStart = [];
    % newCodeEnd = ['afMat = zeros\(1,',num2str(numel(recurrLayerId)),'\);\n'];
    newCodeEnd = ['afMat = zeros\(',num2str(numel(recurrLayerId)),',1\);\n'];
    % Form code to add in function
    for iLayer = 1:numel(recurrLayerId)
        % Code to add at the beginning of the file
        newCodeStart = [...
            newCodeStart,...
            ['ai',num2str(recurrLayerId(iLayer)),' = aiMat\(',num2str(iLayer),'\);\n']];
        % Code to add at the end of the file
        newCodeEnd = [...
            newCodeEnd,...
            ['afMat\(',num2str(iLayer),'\) = af',num2str(recurrLayerId(iLayer)),';\n']];
    end
else
    newCodeStart = 'afMat = aiMat;\n';
end

% Replace ai1,ai2,etc. by aiMat
if nPreviousStates >= 1
    find_and_replace(funPath,',ai(.*?)\)', ',aiMat\)');
else
    find_and_replace(funPath,[funName,'\((.*?)\)'],[funName,'\($1,aiMat\)']);
end

% Replace af1,af2,etc. by afMat
if nPreviousStates >= 1
    find_and_replace(funPath,',af(.*?)\]', ',afMat\]');
else
    find_and_replace(funPath,['\[(.*?)\] = ',funName],['\[$1,afMat\] = ',funName]);
end

% Add code at the beginning
find_and_replace(funPath,...
    '% ===== NEURAL NETWORK CONSTANTS =====\n',...
    ['% ===== NEURAL NETWORK CONSTANTS =====\n\n',newCodeStart]);

% Add code at the end
if nPreviousStates >= 1
    find_and_replace(funPath,...
        'end\n\n% ===',...
        [newCodeEnd,'end\n\n% ===']);
end

% Replace activation function for Kenneth Stanley version
find_and_replace(funPath, 'exp\(-n\)', 'exp(-4.9*n)');

end
