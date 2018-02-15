function closereqmod(varargin)
%CLOSEREQMOD  Figure close request function modified by Maxence Boutet
% This modified function prevents the figure from being closed.

if isempty(gcbf)
    if length(dbstack) == 1
        warning(message('MATLAB:closereq:ObsoleteUsage'));
    end
    fprintf('\nClosing figure while NEAT is running is inhibited.\n');
    % close('force');
else
    fprintf('\nClosing figure while NEAT is running is inhibited.\n');
    % delete(gcbf);
end
