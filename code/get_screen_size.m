function [screenSizeX, screenSizeY, varargout] = get_screen_size()

% Get current units
oldUnits = get(0,'Units');

% Set current units to pixels
set(0,'Units','pixels');

% Get screen size in pixels
screenSize = get(0,'Screensize');
screenSizeX = screenSize(3);
screenSizeY = screenSize(4);

% Set units back to old unit
set(0,'Units',oldUnits);

if nargout == 3
    varargout{1} = oldUnits;
end

end
