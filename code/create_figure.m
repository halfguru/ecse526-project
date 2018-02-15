function figHandle = create_figure(resolutionX, resolutionY, visibility, varargin)

% resolutionX :
% - width in pixel
%
% resolutionY :
% - height in pixel
%
% visibility:
% - 'on'
% - 'off'

% Get screen resolution in pixels
[screenSizeX,screenSizeY,oldUnits] = get_screen_size();

% Set current units to pixels
set(0,'Units','pixels');

% Create figure
if nargin == 4
    figHandle = figure(varargin{1});
else    
    figHandle = figure;
end

% Figure visibility
figHandle.Visible = visibility;

% Figure position and size
positionX = 0.5*screenSizeX-0.5*resolutionX;
positionY = 0.5*screenSizeY-0.5*resolutionY;
figHandle.Position = [positionX,positionY,resolutionX,resolutionY];

% Figure background color
figHandle.Color = 'w'; % none

% Set renderer
% figHandle.AlignVertexCenters = 'on';

% Set units back to old unit
set(0,'Units',oldUnits);

end
