function [figHandle, figHandleError] = plot_uav_sim(Data, showSignal, showReference, ...
    showTolerance, showInput, useSeparatePlot, folderName, fileName, figureVisibility, ...
    exportEn, customLegend, plotError, varargin)
% Required fields in Data (for each UAV)
% Data(iUav).t
% Data(iUav).xTol              (if showTolerance)
% Data(iUav).yTol              (if showTolerance)
% Data(iUav).zTol              (if showTolerance)
% Data(iUav).psiTol            (if showTolerance)
% Data(iUav).thetaTol          (if showTolerance)
% Data(iUav).phiTol            (if showTolerance)
% Data(iUav).xDesired          (if showTolerance/showReference)
% Data(iUav).yDesired          (if showTolerance/showReference)
% Data(iUav).zDesired          (if showTolerance/showReference)
% Data(iUav).psiDesired        (if showTolerance/showReference)
% Data(iUav).thetaDesired      (if showTolerance/showReference)
% Data(iUav).phiDesired        (if showTolerance/showReference)
% Data(iUav).x
% Data(iUav).y
% Data(iUav).z
% Data(iUav).psi
% Data(iUav).theta
% Data(iUav).phi
% Data(iUav).Va1                (if showInput)
% Data(iUav).Va2                (if showInput)
% Data(iUav).Va3                (if showInput)
% Data(iUav).Va4                (if showInput)

%% Init

% Number of UAVs
nUavs = numel(Data);

% Constants
DEG_2_RAD = pi/180;
RAD_2_DEG = 180/pi;
MM_2_M = 1/1000;
M_2_MM = 1000/1;

%% Initialize figures

% Resolution
if showInput
    resolutionX = 800;
    resolutionY = 800;
else
    resolutionX = 800;
    resolutionY = round(800*(2.85/3));
end

% Signal plot
if useSeparatePlot
    for iUav = 1:nUavs
        if figureVisibility
            figHandle(iUav) = create_figure(resolutionX, resolutionY, 'on');
        else
            figHandle(iUav) = create_figure(resolutionX, resolutionY, 'off');
        end
    end
else
    if figureVisibility
        if nargin == 13 || nargin == 14
            figHandle = varargin{1};
        else
            figHandle = create_figure(resolutionX, resolutionY, 'on');
        end
    else
        if nargin == 13 || nargin == 14
            figHandle = varargin{1};
        else
            figHandle = create_figure(resolutionX, resolutionY, 'off');
        end
    end
end

% Error plot
resolutionX = 800;
resolutionY = round(800*(2.85/3));
if plotError
    if useSeparatePlot
        for iUav = 1:nUavs
            if figureVisibility
                figHandleError(iUav) = create_figure(resolutionX, resolutionY, 'on');
            else
                figHandleError(iUav) = create_figure(resolutionX, resolutionY, 'off');
            end
        end
    else
        if figureVisibility
            if nargin == 14
                figHandleError = varargin{2};
            else
                figHandleError = create_figure(resolutionX, resolutionY, 'on');
            end
        else
            if nargin == 14
                figHandleError = varargin{2};
            else
                figHandleError = create_figure(resolutionX, resolutionY, 'off');
            end
        end
    end
end

%% Create axes

if showInput
    subplotConfig = {3,3};
else
    subplotConfig = {2,3};
end

tickFontSize = 8;

marginSubPlot = 0.07;
if useSeparatePlot
    for iUav = 1:nUavs
        % 1st row - x, y, z
        set(0, 'currentfigure', figHandle(iUav));
        AxHandle(iUav).subPlot1 = subplot_tight(subplotConfig{:},1,marginSubPlot);
        AxHandle(iUav).subPlot1.XAxis.FontSize = tickFontSize;
        AxHandle(iUav).subPlot1.YAxis.FontSize = tickFontSize;
        set(0, 'currentfigure', figHandle(iUav));
        AxHandle(iUav).subPlot2 = subplot_tight(subplotConfig{:},2,marginSubPlot);
        AxHandle(iUav).subPlot2.XAxis.FontSize = tickFontSize;
        AxHandle(iUav).subPlot2.YAxis.FontSize = tickFontSize;
        set(0, 'currentfigure', figHandle(iUav));
        AxHandle(iUav).subPlot3 = subplot_tight(subplotConfig{:},3,marginSubPlot);
        AxHandle(iUav).subPlot3.XAxis.FontSize = tickFontSize;
        AxHandle(iUav).subPlot3.YAxis.FontSize = tickFontSize;
        % 2nd row - psi, theta, phi
        set(0, 'currentfigure', figHandle(iUav));
        AxHandle(iUav).subPlot4 = subplot_tight(subplotConfig{:},4,marginSubPlot);
        AxHandle(iUav).subPlot4.XAxis.FontSize = tickFontSize;
        AxHandle(iUav).subPlot4.YAxis.FontSize = tickFontSize;
        set(0, 'currentfigure', figHandle(iUav));
        AxHandle(iUav).subPlot5 = subplot_tight(subplotConfig{:},5,marginSubPlot);
        AxHandle(iUav).subPlot5.XAxis.FontSize = tickFontSize;
        AxHandle(iUav).subPlot5.YAxis.FontSize = tickFontSize;
        set(0, 'currentfigure', figHandle(iUav));
        AxHandle(iUav).subPlot6 = subplot_tight(subplotConfig{:},6,marginSubPlot);
        AxHandle(iUav).subPlot6.XAxis.FontSize = tickFontSize;
        AxHandle(iUav).subPlot6.YAxis.FontSize = tickFontSize;
        % 3rd row - Va
        if showInput
            set(0, 'currentfigure', figHandle(iUav));
            AxHandle(iUav).subPlot7 = subplot_tight(subplotConfig{:},[7,8,9],marginSubPlot);
            AxHandle(iUav).subPlot7.XAxis.FontSize = tickFontSize;
            AxHandle(iUav).subPlot7.YAxis.FontSize = tickFontSize;
        end
    end
else
    % 1st row - x, y, z
    set(0, 'currentfigure', figHandle);
    AxHandle(1).subPlot1 = subplot_tight(subplotConfig{:},1,marginSubPlot);
    AxHandle(1).subPlot1.XAxis.FontSize = tickFontSize;
    AxHandle(1).subPlot1.YAxis.FontSize = tickFontSize;
    set(0, 'currentfigure', figHandle);
    AxHandle(1).subPlot2 = subplot_tight(subplotConfig{:},2,marginSubPlot);
    AxHandle(1).subPlot2.XAxis.FontSize = tickFontSize;
    AxHandle(1).subPlot2.YAxis.FontSize = tickFontSize;
    set(0, 'currentfigure', figHandle);
    AxHandle(1).subPlot3 = subplot_tight(subplotConfig{:},3,marginSubPlot);
    AxHandle(1).subPlot3.XAxis.FontSize = tickFontSize;
    AxHandle(1).subPlot3.YAxis.FontSize = tickFontSize;
    % 2nd row - psi, theta, phi
    set(0, 'currentfigure', figHandle);
    AxHandle(1).subPlot4 = subplot_tight(subplotConfig{:},4,marginSubPlot);
    AxHandle(1).subPlot4.XAxis.FontSize = tickFontSize;
    AxHandle(1).subPlot4.YAxis.FontSize = tickFontSize;
    set(0, 'currentfigure', figHandle);
    AxHandle(1).subPlot5 = subplot_tight(subplotConfig{:},5,marginSubPlot);
    AxHandle(1).subPlot5.XAxis.FontSize = tickFontSize;
    AxHandle(1).subPlot5.YAxis.FontSize = tickFontSize;
    set(0, 'currentfigure', figHandle);
    AxHandle(1).subPlot6 = subplot_tight(subplotConfig{:},6,marginSubPlot);
    AxHandle(1).subPlot6.XAxis.FontSize = tickFontSize;
    AxHandle(1).subPlot6.YAxis.FontSize = tickFontSize;
    % 3rd row - Va
    if showInput
        set(0, 'currentfigure', figHandle);
        AxHandle(1).subPlot7 = subplot_tight(subplotConfig{:},[7,8,9],marginSubPlot);
        AxHandle(1).subPlot7.XAxis.FontSize = tickFontSize;
        AxHandle(1).subPlot7.YAxis.FontSize = tickFontSize;
    end
end

subplotConfig = {2,3};
if useSeparatePlot
    for iUav = 1:nUavs
        % 1st row - Errors x, y, z
        set(0, 'currentfigure', figHandleError(iUav));
        AxHandleError(iUav).subPlot1 = subplot_tight(subplotConfig{:},1,marginSubPlot);
        AxHandleError(iUav).subPlot1.XAxis.FontSize = tickFontSize;
        AxHandleError(iUav).subPlot1.YAxis.FontSize = tickFontSize;
        set(0, 'currentfigure', figHandleError(iUav));
        AxHandleError(iUav).subPlot2 = subplot_tight(subplotConfig{:},2,marginSubPlot);
        AxHandleError(iUav).subPlot2.XAxis.FontSize = tickFontSize;
        AxHandleError(iUav).subPlot2.YAxis.FontSize = tickFontSize;
        set(0, 'currentfigure', figHandleError(iUav));
        AxHandleError(iUav).subPlot3 = subplot_tight(subplotConfig{:},3,marginSubPlot);
        AxHandleError(iUav).subPlot3.XAxis.FontSize = tickFontSize;
        AxHandleError(iUav).subPlot3.YAxis.FontSize = tickFontSize;
        % 2nd row - Errors psi, theta, phi
        set(0, 'currentfigure', figHandleError(iUav));
        AxHandleError(iUav).subPlot4 = subplot_tight(subplotConfig{:},4,marginSubPlot);
        AxHandleError(iUav).subPlot4.XAxis.FontSize = tickFontSize;
        AxHandleError(iUav).subPlot4.YAxis.FontSize = tickFontSize;
        set(0, 'currentfigure', figHandleError(iUav));
        AxHandleError(iUav).subPlot5 = subplot_tight(subplotConfig{:},5,marginSubPlot);
        AxHandleError(iUav).subPlot5.XAxis.FontSize = tickFontSize;
        AxHandleError(iUav).subPlot5.YAxis.FontSize = tickFontSize;
        set(0, 'currentfigure', figHandleError(iUav));
        AxHandleError(iUav).subPlot6 = subplot_tight(subplotConfig{:},6,marginSubPlot);
        AxHandleError(iUav).subPlot6.XAxis.FontSize = tickFontSize;
        AxHandleError(iUav).subPlot6.YAxis.FontSize = tickFontSize;
    end
else
    % 1st row - Errors x, y, z
    set(0, 'currentfigure', figHandleError);
    AxHandleError(1).subPlot1 = subplot_tight(subplotConfig{:},1,marginSubPlot);
    AxHandleError(1).subPlot1.XAxis.FontSize = tickFontSize;
    AxHandleError(1).subPlot1.YAxis.FontSize = tickFontSize;
    set(0, 'currentfigure', figHandleError);
    AxHandleError(1).subPlot2 = subplot_tight(subplotConfig{:},2,marginSubPlot);
    AxHandleError(1).subPlot2.XAxis.FontSize = tickFontSize;
    AxHandleError(1).subPlot2.YAxis.FontSize = tickFontSize;
    set(0, 'currentfigure', figHandleError);
    AxHandleError(1).subPlot3 = subplot_tight(subplotConfig{:},3,marginSubPlot);
    AxHandleError(1).subPlot3.XAxis.FontSize = tickFontSize;
    AxHandleError(1).subPlot3.YAxis.FontSize = tickFontSize;
    % 2nd row - Errors psi, theta, phi
    set(0, 'currentfigure', figHandleError);
    AxHandleError(1).subPlot4 = subplot_tight(subplotConfig{:},4,marginSubPlot);
    AxHandleError(1).subPlot4.XAxis.FontSize = tickFontSize;
    AxHandleError(1).subPlot4.YAxis.FontSize = tickFontSize;
    set(0, 'currentfigure', figHandleError);
    AxHandleError(1).subPlot5 = subplot_tight(subplotConfig{:},5,marginSubPlot);
    AxHandleError(1).subPlot5.XAxis.FontSize = tickFontSize;
    AxHandleError(1).subPlot5.YAxis.FontSize = tickFontSize;
    set(0, 'currentfigure', figHandleError);
    AxHandleError(1).subPlot6 = subplot_tight(subplotConfig{:},6,marginSubPlot);
    AxHandleError(1).subPlot6.XAxis.FontSize = tickFontSize;
    AxHandleError(1).subPlot6.YAxis.FontSize = tickFontSize;
end

%% Compute tolerance curves

if showTolerance
    xLower      = cell(1,nUavs);
    xUpper      = cell(1,nUavs);
    yLower      = cell(1,nUavs);
    yUpper      = cell(1,nUavs);
    zLower      = cell(1,nUavs);
    zUpper      = cell(1,nUavs);
    psiLower    = cell(1,nUavs);
    psiUpper    = cell(1,nUavs);
    thetaLower  = cell(1,nUavs);
    thetaUpper  = cell(1,nUavs);
    phiLower    = cell(1,nUavs);
    phiUpper    = cell(1,nUavs);
    for iUav = 1:nUavs
        [xLower{iUav}, xUpper{iUav}] = compute_lower_upper_bound(...
            Data(iUav).xTol(1), Data(iUav).xTol(2), Data(iUav).t, Data(iUav).xDesired);
        [yLower{iUav}, yUpper{iUav}] = compute_lower_upper_bound(...
            Data(iUav).yTol(1), Data(iUav).yTol(2), Data(iUav).t, Data(iUav).yDesired);
        [zLower{iUav}, zUpper{iUav}] = compute_lower_upper_bound(...
            Data(iUav).zTol(1), Data(iUav).zTol(2), Data(iUav).t, Data(iUav).zDesired);
        [psiLower{iUav}, psiUpper{iUav}] = compute_lower_upper_bound(...
            Data(iUav).psiTol(1), Data(iUav).psiTol(2), Data(iUav).t, Data(iUav).psiDesired);
        [thetaLower{iUav}, thetaUpper{iUav}] = compute_lower_upper_bound(...
            Data(iUav).thetaTol(1), Data(iUav).thetaTol(2), Data(iUav).t, Data(iUav).thetaDesired);
        [phiLower{iUav}, phiUpper{iUav}] = compute_lower_upper_bound(...
            Data(iUav).phiTol(1), Data(iUav).phiTol(2), Data(iUav).t, Data(iUav).phiDesired);
    end
end

%% Compute errors

if plotError
    xError      = cell(1,nUavs);
    yError      = cell(1,nUavs);
    zError      = cell(1,nUavs);
    psiError    = cell(1,nUavs);
    thetaError  = cell(1,nUavs);
    phiError    = cell(1,nUavs);
    for iUav = 1:nUavs
        xError{iUav}     = Data(iUav).x - Data(iUav).xDesired;
        yError{iUav}     = Data(iUav).y - Data(iUav).yDesired;
        zError{iUav}     = Data(iUav).z - Data(iUav).zDesired;
        psiError{iUav}   = Data(iUav).psi - Data(iUav).psiDesired;
        thetaError{iUav} = Data(iUav).theta - Data(iUav).thetaDesired;
        phiError{iUav}   = Data(iUav).phi - Data(iUav).phiDesired;
    end
end

%% Colors settings

fontSize = 10;
fontSizeLegend = 7;
if useSeparatePlot
    lineSpecSignal = {'LineStyle','-','LineWidth',2};
    lineSpecReference = {'LineStyle','--','LineWidth',1.5};
    lineSpecTolerance = {'LineStyle','-.','LineWidth',1.5};
else
    lineSpecSignal = {'LineStyle','-','LineWidth',1};
    lineSpecReference = {'LineStyle','--'};
    lineSpecTolerance = {'LineStyle','-.'};
end
colorListSignal = [...
    [255,0,0];      % red
    [255,128,0];    % orange
    [51,102,0];     % green
    [0,153,153];    % light blue
    [0,0,255];      % dark blue
    [127,0,255];    % violet
    [153,0,76];     % pink
    ]/255;
colorListReference = 1 + (colorListSignal-1)*0.9;
colorListTolerance = 1 + (colorListSignal-1)*0.8;

if nUavs > 7
    for iColor = 1:(nUavs-7)
        colorListSignal = [...
            colorListSignal;
            randi([0,150],[1,3])/255];
        % http://stackoverflow.com/a/12228643
        colorListReference = [...
            colorListReference;
            1 + (colorListSignal(end,:)-1)*0.7];
        colorListTolerance = [...
            colorListTolerance;
            1 + (colorListSignal(end,:)-1)*0.4];
    end
end

%% Plotting

legendStr1 = cell(1,nUavs);
for iUav = 1:nUavs
    
    %% 1st row - x
    
    % Get figure and axis
    if useSeparatePlot
        set(figHandle(iUav), 'CurrentAxes', AxHandle(iUav).subPlot1);
        iUav_ = iUav;
    else
        set(figHandle(1), 'CurrentAxes', AxHandle(1).subPlot1);
        iUav_ = 1;
    end
    
    hold(AxHandle(iUav_).subPlot1,'on');
    
    % Plot signal
    if showSignal
        hLine55(iUav) = plot(AxHandle(iUav_).subPlot1,...
            Data(iUav).t, Data(iUav).x, lineSpecSignal{:}, ...
            'Color', 'black');
        hLine1(iUav) = plot(AxHandle(iUav_).subPlot1,...
            Data(iUav).t, Data(iUav).x, lineSpecSignal{:}, ...
            'Color', colorListSignal(iUav,:));
    end
    
    % Plot tolerance curves
    if showTolerance && showSignal
        hLine56(iUav) = plot(AxHandle(iUav_).subPlot1,...
            Data(iUav).t, xLower{iUav}, lineSpecTolerance{:}, ...
            'Color', 'black');
        plot(AxHandle(iUav_).subPlot1,...
            Data(iUav).t, xLower{iUav}, lineSpecTolerance{:}, ...
            'Color', colorListTolerance(iUav,:));
        plot(AxHandle(iUav_).subPlot1,...
            Data(iUav).t, xUpper{iUav}, lineSpecTolerance{:}, ...
            'Color', colorListTolerance(iUav,:));
    elseif showTolerance && ~showSignal
        hLine56(iUav) = plot(AxHandle(iUav_).subPlot1,...
            Data(iUav).t, xLower{iUav}, lineSpecTolerance{:}, ...
            'Color', 'black');
        hLine1(iUav) = plot(AxHandle(iUav_).subPlot1,...
            Data(iUav).t, xLower{iUav}, lineSpecTolerance{:}, ...
            'Color', colorListTolerance(iUav,:));
        plot(AxHandle(iUav_).subPlot1,...
            Data(iUav).t, xUpper{iUav}, lineSpecTolerance{:}, ...
            'Color', colorListTolerance(iUav,:));
    end
    
    % Plot reference curves
    if showReference
        if  showSignal || showTolerance
            hLine57(iUav) = plot(AxHandle(iUav_).subPlot1,...
                Data(iUav).t, Data(iUav).xDesired, lineSpecReference{:}, ...
                'Color', 'black');
            plot(AxHandle(iUav_).subPlot1,...
                Data(iUav).t, Data(iUav).xDesired, lineSpecReference{:}, ...
                'Color', colorListReference(iUav,:));
        elseif ~showSignal && ~showTolerance
            hLine57(iUav) = plot(AxHandle(iUav_).subPlot1,...
                Data(iUav).t, Data(iUav).xDesired, lineSpecReference{:}, ...
                'Color', 'black');
            hLine1(iUav) = plot(AxHandle(iUav_).subPlot1,...
                Data(iUav).t, Data(iUav).xDesired, lineSpecReference{:}, ...
                'Color', colorListReference(iUav,:));
        end
    end
    
    % Labels and grid
    ylabel(AxHandle(iUav_).subPlot1,'x (m)','FontSize', fontSize);
    xlabel(AxHandle(iUav_).subPlot1,'t (s)','FontSize', fontSize);
    grid(AxHandle(iUav_).subPlot1,'on');
    
    % Legend string
    if isempty(customLegend)
        legendStr1{iUav} = ['UAV ',num2str(iUav)];
    else
        legendStr1{iUav} = customLegend{iUav};
    end
    
    if showSignal
        uistack([hLine55(iUav)], 'bottom');
    end
    if showTolerance
        uistack([hLine56(iUav)], 'bottom');
    end
    if showReference
        uistack([hLine57(iUav)], 'bottom');
    end
    
    hold(AxHandle(iUav_).subPlot1,'off');
    
    %% 1st row - y
    
    % Get figure and axis
    if useSeparatePlot
        set(figHandle(iUav), 'CurrentAxes', AxHandle(iUav).subPlot2);
        iUav_ = iUav;
    else
        set(figHandle(1), 'CurrentAxes', AxHandle(1).subPlot2);
        iUav_ = 1;
    end
    
    hold(AxHandle(iUav_).subPlot2,'on');
    
    % Plot signal
    if showSignal
        plot(AxHandle(iUav_).subPlot2,...
            Data(iUav).t, Data(iUav).y, lineSpecSignal{:}, ...
            'Color', colorListSignal(iUav,:));
    end
    
    % Plot tolerance curves
    if showTolerance
        plot(AxHandle(iUav_).subPlot2,...
            Data(iUav).t, yLower{iUav}, lineSpecTolerance{:}, ...
            'Color', colorListTolerance(iUav,:));
        plot(AxHandle(iUav_).subPlot2,...
            Data(iUav).t, yUpper{iUav}, lineSpecTolerance{:}, ...
            'Color', colorListTolerance(iUav,:));
    end
    
    % Plot reference curves
    if showReference
        plot(AxHandle(iUav_).subPlot2,...
            Data(iUav).t, Data(iUav).yDesired, lineSpecReference{:}, ...
            'Color', colorListReference(iUav,:));
    end
    ylabel(AxHandle(iUav_).subPlot2,'y (m)','FontSize', fontSize);
    xlabel(AxHandle(iUav_).subPlot2,'t (s)','FontSize', fontSize);
    grid(AxHandle(iUav_).subPlot2,'on');
    
    hold(AxHandle(iUav_).subPlot2,'off');
    
    %% 1st row - z
    
    % Get figure and axis
    if useSeparatePlot
        set(figHandle(iUav), 'CurrentAxes', AxHandle(iUav).subPlot3);
        iUav_ = iUav;
    else
        set(figHandle(1), 'CurrentAxes', AxHandle(1).subPlot3);
        iUav_ = 1;
    end
    
    hold(AxHandle(iUav_).subPlot3,'on');
    
    % Plot signal
    if showSignal
        plot(AxHandle(iUav_).subPlot3,...
            Data(iUav).t, Data(iUav).z, lineSpecSignal{:}, ...
            'Color', colorListSignal(iUav,:));
    end
    
    % Plot tolerance curves
    if showTolerance
        plot(AxHandle(iUav_).subPlot3,...
            Data(iUav).t, zLower{iUav}, lineSpecTolerance{:}, ...
            'Color', colorListTolerance(iUav,:));
        plot(AxHandle(iUav_).subPlot3,...
            Data(iUav).t, zUpper{iUav}, lineSpecTolerance{:}, ...
            'Color', colorListTolerance(iUav,:));
    end
    
    % Plot reference curves
    if showReference
        plot(AxHandle(iUav_).subPlot3,...
            Data(iUav).t, Data(iUav).zDesired, lineSpecReference{:}, ...
            'Color', colorListReference(iUav,:));
    end
    ylabel(AxHandle(iUav_).subPlot3,'z (m)','FontSize', fontSize);
    xlabel(AxHandle(iUav_).subPlot3,'t (s)','FontSize', fontSize);
    grid(AxHandle(iUav_).subPlot3,'on');
    
    hold(AxHandle(iUav_).subPlot3,'off');
    
    %% 2nd row - psi
    
    % Get figure and axis
    if useSeparatePlot
        set(figHandle(iUav), 'CurrentAxes', AxHandle(iUav).subPlot4);
        iUav_ = iUav;
    else
        set(figHandle(1), 'CurrentAxes', AxHandle(1).subPlot4);
        iUav_ = 1;
    end
    
    hold(AxHandle(iUav_).subPlot4,'on');
    
    % Plot signal
    if showSignal
        plot(AxHandle(iUav_).subPlot4,...
            Data(iUav).t, Data(iUav).psi*RAD_2_DEG, lineSpecSignal{:}, ...
            'Color', colorListSignal(iUav,:));
    end
    
    % Plot tolerance curves
    if showTolerance
        plot(AxHandle(iUav_).subPlot4,...
            Data(iUav).t, psiLower{iUav}*RAD_2_DEG, lineSpecTolerance{:}, ...
            'Color', colorListTolerance(iUav,:));
        plot(AxHandle(iUav_).subPlot4,...
            Data(iUav).t, psiUpper{iUav}*RAD_2_DEG, lineSpecTolerance{:}, ...
            'Color', colorListTolerance(iUav,:));
    end
    
    % Plot reference curves
    if showReference
        plot(AxHandle(iUav_).subPlot4,...
            Data(iUav).t, Data(iUav).psiDesired*RAD_2_DEG, ...
            lineSpecReference{:}, ...
            'Color', colorListReference(iUav,:));
    end
    ylabel(AxHandle(iUav_).subPlot4,'\psi (deg)','FontSize', fontSize);
    if showInput || useSeparatePlot
        xlabel(AxHandle(iUav_).subPlot4,'t (s)','FontSize', fontSize);
    end
    grid(AxHandle(iUav_).subPlot4,'on');
    
    hold(AxHandle(iUav_).subPlot4,'off');
    
    %% 2nd row - theta
    
    % Get figure and axis
    if useSeparatePlot
        set(figHandle(iUav), 'CurrentAxes', AxHandle(iUav).subPlot5);
        iUav_ = iUav;
    else
        set(figHandle(1), 'CurrentAxes', AxHandle(1).subPlot5);
        iUav_ = 1;
    end
    
    hold(AxHandle(iUav_).subPlot5,'on');
    
    % Plot signal
    if showSignal
        plot(AxHandle(iUav_).subPlot5,...
            Data(iUav).t, Data(iUav).theta*RAD_2_DEG, lineSpecSignal{:}, ...
            'Color', colorListSignal(iUav,:));
    end
    
    % Plot tolerance curves
    if showTolerance
        plot(AxHandle(iUav_).subPlot5,...
            Data(iUav).t, thetaLower{iUav}*RAD_2_DEG, lineSpecTolerance{:}, ...
            'Color', colorListTolerance(iUav,:));
        plot(AxHandle(iUav_).subPlot5,...
            Data(iUav).t, thetaUpper{iUav}*RAD_2_DEG, lineSpecTolerance{:}, ...
            'Color', colorListTolerance(iUav,:));
    end
    
    % Plot reference curves
    if showReference
        plot(AxHandle(iUav_).subPlot5,...
            Data(iUav).t, Data(iUav).thetaDesired*RAD_2_DEG, ...
            lineSpecReference{:}, ...
            'Color', colorListReference(iUav,:));
    end
    ylabel(AxHandle(iUav_).subPlot5,'\theta (deg)','FontSize',fontSize);
    if showInput || useSeparatePlot
        xlabel(AxHandle(iUav_).subPlot5,'t (s)','FontSize',fontSize);
    end
    grid(AxHandle(iUav_).subPlot5,'on');
    
    hold(AxHandle(iUav_).subPlot5,'off');
    
    %% 2nd row - phi
    
    % Get figure and axis
    if useSeparatePlot
        set(figHandle(iUav), 'CurrentAxes', AxHandle(iUav).subPlot6);
        iUav_ = iUav;
    else
        set(figHandle(1), 'CurrentAxes', AxHandle(1).subPlot6);
        iUav_ = 1;
    end
    
    hold(AxHandle(iUav_).subPlot6,'on');
    
    % Plot signal
    if showSignal
        plot(AxHandle(iUav_).subPlot6,...
            Data(iUav).t, Data(iUav).phi*RAD_2_DEG, lineSpecSignal{:}, ...
            'Color', colorListSignal(iUav,:));
    end
    
    % Plot tolerance curves
    if showTolerance
        plot(AxHandle(iUav_).subPlot6,...
            Data(iUav).t, phiLower{iUav}*RAD_2_DEG, lineSpecTolerance{:}, ...
            'Color', colorListTolerance(iUav,:));
        plot(AxHandle(iUav_).subPlot6,...
            Data(iUav).t, phiUpper{iUav}*RAD_2_DEG, lineSpecTolerance{:}, ...
            'Color', colorListTolerance(iUav,:));
    end
    
    % Plot reference curves
    if showReference
        plot(AxHandle(iUav_).subPlot6,...
            Data(iUav).t, Data(iUav).phiDesired*RAD_2_DEG, ...
            lineSpecReference{:}, ...
            'Color', colorListReference(iUav,:));
    end
    ylabel(AxHandle(iUav_).subPlot6,'\phi (deg)','FontSize',fontSize);
    if showInput || useSeparatePlot
        xlabel(AxHandle(iUav_).subPlot6,'t (s)','FontSize',fontSize);
    end
    grid(AxHandle(iUav_).subPlot6,'on');
    
    hold(AxHandle(iUav_).subPlot6,'off');
    
    %% 3rd row - Va
    
    if showInput
        % Marker size
        markerSize = 5;
        
        % Get figure and axis
        if useSeparatePlot
            set(figHandle(iUav), 'CurrentAxes', AxHandle(iUav).subPlot7);
            iUav_ = iUav;
        else
            set(figHandle(1), 'CurrentAxes', AxHandle(1).subPlot7);
            iUav_ = 1;
        end
        
        hold(AxHandle(iUav_).subPlot7,'on');
        
        % Marker decimation
        nPoints = numel(Data(iUav).t);
        markerDecimation = 1:round(nPoints/10):nPoints;
        
        % Plot fake signals (for the legend)
        hLine2 = plot(AxHandle(iUav_).subPlot7,...
            Data(iUav).t(markerDecimation), Data(iUav).Va1(markerDecimation),...
            lineSpecSignal{:}, ...
            'Marker', 'o', 'MarkerSize', markerSize, 'MarkerEdgeColor', 'black',...
            'MarkerFaceColor', 'black', 'Color', 'none');
        hLine3 = plot(AxHandle(iUav_).subPlot7,...
            Data(iUav).t(markerDecimation), Data(iUav).Va2(markerDecimation),...
            lineSpecSignal{:}, ...
            'Marker', 'pentagram', 'MarkerSize', markerSize, 'MarkerEdgeColor', 'black',...
            'MarkerFaceColor', 'black', 'Color', 'none');
        hLine4 = plot(AxHandle(iUav_).subPlot7,...
            Data(iUav).t(markerDecimation), Data(iUav).Va3(markerDecimation),...
            lineSpecSignal{:}, ...
            'Marker', 'diamond', 'MarkerSize', markerSize, 'MarkerEdgeColor', 'black',...
            'MarkerFaceColor', 'black', 'Color', 'none');
        hLine5 = plot(AxHandle(iUav_).subPlot7,...
            Data(iUav).t(markerDecimation), Data(iUav).Va4(markerDecimation),...
            lineSpecSignal{:}, ...
            'Marker', 'square', 'MarkerSize', markerSize, 'MarkerEdgeColor', 'black',...
            'MarkerFaceColor', 'black', 'Color', 'none');
        
        % Plot real signals
        plot(AxHandle(iUav_).subPlot7,...
            Data(iUav).t, Data(iUav).Va1,...
            lineSpecSignal{:}, ...
            'Color', 1 + (colorListSignal(iUav,:)-1)*1.0);
        plot(AxHandle(iUav_).subPlot7,...
            Data(iUav).t(markerDecimation), Data(iUav).Va1(markerDecimation),...
            lineSpecSignal{:}, ...
            'Marker', 'o', 'MarkerSize', markerSize, ...
            'MarkerEdgeColor', 1 + (colorListSignal(iUav,:)-1)*1.0, ...
            'MarkerFaceColor', 1 + (colorListSignal(iUav,:)-1)*1.0, ...
            'Color', 'None');
        
        plot(AxHandle(iUav_).subPlot7,...
            Data(iUav).t, Data(iUav).Va2,...
            lineSpecSignal{:}, ...
            'Color', 1 + (colorListSignal(iUav,:)-1)*0.8);
        plot(AxHandle(iUav_).subPlot7,...
            Data(iUav).t(markerDecimation), Data(iUav).Va2(markerDecimation),...
            lineSpecSignal{:}, ...
            'Marker', 'pentagram', 'MarkerSize', markerSize, ...
            'MarkerEdgeColor', 1 + (colorListSignal(iUav,:)-1)*0.9, ...
            'MarkerFaceColor', 1 + (colorListSignal(iUav,:)-1)*0.9, ...
            'Color', 'None');
        
        plot(AxHandle(iUav_).subPlot7,...
            Data(iUav).t, Data(iUav).Va3,...
            lineSpecSignal{:}, ...
            'Color', 1 + (colorListSignal(iUav,:)-1)*0.6);
        plot(AxHandle(iUav_).subPlot7,...
            Data(iUav).t(markerDecimation), Data(iUav).Va3(markerDecimation),...
            lineSpecSignal{:}, ...
            'Marker', 'diamond', 'MarkerSize', markerSize, ...
            'MarkerEdgeColor', 1 + (colorListSignal(iUav,:)-1)*0.8, ...
            'MarkerFaceColor', 1 + (colorListSignal(iUav,:)-1)*0.8, ...
            'Color', 'None');
        
        plot(AxHandle(iUav_).subPlot7,...
            Data(iUav).t, Data(iUav).Va4,...
            lineSpecSignal{:}, ...
            'Color', 1 + (colorListSignal(iUav,:)-1)*0.4);
        plot(AxHandle(iUav_).subPlot7,...
            Data(iUav).t(markerDecimation), Data(iUav).Va4(markerDecimation),...
            lineSpecSignal{:}, ...
            'Marker', 'square', 'MarkerSize', markerSize, ...
            'MarkerEdgeColor', 1 + (colorListSignal(iUav,:)-1)*0.7, ...
            'MarkerFaceColor', 1 + (colorListSignal(iUav,:)-1)*0.7, ...
            'Color', 'None');
        
        % Move fake plots to bottom
        uistack(hLine2, 'bottom');
        uistack(hLine3, 'bottom');
        uistack(hLine4, 'bottom');
        uistack(hLine5, 'bottom');
        
        % Adjust Y tick marks
        if abs(AxHandle(iUav_).subPlot7.YTick(end) - AxHandle(iUav_).subPlot7.YTick(1)) <= 0.0001
            tmpLabel = round(AxHandle(iUav_).subPlot7.YTick,4);
            tmpLabel = cellfun(@num2str,num2cell(tmpLabel(:)),'uniformoutput',false);
            AxHandle(iUav_).subPlot7.YTickLabel = tmpLabel;
        end
        
        % Labels and grid
        ylabel(AxHandle(iUav_).subPlot7,'V_a (V)','FontSize',fontSize);
        if useSeparatePlot
            xlabel(AxHandle(iUav_).subPlot7,'t (s)','FontSize',fontSize);
        end
        grid(AxHandle(iUav_).subPlot7,'on');
        
        % Legend
        if useSeparatePlot
            hLegend2 = legend(AxHandle(iUav_).subPlot7,[hLine2,hLine3,hLine4,hLine5], ...
                {'V_a_1','V_a_2','V_a_3','V_a_4'});
            hLegend2.FontSize = fontSizeLegend;
            hLegend2.Location = 'best';
        elseif ~useSeparatePlot && iUav == nUavs
            hLegend2 = legend(AxHandle(iUav_).subPlot7,[hLine2,hLine3,hLine4,hLine5], ...
                {'V_a_1','V_a_2','V_a_3','V_a_4'});
            hLegend2.FontSize = fontSizeLegend;
            hLegend2.Location = 'best';
        end
        
        hold(AxHandle(iUav_).subPlot7,'off');
    end
    
    %% Final legend for the 6 first plots
    
    if iUav == nUavs && ~useSeparatePlot
        % Create legend
        if numel(legendStr1) > 12
            hLegend1 = legend(AxHandle(iUav_).subPlot5,[hLine1(1:12)],{legendStr1{1:12}});
            hLegend1.Orientation = 'horizontal';
            hLegend1.FontSize = fontSizeLegend;
            % Programatically move the Legend
            newPosition = [0.5,0.034,0.0,0.0];
            newUnits = 'normalized';
            set(hLegend1, 'Position', newPosition, 'Units', newUnits);
            
            hLegend2 = legend(AxHandle(iUav_).subPlot4,[hLine1(13:end)],{legendStr1{13:end}});
            hLegend2.Orientation = 'horizontal';
            hLegend2.FontSize = fontSizeLegend;
            % Programatically move the Legend
            newPosition = [0.5,0.013,0.0,0.0];
            newUnits = 'normalized';
            set(hLegend2, 'Position', newPosition, 'Units', newUnits);
        else
            hLegend1 = legend(AxHandle(iUav_).subPlot5,[hLine1(:)],legendStr1);
            hLegend1.Orientation = 'horizontal';
            hLegend1.FontSize = fontSizeLegend;
            % Programatically move the Legend
            newPosition = [0.5,0.027,0.0,0.0];
            newUnits = 'normalized';
            set(hLegend1, 'Position', newPosition, 'Units', newUnits);
        end
        if showInput
            uistack(hLine2, 'bottom');
        end
    end
    
    if useSeparatePlot
        set(figHandle(iUav), 'CurrentAxes', AxHandle(iUav).subPlot6);
        iUav_ = iUav;
    else
        set(figHandle(1), 'CurrentAxes', AxHandle(1).subPlot6);
        iUav_ = 1;
    end
    if showSignal && showTolerance && showReference
        hLegend1 = legend(AxHandle(iUav_).subPlot6,...
            [hLine55(1),hLine56(1),hLine57(1)],...
            {'Signal','Tolerance','Reference'});
    elseif ~showSignal && showTolerance && showReference
        hLegend1 = legend(AxHandle(iUav_).subPlot6,...
            [hLine56(1),hLine57(1)],...
            {'Tolerance','Reference'});
    elseif showSignal && showTolerance && ~showReference
        hLegend1 = legend(AxHandle(iUav_).subPlot6,...
            [hLine55(1),hLine56(1)],...
            {'Signal','Tolerance'});
    elseif showSignal && ~showTolerance && showReference
        hLegend1 = legend(AxHandle(iUav_).subPlot6,...
            [hLine55(1),hLine57(1)],...
            {'Signal','Reference'});
    elseif ~showSignal && showTolerance && ~showReference
        hLegend1 = legend(AxHandle(iUav_).subPlot6,...
            [hLine56(1)],...
            {'Tolerance'});
    elseif showSignal && ~showTolerance && ~showReference
        hLegend1 = legend(AxHandle(iUav_).subPlot6,...
            [hLine55(1)],...
            {'Signal'});
    elseif ~showSignal && ~showTolerance && showReference
        hLegend1 = legend(AxHandle(iUav_).subPlot6,...
            [hLine57(1)],...
            {'Reference'});
    end
    hLegend1.Orientation = 'horizontal';
    hLegend1.FontSize = fontSizeLegend;
    newPosition = [0.5,0.97,0.0,0.0];
    newUnits = 'normalized';
    set(hLegend1, 'Position', newPosition, 'Units', newUnits);
    
end

%% Plotting errors

legendStr1 = cell(1,nUavs);
for iUav = 1:nUavs
    
    %% 1st row - x
    
    % Get figure and axis
    if useSeparatePlot
        set(figHandleError(iUav), 'CurrentAxes', AxHandleError(iUav).subPlot1);
        iUav_ = iUav;
    else
        set(figHandleError(1), 'CurrentAxes', AxHandleError(1).subPlot1);
        iUav_ = 1;
    end
    
    hold(AxHandleError(iUav_).subPlot1,'on');
    
    hLine55Error(iUav) = plot(AxHandleError(iUav_).subPlot1,...
        Data(iUav).t, xError{iUav}, lineSpecSignal{:}, ...
        'Color', 'black');
    hLine1Error(iUav) = plot(AxHandleError(iUav_).subPlot1,...
        Data(iUav).t, xError{iUav}, lineSpecSignal{:}, ...
        'Color', colorListSignal(iUav,:));
    
    % Labels and grid
    ylabel(AxHandleError(iUav_).subPlot1,'Error x (m)','FontSize', fontSize);
    xlabel(AxHandleError(iUav_).subPlot1,'t (s)','FontSize', fontSize);
    grid(AxHandleError(iUav_).subPlot1,'on');
    
    % Legend string
    if isempty(customLegend)
        legendStr1Error{iUav} = ['UAV ',num2str(iUav)];
    else
        legendStr1Error{iUav} = customLegend{iUav};
    end
    
    uistack([hLine55Error(iUav)], 'bottom');
    
    hold(AxHandleError(iUav_).subPlot1,'off');
    
    %% 1st row - y
    
    % Get figure and axis
    if useSeparatePlot
        set(figHandleError(iUav), 'CurrentAxes', AxHandleError(iUav).subPlot2);
        iUav_ = iUav;
    else
        set(figHandleError(1), 'CurrentAxes', AxHandleError(1).subPlot2);
        iUav_ = 1;
    end
    
    hold(AxHandleError(iUav_).subPlot2,'on');
    
    % Plot error signal
    plot(AxHandleError(iUav_).subPlot2,...
        Data(iUav).t, yError{iUav}, lineSpecSignal{:}, ...
        'Color', colorListSignal(iUav,:));
    
    ylabel(AxHandleError(iUav_).subPlot2,'Error y (m)','FontSize', fontSize);
    xlabel(AxHandleError(iUav_).subPlot2,'t (s)','FontSize', fontSize);
    grid(AxHandleError(iUav_).subPlot2,'on');
    
    hold(AxHandleError(iUav_).subPlot2,'off');
    
    %% 1st row - z
    
    % Get figure and axis
    if useSeparatePlot
        set(figHandleError(iUav), 'CurrentAxes', AxHandleError(iUav).subPlot3);
        iUav_ = iUav;
    else
        set(figHandleError(1), 'CurrentAxes', AxHandleError(1).subPlot3);
        iUav_ = 1;
    end
    
    hold(AxHandleError(iUav_).subPlot3,'on');
    
    % Plot error signal
    plot(AxHandleError(iUav_).subPlot3,...
        Data(iUav).t, zError{iUav}, lineSpecSignal{:}, ...
        'Color', colorListSignal(iUav,:));
    
    ylabel(AxHandleError(iUav_).subPlot3,'Error z (m)','FontSize', fontSize);
    xlabel(AxHandleError(iUav_).subPlot3,'t (s)','FontSize', fontSize);
    grid(AxHandleError(iUav_).subPlot3,'on');
    
    hold(AxHandleError(iUav_).subPlot3,'off');
    
    %% 2nd row - psi
    
    % Get figure and axis
    if useSeparatePlot
        set(figHandleError(iUav), 'CurrentAxes', AxHandleError(iUav).subPlot4);
        iUav_ = iUav;
    else
        set(figHandleError(1), 'CurrentAxes', AxHandleError(1).subPlot4);
        iUav_ = 1;
    end
    
    hold(AxHandleError(iUav_).subPlot4,'on');
    
    % Plot error signal
    plot(AxHandleError(iUav_).subPlot4,...
        Data(iUav).t, psiError{iUav}*RAD_2_DEG, lineSpecSignal{:}, ...
        'Color', colorListSignal(iUav,:));
    
    ylabel(AxHandleError(iUav_).subPlot4,'Error \psi (deg)','FontSize', fontSize);
    if useSeparatePlot
        xlabel(AxHandleError(iUav_).subPlot4,'t (s)','FontSize', fontSize);
    end
    grid(AxHandleError(iUav_).subPlot4,'on');
    
    hold(AxHandleError(iUav_).subPlot4,'off');
    
    %% 2nd row - theta
    
    % Get figure and axis
    if useSeparatePlot
        set(figHandleError(iUav), 'CurrentAxes', AxHandleError(iUav).subPlot5);
        iUav_ = iUav;
    else
        set(figHandleError(1), 'CurrentAxes', AxHandleError(1).subPlot5);
        iUav_ = 1;
    end
    
    hold(AxHandleError(iUav_).subPlot5,'on');
    
    % Plot error signal
    plot(AxHandleError(iUav_).subPlot5,...
        Data(iUav).t, thetaError{iUav}*RAD_2_DEG, lineSpecSignal{:}, ...
        'Color', colorListSignal(iUav,:));
    
    ylabel(AxHandleError(iUav_).subPlot5,'Error \theta (deg)','FontSize',fontSize);
    if useSeparatePlot
        xlabel(AxHandleError(iUav_).subPlot5,'t (s)','FontSize',fontSize);
    end
    grid(AxHandleError(iUav_).subPlot5,'on');
    
    hold(AxHandleError(iUav_).subPlot5,'off');
    
    %% 2nd row - phi
    
    % Get figure and axis
    if useSeparatePlot
        set(figHandleError(iUav), 'CurrentAxes', AxHandleError(iUav).subPlot6);
        iUav_ = iUav;
    else
        set(figHandleError(1), 'CurrentAxes', AxHandleError(1).subPlot6);
        iUav_ = 1;
    end
    
    hold(AxHandleError(iUav_).subPlot6,'on');
    
    % Plot error signal
    plot(AxHandleError(iUav_).subPlot6,...
        Data(iUav).t, phiError{iUav}*RAD_2_DEG, lineSpecSignal{:}, ...
        'Color', colorListSignal(iUav,:));
    
    ylabel(AxHandleError(iUav_).subPlot6,'Error \phi (deg)','FontSize',fontSize);
    if useSeparatePlot
        xlabel(AxHandleError(iUav_).subPlot6,'t (s)','FontSize',fontSize);
    end
    grid(AxHandleError(iUav_).subPlot6,'on');
    
    hold(AxHandleError(iUav_).subPlot6,'off');
    
    %% Final legends
    
    if iUav == nUavs && ~useSeparatePlot
        % Create legend
        if numel(legendStr1Error) > 12
            hLegend1Error = legend(AxHandleError(iUav_).subPlot5,[hLine1Error(1:12)],{legendStr1Error{1:12}});
            hLegend1Error.Orientation = 'horizontal';
            hLegend1Error.FontSize = fontSizeLegend;
            % Programatically move the Legend
            newPosition = [0.5,0.034,0.0,0.0];
            newUnits = 'normalized';
            set(hLegend1Error, 'Position', newPosition, 'Units', newUnits);
            
            hLegend2Error = legend(AxHandleError(iUav_).subPlot4,[hLine1Error(13:end)],{legendStr1Error{13:end}});
            hLegend2Error.Orientation = 'horizontal';
            hLegend2Error.FontSize = fontSizeLegend;
            % Programatically move the Legend
            newPosition = [0.5,0.013,0.0,0.0];
            newUnits = 'normalized';
            set(hLegend2Error, 'Position', newPosition, 'Units', newUnits);
        else
            hLegend1Error = legend(AxHandleError(iUav_).subPlot5,[hLine1Error(:)],legendStr1Error);
            hLegend1Error.Orientation = 'horizontal';
            hLegend1Error.FontSize = fontSizeLegend;
            % Programatically move the Legend
            newPosition = [0.5,0.027,0.0,0.0];
            newUnits = 'normalized';
            set(hLegend1Error, 'Position', newPosition, 'Units', newUnits);
        end
    end
    
end

%% Export plots

if exportEn
    pngEn = true;
    if useSeparatePlot
        for iUav = 1:nUavs
            export_graph(figHandle(iUav), [fileName,'_',num2str(iUav)], folderName, pngEn);
            if plotError
                export_graph(figHandleError(iUav), [fileName,'_Error_',num2str(iUav)], folderName, pngEn);
            end
        end
    else
        export_graph(figHandle, fileName, folderName, pngEn);
        if plotError
            export_graph(figHandleError, [fileName,'_Error'], folderName, pngEn);
        end
    end
end

end
