function figHandle = plot_neural_network(matFileFolder, ...
    matFilePrefix, folderName, fileName, nGenToPlot, ...
    figureVisibility, exportEn, varargin)

%% Import data

if nargin == 7
    [populationArray, ~, ~, ~, nGenerations] = ...
        import_data(matFileFolder, matFilePrefix);
elseif nargin == 8 || nargin == 9
    populationArray = {varargin{1}};
    nGenerations = 1;
end

%% Create directed graph

nGenToPlot = min([min([nGenToPlot,nGenerations]),12]);
genListTmp = 1:nGenerations;
if nGenToPlot == 1
    genToPlot = nGenerations;
elseif nGenToPlot == 2
    genToPlot = [1,nGenerations];
elseif nGenToPlot > 2
    genToPlot = ...
        [1,genListTmp(1+randperm(nGenerations-2,nGenToPlot-2)),nGenerations];
end
genToPlot = sort(genToPlot);

directedGraph         = cell(1,nGenToPlot);
recurrentConnection   = cell(1,nGenToPlot);
nRecurrentConnections = cell(1,nGenToPlot);
for iGen = 1:nGenToPlot    
    
    % Assign variables
    iGen2 = genToPlot(iGen);
    Population = populationArray{iGen2};
    
    % Directed graph
    [nRecurrentConnections{iGen}, directedGraph{iGen}, recurrentConnection{iGen}] = ...
        find_recurrent_connection(Population);
    
end

%% Initialize figure

resolutionX = 1000;
resolutionY = 940;
if figureVisibility
    if nargin == 9
        figHandle = varargin{2};
    else
        figHandle = create_figure(resolutionX, resolutionY, 'on');
    end
else
    if nargin == 9
        figHandle = varargin{2};
    else
        figHandle = create_figure(resolutionX, resolutionY, 'off');
    end
end

%% Create axes

if nGenToPlot == 1
    subplotConfig = {1,1};
elseif nGenToPlot == 2
    subplotConfig = {1,2};
elseif nGenToPlot == 3
    subplotConfig = {1,3};
elseif nGenToPlot > 3 && nGenToPlot <= 6
    subplotConfig = {2,3};
elseif nGenToPlot > 6 && nGenToPlot <= 9
    subplotConfig = {3,3};
elseif nGenToPlot > 9 && nGenToPlot <= 12
    subplotConfig = {3,4};
end

for iGen = 1:nGenToPlot
    set(0, 'currentfigure', figHandle);
    AxHandle(iGen).subPlot = subplot_tight(subplotConfig{:},iGen,0.00);
    color = get(figHandle,'Color');
    set(AxHandle(iGen).subPlot,'XColor',color,'YColor',color,'TickDir','out')
end

%% Plot

fontSize = 9;
weightMin = 1;
weightMax = 2;
if nGenToPlot <= 2
    lineWidthMin = 0.1;
    lineWidthMax = 2.5;
    markerSize = 4;
    arrowSize = 7;
elseif nGenToPlot <= 6
    lineWidthMin = 0.1/2;
    lineWidthMax = 2.5/2;
    markerSize = 3;
    arrowSize = 4;
elseif nGenToPlot <= 9
    lineWidthMin = 0.1/3;
    lineWidthMax = 2.5/3;
    markerSize = 2;
    arrowSize = 3;
elseif nGenToPlot <= 12
    lineWidthMin = 0.1/4;
    lineWidthMax = 2.5/4;
    markerSize = 1;
    arrowSize = 2;
end

nInputs = 25;
nOutputs = 4;

for iGen = 1:nGenToPlot
 
    % Find best individual of this generation
    iGen2 = genToPlot(iGen);
    Population = populationArray{iGen2};
    tmpFitness = [Population(:).fitness];
    [~, iTopIndividual] = max(tmpFitness);
    
    % Directed graph
    directedGraphTmp = directedGraph{iGen}{iTopIndividual};
    set(figHandle, 'CurrentAxes', AxHandle(iGen).subPlot);
    
    hold(AxHandle(iGen).subPlot,'on');
    
    % Plotting
    pHandle = plot(AxHandle(iGen).subPlot, directedGraphTmp);
    
    % Node and arrow sizes
    pHandle.MarkerSize = markerSize;
    pHandle.ArrowSize  = arrowSize;
    
    % Node label
    if nGenToPlot >= 6
        pHandle.NodeLabel = {};
    end
    
    % Linewidth
    weightTmp = abs(directedGraphTmp.Edges.Weight);
    weightTmp = mapminmax(weightTmp',weightMin,weightMax);
    lineWidthTmp = mapminmax(weightTmp,lineWidthMin,lineWidthMax);
    directedGraphTmp.Edges.LWidths = lineWidthTmp';
    pHandle.LineWidth = directedGraphTmp.Edges.LWidths;
    
    % Remove axes tick lines
    AxHandle(iGen).subPlot.XTick = [];
    AxHandle(iGen).subPlot.XTickLabel = '';
    AxHandle(iGen).subPlot.YTick = [];
    AxHandle(iGen).subPlot.YTickLabel = '';
    
    % Set up layout
    nNodes = size(directedGraphTmp.Nodes,1);
    layout(pHandle,'layered','Sources',1:nInputs,'Sinks',(nNodes-nOutputs+1):(nNodes));
    % layout(pHandle,'layered','Sinks',(nNodes-nOutputs+1):(nNodes));
    
    % Color nodes
    highlight(pHandle,1:nInputs,'NodeColor','r');
    highlight(pHandle,(nNodes-nOutputs+1):(nNodes),'NodeColor','g');
    
    % Highlight recurrent connections
    if nRecurrentConnections{iGen}(iTopIndividual) ~= 0
        ss = recurrentConnection{iGen}{iTopIndividual}(:,1);
        tt = recurrentConnection{iGen}{iTopIndividual}(:,2);        
        highlight(pHandle,ss,tt,'EdgeColor','m');
    end
    
    % Axes position
    axesPos = AxHandle(iGen).subPlot.Position;
    
    % Title
    if nGenToPlot > 1
        textHandle = text(0.5,0.5,['Generation ',num2str(genToPlot(iGen))]);
        textHandle.FontSize = fontSize;
        textHandle.FontWeight = 'bold';
        textHandle.Units = 'normalized';
        textHandle.HorizontalAlignment = 'center';
        textHandle.VerticalAlignment = 'baseline';
        % xPos = axesPos(1) + 0.5*axesPos(3);
        % yPos = axesPos(2) + axesPos(4);
        xPos = 0.5;
        yPos = 0.95;
        textHandle.Position = [xPos,yPos,0];
    end
    
    axis(AxHandle(iGen).subPlot,'tight');
    
    % titleHandle = title(AxHandle(iGen).subPlot,['Generation ',num2str(genToPlot(iGen))],...
    %     'FontSize',fontSize);
    % axesPos = plotboxpos(AxHandle(iGen).subPlot);
    % axesPos = AxHandle(iGen).subPlot.Position;
    % oldUnits = titleHandle.Units;
    % titleHandle.Units = AxHandle(iGen).subPlot.Units;
    % xTitle = titleHandle.Position(1);
    % yTitle = (AxHandle(iGen).subPlot.Position(2) + AxHandle(iGen).subPlot.Position(4));
    % yTitle = (axesPos(2) + axesPos(4));
    % zTitle = titleHandle.Position(3);
    % titleHandle.Position = [xTitle,yTitle,zTitle];
    % titleHandle.Units = oldUnits;
    
    hold(AxHandle(iGen).subPlot,'off');
    
end

%% Export plots

if exportEn
    pngEn = true;
    export_graph(figHandle, fileName, folderName, pngEn);
end

end
