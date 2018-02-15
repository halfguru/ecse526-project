function [nRecurrentConnections, directedGraph, recurrentConnection] = ...
    find_recurrent_connection(Population)

nIndividuals = numel([Population(:).fitness]);
nRecurrentConnections = zeros(1,nIndividuals);
directedGraph = cell(1,nIndividuals);
recurrentConnection = cell(1,nIndividuals);
for iIndividual = 1:nIndividuals
    
    % Individual
    Individual = Population(iIndividual);
    
    % Extract node ID and node type
    nodeId   = Individual.nodeGenes(1,:);
    nodeType = Individual.nodeGenes(2,:);
    
    % Sort nodeId in ascending order
    [nodeId, sortIndex] = sort(nodeId);
    
    % Sort nodeType according to sortIndex
    nodeType = nodeType(sortIndex);
    
    % We put output layer at the end for visualisation purposes (when using
    % the view function)
    nodeId    = [nodeId(nodeType ~= 2),   nodeId(nodeType == 2)];
    nodeType  = [nodeType(nodeType ~= 2), nodeType(nodeType == 2)];
    
    % Extract connections info
    connId      = Individual.connectionGenes(1,:);
    connFrom    = Individual.connectionGenes(2,:);
    connTo      = Individual.connectionGenes(3,:);
    connWeight  = Individual.connectionGenes(4,:);
    connEnabled = Individual.connectionGenes(5,:);
    
    % We remove disabled connections
    indexEnabled = connEnabled == 1;
    connId       = connId(indexEnabled);
    connFrom     = connFrom(indexEnabled);
    connTo       = connTo(indexEnabled);
    connWeight   = connWeight(indexEnabled);
    connEnabled  = connEnabled(indexEnabled);
    
    % Sort connId in ascending order
    [connId, sortIndex] = sort(connId);
    
    % Sort connFrom, connTo, connWeight and connEnabled according to sortIndex
    connFrom    = connFrom(sortIndex);
    connTo      = connTo(sortIndex);
    connWeight  = connWeight(sortIndex);
    connEnabled = connEnabled(sortIndex);
    
    % Inputs' node ID
    inputId = nodeId(nodeType == 1);
    
    % Outputs' node ID
    outputId = nodeId(nodeType == 2);
    
    % Hidden nodes' node ID
    hiddenId = nodeId(nodeType == 3);
    
    % Biases' node ID
    biasId = nodeId(nodeType == 4);
    
    % Number of nodes
    nNodes = numel(nodeId);
    
    % Number of connections (including biases)
    nConns = numel(connId);
    
    % Number of biases
    nBiases = sum(nodeType == 4);
    
    % Number of layers
    numLayers = nNodes - nBiases;
    
    % Array with first column containing the layer ID and second column its
    % corresponding ID
    layerList = [[1:numLayers]',[nodeId(nodeType ~= 4)]'];
    
    % Inputs node's layer ID according to layerList
    inputLayerId = layerList(ismember(layerList(:,2),inputId),1);
    
    % Outputs node's layer ID according to layerList
    outputLayerId = layerList(ismember(layerList(:,2),outputId),1);
    
    % Hidden node's layer ID according to layerList
    hiddenLayerId = layerList(ismember(layerList(:,2),hiddenId),1);
    
    % First, remove "from connections" of biases from connFrom
    connIdNoBias   = connId;
    connFromNoBias = connFrom;
    connToNoBias   = connTo;
    connWeightNoBias = connWeight;
    for iBias = 1:nBiases
        indexOfBias     = connFromNoBias ~= biasId(iBias);
        connIdNoBias    = connIdNoBias(indexOfBias);
        connFromNoBias  = connFromNoBias(indexOfBias);
        connToNoBias    = connToNoBias(indexOfBias);
        connWeightNoBias = connWeightNoBias(indexOfBias);
    end
    
    % Next, create start and target layers arrays
    % Important: This will fail if no input are connected!
    startLayers  = [];
    targetLayers = [];
    layerWeight = [];
    for ii = 1:numel(connFromNoBias)
        startLayers  = [startLayers,...
            layerList(layerList(:,2) == connFromNoBias(ii),1)];
        targetLayers = [targetLayers,...
            layerList(layerList(:,2) == connToNoBias(ii),1)];
        layerWeight = [layerWeight,...
            connWeightNoBias(...
            connFromNoBias == connFromNoBias(ii) & ...
            connToNoBias == connToNoBias(ii) ...
            )...
            ];
    end
    tmpArray = [...
        startLayers;
        targetLayers]';
    [tmpArray,iUnique] = unique(tmpArray,'rows');
    startLayers = tmpArray(:,1)';
    targetLayers = tmpArray(:,2)';
    layerWeight = layerWeight(iUnique);
    
    % Build directed graph object
    % directedGraph = digraph(startLayers, targetLayers);
    directedGraph{iIndividual} = digraph(startLayers, targetLayers, layerWeight);
    
    % Find recurrent connections
    recurrentConnection{iIndividual} = ...
        dfsearch(directedGraph{iIndividual}, 1, 'edgetodiscovered',...
        'Restart', true);
    
    % Find connections that go from one output to another output
    % excluding output that goes into itself
    for ii = 1:numel(connFromNoBias)
        condition1 = nodeType(nodeId == connFromNoBias(ii)) == 2;
        condition2 = nodeType(nodeId == connToNoBias(ii)) == 2;
        condition3 = nodeId(nodeId == connFromNoBias(ii)) ~= nodeId(nodeId == connToNoBias(ii));
        if condition1 && condition2 && condition3
            recurrentConnection{iIndividual} = ...
                [recurrentConnection{iIndividual};
                [layerList(layerList(:,2) == connFromNoBias(ii),1),...
                layerList(layerList(:,2) == connToNoBias(ii),1)]...
                ];
        end
    end
    
    % Number of recurrent connections
    if ~isempty(recurrentConnection{iIndividual})
        nRecurrentConnections(iIndividual) = numel(recurrentConnection{iIndividual}(:,1));
    else
        nRecurrentConnections(iIndividual) = 0;
    end
    
%     if nRecurrentConnections(iIndividual) ~= 0
%         pause
%     end
    
end

end
