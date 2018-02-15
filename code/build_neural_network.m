%% build_neural_network() (Main)
function [neuralNet, nDelays, varargout] = build_neural_network(Individual)
% This function generate a neural network from each individual's codeGenes and
% connectionGenes
% References:   https://goo.gl/FcHGOz
%               https://goo.gl/rdtEzo
%               https://goo.gl/odfpyB (genFunction)
%
% Note: Neural network toolbox uses layers for nodes since a layer can contains
% multiple neurons. However, NEAT generates only nodes with single neuron. Hence,
% in this function, layers correspond to nodes.
%
% Also, NEAT's input nodes (nodeType = 1) are input ports in Neural network
% toolbox's terminology. The Neural network toolbox's input layers are the
% layers directly connected to the input ports. To simplify things, we'll define
% NEAT's input nodes as input layers with weight = 1 and linear activation
% function. NEAT's output nodes corresponds to the toolbox's last layer before
% connection to output.

%% Preprocessing

% Extract node ID and node type
nodeId   = Individual.nodeGenes(1,:);
nodeType = Individual.nodeGenes(2,:);

% Sort nodeId in ascending order
[nodeId, sortIndex] = sort(nodeId);

% Sort nodeType according to sortIndex
nodeType = nodeType(sortIndex);

% We put output layer at the end for visualisation purposes (when using the view
% function)
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

% Biases' node ID. NEAT uses only one bias for all the networks and choose
% whether to connect it to other nodes or not. By assigning a different weight
% to every connection from the bias to a node, it is as if these nodes had their
% own bias.
biasId = nodeId(nodeType == 4);

%% Network architecture

% Initialize network object
neuralNet = network();

% Number of nodes
nNodes = numel(nodeId);

% Number of connections (including biases)
nConns = numel(connId);

% Number of biases. NEAT uses only one bias for all the networks and choose
% whether to connect it to other nodes or not. By assigning a different weight
% to every connection from the bias to a node, it is as if these nodes had their
% own bias.
nBiases = sum(nodeType == 4);

% Number of inputs
neuralNet.numInputs = sum(nodeType == 1);

% Number of layers
% Neural network toolbox does not consider biases as input nodes
neuralNet.numLayers = nNodes - nBiases;

% Note: Neural network toolbox uses layers for nodes since a layer can contains
% multiple neurons. However, NEAT generates only nodes with single neuron. So,
% layers correspond to nodes here.

% Array with first column containing the layer ID and second column its
% corresponding ID. We need to do this since neural network toolbox does not
% consider biases as input nodes
layerList = [[1:neuralNet.numLayers]',[nodeId(nodeType ~= 4)]'];

% Inputs node's layer ID according to layerList
inputLayerId = layerList(ismember(layerList(:,2),inputId),1);

% Outputs node's layer ID according to layerList
outputLayerId = layerList(ismember(layerList(:,2),outputId),1);

% Hidden node's layer ID according to layerList
hiddenLayerId = layerList(ismember(layerList(:,2),hiddenId),1);

%% Bias connections

% We connect the biases to their layers
for jj = biasId
    % Find nodes connected to this bias. connFrom == jj gives a logical indexing
    % vector of the node that receive a connection from the bias.
    % connTo(connFrom == jj) gives the actual node IDs connected to the bias.
    iNodeConnectedToBias = connTo(connFrom == jj);
    % Find layers connected to this bias
    iLayerConnectedToBias = [];
    for kk = 1:numel(iNodeConnectedToBias)
        iLayerConnectedToBias = ...
            [iLayerConnectedToBias, ...
            layerList(layerList(:,2) == iNodeConnectedToBias(kk),1)];
    end
    % Set bias for the above layers (if any)
    if ~isempty(iLayerConnectedToBias)
        for ii = iLayerConnectedToBias
            % neuralNet.biasConnect(TO ith layer)
            neuralNet.biasConnect(ii) = 1;
        end
    end
end

%% Connect inputs, outputs and hidden layers

% Input counter for the following loop
inputCounter = 1;

% Cycle through layers (nodes)
for iNode = 1:nNodes
    
    %% Input connections
    
    % We specify which layers are input layers. Note that this loop assumes that
    % the input nodes' IDs come before the hidden and output layers.
    if nodeType(iNode) == 1 % If is input node
        % Find layer connected to this input
        iLayerConnectedToInput = layerList(layerList(:,2) == nodeId(iNode),1);
        % Set input for the above layers (if any)
        if ~isempty(iLayerConnectedToInput)
            % neuralNet.inputConnect(TO ith layer, FROM jth input)
            neuralNet.inputConnect(iLayerConnectedToInput,inputCounter) = 1;
        end
        inputCounter = inputCounter + 1;
    end
    
    %% Layer connection from input layer and hidden layers
    
    % We connect input layers to hidden layers as well as hidden layers to other
    % hidden layers
    if nodeType(iNode) == 1 || nodeType(iNode) == 2 || nodeType(iNode) == 3 % If is input or hidden node *********
        % Find nodes that receive connection FROM this node
        iNodeConnectedFromCurrLayer = connTo(connFrom == nodeId(iNode));
        % Find layers that receive connection FROM this node
        iLayerConnectedFromCurrLayer = [];
        for kk = 1:numel(iNodeConnectedFromCurrLayer)
            iLayerConnectedFromCurrLayer = ...
                [iLayerConnectedFromCurrLayer, ...
                layerList(layerList(:,2) == iNodeConnectedFromCurrLayer(kk),1)];
        end
        % Set FROM current layer for the above layers (if any)
        if ~isempty(iLayerConnectedFromCurrLayer)
            for ii = iLayerConnectedFromCurrLayer
                % neuralNet.layerConnect(TO ith layer, FROM jth layer)
                iTmp = layerList(layerList(:,2) == nodeId(iNode),1);
                neuralNet.layerConnect(ii,iTmp) = 1;
            end
        end
    end
    
    %% Layer connection to output layers
    
    % We connect hidden layers to the output layers. At this point, the
    % connections FROM other layers (hidden, bias or input) to the output nodes
    % is already done (in the "Layer connection from input layer and hidden
    % layers" section above). What is left to do is to connect the output nodes'
    % output to previous layers to form recurrent connection.
    %     if nodeType(iNode) == 2 % If is output node
    %         % Find nodes that are connected TO this node
    %         iNodeConnectedToCurrLayer = connFrom(connTo == nodeId(iNode));
    %         % Find layers connected TO this node
    %         iLayerConnectedToCurrLayer = [];
    %         for kk = 1:numel(iNodeConnectedToCurrLayer)
    %             % If iNodeConnectedToCurrLayer(kk) is not a bias
    %             if ~any(iNodeConnectedToCurrLayer(kk) == biasId)
    %                 iLayerConnectedToCurrLayer = ...
    %                     [iLayerConnectedToCurrLayer, ...
    %                     layerList(layerList(:,2) == iNodeConnectedToCurrLayer(kk),1)];
    %             end
    %         end
    %         % Set FROM current layer for the above layers (if any)
    %         if ~isempty(iLayerConnectedToCurrLayer)
    %             for ii = iLayerConnectedToCurrLayer
    %                 % neuralNet.layerConnect(TO ith layer, FROM jth layer)
    %                 iTmp = layerList(layerList(:,2) == nodeId(iNode),1);
    %                 neuralNet.layerConnect(iTmp,ii) = 1;
    %             end
    %         end
    %     end
    
    %     % We connect hidden layers to the output layers. At this point, the
    %     % connections FROM other layers (hidden, bias or input) to the output nodes
    %     % is already done (in the "Layer connection from input layer and hidden
    %     % layers" section above). What is left to do is to connect the output nodes'
    %     % output to previous layers to form recurrent connection.
    %     if nodeType(iNode) == 2 % If is output node
    %         % Find nodes that are connected TO this node
    %         iNodeConnectedToCurrLayer = connFrom(connTo == nodeId(iNode));
    %         % Find layers connected TO this node
    %         iLayerConnectedToCurrLayer = [];
    %         for kk = 1:numel(iNodeConnectedToCurrLayer)
    %             % If iNodeConnectedToCurrLayer(kk) is not a bias
    %             if ~any(iNodeConnectedToCurrLayer(kk) == biasId)
    %                 iLayerConnectedToCurrLayer = ...
    %                     [iLayerConnectedToCurrLayer, ...
    %                     layerList(layerList(:,2) == iNodeConnectedToCurrLayer(kk),1)];
    %             end
    %         end
    %         % Set FROM current layer for the above layers (if any)
    %         if ~isempty(iLayerConnectedToCurrLayer)
    %             for ii = iLayerConnectedToCurrLayer
    %                 % neuralNet.layerConnect(TO ith layer, FROM jth layer)
    %                 iTmp = layerList(layerList(:,2) == nodeId(iNode),1);
    %                 neuralNet.layerConnect(iTmp,ii) = 1;
    %             end
    %         end
    %     end
    
    %% Output connections
    
    % We specify which layers are output layers
    if nodeType(iNode) == 2 % If is output node
        % Find layer ID corresponding to current node ID
        iLayerConnectedToOutput = layerList(layerList(:,2) == nodeId(iNode),1);
        % Set FROM for this output
        if ~isempty(iLayerConnectedToOutput)
            % neuralNet.outputConnect(FROM ith layer)
            neuralNet.outputConnect(iLayerConnectedToOutput) = 1;
        end
    end
    
end

%% Set layer connection and delays for recurrent connections
% We use MATLAB's graph theory toolbox (https://goo.gl/O0cbX5) to find
% recurrences in the network

% First, remove "from connections" of biases from connFrom
connFromNoBias = connFrom;
connToNoBias   = connTo;
for iBias = 1:nBiases
    indexOfBias     = connFromNoBias ~= biasId(iBias);
    connFromNoBias  = connFromNoBias(indexOfBias);
    connToNoBias    = connToNoBias(indexOfBias);
end

% Next, create start and target layers arrays
% Important: This will fail if no input are connected!
startLayers  = [];
targetLayers = [];
for ii = 1:numel(connFromNoBias)
    startLayers  = [startLayers,layerList(layerList(:,2) == connFromNoBias(ii),1)];
    targetLayers = [targetLayers,layerList(layerList(:,2) == connToNoBias(ii),1)];
end
tmpArray = [...
    startLayers;
    targetLayers]';
tmpArray = unique(tmpArray,'rows');
startLayers = tmpArray(:,1)';
targetLayers = tmpArray(:,2)';

% Build directed graph object
directedGraph = digraph(startLayers, targetLayers);

% Find recurrent connections
recurrentConnection = dfsearch(directedGraph, 1, 'edgetodiscovered', 'Restart', true);
% plot(directedGraph, 'Layout', 'force');

% Find connections that go from one output to another output
% excluding output that goes into itself
for ii = 1:numel(connFromNoBias)
    condition1 = nodeType(nodeId == connFromNoBias(ii)) == 2;
    condition2 = nodeType(nodeId == connToNoBias(ii)) == 2;
    condition3 = nodeId(nodeId == connFromNoBias(ii)) ~= nodeId(nodeId == connToNoBias(ii));
    if condition1 && condition2 && condition3
        recurrentConnection = ...
            [recurrentConnection;
            [layerList(layerList(:,2) == connFromNoBias(ii),1),...
            layerList(layerList(:,2) == connToNoBias(ii),1)]...
            ];
    end
end

% Finally, assign delays to recurrent connections (layers)
for ii = 1:numel(recurrentConnection(:,1))
    % neuralNet.layerWeights{TO ith layer, FROM jth layer}.delays
    iLayer = recurrentConnection(ii,2);
    jLayer = recurrentConnection(ii,1);
    neuralNet.layerConnect(iLayer,jLayer) = 1; % Is this necessary ?
    neuralNet.layerWeights{iLayer,jLayer}.delays = 1;
end

%% Inputs and hidden layers size and activation function

% Set input sizes (all inputs are size 1)
neuralNet.inputs{:}.size = 1;

% Set layers size, activation function and weight initialization function
% We don't use the weight initialization function but the toolbox requires it
neuralNet.layers{:}.size = 1;
neuralNet.layers{hiddenLayerId}.transferFcn = 'logsig';
neuralNet.layers{outputLayerId}.transferFcn = 'logsig';
neuralNet.layers{inputLayerId}.transferFcn  = 'purelin';
neuralNet.layers{:}.initFcn                 = ''; % initnw

%% Set network functions

% We don't need those but Neural network toolbox requires them in order to
% initialize the network
neuralNet.initFcn = ''; % initlay

neuralNet.performFcn = ''; % mse
neuralNet.trainFcn = ''; % trainlm

neuralNet.divideFcn = ''; % dividerand

% neuralNet.plotFcns = {'plotperform','plottrainstate'};
neuralNet.plotFcns = {};

%% Initialize network

neuralNet = init(neuralNet);

%% Set biases weights

for jj = biasId
    % Find nodes connected to this bias (i.e. nodes that receive connection from
    % the bias)
    iNodeConnectedToBias = connTo(connFrom == jj);
    % Find layers connected to this bias
    iLayerConnectedToBias = [];
    for kk = 1:numel(iNodeConnectedToBias)
        iLayerConnectedToBias = ...
            [iLayerConnectedToBias, ...
            layerList(layerList(:,2) == iNodeConnectedToBias(kk),1)];
    end
    % Set bias
    if ~isempty(iLayerConnectedToBias)
        for kk = 1:numel(iLayerConnectedToBias)
            ii = iLayerConnectedToBias(kk);
            biasWeight = connWeight(connFrom == jj & connTo == iNodeConnectedToBias(kk));
            neuralNet.b{ii} = biasWeight;
        end
    end
end

%% Set layers weights

% Input counter for the following loop
inputCounter = 1;

% Cycle through layers (nodes)
for iNode = 1:nNodes
    
    %% Weight input layer
    
    if nodeType(iNode) == 1 % If is input node
        % Find input's layer Id connected to this input
        iLayerConnectedToInput = layerList(layerList(:,2) == nodeId(iNode),1);
        % Set weight
        if ~isempty(iLayerConnectedToInput)
            % neuralNet.IW{TO ith layer, FROM jth input}
            neuralNet.IW{iLayerConnectedToInput,inputCounter} = 1;
        end
        inputCounter = inputCounter + 1;
    end
    
    %% Connections to hidden and output layers
    
    if nodeType(iNode) == 2 || nodeType(iNode) == 3 % If output or hidden node
        % Find nodes that are connected FROM this node
        iNodeConnectedFromCurrLayer = connFrom(connTo == nodeId(iNode));
        % Find layers connected from this node
        iLayerConnectedFromCurrLayer = [];
        for kk = 1:numel(iNodeConnectedFromCurrLayer)
            iLayerConnectedFromCurrLayer = ...
                [iLayerConnectedFromCurrLayer, ...
                layerList(layerList(:,2) == iNodeConnectedFromCurrLayer(kk),1)];
        end
        % Set FROM current layer for the above layers (if any)
        if ~isempty(iLayerConnectedFromCurrLayer)
            for jj = iLayerConnectedFromCurrLayer
                % neuralNet.LW{TO ith layer, FROM jth layer}
                iConnId = ...
                    (connFrom == layerList(jj,2)) & (connTo == nodeId(iNode));
                ii = layerList(layerList(:,2) == nodeId(iNode),1);
                
                % neuralNet.LW{ii,jj}
                % connWeight(iConnId)
                % if isempty(neuralNet.LW{ii,jj})
                % 	neuralNet.LW{ii,jj} = [];
                % end
                % neuralNet.LW{ii,jj}
                % if isempty(neuralNet.LW{ii,jj})
                %   disp('debug');
                % end
                % neuralNet.LW{ii,jj} = connWeight(iConnId);
                
                % if jj == 18 && ii == 24
                %    disp(' ')
                % end
                
                % Sometimes, NEAT will create two same connections but with
                % different weights. When this happens, neuralNet.LW{ii,jj}
                % fails because there are more than one weight. A solution could
                % be to add the two weights.
                neuralNet.LW{ii,jj} = sum(connWeight(iConnId));
            end
        end
    end
    
end

% Weights for reccurent connections
nDelays = numel(recurrentConnection(:,1));
for ii = 1:nDelays
    % neuralNet.LW{TO ith layer, FROM jth layer}
    iLayer = recurrentConnection(ii,2);
    jLayer = recurrentConnection(ii,1);
    iConnId = (connFrom == layerList(jLayer,2)) & (connTo == layerList(iLayer,2));
    % We sum in case NEAT decides to duplicate the recurrent connection
    neuralNet.LW{iLayer,jLayer} = sum(connWeight(iConnId));
end

%% Old

% % Find recurrent connections
% for iNode = 1:nNodes
%     if nodeType(iNode) ~= 4 % If is not bias node
%
%         % We use a recursive call to find_recurrence() to find if there's a
%         % reccurence involded with the current node
%         recurrenceFlag = false;
%         fromNode = -1;
%         previousINodeId = -1;
%         alreadyVisited = []; % List of already visited nodes
%         parentNodes = []; % List of parents nodes
%         [recurrenceFlag, fromNode, alreadyVisited] = ...
%             find_recurrence(nodeId(iNode), nodeId(iNode), previousINodeId, ...
%             nodeId, nodeType, layerList, connFrom, connTo, alreadyVisited, parentNodes,...
%             recurrenceFlag, fromNode, true);
%
%         % Set layers delays for reccurent connections
%         % neuralNet.layerWeights{TO ith layer, FROM jth layer}.delays
%         if fromNode(1) ~= -1
%             for kk = 1:numel(fromNode)
%                 jj = layerList(layerList(:,2) == nodeId(iNode),1);
%                 ii = layerList(layerList(:,2) == fromNode(kk),1);
%                 neuralNet.layerWeights{ii,jj}.delays = 1;
%             end
%         end
%
%     end
% end
%
% % Initialize neural net
% neuralNet = init(neuralNet);
%
% end
%
% %% find_recurrence()
% function [recurrenceFlag, fromNode, alreadyVisited] = ...
%     find_recurrence(iNodeIdOrigin, iNodeId, previousINodeId, nodeId, nodeType, layerList, connFrom, ...
%     connTo, alreadyVisited, parentNodes, recurrenceFlag, fromNode, firstCall)
%
% % Check if the node is the origin node which means that it is reccurent
% if recurrenceFlag == false
%     if (iNodeId == iNodeIdOrigin) && ...
%             (firstCall == false) && ...
%             (nodeType(nodeId == previousINodeId) ~= 2)
%         recurrenceFlag = false;
%         fromNode = [fromNode(fromNode ~= -1),previousINodeId];
%         return;
%     elseif (iNodeId == iNodeIdOrigin) && ...
%             (firstCall == false) && ...
%             (nodeType(nodeId == previousINodeId) == 2)
%         recurrenceFlag = true;
%         fromNode = [fromNode(fromNode ~= -1),previousINodeId];
%         return;
%     else
%         recurrenceFlag = false;
%         if fromNode(1) == -1
%             fromNode = -1;
%         end
%     end
% else
%     return;
% end
%
% % ID of nodes having connection from this node
% iFromConn = connTo(connFrom == iNodeId);
%
% % If iNodeId is contained in iFromConn, we remove it from iFromConn unless
% % iNodeId is equal to iNodeIdOrigin. This is to avoid infinite recursive calls
% if (iNodeId ~= iNodeIdOrigin) && (sum(iFromConn == iNodeId) > 0)
%     iFromConn = iFromConn(iFromConn ~= iNodeId);
% end
%
% % We also remove parent nodes contained in iFromConn to avoid infinite recursive
% % calls
% if ~isempty(parentNodes)
%     sumTmp = 0;
%     for iParentNode = 1:numel(parentNodes)
%         sumTmp = sumTmp + sum(parentNodes(iParentNode) == iFromConn);
%     end
%     if sumTmp > 0
%         for iParentNode = 1:numel(parentNodes)
%             iFromConn = ...
%                 iFromConn(...
%                 (parentNodes(iParentNode) ~= iFromConn) | ...
%                 (iNodeIdOrigin == iFromConn) ...
%                 );
%         end
%     end
% end
%
% % If iFromConn is empty, then it means that we reached an output without
% % connection to previous node. We only return to previous call of
% % find_recurrence().
% if isempty(iFromConn)
%     return;
% end
%
% % Store parent node
% parentNodes = [parentNodes, iNodeId];
%
% % Cycle through FROM nodes
% for ii = iFromConn
%     previousINodeId = iNodeId;
%     [recurrenceFlag, fromNode, alreadyVisited] = ...
%         find_recurrence(iNodeIdOrigin, ii, previousINodeId, nodeId, nodeType, layerList, connFrom, ...
%         connTo, alreadyVisited, parentNodes, recurrenceFlag, fromNode, false);
% end
%
% % Update list of already visited nodes
% alreadyVisited = [alreadyVisited, iNodeId];
%

end
