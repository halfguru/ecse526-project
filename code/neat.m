%% neat() (main)
function [Population, bestIndividual, varargout] = neat(...
    nMaxGeneration,...
    loadGene,...
    geneFile,...
    loadFlag,...
    saveFlag,...
    populationSize,...
    nInputNodes,...
    nOutputNodes,...
    initialConnectedInputNodes,...
    fitnessLimit,...
    Speciation,...
    Initial,...
    Stagnation, ...
    Refocus, ...
    Selection, ...
    Crossover, ...
    Mutation, ...
    evalFunHandle, ...
    EvalFunParameters, ...
    dateAndTime, ...
    varargin)

%% Initialization

% Timing
tic;

% Initialize vectors
nAverageNonDisabledConnections = [];
nAverageHiddenNodes = [];
maxOverallFitness = [];
nAverageDisabledConnections = [];
nSpeciesArray = [];
meanOverallFitness = [];
meanRecurrentConnection = [];
nIndPerSpecie = [];

% Initialize Speciation structure. It will contain various information on
% single species. This data will be used for fitness sharing, reproduction,
% and for visualisation purposes.
SpeciesRecord(1).id = 0; % Consecutive species ID's
SpeciesRecord(1).nIndividuals = 0; % Number of individuals in species
% generationRecord matrix will be 4 rows by (number of generations
% existent) columns, will contain (from top to bottom):
% - Number of generation
% - Mean raw fitness
% - Max raw fitness
% - Index of individual in population which has produced max raw fitness
SpeciesRecord(1).generationRecord = [];

if loadFlag == 0 && loadGene == 0
    
    % Call function to create Initial population
    % for information about the make-up of the population structure and the
    % innovationRecord, look at initialize_population()
    [Population, innovationRecord] = initialize_population(...
        populationSize, nInputNodes, nOutputNodes, initialConnectedInputNodes);
    
    % Initial speciation
    [Population, SpeciesRecord] = ...
        initial_speciation(initialConnectedInputNodes, Population, nOutputNodes, ...
        SpeciesRecord, Speciation, loadGene);
    
    % Generation counter
    iGeneration = 1;
    
elseif loadFlag == 0 && loadGene == 1
    
    % Load gene
    TmpStruct = load(geneFile,'Population');
    PopulationTmp = TmpStruct.Population;
    if numel([PopulationTmp(:).fitness]) == 1
        startNodeGenes       = PopulationTmp.nodeGenes;
        startConnectionGenes = PopulationTmp.connectionGenes;
    else
        tmpFitness           = [PopulationTmp(:).fitness];
        [~, iTopIndividual]  = max(tmpFitness);
        bestIndividual       = PopulationTmp(iTopIndividual);
        startNodeGenes       = bestIndividual.nodeGenes;
        startConnectionGenes = bestIndividual.connectionGenes;
    end
    
    % Initialize population with gene
    [Population,innovationRecord] = initialize_population_custom_gene(...
        populationSize,...
        startNodeGenes,...
        startConnectionGenes);
    
    % Initial speciation
    [Population, SpeciesRecord] = ...
        initial_speciation(initialConnectedInputNodes, Population, nOutputNodes, ...
        SpeciesRecord, Speciation, loadGene);
    
    % Generation counter
    iGeneration = 1;
    
else
    
    % Backup mutation schedulings
    probabilityAddNodeSchedulingTmp       = Mutation.probabilityAddNodeScheduling;
    probabilityAddConnectionSchedulingTmp = Mutation.probabilityAddConnectionScheduling;
    if true
        weightCapTmp   = Mutation.weightCap;
        weightRangeTmp = Mutation.weightRange;
    end
    
    % Start with saved version of evolution
    load('neatsave.mat');
    
    % Restore mutation schedulings
    Mutation.probabilityAddNodeScheduling       = probabilityAddNodeSchedulingTmp;
    Mutation.probabilityAddConnectionScheduling = probabilityAddConnectionSchedulingTmp;
    if true
        Mutation.weightCap   = weightCapTmp;
        Mutation.weightRange = weightRangeTmp;
    end
    
end

if Mutation.probabilityAddNodeScheduling(1) ~= -1
    Mutation.probabilityAddNode = Mutation.probabilityAddNodeScheduling(2,1);
end

if Mutation.probabilityAddConnectionScheduling(1) ~= -1
    Mutation.probabilityAddConnection = Mutation.probabilityAddConnectionScheduling(2,1);
end

%% Generational Loop

print_generation(iGeneration);
flagSolution = 0;
while (iGeneration < nMaxGeneration) && (flagSolution == 0)
    
    % Backup copies of current generation
    if saveFlag == 1
        save('neatsave.mat','Population','iGeneration','innovationRecord',...
            'SpeciesRecord','nAverageNonDisabledConnections','nAverageHiddenNodes',...
            'maxOverallFitness','dateAndTime','Mutation','nAverageDisabledConnections',...
            'nSpeciesArray','meanOverallFitness','meanRecurrentConnection','nIndPerSpecie');
        if exist('figHandle','var')
            if ishandle(figHandle)
                savefig(figHandle,'neatfig.fig');
            end
        end
    end
    
    % Call evaluation function (in this case XOR), fitnesses of individuals
    % will be stored in Population(:).fitness.
    % IMPORTANT: reproduction assumes an (all positive!) evaluation function
    % where a higher value means better fitness (in other words, the
    % algorithm is geared towards maximizing a fitness function which can
    % only assume values between 0 and +Inf)
    % Population = xor_experiment(Population);
    % Population = fulladder_experiment(Population);
    if isempty(EvalFunParameters)
        Population = evalFunHandle(Population, fitnessLimit);
    else
        nSimPoints = 5000;
        SimData = initialize_sim_data_structure(EvalFunParameters, populationSize, ...
            nSimPoints);
        EvalFunParameters.iGeneration = iGeneration;
        [Population, SimData] = evalFunHandle(Population, EvalFunParameters, ...
            fitnessLimit, SimData);
    end
    
    % Check for stagnation
    [maxFitnessesCurrentGeneration, SpeciesRecord] = ...
        stagnation(SpeciesRecord, Population, Stagnation, iGeneration);
    
    % Check for refocus
    [SpeciesRecord] = ...
        refocus(maxFitnessesCurrentGeneration, SpeciesRecord, Refocus);
    
    % Save current generation
    folderName = 'neat_data';
    if ~exist(folderName,'dir')
        mkdir(folderName);
    end
    save([folderName,'/neatsave_',dateAndTime,'_',num2str(iGeneration,'%03.0f'),'.mat'],...
        'Population','iGeneration','innovationRecord','SpeciesRecord');
    
    % Compute stats
    [nAverageNonDisabledConnections, nAverageHiddenNodes, maxFitness,...
        maxOverallFitness, nAverageDisabledConnections, ...
        nSpeciesArray, meanOverallFitness, meanRecurrentConnection, nIndPerSpecie] = ...
        compute_stats(Population, populationSize, iGeneration, ...
        SpeciesRecord, nAverageNonDisabledConnections, nAverageHiddenNodes, ...
        maxOverallFitness, nAverageDisabledConnections, ...
        nSpeciesArray, meanOverallFitness, meanRecurrentConnection, nIndPerSpecie);
    
    % Create figure handle and axes if don't exist or have been deleted
    if exist('figHandle','var')
        [figHandle, AxHandle] = neat_plot_helper(loadFlag, figHandle, AxHandle);
        figHandle.CloseRequestFcn = @closereqmod;
        figHandle.Resize = 'off';
    else
        [figHandle, AxHandle] = neat_plot_helper(loadFlag);
        figHandle.CloseRequestFcn = @closereqmod;
        figHandle.Resize = 'off';
    end
    
    % Visualisation fitness & species
    visualization(figHandle, AxHandle, nAverageNonDisabledConnections, ...
        nAverageHiddenNodes, maxOverallFitness, nAverageDisabledConnections, ...
        nSpeciesArray, meanOverallFitness, meanRecurrentConnection, nIndPerSpecie);
    
    % Find best individual of this generation
    tmpFitness = [Population(:).fitness];
    [~, iTopIndividual] = max(tmpFitness);
    bestIndividual = Population(iTopIndividual);
    
    % Visualisation of simulation results of the best individual
    if exist('SimData','var')
        if ~exist('figHandleSim','var'); figHandleSim = []; end
        if ~exist('figHandleSimError','var'); figHandleSimError = []; end
        [figHandleSim, figHandleSimError] = ...
            sim_visualization(figHandleSim, figHandleSimError, SimData, ...
            iTopIndividual, EvalFunParameters);
        figHandleSim.CloseRequestFcn = @closereqmod;
        figHandleSimError.CloseRequestFcn = @closereqmod;
        figHandleSim.Resize = 'off';
        figHandleSimError.Resize = 'off';
    end
    
    % Visualisation of neural network of the best individual
    if exist('figHandleNeuralNet','var')
        figHandleNeuralNet = nn_visualization(Population, figHandleNeuralNet);
        figHandleNeuralNet.CloseRequestFcn = @closereqmod;
        figHandleNeuralNet.Resize = 'off';
    else
        figHandleNeuralNet = nn_visualization(Population);
        figHandleNeuralNet.CloseRequestFcn = @closereqmod;
        figHandleNeuralNet.Resize = 'off';
    end
    
    % Check for solution
    flagSolution = solution_check(maxFitness, fitnessLimit, flagSolution);
    
    % We call reproduce if there's not solution yet
    if flagSolution == 0
        % Call reproduction function with parameters, current population
        % and species record, returns new population, new species record
        % and new innovation record
        [Population, SpeciesRecord, innovationRecord] = ...
            reproduce(Population, SpeciesRecord, innovationRecord, Initial,...
            Selection, Crossover, Mutation, Speciation, iGeneration, populationSize);
        % Display time
        display_elapsed_time(toc());
    end
    
    % Increment generational counter
    iGeneration = iGeneration + 1;
    
    % Print generation to command window
    print_generation(iGeneration);
    
    % Reduce mutation probabilities (if shrink rate > 0)
    Mutation = shrink_mutation(Mutation, iGeneration);
    
end

end

%% initialize_sim_data_structure()
function SimData = ...
    initialize_sim_data_structure(EvalFunParameters, populationSize, nSimPoints)

cellArrayTmp = cell(1,numel(EvalFunParameters.Trajectory));
[cellArrayTmp{:}]                           = deal(zeros(nSimPoints,1));
[SimData(1:populationSize).t]               = deal(cellArrayTmp);
[SimData(1:populationSize).xDesired]        = deal(cellArrayTmp);
[SimData(1:populationSize).yDesired]        = deal(cellArrayTmp);
[SimData(1:populationSize).zDesired]        = deal(cellArrayTmp);
[SimData(1:populationSize).psiDesired]      = deal(cellArrayTmp);
[SimData(1:populationSize).thetaDesired]    = deal(cellArrayTmp);
[SimData(1:populationSize).phiDesired]      = deal(cellArrayTmp);
[SimData(1:populationSize).x]               = deal(cellArrayTmp);
[SimData(1:populationSize).y]               = deal(cellArrayTmp);
[SimData(1:populationSize).z]               = deal(cellArrayTmp);
[SimData(1:populationSize).psi]             = deal(cellArrayTmp);
[SimData(1:populationSize).theta]           = deal(cellArrayTmp);
[SimData(1:populationSize).phi]             = deal(cellArrayTmp);
[SimData(1:populationSize).Va1]             = deal(cellArrayTmp);
[SimData(1:populationSize).Va2]             = deal(cellArrayTmp);
[SimData(1:populationSize).Va3]             = deal(cellArrayTmp);
[SimData(1:populationSize).Va4]             = deal(cellArrayTmp);

end

%% display_elapsed_time()
function display_elapsed_time(elapsedTime)

% elapsedTime = toc;
elapsedTimeHours = floor(elapsedTime/3600);
elapsedTimeMins  = floor((elapsedTime - 3600*elapsedTimeHours)/60);
elapsedTimeSecs  = floor((elapsedTime - 3600*elapsedTimeHours - 60*elapsedTimeMins));
fprintf('Elapsed time is %02.0fh%02.0fm%02.0fs\n\n',...
    elapsedTimeHours, elapsedTimeMins, elapsedTimeSecs);

end

%% sim_visualization()
function [figHandleSim, figHandleSimError] = ...
    sim_visualization(figHandleSim, figHandleSimError, SimData, iTopIndividual, ...
    EvalFunParameters)

if exist('SimData','var')
    
    if ~isempty(figHandleSim) && ~isempty(figHandleSimError)
        if ishandle(figHandleSim) && ishandle(figHandleSimError)
            clf(figHandleSim);
            clf(figHandleSimError);
        elseif ~ishandle(figHandleSim) || ~ishandle(figHandleSimError)
            if ishandle(figHandleSim); close(figHandleSim); end;
            if ishandle(figHandleSimError); close(figHandleSimError); end;
            delete figHandleSim;
            delete figHandleSimError;
        end
    end
    
    nTraj = numel(EvalFunParameters.Trajectory);
    [Data(1:nTraj).t]               = deal(SimData(iTopIndividual).t{:});
    [Data(1:nTraj).xDesired]        = deal(SimData(iTopIndividual).xDesired{:});
    [Data(1:nTraj).yDesired]        = deal(SimData(iTopIndividual).yDesired{:});
    [Data(1:nTraj).zDesired]        = deal(SimData(iTopIndividual).zDesired{:});
    [Data(1:nTraj).psiDesired]      = deal(SimData(iTopIndividual).psiDesired{:});
    [Data(1:nTraj).thetaDesired]    = deal(SimData(iTopIndividual).thetaDesired{:});
    [Data(1:nTraj).phiDesired]      = deal(SimData(iTopIndividual).phiDesired{:});
    [Data(1:nTraj).x]               = deal(SimData(iTopIndividual).x{:});
    [Data(1:nTraj).y]               = deal(SimData(iTopIndividual).y{:});
    [Data(1:nTraj).z]               = deal(SimData(iTopIndividual).z{:});
    [Data(1:nTraj).psi]             = deal(SimData(iTopIndividual).psi{:});
    [Data(1:nTraj).theta]           = deal(SimData(iTopIndividual).theta{:});
    [Data(1:nTraj).phi]             = deal(SimData(iTopIndividual).phi{:});
    [Data(1:nTraj).Va1]             = deal(SimData(iTopIndividual).Va1{:});
    [Data(1:nTraj).Va2]             = deal(SimData(iTopIndividual).Va2{:});
    [Data(1:nTraj).Va3]             = deal(SimData(iTopIndividual).Va3{:});
    [Data(1:nTraj).Va4]             = deal(SimData(iTopIndividual).Va4{:});
    
    customLegend = cell(1,nTraj);
    for jTraj = 1:nTraj; customLegend{jTraj} = ['Traj ',num2str(jTraj)]; end;
    
    if ~isempty(figHandleSim) && ~isempty(figHandleSimError) && ...
            ishandle(figHandleSim) && ishandle(figHandleSimError)
        [figHandleSim, figHandleSimError] = ...
            plot_uav_sim(Data, true, false, false, true, false, '', '', true, ...
            false, customLegend, true, figHandleSim, figHandleSimError);
        drawnow;
    else
        [figHandleSim, figHandleSimError] = ...
            plot_uav_sim(Data, true, false, false, true, false, '', '', true, ...
            false, customLegend, true);
        drawnow;
    end
    
end

end

%% nn_visualization()
function figHandleNeuralNet = nn_visualization(Population, varargin)

if nargin == 2
    figHandleNeuralNet = varargin{1};
end

if exist('figHandleNeuralNet','var') && ishandle(figHandleNeuralNet)
    clf(figHandleNeuralNet);
end

if exist('figHandleNeuralNet','var') && ishandle(figHandleNeuralNet)
    figHandleNeuralNet = plot_neural_network([], ...
        [], [], [], 1, true, false, Population, figHandleNeuralNet);
    drawnow;
else
    figHandleNeuralNet = plot_neural_network([], ...
        [], [], [], 1, true, false, Population);
    drawnow;
end

end

%% neat_plot_helper()
function [figHandle, AxHandle] = neat_plot_helper(loadFlag, varargin)

% This function simply manages the figure and axes of neat plot.

if nargin == 3
    figHandle = varargin{1};
    AxHandle = varargin{2};
end

marginSubPlot = 0.10;
subplotSettings = {marginSubPlot, 'NextPlot', 'replacechildren', 'Tag'};

if ~exist('figHandle','var') && loadFlag ~= 1
    
    figHandle = create_figure(800, 700, 'on');
    AxHandle.subPlot1 = subplot_tight(3,3,[1,2], subplotSettings{:}, 'subPlot1');
    AxHandle.subPlot2 = subplot_tight(3,3,3, subplotSettings{:}, 'subPlot2');
    AxHandle.subPlot3 = subplot_tight(3,3,4, subplotSettings{:}, 'subPlot3');
    AxHandle.subPlot4 = subplot_tight(3,3,5, subplotSettings{:}, 'subPlot4');
    AxHandle.subPlot5 = subplot_tight(3,3,6, subplotSettings{:}, 'subPlot5');
    AxHandle.subPlot6 = subplot_tight(3,3,7, subplotSettings{:}, 'subPlot6');
    AxHandle.subPlot7 = subplot_tight(3,3,8, subplotSettings{:}, 'subPlot7');
    AxHandle.subPlot8 = subplot_tight(3,3,9, subplotSettings{:}, 'subPlot8');
    
elseif ~exist('figHandle','var') && loadFlag == 1
    
    figHandle = openfig('neatfig.fig');
    AxHandle.subPlot1 = findobj(figHandle, 'Type', 'axes', 'Tag', 'subPlot1');
    AxHandle.subPlot2 = findobj(figHandle, 'Type', 'axes', 'Tag', 'subPlot2');
    AxHandle.subPlot3 = findobj(figHandle, 'Type', 'axes', 'Tag', 'subPlot3');
    AxHandle.subPlot4 = findobj(figHandle, 'Type', 'axes', 'Tag', 'subPlot4');
    AxHandle.subPlot5 = findobj(figHandle, 'Type', 'axes', 'Tag', 'subPlot5');
    AxHandle.subPlot6 = findobj(figHandle, 'Type', 'axes', 'Tag', 'subPlot6');
    AxHandle.subPlot7 = findobj(figHandle, 'Type', 'axes', 'Tag', 'subPlot7');
    AxHandle.subPlot8 = findobj(figHandle, 'Type', 'axes', 'Tag', 'subPlot8');
    
end

% In case the figure has been closed
if exist('figHandle','var') && ~ishandle(figHandle)
    
    fprintf('\nNEAT plot window has been closed. It will reopened the next generation.\n');
    delete AxHandle; clear AxHandle;
    delete figHandle; clear figHandle;
    pause(0.1);
    figHandle = create_figure(800, 700, 'on');
    drawnow;
    set(0, 'currentfigure', figHandle);
    AxHandle.subPlot1 = subplot_tight(3,3,[1,2], subplotSettings{:}, 'subPlot1');
    set(0, 'currentfigure', figHandle);
    AxHandle.subPlot2 = subplot_tight(3,3,3, subplotSettings{:}, 'subPlot2');
    set(0, 'currentfigure', figHandle);
    AxHandle.subPlot3 = subplot_tight(3,3,4, subplotSettings{:}, 'subPlot3');
    set(0, 'currentfigure', figHandle);
    AxHandle.subPlot4 = subplot_tight(3,3,5, subplotSettings{:}, 'subPlot4');
    set(0, 'currentfigure', figHandle);
    AxHandle.subPlot5 = subplot_tight(3,3,6, subplotSettings{:}, 'subPlot5');
    set(0, 'currentfigure', figHandle);
    AxHandle.subPlot6 = subplot_tight(3,3,7, subplotSettings{:}, 'subPlot6');
    set(0, 'currentfigure', figHandle);
    AxHandle.subPlot7 = subplot_tight(3,3,8, subplotSettings{:}, 'subPlot7');
    set(0, 'currentfigure', figHandle);
    AxHandle.subPlot8 = subplot_tight(3,3,9, subplotSettings{:}, 'subPlot8');
    
end

end

%% initialize_population()
function [Population,innovationRecord] = initialize_population(...
    populationSize,...
    nInputNodes,...
    nOutputNodes,...
    initialConnectedInputNodes)

% nodeGenes is array 4 rows * (nInputNodes + nOutputNodes +
% hidden-nodes (not existent in initial population) + 1 (bias-node)) columns
% nodeGenes contains:
% - Consecutive node ID's (upper row)
% - Node type (lower row) 1 = input 2 = output 3 = hidden 4 = bias
% - Node input state
% - Node output state (used for evaluation, all input states zero
%   initially, except bias node, which is always 1)
%
% connectionGenes is array 5 rows * nConnections columns
% from top to bottom, those five rows contain:
% - Innovation number
% - Connection from
% - Connection to
% - Weight
% - Enable bit
% The rest of the elements in the structure for an individual should be
% self-explanatory
%
% innovationRecord tracks innovations in a 5 rows by (number of
% innovations) columns matrix, contains:
% - innovation number
% - connectFromNode for this innovation
% - connectToNode for this innovation
% - The new node (if it is a new node mutation, then this node will appear in
%   the 4th row when it is first connected. There will always be two
%   innovations with one node mutation, since there is a connection to and
%   from the new node. In the initial population, this will be abbreviated to
%   the node with the highest number appearing in the last column of the record,
%   since only this is needed as starting point for the rest of the algorithm)
% - 5th row is generation this innovation occured (generation is assumed
%   to be zero for the innovations in the initial population)

% Compute number and matrix of initial connections (all connections between
% output nodes and the nodes listed in initialConnectedInputNodes). In
% nConnections, we add 1 since the bias is connected. vectorConnectionFrom is a
% row vector (i.e. one row, multiple columns) consists of a nOutputNodes times
% repeated row vector containing the ID of the initial connected input plus the
% bias node ID. For example, if there are 3 inputs (all connected) and 3
% outputs, the vector will be:
% vectorConnectionFrom = [[1,2,3,4],[1,2,3,4],[1,2,3,4]]
% where 4 is the bias' node ID.
nConnections         = (length(initialConnectedInputNodes)+1)*nOutputNodes;
vectorConnectionFrom = repmat([initialConnectedInputNodes,nInputNodes+1],[1,nOutputNodes]);
% vectorConnectionTo contains the node ID to which are connected every node
% listed in vectorConnectionFrom. For instance, if we take the example above
% with 3 inputs (all connected) and 3 outputs, we'd have:
% vectorConnectionTo = [[5,5,5,5],[6,6,6,6],[7,7,7,7]]
vectorConnectionTo   = [];
% Outputs' node IDs begin at the ID after the bias' ID which is why we take
% nInputNodes + 2. For instance, if we have 3 inputs + 1 bias, the IDs 1,2,3,4
% are taken, so we start at 3 + 2 = 5 for the first output ID.
iOutputNodeStart = nInputNodes + 2;
iOutputNodeEnd   = nInputNodes + 1 + nOutputNodes;
for iOutputNode  = iOutputNodeStart:iOutputNodeEnd
    % There is a +1 because the bias is connected.
    vectorConnectionTo = ...
        [vectorConnectionTo, ...
        iOutputNode*ones(1,length(initialConnectedInputNodes)+1)];
end
% Connection matrix is a matrix of two rows by (number of connections) columns
connectionMatrix = [...
    vectorConnectionFrom;
    vectorConnectionTo];

% We build the initial population nodeGenes and codeGenes
for iIndividual = 1:populationSize
    Population(iIndividual).nodeGenes = [...
        1:(nInputNodes+1+nOutputNodes);
        ones(1,nInputNodes),4,2*ones(1,nOutputNodes);
        zeros(1,nInputNodes),1,zeros(1,nOutputNodes);
        zeros(1,nInputNodes+1+nOutputNodes)];
    % All weights uniformly distributed in [-1,+1], all connections enabled
    Population(iIndividual).connectionGenes = [...
        1:nConnections;
        connectionMatrix;
        rand(1,nConnections)*2-1;
        ones(1,nConnections)];
    Population(iIndividual).fitness = 0;
    Population(iIndividual).species = 0;
end
innovationRecord = [...
    Population(populationSize).connectionGenes(1:3,:);
    zeros(size(Population(populationSize).connectionGenes(1:2,:)))];
% Highest node ID for initial population
innovationRecord(4,size(innovationRecord,2)) = max(Population(1).nodeGenes(1,:));

end

%% initialize_population_custom_gene()
function [Population,innovationRecord] = initialize_population_custom_gene(...
    populationSize,...
    startNodeGenes,...
    startConnectionGenes)

% See initialize_population() for documentation

% Number of connections
nConnections     = numel(startConnectionGenes(1,:));

% We build the initial population nodeGenes and codeGenes
for iIndividual = 1:populationSize
    Population(iIndividual).nodeGenes = startNodeGenes;
    % All weights uniformly distributed in [-1,+1]
    Population(iIndividual).connectionGenes = [...
        startConnectionGenes(1:3,:);
        rand(1,nConnections)*2-1;
        startConnectionGenes(5,:);];
    Population(iIndividual).fitness = 0;
    Population(iIndividual).species = 0;
end
innovationRecord = [...
    Population(populationSize).connectionGenes(1:3,:);
    zeros(size(Population(populationSize).connectionGenes(1:2,:)))];
% Highest node ID for initial population
innovationRecord(4,size(innovationRecord,2)) = max(Population(1).nodeGenes(1,:));

end

%% initial_speciation()
function [Population, SpeciesRecord] = ...
    initial_speciation(initialConnectedInputNodes, Population, nOutputNodes, ...
    SpeciesRecord, Speciation, loadGene)

% This function performs initial speciation. That is, it assigns each
% individual in a specie according to their compatibility distance.

% Number of connections
if loadGene == 0
    nConnections = (length(initialConnectedInputNodes)+1)*nOutputNodes;
else
    nConnections = numel(Population(1).connectionGenes(1,:));
end
% Put first individual in species one and update SpeciesRecord
Population(1).species = 1;
% Species reference matrix (abbreviated, only weights, since there are
% no topology differences in initial population)
matrixReferenceIndividuals = Population(1).connectionGenes(4,:);
SpeciesRecord(1).id = 1;
SpeciesRecord(1).nIndividuals = 1;

% Loop through rest of individuals and either assign to existing species
% or create new species and use first individual of new species as
% reference
for iIndividual = 2:size(Population,2);
    assignedExistingSpeciesFlag = 0;
    newSpeciesFlag = 0;
    iSpecie = 1;
    % Loops through the existing species, terminates when either the
    % individual is assigned to existing species or there are no more
    % species to test it against, which means it is a new species
    while (assignedExistingSpeciesFlag == 0) && (newSpeciesFlag == 0)
        % Computes compatibility distance, abbreviated, only average
        % weight distance considered
        averageWeightDiffMatchingGenes = ...
            sum(...
            abs(Population(iIndividual).connectionGenes(4,:) - ...
            matrixReferenceIndividuals(iSpecie,:)) ...
            );
        distance = Speciation.c3*averageWeightDiffMatchingGenes/nConnections;
        % If within threshold, assign to the existing species
        if distance < Speciation.threshold
            Population(iIndividual).species = iSpecie;
            assignedExistingSpeciesFlag = 1;
            SpeciesRecord(iSpecie).nIndividuals = SpeciesRecord(iSpecie).nIndividuals+1;
        end
        iSpecie = iSpecie + 1;
        % Outside of species references, must be new species
        if (iSpecie > size(matrixReferenceIndividuals,1)) && (assignedExistingSpeciesFlag == 0)
            newSpeciesFlag = 1;
        end
    end
    % Check for new species, if it is, update the SpeciesRecord and
    % use individual as reference for new species
    if newSpeciesFlag == 1
        Population(iIndividual).species = iSpecie;
        matrixReferenceIndividuals = [...
            matrixReferenceIndividuals;
            Population(iIndividual).connectionGenes(4,:)];
        SpeciesRecord(iSpecie).id = iSpecie;
        % If number individuals in a species is zero, that species is extinct
        SpeciesRecord(iSpecie).nIndividuals = 1;
    end
end

end

%% stagnation()
function [maxFitnessesCurrentGeneration, SpeciesRecord] = ...
    stagnation(SpeciesRecord, Population, Stagnation, iGeneration)

% This function checks for stagnation. We also compute mean and max raw
% fitnesses in each species and store in SpeciesRecord.generationRecord

% Number of species in current generation
nSpecies = size(SpeciesRecord,2);

% Initialize maxFitnessesCurrentGeneration vector
maxFitnessesCurrentGeneration = zeros(1,nSpecies);

% The vector species contains the species index of each individual in the
% population
popSpecies = [Population(:).species];

% The vector fitness contains the fitness of each individual in the
% population
popFitness = [Population(:).fitness];

% Cycle through each specie
for iSpecie = 1:nSpecies
    if SpeciesRecord(iSpecie).nIndividuals > 0
        % Max fitness of current specie
        [maxFitness, iMaxFitness] = max((popSpecies == iSpecie).*popFitness);
        % Mean fitness of current specie
        meanFitness = ...
            sum((popSpecies == iSpecie).*popFitness)/SpeciesRecord(iSpecie).nIndividuals;
        % Number of generation in current specie
        nGenerations = size(SpeciesRecord(iSpecie).generationRecord,2);
        % Compute stagnation vector (last Stagnation.nGenerations-1 max
        % fitnesses plus current fitness)
        if nGenerations > Stagnation.nGenerations-2
            iGenerationStart = nGenerations - Stagnation.nGenerations + 2;
            lastMaxFitnesses = ...
                SpeciesRecord(iSpecie).generationRecord(3,iGenerationStart:nGenerations);
            stagnationVector = [lastMaxFitnesses,maxFitness];
            % Average max fitnesses change for the last Stagnation.nGenerations generations
            averageMaxFitnessesChange = abs(stagnationVector-mean(stagnationVector));
            % Number of stagned generation in current specie
            nStagnedGenerations = sum(averageMaxFitnessesChange < Stagnation.threshold);
            % Check for stagnation
            if nStagnedGenerations == Stagnation.nGenerations
                % Set mean fitness to small value to eliminate species
                % (cannot be set to 0, if only one species is present,
                % we would have divide by zero in fitness sharing.
                % anyways, with only one species present, we have to keep it)
                meanFitness = 0.01;
            end
        end
        % We add new generation to current specie
        SpeciesRecord(iSpecie).generationRecord = [...
            SpeciesRecord(iSpecie).generationRecord,...
            [iGeneration;meanFitness;maxFitness;iMaxFitness]];
        maxFitnessesCurrentGeneration(1,iSpecie) = maxFitness;
    end
end

end

%% refocus()
function [SpeciesRecord] = ...
    refocus(maxFitnessesCurrentGeneration, SpeciesRecord, Refocus)

% In rare cases when the fitness of the entire population does not
% improve for more than Refocus.nGenerations generations, only the top
% two species are allowed to reproduce, refocusing the search into the
% most promising spaces

% Index of the top specie
[~, iTopSpecie] = max(maxFitnessesCurrentGeneration);

% Number of generations in the top specie
nGenerationsTopSpecie = size(SpeciesRecord(iTopSpecie).generationRecord,2);

% Check for refocus
if nGenerationsTopSpecie > Refocus.nGenerations
    index1 = nGenerationsTopSpecie-Refocus.nGenerations;
    index2 = nGenerationsTopSpecie;
    % Max fitnesses from (nGenerationsTopSpecie-Refocus.nGenerations) to
    % nGenerationsTopSpecie of top specie
    maxFitnessesTopSpecie = SpeciesRecord(iTopSpecie).generationRecord(3,index1:index2);
    % Average max fitnesses change for the top specie over the last
    % nGenerationsTopSpecie generations
    averageMaxFitnessesChangeTopSpecie = abs(maxFitnessesTopSpecie - mean(maxFitnessesTopSpecie));
    % Number of stagned generation in current specie
    nStagnedGenerationsTopSpecie = sum(averageMaxFitnessesChangeTopSpecie < Refocus.threshold);
    % Check if number of stagned generations in the top specie is equal
    % to the threshold refocus generations number
    if nStagnedGenerationsTopSpecie == Refocus.nGenerations
        % Sort max species' fitnesses from high fitnesses to low fitnesses
        [~, vectorCull] = sort(-maxFitnessesCurrentGeneration);
        % Index of the species to discard
        vectorCull = vectorCull(1,3:sum(maxFitnessesCurrentGeneration > 0));
        for iSpecie = 1:size(vectorCull,2)
            % Index of specie to discard
            iCull = vectorCull(1,iSpecie);
            % Index of the the last generation of the current specie to
            % discard
            iGenerationCull = size(SpeciesRecord(iCull).generationRecord,2);
            % We assign a mean raw fitness of 0.01 to the last
            % generation of the current specie to discard
            SpeciesRecord(iCull).generationRecord(2,iGenerationCull) = 0.01;
        end
    end
end

end

%% reproduce()
function [NewPopulation, updatedSpeciesRecord, updatedInnovationRecord] = ...
    reproduce(Population, SpeciesRecord, innovationRecord, Initial, Selection,...
    Crossover, Mutation, Speciation, iGeneration, populationSize)

%% Initial selection

% Initial selection
[matExistingAndPropagSpecies, NewPopulation, Population, iIndividual] = ...
    intial_selection(SpeciesRecord, Initial, Population, populationSize);

% Generate reference population
populationRef = initialize_ref_population(Population, SpeciesRecord);

fprintf('\nExisting and propagating matrix:\n');
tmpMatExistingAndPropagSpecies = num2str(matExistingAndPropagSpecies,' %4.0f');
fprintf([tmpMatExistingAndPropagSpecies(1,:),'\n']);
fprintf([tmpMatExistingAndPropagSpecies(2,:),'\n']);
fprintf([tmpMatExistingAndPropagSpecies(3,:),'\n\n']);

%% Reproduction

% Cycle through all existing species
for iSpecie = 1:size(matExistingAndPropagSpecies,2)
    
    % Number of individuals in specie
    countIndividualsSpecies = 0;
    
    % This is the ID of species which will be reproduced this cycle.
    % IMPORTANT: iSpecie only has relevance to
    % matExistingAndPropagSpecies, all other mechanisms using
    % species in some way must use speciesId
    speciesId = matExistingAndPropagSpecies(1,iSpecie);
    
    % Linear Ranking and stochastic universal sampling Ranking with
    % Selection.pressure
    [NewChrIx, numberCrossover, numberMutate] = ...
        ranking_and_sampling(Population, Selection, speciesId, iSpecie, ...
        matExistingAndPropagSpecies, Crossover);
    
    % Cycle until actual number of offspring has reached allotted number of offspring
    while matExistingAndPropagSpecies(3,iSpecie) < matExistingAndPropagSpecies(2,iSpecie)
        
        % Increment number of individuals in population
        iIndividual = iIndividual + 1;
        
        % Increment number of individuals in current specie
        countIndividualsSpecies = countIndividualsSpecies + 1;
        
        % Crossover
        if countIndividualsSpecies <= numberCrossover
            % OK we are doing crossover
            NewIndividual = crossover(...
                Population, NewChrIx, countIndividualsSpecies, Crossover, ...
                matExistingAndPropagSpecies);
        else
            % No crossover, just copy a individual of the species and mutate in
            % subsequent steps
            NewIndividual = Population(NewChrIx(numberCrossover+countIndividualsSpecies));
        end
        
        % Hidden nodes culling
        NewIndividual = node_culling(NewIndividual);
        
        % Disabled Genes Mutation
        NewIndividual = enable_gene_mutation(NewIndividual, Mutation);
        
        % Weight Mutation
        NewIndividual = weight_mutation(NewIndividual, Mutation);
        
        % IMPORTANT: The checks for duplicate innovations in the following
        % two types of mutation can only check in the current generation
        
        flag1 = rand() < Mutation.probabilityAddNode;
        if flag1 == 0
            % Add Connection Mutation
            [NewIndividual, innovationRecord] = ...
                add_connection_mutation(NewIndividual, innovationRecord,...
                Mutation, iGeneration, flag1);
        elseif flag1 == 1
            % Add (insert) Node Mutation
            [NewIndividual, innovationRecord] = ...
                add_node_mutation(innovationRecord, iGeneration, NewIndividual);
        end
        
        % Speciation
        [NewIndividual, SpeciesRecord, populationRef] = ...
            speciation(NewIndividual, SpeciesRecord, populationRef, Speciation);
        
        % Add NewIndividual to NewPopulation
        NewPopulation(iIndividual) = NewIndividual;
        
        % Increment species
        matExistingAndPropagSpecies(3,iSpecie) = ...
            matExistingAndPropagSpecies(3,iSpecie) + 1;
        
    end
end

% Final update of SpeciesRecord (can only be done now since old population
% sizes were needed during reproduction cycle)
for iSpecie = 1:size(SpeciesRecord,2)
    SpeciesRecord(iSpecie).nIndividuals = sum([NewPopulation(:).species] == iSpecie);
end

% Assign updated SpeciesRecord to output
updatedSpeciesRecord = SpeciesRecord;

% Assign updated innovationRecord to output
updatedInnovationRecord = innovationRecord;

end

%% compute_sum_average_fitnesses()
function sumAverageFitnesses = compute_sum_average_fitnesses(SpeciesRecord)

% Compute the sum of average fitenesses
sumAverageFitnesses = 0;
for iSpecie = 1:size(SpeciesRecord,2)
    % Number of generation in current specie
    nGenerations = size(SpeciesRecord(iSpecie).generationRecord,2);
    % Mean raw fitness of the last generation of current specie "iSpecie"
    meanRawFitness = SpeciesRecord(iSpecie).generationRecord(2,nGenerations);
    % Sum average fitnesses of specie if number of individuals in this
    % specie is > 0
    sumAverageFitnesses = sumAverageFitnesses + ...
        meanRawFitness*(SpeciesRecord(iSpecie).nIndividuals > 0);
end

end

%% initial_selection()
function [matExistingAndPropagSpecies, NewPopulation, Population, iIndividual] = ...
    intial_selection(SpeciesRecord, Initial, Population, populationSize)

% This function implements this part of the paper:
% Every species is assigned a potentially different number of offspring in
% proportion to the sum of adjusted fitnesses of its member organisms.
% Species then reproduce by first eliminating the lowest performing members
% from the population.
%
% The following 'for' loop has these three objectives:
%
% 1. Compute matrix of existing and propagating species from SpeciesRecord
% (first row), assign their alloted number of offspring from the shared
% fitness (second row), and set their actual number of offspring to zero
% (third row) (will be incremented as new individuals are created from this
% species)
%
% 2. Copy most fit individual in every species with more than
% Initial.numberCopy individuals unchanged into new generation (elitism)
% (But only if species is not dying out, i.e. has at least one individual
% allotted to itself in the new generation) utilizes data stored in
% SpeciesRecord.generationRecord (index of individual in population
% having the highest fitness)
%
% 3. Erase lowest percentage (Initial.killPercentage) in species with more
% than Initial.numberForKill individuals to keep them from taking part in
% reproduction

% Compute sum of average fitenesses
sumAverageFitnesses = compute_sum_average_fitnesses(SpeciesRecord);

% The following two lines only initialize the NewPopulation structure.
% Since its species is set to 0, the rest of the algorithm will not bother
% with it. It gets overwritten as soon as the first really new individual
% is created
NewPopulation(1) = Population(1);
NewPopulation(1).species = 0;

% Cycle through all existing species
overflow = 0;
iIndividual = 0;
matExistingAndPropagSpecies = [];
for iSpecie = 1:size(SpeciesRecord,2)
    % Test if species existed in old generation
    if SpeciesRecord(iSpecie).nIndividuals > 0
        
        % Number of generation in current specie
        nGenerations = size(SpeciesRecord(iSpecie).generationRecord,2);
        % Mean raw fiteness of current generation
        meanRawFitness = SpeciesRecord(iSpecie).generationRecord(2,nGenerations);
        % Compute number of offspring in new generation
        nOffsprings = (meanRawFitness/sumAverageFitnesses)*populationSize;
        overflow = overflow + nOffsprings - floor(nOffsprings);
        % Since new species sizes are fractions, overflow sums up the
        % difference between size and floor(size), and everytime this
        % overflow is above 1, the species gets one additional individual
        % alotted to it
        if overflow >= 1
            nOffsprings = ceil(nOffsprings);
            overflow = overflow - 1;
        else
            nOffsprings = floor(nOffsprings);
        end
        
        % Check to see if species is dying out, only add those species to
        % matExistingAndPropagSpecies which have offspring in the
        % new generation
        if nOffsprings > 0
            % Matrix (objective 1)
            matExistingAndPropagSpecies = ...
                [matExistingAndPropagSpecies,...
                [SpeciesRecord(iSpecie).id;
                nOffsprings;
                0]];
            % Check for condition for objective 2
            if SpeciesRecord(iSpecie).nIndividuals >= Initial.numberCopy
                iIndividual = iIndividual + 1;
                % Objective 2
                NewPopulation(iIndividual) = ...
                    Population(SpeciesRecord(iSpecie).generationRecord(4,nGenerations));
                % Update matExistingAndPropagSpecies
                matExistingAndPropagSpecies(3,size(matExistingAndPropagSpecies,2)) = 1;
            end
        end
        
        % condition1 checks if specie has more individuals than numberForKill
        condition1 = SpeciesRecord(iSpecie).nIndividuals > Initial.killPercentage;
        % condition2 checks if after killing killPercentage, there's at least 2
        % individuals to be able to cross over or at least the single
        % individual in a species with only one individual
        condition2 = ceil(SpeciesRecord(iSpecie).nIndividuals*(1-Initial.killPercentage)) > 2;
        % Check condition for objective 3
        if condition1 && condition2
            % The vector popSpecies contains the species index of each
            % individual in the population
            popSpecies = [Population(:).species];
            % Index in population of individuals in iSpecie
            iIndividualsToKill = find(popSpecies == iSpecie);
            % Matrix containing on first row the index of individual in
            % current specie to kill and their fitness on the second row
            matrixIndividualsSpecies = [...
                iIndividualsToKill;
                [Population(iIndividualsToKill).fitness]];
            % We sort the individual in the specie according to their fitness
            [~, sortingVector] = sort(matrixIndividualsSpecies(2,:));
            % Sorted matrixIndividualsSpecies
            matrixIndividualsSpecies = matrixIndividualsSpecies(:,sortingVector);
            % Index of sorted individual in specie
            sortingVector = matrixIndividualsSpecies(1,:);
            % MATLAB actually doesn't offer a facility for redirecting the
            % pointers which link one element in a structure with the next,
            % so essentially this individual entry in the population
            % structure cannot be erased. Rather, it will have its species'
            % ID set to zero, which has the same effect of removing it from
            % the rest of the reproduction cycle, since all reproduction
            % functions access the population structure through the
            % species' ID and no species has an ID of zero.
            % We kill the Initial.killPercentage first individual in the specie
            nKills = floor(SpeciesRecord(iSpecie).nIndividuals*Initial.killPercentage);
            for ikill = 1:nKills
                % Objective 3
                Population(sortingVector(ikill)).species = 0;
            end
        end
        
    end
end

end

%% initialize_ref_population()
function populationRef = initialize_ref_population(Population, SpeciesRecord)

% Generate reference population of random individuals from every species from
% old population. Cycle through species ID's, and add reference individuals from
% old population. New species of new population will be added during reproduction.

% The vector popSpecies contains the species' index of each individual in
% the population (the index of the specie in which each individual belongs to)
popSpecies = [Population(:).species];

% Number of species in current generation
nSpecies = size(SpeciesRecord,2);

% Create reference population
iRef = 0;
for iSpecieRef = 1:nSpecies
    % Index of individuals in current specie (iSpecieRef) that exists in
    % old population. By exist, we mean the individuals that haven't been
    % killed during intial_selection()
    iExistingIndividual = popSpecies == iSpecieRef;
    % Check if species exists in old population
    if sum(iExistingIndividual) > 0
        iRef = iRef + 1;
        % Number of individuals in population
        nIndividuals = size(Population,2);
        % iRefOld is the index of the individual in the old population that
        % is assigned to the new population as the representative of the
        % specie
        [~, iRefOld] = max(iExistingIndividual.*rand(1,nIndividuals));
        % We assign the chosen individual of the old population in the
        % reference population
        populationRef(iRef) = Population(iRefOld);
    end
end

end

%% ranking_and_sampling()
function [NewChrIx, numberCrossover, numberMutate] = ...
    ranking_and_sampling(Population, Selection, speciesId, iSpecie, ...
    matExistingAndPropagSpecies, Crossover)

% This function compute select the individuals within the current specie
% for reproduction (NewChrIx). It also compute the number of crossovers and
% mutations to perform within the selected individuals in the specie
%
% See http://www.geatbx.com/ver_3_7/selsus.html and
% http://www.geatbx.com/docu/algindex-02.html for explaination

%%  Linear Ranking with Selection.pressure

% Index of individual with specie's ID equal to speciesId
iFitnessesSpecies      = find([Population(:).species] == speciesId);
% Fitnesses of individual with specie's ID equal to speciesId
fitnessesSpecies       = [Population(iFitnessesSpecies).fitness];
% Sorting fitnesses in ascending order
[~, sortedFitnesses]   = sort(fitnessesSpecies);
% Number of individuals in speciesId
nIndividualsInSpecie     = size(fitnessesSpecies,2);
% Ranking according to their fitnesses
ranking                  = zeros(1,nIndividualsInSpecie);
ranking(sortedFitnesses) = 1:nIndividualsInSpecie;
if nIndividualsInSpecie > 1
    FitnV = 2 - Selection.pressure + ...
        2*(ranking-1)*(Selection.pressure-1)/(nIndividualsInSpecie - 1);
    FitnV = FitnV';
else
    FitnV = 2;
end

% Compute number of individuals to be selected (two parents
% required for every offspring through crossover, one for mutation)
numberOverall   = matExistingAndPropagSpecies(2,iSpecie) - ...
    matExistingAndPropagSpecies(3,iSpecie);
numberCrossover = round(Crossover.percentage*numberOverall);
numberMutate    = numberOverall - numberCrossover;
nInd            = nIndividualsInSpecie;
nSel            = 2*numberCrossover + numberMutate;

% Rare case, in which a species with at least Initial.numberCopy
% individuals gets individual copied, but compares poorly to new
% generation, which results in this copied individual being the only
% individual of this species in new generation, so no crossover or
% mutation takes place. setting nSel to 1 will prevent division by zero error,
% but will have no other effect since the propagation loop is governed by
% matExistingAndPropagSpecies, not by nSel
if nSel == 0
    nSel = 1;
end

%% Perform stochastic universal sampling

% (Code Snippet from Genetic Algorithm toolbox 1.2 by Chipperfield et al)

cumFit = cumsum(FitnV);
trials = cumFit(nInd) / nSel * (rand() + (0:nSel-1)');
mF = cumFit(:, ones(1, nSel));
mT = trials(:, ones(1, nInd))';
[NewChrIx, ~] = find(mT < mF & [zeros(1, nSel); mF(1:nInd-1, :)] <= mT);
% Shuffle selected Individuals
[~, shuf] = sort(rand(nSel, 1));
% NewChrIx is a column vector containing the indexes of the selected
% individuals relative to the original population, shuffeld
NewChrIx = NewChrIx(shuf);
% Relate to indexes of individuals in population
NewChrIx = iFitnessesSpecies(NewChrIx);

end

%% crossover()
function NewIndividual = crossover(...
    Population, NewChrIx, countIndividualsSpecies, Crossover, ...
    matExistingAndPropagSpecies)

% Select Parent1
Parent1 = Population(NewChrIx(2*countIndividualsSpecies-1));

% Select Parent2
foundParent2 = 0;
% Select Parent2 from other species (can only be done if there
% is more than one species in old population)
condition1 = rand() < Crossover.probabilityInterspecies;
condition2 = size(matExistingAndPropagSpecies,2) > 1;
if condition1 && condition2
    while foundParent2 == 0
        [~, iParent2] = max(rand(1,size(Population,2)));
        Parent2 = Population(iParent2);
        % Check if Parent2.species is not species 0 (deleted
        % individual) or species of Parent1
        foundParent2 = ...
            ((Parent2.species ~= 0) & (Parent2.species ~= Parent1.species));
    end
    % Set fitnesses to same to ensure that disjoint and excess
    % genes are inherited fully from both parents (tip from ken)
    Parent2.fitness = Parent1.fitness;
else
    % OK we take Parent2 from same species as Parent1
    Parent2 = Population(NewChrIx(2*countIndividualsSpecies));
end

% Inherit nodes from both parents
NewIndividual.nodeGenes = [];
nGenesP1 = size(Parent1.nodeGenes,2);
nGenesP2 = size(Parent2.nodeGenes,2);
matrixNodeLineup   = [...
    [Parent1.nodeGenes(1,:);1:nGenesP1;zeros(1,nGenesP1)],...
    [Parent2.nodeGenes(1,:);zeros(1,nGenesP2);1:nGenesP2]];
[~, sortNodeVec] = sort(matrixNodeLineup(1,:));
matrixNodeLineup = matrixNodeLineup(:,sortNodeVec);
nodeNumber       = 0;
for iNodeSort = 1:size(matrixNodeLineup,2)
    if nodeNumber ~= matrixNodeLineup(1,iNodeSort)
        if matrixNodeLineup(2,iNodeSort) > 0
            NewIndividual.nodeGenes = ...
                [NewIndividual.nodeGenes,...
                Parent1.nodeGenes(:,matrixNodeLineup(2,iNodeSort))];
        else
            NewIndividual.nodeGenes = ...
                [NewIndividual.nodeGenes,...
                Parent2.nodeGenes(:,matrixNodeLineup(3,iNodeSort))];
        end
        nodeNumber = matrixNodeLineup(1,iNodeSort);
    end
end

% Crossover connection genes
% First do lineup of connection genes
nConnGenesP1 = size(Parent1.connectionGenes,2);
nConnGenesP2 = size(Parent2.connectionGenes,2);
matrixLineup = [...
    [Parent1.connectionGenes(1,:);1:nConnGenesP1;zeros(1,nConnGenesP1)],...
    [Parent2.connectionGenes(1,:);zeros(1,nConnGenesP2);1:nConnGenesP2]];
[~, sortVec] = sort(matrixLineup(1,:));
matrixLineup = matrixLineup(:,sortVec);
finalMatrixLineup = [];
innovationNumber = 0;
for iSort = 1:size(matrixLineup,2)
    if innovationNumber ~= matrixLineup(1,iSort)
        finalMatrixLineup = [finalMatrixLineup, matrixLineup(:,iSort)];
        innovationNumber = matrixLineup(1,iSort);
    else
        finalMatrixLineup(2:3,size(finalMatrixLineup,2)) = ...
            finalMatrixLineup(2:3,size(finalMatrixLineup,2)) + ...
            matrixLineup(2:3,iSort);
    end
end

% OK Connection Genes are lined up, start with crossover
NewIndividual.connectionGenes = [];
for iLineup = 1:size(finalMatrixLineup,2)
    % Check for matching genes, do crossover
    if (finalMatrixLineup(2,iLineup) > 0) && (finalMatrixLineup(3,iLineup) > 0)
        % Random crossover for matching genes
        if rand() < 0.5
            NewIndividual.connectionGenes = ...
                [NewIndividual.connectionGenes,...
                Parent1.connectionGenes(:,finalMatrixLineup(2,iLineup))];
        else
            NewIndividual.connectionGenes = ...
                [NewIndividual.connectionGenes,...
                Parent2.connectionGenes(:,finalMatrixLineup(3,iLineup))];
        end
        % Weight averaging for offspring, otherwise the randomly inherited
        % weights are left undisturbed
        if rand() > Crossover.probabilityMultipoint
            NewIndividual.connectionGenes(4,size(NewIndividual.connectionGenes,2)) = ...
                (...
                Parent1.connectionGenes(4,finalMatrixLineup(2,iLineup)) + ...
                Parent2.connectionGenes(4,finalMatrixLineup(3,iLineup))...
                )/2;
        end
    end
    % Test if there exist further connection genes from
    % iLineup+1 to end of finalMatrixLineup for Parent1 (to
    % detect excess)
    parent1Flag = sum(finalMatrixLineup(2,iLineup+1:size(finalMatrixLineup,2)));
    % Test if there exist further connection genes from
    % iLineup+1 to end of finalMatrixLineup for Parent1 (to
    % detect excess)
    parent2Flag = sum(finalMatrixLineup(3,iLineup+1:size(finalMatrixLineup,2)));
    % Two cases to check (excess is taken care of in the disjoint gene checks)
    % Disjoint Parent1
    if (finalMatrixLineup(2,iLineup) > 0) && (finalMatrixLineup(3,iLineup) == 0)
        if Parent1.fitness >= Parent2.fitness
            NewIndividual.connectionGenes = ...
                [NewIndividual.connectionGenes,...
                Parent1.connectionGenes(:,finalMatrixLineup(2,iLineup))];
        end
    end
    % Disjoint Parent2
    if (finalMatrixLineup(2,iLineup) == 0) && (finalMatrixLineup(3,iLineup) > 0)
        if Parent2.fitness >= Parent1.fitness
            NewIndividual.connectionGenes = ...
                [NewIndividual.connectionGenes,...
                Parent2.connectionGenes(:,finalMatrixLineup(3,iLineup))];
        end
    end
end

% Has no impact on algorithm, only required for assignment to new population
NewIndividual.fitness = 0;

% Will be species hint for speciation
NewIndividual.species = Parent1.species;

end

%% node_culling()
function NewIndividual = node_culling(NewIndividual)

% Hidden nodes culling (remove any hidden nodes where there is no
% corresponding connection gene in the new individual)
connectedNodes = [];
for iNodeCulling = 1:size(NewIndividual.nodeGenes,2)
    nodeConnectedFlag = ...
        sum(NewIndividual.connectionGenes(2,:) == NewIndividual.nodeGenes(1,iNodeCulling)) + ...
        sum(NewIndividual.connectionGenes(3,:) == NewIndividual.nodeGenes(1,iNodeCulling));
    if (nodeConnectedFlag > 0) || (NewIndividual.nodeGenes(2,iNodeCulling) ~= 3)
        connectedNodes = [connectedNodes,NewIndividual.nodeGenes(:,iNodeCulling)];
    end
end
NewIndividual.nodeGenes = connectedNodes;

end

%% enable_gene_mutation()
function NewIndividual = enable_gene_mutation(NewIndividual, Mutation)

% Disabled Genes Mutation
% Run through all connection genes in a NewIndividual, find disabled
% connection genes, enable again with Crossover.probabilityGeneReenabled
% probability
for iConnectionGene = 1:size(NewIndividual.connectionGenes,2)
    condition1 = NewIndividual.connectionGenes(5,iConnectionGene) == 0;
    condition2 = rand() < Mutation.probabilityGeneReenabled;
    if condition1 && condition2
        NewIndividual.connectionGenes(5,iConnectionGene) = 1;
    end
end

end

%% weight_mutation()
function NewIndividual = weight_mutation(NewIndividual, Mutation)

% Weight Mutation
% Run through all connection genes in a NewIndividual, decide on
% mutating or not
for iConnectionGene = 1:size(NewIndividual.connectionGenes,2)
    % *iConnectionGene/size(NewIndividual.connectionGenes,2)
    % linearly biased towards higher probability of mutation at end of
    % connection genes
    if rand() < Mutation.probabilityMutateWeight
        NewIndividual.connectionGenes(4,iConnectionGene) = ...
            NewIndividual.connectionGenes(4,iConnectionGene) + ...
            Mutation.weightRange*(rand()-0.5);
    end
    % Weight of connection genes of NewIndividual
    weightConnGenes = NewIndividual.connectionGenes(4,iConnectionGene);
    % Weight capping
    NewIndividual.connectionGenes(4,iConnectionGene) = ...
        weightConnGenes*(abs(weightConnGenes) <= Mutation.weightCap) + ...
        (sign(weightConnGenes)*Mutation.weightCap)*(abs(weightConnGenes) > Mutation.weightCap);
end

end

%% add_connection_mutation()
function [NewIndividual, innovationRecord] = ...
    add_connection_mutation(NewIndividual, innovationRecord, Mutation,...
    iGeneration, flag1)

% Add Connection Mutation
flagRecurrencyEnabled = rand() < Mutation.probabilityRecurrency;

% Connections can run from every node
vectorPossibleConnectFromNodes = NewIndividual.nodeGenes(1,:);

% Connections can only run into hidden and output nodes
vectorPossibleConnectToNodes = ...
    NewIndividual.nodeGenes(...
    1,find((NewIndividual.nodeGenes(2,:) == 2) + (NewIndividual.nodeGenes(2,:) == 3))...
    );
numberPossibleConnection = ...
    length(vectorPossibleConnectFromNodes)*length(vectorPossibleConnectToNodes) - ...
    size(NewIndividual.connectionGenes,2);

% Check if new connections can be added to genes (if there are any
% possible connections which are not already existing in genes of
% new individual)
condition1 = rand() < Mutation.probabilityAddConnection;
condition2 = numberPossibleConnection > 0;
condition3 = flag1 == 0;
if condition1 && condition2 && condition3
    % First build matrix containing all possible new connection for
    % nodegene of new individual
    newConnectionMatrix = [];
    for iConnectFrom = 1:length(vectorPossibleConnectFromNodes)
        for iConnectTo = 1:length(vectorPossibleConnectToNodes)
            possibleConnection = ...
                [...
                vectorPossibleConnectFromNodes(iConnectFrom);
                vectorPossibleConnectToNodes(iConnectTo)];
            % Check if proposed connection is not already contained in gene
            condition4 = sum(...
                (NewIndividual.connectionGenes(2,:) == possibleConnection(1)).* ...
                (NewIndividual.connectionGenes(3,:) == possibleConnection(2)) ...
                ) == 0;
            if condition4
                newConnectionMatrix = [newConnectionMatrix, possibleConnection];
            end
        end
    end
    % Shuffle possible new connections randomly
    [~, shuffle] = sort(rand(1,size(newConnectionMatrix,2)));
    newConnectionMatrix = newConnectionMatrix(:,shuffle);
    
    iNewConnection = 0;
    flagConnectionOk = 0;
    % Check if connection is ok. (meaning either non-recurrent or
    % recurrent and flagRecurrencyEnabled set to 1) if no
    % connection is found which is ok, no connection will be added
    % to connection genes of new individual
    while (flagConnectionOk == 0) && (iNewConnection < size(newConnectionMatrix,2))
        iNewConnection = iNewConnection + 1;
        newConnection = newConnectionMatrix(:,iNewConnection);
        % Test new connection if it is recurrent (i.e. at least one
        % of the possibles path starting from connect_to node in
        % the network leads back to the connect_from node
        flagRecurrent = 0;
        % Trivial recurrency
        if newConnection(1) == newConnection(2)
            flagRecurrent = 1;
        end
        nodesCurrentLevel = newConnection(2);
        depth = 0;
        while (flagRecurrent == 0) && ...
                (depth < size(NewIndividual.connectionGenes,2)) && ...
                ~isempty(nodesCurrentLevel)
            depth = depth + 1;
            nodesNextLevel = [];
            for indexCheck = 1:size(nodesCurrentLevel);
                indexTemp = find(...
                    NewIndividual.connectionGenes(2,:) == ...
                    nodesCurrentLevel(indexCheck) ...
                    );
                nodesNextLevel = ...
                    [nodesNextLevel,...
                    NewIndividual.connectionGenes(3,indexTemp)];
            end
            if sum(nodesNextLevel(:) == newConnection(1)) > 0
                flagRecurrent = 1;
            end
            nodesCurrentLevel = nodesNextLevel;
        end
        if flagRecurrent == 0
            flagConnectionOk = 1;
        elseif flagRecurrencyEnabled
            flagConnectionOk = 1;
        end
    end
    
    % Now we test if it is a true innovation (i.e. hasn't already
    % happened in this generation) we can only do this if a valid
    % new connection has been found
    if flagConnectionOk
        % Set flag signifying new innovation (connection not
        % contained in innovationRecord of this generation)
        iAlreadyHappened = find(...
            (innovationRecord(5,:) == iGeneration).* ...
            (innovationRecord(2,:) == newConnection(1)).* ...
            (innovationRecord(3,:) == newConnection(2)) ...
            );
        newInnovation = not(sum(iAlreadyHappened));
        if newInnovation == 1 % OK is new innovation
            % Update the new connection with its innovation number
            newConnection = [max(innovationRecord(1,:))+1;newConnection];
            % Update connection_genes. Weight is uniformly distributed between
            % [-1,+1]
            NewIndividual.connectionGenes = [...
                NewIndividual.connectionGenes,...
                [newConnection;rand()*2-1;1]];
            % Update innovationRecord
            innovationRecord = [innovationRecord,[newConnection;0;iGeneration]];
        else % Connection gene already exists in innovationRecord of this generation
            % Update connectionGenes. Weight is uniformly distributed between
            % [-1,+1]
            NewIndividual.connectionGenes = [...
                NewIndividual.connectionGenes,...
                [innovationRecord(1:3,iAlreadyHappened);rand()*2-1;1]];
        end
    end
    
end

end

%% add_node_mutation()
function [NewIndividual, innovationRecord] = ...
    add_node_mutation(innovationRecord, iGeneration, NewIndividual)

% Add (Insert) Node Mutation
newInnovation = 0;
% Highest innovation number from last generation (to ensure
% that only connections from from last generation or older are
% chosen for add node mutation, otherwise a new connection
% added in the last mutation might instantly be disabled)
maxOldInnovationNumber = max(...
    (innovationRecord(5,:) < iGeneration).*innovationRecord(1,:) ...
    );
% Compute vector of connections into which a new node could be
% inserted and their positions in the connectionGene matrix.
% This vector is composed of all nondisabled connections which
% stem at least from the last generation or older
indexTemp = find(...
    (NewIndividual.connectionGenes(5,:) == 1) & ...
    (NewIndividual.connectionGenes(1,:) <= maxOldInnovationNumber) ...
    );
vectorPossibleConnections = [...
    NewIndividual.connectionGenes(2:3,indexTemp);
    indexTemp];
insertNodeConnection = ...
    vectorPossibleConnections(:,round(rand()*size(vectorPossibleConnections,2) + 0.5));
% Set provisionally to 1, will be checked
newInnovation = 1;
% Beginning of check innovation record to test for real
% innovation. existInnovation contains vector of index of
% elements in innovation record which fulfil three things:
% current generation, add node mutation and same connect from
% as current innovation
existInnovation = find((innovationRecord(5,:) == ...
    iGeneration).*...
    (innovationRecord(4,:) > 0).*...
    (innovationRecord(2,:) == insertNodeConnection(1)));
% If these are fulfilled, we have to test for connect_to node
% to see if innovation really is the same
if sum(existInnovation) > 0
    for indexCheck = 1:length(existInnovation)
        if innovationRecord(3,existInnovation(indexCheck)+1) == insertNodeConnection(2)
            newInnovation = 0;
            iAlreadyExistentThisGeneration = existInnovation(indexCheck);
        end
    end
end
if newInnovation == 1 % OK is true innovation for current generation
    % Update node genes. By [newNodeNumber;3;0;0], we see that only hidden nodes
    % can be added by add_node_mutation.
    newNodeNumber = max(innovationRecord(4,:)) + 1;
    NewIndividual.nodeGenes = [NewIndividual.nodeGenes,[newNodeNumber;3;0;0]];
    % Disable old connection gene
    NewIndividual.connectionGenes(5,insertNodeConnection(3)) = 0;
    % Update connectionGenes. newConnections is a 5 rows by two columns matrix
    % where the first column is the new connection TO the new node whereas the
    % second column is the connection FROM the new node
    newConnections = [...
        [max(innovationRecord(1,:)) + 1;
        insertNodeConnection(1);
        newNodeNumber;
        1;
        1],...
        [max(innovationRecord(1,:)) + 2;
        newNodeNumber;
        insertNodeConnection(2);
        NewIndividual.connectionGenes(4,insertNodeConnection(3));
        1]];
    % Extend connection genes by the two new connections
    NewIndividual.connectionGenes = [NewIndividual.connectionGenes, newConnections];
    % Update innovationRecord
    innovationRecord = [...
        innovationRecord,...
        [newConnections(1:3,:);newNodeNumber,0;iGeneration,iGeneration]];
else % No new innovation, has already happened at least once in this generation
    % Update node genes
    nodeNumber = innovationRecord(4,iAlreadyExistentThisGeneration);
    NewIndividual.nodeGenes = [NewIndividual.nodeGenes,[nodeNumber;3;0;0]];
    % Update connectionGenes
    % Disable old connection gene
    NewIndividual.connectionGenes(5,insertNodeConnection(3)) = 0;
    newConnections = [...
        innovationRecord(1:3,iAlreadyExistentThisGeneration:iAlreadyExistentThisGeneration+1);
        1,NewIndividual.connectionGenes(4,insertNodeConnection(3));1,1];
    % Length of the connection genes of current NewIndividual
    lengthConGen = size(NewIndividual.connectionGenes,2);
    % Check if there was an add_connection_mutation to current
    % NewIndividual which has a higher innovation number than
    % current add_node_mutation
    if NewIndividual.connectionGenes(1,lengthConGen) > newConnections(1,2)
        NewIndividual.connectionGenes = [...
            NewIndividual.connectionGenes(:,1:lengthConGen-1), ...
            newConnections, ...
            NewIndividual.connectionGenes(:,lengthConGen)];
    else
        NewIndividual.connectionGenes = [...
            NewIndividual.connectionGenes,...
            newConnections];
    end
end

end

%% speciation()
function [NewIndividual, SpeciesRecord, populationRef] = ...
    speciation(NewIndividual, SpeciesRecord, populationRef, Speciation)

% This function performs speciation to assign the new individual into an
% existing or new specie

% Loop through comparison vector
speciesAssigned = 0;
iPopulationRef = 0;
while (speciesAssigned == 0) && (iPopulationRef < size(populationRef,2))
    % Extract referenceIndividual from reference population
    iPopulationRef = iPopulationRef + 1;
    referenceIndividual = populationRef(iPopulationRef);
    % Run through both connection genes, compute disjoint, excess,
    % and average weight difference.
    % Maximum number of genes between new individual and refenrence individual
    maxNumGenes = max(...
        [size(NewIndividual.connectionGenes,2),...
        size(referenceIndividual.connectionGenes,2)]...
        );
    % Maximum innovation number between new individual and refenrence individual
    maxNumInnovation = max(...
        [NewIndividual.connectionGenes(1,:),...
        referenceIndividual.connectionGenes(1,:)]...
        );
    % New innovation vector
    vectorInnovationNew = [...
        zeros(1,max(NewIndividual.connectionGenes(1,:))),...
        ones(1,maxNumInnovation - max(NewIndividual.connectionGenes(1,:)))];
    vectorInnovationNew(NewIndividual.connectionGenes(1,:)) = 2;
    % New weight vector
    vectorWeightNew = zeros(1,maxNumInnovation);
    vectorWeightNew(NewIndividual.connectionGenes(1,:)) = ...
        NewIndividual.connectionGenes(4,:);
    % Innovation vector for refenrence individual
    vectorInnovationRef = [...
        4*ones(1,max(referenceIndividual.connectionGenes(1,:))),...
        8*ones(1,maxNumInnovation - max(referenceIndividual.connectionGenes(1,:)))];
    vectorInnovationRef(referenceIndividual.connectionGenes(1,:)) = 16;
    % Weight vector for refenrence individual
    vectorWeightRef = zeros(1,maxNumInnovation);
    vectorWeightRef(referenceIndividual.connectionGenes(1,:)) = ...
        referenceIndividual.connectionGenes(4,:);
    % Lineup vector between new individual and reference individual
    vectorLineup = vectorInnovationNew + vectorInnovationRef;
    % Number of excess genes between new and refrence individuals
    excess = sum(vectorLineup == 10) + sum(vectorLineup == 17);
    % Number of disjoint genes between new and refrence individuals
    disjoint = sum(vectorLineup == 6) + sum(vectorLineup == 16);
    % Matching genes
    vectorMatching = find(vectorLineup == 18);
    % Average weight difference
    averageWeightDifference = ...
        sum(...
        abs(vectorWeightNew(vectorMatching) - vectorWeightRef(vectorMatching))...
        )/length(vectorMatching);
    % Compute compatbility distance between new and refenrence individuals
    maxNumGenes = 1;
    distance = Speciation.c1*excess/maxNumGenes + ...
        Speciation.c2*disjoint/maxNumGenes + ...
        Speciation.c3*averageWeightDifference;
    if distance < Speciation.threshold
        % Assign individual to same species as current reference individual
        NewIndividual.species = referenceIndividual.species;
        % Set flag indicating NewIndividual has been assigned to species
        speciesAssigned = 1;
    end
end

% Not compatible with any? well, then create new species
if speciesAssigned == 0
    newSpeciesId = size(SpeciesRecord,2) + 1;
    % Assign individual to new species
    NewIndividual.species = newSpeciesId;
    % Update SpeciesRecord
    SpeciesRecord(newSpeciesId).id = newSpeciesId;
    SpeciesRecord(newSpeciesId).nIndividuals = 1;
    SpeciesRecord(newSpeciesId).generationRecord = [];
    % Update population reference
    populationRef(size(populationRef,2)+1) = NewIndividual;
end

end

%% compute_stats()
function [nAverageNonDisabledConnections, nAverageHiddenNodes, maxFitness,...
    maxOverallFitness, nAverageDisabledConnections, ...
    nSpeciesArray, meanOverallFitness, meanRecurrentConnection, nIndPerSpecie] = ...
    compute_stats(Population, populationSize, iGeneration, ...
    SpeciesRecord, nAverageNonDisabledConnections, nAverageHiddenNodes, ...
    maxOverallFitness, nAverageDisabledConnections, ...
    nSpeciesArray, meanOverallFitness, meanRecurrentConnection, nIndPerSpecie)

% Number of individuals
nIndividuals = size(Population,2);

% Compute the average number of enabled connections and hidden nodes in current
% population
a = 0; % Sum of enabled connection in current population
b = 0; % Sum of hidden nodes in current population
c = 0; % Sum of disabled connection in current population
for iIndividual = 1:nIndividuals
    a = a + sum(Population(iIndividual).connectionGenes(5,:) == 1);
    b = b + sum(Population(iIndividual).nodeGenes(2,:) == 3);
    c = c + sum(Population(iIndividual).connectionGenes(5,:) == 0);
end
nAverageNonDisabledConnections = [...
    nAverageNonDisabledConnections,...
    [a/populationSize;iGeneration]];
nAverageDisabledConnections = [...
    nAverageDisabledConnections,...
    [c/populationSize;iGeneration]];
nAverageHiddenNodes = [...
    nAverageHiddenNodes,...
    [b/populationSize;iGeneration]];

% Number of species in current generation
nSpecies = size(SpeciesRecord,2);
% nSpeciesArray = [nSpeciesArray,nSpecies];
nSpeciesArray = [nSpeciesArray,numel(unique([Population(:).species]))];

% Average number of individuals per specie
uniqueSpecieId = unique([Population(:).species]);
nIndPerSpecieTmp = zeros(1,numel(uniqueSpecieId));
for k = 1:numel(uniqueSpecieId)
    nIndPerSpecieTmp(k) = sum([Population(:).species] == uniqueSpecieId(k));
end
nIndPerSpecie = [nIndPerSpecie,mean(nIndPerSpecieTmp)];

% We build the fitCurrGenSpecie matrix with 3 rows and nSpecies columns. The
% first row is the index of the current generation of the specie. The second row
% is mean raw fitness of the current generation for the specie and the third row
% is the max raw fitness of the current generation for the specie.
fitCurrGenSpecie = [];
for iSpecie = 1:nSpecies
    % Number of generations in current specie
    nGenerations = size(SpeciesRecord(iSpecie).generationRecord,2);
    fitCurrGenSpecie = [...
        fitCurrGenSpecie,...
        SpeciesRecord(iSpecie).generationRecord(1:3,nGenerations)];
end

% Maximal fitness for current generation
maxFitness = max(fitCurrGenSpecie(3,:).*(fitCurrGenSpecie(1,:) == iGeneration));
fprintf('Max fitness is %8.6f\n', maxFitness);

% Maximum overall fitness. It is a matrix with 2 rows and nGenerations columns.
% The first row is the maximum overall fitness of the current generation and the
% second row is the generation number.
maxOverallFitness = [...
    maxOverallFitness,...
    [maxFitness;iGeneration]];

% Mean fitness for current generation
meanFitness = mean([Population(:).fitness]);
% meanFitness = mean(fitCurrGenSpecie(3,:).*(fitCurrGenSpecie(1,:) == iGeneration));
meanOverallFitness = [...
    meanOverallFitness,...
    [meanFitness;iGeneration]];

% Number of recurrent connections
nRecurrentConnections = ...
    find_recurrent_connection(Population);
meanRecurrentConnection = [meanRecurrentConnection,mean(nRecurrentConnections)];

end

%% solution_check()
function flagSolution = solution_check(maxFitness, fitnessLimit, flagSolution)

% Check for solution
if maxFitness > fitnessLimit
    flagSolution = 1;
end

end

%% visualization()
function visualization(figHandle, AxHandle, nAverageNonDisabledConnections, ...
    nAverageHiddenNodes, maxOverallFitness, nAverageDisabledConnections, ...
    nSpeciesArray, meanOverallFitness, meanRecurrentConnection, nIndPerSpecie)

nGenerations = numel(nSpeciesArray);

lineSpec = {'LineWidth',1};
fontSize = 8;
nYTick = 9;
if 1 <= nGenerations && nGenerations < 6
    XTick1 = unique(round(linspace(1,5,5)));
    XTick2 = unique(round(linspace(1,5,5)));
elseif 6 <= nGenerations && nGenerations < 10
    XTick1 = unique(round(linspace(1,nGenerations,9)));
    XTick2 = unique(round(linspace(1,nGenerations,6)));
elseif 10 <= nGenerations && nGenerations < Inf
    XTick1 = unique(round(linspace(1,nGenerations,15)));
    XTick2 = unique(round(linspace(1,nGenerations,6)));
end
XTickLabel1 = cellfun(@num2str,num2cell(XTick1(:)),'uniformoutput',false);
XTickLabel2 = cellfun(@num2str,num2cell(XTick2(:)),'uniformoutput',false);

% Maximum fitness
set(figHandle, 'CurrentAxes', AxHandle.subPlot1);
hold(AxHandle.subPlot1,'on');
plot(AxHandle.subPlot1, maxOverallFitness(2,:), maxOverallFitness(1,:), lineSpec{:});
% AxHandle.subPlot1.XTick = XTick1;
% AxHandle.subPlot1.XTickLabel = XTickLabel1;
% AxHandle.subPlot1.YTick = linspace(min(AxHandle.subPlot1.YLim),max(AxHandle.subPlot1.YLim),nYTick);
xlabel(AxHandle.subPlot1,'Generation','FontSize', fontSize);
ylabel(AxHandle.subPlot1,'Maximum fitness','FontSize', fontSize);
grid(AxHandle.subPlot1,'on');
hold(AxHandle.subPlot1,'off');

% Average fitness
set(figHandle, 'CurrentAxes', AxHandle.subPlot2);
hold(AxHandle.subPlot2,'on');
plot(AxHandle.subPlot2, meanOverallFitness(2,:), meanOverallFitness(1,:), lineSpec{:});
% AxHandle.subPlot2.XTick = XTick2;
% AxHandle.subPlot2.XTickLabel = XTickLabel2;
% AxHandle.subPlot2.YTick = linspace(min(AxHandle.subPlot2.YLim),max(AxHandle.subPlot2.YLim),nYTick);
xlabel(AxHandle.subPlot2,'Generation','FontSize', fontSize);
ylabel(AxHandle.subPlot2,'Average fitness','FontSize', fontSize);
grid(AxHandle.subPlot2,'on');
hold(AxHandle.subPlot2,'off');

% Average number of hidden nodes
set(figHandle, 'CurrentAxes', AxHandle.subPlot3);
hold(AxHandle.subPlot3,'on');
plot(AxHandle.subPlot3, nAverageHiddenNodes(2,:), nAverageHiddenNodes(1,:), lineSpec{:});
% AxHandle.subPlot3.XTick = XTick2;
% AxHandle.subPlot3.XTickLabel = XTickLabel2;
% AxHandle.subPlot3.YTick = linspace(min(AxHandle.subPlot3.YLim),max(AxHandle.subPlot3.YLim),nYTick);
xlabel(AxHandle.subPlot3,'Generation','FontSize', fontSize);
ylabel(AxHandle.subPlot3,{'Avg. number of';'hidden nodes'},'FontSize', fontSize);
grid(AxHandle.subPlot3,'on');
hold(AxHandle.subPlot3,'off');

% Average number of non disabled connections
set(figHandle, 'CurrentAxes', AxHandle.subPlot4);
hold(AxHandle.subPlot4,'on');
plot(AxHandle.subPlot4, nAverageNonDisabledConnections(2,:),...
    nAverageNonDisabledConnections(1,:), lineSpec{:});
% AxHandle.subPlot4.XTick = XTick2;
% AxHandle.subPlot4.XTickLabel = XTickLabel2;
% AxHandle.subPlot4.YTick = linspace(min(AxHandle.subPlot4.YLim),max(AxHandle.subPlot4.YLim),nYTick);
xlabel(AxHandle.subPlot4,'Generation','FontSize', fontSize);
ylabel(AxHandle.subPlot4,{'Avg. number of';'non disabled connections'}, ...
    'FontSize', fontSize);
grid(AxHandle.subPlot4,'on');
hold(AxHandle.subPlot4,'off');

% Average number of disabled connections
set(figHandle, 'CurrentAxes', AxHandle.subPlot5);
hold(AxHandle.subPlot5,'on');
plot(AxHandle.subPlot5, nAverageDisabledConnections(2,:), ...
    nAverageDisabledConnections(1,:), lineSpec{:});
% AxHandle.subPlot5.XTick = XTick2;
% AxHandle.subPlot5.XTickLabel = XTickLabel2;
% AxHandle.subPlot5.YTick = linspace(min(AxHandle.subPlot5.YLim),max(AxHandle.subPlot5.YLim),nYTick);
xlabel(AxHandle.subPlot5,'Generation','FontSize', fontSize);
ylabel(AxHandle.subPlot5,{'Avg. number of';'disabled connections'}, ...
    'FontSize', fontSize);
grid(AxHandle.subPlot5,'on');
hold(AxHandle.subPlot5,'off');

% Average number of recurrent connections
set(figHandle, 'CurrentAxes', AxHandle.subPlot6);
hold(AxHandle.subPlot6,'on');
plot(AxHandle.subPlot6, nAverageDisabledConnections(2,:), ...
    meanRecurrentConnection, lineSpec{:});
% AxHandle.subPlot6.XTick = XTick2;
% AxHandle.subPlot6.XTickLabel = XTickLabel2;
% AxHandle.subPlot6.YTick = linspace(min(AxHandle.subPlot6.YLim),max(AxHandle.subPlot6.YLim),nYTick);
xlabel(AxHandle.subPlot6,'Generation','FontSize', fontSize);
ylabel(AxHandle.subPlot6,{'Avg. number of';'recurrent connections'}, ...
    'FontSize', fontSize);
grid(AxHandle.subPlot6,'on');
hold(AxHandle.subPlot6,'off');

% Number of species
set(figHandle, 'CurrentAxes', AxHandle.subPlot7);
hold(AxHandle.subPlot7,'on');
plot(AxHandle.subPlot7, nAverageDisabledConnections(2,:), nSpeciesArray, lineSpec{:});
% AxHandle.subPlot7.XTick = XTick2;
% AxHandle.subPlot7.XTickLabel = XTickLabel2;
% AxHandle.subPlot7.YTick = linspace(min(AxHandle.subPlot7.YLim),max(AxHandle.subPlot7.YLim),nYTick);
xlabel(AxHandle.subPlot7,'Generation','FontSize', fontSize);
ylabel(AxHandle.subPlot7,{'Number of';'species'}, ...
    'FontSize', fontSize);
grid(AxHandle.subPlot7,'on');
hold(AxHandle.subPlot7,'off');

% Average number of individuals per specie
% if numel(nIndPerSpecie) >= 3 && numel(unique(nIndPerSpecie)) >= 3
%     errorBreak = false;
%     try
%         firstMaxVal     = max(nIndPerSpecie);
%         secondMaxVal    = max(nIndPerSpecie(nIndPerSpecie ~= firstMaxVal));
%         thirdMaxVal     = max(nIndPerSpecie(nIndPerSpecie ~= firstMaxVal & nIndPerSpecie ~= secondMaxVal));
%         if abs(firstMaxVal-secondMaxVal) < 10
%             lowerBreak = round(thirdMaxVal+10,0);
%             upperBreak = round(secondMaxVal-10,0);
%         else
%             lowerBreak = round(secondMaxVal+10,0);
%             upperBreak = round(firstMaxVal-10,0);
%         end
%     catch
%         errorBreak = true;
%     end
% end
set(figHandle, 'CurrentAxes', AxHandle.subPlot8);
hold(AxHandle.subPlot8,'on');
% if numel(nIndPerSpecie) >= 3 && numel(unique(nIndPerSpecie)) >= 3
%     if ~errorBreak
%         try
%             %breakInfo = breakyaxis(AxHandle.subPlot8,[lowerBreak,upperBreak]);
%         catch
%
%         end
%     end
% end
plot(AxHandle.subPlot8, nAverageDisabledConnections(2,:), nIndPerSpecie, lineSpec{:});
% AxHandle.subPlot8.XTick = XTick2;
% AxHandle.subPlot8.XTickLabel = XTickLabel2;
% AxHandle.subPlot8.YTick = linspace(min(AxHandle.subPlot8.YLim),max(AxHandle.subPlot8.YLim),nYTick);
xlabel(AxHandle.subPlot8,'Generation','FontSize', fontSize);
ylabel(AxHandle.subPlot8,{'Avg. number of';'individuals per specie'}, ...
    'FontSize', fontSize);
grid(AxHandle.subPlot8,'on');
hold(AxHandle.subPlot8,'off');

drawnow;

hold off

end

%% print_generation()
function print_generation(iGeneration)

fprintf([...
    '\n\n---------------------------------------------------------------------------\n',...
    '|                             Generation %4.f                             |\n',...
    '---------------------------------------------------------------------------\n\n'...
    ],iGeneration);

end

%% shrink_mutation()
function Mutation = shrink_mutation(Mutation, iGeneration)

% Add node
if Mutation.probabilityAddNodeScheduling(1) ~= -1
    for iGen = 1:numel(Mutation.probabilityAddNodeScheduling(1,:))
        if iGeneration <= Mutation.probabilityAddNodeScheduling(1,iGen)
            Mutation.probabilityAddNode = Mutation.probabilityAddNodeScheduling(2,iGen);
            break;
        end
    end
else
    Mutation.probabilityAddNode = max([...
        Mutation.probabilityAddNode - Mutation.probabilityAddNodeShrinkRate,...
        Mutation.probabilityAddNodeMin]);
end

% Add connection
if Mutation.probabilityAddConnectionScheduling(1) ~= -1
    for iGen = 1:numel(Mutation.probabilityAddConnectionScheduling(1,:))
        if iGeneration <= Mutation.probabilityAddConnectionScheduling(1,iGen)
            Mutation.probabilityAddConnection = Mutation.probabilityAddConnectionScheduling(2,iGen);
            break;
        end
    end
else
    Mutation.probabilityAddConnection = max([...
        Mutation.probabilityAddConnection - Mutation.probabilityAddConnectionShrinkRate,...
        Mutation.probabilityAddConnectionMin]);
end

% Recurrent connection
Mutation.probabilityRecurrency = max([...
    Mutation.probabilityRecurrency - Mutation.probabilityRecurrencyShrinkRate,...
    Mutation.probabilityRecurrencyMin]);

% Weight mutation
Mutation.probabilityMutateWeight = max([...
    Mutation.probabilityMutateWeight - Mutation.probabilityMutateWeightShrinkRate,...
    Mutation.probabilityMutateWeightMin]);

% Display new mutation probabilities
fprintf([...
    'Mutation probabilities update:\n',...
    'Mutation.probabilityAddNode       = %4.5f\n',...
    'Mutation.probabilityAddConnection = %4.5f\n',...
    'Mutation.probabilityRecurrency    = %4.5f\n',...
    'Mutation.probabilityMutateWeight  = %4.5f\n\n',...
    ],...
    Mutation.probabilityAddNode,...
    Mutation.probabilityAddConnection,...
    Mutation.probabilityRecurrency,...
    Mutation.probabilityMutateWeight...
    );

end
