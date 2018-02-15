clear all
close all
clc

% Uncomment for repetability between runs
% rng(6541651)
rng('shuffle');

% Constants
DEG_2_RAD = pi/180;
RAD_2_DEG = 180/pi;
MM_2_M = 1/1000;
M_2_MM = 1000/1;

% Date and time
dateAndTime = datestr(datetime('now','TimeZone','local',...
    'Format','d-MMM-y HH:mm:ss Z'));
dateAndTime = strrep(dateAndTime,':','_');
dateAndTime = strrep(dateAndTime,'-','_');
dateAndTime = strrep(dateAndTime,' ','_');

% Set up diary
if ~exist('diary','dir')
    mkdir('diary');
end
diary(['diary/','design_controller','_',dateAndTime,'.txt']);
fprintf('Diary stored in %s\n\n',...
    ['diary/','design_controller','_',dateAndTime,'.txt']);

% Close all figure handles
figAllHandle = findall(0,'Type','figure');
if ~isempty(figAllHandle)
    for iFigHandle = 1:numel(figAllHandle)
        % figAllHandle(iFigHandle).CloseRequestFcn = @closereq;
        set(0, 'currentfigure', figAllHandle(iFigHandle));
        close('force');
        % close(figAllHandle(iFigHandle));
        delete(figAllHandle(iFigHandle));
    end
    clear figAllHandle;
end

%% Delete paths from MATLAB search path

% Code snippet from https://goo.gl/0OwRMn

% Turn off warnings temporarily
% warning off %#ok<WNOFF>

% Get a list of all files and folders in this folder.
if exist('temp_uav_nn_controller','dir')
    files = dir('temp_uav_nn_controller');
    files(1:2) = [];
end

% Get a logical vector that tells which is a directory.
if exist('temp_uav_nn_controller','dir')
    dirFlags = [files.isdir];
end

% Extract only those that are directories.
if exist('temp_uav_nn_controller','dir')
    subFolders = files(dirFlags);
end

% Remove path
if exist('temp_uav_nn_controller','dir')
    for k = 1:length(subFolders)
        rmpath(['temp_uav_nn_controller/',subFolders(k).name])
    end
end
rehash();

% Turn on warnings again
% warning on %#ok<WNON>

%% Parameters - General

% If set to 1, will load population, generation, innovation_record and
% species_record from neatsave.mat at start of algorithm, if set to 0,
% algorithm will start with initial population, new species record and new
% innovation record, at generation = 1 (default option)
loadFlag = 0;

% If set to 1, will save population, generation, innovation_record and
% species_record to neatsave.mat at every generation (default option)
saveFlag = 1;

% If set to 1, NEAT will start with the gene from geneFile and intialize the
% weight randomnly between [-1,+1]
loadGene = 0;
geneFile = 'neat_data/neatsave_27_Nov_2016_20_27_44_043.mat';

% Maximum number of generations for generational loop
nMaxGeneration = 200; % 200

% Population size
populationSize = 10; % 150, 500, 300

% Input and output nodes numbers
nInputNodes = 25;
nOutputNodes = 4;

% Vector of initially connected input nodes out of complete number of input nodes
% (if you want to start with a subset and let evolution decide which ones are necessary)
% for a fully connected initial population, list all the inputs. For each input
% node listed in initialConnectedInputNodes, it will be connected to every
% output node during population's intialization.
% initialConnectedInputNodes = 1:nInputNodes;
initialConnectedInputNodes = [1,2,3,4,5,6,7,8,9,16,17,18,22,23,24,25];
% initialConnectedInputNodes = 1;

% Fitness limit
fitnessLimit = Inf;

%% Parameters - Speciation

% Speciation parameters as published by Ken Stanley
Speciation.c1 = 1.0;
Speciation.c2 = 1.0;
Speciation.c3 = 0.4; % 0.4
Speciation.threshold = 25.0; % 3 (10.0,12.0,14.0)

%% Parameters - Reproduction - Stagnation + Refocuse

% Threshold to judge if a species is in stagnation (max fitness of species
% varies below threshold) this threshold is of course dependent on your
% fitness function, if you have a fitness function which has a large
% spread, you might want to increase this threshold
Stagnation.threshold = 1e-2;

% If max fitness of species has stayed within Stagnation.threshold in the
% last Stagnation.nGenerations generations, all its fitnesses will be
% reduced to 0, so it will die out
Stagnation.nGenerations = 15; % 15

% Computation is done the following way: the absolute difference between
% the average max fitness of the last Stagnation.nGenerations
% generations and the max fitness of each of these generations is computed
% and compared to Stagnation.threshold.
% if it stays within this threshold for the indicated number of
% generations, the species is eliminated
Refocus.threshold = 1e-2;

% If maximum overall fitness of population doesn't change within threshold
% for this number of generations, only the top two species are allowed to
% reproduce
Refocus.nGenerations = 20; % 20

%% Parameters - Reproduction - Initial

% The percentage of each species which will be eliminated (lowest
% performing individuals). This percentage for eliminating individuals will
% only be used in species which have more individuals than numberForKill
Initial.killPercentage = 0.2;

% Please note that whatever the above settings, the code always ensures
% that at least 2 individuals are kept to be able to cross over, or at
% least the single individual in a species with only one individual
Initial.numberForKill = 5;

% Species which have equal or greater than numberCopy individuals will
% have their best individual copied unchanged into the next generation
Initial.numberCopy = 5;

%% Parameters - Reproduction - Selection

% Selection (ranking and stochastic universal sampling)
% Number between 1.1 and 2.0, determines selective pressure towards most
% fit individual of species
Selection.pressure = 2;

%% Parameters - Reproduction - Crossover

% Percentage governs the way in which new population will be composed from
% old population. exception: species with just one individual can only use
% mutation
Crossover.percentage = 0.8;

% If crossover has been selected, this probability governs the
% intra/interspecies parent composition being used for the
% standard-crossover in which matching connection genes are inherited
% randomly from both parents. In the (1-Crossover.probabilityMultipoint)
% cases, weights of the new connection genes are the mean of the
% corresponding parent genes
Crossover.probabilityInterspecies = 0.001;
Crossover.probabilityMultipoint = 0.6;

%% Parameters - Reproduction - Mutation

% Adding node probability
Mutation.probabilityAddNode                 = 0.60; % 0.03
Mutation.probabilityAddNodeMin              = 0.05; % 0.03
Mutation.probabilityAddNodeShrinkRate       = ...
    (Mutation.probabilityAddNode-Mutation.probabilityAddNodeMin)/60;
if true
    Mutation.probabilityAddNodeScheduling = [...
        0005,  10,  15,  20,  30,  40,  70,  75,  90,  100, Inf;
        0.90,0.50,0.40,0.30,0.20,0.10,0.90,0.85,0.80,0.10,0.05];
else
    Mutation.probabilityAddNodeScheduling = -1;
end

% Adding connection probability
Mutation.probabilityAddConnection           = 0.80; % 0.05
Mutation.probabilityAddConnectionMin        = 0.20;
Mutation.probabilityAddConnectionShrinkRate = ...
    (Mutation.probabilityAddConnection-Mutation.probabilityAddConnectionMin)/60;
if true
    Mutation.probabilityAddConnectionScheduling = [...
        0005,  10,  15,  20,  30,  40,  70,  75,  90,  Inf;
        0.90,0.70,0.60,0.45,0.40,0.30,0.90,0.85,0.80,0.20];
else
    Mutation.probabilityAddNodeScheduling = -1;
end

% If we are in add_connection_mutation, this governs if a recurrent
% connection is allowed. Note: this will only activate if the random
% connection is a recurrent one, otherwise the connection is simply
% accepted. If no possible non-recurrent connections exist for the
% current node genes, then for e.g. a probability of 0.1, 9 times out
% of 10 no connection is added.
Mutation.probabilityRecurrency              = 0.05; % 0.05
Mutation.probabilityRecurrencyMin           = Mutation.probabilityRecurrency;
Mutation.probabilityRecurrencyShrinkRate    = ...
    (Mutation.probabilityRecurrency-Mutation.probabilityRecurrencyMin)/60;
Mutation.probabilityMutateWeight            = 0.9;  % 0.9
Mutation.probabilityMutateWeightMin         = 0.9;
Mutation.probabilityMutateWeightShrinkRate  = ...
    (Mutation.probabilityMutateWeight-Mutation.probabilityMutateWeightMin)/60;

% Weights will be restricted from -Mutation.weightCap to Mutation.weightCap
Mutation.weightCap = 0.5; % 8

% Random distribution with width Mutation.weightRange, centered on 0.
% mutation range of 5 will give random distribution from -2.5 to 2.5
Mutation.weightRange = 0.25; % 5

% Probability of a connection gene being reenabled in offspring if it was
% inherited disabled
Mutation.probabilityGeneReenabled = 0.10; % 0.20

%% Display parameters

if Mutation.probabilityAddNodeScheduling(1) == -1
    mutationAddNode1Str = 'Mutation.probabilityAddNode                 = %4.5f\n';
    mutationAddNode2Str = 'Mutation.probabilityAddNodeShrinkRate       = %4.5f\n';
    mutationAddNode3Str = 'Mutation.probabilityAddNodeMin              = %4.5f\n';
    mutationAddNode1Val = {Mutation.probabilityAddNode};
    mutationAddNode2Val = Mutation.probabilityAddNodeShrinkRate;
    mutationAddNode3Val = Mutation.probabilityAddNodeMin;
else
    mutationAddNode1Str = [...
        'Mutation.probabilityAddNode                 = [%s]\n',...
        '                                              [%s]\n'];
    mutationAddNode2Str = 'Mutation.probabilityAddNodeShrinkRate       = %s\n';
    mutationAddNode3Str = 'Mutation.probabilityAddNodeMin              = %s\n';
    mutationAddNode1Val = {...
        num2str(Mutation.probabilityAddNodeScheduling(1,:),'\t\t%2.2f'), ...
        num2str(Mutation.probabilityAddNodeScheduling(2,:),'\t\t%2.2f')};
    mutationAddNode2Val = 'N/A';
    mutationAddNode3Val = 'N/A';
end

if Mutation.probabilityAddConnectionScheduling(1) == -1
    mutationAddConn1 = 'Mutation.probabilityAddConnection           = %4.5f\n';
    mutationAddConn2 = 'Mutation.probabilityAddConnectionShrinkRate = %4.5f\n';
    mutationAddConn3 = 'Mutation.probabilityAddConnectionMin        = %4.5f\n';
    mutationAddConn1Val = {Mutation.probabilityAddConnection};
    mutationAddConn2Val = Mutation.probabilityAddConnectionShrinkRate;
    mutationAddConn3Val = Mutation.probabilityAddConnectionMin;
else
    mutationAddConn1Str = [...
        'Mutation.probabilityAddConnection           = [%s]\n',...
        '                                              [%s]\n'];
    mutationAddConn2Str = 'Mutation.probabilityAddConnectionShrinkRate = %s\n';
    mutationAddConn3Str = 'Mutation.probabilityAddConnectionMin        = %s\n';
    mutationAddConn1Val = {...
        num2str(Mutation.probabilityAddConnectionScheduling(1,:),'\t\t%2.2f'), ...
        num2str(Mutation.probabilityAddConnectionScheduling(2,:),'\t\t%2.2f')};
    mutationAddConn2Val = 'N/A';
    mutationAddConn3Val = 'N/A';
end

fprintf([...
    'Using the following parameters:\n',...
    'loadFlag                                    = %4.f\n',...
    'saveFlag                                    = %4.f\n',...
    'nMaxGeneration                              = %4.f\n',...
    'populationSize                              = %4.f\n',...
    'nInputNodes                                 = %4.f\n',...
    'nOutputNodes                                = %4.f\n',...
    'initialConnectedInputNodes                  = [%s]\n',...
    'fitnessLimit                                = %4.5f\n',...
    'Speciation.c1                               = %4.5f\n',...
    'Speciation.c2                               = %4.5f\n',...
    'Speciation.c3                               = %4.5f\n',...
    'Speciation.threshold                        = %4.5f\n',...
    'Stagnation.threshold                        = %4.5f\n',...
    'Stagnation.nGenerations                     = %4.f\n',...
    'Refocus.threshold                           = %4.5f\n',...
    'Refocus.nGenerations                        = %4.f\n',...
    'Initial.killPercentage                      = %4.5f\n',...
    'Initial.numberForKill                       = %4.5f\n',...
    'Initial.numberCopy                          = %4.f\n',...
    'Selection.pressure                          = %4.5f\n',...
    'Crossover.percentage                        = %4.5f\n',...
    'Crossover.probabilityInterspecies           = %4.5f\n',...
    'Crossover.probabilityMultipoint             = %4.5f\n',...
    mutationAddNode1Str,...
    mutationAddNode2Str,...
    mutationAddNode3Str,...
    mutationAddConn1Str,...
    mutationAddConn2Str,...
    mutationAddConn3Str,...
    'Mutation.probabilityRecurrency              = %4.5f\n',...
    'Mutation.probabilityRecurrencyShrinkRate    = %4.5f\n',...
    'Mutation.probabilityRecurrencyMin           = %4.5f\n',...
    'Mutation.probabilityMutateWeight            = %4.5f\n',...
    'Mutation.probabilityMutateWeightShrinkRate  = %4.5f\n',...
    'Mutation.probabilityMutateWeightMin         = %4.5f\n',...
    'Mutation.weightCap                          = %4.5f\n',...
    'Mutation.weightRange                        = %4.5f\n',...
    'Mutation.probabilityGeneReenabled           = %4.5f\n\n',...
    ],...
    loadFlag,...
    saveFlag,...
    nMaxGeneration,...
    populationSize,...
    nInputNodes,...
    nOutputNodes,...
    num2str(initialConnectedInputNodes),...
    fitnessLimit,...
    Speciation.c1,...
    Speciation.c2,...
    Speciation.c3,...
    Speciation.threshold,...
    Stagnation.threshold,...
    Stagnation.nGenerations,...
    Refocus.threshold,...
    Refocus.nGenerations,...
    Initial.killPercentage,...
    Initial.numberForKill,...
    Initial.numberCopy,...
    Selection.pressure,...
    Crossover.percentage,...
    Crossover.probabilityInterspecies,...
    Crossover.probabilityMultipoint,...
    mutationAddNode1Val{:},...
    mutationAddNode2Val,...
    mutationAddNode3Val,...
    mutationAddConn1Val{:},...
    mutationAddConn2Val,...
    mutationAddConn3Val,...
    Mutation.probabilityRecurrency,...
    Mutation.probabilityRecurrencyShrinkRate,...
    Mutation.probabilityRecurrencyMin,...
    Mutation.probabilityMutateWeight,...
    Mutation.probabilityMutateWeightShrinkRate,...
    Mutation.probabilityMutateWeightMin,...
    Mutation.weightCap,...
    Mutation.weightRange,...
    Mutation.probabilityGeneReenabled...
    );

%% Start pool

clusterName = 'local';
nWorkers = 8;
create_parpool(clusterName, nWorkers, false);

%% Evaluation function parameters

% Handle to evaluation function
evalFunHandle = @uav_evaluation_function;

% Use parallel evaluation
EvalFunParameters.useParallel = true;

% Default simulation parameters
EvalFunParameters.SimParam.absTol         = 'auto'; % 'auto'
EvalFunParameters.SimParam.relTol         = 'auto'; % 'auto'
EvalFunParameters.SimParam.stopTime       = '10';
EvalFunParameters.SimParam.solverType     = 'Variable-step'; % 'Variable-step', 'Fixed-step'
EvalFunParameters.SimParam.solver         = 'VariableStepAuto'; % 'VariableStepAuto'
EvalFunParameters.SimParam.minStep        = '0.000001'; % 'auto', '0.000001'
EvalFunParameters.SimParam.maxStep        = 'auto'; % 'auto'
EvalFunParameters.SimParam.maxConsecutiveMinStep = '100'; % '1'
EvalFunParameters.SimParam.fixedStep      = EvalFunParameters.SimParam.maxStep;
EvalFunParameters.SimParam.initialStep    = 'auto'; % 'auto'
EvalFunParameters.SimParam.simulationMode = 'accelerator'; % 'normal', 'accelerator', 'rapid-accelerator' 
EvalFunParameters.SimParam.fastRestart    = 'on';

% UAV parameters
EvalFunParameters.UavParam.g        = 9.80665;
EvalFunParameters.UavParam.nMotors  = 4;
EvalFunParameters.UavParam.m        = 0.312;
EvalFunParameters.UavParam.I        = [0.0032373,0,0;0,0.0032373,0;0,0,0.0058673];
EvalFunParameters.UavParam.PP       = 0.1591*[1,1,-1,-1;1,-1,-1,1;0,0,0,0];
EvalFunParameters.UavParam.Ra       = 0.5033*ones(4,1);
EvalFunParameters.UavParam.La       = 0*1e-3*ones(4,1);
EvalFunParameters.UavParam.ke       = 7.8e-3*ones(4,1);
EvalFunParameters.UavParam.km       = 3.7e-3*ones(4,1);
EvalFunParameters.UavParam.kD       = 1.140e-7*ones(4,1);
EvalFunParameters.UavParam.Jm       = 3.357e-5*ones(4,1);
EvalFunParameters.UavParam.Bm       = 1.23e-6*ones(4,1);
EvalFunParameters.UavParam.kT       = 2.980e-6*ones(4,1);
EvalFunParameters.UavParam.VaMin    = zeros(4,1);
EvalFunParameters.UavParam.VaMax    = 10*ones(4,1);

% Input limits
EvalFunParameters.UavParam.xDesiredMin      = -100;
EvalFunParameters.UavParam.xDesiredMax      = +100;
EvalFunParameters.UavParam.yDesiredMin      = -100;
EvalFunParameters.UavParam.yDesiredMax      = +100;
EvalFunParameters.UavParam.zDesiredMin      = -100;
EvalFunParameters.UavParam.zDesiredMax      = +100;
EvalFunParameters.UavParam.psiDesiredMin    = -pi;
EvalFunParameters.UavParam.psiDesiredMax    = +pi;
EvalFunParameters.UavParam.thetaDesiredMin  = -pi/2;
EvalFunParameters.UavParam.thetaDesiredMax  = +pi/2;
EvalFunParameters.UavParam.phiDesiredMin    = -pi;
EvalFunParameters.UavParam.phiDesiredMax    = +pi;
EvalFunParameters.UavParam.xMin             = -100;
EvalFunParameters.UavParam.xMax             = +100;
EvalFunParameters.UavParam.yMin             = -100;
EvalFunParameters.UavParam.yMax             = +100;
EvalFunParameters.UavParam.zMin             = -100;
EvalFunParameters.UavParam.zMax             = +100;
EvalFunParameters.UavParam.uMin             = -30;
EvalFunParameters.UavParam.uMax             = +30;
EvalFunParameters.UavParam.vMin             = -30;
EvalFunParameters.UavParam.vMax             = +30;
EvalFunParameters.UavParam.wMin             = -30;
EvalFunParameters.UavParam.wMax             = +30;
EvalFunParameters.UavParam.udotMin          = -60;
EvalFunParameters.UavParam.udotMax          = +60;
EvalFunParameters.UavParam.vdotMin          = -60;
EvalFunParameters.UavParam.vdotMax          = +60;
EvalFunParameters.UavParam.wdotMin          = -60;
EvalFunParameters.UavParam.wdotMax          = +60;
EvalFunParameters.UavParam.psiMin           = -pi;
EvalFunParameters.UavParam.psiMax           = +pi;
EvalFunParameters.UavParam.thetaMin         = -pi/2;
EvalFunParameters.UavParam.thetaMax         = +pi/2;
EvalFunParameters.UavParam.phiMin           = -pi;
EvalFunParameters.UavParam.phiMax           = +pi;
EvalFunParameters.UavParam.pMin             = -2.4;
EvalFunParameters.UavParam.pMax             = +2.4;
EvalFunParameters.UavParam.qMin             = -2.4;
EvalFunParameters.UavParam.qMax             = +2.4;
EvalFunParameters.UavParam.rMin             = -4.72;
EvalFunParameters.UavParam.rMax             = +4.72;
EvalFunParameters.UavParam.omega1Min        = 0;
EvalFunParameters.UavParam.omega1Max        = +800;
EvalFunParameters.UavParam.omega2Min        = 0;
EvalFunParameters.UavParam.omega2Max        = +800;
EvalFunParameters.UavParam.omega3Min        = 0;
EvalFunParameters.UavParam.omega3Max        = +800;
EvalFunParameters.UavParam.omega4Min        = 0;
EvalFunParameters.UavParam.omega4Max        = +800;

% Load trajectories from trajectory_definition.m
trajectory_definition();

% Save mat file with EvalFunParameters
folderName = 'eval_fun_parameters';
if ~exist(folderName,'dir')
    mkdir(folderName);
end
save([folderName,'/eval_fun_parameters_',dateAndTime,'.mat'],'EvalFunParameters');

%% Run NEAT

% Call NEAT
[Population, bestIndividual] = neat(...
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
    Stagnation,...
    Refocus,...
    Selection,...
    Crossover,...
    Mutation,...
    evalFunHandle, ...
    EvalFunParameters, ...
    dateAndTime)

% Build best neural network
[neuralNet, nDelays] = build_neural_network(bestIndividual);

% Generate function
genFunctionMod(neuralNet, 'uav_best_nn_controller', ...
    'MatrixOnly', 'yes', 'ShowLinks', 'yes');
pause(0.05);

% Number of previous states
nPreviousStates = nargin(funHandle) - nInputs;

% Post-process generated function
generated_nn_processing(...
    [strrep(pwd,'\','/'),'/uav_best_nn_controller.m'], ...
    'uav_best_nn_controller', nPreviousStates);

% View net
view(neuralNet);

% Terminate diary
diary('off');
