function [PopulationWithFitnesses, SimDataOut] = ...
    uav_evaluation_function(Population, EvalFunParameters, fitnessLimit, SimData)

%% Init

% It seems that initializing the SimData structure inside
% uav_evaluation_function() causes problem in parfor simulation (i.e. fitness
% are all the same)
SimDataOut = SimData;

% We store Population in PopulationWithFitnesses
PopulationWithFitnesses = Population;

% Number of individuals
nIndividuals = size(Population, 2);

% Number of inputs/outputs excluding previous layers' states
nInputs = 25;
nOutputs = 4;

% Create temporary directory to put the generated neural networks if it doesn't
% already exist
folderName = 'temp_uav_nn_controller';
if ~exist(folderName,'dir')
    mkdir(folderName);
end

% Number of trajectories
nTrajectories = numel(EvalFunParameters.Trajectory);

% Current generation number
iGeneration = EvalFunParameters.iGeneration;

%% Generate neural network functions for each individual

% Function name string
funNameString = 'uav_nn_controller';

% Folders name initialization
subFolderNameArray                  = cell(1,nIndividuals);
folderWithSubFolderName             = cell(1,nIndividuals);
folderWithSubFolderAndFileWOExtName = cell(1,nIndividuals);
folderWithSubFolderAndFileWExtName  = cell(1,nIndividuals);
workerDir                           = cell(1,nIndividuals);

% Main directory
mainDir = strrep(pwd,'\','/');

% Temporary subfolder. This helps avoid fuck up when using parfor
for iIndividual = 1:nIndividuals
    subFolderNameArray{iIndividual} = gen_temp_string(30);
end

% We create an array to log which individual failed to have their network built.
% This is to make the code more robust to errors. A 1 means the individual is ok
% and a zero means it is not.
statusIndividualRecord = ones(1,nIndividuals);

% We create the temporary directories and generate the neural network functions
parfor iIndividual = 1:nIndividuals
    
    % Suffle random number generator seed
    rng('shuffle');
    
    % Worker directory
    workerDir{iIndividual} = strrep(pwd,'\','/');
    
    % cd to main directory
    cd(mainDir);
    
    % Build neural network of current individual
    try
        [neuralNet, nDelays] = build_neural_network(Population(iIndividual));
    catch
        statusIndividualRecord(iIndividual) = 0;
        fprintf('Neural network build failed for individual % 4.f in generation % 4.f\n', ...
            iIndividual, iGeneration);
    end
    
    % Folder/subfolder path
    folderWithSubFolderName{iIndividual} = ...
        [folderName,'/',subFolderNameArray{iIndividual}];
    
    % Folder/subfolder/file without extension path
    folderWithSubFolderAndFileWOExtName{iIndividual} = ...
        [folderWithSubFolderName{iIndividual},'/',funNameString];
    
    % Folder/subfolder/file with extension path
    folderWithSubFolderAndFileWExtName{iIndividual} = ...
        [folderWithSubFolderAndFileWOExtName{iIndividual},'.m'];
    
    % Add temporary folder to path
    cd([mainDir,'/',folderName]);
    mkdir(subFolderNameArray{iIndividual});
    cd(mainDir);
    addpath(folderWithSubFolderName{iIndividual});
    
    % Generate function
    if statusIndividualRecord(iIndividual)
        genFunctionMod(neuralNet, ...
            folderWithSubFolderAndFileWOExtName{iIndividual}, ...
            'MatrixOnly', 'yes', 'ShowLinks', 'no');
    end
    
    % Rehash (https://goo.gl/IDtpzk)
    exist([funNameString,'.m'],'file');
    rehash();
    
    % Remove path
    rmpath(folderWithSubFolderName{iIndividual});
    rehash();
    
    % cd back to worker directory
    cd(workerDir{iIndividual});
    
end

% We add a little pause in order for the filesystem to record to new created files
pause(1);

% We then modify the generated function. Initially, this was in the same parfor
% loop above but it was not reliable (don't know why). Moving this step in a
% seperate parfor loop seems to have helped avoid problem.
nPreviousStatesArray = zeros(1,nIndividuals);
parfor iIndividual = 1:nIndividuals
    
    if statusIndividualRecord(iIndividual)
        % Suffle random number generator seed
        rng('shuffle');
        
        % Worker directory
        % workerDir{iIndividual} = strrep(pwd,'\','/');
        
        % cd to main directory
        cd(mainDir);
        
        % Add temporary folder to path
        addpath(folderWithSubFolderName{iIndividual});
        
        % Create function handle
        funHandle = str2func(funNameString);
        
        % https://goo.gl/IDtpzk
        exist([funNameString,'.m'],'file');
        rehash();
        
        % Number of previous states
        nPreviousStates = nargin(funHandle) - nInputs;
        if nPreviousStates > 0
            fprintf('nPreviousStates = % 4.f for individual % 4.f in generation % 4.f (%s)\n', ...
                nPreviousStates, iIndividual, iGeneration, subFolderNameArray{iIndividual});
        end
        nPreviousStatesArray(iIndividual) = nPreviousStates;
        
        % Post-process generated function
        generated_nn_processing(...
            folderWithSubFolderAndFileWExtName{iIndividual}, ...
            funNameString, nPreviousStates);
        
        % Remove path
        rmpath(folderWithSubFolderName{iIndividual});
        rehash();
        
        % cd back to worker directory
        cd(workerDir{iIndividual});
    end
    
end

pause(1);

%% Evaluate networks

% To avoid broadcasting the entire structure into the parfor loop
Trajectory = EvalFunParameters.Trajectory;
SimParam = EvalFunParameters.SimParam;
UavParam = EvalFunParameters.UavParam;
xMin        = UavParam.xMin;
xMax        = UavParam.xMax;
yMin        = UavParam.yMin;
yMax        = UavParam.yMax;
zMin        = UavParam.zMin;
zMax        = UavParam.zMax;
psiMin      = UavParam.psiMin;
psiMax      = UavParam.psiMax;
thetaMin    = UavParam.thetaMin;
thetaMax    = UavParam.thetaMax;
phiMin      = UavParam.phiMin;
phiMax      = UavParam.phiMax;

% Compute maximum time fitness
timeCellArray = {Trajectory(:).XDesired};
timeVector = zeros(1,nTrajectories);
for iTraj = 1:nTrajectories
    timeVector(iTraj) = timeCellArray{iTraj}.t(end);
end
maxTimeFitness = sum(timeVector);
% % maxTimeFitness = 100;

% Model name
modelName = 'uav_closed_loop_v5';

% https://goo.gl/uNCs79 (In parfor, how to create a temp folder for every worker)
if EvalFunParameters.useParallel
    
    %% Parallel evaluation
    
    % Load model
    load_system(modelName);
    
    % Shuffled individual vector
    % iIndividualShuffled = randperm(nIndividuals);
    
    % Cycle through population
    parfor iIndividual = 1:nIndividuals % iIndividualShuffled
        
        %iIndividualSh = iIndividualShuffled(iIndividual);
        
        % We only simulate if the individual is ok
        if statusIndividualRecord(iIndividual)
            
            % Suffle random number generator seed
            rng('shuffle');
            
            % Current work directory
            cwd = strrep(pwd,'\','/');
            
            % Add current work directory in current session
            addpath(cwd);
            
            % Create temporary folder and change current directory to this folder
            tmpdir = strrep(tempname,'\','/');
            mkdir(tmpdir);
            cd(tmpdir);
            
            % Add path to uav_nn_controller
            addpath([mainDir,'/',folderWithSubFolderName{iIndividual}]);
            rehash();
            
            % Load system
            load_system(modelName);
            % load_system(modelName);
            
            % Create function handle
            funHandle = str2func(funNameString);
            
            % Initial state of recurrent layer
            if nPreviousStatesArray(iIndividual) == 0
                initialState = 0;
            else
                initialState = zeros(1,nPreviousStatesArray(iIndividual));
            end
            
            % Neural network Id
            neuralNetworkId = iIndividual;
            
            % Initialize cell array that will contain simulation results
            t            = cell(1,nTrajectories);
            xDesired     = cell(1,nTrajectories);
            yDesired     = cell(1,nTrajectories);
            zDesired     = cell(1,nTrajectories);
            psiDesired   = cell(1,nTrajectories);
            thetaDesired = cell(1,nTrajectories);
            phiDesired   = cell(1,nTrajectories);
            x            = cell(1,nTrajectories);
            y            = cell(1,nTrajectories);
            z            = cell(1,nTrajectories);
            psi          = cell(1,nTrajectories);
            theta        = cell(1,nTrajectories);
            phi          = cell(1,nTrajectories);
            Va1          = cell(1,nTrajectories);
            Va2          = cell(1,nTrajectories);
            Va3          = cell(1,nTrajectories);
            Va4          = cell(1,nTrajectories);
            
            % Simulate trajectories
            rapidAccEn = false;
            if strcmp(SimParam.solverType,'Variable-step')
                set_param(modelName,...
                    'SolverType',SimParam.solverType,...
                    'Solver',SimParam.solver,...
                    'SimulationMode',SimParam.simulationMode,...
                    'InitialStep',SimParam.initialStep,...
                    'SaveFormat','StructureWithTime',...
                    'RelTol',SimParam.relTol,...
                    'AbsTol',SimParam.absTol,...
                    'MinStep',SimParam.minStep,...
                    'MaxStep',SimParam.maxStep,...
                    'MaxConsecutiveMinStep',SimParam.maxConsecutiveMinStep,...
                    'SFSimEcho','off',...
                    'SimIntegrity','off');
            else
                set_param(modelName,...
                    'SolverType',SimParam.solverType,...
                    'Solver',SimParam.solver,...
                    'SimulationMode',SimParam.simulationMode,...
                    'FixedStep',SimParam.fixedStep,...
                    'SaveFormat','StructureWithTime',...
                    'RelTol',SimParam.relTol,...
                    'AbsTol',SimParam.absTol,...
                    'SFSimEcho','off',...
                    'SimIntegrity','off');
            end
            if nPreviousStatesArray(iIndividual) > 0 && rapidAccEn
                set_param(modelName, ...
                    'SimulationMode','rapid-accelerator');
            end
            % set_param(modelName,'SimCtrlC','on');
            for iTraj = 1:nTrajectories
                
                %% Assign simulation variables to base workspace
                assignin('base', 'g', UavParam.g); %#ok
                assignin('base', 'nMotors', UavParam.nMotors); %#ok
                assignin('base', 'm', UavParam.m); %#ok
                assignin('base', 'I', UavParam.I); %#ok
                assignin('base', 'PP', UavParam.PP); %#ok
                assignin('base', 'Ra', UavParam.Ra); %#ok
                assignin('base', 'La', UavParam.La); %#ok
                assignin('base', 'ke', UavParam.ke); %#ok
                assignin('base', 'km', UavParam.km); %#ok
                assignin('base', 'kD', UavParam.kD); %#ok
                assignin('base', 'Jm', UavParam.Jm); %#ok
                assignin('base', 'Bm', UavParam.Bm); %#ok
                assignin('base', 'kT', UavParam.kT); %#ok
                assignin('base', 'VaMin', UavParam.VaMin); %#ok
                assignin('base', 'VaMax', UavParam.VaMax); %#ok
                assignin('base', 'omegaInitial', Trajectory(iTraj).IniCond.omegaInitial); %#ok
                assignin('base', 'omegaDelay', Trajectory(iTraj).IniCond.omegaDelay); %#ok
                assignin('base', 'iaInitial', Trajectory(iTraj).IniCond.iaInitial); %#ok
                assignin('base', 'iaDelay', Trajectory(iTraj).IniCond.iaDelay); %#ok
                assignin('base', 'vPtInitial', Trajectory(iTraj).IniCond.vPtInitial); %#ok
                assignin('base', 'vPtDelay', Trajectory(iTraj).IniCond.vPtDelay); %#ok
                assignin('base', 'xPInitial', Trajectory(iTraj).IniCond.xPInitial); %#ok
                assignin('base', 'xPDelay', Trajectory(iTraj).IniCond.xPDelay); %#ok
                assignin('base', 'omegabtInitial', Trajectory(iTraj).IniCond.omegabtInitial); %#ok
                assignin('base', 'omegabtDelay', Trajectory(iTraj).IniCond.omegabtDelay); %#ok
                assignin('base', 'eulerZYXInitial', Trajectory(iTraj).IniCond.eulerZYXInitial); %#ok
                assignin('base', 'eulerZYXDelay', Trajectory(iTraj).IniCond.eulerZYXDelay); %#ok
                assignin('base', 'XDesired', Trajectory(iTraj).XDesired); %#ok
                assignin('base', 'YDesired', Trajectory(iTraj).YDesired); %#ok
                assignin('base', 'ZDesired', Trajectory(iTraj).ZDesired); %#ok
                assignin('base', 'PsiDesired', Trajectory(iTraj).PsiDesired); %#ok
                assignin('base', 'ThetaDesired', Trajectory(iTraj).ThetaDesired); %#ok
                assignin('base', 'PhiDesired', Trajectory(iTraj).PhiDesired); %#ok
                assignin('base', 'neuralNetworkId', neuralNetworkId); %#ok
                assignin('base', 'nPreviousStates', nPreviousStatesArray(iIndividual)); %#ok
                assignin('base', 'initialState', initialState); %#ok
                assignin('base', 'funNameString', double(funNameString)); %#ok
                assignin('base', 'lowerXTol', Trajectory(iTraj).XDesired.tol(1)); %#ok
                assignin('base', 'upperXTol', Trajectory(iTraj).XDesired.tol(2)); %#ok
                assignin('base', 'lowerYTol', Trajectory(iTraj).YDesired.tol(1)); %#ok
                assignin('base', 'upperYTol', Trajectory(iTraj).YDesired.tol(2)); %#ok
                assignin('base', 'lowerZTol', Trajectory(iTraj).ZDesired.tol(1)); %#ok
                assignin('base', 'upperZTol', Trajectory(iTraj).ZDesired.tol(2)); %#ok
                assignin('base', 'lowerPsiTol', Trajectory(iTraj).PsiDesired.tol(1)); %#ok
                assignin('base', 'upperPsiTol', Trajectory(iTraj).PsiDesired.tol(2)); %#ok
                assignin('base', 'lowerThetaTol', Trajectory(iTraj).ThetaDesired.tol(1)); %#ok
                assignin('base', 'upperThetaTol', Trajectory(iTraj).ThetaDesired.tol(2)); %#ok
                assignin('base', 'lowerPhiTol', Trajectory(iTraj).PhiDesired.tol(1)); %#ok
                assignin('base', 'upperPhiTol', Trajectory(iTraj).PhiDesired.tol(2)); %#ok
                assignin('base', 'maxViolationX', Trajectory(iTraj).maxViolationX); %#ok
                assignin('base', 'maxViolationY', Trajectory(iTraj).maxViolationY); %#ok
                assignin('base', 'maxViolationZ', Trajectory(iTraj).maxViolationZ); %#ok
                assignin('base', 'maxViolationPsi', Trajectory(iTraj).maxViolationPsi); %#ok
                assignin('base', 'maxViolationTheta', Trajectory(iTraj).maxViolationTheta); %#ok
                assignin('base', 'maxViolationPhi', Trajectory(iTraj).maxViolationPhi); %#ok
                assignin('base', 'violationFactorX', Trajectory(iTraj).violationFactorX); %#ok
                assignin('base', 'violationFactorY', Trajectory(iTraj).violationFactorY); %#ok
                assignin('base', 'violationFactorZ', Trajectory(iTraj).violationFactorZ); %#ok
                assignin('base', 'violationFactorPsi', Trajectory(iTraj).violationFactorPsi); %#ok
                assignin('base', 'violationFactorTheta', Trajectory(iTraj).violationFactorTheta); %#ok
                assignin('base', 'violationFactorPhi', Trajectory(iTraj).violationFactorPhi); %#ok
                assignin('base', 'currentTimeThresholdX', Trajectory(iTraj).violationTimeThresholdX); %#ok
                assignin('base', 'currentTimeThresholdY', Trajectory(iTraj).violationTimeThresholdY); %#ok
                assignin('base', 'currentTimeThresholdZ', Trajectory(iTraj).violationTimeThresholdZ); %#ok
                assignin('base', 'currentTimeThresholdPsi', Trajectory(iTraj).violationTimeThresholdPsi); %#ok
                assignin('base', 'currentTimeThresholdTheta', Trajectory(iTraj).violationTimeThresholdTheta); %#ok
                assignin('base', 'currentTimeThresholdPhi', Trajectory(iTraj).violationTimeThresholdPhi); %#ok
                assignin('base', 'fullLogging', -1); %#ok
                assignin('base', 'loggingSampleTime', 0.01); %#ok
                assignin('base', 'xDesiredMin', UavParam.xDesiredMin); %#ok
                assignin('base', 'xDesiredMax', UavParam.xDesiredMax); %#ok
                assignin('base', 'yDesiredMin', UavParam.yDesiredMin); %#ok
                assignin('base', 'yDesiredMax', UavParam.yDesiredMax); %#ok
                assignin('base', 'zDesiredMin', UavParam.zDesiredMin); %#ok
                assignin('base', 'zDesiredMax', UavParam.zDesiredMax); %#ok
                assignin('base', 'psiDesiredMin', UavParam.psiDesiredMin); %#ok
                assignin('base', 'psiDesiredMax', UavParam.psiDesiredMax); %#ok
                assignin('base', 'thetaDesiredMin', UavParam.thetaDesiredMin); %#ok
                assignin('base', 'thetaDesiredMax', UavParam.thetaDesiredMax); %#ok
                assignin('base', 'phiDesiredMin', UavParam.phiDesiredMin); %#ok
                assignin('base', 'phiDesiredMax', UavParam.phiDesiredMax); %#ok
                assignin('base', 'xMin', UavParam.xMin); %#ok
                assignin('base', 'xMax', UavParam.xMax); %#ok
                assignin('base', 'yMin', UavParam.yMin); %#ok
                assignin('base', 'yMax', UavParam.yMax); %#ok
                assignin('base', 'zMin', UavParam.zMin); %#ok
                assignin('base', 'zMax', UavParam.zMax); %#ok
                assignin('base', 'uMin', UavParam.uMin); %#ok
                assignin('base', 'uMax', UavParam.uMax); %#ok
                assignin('base', 'vMin', UavParam.vMin); %#ok
                assignin('base', 'vMax', UavParam.vMax); %#ok
                assignin('base', 'wMin', UavParam.wMin); %#ok
                assignin('base', 'wMax', UavParam.wMax); %#ok
                assignin('base', 'udotMin', UavParam.udotMin); %#ok
                assignin('base', 'udotMax', UavParam.udotMax); %#ok
                assignin('base', 'vdotMin', UavParam.vdotMin); %#ok
                assignin('base', 'vdotMax', UavParam.vdotMax); %#ok
                assignin('base', 'wdotMin', UavParam.wdotMin); %#ok
                assignin('base', 'wdotMax', UavParam.wdotMax); %#ok
                assignin('base', 'psiMin', UavParam.psiMin); %#ok
                assignin('base', 'psiMax', UavParam.psiMax); %#ok
                assignin('base', 'thetaMin', UavParam.thetaMin); %#ok
                assignin('base', 'thetaMax', UavParam.thetaMax); %#ok
                assignin('base', 'phiMin', UavParam.phiMin); %#ok
                assignin('base', 'phiMax', UavParam.phiMax); %#ok
                assignin('base', 'pMin', UavParam.pMin); %#ok
                assignin('base', 'pMax', UavParam.pMax); %#ok
                assignin('base', 'qMin', UavParam.qMin); %#ok
                assignin('base', 'qMax', UavParam.qMax); %#ok
                assignin('base', 'rMin', UavParam.rMin); %#ok
                assignin('base', 'rMax', UavParam.rMax); %#ok
                assignin('base', 'omega1Min', UavParam.omega1Min); %#ok
                assignin('base', 'omega1Max', UavParam.omega1Max); %#ok
                assignin('base', 'omega2Min', UavParam.omega2Min); %#ok
                assignin('base', 'omega2Max', UavParam.omega2Max); %#ok
                assignin('base', 'omega3Min', UavParam.omega3Min); %#ok
                assignin('base', 'omega3Max', UavParam.omega3Max); %#ok
                assignin('base', 'omega4Min', UavParam.omega4Min); %#ok
                assignin('base', 'omega4Max', UavParam.omega4Max); %#ok
                
                %% Launch simulation
                try
                    if nPreviousStatesArray(iIndividual) == 0 || ~rapidAccEn
                        simout = sim(modelName,...
                            'StopTime',num2str(max(Trajectory(iTraj).XDesired.t)),...
                            'FastRestart',SimParam.fastRestart);
                    else
                        simout = sim(modelName,...
                            'StopTime',num2str(max(Trajectory(iTraj).XDesired.t)));
                    end
                catch ME
                    statusIndividualRecord(iIndividual) = 0;
                    fprintf('Simulation failed for individual % 8.f     in generation % 8.f. Error is:\n', ...
                        iIndividual, iGeneration);
                    fprintf('identifier: %s    message: %s\n', ME.identifier, ME.message);
                    break;
                end
                
                %% Extract simulation results
                simout1              = simout.get('simout1');
                t{iTraj}             = simout1.time;
                xDesired{iTraj}      = simout1.signals.values(:,1);
                yDesired{iTraj}      = simout1.signals.values(:,2);
                zDesired{iTraj}      = simout1.signals.values(:,3);
                psiDesired{iTraj}    = simout1.signals.values(:,4);
                thetaDesired{iTraj}  = simout1.signals.values(:,5);
                phiDesired{iTraj}    = simout1.signals.values(:,6);
                x{iTraj}             = simout1.signals.values(:,7);
                y{iTraj}             = simout1.signals.values(:,8);
                z{iTraj}             = simout1.signals.values(:,9);
                psi{iTraj}           = simout1.signals.values(:,16);
                theta{iTraj}         = simout1.signals.values(:,17);
                phi{iTraj}           = simout1.signals.values(:,18);
                Va1{iTraj}           = simout1.signals.values(:,41);
                Va2{iTraj}           = simout1.signals.values(:,42);
                Va3{iTraj}           = simout1.signals.values(:,43);
                Va4{iTraj}           = simout1.signals.values(:,44);
                
            end
            if nPreviousStatesArray(iIndividual) == 0 || ~rapidAccEn
                set_param(modelName,'FastRestart','off');
            end
            
            % Close system
            close_system(modelName,0);
            
            % Remove path and folder to uav_nn_controller
            rmpath([mainDir,'/',folderWithSubFolderName{iIndividual}]);
            try
                rmdir([mainDir,'/',folderWithSubFolderName{iIndividual}], 's');
            catch
                try
                    rmdir([mainDir,'/',folderWithSubFolderName{iIndividual}]);
                catch
                    fprintf('Unable to remove %s for individual % 8.f     in generation % 8.f\n', ...
                        [mainDir,'/',folderWithSubFolderName{iIndividual}], iIndividual, iGeneration);
                end
            end
            rehash();
            
            % Change back to work directory and remove temporary directory
            cd(cwd);
            try
                rmdir(tmpdir, 's');
            catch
                try
                    rmdir(tmpdir);
                catch
                    %                     fprintf('Unable to remove %s for individual % 8.f     in generation % 8.f\n', ...
                    %                         tmpdir, iIndividual, iGeneration);
                end
            end
            rmpath(cwd);
            
        end
        
        % Compute fitness for this individual
        if statusIndividualRecord(iIndividual)
            [fitnessSignal,rewardTime] = compute_signal_error_fitness(nTrajectories, maxTimeFitness, ...
                xMin, xMax, yMin, yMax, zMin, zMax, psiMin, psiMax, thetaMin, thetaMax, phiMin, phiMax, ...
                t, xDesired, yDesired, zDesired, psiDesired, thetaDesired, phiDesired, ...
                x, y, z, psi, theta, phi, Va1, Va2, Va3, Va4);
            fitness = compute_fitness(nTrajectories, Trajectory, ...
                t, xDesired, yDesired, zDesired, psiDesired, thetaDesired, phiDesired, ...
                x, y, z, psi, theta, phi, rewardTime);
            PopulationWithFitnesses(iIndividual).fitness = fitness + fitnessSignal;
        else
            PopulationWithFitnesses(iIndividual).fitness = 0.0;
        end
        
        fprintf('fitness = % 8.2f     for individual % 8.f     in generation % 8.f\n', ...
            PopulationWithFitnesses(iIndividual).fitness, iIndividual, iGeneration);
        
        % Store simulation results
        if statusIndividualRecord(iIndividual)
            for iTraj = 1:nTrajectories
                nPtsTmp = numel(t{iTraj});
                SimDataOut(iIndividual).t{iTraj}(1:nPtsTmp)           = t{iTraj};
                SimDataOut(iIndividual).xDesired{iTraj}(1:nPtsTmp)    = xDesired{iTraj};
                SimDataOut(iIndividual).yDesired{iTraj}(1:nPtsTmp)    = yDesired{iTraj};
                SimDataOut(iIndividual).zDesired{iTraj}(1:nPtsTmp)    = zDesired{iTraj};
                SimDataOut(iIndividual).psiDesired{iTraj}(1:nPtsTmp)  = psiDesired{iTraj};
                SimDataOut(iIndividual).thetaDesired{iTraj}(1:nPtsTmp) = thetaDesired{iTraj};
                SimDataOut(iIndividual).phiDesired{iTraj}(1:nPtsTmp)  = phiDesired{iTraj};
                SimDataOut(iIndividual).x{iTraj}(1:nPtsTmp)           = x{iTraj};
                SimDataOut(iIndividual).y{iTraj}(1:nPtsTmp)           = y{iTraj};
                SimDataOut(iIndividual).z{iTraj}(1:nPtsTmp)           = z{iTraj};
                SimDataOut(iIndividual).psi{iTraj}(1:nPtsTmp)         = psi{iTraj};
                SimDataOut(iIndividual).theta{iTraj}(1:nPtsTmp)       = theta{iTraj};
                SimDataOut(iIndividual).phi{iTraj}(1:nPtsTmp)         = phi{iTraj};
                SimDataOut(iIndividual).Va1{iTraj}(1:nPtsTmp)         = Va1{iTraj};
                SimDataOut(iIndividual).Va2{iTraj}(1:nPtsTmp)         = Va2{iTraj};
                SimDataOut(iIndividual).Va3{iTraj}(1:nPtsTmp)         = Va3{iTraj};
                SimDataOut(iIndividual).Va4{iTraj}(1:nPtsTmp)         = Va4{iTraj};
                
                SimDataOut(iIndividual).t{iTraj}           = SimDataOut(iIndividual).t{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).xDesired{iTraj}    = SimDataOut(iIndividual).xDesired{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).yDesired{iTraj}    = SimDataOut(iIndividual).yDesired{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).zDesired{iTraj}    = SimDataOut(iIndividual).zDesired{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).psiDesired{iTraj}  = SimDataOut(iIndividual).psiDesired{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).thetaDesired{iTraj} = SimDataOut(iIndividual).thetaDesired{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).phiDesired{iTraj}  = SimDataOut(iIndividual).phiDesired{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).x{iTraj}           = SimDataOut(iIndividual).x{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).y{iTraj}           = SimDataOut(iIndividual).y{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).z{iTraj}           = SimDataOut(iIndividual).z{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).psi{iTraj}         = SimDataOut(iIndividual).psi{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).theta{iTraj}       = SimDataOut(iIndividual).theta{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).phi{iTraj}         = SimDataOut(iIndividual).phi{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).Va1{iTraj}         = SimDataOut(iIndividual).Va1{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).Va2{iTraj}         = SimDataOut(iIndividual).Va2{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).Va3{iTraj}         = SimDataOut(iIndividual).Va3{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).Va4{iTraj}         = SimDataOut(iIndividual).Va4{iTraj}(1:nPtsTmp);
            end
        end
        
    end
    
else % Don't use parallel
    
    %% Serial evaluation
    
    % Load system (just once)
    load_system(modelName);
    
    % Cycle through population
    for iIndividual = 1:nIndividuals
        
        % We only simulate if the individual is ok
        if statusIndividualRecord(iIndividual)
            
            % Add temporary folder in path
            addpath(folderWithSubFolderName{iIndividual});
            rehash();
            
            % Create function handle
            funHandle = str2func(funNameString);
            
            % Initial state of recurrent layer
            if nPreviousStatesArray(iIndividual) == 0
                initialState = 0;
            else
                initialState = zeros(1,nPreviousStatesArray(iIndividual));
            end
            
            % Neural network Id
            neuralNetworkId = iIndividual;
            
            % Initialize cell array that will contain simulation results
            t            = cell(1,nTrajectories);
            xDesired     = cell(1,nTrajectories);
            yDesired     = cell(1,nTrajectories);
            zDesired     = cell(1,nTrajectories);
            psiDesired   = cell(1,nTrajectories);
            thetaDesired = cell(1,nTrajectories);
            phiDesired   = cell(1,nTrajectories);
            x            = cell(1,nTrajectories);
            y            = cell(1,nTrajectories);
            z            = cell(1,nTrajectories);
            psi          = cell(1,nTrajectories);
            theta        = cell(1,nTrajectories);
            phi          = cell(1,nTrajectories);
            Va1          = cell(1,nTrajectories);
            Va2          = cell(1,nTrajectories);
            Va3          = cell(1,nTrajectories);
            Va4          = cell(1,nTrajectories);
            
            % Simulate trajectories
            rapidAccEn = true;
            set_param(modelName,'FastRestart','off');
            if strcmp(SimParam.solverType,'Variable-step')
                set_param(modelName,...
                    'SolverType',SimParam.solverType,...
                    'Solver',SimParam.solver,...
                    'SimulationMode',SimParam.simulationMode,...
                    'InitialStep',SimParam.initialStep,...
                    'SaveFormat','StructureWithTime',...
                    'RelTol',SimParam.relTol,...
                    'AbsTol',SimParam.absTol,...
                    'MinStep',SimParam.minStep,...
                    'MaxStep',SimParam.maxStep,...
                    'MaxConsecutiveMinStep',SimParam.maxConsecutiveMinStep,...
                    'SFSimEcho','off',...
                    'SimIntegrity','off');
            else
                set_param(modelName,...
                    'SolverType',SimParam.solverType,...
                    'Solver',SimParam.solver,...
                    'SimulationMode',SimParam.simulationMode,...
                    'FixedStep',SimParam.fixedStep,...
                    'SaveFormat','StructureWithTime',...
                    'RelTol',SimParam.relTol,...
                    'AbsTol',SimParam.absTol,...
                    'SFSimEcho','off',...
                    'SimIntegrity','off');
            end
            if nPreviousStatesArray(iIndividual) > 0 && rapidAccEn
                set_param(modelName, ...
                    'SimulationMode','rapid-accelerator');
            end
            if nPreviousStatesArray(iIndividual) == 0 || ~rapidAccEn
                set_param(modelName,'FastRestart',SimParam.fastRestart);
            end
            for iTraj = 1:nTrajectories
                
                %% Assign simulation variables to base workspace
                assignin('base', 'g', UavParam.g);
                assignin('base', 'nMotors', UavParam.nMotors);
                assignin('base', 'm', UavParam.m);
                assignin('base', 'I', UavParam.I);
                assignin('base', 'PP', UavParam.PP);
                assignin('base', 'Ra', UavParam.Ra);
                assignin('base', 'La', UavParam.La);
                assignin('base', 'ke', UavParam.ke);
                assignin('base', 'km', UavParam.km);
                assignin('base', 'kD', UavParam.kD);
                assignin('base', 'Jm', UavParam.Jm);
                assignin('base', 'Bm', UavParam.Bm);
                assignin('base', 'kT', UavParam.kT);
                assignin('base', 'VaMin', UavParam.VaMin);
                assignin('base', 'VaMax', UavParam.VaMax);
                assignin('base', 'omegaInitial', Trajectory(iTraj).IniCond.omegaInitial);
                assignin('base', 'omegaDelay', Trajectory(iTraj).IniCond.omegaDelay);
                assignin('base', 'iaInitial', Trajectory(iTraj).IniCond.iaInitial);
                assignin('base', 'iaDelay', Trajectory(iTraj).IniCond.iaDelay);
                assignin('base', 'vPtInitial', Trajectory(iTraj).IniCond.vPtInitial);
                assignin('base', 'vPtDelay', Trajectory(iTraj).IniCond.vPtDelay)
                assignin('base', 'xPInitial', Trajectory(iTraj).IniCond.xPInitial);
                assignin('base', 'xPDelay', Trajectory(iTraj).IniCond.xPDelay);
                assignin('base', 'omegabtInitial', Trajectory(iTraj).IniCond.omegabtInitial);
                assignin('base', 'omegabtDelay', Trajectory(iTraj).IniCond.omegabtDelay);
                assignin('base', 'eulerZYXInitial', Trajectory(iTraj).IniCond.eulerZYXInitial);
                assignin('base', 'eulerZYXDelay', Trajectory(iTraj).IniCond.eulerZYXDelay);
                assignin('base', 'XDesired', Trajectory(iTraj).XDesired);
                assignin('base', 'YDesired', Trajectory(iTraj).YDesired);
                assignin('base', 'ZDesired', Trajectory(iTraj).ZDesired);
                assignin('base', 'PsiDesired', Trajectory(iTraj).PsiDesired);
                assignin('base', 'ThetaDesired', Trajectory(iTraj).ThetaDesired);
                assignin('base', 'PhiDesired', Trajectory(iTraj).PhiDesired);
                assignin('base', 'neuralNetworkId', neuralNetworkId);
                assignin('base', 'nPreviousStates', nPreviousStatesArray(iIndividual));
                assignin('base', 'initialState', initialState);
                assignin('base', 'funNameString', double(funNameString));
                assignin('base', 'lowerXTol', Trajectory(iTraj).XDesired.tol(1));
                assignin('base', 'upperXTol', Trajectory(iTraj).XDesired.tol(2));
                assignin('base', 'lowerYTol', Trajectory(iTraj).YDesired.tol(1));
                assignin('base', 'upperYTol', Trajectory(iTraj).YDesired.tol(2));
                assignin('base', 'lowerZTol', Trajectory(iTraj).ZDesired.tol(1));
                assignin('base', 'upperZTol', Trajectory(iTraj).ZDesired.tol(2));
                assignin('base', 'lowerPsiTol', Trajectory(iTraj).PsiDesired.tol(1));
                assignin('base', 'upperPsiTol', Trajectory(iTraj).PsiDesired.tol(2));
                assignin('base', 'lowerThetaTol', Trajectory(iTraj).ThetaDesired.tol(1));
                assignin('base', 'upperThetaTol', Trajectory(iTraj).ThetaDesired.tol(2));
                assignin('base', 'lowerPhiTol', Trajectory(iTraj).PhiDesired.tol(1));
                assignin('base', 'upperPhiTol', Trajectory(iTraj).PhiDesired.tol(2));
                assignin('base', 'maxViolationX', Trajectory(iTraj).maxViolationX);
                assignin('base', 'maxViolationY', Trajectory(iTraj).maxViolationY);
                assignin('base', 'maxViolationZ', Trajectory(iTraj).maxViolationZ);
                assignin('base', 'maxViolationPsi', Trajectory(iTraj).maxViolationPsi);
                assignin('base', 'maxViolationTheta', Trajectory(iTraj).maxViolationTheta);
                assignin('base', 'maxViolationPhi', Trajectory(iTraj).maxViolationPhi);
                assignin('base', 'violationFactorX', Trajectory(iTraj).violationFactorX);
                assignin('base', 'violationFactorY', Trajectory(iTraj).violationFactorY);
                assignin('base', 'violationFactorZ', Trajectory(iTraj).violationFactorZ);
                assignin('base', 'violationFactorPsi', Trajectory(iTraj).violationFactorPsi);
                assignin('base', 'violationFactorTheta', Trajectory(iTraj).violationFactorTheta);
                assignin('base', 'violationFactorPhi', Trajectory(iTraj).violationFactorPhi);
                assignin('base', 'currentTimeThresholdX', Trajectory(iTraj).violationTimeThresholdX);
                assignin('base', 'currentTimeThresholdY', Trajectory(iTraj).violationTimeThresholdY);
                assignin('base', 'currentTimeThresholdZ', Trajectory(iTraj).violationTimeThresholdZ);
                assignin('base', 'currentTimeThresholdPsi', Trajectory(iTraj).violationTimeThresholdPsi);
                assignin('base', 'currentTimeThresholdTheta', Trajectory(iTraj).violationTimeThresholdTheta);
                assignin('base', 'currentTimeThresholdPhi', Trajectory(iTraj).violationTimeThresholdPhi);
                assignin('base', 'fullLogging', -1);
                assignin('base', 'loggingSampleTime', 0.01);
                assignin('base', 'xDesiredMin', UavParam.xDesiredMin);
                assignin('base', 'xDesiredMax', UavParam.xDesiredMax);
                assignin('base', 'yDesiredMin', UavParam.yDesiredMin);
                assignin('base', 'yDesiredMax', UavParam.yDesiredMax);
                assignin('base', 'zDesiredMin', UavParam.zDesiredMin);
                assignin('base', 'zDesiredMax', UavParam.zDesiredMax);
                assignin('base', 'psiDesiredMin', UavParam.psiDesiredMin);
                assignin('base', 'psiDesiredMax', UavParam.psiDesiredMax);
                assignin('base', 'thetaDesiredMin', UavParam.thetaDesiredMin);
                assignin('base', 'thetaDesiredMax', UavParam.thetaDesiredMax);
                assignin('base', 'phiDesiredMin', UavParam.phiDesiredMin);
                assignin('base', 'phiDesiredMax', UavParam.phiDesiredMax);
                assignin('base', 'xMin', UavParam.xMin);
                assignin('base', 'xMax', UavParam.xMax);
                assignin('base', 'yMin', UavParam.yMin);
                assignin('base', 'yMax', UavParam.yMax);
                assignin('base', 'zMin', UavParam.zMin);
                assignin('base', 'zMax', UavParam.zMax);
                assignin('base', 'uMin', UavParam.uMin);
                assignin('base', 'uMax', UavParam.uMax);
                assignin('base', 'vMin', UavParam.vMin);
                assignin('base', 'vMax', UavParam.vMax);
                assignin('base', 'wMin', UavParam.wMin);
                assignin('base', 'wMax', UavParam.wMax);
                assignin('base', 'udotMin', UavParam.udotMin);
                assignin('base', 'udotMax', UavParam.udotMax);
                assignin('base', 'vdotMin', UavParam.vdotMin);
                assignin('base', 'vdotMax', UavParam.vdotMax);
                assignin('base', 'wdotMin', UavParam.wdotMin);
                assignin('base', 'wdotMax', UavParam.wdotMax);
                assignin('base', 'psiMin', UavParam.psiMin);
                assignin('base', 'psiMax', UavParam.psiMax);
                assignin('base', 'thetaMin', UavParam.thetaMin);
                assignin('base', 'thetaMax', UavParam.thetaMax);
                assignin('base', 'phiMin', UavParam.phiMin);
                assignin('base', 'phiMax', UavParam.phiMax);
                assignin('base', 'pMin', UavParam.pMin);
                assignin('base', 'pMax', UavParam.pMax);
                assignin('base', 'qMin', UavParam.qMin);
                assignin('base', 'qMax', UavParam.qMax);
                assignin('base', 'rMin', UavParam.rMin);
                assignin('base', 'rMax', UavParam.rMax);
                assignin('base', 'omega1Min', UavParam.omega1Min);
                assignin('base', 'omega1Max', UavParam.omega1Max);
                assignin('base', 'omega2Min', UavParam.omega2Min);
                assignin('base', 'omega2Max', UavParam.omega2Max);
                assignin('base', 'omega3Min', UavParam.omega3Min);
                assignin('base', 'omega3Max', UavParam.omega3Max);
                assignin('base', 'omega4Min', UavParam.omega4Min);
                assignin('base', 'omega4Max', UavParam.omega4Max);
                
                %% Launch simulation
                try
                    simout = sim(modelName,...
                        'StopTime',num2str(max(Trajectory(iTraj).XDesired.t)));
                catch ME
                    statusIndividualRecord(iIndividual) = 0;
                    fprintf('Simulation failed for individual % 8.f     in generation % 8.f. Error is:\n', ...
                        iIndividual, iGeneration);
                    fprintf('identifier: %s    message: %s\n', ME.identifier, ME.message);
                    break;
                end
                
                %% Extract simulation results
                simout1              = simout.get('simout1');
                t{iTraj}             = simout1.time;
                xDesired{iTraj}      = simout1.signals.values(:,1);
                yDesired{iTraj}      = simout1.signals.values(:,2);
                zDesired{iTraj}      = simout1.signals.values(:,3);
                psiDesired{iTraj}    = simout1.signals.values(:,4);
                thetaDesired{iTraj}  = simout1.signals.values(:,5);
                phiDesired{iTraj}    = simout1.signals.values(:,6);
                x{iTraj}             = simout1.signals.values(:,7);
                y{iTraj}             = simout1.signals.values(:,8);
                z{iTraj}             = simout1.signals.values(:,9);
                psi{iTraj}           = simout1.signals.values(:,16);
                theta{iTraj}         = simout1.signals.values(:,17);
                phi{iTraj}           = simout1.signals.values(:,18);
                Va1{iTraj}           = simout1.signals.values(:,41);
                Va2{iTraj}           = simout1.signals.values(:,42);
                Va3{iTraj}           = simout1.signals.values(:,43);
                Va4{iTraj}           = simout1.signals.values(:,44);
                
            end
            if nPreviousStatesArray(iIndividual) == 0 || ~rapidAccEn
                set_param(modelName,'FastRestart','off');
            end
            
            % Remove temporary folder in path
            rmpath(folderWithSubFolderName{iIndividual});
            try
                rmdir(folderWithSubFolderName{iIndividual}, 's');
            catch
                try
                    rmdir(folderWithSubFolderName{iIndividual});
                catch
                    fprintf('Unable to remove %s for individual % 8.f     in generation % 8.f\n', ...
                        folderWithSubFolderName{iIndividual}, iIndividual, iGeneration);
                end
            end
            rehash();
            
        end
        
        % Compute fitness for this individual
        if statusIndividualRecord(iIndividual)
            [fitnessSignal,rewardTime] = compute_signal_error_fitness(nTrajectories, maxTimeFitness, ...
                xMin, xMax, yMin, yMax, zMin, zMax, psiMin, psiMax, thetaMin, thetaMax, phiMin, phiMax, ...
                t, xDesired, yDesired, zDesired, psiDesired, thetaDesired, phiDesired, ...
                x, y, z, psi, theta, phi, Va1, Va2, Va3, Va4);
            fitness = compute_fitness(nTrajectories, Trajectory, ...
                t, xDesired, yDesired, zDesired, psiDesired, thetaDesired, phiDesired, ...
                x, y, z, psi, theta, phi, rewardTime);
            PopulationWithFitnesses(iIndividual).fitness = fitness + fitnessSignal;
        else
            PopulationWithFitnesses(iIndividual).fitness = 0.0;
        end
        
        fprintf('fitness = % 8.2f     for individual % 8.f     in generation % 8.f\n', ...
            PopulationWithFitnesses(iIndividual).fitness, iIndividual, iGeneration);
        
        % Store simulation results
        if statusIndividualRecord(iIndividual)
            for iTraj = 1:nTrajectories
                nPtsTmp = numel(t{iTraj});
                SimDataOut(iIndividual).t{iTraj}(1:nPtsTmp)           = t{iTraj};
                SimDataOut(iIndividual).xDesired{iTraj}(1:nPtsTmp)    = xDesired{iTraj};
                SimDataOut(iIndividual).yDesired{iTraj}(1:nPtsTmp)    = yDesired{iTraj};
                SimDataOut(iIndividual).zDesired{iTraj}(1:nPtsTmp)    = zDesired{iTraj};
                SimDataOut(iIndividual).psiDesired{iTraj}(1:nPtsTmp)  = psiDesired{iTraj};
                SimDataOut(iIndividual).thetaDesired{iTraj}(1:nPtsTmp) = thetaDesired{iTraj};
                SimDataOut(iIndividual).phiDesired{iTraj}(1:nPtsTmp)  = phiDesired{iTraj};
                SimDataOut(iIndividual).x{iTraj}(1:nPtsTmp)           = x{iTraj};
                SimDataOut(iIndividual).y{iTraj}(1:nPtsTmp)           = y{iTraj};
                SimDataOut(iIndividual).z{iTraj}(1:nPtsTmp)           = z{iTraj};
                SimDataOut(iIndividual).psi{iTraj}(1:nPtsTmp)         = psi{iTraj};
                SimDataOut(iIndividual).theta{iTraj}(1:nPtsTmp)       = theta{iTraj};
                SimDataOut(iIndividual).phi{iTraj}(1:nPtsTmp)         = phi{iTraj};
                SimDataOut(iIndividual).Va1{iTraj}(1:nPtsTmp)         = Va1{iTraj};
                SimDataOut(iIndividual).Va2{iTraj}(1:nPtsTmp)         = Va2{iTraj};
                SimDataOut(iIndividual).Va3{iTraj}(1:nPtsTmp)         = Va3{iTraj};
                SimDataOut(iIndividual).Va4{iTraj}(1:nPtsTmp)         = Va4{iTraj};
                
                SimDataOut(iIndividual).t{iTraj}           = SimDataOut(iIndividual).t{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).xDesired{iTraj}    = SimDataOut(iIndividual).xDesired{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).yDesired{iTraj}    = SimDataOut(iIndividual).yDesired{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).zDesired{iTraj}    = SimDataOut(iIndividual).zDesired{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).psiDesired{iTraj}  = SimDataOut(iIndividual).psiDesired{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).thetaDesired{iTraj} = SimDataOut(iIndividual).thetaDesired{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).phiDesired{iTraj}  = SimDataOut(iIndividual).phiDesired{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).x{iTraj}           = SimDataOut(iIndividual).x{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).y{iTraj}           = SimDataOut(iIndividual).y{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).z{iTraj}           = SimDataOut(iIndividual).z{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).psi{iTraj}         = SimDataOut(iIndividual).psi{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).theta{iTraj}       = SimDataOut(iIndividual).theta{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).phi{iTraj}         = SimDataOut(iIndividual).phi{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).Va1{iTraj}         = SimDataOut(iIndividual).Va1{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).Va2{iTraj}         = SimDataOut(iIndividual).Va2{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).Va3{iTraj}         = SimDataOut(iIndividual).Va3{iTraj}(1:nPtsTmp);
                SimDataOut(iIndividual).Va4{iTraj}         = SimDataOut(iIndividual).Va4{iTraj}(1:nPtsTmp);
            end
        end
        
    end
    
    % Close system
    close_system(modelName,0);
    
end

end

%% Function to compute the fitness of an individual
function f = compute_fitness(nTrajectories, Trajectory, ...
    t, xDesired, yDesired, zDesired, psiDesired, thetaDesired, phiDesired, ...
    x, y, z, psi, theta, phi, rewardTime)

% Initialize total fitness to zero
f = 0;

% Compute fitness for each trajectory
for iTraj = 1:nTrajectories
    
    % Compute the sample time vector (i.e. t_i - t_i-1)
    deltaT = [0,diff(t{iTraj})'];
    
    % Compute the parallel trajectory curves according to their tolerance
    [xLower, xUpper] = compute_lower_upper_bound(Trajectory(iTraj).XDesired.tol(1), ...
        Trajectory(iTraj).XDesired.tol(2), t{iTraj}, xDesired{iTraj});
    
    [yLower, yUpper] = compute_lower_upper_bound(Trajectory(iTraj).YDesired.tol(1), ...
        Trajectory(iTraj).YDesired.tol(2), t{iTraj}, yDesired{iTraj});
    
    [zLower, zUpper] = compute_lower_upper_bound(Trajectory(iTraj).ZDesired.tol(1), ...
        Trajectory(iTraj).ZDesired.tol(2), t{iTraj}, zDesired{iTraj});
    
    [psiLower, psiUpper] = compute_lower_upper_bound(Trajectory(iTraj).PsiDesired.tol(1), ...
        Trajectory(iTraj).PsiDesired.tol(2), t{iTraj}, psiDesired{iTraj});
    
    [thetaLower, thetaUpper] = compute_lower_upper_bound(Trajectory(iTraj).ThetaDesired.tol(1), ...
        Trajectory(iTraj).ThetaDesired.tol(2), t{iTraj}, thetaDesired{iTraj});
    
    [phiLower, phiUpper] = compute_lower_upper_bound(Trajectory(iTraj).PhiDesired.tol(1), ...
        Trajectory(iTraj).PhiDesired.tol(2), t{iTraj}, phiDesired{iTraj});
    
    % Compute the total duration that each trajectory is within the tolerance
    durationX       = sum(deltaT((xLower     <= x{iTraj})     &   (x{iTraj}     <= xUpper)));
    durationY       = sum(deltaT((yLower     <= y{iTraj})     &   (y{iTraj}     <= yUpper)));
    durationZ       = sum(deltaT((zLower     <= z{iTraj})     &   (z{iTraj}     <= zUpper)));
    durationPsi     = sum(deltaT((psiLower   <= psi{iTraj})   &   (psi{iTraj}   <= psiUpper)));
    durationTheta   = sum(deltaT((thetaLower <= theta{iTraj}) &   (theta{iTraj} <= thetaUpper)));
    durationPhi     = sum(deltaT((phiLower   <= phi{iTraj})   &   (phi{iTraj}   <= phiUpper)));
    
    % Compute total fitness
    f = f + rewardTime(iTraj)*(durationX + durationY + durationZ + durationPsi + durationTheta + durationPhi);
    
end

end

%% Function to compute the fitness associated to the trajectory error
function [fitnessSignal,rewardTime] = compute_signal_error_fitness(nTrajectories, maxTimeFitness, ...
    xMin, xMax, yMin, yMax, zMin, zMax, psiMin, psiMax, thetaMin, thetaMax, phiMin, phiMax, ...
    t, xDesired, yDesired, zDesired, psiDesired, thetaDesired, phiDesired, ...
    x, y, z, psi, theta, phi, Va1, Va2, Va3, Va4)

% Initialize total fitness to zero
fitnessSignal = 0;
fitnessSignalMax = 0;
rewardTime = zeros(1,nTrajectories);

% Compute fitness for each trajectory
for iTraj = 1:nTrajectories
    
    % Compute fitnesses
    [fitnessSignalX,fitnessMaxX] = sub_fun_error_fitness(xMin, xMax, xDesired{iTraj}, x{iTraj});
    [fitnessSignalY,fitnessMaxY] = sub_fun_error_fitness(yMin, yMax, yDesired{iTraj}, y{iTraj});
    [fitnessSignalZ,fitnessMaxZ,signalErrorSignedZ] = ...
        sub_fun_error_fitness(zMin, zMax, zDesired{iTraj}, z{iTraj});
    [fitnessSignalPsi,fitnessMaxPsi,signalErrorSignedPsi] = ...
        sub_fun_error_fitness(psiMin, psiMax, psiDesired{iTraj}, psi{iTraj});
    [fitnessSignalTheta,fitnessMaxTheta] = sub_fun_error_fitness(thetaMin, thetaMax, thetaDesired{iTraj}, theta{iTraj});
    [fitnessSignalPhi,fitnessMaxPhi]     = sub_fun_error_fitness(phiMin, phiMax, phiDesired{iTraj}, phi{iTraj});
    
    % Conditions to penalize controllers that send very low or very high voltages
    VaMinThreshold      = 1;
    VaMinRatioThreshold = 0.80;
    VaMinTimeThreshold  = 2;
    VaMaxThreshold      = 8.8;
    VaMaxRatioThreshold = 0.80;
    VaMaxTimeThreshold  = 2;
    [condVaMin,condVaMax] = check_cond_Va(t{iTraj},Va1{iTraj},Va2{iTraj},Va3{iTraj},Va4{iTraj}, ...
        VaMinThreshold,VaMinRatioThreshold,VaMinTimeThreshold, ...
        VaMaxThreshold,VaMaxRatioThreshold,VaMaxTimeThreshold);
    % VaMinThreshold = 1;
    % condVaMin1 = ((max(Va1{iTraj}) - min(Va1{iTraj}))) < VaMinThreshold && (max(Va1{iTraj}) < VaMinThreshold);
    % condVaMin2 = ((max(Va2{iTraj}) - min(Va2{iTraj}))) < VaMinThreshold && (max(Va2{iTraj}) < VaMinThreshold);
    % condVaMin3 = ((max(Va3{iTraj}) - min(Va3{iTraj}))) < VaMinThreshold && (max(Va3{iTraj}) < VaMinThreshold);
    % condVaMin4 = ((max(Va4{iTraj}) - min(Va4{iTraj}))) < VaMinThreshold && (max(Va4{iTraj}) < VaMinThreshold);
    
    % Conditions to penalize controllers that let the altitude deacrease or increase
    % sharply in Z and/or psi. Compute gradient of signalErrorZ and signalErrorPsi:
    minErrRatioThrshld = 0.80;
    errZThrshld        = 0.5;
    maxErrZThrshld     = 2;
    errPsiThrshld      = 4*pi/180;
    maxErrPsiThrshld   = 8*pi/180;
    gradErrZ    = gradient(signalErrorSignedZ);
    gradErrPsi  = gradient(signalErrorSignedPsi);
    [condGradPosZ,condGradNegZ] = check_cond_grad(gradErrZ, signalErrorSignedZ, ...
        errZThrshld, maxErrZThrshld, minErrRatioThrshld);
    [condGradPosPsi,condGradNegPsi] = check_cond_grad(gradErrPsi, signalErrorSignedPsi, ...
        errPsiThrshld, maxErrPsiThrshld, minErrRatioThrshld);
    
    % Compute rewards
    overallCondGrad = condGradPosZ || condGradNegZ || condGradPosPsi || condGradNegPsi;
    [rewardTime(iTraj),rewardZ,rewardPsi] = compute_reward(overallCondGrad,...
        condGradPosZ,condGradNegZ,condGradPosPsi,condGradNegPsi,t{iTraj});
    
    % Compute total fitness
    if overallCondGrad && ~condVaMin && ~condVaMax
        fitnessSignal = fitnessSignal + fitnessSignalX + fitnessSignalY + ...
            rewardZ*fitnessSignalZ + rewardPsi*fitnessSignalPsi + ...
            fitnessSignalTheta + fitnessSignalPhi;
        fitnessSignalMax = fitnessSignalMax + fitnessMaxX + fitnessMaxY + ...
            fitnessMaxZ + fitnessMaxPsi + fitnessMaxTheta + fitnessMaxPhi;
    elseif condVaMin || condVaMax
        fitnessSignalMax = fitnessSignalMax + fitnessMaxX + fitnessMaxY + ...
            fitnessMaxZ + fitnessMaxPsi + fitnessMaxTheta + fitnessMaxPhi;
    else
        fitnessSignal = fitnessSignal + fitnessSignalX + fitnessSignalY + ...
            rewardZ*fitnessSignalZ + rewardPsi*fitnessSignalPsi + ...
            fitnessSignalTheta + fitnessSignalPhi;
        fitnessSignalMax = fitnessSignalMax + fitnessMaxX + fitnessMaxY + ...
            fitnessMaxZ + fitnessMaxPsi + fitnessMaxTheta + fitnessMaxPhi;
    end
    
end

% Fitness scaling
fitnessSignal = maxTimeFitness*fitnessSignal/fitnessSignalMax;

    function [fitnessSignal,fitnessMax,signalErrorSigned] = sub_fun_error_fitness(qMin, qMax, qDesired, q)
        % Compute max possible error
        sumSignalError1    = sum(abs(qDesired - qMin*ones(size(qDesired))));
        sumSignalError2    = sum(abs(qDesired - qMax*ones(size(qDesired))));
        sumSignalErrorMax  = max([sumSignalError1,sumSignalError2]);
        % Compute actual error
        signalErrorSigned = q - qDesired;
        signalError     = abs(signalErrorSigned);
        sumSignalError  = sum(signalError);
        % Fitness
        fitnessSignal   = abs(sumSignalErrorMax - sumSignalError);
        fitnessMax      = abs(sumSignalErrorMax);
    end

    function [condVaMin,condVaMax] = check_cond_Va(t,Va1,Va2,Va3,Va4, ...
            VaMinThreshold,VaMinRatioThreshold,VaMinTimeThreshold, ...
            VaMaxThreshold,VaMaxRatioThreshold,VaMaxTimeThreshold)
        % Conditions to penalize controllers that send very low voltages
        condVaMin1 = sum(Va1 < VaMinThreshold) >= VaMinRatioThreshold*numel(Va1);
        condVaMin2 = sum(Va2 < VaMinThreshold) >= VaMinRatioThreshold*numel(Va2);
        condVaMin3 = sum(Va3 < VaMinThreshold) >= VaMinRatioThreshold*numel(Va3);
        condVaMin4 = sum(Va4 < VaMinThreshold) >= VaMinRatioThreshold*numel(Va4);
        % Conditions to penalize controllers that send very high voltages
        condVaMax1 = sum(Va1 > VaMaxThreshold) >= VaMaxRatioThreshold*numel(Va1);
        condVaMax2 = sum(Va2 > VaMaxThreshold) >= VaMaxRatioThreshold*numel(Va2);
        condVaMax3 = sum(Va3 > VaMaxThreshold) >= VaMaxRatioThreshold*numel(Va3);
        condVaMax4 = sum(Va4 > VaMaxThreshold) >= VaMaxRatioThreshold*numel(Va4);
        % Compute final condition
        if t(end) >= VaMinTimeThreshold
            condVaMin = condVaMin1 || condVaMin2 || condVaMin3 || condVaMin4;
        else
            condVaMin = false;
        end
        if t(end) >= VaMaxTimeThreshold
            condVaMax = condVaMax1 || condVaMax2 || condVaMax3 || condVaMax4;
        else
            condVaMax = false;
        end
    end

    function [condGradPos,condGradNeg] = check_cond_grad(gradErr, signalErrorSigned, ...
            errThrshld, maxErrThrshld, minErrRatioThrshld)
        % Positive gradient
        condGradPosa = sum(gradErr > 0 & abs(signalErrorSigned) > errThrshld) >= minErrRatioThrshld*numel(gradErr);
        condGradPosb = max(signalErrorSigned) > maxErrThrshld;
        condGradPos  = condGradPosa && condGradPosb;
        % Negative gradient
        condGradNega = sum(gradErr < 0 & abs(signalErrorSigned) > errThrshld) >= minErrRatioThrshld*numel(gradErr);
        condGradNegb = min(signalErrorSigned) < -maxErrThrshld;
        condGradNeg  = condGradNega && condGradNegb;
    end

    function [rewardTime,rewardZ,rewardPsi] = compute_reward(overallCondGrad,...
            condGradPosZ,condGradNegZ,condGradPosPsi,condGradNegPsi,t)
        if overallCondGrad
            if t(end) >= 0 && t(end) < 3
                rewardTime = 0.5;
            elseif t(end) >= 3 && t(end) < 6
                rewardTime = 0.25;
            elseif t(end) >= 6
                rewardTime = 0.125;
            end
            if t(end) >= 1 && t(end) < 3
                rewardZ    = 1;
                rewardPsi  = 1;
            elseif t(end) >= 3 && t(end) < 6
                rewardZ    = 1.25;
                rewardPsi  = 1.25;
            elseif t(end) >= 6 && t(end) < 8
                rewardZ    = 1.5;
                rewardPsi  = 1.5;
            elseif t(end) >= 8
                rewardZ    = 1.75;
                rewardPsi  = 1.75;
            else
                rewardZ    = 0.5;
                rewardPsi  = 0.5;
            end
            if condGradPosZ || condGradNegZ
                rewardZ = 0;
            end
            if condGradPosPsi || condGradNegPsi
                rewardPsi = 0;
            end
        else
            if t(end) >= 1 && t(end) < 3
                rewardTime = 1.5;
                rewardZ    = 1.5;
                rewardPsi  = 1.5;
            elseif t(end) >= 3 && t(end) < 6
                rewardTime = 2;
                rewardZ    = 2;
                rewardPsi  = 2;
            elseif t(end) >= 6 && t(end) < 8
                rewardTime = 2.5;
                rewardZ    = 2.5;
                rewardPsi  = 2.5;
            elseif t(end) >= 8
                rewardTime = 3;
                rewardZ    = 3;
                rewardPsi  = 3;
            else
                rewardTime = 1;
                rewardZ    = 1;
                rewardPsi  = 1;
            end
        end
    end

end
