classdef SimulationEngine < handle
    %SIMULATIONENGINE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %simualtion environment, contain all calibrated IR/FX models
        environment
        correlationMatrix
        timeGrid
        nPath
        %vector of the handle name of IR models that needs to be simulated
        %eg. ['USD_G2', 'GBP_G2']
        simulationTarget
        %vector of number of simulated factor for each model
        %ag. [2, 3]
        nSimulationFactor
        %simulated random number
        %simulationResult
        nonSimulatedModels
    end
    
    methods
        
        %constructor
        function obj = SimulationEngine(env, corrMatrix, timeGrid, nPath, nSimulationFactor, simulationTarget,nonSimulatedModels)
            obj.environment = env;
            obj.correlationMatrix = corrMatrix;
            obj.timeGrid = timeGrid;
            obj.nPath = nPath;  
            obj.nSimulationFactor =  nSimulationFactor; 
            obj.simulationTarget =  simulationTarget; 
            obj.nonSimulatedModels = nonSimulatedModels;
            
            if size(corrMatrix, 1) ~= size(corrMatrix, 2) || size(corrMatrix, 1) ~= sum(nSimulationFactor)
                error(message('correlation matrix size doesnt match the factors to be simulated'));
            end
        end
        
        %simualte all required random numbers
        function obj = simulate(obj)
            dw = zeros(length(obj.timeGrid), obj.nPath, sum(obj.nSimulationFactor));
            
            %simulate all random numbers
            s = RandStream('mt19937ar','Seed',0);
            RandStream.setGlobalStream(s);
            r = mvnrnd(zeros(1, sum(obj.nSimulationFactor)), obj.correlationMatrix, length(obj.timeGrid) * obj.nPath);
             
            %reshape to dw
            for iFactor = 1 : sum(obj.nSimulationFactor)
                series = r(:, iFactor);
                series = reshape(series, length(obj.timeGrid), obj.nPath);
                dw(:, :,iFactor) = series;
            end
            
            %simulate each model
            factorIndex = 1;
            for iModel = 1 : length(obj.simulationTarget)
                %get model handle fro environment
                irModel = getIRModelHandle(obj.environment, obj.simulationTarget{iModel});
                generateScenario(irModel, obj.timeGrid, dw(:, :, factorIndex: factorIndex + obj.nSimulationFactor(iModel)-1));
                factorIndex = factorIndex + obj.nSimulationFactor(iModel);
            end
            if ~isempty(obj.nonSimulatedModels)
                for iModel = 1 : length(obj.nonSimulatedModels)
                    irModel = getIRModelHandle(obj.environment, char(obj.nonSimulatedModels{iModel}));
                    generateScenario(irModel, obj.timeGrid, []);
                end
            end
        end
        
    end
    
end

