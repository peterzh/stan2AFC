classdef behavModel
    properties
        modelName;
    end
    
    properties(Access=private)
        stanModelObj;
    end
    
    methods
        function obj = behavModel(modelName)
            obj.modelName = modelName;
            
            %Define directory containing the compiled model
            models_dir = fullfile( fileparts(mfilename('fullpath')) , 'models');
            
            %If compiled model .mat exists, load it
            modelFile = fullfile(models_dir, [modelName '.mat']);
            if exist( modelFile, 'file')
                load(modelFile,'sm');
                fprintf('Loaded model %s\n', modelName);
            else %If it doesn't exist, compile and save it
                fprintf('Compiled model %s not found, compiling now...\n', modelName);
                
                %Load stan model string
                stanFile = strrep(modelFile, '.mat', '.stan');
                
                %Compile model
                sm = StanModel('file',stanFile,'working_dir',tempdir);
                sm.compile();
                
                %Save copy of compiled model .mat object
                save(modelFile,'sm');
            end

        end
        
        function simulateAndFit
            %Simulate data and fit model
        end
        
        function fit
            %Fit model on data
        end
        
        function plotFit
            %Plot psych data fit with uncertainty
        end
    end
    
    
end