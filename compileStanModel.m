function compileStanModel(stanFile)
    baseDir = fileparts(mfilename('fullpath'));
    binDir = fullfile(baseDir,'bin');
    if ~exist(binDir,'dir')
        mkdir(binDir);
    end
    
    %Copy stanFile to binary directory
    copyfile(fullfile(baseDir,stanFile), fullfile(binDir,stanFile));
    
    %Compile model
    sm = StanModel('file',fullfile(binDir,stanFile));
    sm.compile();
    
    %Save copy of stanmodel object
    save(sm.model_name,'sm');
end