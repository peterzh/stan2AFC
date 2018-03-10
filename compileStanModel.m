function compileStanModel(stanFile)
    sm = StanModel('file',stanFile);
    sm.compile();
    save(sm.model_name,'sm');
end