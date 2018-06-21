function resFull = runRealworld(dataFull,resFull,nPart,options,is,idataset)

for ip = 1:nPart
    fprintf('Dataset %d\t Part %d\t Sample %d\n',idataset,ip,is);
    data = dataFull{ip};
    data.nU = data.nTrain;
    %% censored data with GP
    nMethod = 1;
    options.inferType = 3;
    % assume the same intensity function
    options.wModel = 0;
    [model,nTime] = varTrainIntCen(data,options);
    res = testmodel(data,model,options);
    res.nTime = nTime;
    displayResult(res);
    resFull{ip}{nMethod,is} = res;
    %% censored data with GP and adding additional weight
    nMethod = 2;
    options.wModel = 1;
    options.wLow = 1e-6;
    options.inferType = 3;
    [model,nTime] = varTrainIntCen(data,options);
    res = testmodel(data,model,options);
    res.nTime = nTime;
    displayResult(res);
    resFull{ip}{nMethod,is} = res;
    %% Benchmark -- localEM
    nMethod = 3;
    [model,nTime] = localEM(data,options);
    res = testmodel(data,model,options);
    res.nTime = nTime;
    displayResult(res);
    resFull{ip}{nMethod,is} = res;
end

end