function data_rev = reverseSplitting(data)
%% Reverse the train-test splitting
data_rev = data;

npart = length(data_rev);
for ip = 1:npart
    centmp = data_rev{ip}.centest;
    data_rev{ip}.centest = data_rev{ip}.censor;
    data_rev{ip}.censor = centmp;
    data_rev{ip}.nTrain = length(data_rev{ip}.censor);
    data_rev{ip}.nTest = length(data_rev{ip}.centest);
end


end