function displayResult(resTest)
% Display the result according to whether it is GP or Benchmark

if isfield(resTest,'nSampLike')     % variational method
    if isfield(resTest,'nSampMISE') % synthetic set
        fprintf('MISE %.4f\tTestLike %.4f\tSampleMISE %.4f\tSampleLike %.4f\tTime %.4f\n',...
            resTest.nMISE,resTest.nTestLike,resTest.nSampMISE,resTest.nSampLike,resTest.nTime);
    else                            % realworld set
        fprintf('TestLike %.4f\tSampleLike %.4f\tTime %.4f\n',resTest.nTestLike,resTest.nSampLike,resTest.nTime);
    end
else                                % localEM
    if isfield(resTest,'nMISE')     % synthetic set
        fprintf('MISE %.4f\tTestLike %.4f\tTime %.4f\n',resTest.nMISE,resTest.nTestLike,resTest.nTime);
    else                            % realworld set
        fprintf('TestLike %.4f\tTime %.4f\n',resTest.nTestLike,resTest.nTime);
    end
end
end