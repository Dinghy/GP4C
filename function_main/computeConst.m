function nA = computeConst(vPoint)
%% compute the constant part in the likelihood of a Poisson likelihood
% nA = \sum_i \ln (m_i)!

nA = sum(arrayfun(@(x)sum(log(1:x)),vPoint));
end