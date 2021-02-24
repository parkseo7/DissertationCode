function output = detEstimatesExp(z, taum, Omega, num_trials, parameters)
% Provides estimates of determinant statistics from the i.i.d. sampled
% exponential matrix.
% Returns:
% - Mean of det
% - Variance of det
% - Normality measure using lillietest

% Fixed parameters
g = parameters.g;
omega0 = parameters.omega0;
N = parameters.N;

% Theoretical values
mu = 1 / taum;
mu_0 = mu^2 / (mu^2 + Omega^2);
beta = - g / (z + g*mu_0);

% Set up determinant trials
D = eye(N);

% Arrays
arrDet = zeros(num_trials,1);
arrAvgCos = zeros(num_trials,1);
arrVarCosExp = zeros(num_trials,1);

for i=1:num_trials
    tau0 = exprnd(taum, N);
    X = cos(Omega*tau0).*exp(-z*tau0);
    Y = D + beta*(X/N);
    
    arrDet(i) = det(Y);
    arrAvgCos(i) = mean(cos(Omega*tau0), 'all');
    arrVarCosExp(i) = var(X(:));
    
end

% Compute statistics
[h,p,k,c] = lillietest(log(abs(arrDet)));

output = struct;
output.pvalue = p;
output.detMean = mean(arrDet);
output.detVar = var(arrDet);
output.Omega = Omega;

% Other statistics
output.avgCosVar = var(arrAvgCos);
output.avgVarCosExp = mean(arrVarCosExp);

end

