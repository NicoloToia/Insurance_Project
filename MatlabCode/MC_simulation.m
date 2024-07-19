function MC_value = MC_simulation(S0 ,deltas ,N ,rates, dividend ,sigma, RD)
% The MC_simulation function run a Monte carlo simulation for an asset with
% initial value S0, that follows a GBM dynamics over time

% INPUTS
% S0            : Initial underlying value
% deltas        : delta time for each year in act 365
% N             : number of MC simulations
% rates         : forward rates
% dividend      : dividend yield (per year)
% sigma         : volatility
% RD            : regular deduction

% OUTPUTS
% MC_value      : matrix of the simulated Monte Carlo GBM value

% number of time steps
n_steps = length(deltas);

% extrapolate N random numbers from a normal distribution (standard normal) for n_steps 
g = randn(N,n_steps);

% initialize the matrix of undelying values:
%   -> rows are the paths
%   -> columns are the time step
S = zeros(N,n_steps+1);
S(:,1) = S0*ones(N,1);

% simulate the GBM process
for i = 1:n_steps

    S(:,i+1) = S(:,i) * (1 - RD) .* exp( ( (rates(i+1) - dividend) - 0.5*sigma^2)*...
        deltas(i) + sigma*sqrt(deltas(i))*g(:,i) );

end

% return a matrix of the simulated GBM values
MC_value = S;

end