function martingaleTest = mtgTest(underlying ,delta_t0, deltas ,N ,rates, ...
    fwd_rates, dividend ,sigma, RD, index, fig_number)
% This function performs the Martingale test for the underlying asset:
%   - check if the selected parameters are consistent with the martingale property
%
% INPUTS:
% underlying            : initial value of the underlying asset
% delta_t0              : time intevals between t0 and ti
% deltas                : time steps
% N                     : number of simulations
% rates                 : interest rates
% fwd_rates             : forward rates
% dividend              : dividend yield
% sigma                 : volatility
% RD                    : regular deduction
% index                 : index of the subplot
% fig_number            : number of the figure
%
% OUTPUTS:
% martingaleTest        : vector of the martingale test

% Rates vector
rates = [1; rates];

% Simulate the underlying asset
simulated_UL = MC_simulation(underlying ,deltas ,N ,fwd_rates, dividend ,sigma, RD);
% Compute the mean of the simulated underlying asset
mean_UL = mean(simulated_UL,1)';
% Compute the flat underlying asset (constant initial value)
flat_UL = underlying * ones(length(delta_t0),1);
% Compute the discounted mean of the simulated underlying asset
DF_UL = exp(-rates.*delta_t0).*mean_UL;

% to check the martingale property we compute the difference between the
% flat underlying asset and the discounted mean of the simulated underlying asset
% hence we accept the martingale property if the difference is smaller than a trheshold
martingaleTest = abs(flat_UL - DF_UL);

% plot Martingale test in a subplot with index i
figure(fig_number)
subplot(2,3,index)
plot(delta_t0, flat_UL,'--r')
hold on
plot(delta_t0, mean_UL, 'k','LineWidth',1)
hold on
plot(delta_t0,DF_UL,'+b')
% tile with number of simulations
title(['Martingale Test with N = ', num2str(N)])
hold off

end