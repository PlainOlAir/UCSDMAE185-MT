function mu = sutherland(T, mu0, T0, S1)
%% Description
% Calculates dynamic viscoity using Sutherland's law (PRESET AIR)
% Assumes calorically perfect gas

% INPUTS
% T - temperature with size Nx x Ny
% mu0 - dynamic viscosity at sea-level (OPTIONAL)
% T0 - temperature at sea-level (OPTIONAL)
% S1 - Sutherland temperature at sea-level (OPTIONAL)

% OUTPUTS
% mu - dynamic viscosity with size Nx x Ny

%% Setup
    % process input variables
    switch nargin
        case 1
            mu0 = 1.735e-5; %N-s/m^2
            T0 = 288.15; %K
            S1 = 110.4; %K
        case 2
            T0 = 288.15; %K
            S1 = 110.4; %K
        case 3
            S1 = 110.4; %K
    end

%% Calculate
    % Sutherland's law
    mu = mu0 .* (T ./ T0) .^ (3/2) .* (T0 + S1) ./ (T + S1);

end