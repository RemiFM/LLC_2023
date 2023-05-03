% Remi De Coster, May 2023
% https://www.infineon.com/dgdl/Design+Guide+for+LLC+Converter+with+ICE2HS01G+05072011.pdf?fileId=db3a304330f68606013103ebd94f3e98

clear; close all; clc;

%% INPUTS
m = [5 6 7];                % Resonant factors
Q = [4 2 1 0.5 0.4 0.3];    % Quality factors

V_in = 156.8;       % Nominal input voltage (V)
V_in_min = 98.0;      % Minimum input voltage (V)
V_in_max = 178.85;  % Maximum input voltage (V)

V_out = 12;     % Output voltage (V)
I_out = 2;      % Maximum output current (A)

f_r = 100e3;    % Resonant frequency (Hz)
V_f = 0.2;      % Rectification voltage drop (V)

%% CALCULATIONS
% Calculate transformer turns ratio
n = V_in / (2 * (V_out + 2*V_f)); % Division by two due to half-bridge switching

% Calculate worst-case gain
M_max = n * ((V_out + 2*V_f) / (V_in_min/2));
M_min = n * ((V_out + 2*V_f) / (V_in_max/2));

% Calculate effective load resistance
R_eff = (8 / pi^2) * n^2 * (V_out/I_out);

% Calculate normalised frequency
f = 20e3:1e2:180e3;
F = f/f_r;

for j = 1:length(m)
    % Calculate and plot gain magnitude for all Q and m
    Mi = zeros(length(Q), length(F));
    Gi = zeros(length(Q), length(F));
    
    for i = 1:length(Q)
        Mi(i, 1:length(F)) = ((m(j)-1) .* (F.^2)) ./ ( ((F.^2).*m(j) - 1) + ((1i.*F).*(F.^2 - 1).*(m(j)-1).*Q(i)));
        Gi(i, 1:length(F)) = sqrt(real(Mi(i, 1:length(F))).^2 + imag(Mi(i, 1:length(F))).^2);
    end
    figure();
    plot(F, Gi)
    
    grid on;
    yline(M_max, '--', 'Mmax = ' + string(M_max))
    yline(M_min, '--', 'Mmin = ' + string(M_min))
    axis([0.2 1.8 0.2 1.8])
    title('Magnitude vs Frequency for m = ' + string(m(j)))
    legend("Q = " + string(Q))
    xlabel('Normalized frequency F');
    ylabel('Voltage gain G')
end

% Calculate resonant tank for selected Q & m values:
% [C_r, L_r, L_p, L_m] = LLC_Calculate_Components(Q, m, f_r, R_eff)
Q = 0.3;
m = 6;
[C_r, L_r, L_p, L_m] = LLC_Calculate_Components(Q, m, f_r, R_eff)

% % Calculate filter capacitor
Th = 20E-3;     % Hold-up time
% syms C_out
% P_in = 27;  % Input power
% C_out = vpa(solve(V_in_min == sqrt(V_in^2 - ((2*P_in*Th)/C_out)), C_out))

function [C_r, L_r, L_p, L_m] = LLC_Calculate_Components(Q, m, f_r, R_eff)
    C_r = 1 / (2 * pi * Q * f_r * R_eff);
    L_r = 1 / ((2 * pi * f_r)^2 * C_r);
    L_p = m * L_r;
    L_m = L_p - L_r;
end