function [Ynm, Ynm_theta, Ynm_phi, Ynm_theta_theta, Ynm_theta_phi,...
    Ynm_phi_phi] = sphericalHarmonicY(n, m, theta, phi)
% Computes spherical harmonics and their derivatives
% just like in mathematica
% m is assumed always positive
%
%  Ynm(theta, phi) = sqrt((2*n+1) / (4*pi)) * sqrt((n-m)! / (n+m)!) ...
%                     * Pnm(cos(theta)) * exp(i*m*phi)
%
% First order derivatives 
%
%  Ynm_theta(theta, phi) = i * m * Ynm(theta, phi)
%  Ynm_phi(theta, phi) = i * m * Ynm(theta, phi)
%
%

%% Check input
if ~isscalar(n) || ~isreal(n) || ~isfinite(n) || n < 0 || n ~= round(n)
    error(message('sphericalHarmonicY:InvalidN'));
end
if ~isscalar(m) || ~isreal(m) || ~isfinite(m) || m < 0 || m > n || m ~= round(m)
    error(message('sphericalHarmonicY:InvalidM'));
end

%% Preallocate
N = length(theta);
Ynm = zeros(N,1);
Ynm_theta = zeros(N,1);
Ynm_phi = zeros(N,1);
Ynm_theta_theta = zeros(N,1);
Ynm_theta_phi = zeros(N,1);
Ynm_phi_phi = zeros(N,1);


%% Trig function
cosTheta = cos(theta);
cotTheta = cot(theta);
cscTheta = csc(theta);


%% Distinguish singular points
iZero = (theta == 0);
iPi = (theta == pi);
iOther = ~iZero & ~iPi;


%% Compute Spherical Harmonic Y
Pn = legendre(n, cosTheta, 'sch')';  % (N x n)
if m == 0
    coeff = sqrt((2 * n + 1) / (4 * pi));
else
    coeff = sqrt((n + 0.5) / (4 * pi)) / (-1)^m;    
end
Pnm = Pn(:, m + 1);
Ynm = coeff * Pnm .* exp(1i * m * phi);

        
%% Compute First Derivatives
if n == 0
    % zero result - Ynm_theta = zeros(N,1);
else
    if n == m
        PnmMinusOne = zeros(size(Ynm));
    else
        PnMinusOne = legendre(n - 1, cosTheta, 'sch')'; % (N x n - 1)
        PnmMinusOne = PnMinusOne(:, m + 1) * sqrt((n - m) / (n + m));
    end
    
    iZero = (theta == 0);
    iPi = (theta == pi);
    iOther = ~iZero & ~iPi;
    
    Ynm_theta(iOther) = coeff .* exp(1i * m * phi(iOther)) .* ...
        (n * cotTheta(iOther) .* Pnm(iOther) - (n + m) * cscTheta(iOther) .* PnmMinusOne(iOther));
    
    if m == 1
        Ynm_theta(iZero) = -exp(1i * phi(iZero)) * sqrt(n * (1 + 2*n) * (1 + n)) / sqrt(16*pi);
        Ynm_theta(iPi) = -(-1)^n * exp(1i * phi(iPi)) * sqrt(n * (1 + 2*n) * (1 + n)) / sqrt(16*pi);
    end % else it's zero    
end
Ynm_phi = 1i * m * Ynm;


%% Compute Second Derivatives
if n == 0
    % zero result - Ynm_theta_theta = zeros(N,1);
else
    if n <= m + 1
        PnmMinusTwo = zeros(size(Ynm));
    else
        PnMinusTwo = legendre(n - 2, cosTheta, 'sch')'; % (N x n - 2)
        PnmMinusTwo = PnMinusTwo(:, m + 1) * sqrt(((n - m)*(n - 1 - m)) / ((n + m)*(n - 1 + m)));
    end
    
    Ynm_theta_theta(iOther) = coeff .* exp(1i * m * phi(iOther)) .* ...
        (n * (n * cotTheta(iOther).^2 - cscTheta(iOther).^2) .* Pnm(iOther) -...
        2 * (n - 1) * (m + n) * cotTheta(iOther) .* cscTheta(iOther) .* PnmMinusOne(iOther) + ...
        (n - 1 + m) * (m + n) * cscTheta(iOther).^2 .* PnmMinusTwo(iOther));
     
    if m == 0
        Ynm_theta_theta(iZero) = - sqrt(n^2 * (1 + 2*n) * (1 + n)^2 / (16*pi));
        Ynm_theta_theta(iPi) = -(-1)^n * sqrt(n^2 * (1 + 2*n) * (1 + n)^2 / (16*pi));
    elseif m == 2
        Ynm_theta_theta(iZero) = -exp(2 * 1i * phi(iZero)) * sqrt(n * (1 + 2*n) * (1 + n) * (2 + n) * (n - 1) / (64*pi));
        Ynm_theta_theta(iPi) = (-1)^n * exp(2 * 1i * phi(iPi)) * sqrt(n * (1 + 2*n) * (1 + n) * (2 + n) * (n - 1) / (64*pi));        
    end % else it's zero
end
Ynm_theta_phi = 1i * m * Ynm_theta;
Ynm_phi_phi = -m^2 * Ynm;
