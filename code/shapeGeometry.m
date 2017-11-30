function output = shapeGeometry(input)
%% output = shapeGeometry(input)
% 
% Computes geometric quantities for a particle. 
% 
% input structure contains:
%   order           - order of Lebedev quadrature
%   c               - spheroid aspect ratio (c = height / width)
%
% output structure contains:
%   N               - number of Lebedev points (1x1)
%   weight          - Lebedev weights (Nx1)
%   theta           - polar angle of Lebedev points (Nx1)
%   phi             - azimuthal angle of Lebedev point (Nx1)
%   r               - position of Lebedev points (Nx3)
%   areaElement     - area element dS of Lebedev points (Nx1)
%   surfaceArea     - total surface area (1x1)
%   hTheta          - scale factor in theta-direction (Nx1)
%   hPhi            - scale factor in phi-direction (Nx1)
%   meanCurvature   - mean curvature (Nx1)
%   n               - unit normal vector at mesh points (Nx3)
%   tPhi            - unit tangent vector (phi-direction) (Nx3)
%   tTheta          - unit tangent vector (theta-direction) (Nx3)
%   rpq             - |rq - rp|  distance matrix (NxN)
%   xpq             - xq - xp  (NxN)
%   ypq             - yq - yp  (NxN)
%   zpq             - zq - zp  (NxN)
%   nqKrpqK         - nqk * rpqk (NxN)
%   t1pKnqK         -
%   t1pKrpqK        -
%   t2pKnqK         -   
%   t2pKrpqK        - 


%% Parse Input
order = 95;    % order of Lebedev quadrature
shapeA = [1,1,1];  % shape parameters, A (3x1)
shapeB = input.shapeB;  % shape parameters, B (sparse)
rootPath = input.rootPath;


%% Load points and weights
filename = sprintf('%s/data/lebedev/lebedev_%03d.txt', rootPath, order);
lebedev = load(filename);

theta = lebedev(:,2) * pi / 180;
phi = lebedev(:,1) * pi / 180;
weight = lebedev(:,3) * 4 * pi;
N = length(theta);


%% Evaluate Surface Points and Derivatives
[n,m,Bnm] = find(shapeB); m = m - 1;

% Spherical harmonics
Ynm = zeros(N,length(n));
Ynm_theta = zeros(N,length(n));
Ynm_phi = zeros(N,length(n));

Ynm_theta_theta = zeros(N,length(n));
Ynm_theta_phi = zeros(N,length(n));
Ynm_phi_phi = zeros(N,length(n));

% Shape functions
f = ones(N,1);
f_theta = zeros(N,1);
f_phi = zeros(N,1);

f_theta_theta = zeros(N,1);
f_theta_phi = zeros(N,1);
f_phi_phi = zeros(N,1);

for i = 1:length(n)    
    [Ynm(:,i), Ynm_theta(:,i), Ynm_phi(:,i), Ynm_theta_theta(:,i), ...
        Ynm_theta_phi(:,i), Ynm_phi_phi(:,i)] = ...
        sphericalHarmonicY(n(i), m(i), theta, phi);

    f = f + real(Bnm(i) * Ynm(:,i));
    f_theta = f_theta + real(Bnm(i) * Ynm_theta(:,i));
    f_phi = f_phi + real(Bnm(i) * Ynm_phi(:,i));
    f_theta_theta = f_theta_theta + real(Bnm(i) * Ynm_theta_theta(:,i));
    f_theta_phi = f_theta_phi + real(Bnm(i) * Ynm_theta_phi(:,i));
    f_phi_phi = f_phi_phi + real(Bnm(i) * Ynm_phi_phi(:,i));
end

% Lebedev points
r = [...
    shapeA(1) * f .* cos(phi) .* sin(theta), ...
    shapeA(2) * f .* sin(phi) .* sin(theta), ...
    shapeA(3) * f .* cos(theta)];
rMag = sqrt(dot(r,r,2));

% First order derivatives at Lebedev points
r_theta = [...
    shapeA(1) * cos(phi) .* (f .* cos(theta) + f_theta .* sin(theta)), ...
    shapeA(2) * sin(phi) .* (f .* cos(theta) + f_theta .* sin(theta)), ...
    shapeA(3) * (-f .* sin(theta) + f_theta .* cos(theta))];
r_phi = [...
    shapeA(1) * sin(theta) .* (-f .* sin(phi) + f_phi .* cos(phi)), ...
    shapeA(2) * sin(theta) .* (f .* cos(phi) + f_phi .* sin(phi)), ...
    shapeA(3) * cos(theta) .* f_phi];

% Second order derivatives at Lebedev points
r_theta_theta = [...
    shapeA(1) * cos(phi) .* (-f .* sin(theta) + 2 * f_theta .* cos(theta) + f_theta_theta .* sin(theta)), ...
    shapeA(2) * sin(phi) .* (-f .* sin(theta) + 2 * f_theta .* cos(theta) + f_theta_theta .* sin(theta)), ...
    shapeA(3) * (-f .* cos(theta) - 2 * f_theta .* sin(theta) + f_theta_theta .* cos(theta))];
r_theta_phi = [...
    shapeA(1) * (-sin(phi) .* (f .* cos(theta) + f_theta .* sin(theta)) + cos(phi) .* (f_phi .* cos(theta) + f_theta_phi .* sin(theta))), ...
    shapeA(2) * (cos(phi) .* (f .* cos(theta) + f_theta .* sin(theta)) + sin(phi) .* (f_phi .* cos(theta) + f_theta_phi .* sin(theta))), ...
    shapeA(3) * (-f_phi .* sin(theta) + f_theta_phi .* cos(theta))];
r_phi_phi = [...
    shapeA(1) * sin(theta) .* (-f .* cos(phi) - 2 * f_phi .* sin(phi) + f_phi_phi .* cos(phi)), ...
    shapeA(2) * sin(theta) .* (-f .* sin(phi) + 2 * f_phi .* cos(phi) + f_phi_phi .* sin(phi)), ...
    shapeA(3) * cos(theta) .* f_phi_phi];

       
%% Correct for singularities at theta = 0, pi
iZero = (theta == 0);
iPi = (theta == pi);
phiZero = phi(iZero) + pi/2;
phiPi = phi(iPi) + pi/2;

YnmZero = zeros(1, length(n));
YnmZero_theta = zeros(1,length(n));
YnmZero_phi = zeros(1,length(n));
YnmZero_theta_theta = zeros(1,length(n));
YnmZero_theta_phi = zeros(1,length(n));
YnmZero_phi_phi = zeros(1,length(n));

YnmPi = zeros(1, length(n));
YnmPi_theta = zeros(1,length(n));
YnmPi_phi = zeros(1,length(n));
YnmPi_theta_theta = zeros(1,length(n));
YnmPi_theta_phi = zeros(1,length(n));
YnmPi_phi_phi = zeros(1,length(n));

fZero = 1;
fZero_theta = 0;
fZero_phi = 0;
fZero_theta_theta = 0;
fZero_theta_phi = 0;
fZero_phi_phi = 0;

fPi = 1;
fPi_theta = 0;
fPi_phi = 0;
fPi_theta_theta = 0;
fPi_theta_phi = 0;
fPi_phi_phi = 0;

for i = 1:length(n)    
    [YnmZero(:,i), YnmZero_theta(:,i), YnmZero_phi(:,i), ...
         YnmZero_theta_theta(:,i), YnmZero_theta_phi(:,i), ...
         YnmZero_phi_phi(:,i)] = sphericalHarmonicY(n(i), m(i), ...
         theta(iZero), phiZero);
    [YnmPi(:,i), YnmPi_theta(:,i), YnmPi_phi(:,i), ...
        YnmPi_theta_theta(:,i), YnmPi_theta_phi(:,i), ...
        YnmPi_phi_phi(:,i)] = sphericalHarmonicY(n(i), m(i), ...
        theta(iPi), phiPi);

    fZero = fZero + real(Bnm(i) * YnmZero(:,i));
    fZero_theta = fZero_theta + real(Bnm(i) * YnmZero_theta(:,i));
    fZero_phi = fZero_phi + real(Bnm(i) * YnmZero_phi(:,i));
    fZero_theta_theta = fZero_theta_theta + real(Bnm(i) * YnmZero_theta_theta(:,i));
    fZero_theta_phi = fZero_theta_phi + real(Bnm(i) * YnmZero_theta_phi(:,i));
    fZero_phi_phi = fZero_phi_phi + real(Bnm(i) * YnmZero_phi_phi(:,i));
    
    fPi = fPi + real(Bnm(i) * YnmPi(:,i));
    fPi_theta = fPi_theta + real(Bnm(i) * YnmPi_theta(:,i));
    fPi_phi = fPi_phi + real(Bnm(i) * YnmPi_phi(:,i));
    fPi_theta_theta = fPi_theta_theta + real(Bnm(i) * YnmPi_theta_theta(:,i));
    fPi_theta_phi = fPi_theta_phi + real(Bnm(i) * YnmPi_theta_phi(:,i));
    fPi_phi_phi = fPi_phi_phi + real(Bnm(i) * YnmPi_phi_phi(:,i));
end

% First order derivatives at Lebedev points
r_phi(iZero,:) = [...
    shapeA(1) * cos(phiZero) * fZero, ...
    shapeA(2) * sin(phiZero) * fZero, ...
    shapeA(3) * fZero_theta];
r_phi(iPi,:) = [...
    shapeA(1) * cos(phiPi) * fPi, ...
    shapeA(2) * sin(phiPi) * fPi, ...
    shapeA(3) * fPi_theta];

% Second order derivatives at Lebedev points
r_theta_phi(iZero,:) = [0,0,0];
r_theta_phi(iPi,:) = [0,0,0];

r_phi_phi(iZero,:) = [ ...
    shapeA(1) * cos(phiZero) * (2 * fZero_theta), ...
    shapeA(2) * sin(phiZero) * (2 * fZero_theta), ...
    shapeA(3) * (-fZero + fZero_theta_theta)];
r_phi_phi(iPi,:) = [ ...
    shapeA(1) * cos(phiPi) * (-2 * fPi_theta), ...
    shapeA(2) * sin(phiPi) * (-2 * fPi_theta), ...
    shapeA(3) * (fPi - fPi_theta_theta)];


%% Scale factors (N x 1) and Surface tangent vectors (N x 3)
hTheta = sqrt(dot(r_theta, r_theta, 2));
hPhi = sqrt(dot(r_phi, r_phi, 2));

tTheta = bsxfun(@rdivide, r_theta, hTheta);
tPhi = bsxfun(@rdivide, r_phi, hPhi);


%% Surface normal vector (N x 3) and area element (N x 1)
n = cross(r_theta, r_phi, 2);
areaElement = sqrt(dot(n, n, 2));
n = bsxfun(@rdivide, n, areaElement);
areaElement = areaElement ./ sin(theta);

% Correct area elements at singular points
r2Zero = r(iZero,:) * r(iZero,:)';
dr2Zero = r_phi(iZero,3).^2 + r_theta(iZero,3)^2;
areaElement(iZero) = r2Zero * sqrt(1 + dr2Zero/r2Zero);

r2Pi = r(iPi,:) * r(iPi,:)';
dr2Pi = r_phi(iPi,3)^2 + r_theta(iPi,3)^2;
areaElement(iPi) = r2Pi * sqrt(1 + dr2Pi/r2Pi);


%% Mean Curvature
% First fundamental form
EE = dot(r_theta, r_theta, 2);
FF = dot(r_theta, r_phi, 2);
GG = dot(r_phi, r_phi, 2);

% Second fundamental form
LL = dot(r_theta_theta, n, 2);
MM = dot(r_theta_phi, n, 2);
NN = dot(r_phi_phi, n, 2);

meanCurvature = zeros(N,1);
principalCurvature = zeros(N,2);
for i = 1:N    
    F1 = [EE(i), FF(i); FF(i), GG(i)];
    F2 = [LL(i), MM(i); MM(i), NN(i)];
    A = inv(F1) * F2;
    meanCurvature(i) = 0.5 * trace(A);
    principalCurvature(i,:) = eig(A);
end


%% Calculate surface area
surfaceArea = sum(areaElement .* weight);


%% Distance matrix (N x N)
% p,q indices refer to pairs of points on the particle surface
xpq = bsxfun(@minus, r(:,1)', r(:,1)); % xq - xp
ypq = bsxfun(@minus, r(:,2)', r(:,2)); % yq - yp
zpq = bsxfun(@minus, r(:,3)', r(:,3)); % zq - zp

rpq = sqrt(xpq.^2 + ypq.^2 + zpq.^2); % scalar distance between points p and q
rpq(logical(speye(N))) = Inf; % set diagonals (p,p) to infinity

% normalize
xpq = xpq ./ rpq;
ypq = ypq ./ rpq;
zpq = zpq ./ rpq;

% nqK * rpqK (contracting over capital index)
nqKrpqK = bsxfun(@times, n(:,1)', xpq) ...
    + bsxfun(@times, n(:,2)', ypq) ...
    + bsxfun(@times, n(:,3)', zpq);

% t1pK * nqK (theta direction)
t1pKnqK = bsxfun(@times, tTheta(:,1), n(:,1)') ...
    + bsxfun(@times, tTheta(:,2), n(:,2)') ...
    + bsxfun(@times, tTheta(:,3), n(:,3)');

% t1pK * rpqK (theta direction)
t1pKrpqK = bsxfun(@times, tTheta(:,1), xpq) ...
    + bsxfun(@times, tTheta(:,2), ypq) ...
    + bsxfun(@times, tTheta(:,3), zpq);

% t2pK * nqK (phi direction)
t2pKnqK = bsxfun(@times, tPhi(:,1), n(:,1)') ...
    + bsxfun(@times, tPhi(:,2), n(:,2)') ...
    + bsxfun(@times, tPhi(:,3), n(:,3)');

% t2pk * rpqk (phi direction)
t2pKrpqK = bsxfun(@times, tPhi(:,1), xpq) ...
    + bsxfun(@times, tPhi(:,2), ypq) ...
    + bsxfun(@times, tPhi(:,3), zpq);


%% Triangulation (used for visualization with trisurf)
thetaExtended = [theta; theta; theta];
phiExtended = [phi; phi + 2*pi; phi - 2*pi];
DT = delaunayTriangulation([thetaExtended, phiExtended]);
idx = (DT.ConnectivityList(:,1) <= N) ...
    | (DT.ConnectivityList(:,2) <= N) ...
    | (DT.ConnectivityList(:,3) <= N);
triConnect = 1 + mod(DT.ConnectivityList(idx,:) - 1, N);


%% Prepare Output
output = input;
output.N = N;           % Number of mesh points (1x1)
output.weight = weight; % Lebedev weights (Nx1)
output.theta = theta;
output.phi = phi;
output.r = r;           % Position of mesh points (Nx3)
output.areaElement = areaElement; % dS (Nx1)
output.surfaceArea = surfaceArea; % surface area (1x1)
output.meanCurvature = meanCurvature; % (Nx1)
output.principalCurvature = principalCurvature; % (Nx1)
output.n = n;           % Unit normal vector at mesh points (Nx3)
output.tPhi = tPhi;     % Unit tangent vector (phi-direction) (Nx3)
output.tTheta = tTheta; % Unit tangent vector (theta-direction) (Nx3)
output.rpq = rpq;       % |rq - rp|  distance matrix (NxN)
output.xpq = xpq;       % xq - xp  (NxN)
output.ypq = ypq;       % yq - yp  (NxN)
output.zpq = zpq;       % zq - zp  (NxN)
output.nqKrpqK = nqKrpqK;     % nqk * rpqk (NxN)
output.t1pKnqK = t1pKnqK;  
output.t1pKrpqK = t1pKrpqK;  
output.t2pKnqK = t2pKnqK;   
output.t2pKrpqK = t2pKrpqK; 
output.triConnect = triConnect;
output.rMag = rMag;

% output.hTheta = hTheta; % scale factor in theta-direction (Nx1)
% output.hPhi = hPhi;     % scale factor in phi-direction (Nx1)

