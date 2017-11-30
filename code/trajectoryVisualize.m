function trajectoryVisualize(sol, geometry, varargin)

if nargin == 2
    %% Input Parameters
    axisLength = 2;
    axesWidth = 2;
    ribbonPoint = 10;
    ribbonFrac = 0.3;
    timeStep = 0.5;
    transparency = 0.6;
    startTransparency = 0.3;
    ambientStrength = 0.5;
else
    %% Input Parameters
    axisLength = varargin{1}.axisLength;
    axesWidth = varargin{1}.axesWidth;
    ribbonPoint = varargin{1}.ribbonPoint;
    ribbonFrac = varargin{1}.ribbonFrac;
    timeStep = varargin{1}.timeStep;
    transparency = varargin{1}.transparency;
    startTransparency = varargin{1}.startTransparency;
    ambientStrength = varargin{1}.ambientStrength;
end

%% Colors
axisXColor = [0.9,0,0];
axisYColor = [0,0.9,0];
axisZColor = [0,0,0.9];
startColor = 1 * ones(1,3);


%% Generate equally spaced points in time
t = linspace(0, sol.x(end), sol.x(end)/timeStep);
y = interp1(sol.x, sol.y', t, 'linear')';
nPoint = size(y,2);


%% Rotation matrices
rotationMatrix = quaternionToRotationMatrix( y(4:7,:)' );


%% Particle Axes Ribbons
axisX = zeros(3, ribbonPoint * nPoint);
axisY = zeros(3, ribbonPoint * nPoint);
axisZ = zeros(3, ribbonPoint * nPoint);

% ribbon point connectivity
[X,Y] = meshgrid(1:nPoint, 1:ribbonPoint);
DT = delaunayTriangulation(X(:), Y(:));

for iStep = 1:nPoint
    Rq = reshape(rotationMatrix(iStep,:), 3, 3)';
    
    for jStep = 1:ribbonPoint
        idx = ribbonPoint * (iStep-1) + jStep;
        axisX(:,idx) = y(1:3,iStep) ...
            + ribbonFrac * axisLength * Rq(:,1) * (jStep - 1) / (ribbonPoint - 1);
        axisY(:,idx) = y(1:3,iStep) ...
            + ribbonFrac * axisLength * Rq(:,2) * (jStep - 1) / (ribbonPoint - 1);
        axisZ(:,idx) = y(1:3,iStep) ...
            + ribbonFrac * axisLength * Rq(:,3) * (jStep - 1) / (ribbonPoint - 1);
    end
end

hAxisX = trisurf(DT.ConnectivityList, axisX(1,:), axisX(2,:), axisX(3,:)); hold on;
hAxisY = trisurf(DT.ConnectivityList, axisY(1,:), axisY(2,:), axisY(3,:)); 
hAxisZ = trisurf(DT.ConnectivityList, axisZ(1,:), axisZ(2,:), axisZ(3,:)); 


%% Plot the trajectory of the particle center
idx = [ribbonPoint : ribbonPoint : ribbonPoint * nPoint];
hLines = plot3(...
    [y(1,1), axisX(1,idx)], [y(1,1), axisX(2,idx)], [y(1,1),axisX(3,idx)],...
    [y(1,1), axisY(1,idx)], [y(1,1), axisY(2,idx)], [y(1,1),axisY(3,idx)],...
    [y(1,1), axisZ(1,idx)], [y(1,1), axisZ(2,idx)], [y(1,1),axisZ(3,idx)]);


%% Plot particle axes
axisRatio = 1 / ribbonFrac;
lAxisX = axisRatio * (axisX(1:3,end) - y(1:3,end)) + y(1:3,end);
lAxisY = axisRatio * (axisY(1:3,end) - y(1:3,end)) + y(1:3,end);
lAxisZ = axisRatio * (axisZ(1:3,end) - y(1:3,end)) + y(1:3,end);

hAxesStart = plot3(...
    [lAxisX(1),y(1,end)], [lAxisX(2),y(2,end)], [lAxisX(3),y(3,end)],...
    [lAxisY(1),y(1,end)], [lAxisY(2),y(2,end)], [lAxisY(3),y(3,end)],...
    [lAxisZ(1),y(1,end)], [lAxisZ(2),y(2,end)], [lAxisZ(3),y(3,end)]);


%% Plot the particle at the starting location
geo = geometry;
Rq = reshape(rotationMatrix(1,:), 3, 3)';
geo.r = bsxfun(@plus,(Rq * geo.r'), sol.y(1:3,1))';
hStart = shapeVisualize(geo);
caxis( [min(geo.rMag(:)), max(geo.rMag(:))] );


%% Plot the particle at the final location
geo = geometry;
Rq = reshape(rotationMatrix(end,:), 3, 3)';
geo.r = bsxfun(@plus,(Rq * geo.r'), sol.y(1:3,end))';
hFinish = shapeVisualize(geo);
caxis( [min(geo.rMag(:)), max(geo.rMag(:))] );


%% Axes and View
axis equal
set(gcf,'Color','White');
set(gca,...
    'Color','None',...
    'XTick',[],...
    'YTick',[],...
    'ZTick',[]);
box on;


%% Lighting and View
view(120,30);
camlight('right')
lighting gouraud
material dull
shading interp


%% Color objects
hAxisX.FaceColor = axisXColor;
hAxisY.FaceColor = axisYColor;
hAxisZ.FaceColor = axisZColor;

hAxisX.FaceAlpha = transparency;
hAxisY.FaceAlpha = transparency;
hAxisZ.FaceAlpha = transparency;

hAxisX.EdgeColor = 'none';
hAxisY.EdgeColor = 'none';
hAxisZ.EdgeColor = 'none';

hAxisX.AmbientStrength = ambientStrength;
hAxisY.AmbientStrength = ambientStrength;
hAxisZ.AmbientStrength = ambientStrength;

hLines(1).Color = axisXColor;
hLines(2).Color = axisYColor;
hLines(3).Color = axisZColor;

hAxesStart(1).LineWidth = axesWidth;
hAxesStart(2).LineWidth = axesWidth;
hAxesStart(3).LineWidth = axesWidth;

hAxesStart(1).Color = axisXColor;
hAxesStart(2).Color = axisYColor;
hAxesStart(3).Color = axisZColor;

hStart.FaceColor = startColor;
hStart.FaceAlpha = startTransparency;
hStart.AmbientStrength = ambientStrength;

hFinish.AmbientStrength = ambientStrength;

hold off;

