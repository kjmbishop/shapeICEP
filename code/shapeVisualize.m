function hParticle = shapeVisualize(geometry)

%% Colors

%% Triangulation
tri = triangulation(geometry.triConnect, ...
    geometry.r(:,1), geometry.r(:,2), geometry.r(:,3));
hParticle = trisurf(tri, geometry.rMag);

%% Axes and Figure
axis equal;
axis on;
box on;
colormap default; 
caxis( [min(geometry.rMag(:)), max(geometry.rMag(:))] );
set(gcf,'Color','w');


%% Lighting and View
view(120,30);
camlight('right')
lighting gouraud
material dull
shading interp

hParticle.AmbientStrength = 0.5;

