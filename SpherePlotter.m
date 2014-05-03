nParticles=1000;

boundingBoxEdge = 100;

% Possible sphere surface colors. The number of elements in this cell array
% matches the number of possible colors. 
my_colors = {'r' 'b'};
% Colors can also be represented in RGB format.
%my_colors = {[1 0 0] [0 0 1]};



% Create random matrix for color selection of each sphere. The vector "c"
% contains a random number of 1's and 2's with the same number of vector
% elements as x
numColors = length(my_colors); %number of possible colors
c = ceil(numColors * rand(size(z)));

% Create unit sphere
[x1, y1, z1] = sphere(10);
% Scale factor for unit sphere
scalefactor = 0.2;

% Apply scale factor to sphere surface data values
x1 = scalefactor*x1;
y1 = scalefactor*y1;
z1 = scalefactor*z1;

% Create a new figure window with a black background

fig = figure('Color', 'k');

% Create a new axes, which is invisible.

ax = axes('Visible', 'off');

hold(ax, 'on')

% Loop through each data element.

for ind = 1:length(x)
    % Create a new surface object. The surface data is the scaled unit
    % sphere where the center has been translated by the coordinate value
    % of the data point the sphere represents.
    surf(ax, x(ind) + x1, y(ind) + y1, z(ind) + z1, 'FaceColor', my_colors{c(ind)}, 'EdgeColor', 'none')
    % The surface object's FaceColor has been set to a scalar value, and
    % the EdgeColor is transparent ('none').
end

% Locate the handles to all of the surface objects in the current figure.
h = findall(ax, 'Type', 'surface');

% Set the surface objects' lighting properties
set(h,'FaceLighting','phong', 'AmbientStrength',0.5);

% Create a new light to illuminate the spheres.
light('Position',[1 0 0],'Style','infinite');

     
Box([0,boundingBoxEdge,0,boundingBoxEdge,0,boundingBoxEdge,0],Filled:FALSE, LineColor:White);
Axes=None;
Scaling=Constrained;

axis([0,100,0,100,0,100,0,255]);

% Set view angle to default 3-D view.
view(3)
