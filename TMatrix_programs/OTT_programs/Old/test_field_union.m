% Set up constants.
n_medium = 1.33;        % Water
n_particle = 1.59;      % Polystyrene
wavelength0 = 1064e-9;  % Vacuum wavelength
wavelength_medium = wavelength0 / n_medium;
NA = 1.02;              % Numerical aperture of beam
nrel = n_particle/n_medium;
plotSize = 4e-6;
close all;

% Make beam.
beamPos = [0; 0; 0];
beam = ott.BscPlane(0,0, ...
     'polarisation', [ 1 1i ], ...
     'index_medium', n_medium, 'wavelength0', wavelength0);
beam.basis = 'regular';


radius = 10e-6;
pos1 = [-4*radius; 0; 0];
pos2 = [-2*radius; 0; 0];
pos3 = [ 0*radius; 0; 0];
pos4 = [ 2*radius; 0; 0];
pos5 = [ 4*radius; 0; 0];

% make shapes
emptyShape = ott.shapes.Shape.simple("sphere", 0.0);

shape1 = ott.shapes.Shape.simple("sphere", radius);
shape1.position = pos1;
shape2 = ott.shapes.Shape.simple("sphere", radius);
shape2.position = pos2;
shape3 = ott.shapes.Shape.simple("sphere", radius);
shape3.position = pos3;
shape4 = ott.shapes.Shape.simple("sphere", radius);
shape4.position = pos4;
shape5 = ott.shapes.Shape.simple("sphere", radius);
shape5.position = pos5;

%{
unionL = ott.shapes.Union([shape1, shape3, shape5]);
unionR = ott.shapes.Union([shape1, shape2, shape5]);
%}

unionL = ott.shapes.Union([shape1]);
unionR = ott.shapes.Union([shape5]);

unionT = ott.shapes.Union([unionL, unionR]);

% I want to know if the union is different from the sum of the individuals.

z_plane = 0e-6;
plotSize = 4e-6;
plotN = 500;
shouldLogarithm = false;
shouldPlot = true;



% make T matrices
Tempty = ott.TmatrixMie.simple(emptyShape, 'index_relative', nrel, ...
    'index_medium', 1.0, 'wavelength0', wavelength0);
TL = ott.TmatrixMie.simple(unionL, 'index_relative', nrel, ...
    'index_medium', 1.0, 'wavelength0', wavelength0);
TR = ott.TmatrixMie.simple(unionR, 'index_relative', nrel, ...
    'index_medium', 1.0, 'wavelength0', wavelength0);
TT = ott.TmatrixMie.simple(unionT, 'index_relative', nrel, ...
    'index_medium', 1.0, 'wavelength0', wavelength0);

% TT = TL + TR;


sempty = Tempty * beam;
sbeamL = TL * beam;
sbeamR = TR * beam;
sbeamT = TT * beam;

% ZcolEmpty is just always 0.
% [~, ~, ZcolEmpty] = PlotField(sempty.scatteredField, z_plane, plotSize, plotN, shouldLogarithm, shouldPlot);
[X, Y, ZcolL] = PlotField(sbeamL.scatteredField, z_plane, plotSize, plotN, shouldLogarithm, shouldPlot);
[~, ~, ZcolR] = PlotField(sbeamR.scatteredField, z_plane, plotSize, plotN, shouldLogarithm, shouldPlot);
[~, ~, ZcolT] = PlotField(sbeamT.scatteredField, z_plane, plotSize, plotN, shouldLogarithm, shouldPlot);

% All equal for all spheres of the same radius, different if one of the
% radii is changed.
if ZcolR == ZcolT
    disp("RT equal")
end
if ZcolR == ZcolL
    disp("RL equal")
end
if ZcolT == ZcolL
    disp("LT equal")
end




%{
figure
surf(X, Y, Zcol)
shading interp
% view(0,90)
%}

%{
T = ott.TmatrixMie.simple(union, 'index_relative', nrel, ...
    'index_medium', 1.0, 'wavelength0', wavelength0);

beam = ott.BscPmGauss(...
     'polarisation', [ 1 1i ], 'NA', NA, ...
     'index_medium', n_medium, 'wavelength0', wavelength0);
beam.basis = 'regular';

sbeam = T * beam;
field = sbeam.scatteredField;
[X, Y, Zcol] = PlotField(field, z_plane, plotSize, plotN, shouldLogarithm, shouldPlot);
%}


%{
shouldVisualise = false;
if shouldVisualise
    % Plot the field. | (vecnorm(xyz) < radius)
    figure();
    subplot(1, 2, 1);
    sbeam.basis = 'regular';
    sbeam.visualise('axis', 'y', ...
       'mask', @(xyz) (vecnorm(xyz - pos1) < radius) | (vecnorm(xyz - pos2) < radius) , ...
       'range', [1,1]*plotSize)
    title('Scattered field');

    subplot(1, 2, 2);
    tbeam = sbeam.totalField(beam);
    tbeam.basis = 'regular';
    tbeam.visualise('axis', 'y', ...
        'mask', @(xyz) (vecnorm(xyz - pos1) < radius) | (vecnorm(xyz - pos2) < radius), ...
        'range', [1,1]*plotSize)
    title('Total field');
end




plotBeam = false;
if plotBeam
    beam.basis = 'regular';
    figure();
    subplot(1, 2, 1);
    beam.visualise('axis', 'y');
    subplot(1, 2, 2);
    beam.visualise('axis', 'z');
end
%}