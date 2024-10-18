n_medium = 1.33;        % Water
n_particle = 1.59;      % Polystyrene
wavelength0 = 1064e-9;  % Vacuum wavelength
wavelength_medium = wavelength0 / n_medium;
NA = 1.02;              % Numerical aperture of beam
nrel = n_particle/n_medium;
beamPos = [0; 0; 0];
plotSize = 4e-6;
radius = 1e-6;

close all;
pos1 = [ radius*0.9999999; 0; 0];

shape1 = ott.shapes.Shape.simple("sphere", radius);
shape1.position = pos1;
T = ott.TmatrixMie.simple(shape1, 'index_relative', nrel, ...
    'index_medium', 1.0, 'wavelength0', wavelength0);

beam = ott.BscPlane(0,0, ...
     'polarisation', [ 1 1i ], ...
     'index_medium', n_medium, 'wavelength0', wavelength0);
beam.basis = 'regular';

sbeam = T * beam;
field = sbeam.scatteredField;

for iter = 0:5
    samples = 100+iter;
    Ftot = ForceCalc_modified(pos1, radius, field, samples);
    disp("=== "+iter+" ===");
    disp(Ftot);
end