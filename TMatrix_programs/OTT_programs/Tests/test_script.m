% Wavelength in medium/vacuum [m]
n_medium = 1.33;        % Water
n_particle = 1.59;      % Polystyrene
wavelength = 1064e-9;  % Vacuum wavelength
wavelength_medium = wavelength / n_medium;
radius = 1.0*wavelength_medium;
NA = 1.02;              % Numerical aperture of beam
% Setup the T-matrix object
T = ott.Tmatrix.simple('sphere', wavelength, 'index_medium', 1.0, ...
    'index_particle', 1.2, 'wavelength0', wavelength);
% Leg beam
beam = ott.BscPmGauss('hg', [ 0 3 ], ...
     'polarisation', [ 1 1i ], 'NA', NA, ...
     'index_medium', n_medium, 'wavelength0', wavelength);

shape1 = ott.shapes.Sphere(1.0, [0, 0, -2]);
methods(T);
disp("========");
methods(beam);

figure();
beam.visualise('axis', 'z');
title('Single Gaussian beam');

new_beam = translateXyz(beam, [2*wavelength;0;0]);

figure();
new_beam.visualise('axis', 'z');
title("NEW BEAM");

scattered_beam = T*beam;

figure();
scattered_beam.visualise('axis', 'z');
title("SCATTERED BEAM");