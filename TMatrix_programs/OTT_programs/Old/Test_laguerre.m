n_medium = 1.33;        % Water
n_particle = 1.59;      % Polystyrene
wavelength0 = 1064e-9;  % Vacuum wavelength
wavelength_medium = wavelength0 / n_medium;
radius = 1.0*wavelength_medium;
NA = 1.02;              % Numerical aperture of beam
beam = ott.BscPmGauss('lg', [ 0 3 ], ...
     'polarisation', [ 1 1i ], 'NA', NA, ...
     'index_medium', n_medium, 'wavelength0', wavelength0);
beam.power = 1.0;
beam.basis = 'regular';
figure();
subplot(1, 2, 1);
beam.visualise('axis', 'y');
subplot(1, 2, 2);
beam.visualise('axis', 'z');