% Example calculation of forces on a sphere in a optical beam
%
% This example includes a couple of example beams (Gaussian, LG and HG),
% combining the separate example scripts found in ottv1.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Add the toolbox to the path (assuming we are in ott/examples)
addpath('../');

% Close open figures
close all;

%% Describe the particle, beam and surrounding medium

% Make warnings less obtrusive
ott.warning('once');
ott.change_warnings('off');

% Refractive index of particle and medium
n_medium = 1.33;
n_particle = 1.59;
nrel = n_particle/n_medium;

% Wavelength of light in vacuum [m]
wavelength0 = 1064e-9;

% Calculate the wavelength in the medium
wavelength_medium = wavelength0/n_medium;

% Radius of particle
radius = 1.0*wavelength_medium;

% Set the beam type (must be one of lg, gaussian or hg)
% This is used bellow to select from the beam presets
beam_type = 'lg';

% Specify the numerical aparture of the beam/microscope
NA = 1.02;

%% Setup the T-matrix for the particle

tic

% Create a T-matrix for a sphere
%T = ott.Tmatrix.simple('sphere', radius, 'wavelength0', wavelength0, ...
%    'index_medium', n_medium, 'index_particle', n_particle);

%{
. Notes;
--> Tmatrix seems to only like
--> Mie vs regualr
%}
%shape1 = ott.shapes.Sphere(1.0*wavelength0);
particle_radius = 2e-7;
particle_offset = 1e-6;
shape_a = ott.shapes.Shape.simple('sphere',particle_radius);
shape_b = ott.shapes.Shape.simple('sphere',particle_radius);
shape_a.position = [-particle_offset; 0; 0];
shape_b.position = [ particle_offset; 0; 0];
shape_union = ott.shapes.Union([shape_a, shape_b]);
%shape1 = ott.shapes.Shape.simple(shape1, p.Results.parameters);
T = ott.TmatrixMie.simple(shape_union, 'index_relative', nrel, ...
    'index_medium', n_medium, 'wavelength0', wavelength0);

newPos = [1;2;3];
disp( vecnorm(newPos) );
disp(['T-matrix calculation took ' num2str(toc) ' seconds']);

%% Calculate the beam shape coefficients

tic

switch beam_type
  case 'gaussian'

    % Create a simple Gaussian beam with circular polarisation
    % We could also calculate the beam angle ourseves and specify that.
    beam = ott.BscPmGauss('NA', NA, 'polarisation', [ 1 1i ], ...
        'index_medium', n_medium, 'wavelength0', wavelength0);

  case 'lg'

    % Create a LG03 beam with circular polarisation
    beam = ott.BscPmGauss('lg', [ 0 3 ], ...
        'polarisation', [ 1 1i ], 'NA', NA, ...
        'index_medium', n_medium, 'wavelength0', wavelength0);

  case 'hg'

    % Create a HG23 beam with circular polarisation
    beam = ott.BscPmGauss('hg', [ 2 3 ], ...
        'polarisation', [ 1 1i ], 'NA', NA, ...
        'index_medium', n_medium, 'wavelength0', wavelength0);

  otherwise
    error('Unsupported beam type');
end

% Normalize the beam power
beam.power = 1.0;

disp(['Beam calculation took ' num2str(toc) ' seconds']);

%% Generate force/position graphs

tic

% Calculate the force along z
%fXyz = ott.forcetorque(beam, T, 'position', xyz);

%%
%%
% Define ranges for x and y
width = 8e-6;
sampleNumber = 30;
jumpWidth = width/sampleNumber;
x = [1;0;0] .* linspace(-width/2.0, width/2.0, sampleNumber);
y = [0;1;0] .* linspace(-width/2.0, width/2.0, sampleNumber);
z_sample_section = 0;

[X, Y] = meshgrid(x(1,:),y(2,:));
Z = [];
for i = 1: sampleNumber;
    for j = 1: sampleNumber;
        fXyz = ott.forcetorque(beam, T, 'position', [x(1, i); y(2, j); z_sample_section]);
        Z(i, j) = vecnorm(fXyz);
    end
end
% Create a 3D surface plot
figure();
surf(X, Y, Z);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Force Norm');
%%
%%

Z = [];
for i = 1: sampleNumber;
    for j = 1: sampleNumber;
        fXyz = ott.forcetorque(beam, T, 'position', [x(1, i); y(2, j); z_sample_section]);
        Z(i, j) = fXyz(1);
    end
end
% Create a 3D surface plot
figure();
surf(X, Y, Z);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Force X');


Z = [];
for i = 1: sampleNumber;
    for j = 1: sampleNumber;
        fXyz = ott.forcetorque(beam, T, 'position', [x(1, i); y(2, j); z_sample_section]);
        Z(i, j) = fXyz(2);
    end
end
% Create a 3D surface plot
figure();
surf(X, Y, Z);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Force Y');


Z = [];
for i = 1: sampleNumber;
    for j = 1: sampleNumber;
        fXyz = ott.forcetorque(beam, T, 'position', [x(1, i); y(2, j); z_sample_section]);
        Z(i, j) = fXyz(3);
    end
end
% Create a 3D surface plot
figure();
surf(X, Y, Z);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Force Z');


disp(['Force calculation took ' num2str(toc) ' seconds']);

%% Generate the plots

figure();
nPc = n_medium / 3e8;  % n / vacuum_speed
plot(xyz(3, :), fXyz .* nPc);
xlabel('Z position [m]');
ylabel('Force [N/W]');
legend({'Fx', 'Fy', 'Fz'});

figure();
beam.visualise('axis', 'z', ...
    'mask', @(xyz) ((vecnorm(xyz -shape_a.position) < shape_a.radius) | (vecnorm(xyz -shape_b.position) < shape_b.radius)) );
title("Current beam");

sbeam = T*beam;
figure();
subplot(1, 2, 1);
sbeam.basis = 'outgoing';
sbeam.visualise('axis', 'y', ...
   'mask', @(xyz) ((vecnorm(xyz -shape_a.position) < shape_a.radius) | (vecnorm(xyz -shape_b.position) < shape_b.radius)), ...
   'range', [1,1]*2e-6)
title('Scattered field');

subplot(1, 2, 2);
tbeam = sbeam.totalField(beam);
tbeam.basis = 'regular';
tbeam.visualise('axis', 'y', ...
   'mask', @(xyz) ((vecnorm(xyz -shape_a.position) < shape_a.radius) | (vecnorm(xyz -shape_b.position) < shape_b.radius)), ...
   'range', [1,1]*2e-6)
title('Total field');

%'mask', @(xyz) vecnorm([xyz(1)+wavelength0; xyz(2)+wavelength0; xyz(3)+wavelength0]) < 1.0*wavelength0, 
%'mask', @(xyz) vecnorm(xyz) < 1.0*wavelength0, 