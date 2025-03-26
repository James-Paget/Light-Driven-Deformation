n_medium = 1.33;        % Water
n_particle = 1.59;      % Polystyrene
wavelength0 = 1064e-9;  % Vacuum wavelength
wavelength_medium = wavelength0 / n_medium;
NA = 1.02;              % Numerical aperture of beam
nrel = n_particle/n_medium;
beamPos = [0; 0; 0];
plotSize = 4e-6;

radius = 1e-6;
pos1 = [ radius; 0; 0];

close all;

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

%{
disp("trapezium 10")
CubeN = 25;
Ftot = ForceCalc(pos1, radius, field, CubeN);
disp(Ftot)
disp("simpson 10")
Ftot = ForceCalc_simpson(pos1, radius, field, CubeN);
disp(Ftot)
%}

%{
vals = [2,6,10];
Fs = zeros(length(vals));

for i = vals
    disp("START "+i)
    Ftot = ForceCalc(pos1, radius, field, i);
    disp(Ftot)
    % disp("simpson 10")
    Ftot = ForceCalc_simpson(pos1, radius, field, i);
    disp(Ftot)
end

ott.forcetorque(beam, T)
% plot(log(vals), Fs)
%}


tic
resultant_force = ForceCalc_quad(pos1, radius, field);
disp("resultant force is ")
disp(resultant_force)
toc

%{

num = 20;
positions = linspace(0.8*radius, 1.2*radius, num);
forcesX = linspace(0,0,num);
forcesY = linspace(0,0,num);
forcesZ = linspace(0,0,num);
for i = 1:num
    pos = [positions(i);0;0];
    resultant_force = ForceCalc_quad(pos, radius, field);
    forcesX(i) = resultant_force(1);
    forcesY(i) = resultant_force(2);
    forcesZ(i) = resultant_force(3);
end

f1 = figure;
hold on;
plot(positions, log10(forcesX))
plot(positions, log10(forcesY))
plot(positions, log10(forcesZ))

%}