n_medium = 1.33;        % Water
n_particle = 1.59;      % Polystyrene
wavelength0 = 1064e-9;  % Vacuum wavelength
nrel = n_particle/n_medium;
beamPos = [0; 0; 0];
plotSize = 4e-6;

radius = 0.05e-6;
pos1 = [-4*radius; 0; 0];
pos2 = [-2*radius; 0; 0];
pos3 = [ 0*radius; 0; 0];
pos4 = [ 2*radius; 0; 0];
pos5 = [ 4*radius; 0; 0];

% make shapes
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

unionL = ott.shapes.Union([shape1, shape5]);
unionR = ott.shapes.Union([shape1, shape2,shape3, shape5]);
unionT = ott.shapes.Union([unionL, unionR]);

% I want to know if the union is different from the sum of the individuals.

% make T matrices
TL = ott.TmatrixMie.simple(unionL, 'index_relative', nrel, ...
    'index_medium', 1.0, 'wavelength0', wavelength0);
TR = ott.TmatrixMie.simple(unionR, 'index_relative', nrel, ...
    'index_medium', 1.0, 'wavelength0', wavelength0);
TT = ott.TmatrixMie.simple(unionT, 'index_relative', nrel, ...
    'index_medium', 1.0, 'wavelength0', wavelength0);

% convert sparse matrix to full
L = TL.data;
R = TR.data;
T = TT.data;
TLout = full(TL.data);
TRout = full(TR.data);
TTout = full(TT.data);

% print all the diagonal elements for union, shape1 and shape2 (off
% diagonal all zero)

disp("Each line shows the diagonal elements of the left, right and total T-matrices. They should be different.")
for i = 1:length(L)
    Z1 = TLout(i,i);
    Z2 = TRout(i,i);
    Z3 = TTout(i,i);
    fprintf('\n %f + %fi, \t%f + %fi, \t\t%f + %fi ',real(Z1),imag(Z1),real(Z2),imag(Z2),real(Z3),imag(Z3))
end
fprintf('\n')
