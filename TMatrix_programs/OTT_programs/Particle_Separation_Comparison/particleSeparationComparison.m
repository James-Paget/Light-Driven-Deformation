%{
%}
disp("Program start");

close all;
ott.warning('once');
ott.change_warnings('off');

n_medium = 1.33;        % Water
n_particle = 1.59;      % Polystyrene
wavelength0 = 1064e-9;  % Vacuum wavelength
wavelength_medium = wavelength0 / n_medium;
%% radius=0.1* wavelength of medium is TOO SMALL for good scattering
radius = 0.5*wavelength_medium;
view_range = 2e-6;

iterations = 4;
factor = 0.01;
sample_number = 50;

beam_inc = ott.BscPlane(0, 0, 'polarisation', [ 1 0 ], ...
     'index_medium', n_medium, 'wavelength0', wavelength0);

for iter = 0:iterations-1
    shape_1 = ott.shapes.Sphere(radius, [-radius/2.0-iter*factor*radius;0;0]);
    shape_2 = ott.shapes.Sphere(radius, [ radius/2.0+iter*factor*radius;0;0]);
    shape_union = ott.shapes.Union([shape_1, shape_2]);
    
    T_union = ott.TmatrixMie.simple(shape_union, 'wavelength0', wavelength0, ...
       'index_medium', n_medium, 'index_particle', n_particle);
    beam_scat_union = T_union * beam_inc;
    
    xrange = linspace(-1, 1, sample_number)*view_range;
    yrange = linspace(-1, 1, sample_number)*view_range;
    [xx, yy] = meshgrid(xrange, yrange);
    xyz = [xx(:) yy(:) zeros(size(xx(:)))].';
    [E,H] = emFieldXyz(beam_scat_union, xyz);
    E_mag=reshape(sqrt(sum(real(E).^2,1)),[sample_number, sample_number]);

    %E_mag = log(E_mag);

    figure();
    surf(xx, yy, E_mag);
    view(0, 90);
    title('Scattered E field at sep. iter= '+iter);
end

disp("Program end");