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

iterations = 10;
factor = 0.4;
sample_number = 50;

beam_inc = ott.BscPlane(0, 0, 'polarisation', [ 1 0 ], ...
     'index_medium', n_medium, 'wavelength0', wavelength0);

figure_E_x_plots = figure;
figure_F_x_plots = figure;
hold on;
for iter = 0:iterations-1
    view_width  = radius/2.0+iter*factor*radius; %Both sides of the axis => view_range = one half of the width covered
    view_height = radius;
    shape_1 = ott.shapes.Sphere(radius, [-radius/2.0-iter*factor*radius;0;0]);
    shape_2 = ott.shapes.Sphere(radius, [ radius/2.0+iter*factor*radius;0;0]);
    shape_union = ott.shapes.Union([shape_1, shape_2]);

    %% Unhash to see interaction of beam on a sphere at the origin, and 2
    %% unioned at origin --> should prove that translations are being applied
    %shape_test_1 = ott.shapes.Sphere(radius, [0;0;0]);
    %shape_test_2 = ott.shapes.Sphere(radius, [0;0;0]);
    %shape_union = ott.shapes.Union([shape_test_1, shape_test_2]);
    %shape_union = shape_test_1
    
    T_union = ott.TmatrixMie.simple(shape_union, 'wavelength0', wavelength0, ...
       'index_medium', n_medium, 'index_particle', n_particle);
    beam_scat_union = T_union * beam_inc;

    %%##
    %%## HERE CAN CONSIDER JUST SCATTERED BEAM, OR THE TOTAL BEAM -> WE WANT OT JUST CONSIDER SCATTERED (IDEALLY JUST INTERACTION BETWEEN THE 2 PARTICLES)
    %%##
    
    %Get E field for this scattering
    xrange = linspace(-1, 1, sample_number)*view_width;
    yrange = linspace(-1, 1, sample_number)*view_height;
    [xx, yy] = meshgrid(xrange, yrange);
    xyz = [xx(:) yy(:) zeros(size(xx(:)))].';
    [E_scat_uniform, H_scat_uniform] = beam_scat_union.emFieldXyz(xyz);
    [E_inc, H_inc] = beam_inc.emFieldXyz(xyz);
    %E_mag=reshape(sqrt(sum(real(E).^2,1)),[sample_number, sample_number]);
    %--> Total magnitude of E at each grid space

    %% Get from this the grid of E field in the X direction -> proportional
    %% to force in X direction from lorentz force eq.
    %% Also, we just want to see the components for field scattered from one
    %% sphere to the other
    E_scat_uniform_x = real(E_scat_uniform(1,:));
    E_scat_uniform_grid_x = reshape(E_scat_uniform_x,[sample_number, sample_number]);

    %Consider the log of the grid instead, as is a bright spot at center
    %%##
    %%## MAY BREAK FOR NEGATIVES
    %%##
    E_scat_uniform_grid_x = log(E_scat_uniform_grid_x);

    %Average values between radius in +- y axis for net E in x direction
    %Note; view_range always set between the 2 particle centres =>
    E_scat_uniform_x_averaged_list = zeros(sample_number); %Holds averaged E_x at each sample pos in x axis between particles
    for i = 1:sample_number
        %This will go from p1_centre to p2_centre, due to how view_range is
        %fixed
        E_scat_uniform_x_averaged = 0.0;
        for j = 1:sample_number
            %This will cover from bottom of p1 & 2 to top of p1 & 2
            E_scat_uniform_x_averaged = E_scat_uniform_x_averaged +E_scat_uniform_grid_x(j,i);
        end
        E_scat_uniform_x_averaged_list(i) = E_scat_uniform_x_averaged;
    end

    %Get force acting on each particle from JUST the SCATTERED field
    force_samples = 53; %Test for varying valeus, sometimes the values land exactly on the central intensity peak and break the value
    Ftot_shape_1 = ForceCalc_modified(shape_1.position, shape_1.radius, beam_scat_union, force_samples);
    %Ftot_shape_2 = ForceCalc_modified(shape_2.position, shape_2.radius, beam_scat_union, force_samples);

    %%##
    %%## VARY NMAX, SEE HOW IT AFFECTS RESULTS
    %%##

    %Plot E_x values across X-axis for each separation
    figure(figure_E_x_plots);
    plotLabel = "iter; "+iter;
    plot(xrange, E_scat_uniform_x_averaged_list);
    title("E_x averaged between 2 spheres of varying separation");
    xlabel("Separation (m)");
    ylabel("Averaged E_x from scattered field");
    %legend(plotLabel);

    %Plot log(F_x) values across X-axis for each separation
    figure(figure_F_x_plots);
    x_separation = shape_2.position(1) -shape_1.position(1);
    scatter(x_separation, log(Ftot_shape_1(1)));
    %scatter(x_separation, log(Ftot_shape_2(1)));
    title("log(F_x) on left particle for varying separation of 2 particles");
    xlabel("Separation (m)");
    ylabel("F_x from scattered field");

    %Unhash to view surf plot every N iterations
    %{
    if(mod(iter, 2) == 0)
        figure();
        surf(xx, yy, abs(E_scat_uniform_grid_x));
        view(0, 90);
        title("Beam E Field, Scattered Only, iter= "+iter);
    end
    %}
end
hold off;

disp("Program end");