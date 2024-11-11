%{
%}
disp("Program start");

close all;
ott.warning('once');
ott.change_warnings('off');

n_medium = 1.33;        % Water
n_particle = 1.59;      % Polystyrene

view_range = 2e-6;

Nmax = 50;
iterations = 1000;
iterations_wavelength = 1;
%force_samples = 100;
factor = 0.05;  %Fractions of a particle radius to jump up by

wavelength0_base = 1064e-9;  % Vacuum wavelength
wavelength_medium_base = wavelength0_base / n_medium;
radius = 0.5*wavelength_medium_base;

hold on;
for iter_wave = 1:iterations_wavelength
    force_values = zeros(2,iterations);
    wavelength0 = wavelength0_base*iter_wave;

    %Generate beam, centered at origin, moving in +Z direction
    beam = ott.BscPlane(0, 0, ...
        'polarisation', [ 1 0 ], ...
        'index_medium', n_medium, ...
        'wavelength0', wavelength0, ...
        'Nmax', Nmax ...
    );
    for iter = 1:iterations
        %Setup no particle positions
        p1_pos = [-radius-(iter)*factor*radius; 0; 0];
        p2_pos = [ radius+(iter)*factor*radius; 0; 0];
        p1_shape = ott.shapes.Sphere(radius, p1_pos);
        p2_shape = ott.shapes.Sphere(radius, p2_pos);
        p_union = ott.shapes.Union([p1_shape, p2_shape]);
        %Generate T-Matrix for setup
        p1_T = ott.TmatrixMie.simple(p1_shape, ...
           'wavelength0', wavelength0, ...
           'index_medium', n_medium, ...
           'index_particle', n_particle ...
        );
        p2_T = ott.TmatrixMie.simple(p2_shape, ...
           'wavelength0', wavelength0, ...
           'index_medium', n_medium, ...
           'index_particle', n_particle ...
        );
        T = ott.TmatrixMie.simple(p_union, ...
           'wavelength0', wavelength0, ...
           'index_medium', n_medium, ...
           'index_particle', n_particle ...
        );
        %Get scattered and total field
        scattered_p1_beam = p1_T*translateXyz(beam, -p1_pos);
        scattered_p2_beam = p2_T*translateXyz(beam, -p2_pos);
        scattered_union_beam = T*beam;
        %Find force along X axis
        [p1_force, p1_torque] = ott.forcetorque(beam, scattered_p1_beam);
        [p2_force, p2_torque] = ott.forcetorque(beam, scattered_p2_beam);
        [union_force, union_torque] = ott.forcetorque(beam, scattered_union_beam);
        interparticle_force = union_torque -p1_force -p2_force;
        force_values(1, iter) = 2*(iter)*factor*radius;
        force_values(2, iter) = interparticle_force(1);%vecnorm(force);%force(1);
    
        %{
        if(interparticle_force(1) ~= 0)
            %View forces at some moment in time to observe what is occurring
            disp("Viewing interparticle force");
            disp(interparticle_force);
        end
        %}
    end
    figure();
    plot(force_values(1,:), force_values(2,:));
    ylim([-1e-15 1e-15]);
    xlabel("Separation (m)");
    ylabel("Resultant X force from inter-particle scattering");
end
hold off;

disp("Program end");