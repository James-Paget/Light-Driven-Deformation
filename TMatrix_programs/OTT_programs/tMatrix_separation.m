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

view_range = 2e-6;

radius = 0.5*wavelength_medium;
iterations = 1000;
%force_samples = 100;
factor = 0.01;  %Fractions of a particle radius to jump up by

%Generate beam, centered at origin, moving in +Z direction
beam = ott.BscPlane(0, 0, ...
    'polarisation', [ 1 0 ], ...
    'index_medium', n_medium, ...
    'wavelength0', wavelength0 ...
);

force_values = zeros(2,iterations);

for iter = 1:iterations
    %Setup no particle positions
    p1_pos = [-radius-(iter-1)*factor*radius; 0; 0];
    p2_pos = [ radius+(iter-1)*factor*radius; 0; 0];
    p1_shape = ott.shapes.Sphere(radius, p1_pos);
    p2_shape = ott.shapes.Sphere(radius, p2_pos);
    p_union = ott.shapes.Union([p1_shape, p2_shape]);
    %Generate T-Matrix for setup
    T = ott.TmatrixMie.simple(p_union, ...
       'wavelength0', wavelength0, ...
       'index_medium', n_medium, ...
       'index_particle', n_particle ...
    );
    %Get scattered and total field
    scattered_beam = T*beam;
    total_beam     = scattered_beam.totalField(beam);
    %Find force along X axis
    %xrange = [1;0;0] .* linspace(p1_pos(1), p2_pos(1), force_samples);
    %fxyz = ott.forcetorque(beam, T, xrange);
    %% USING VECTOR NORM INSTEAD OF X DIRECTION AS NET MOMENTUM CANCELLED
    [force, torque] = ott.forcetorque(beam, scattered_beam);
    force_values(1, iter) = (iter-1)*factor;
    force_values(2, iter) = vecnorm(force);%force(1);
end

figure();
plot(force_values(1,:), force_values(2,:));
xlabel("Separation in radi ("+radius+"m)");
ylabel("Relative joint force magnitude");


%% OSCILLATION SEEN, UNUSUAL MAX AND MIN SEEN FOR LARGE SEPARATION
disp("Program end");