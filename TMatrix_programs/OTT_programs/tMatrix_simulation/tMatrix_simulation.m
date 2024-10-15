%{
. This program will run a T-Matrix simulation with given parameters,
calculate forces involved and then timestep the particles.
. The particles when then be placed back into the simulation for the
process to be repeated again.
%}

disp("Program start");

%Setup for a fresh run
close all;

n_medium = 1.3;
n_particle = 1.6;
nrel = n_particle/n_medium;
wavelength0 = 1.5*1064e-9;                      %Wavelength in free space
wavelength_medium = wavelength0/n_medium;   %Wavelength in medium
particle_radius = 2e-7;
scale_length = 0.8e-6;
regionWidth = 8e-6;         %The width to test force along (draw plot between)

beam_type = 'leguerre3_counterProp';
numerical_aperture = 1.02;  %Some quality of the beam

particle_number = 4;

%Generate a set of particles to work with
particles = cell(1, particle_number);           %Particle shapes stored here
particle_velocities = cell(1, particle_number); %Shapes do not store a velocity => stored externally like this
for index=1:particle_number
    %Specify shape
    particle_shape = ott.shapes.Shape.simple('sphere',particle_radius);
    particle_shape.position = FetchParticleOrigin("circle", scale_length, particle_number, index, 0);
    %Add shape to list
    particles{index} = particle_shape;
    particle_velocities{index} = [0;0;0];
end

iterations = 20;    %How many cycles of the dynamics of the motion you want to compute
timestep = 1.6e-3;    %How much time to simulate for each cycle step
frames(iterations) = struct('cdata',[],'colormap',[]);

figure_scatter  = figure;
figure_scatter2d= figure;
figure_beam     = figure;
figure_surface  = figure;
figure_movie    = figure;

for iter=1:iterations
    %**Complete a cycle of the dynamics calculation
    %Save data for animation at the end if required
    %pass
    %Generate a T-matrix for all particles
    particle_union = UnionAllParticles(particles);
    T = ott.TmatrixMie.simple(particle_union, 'index_relative', nrel, ...
    'index_medium', n_medium, 'wavelength0', wavelength0);
    %Create the beam present (fetch beam shape coefficient)
    switch beam_type
        case 'gaussian'
            %Circular polarisation Gaussian
            beam = ott.BscPmGauss('NA', NA, 'polarisation', [ 1 1i ], ...
            'index_medium', n_medium, 'wavelength0', wavelength0);
        case 'leguerre3'
            %Circular polarisation Leguerre Gaussian 3
            beam = ott.BscPmGauss('lg', [ 0 3 ], ...
            'polarisation', [ 1 1i ], 'NA', numerical_aperture, ...
            'index_medium', n_medium, 'wavelength0', wavelength0);
        case 'leguerre3_counterProp'
            %Circular polarisation Leguerre Gaussian 3
            beam_upper = ott.BscPmGauss('lg', [ 0 3 ], ...
            'polarisation', [ 1 1i ], 'NA', numerical_aperture, ...
            'index_medium', n_medium, 'wavelength0', wavelength0);

            beam_lower = ott.BscPmGauss('lg', [ 0 3 ], ...
            'polarisation', [ 1 -1i ], 'NA', numerical_aperture, ...
            'index_medium', n_medium, 'wavelength0', wavelength0);
            beam_lower = beam_lower.rotateXyz(pi, 0, 0);

            beam = beam_upper + beam_lower;
        otherwise
            %Circular polarisation Plane-Wave
            disp("Unrecognised beam, defaulting to plane-wave...");
            beam = ott.BscPmGauss('NA', numerical_aperture, 'polarisation', [ 1 1i ], ...
            'index_medium', n_medium, 'wavelength0', wavelength0);
    end
    %beam.power = 1.0;   %NOT recommended anymoreaccording to doc

    %Step all particles (rigidly) according to the forces they are experiencing
    for index = 1:particle_number
        force_net = FetchParticleNetForce(particles{index}, beam, T, 10);
        [stepped_position, stepped_velocity] = StepByForce_Euler(particles{index}, particle_velocities{index}, force_net, timestep); %## CAN DO A MORE ACCURATE RUNGE-KUTTA STEP TOO ##
        particles{index}.position  = stepped_position;
        particle_velocities{index} = stepped_velocity;
    end

    figure(figure_beam);
    beam.visualise('axis', 'z', ...
        'mask', @(xyz) CheckParticlesMaskCondition(xyz, particles), ...
        'range', [1,1]*regionWidth/2.0 );
    title("Incident Beam, iteration="+iter+".png");
    %frames(iter) = getframe(gcf);
    name = "frame_"+iter;
    
    figure(figure_scatter);
    pos_X_list = zeros(particle_number);
    pos_Y_list = zeros(particle_number);
    pos_Z_list = zeros(particle_number);
    for i = 1:particle_number
        pos_X_list(i) = particles{i}.position(1);
        pos_Y_list(i) = particles{i}.position(2);
        pos_Z_list(i) = particles{i}.position(3);
    end
    scatter3(pos_X_list, pos_Y_list, pos_Z_list);
    xlim([-regionWidth/2.0, regionWidth/2.0]);
    ylim([-regionWidth/2.0, regionWidth/2.0]);
    zlim([-regionWidth/2.0, regionWidth/2.0]);
    xlabel("X");
    ylabel("Y");
    zlabel("Z");
    frames(iter) = getframe(gcf);
    %saveas(gcf, name);
    %}

    figure(figure_scatter2d);
    pos_X_list = zeros(particle_number);
    pos_Y_list = zeros(particle_number);
    for i = 1:particle_number
        pos_X_list(i) = particles{i}.position(1);
        pos_Y_list(i) = particles{i}.position(2);
    end
    scatter(pos_X_list, pos_Y_list);
    xlim([-regionWidth/2.0, regionWidth/2.0]);
    ylim([-regionWidth/2.0, regionWidth/2.0]);
    xlabel("X");
    ylabel("Y");
    frames(iter) = getframe(gcf);
    %saveas(gcf, name);

    %Plotting results for given iterations
    %%{
    if(iter == 1)
        %3D plot of force norms
        %##### TRY GET MASK WORKING HERE TOO MAYBE ######
        force_samples = 20;         %How many sections tested on each axis
        force_sampleWidth = regionWidth/force_samples;
        force_X_set = [1;0;0] .* linspace(-regionWidth/2.0, regionWidth/2.0, force_samples);
        force_Y_set = [0;1;0] .* linspace(-regionWidth/2.0, regionWidth/2.0, force_samples);
        force_Z_plane = 0;          %Z=<?> plane that will be tested
        force_norms = [];  %########### FIX 2D SIZE HERE ######## -> NxN
        for i = 1: force_samples
            for j = 1: force_samples
                fXyz = ott.forcetorque(beam, T, 'position', [force_X_set(1, i); force_Y_set(2, j); force_Z_plane]);
                force_norms(i, j) = vecnorm(fXyz);
            end
        end

        [force_X_mesh, force_Y_mesh] = meshgrid(force_X_set(1,:),force_Y_set(2,:));
        figure(figure_surface);
        surf(force_X_mesh, force_Y_mesh, force_norms);
        xlabel('X-axis');
        ylabel('Y-axis');
        zlabel('Z-axis');
        title('Force Norm');
    end
    %%}
end

disp("Program end, showing manimation...");

figure(figure_movie);
movie(frames, 1000);