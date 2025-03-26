%{
%}
disp("Program start");

close all;
ott.warning('once');
ott.change_warnings('off');

%%      NEED TO FIRST FIND THE OPTICAL TRAPPING PARAMETERS SO YOU WOULD EXPECT IT TO CICLE AROUND
%% DO PLOTS OF WHERE PARTICLES ARE LOCATED TOO --> FIND SAME THING FOR OTHER PARTICLES AROUND RING TO PROVE IT IS ACTUALLY WORKING
%%      FIND WHY FORCE IS SO LOW --> COULD BE TO DO WITH TENSO LOW VALUE PROBLEMS --> TRY VARY SAMPLES IT USES
%%          GET TRANSLATION MATRIX WORKING FOR TRUE CONFIRMATION OF SOME OF THESE IDEAS
%%

displayBeams = true;
Nmax = 4;

n_medium = 1.3;
n_particle = 1.6;
iterations = 1;        %Max number of particles in the ring
factor = 0.05;  %Fractions of a particle radius to jump up by

wavelength0_base = 1064e-9;  % Vacuum wavelength
wavelength_medium_base = wavelength0_base / n_medium;
wavelength0 = wavelength0_base;

%Outer part of ring = 0.75*wavelength0
%Middle part of ring = 0.65*wavelength0
%Inner part of ring = 0.55*wavelength0

scale_length = 0.4*wavelength0;
radius = 0.9*(pi*scale_length)/(1+iterations);%0.25*wavelength0;
forceType = "tensor";%forcetorque

force_figure = figure;

%Generate beam, centered at origin, moving in +Z direction
beam = ott.BscPmGauss('lg', [ 0 6 ], ...
    'index_medium', n_medium, ...
    'wavelength0', wavelength0, ...
    'NA', 0.8, ...
    'polarisation', [ 1 1i ], ...
    'Nmax', Nmax ...
);
%{
beam = ott.BscPlane(0, 0, ...
    'polarisation', [ 1 0 ], ...
    'index_medium', n_medium, ...
    'wavelength0', wavelength0 ...
);
%}

for focus_value = 0:0
    focus_iteration = focus_value;    %Which particle to get forces for (starting at particle 0)

    hold on
    for iter_particle = focus_iteration:iterations
        disp("iter_particle: Calc Force= "+iter_particle);
        force_values = zeros(2,iterations);
    
        max_particles = 1+iter_particle;
        particles = [];
        for index_particle = 0:max_particles-1
            particle_pos = FetchParticleOrigin("circle", scale_length, max_particles, index_particle, 0.0);
            particles{1+index_particle} = ott.shapes.Sphere(radius, particle_pos);
        end
    
        force = [0;0;0];
        if(forceType == "tensor")
            %Store matrix for union
            %% USE DDA FORMULTION FOR ACCURATE RESULTS %% --> TO ESNURE CROSS-SCATTERING
            shape_union = UnionAllParticles(particles);
            T_union = ott.TmatrixMie.simple(shape_union, ...
               'wavelength0', wavelength0, ...
               'index_medium', n_medium, ...
               'index_particle', n_particle ...
            );
    
            %Find scattered field
            beam_total = totalField(T_union*beam, beam);
        
            %Find force on specific particle in radial direction
            force_samples = 32;%53; %Test for varying valeus, sometimes the values land exactly on the central intensity peak and break the value
            force = ForceCalc_modified(particles{1+focus_iteration}.position, particles{1+focus_iteration}.radius, beam_total, force_samples);
    
        elseif(forceType == "forcetorque")
            %Store matrix for WHOLE system and system-particle of interest (focus_iteration)
            %% USE DDA FORMULTION FOR ACCURATE RESULTS %%
            shape_unionFull = UnionAllParticles(particles);
            T_unionFull = ott.TmatrixMie.simple(shape_unionFull, ...
               'wavelength0', wavelength0, ...
               'index_medium', n_medium, ...
               'index_particle', n_particle ...
            );
            %particlesReduced = particles(2:size(particles,2));   %Leaves the 0th particle --> particle of interest
            particlesReduced = [];
            for index_particle = 0:max_particles-1
                if(index_particle ~= focus_iteration)   %Leaves the 0th particle --> particle of interest
                    particle_pos = FetchParticleOrigin("circle", scale_length, max_particles, index_particle, 0.0);
                    particlesReduced{1+index_particle} = ott.shapes.Sphere(radius, particle_pos);
                end
            end
            shape_unionReduced = UnionAllParticles(particlesReduced);
            T_unionReduced = ott.TmatrixMie.simple(shape_unionReduced, ...
               'wavelength0', wavelength0, ...
               'index_medium', n_medium, ...
               'index_particle', n_particle ...
            );
    
            %Get force on specific particle from difference ith other union
            beam_scattered_full    = (T_unionFull*beam);
            beam_scattered_reduced = (T_unionReduced*beam);
            beam_scattered_difference = beam_scattered_reduced - beam_scattered_reduced;
            [force, torque] = ott.forcetorque(beam, beam_scattered_difference);
        end
    
        %Plot force in radial direction
        figure(force_figure);
        pos_mag = vecnorm(particles{1+focus_iteration}.position);
        theta_hat = [-particles{1+focus_iteration}.position(2); particles{1+focus_iteration}.position(1); 0]*(1.0/pos_mag);
        force_theta = theta_hat(1)*force(1) + theta_hat(2)*force(2);     %Component of force in theta-hat direction (CCW around ring) => F.ThetaHat
        disp("F      = "+force);
        disp("F_theta= "+force_theta);
        scatter(max_particles, force_theta);
        %ylim([-1e-15 1e-15]);
        title("Spheres in Leguerre Gaussian Order 2 Beam");
        xlabel("Particle Number");
        ylabel("Resultant force in theta-hat dir on particle "+focus_iteration);
    end
    hold off;
    
    if(displayBeams)
        figure();
        beam.visualise('axis', 'x', ...
            'mask', @(xyz) CheckParticlesMaskCondition(xyz, particles), ...
            'range', [-1,1]*2.0*scale_length );
        %Do plots for particle positions now
        for iter_particle = focus_iteration:iterations
            disp("iter_particle: Display Beam= "+iter_particle);
            force_values = zeros(2,iterations);
        
            max_particles = 1+iter_particle;
            particles = [];
            for index_particle = 0:max_particles-1
                particle_pos = FetchParticleOrigin("circle", scale_length, max_particles, index_particle, 0.0);
                particles{1+index_particle} = ott.shapes.Sphere(radius, particle_pos);
            end
    
            %Draw plot for each
            figure();
            beam.visualise('axis', 'z', ...
                'mask', @(xyz) CheckParticlesMaskCondition(xyz, particles), ...
                'range', [-1,1]*2.0*scale_length );
        end
    end
end

disp("Program end");