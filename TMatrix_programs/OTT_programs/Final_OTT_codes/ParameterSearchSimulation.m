%{
This program has several options in a process you can follow in order to
help you find parameters suitable for trapping.

.(1) "View_Beam" = view different orders of laguerre beam, see where fringes
occur, swept through beam orders
.(2) "Sweep_Particle_Dist" = Fixes all other values, sweeps pos only, gives
force on 0th particle in t=his configuration, manually check for trapping
.(3) "Sweep_Particle_Rad" = Fixes all other values, sweeps rad only, gives
force on 0th particle in this configuration, manually check for trapping
.(4) "Sweep_Particle_Dist_Rad" = Fixes all other values, sweeps pos AND rad only, gives
force on 0th particle in this configuration, manually check for trapping
.(5) "All_Particles_R_Theta" = Plots the r and theta forces on each
particle in a circle formation in a Laguerre beam.
%}
disp("Program start");


%Startup
close all;
ott.warning('once');
ott.change_warnings('off');

%Queue up several tests in one run
mode_sequence = ["All_Particles_R_Theta"];
%%MODES
for mode = mode_sequence
    if(mode == "View_Beam")
        disp("Performing View_Beam");
        view_range = 2.0e-6;
        View_Beam(view_range, 3, 6);
    
    elseif(mode == "Sweep_Particle_Dist")
        disp("Performing Sweep_Particle_Dist");
        view_range      = 2.0e-6;
        laguerreOrder   = 3;
        target_dist     = 7.0e-7;
        Sweep_Particle_Dist(view_range, laguerreOrder, target_dist, 0.2, 10);
    
    elseif(mode == "Sweep_Particle_Rad")
        disp("Performing Sweep_Particle_Rad");
        view_range      = 2.0e-6;
        laguerreOrder   = 3;
        distance        = 7.0e-7;
        target_rad      = 2.4e-7;
        Sweep_Particle_Rad(view_range, laguerreOrder, distance, target_rad, 0.2, 16);
    
    elseif(mode == "Sweep_Particle_Dist_Rad")
        disp("Performing Sweep_Particle_Dist_Rad");
        laguerreOrder   = 3;
        target_dist     = 7.0e-7;
        target_rad      = 2.4e-7;
        Sweep_Particle_Dist_Rad(laguerreOrder, target_dist, 0.5, 10, target_rad, 0.5, 10, "X");
        Sweep_Particle_Dist_Rad(laguerreOrder, target_dist, 0.5, 10, target_rad, 0.5, 10, "Y");
    
    elseif(mode == "All_Particles_R_Theta")
        disp("Performing All_Particles_R_Theta");
        particle_rad = 0.5e-6;
        particle_dist = 1e-6;
        num_particles = 5;
        min_gap_factor = 0;
        All_Particles_R_Theta(num_particles, particle_dist, particle_rad, min_gap_factor);
    end
end


function output = View_Beam(view_range, laguerreOrder_target_lower, laguerreOrder_target_upper)
    Nmax = 70;

    n_medium = 1.3;
    n_particle = 1.6;
    wavelength0_base = 1064e-9;         %Vacuum wavelength
    wavelength_medium_base = wavelength0_base / n_medium;
    wavelength0 = wavelength0_base;
    beam_polarisation = [ 1 1i ];
    NA = 1.4;

    beam_order_comp_figure = figure;
    subplot_width  = 1+ laguerreOrder_target_upper-laguerreOrder_target_lower;
    subplot_height = 2;

    figure(beam_order_comp_figure);
    for laguerreOrder_iter = laguerreOrder_target_lower:laguerreOrder_target_upper
        laguerreOrder = laguerreOrder_iter;
        subplot_index = 2.0*(laguerreOrder_iter-laguerreOrder_target_lower)+1;

        beam = Get_Beam(laguerreOrder, n_medium, wavelength0, NA, beam_polarisation, Nmax);
        
        disp("---");
        disp("W; "+ subplot_width);
        disp("H; "+ subplot_height);
        disp("I: "+ subplot_index);
        subplot(subplot_width, subplot_height, subplot_index);
        beam.visualise('axis', 'z', ...
            'mask', @(xyz) vecnorm(xyz)<0.05*view_range, ...
            'range', [-1,1]*view_range );
        title("Laguerre beam order "+laguerreOrder_iter);
        subplot(subplot_width, subplot_height, subplot_index+1);
        beam.visualise('axis', 'x', ...
            'mask', @(xyz) vecnorm(xyz)<0.05*view_range, ...
            'range', [-1,1]*view_range );
        title("Laguerre beam order "+laguerreOrder_iter);
    end
    output=1;
end

function output = All_Particles_R_Theta(particle_iter_max, distance, radius, min_gap_factor)
    % Check if particles will overlap and make radius smaller if needed.
    r_max = distance * sin(pi/particle_iter_max) * (1-min_gap_factor);
    if radius > r_max && particle_iter_max > 1
        radius = r_max;
        disp("All_Particles_R_Theta: changed radius to "+radius);
    end

    displayBeams = true;
    Nmax = 4;
    laguerreOrder = 3;
    forceType = "tensor";

    n_medium = 1.3;
    n_particle = 1.6;          
    wavelength0_base = 1064e-9;         %Vacuum wavelength
    wavelength_medium_base = wavelength0_base / n_medium;
    wavelength0 = wavelength0_base;
    beam_polarisation = [ 1 1i ];
    NA = 1.4;
    force_samples = 10;%32;
    R_Theta_Figure = figure;

    focusses = linspace(1, particle_iter_max, particle_iter_max);
    rs = linspace(0, 0, particle_iter_max);
    thetas = linspace(0, 0, particle_iter_max);
    beam = Get_Beam(laguerreOrder, n_medium, wavelength0, NA, beam_polarisation, Nmax);
    particles = Get_Particles("circle", particle_iter_max, distance, radius, 0.0); %Get a set of particles

    for focus_particle_iter = 1:particle_iter_max %Which particle index to observe in measurements (1->N)
        disp("Focussing particle " + focus_particle_iter);                   
        force = Perform_MST_Calc(particles, focus_particle_iter, beam, wavelength0, n_medium, n_particle);   %Calculate forces on ith particle in set
        force_rhat = Get_Specific_ForceScalar("r_hat", force, particles{focus_particle_iter});
        force_scalar = Get_Specific_ForceScalar("theta_hat", force, particles{focus_particle_iter});
        
        focusses(focus_particle_iter) = focus_particle_iter;
        rs(focus_particle_iter) = force_rhat;
        thetas(focus_particle_iter) = force_scalar;
    end
    
    % make graph line (scatter for one particle). (change to scatter only?)
    hold on;
    if particle_iter_max == 1
        scatter(focusses, rs)
        scatter(focusses, thetas)
    else
        plot(focusses, rs)
        plot(focusses, thetas)
    end
    legend('r components','theta components')
    xlabel("Focus particle")
    ylabel("Polar components")
end

function output = Sweep_Particle_Dist(view_range, laguerreOrder, target_dist, variation_dist, dist_iter_max)
    Nmax = 4;

    n_medium = 1.3;
    n_particle = 1.6;
    wavelength0_base = 1064e-9;
    wavelength_medium_base = wavelength0_base / n_medium;
    wavelength0 = wavelength0_base;
    beam_polarisation = [ 1 1i ];
    NA = 1.4;

    beam = Get_Beam(laguerreOrder, n_medium, wavelength0, NA, beam_polarisation, Nmax);

    dist_sweep_figure = figure;
    figure(dist_sweep_figure);
    subplot_width  = ceil(sqrt(dist_iter_max));
    subplot_height = ceil(sqrt(dist_iter_max));
    for dist_iter = 1:dist_iter_max
        distance = target_dist +( dist_iter-(dist_iter_max/2.0) )*(variation_dist*target_dist);
        particle_pos    = [distance; 0; 0];
        particle_radius = 0.05*view_range;

        subplot(subplot_width, subplot_height, dist_iter);
        beam.visualise('axis', 'z', ...
            'mask', @(xyz) vecnorm(xyz-particle_pos)<particle_radius, ...
            'range', [-1,1]*view_range );
        title("Distance= "+distance);
    end
    output=1;
end

function output = Sweep_Particle_Rad(view_range, laguerreOrder, distance, target_rad, variation_rad, rad_iter_max)
    Nmax = 4;

    n_medium = 1.3;
    n_particle = 1.6;
    wavelength0_base = 1064e-9;
    wavelength_medium_base = wavelength0_base / n_medium;
    wavelength0 = wavelength0_base;
    beam_polarisation = [ 1 1i ];
    NA = 1.4;

    beam = Get_Beam(laguerreOrder, n_medium, wavelength0, NA, beam_polarisation, Nmax);

    rad_sweep_figure = figure;
    figure(rad_sweep_figure);
    subplot_width  = ceil(sqrt(rad_iter_max));
    subplot_height = ceil(sqrt(rad_iter_max));
    for rad_iter = 1:rad_iter_max
        radius = target_rad +( rad_iter-(rad_iter_max/2.0) )*(variation_rad*target_rad);
        particle_pos    = [distance; 0; 0];
        particle_radius = radius;

        subplot(subplot_width, subplot_height, rad_iter);
        beam.visualise('axis', 'z', ...
            'mask', @(xyz) vecnorm(xyz-particle_pos)<particle_radius, ...
            'range', [-1,1]*view_range );
        title("Radius= "+radius+"   (With Distance= "+distance+")");
    end
    output=1;
end

function output = Sweep_Particle_Dist_Rad(laguerreOrder, target_dist, variation_dist, dist_iter_max, target_rad, variation_rad, rad_iter_max, force_format)
    Nmax = 4;

    n_medium = 1.3;
    n_particle = 1.6;
    particle_iter_max = 1;              %Max number of particles in the ring (1->N)
    focus_particle_iter = 1;            %Which particle index to observe in measurements (1->N)
    wavelength0_base = 1064e-9;         %Vacuum wavelength
    wavelength_medium_base = wavelength0_base / n_medium;
    wavelength0 = wavelength0_base;
    beam_polarisation = [ 1 1i ];
    NA = 1.4;

    sweepDistRad_force_figure = figure;

    sweepData = {};     %Holds objects that are vectors [distance, radius, force on ith]
    beam = Get_Beam(laguerreOrder, n_medium, wavelength0, NA, beam_polarisation, Nmax);
    for dist_iter = 1:dist_iter_max
        disp("dist_iter= "+dist_iter+"/"+dist_iter_max);
        distance = target_dist +( dist_iter-(dist_iter_max/2.0) )*(variation_dist*target_dist);
        for rad_iter = 1:rad_iter_max
            disp("    rad_iter= "+rad_iter+"/"+rad_iter_max);
            radius = target_rad +( rad_iter-(rad_iter_max/2.0) )*(variation_rad*target_rad);

            particles = Get_Particles("circle", particle_iter_max, distance, radius, 0.0);                       %Get a set of particles
            force = Perform_MST_Calc(particles, focus_particle_iter, beam, wavelength0, n_medium, n_particle);   %Calculate forces on ith particle in set
            force_scalar = Get_Specific_ForceScalar(force_format, force, particles{focus_particle_iter});
            sweepData{rad_iter, dist_iter} = [distance, radius, force_scalar];
        end
    end
    Display_3D_Sweep(sweepDistRad_force_figure, target_dist, variation_dist, dist_iter_max, target_rad, variation_rad, rad_iter_max, "Force "+force_format+" for varying parameters", "distance (m)", "radius (m)", "Force "+force_format+" on "+focus_particle_iter, sweepData);     %Display force with DIST and RAD sweeping
    output=1;
end

function output = Display_3D_Sweep(spec_figure, target_dist, variation_dist, dist_iter_max, target_rad, variation_rad, rad_iter_max, title_labelling, xaxis_labelling, yaxis_labelling, zaxis_labelling, data)
    %{
    .spec_figure = figure to display this data on
    .*_labelling = axis labels
    .data = 3D grid of data, holds set of vectors with [x,y,z] values
    %}
    %assignin("base", "data", data);
    x = [1;0;0] .* linspace( ...
        target_dist -target_dist*variation_dist, ...
        target_dist +target_dist*variation_dist, ...
        dist_iter_max ...
    );
    y = [0;1;0] .* linspace( ...
        target_rad -target_rad*variation_rad, ...
        target_rad +target_rad*variation_rad, ...
        rad_iter_max ...
    );
    [X, Y] = meshgrid(x(1,:),y(2,:));

    Z = zeros( size(data,1), size(data,2) );
    for i = 1: size(Z,1);
        for j = 1: size(Z,2);
            Z(i, j) = data{i, j}(3);
        end
    end

    figure(spec_figure);
    surf(X,Y,Z);
    title(title_labelling);
    xlabel(xaxis_labelling);
    ylabel(yaxis_labelling);
    zlabel(zaxis_labelling);

    output=1;
end

function beam = Get_Beam(laguerreOrder, n_medium, wavelength0, NA, polarisation, Nmax)
    %{
    .laguereeOrder = order of beam

    . Add option to vary beams
    %}
    beam = ott.BscPmGauss('lg', [ 0 laguerreOrder ], ...
        'index_medium', n_medium, ...
        'wavelength0', wavelength0, ...
        'NA', NA, ...
        'polarisation', polarisation, ...
        'Nmax', Nmax ...
    );
end

function force_scalar = Get_Specific_ForceScalar(force_format, force, particle)
    %{
    .force_format = the specific type of component you want to be extracted
    from this net force
    .force = vector force that you want to pull a scalar from
    %}
    force_scalar = 0.0;
    if(force_format == "theta_hat")
        force_angle = atan2(force(2), force(1));
        pos_angle = atan2(particle.position(2), particle.position(1));
        force_scalar = sin(force_angle - pos_angle);
    elseif(force_format == "r_hat")
        force_angle = atan2(force(2), force(1));
        pos_angle = atan2(particle.position(2), particle.position(1));
        force_scalar = cos(force_angle - pos_angle);
    elseif(force_format == "X")
        force_scalar = force(1);
    elseif(force_format == "Y")
        force_scalar = force(2);
    else
        %Do norm if nothing else specified
        disp("Invalid format, force norm being displayed...");
        force_scalar = vecnorm(force);
    end
end

function particles = Get_Particles(layout_type, particle_number, distance, radius, Z_plane)
    %{
    .layout_type = how to align particles in the system (e.g in a circle or on a line)
    .particle_number = number of particles to put into formation
    .distance = scale distance for each formation (e.g for circles, is the distance from origin for each particle)
    .radius = radius of spherical particles
    .Z_plane = Z plane that all particle sit in
    
    (1) For each particle wanted
    (2) Find its parameters and add to particle set
    %}
    particles = {};
    for particle_iter = 1:particle_number
        if(layout_type == "line")
            %Spreads particle along the x-axis, from 0 to width
            width = distance;
            spacing = width/particle_number;
            pos = [(particle_iter-1)*spacing; 0; Z_plane];
            particles{particle_iter} = ott.shapes.Sphere(radius, pos);
        elseif(layout_type == "circle")
            %Spread evenly along the circumference of a circle in the Z=0 plane
            theta = 2*pi/particle_number;
            pos = [distance*cos(theta*(particle_iter-1));distance*sin(theta*(particle_iter-1));Z_plane];
            particles{particle_iter} = ott.shapes.Sphere(radius, pos);
        else
            %Default to coord origin if invalid layout_type given
            %pass
        end
    end
end

function force = Perform_MST_Calc(particles, focus_particle_index, ibeam, wavelength0, n_medium, n_particle)
    %{
    .Perform a calculation of the Maxwell Stress-Tensor

    .particles = set of particles in the ring when this calcualtion is
    performed
    .focus_particle_index = the index of the particle in particles to
    consider (from 1->N)
    .ibeam = incident beam
    .force samples = number of points to calculate integral over in MST calc

    (1) Combines all particles into 1, gte T matrix for this combination
    (2) Find scattered field from set using incident beam
    %}
    shape_union = UnionAllParticles(particles);
    T_union = ott.TmatrixMie.simple(shape_union, ...
       'wavelength0', wavelength0, ...
       'index_medium', n_medium, ...
       'index_particle', n_particle ...
    );

    %Find scattered field
    beam_total = totalField(T_union*ibeam, ibeam);

    %Find force on specific particle in radial direction
    %force = ForceCalc_modified(particles{focus_particle_index}.position,
    force = ForceCalc_quad(particles{focus_particle_index}.position, particles{focus_particle_index}.radius, beam_total);
end


disp("Program end");