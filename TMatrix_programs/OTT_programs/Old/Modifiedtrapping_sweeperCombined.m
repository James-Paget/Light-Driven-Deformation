%{
This program has several options in a process you can follow in order to
help you find parameters suitable for trapping.

%}
disp("Program start");


%Startup
close all;
ott.warning('once');
ott.change_warnings('off');

%Queue up several tests in one run
mode_sequence = ["Time_Simulation_Multiple"]; % Sweep_Particle_In_Planes_Forces, Time_Simulation, Time_Simulation_Multiple
%%MODES
for mode = mode_sequence
    if(mode == "View_Beam")
        disp("Performing View_Beam");
        view_range = 2.0e-6;
        minLg = 8;
        maxLg = 8;
        % Nmaxs = 10;
        Nmaxs = [10, 15, 20, 30];
        for Nmax = Nmaxs
            View_Beam(view_range, minLg, maxLg, Nmax);
        end
    
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
        laguerreOrder   = 8;
        target_dist     = 10.0e-7;
        target_rad      = 4e-7;
        Sweep_Particle_Dist_Rad(laguerreOrder, target_dist, 0.2, 1, target_rad, 0.2, 3, "theta_hat");
        % Sweep_Particle_Dist_Rad(laguerreOrder, target_dist, 0.5, 10, target_rad, 0.5, 10, "Y");
    
    elseif(mode == "All_Particles_R_Theta")
        disp("Performing All_Particles_R_Theta");
        particle_rad = 0.5e-6;
        particle_dist = 1e-6;
        num_particles = 3;
        min_gap_factor = 0;
        shouldPlot = true;
        All_Particles_R_Theta(num_particles, particle_dist, particle_rad, min_gap_factor, shouldPlot);

    elseif(mode == "Show_Particle_Forces")
        disp("Performing Show_Particle_Forces");
        laguerreOrder     = 8;
        target_dist       = 7.0e-7;
        target_rad        = 4e-7;
        particle_iter_max = 4;
        force_scaling     = 7e9;
        Show_Particle_Forces(laguerreOrder, target_dist, target_rad, particle_iter_max);

    elseif(mode == "Sweep_Particle_In_Planes_Forces")
        disp("Performing Sweep_Particle_In_Planes_Forces");
        laguerreOrder = 8;
        % Num angles
        ang_max = pi*2;
        ang_iter_max = 10;
        
        % Particle Radius
        target_rad = 0.6e-6;
        variation_rad = 0.3;
        rad_iter_max = 1;
        
        % Distance from the origin
        target_dist = 1.5e-6;
        variation_dist = 0.3;
        dist_iter_max = 3;
        
        % Z value
        z_min = 0.8e-6;
        z_max = 1.1e-6;
        z_iter_max = 1;
        
        Sweep_Particle_In_Planes_Forces(laguerreOrder, target_rad, variation_rad, rad_iter_max, target_dist, variation_dist, dist_iter_max, ang_iter_max, ang_max, z_min, z_max, z_iter_max)
    
    elseif(mode == "Time_Simulation")
        disp("Performing Time_Simulation");
        % d = 1.25e-6;
        % angle = pi/6;
        d = 0.1e-6;
        angle = pi/6;
        initial_position = [d*cos(angle);d*sin(angle);0]; % [d;0;0]
        initial_velocity = [0;0;0];

        dt = 0.7e-4;
        t_max = 2*dt; % 150e-4
        radius = 0.5e-6;
        Nmax = 15;
        laguerreOrder = 8;
        figTitle = "radius = "+radius+ ", distance = "+ distance;

        Time_Simulation(initial_position, initial_velocity, t_max, dt, radius, laguerreOrder, figTitle, Nmax)

    elseif(mode == "Time_Simulation_Multiple")
        disp("Performing Time_Simulation_Multiple");

        t_max = 5e-4;
        dt = 1e-4;
        laguerreOrder = 8;
        pos_angle = pi/6;
        initial_velocity = [0;0;0];
        Nmax = 10;

        distMin = 1.2e-6;
        distMax = 1.3e-6;
        distNum = 4;

        radMin = 0.4e-6;
        radMax = 0.6e-6;
        radNum = 3;
    
        dispStr = "";
        for dist_iter = 1:distNum
            distance = distMin + (dist_iter-1)/distNum * (distMax - distMin);
            initial_position = [distance * cos(pos_angle); distance * sin(pos_angle); 0];

            for rad_iter = 1:radNum
                radius = radMin + (rad_iter-1)/radNum * (radMax - radMin);
                figNum = (dist_iter-1) * radNum + rad_iter;
                
                titleStr = "Figure "+ figNum + " : radius = "+radius+ ", distance = "+ distance;
                disp(titleStr)
                dispStr = dispStr +"Figure "+ figNum + " : radius = "+radius+ ", distance = "+ distance + newline;
                Time_Simulation(initial_position, initial_velocity, t_max, dt, radius, laguerreOrder, figTitle, Nmax)
            end
        end
        disp(dispStr)



    end
end

function output = Time_Simulation(initial_position, initial_velocity, t_max, dt, radius, laguerreOrder, figTitle, Nmax)
    n_medium = 1.0;
    n_particle = 1.55;
    focus_particle_iter = 1;
    particle_iter_max = 1;
    wavelength0_base = 1064e-9;         %Vacuum wavelength
    wavelength_medium_base = wavelength0_base / n_medium;
    wavelength0 = wavelength0_base;
    beam_polarisation = [ 1 1i ];
    NA = 1.3;
    beam = Get_Beam(laguerreOrder, n_medium, wavelength0, NA, beam_polarisation, Nmax);
    
    particles = {};
    num_times = floor(t_max/dt);
    times = linspace(0, t_max, num_times);
    
    coords = zeros(num_times, 3);
    velocities = zeros(num_times, 3);
    forces = zeros(num_times, 3);
    position = initial_position;
    velocity = initial_velocity;

    for iter = 1:num_times
        t = times(iter);
        particles{1} = ott.shapes.Sphere(radius, position);
        coords(iter, :) = particles{focus_particle_iter}.position;
        velocities(iter, :) = velocity;


        force = Perform_MST_Calc(particles, focus_particle_iter, beam, wavelength0, n_medium, n_particle);   %Calculate forces on ith particle in set

        force_scalar_x = Get_Specific_ForceScalar("X", force, particles{focus_particle_iter});
        force_scalar_y = Get_Specific_ForceScalar("Y", force, particles{focus_particle_iter});
        force_scalar_z = Get_Specific_ForceScalar("Z", force, particles{focus_particle_iter});
        
        forcexyz = [force_scalar_x; force_scalar_y; force_scalar_z];

        forces(iter, :) = forcexyz;

        % Time step
        % calcite density = 2.7g/cm^3 -> 1 um^3 -> 2.7e-3 * 10^-12 =
        % 2.7e-15 kg
        mass = 2.7e-15;
        velocity = velocity + forces(iter, :)' * dt / mass;
        position = position + velocity * dt;

    end
    
    force_scaling_type = "norm"; % scale, norm, xy
    Plot_Force_Arrows(coords, forces, force_scaling_type, radius, figTitle)
    save('coordMatrix.mat', 'coords');
    save('velMatrix.mat', 'velocities');
end



function output = View_Beam(view_range, laguerreOrder_target_lower, laguerreOrder_target_upper, Nmax)
    % Nmax = 100;

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

function output = All_Particles_R_Theta(particle_iter_max, distance, radius, min_gap_factor, shouldPlot)
    % Check if particles will overlap and make radius smaller if needed.
    r_max = distance * sin(pi/particle_iter_max) * (1-min_gap_factor);
    if radius > r_max && particle_iter_max > 1
        radius = r_max;
        disp("All_Particles_R_Theta: changed radius to "+radius);
    end

    displayBeams = true;
    Nmax = 15;
    laguerreOrder = 8;
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

    if shouldPlot
        % make graph line (scatter for one particle).
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
    output = 1;
    % output = { focusses, rs, thetas};
end

function output = Plot_Force_Arrows(coords, forces, force_scaling_type, radius, figTitle)
    % makes quiver plot of the forces at the coords.
    % force_scaling_type = "scale" scales the force arrows to a good size.
    % force_scaling_type == "norm" normalises all the force arrows.
    % radius is the particle radius and is only used for titling the plot.
    % iter is used for naming the graph when saving.
    fig = figure();

    if length(coords(1,:)) == 2
        hold on
        drawArrow = @(pos,force) quiver( pos(1),pos(2),force(1),force(2) );
        scale = sum(abs(coords))/sum(abs(forces));
        particle_iter_max = length(coords(:,1));
        for i = 1:particle_iter_max
            if force_scaling_type == "scale"
                scaled_forces = forces(i,:) .* scale * 0.2;
            elseif force_scaling_type == "norm"
                scaled_forces = forces(i,:) .* (0.2 * 1/sqrt(forces(i,1)^2+forces(i,2)^2) * sqrt(sum(coords(1,:).^2)) );
            end
    
            disp("Force on particle "+i+" is ("+forces(i,1)+", "+forces(i,2)+")")
            drawArrow(coords(i,:), scaled_forces);
        end
        scatter(coords(:,1), coords(:,2))
        

    elseif length(coords(1,:)) == 3
        drawArrow = @(pos,force) quiver3( pos(1),pos(2),pos(3),force(1),force(2),force(3));
        
        scale = sum(abs(coords))/sum(abs(forces));
        particle_iter_max = length(coords(:,1));
        for i = 1:particle_iter_max
            if force_scaling_type == "scale"
                scaled_forces = forces(i,:) .* scale * 0.2;
            elseif force_scaling_type == "norm"
                scaled_forces = forces(i,:) .* (0.2 * 1/sqrt(forces(i,1)^2+forces(i,2)^2+forces(i,3)^2) * sqrt(sum(coords(1,:).^2)) );
            elseif force_scaling_type == "xy"
                forces(i,3) = 0;
                scaled_forces = forces(i,:) .* (0.2 * 1/sqrt(forces(i,1)^2+forces(i,2)^2+forces(i,3)^2) * sqrt(sum(coords(1,:).^2)) );
            
            end
            
            disp("Position of particle "+i+" is ("+coords(i,1)+", "+coords(i,2)+", "+coords(i,3)+")")
            % disp("Force on particle "+i+" is ("+forces(,1)+", "+forces(i,2)+", "+forces(i,3)+")")

            col = 10 * i/particle_iter_max;
            drawArrow(coords(i,:), scaled_forces);
            hold on
            scatter3(coords(i,1), coords(i,2), coords(i,3), [], col)
            
        end
        % scatter3(coords)
        zlim([min(coords(:,3))-2e-7, max(coords(:,3))+2e-7])
        zlabel("z /m")

    end

    xlim([min(coords(:,1))-2e-7, max(coords(:,1))+2e-7])
    ylim([min(coords(:,2))-2e-7, max(coords(:,2))+2e-7])
    xlabel("x /m")
    ylabel("y /m")
    title(figTitle)
    uuid = char(matlab.lang.internal.uuid());
    savefig("Figure"+uuid+".fig")
    hold off;
    
end

function output = Show_Particle_Forces(laguerreOrder, distance, radius, particle_iter_max)
    Nmax = 10;

    n_medium = 1.3;
    n_particle = 1.6;
    wavelength0_base = 1064e-9;         %Vacuum wavelength
    wavelength_medium_base = wavelength0_base / n_medium;
    wavelength0 = wavelength0_base;
    beam_polarisation = [ 1 1i ];
    NA = 1.4;
    beam = Get_Beam(laguerreOrder, n_medium, wavelength0, NA, beam_polarisation, Nmax);
   

    particles = Get_Particles("circle", particle_iter_max, distance, radius, 0.0, 0.0);    
    coords = zeros(particle_iter_max, 2);
    forces = zeros(particle_iter_max, 2);

    for focus_particle_iter = 1:particle_iter_max
        coords(focus_particle_iter, 1) = particles{focus_particle_iter}.position(1);
        coords(focus_particle_iter, 2) = particles{focus_particle_iter}.position(2);

        force = Perform_MST_Calc(particles, focus_particle_iter, beam, wavelength0, n_medium, n_particle);   %Calculate forces on ith particle in set
        force_scalar_x = Get_Specific_ForceScalar("X", force, particles{focus_particle_iter});
        force_scalar_y = Get_Specific_ForceScalar("Y", force, particles{focus_particle_iter});
        forces(focus_particle_iter, 1) = force_scalar_x;
        forces(focus_particle_iter, 2) = force_scalar_y;
    end
    
    Plot_Force_Arrows(coords, forces, "scale", radius, "Radius = "+radius+"m")
end

function output = Sweep_Particle_In_Planes_Forces(laguerreOrder, target_rad, variation_rad, rad_iter_max, target_dist, variation_dist, dist_iter_max, ang_iter_max, ang_max, z_min, z_max, z_iter_max)
    % Setup variables.
    Nmax = 10;
    n_medium = 1.0;
    n_particle = 1.55;
    focus_particle_iter = 1;
    particle_iter_max = 1;
    wavelength0_base = 1064e-9;         %Vacuum wavelength
    wavelength_medium_base = wavelength0_base / n_medium;
    wavelength0 = wavelength0_base;
    beam_polarisation = [ 1 1i ];
    NA = 1.3;
    beam = Get_Beam(laguerreOrder, n_medium, wavelength0, NA, beam_polarisation, Nmax);

        % %{
    radii = linspace(target_rad*(1-variation_rad), target_rad*(1+variation_rad),  rad_iter_max);
    distances = linspace(target_dist*(1-variation_dist), target_dist*(1+variation_dist), dist_iter_max);

    for rad_iter = 1:rad_iter_max
        if rad_iter_max == 1
            radius = target_rad;
        else
            radius = radii(rad_iter);
        end
        % radius = target_rad + ( rad_iter/rad_iter_max - 1/2 )*(variation_rad*target_rad);
        disp("Calculating for radius = "+radius)


        coords = zeros(dist_iter_max*ang_iter_max*z_iter_max, 3);
        forces = zeros(dist_iter_max*ang_iter_max*z_iter_max, 3);
        for z_iter = 1:z_iter_max
            z = z_min + (z_max-z_min)*(z_iter-1)/(z_iter_max);
    
            for dist_iter = 1:dist_iter_max
                if dist_iter_max == 1
                    distance = target_dist;
                else
                    distance = distances(dist_iter);
                end
    
                if radius > distance
                    disp("r, d = "+radius + ", " + distance)
                    radius = 0.9*distance;
                    disp("Radius > Distance, setting r = 0.9*d = "+radius)
                end
                    
                for offset_iter = 1:ang_iter_max
                    iter = (z_iter-1) * ang_iter_max * dist_iter_max + (dist_iter-1) * ang_iter_max + offset_iter;
                    offset = (offset_iter-1)/ang_iter_max * ang_max;
                    particles = Get_Particles("circle", particle_iter_max, distance, radius, z, offset);
                    assignin("base", "particles", particles)
                    coords(iter, 1) = particles{focus_particle_iter}.position(1);
                    coords(iter, 2) = particles{focus_particle_iter}.position(2);

                    force = Perform_MST_Calc(particles, focus_particle_iter, beam, wavelength0, n_medium, n_particle);   %Calculate forces on ith particle in set

                    force_scalar_x = Get_Specific_ForceScalar("X", force, particles{focus_particle_iter});
                    force_scalar_y = Get_Specific_ForceScalar("Y", force, particles{focus_particle_iter});

                    forces(iter, 1) = force_scalar_x;
                    forces(iter, 2) = force_scalar_y;

                    coords(iter, 3) = particles{focus_particle_iter}.position(3);
                    force_scalar_z = Get_Specific_ForceScalar("Z", force, particles{focus_particle_iter});
                    forces(iter, 3) = force_scalar_z;
                    
                end
            end
        end
        Plot_Force_Arrows(coords, forces, "norm", radius, "Radius = "+radius+"m") % scale, norm, xy
        
    end
end


function output = Sweep_Particle_Dist(view_range, laguerreOrder, target_dist, variation_dist, dist_iter_max)
    Nmax = 15;

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
    Nmax = 15;

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
    Nmax = 10;

    n_medium = 1.3;
    n_particle = 1.6;
    particle_iter_max = 1;              %Max number of particles in the ring (1->N)
    focus_particle_iter = 1;            %Which particle index to observe in measurements (1->N)
    wavelength0_base = 1064e-9;         %Vacuum wavelength
    wavelength_medium_base = wavelength0_base / n_medium;
    wavelength0 = wavelength0_base;
    beam_polarisation = [ 1 1i ];
    NA = 1.4;

    sweepData = {};     %Holds objects that are vectors [distance, radius, force on ith]
    beam = Get_Beam(laguerreOrder, n_medium, wavelength0, NA, beam_polarisation, Nmax);
    
    % show the beam
    figure();
    subplot(1, 1, 1);
    beam.visualise('axis', 'z');

    sweepDistRad_force_figure = figure;

    for dist_iter = 1:dist_iter_max
        disp("dist_iter= "+dist_iter+"/"+dist_iter_max);
        distance = target_dist +( dist_iter-(dist_iter_max/2.0) )*(variation_dist*target_dist);
        for rad_iter = 1:rad_iter_max
            disp("    rad_iter= "+rad_iter+"/"+rad_iter_max);
            radius = target_rad +( rad_iter-(rad_iter_max/2.0) )*(variation_rad*target_rad);

            particles = Get_Particles("circle", particle_iter_max, distance, radius, 0.0, 0.0);                       %Get a set of particles
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
    elseif(force_format == "Z")
        force_scalar = force(3);
    else
        %Do norm if nothing else specified
        disp("Invalid format, force norm being displayed...");
        force_scalar = vecnorm(force);
    end
end

function particles = Get_Particles(layout_type, particle_number, distance, radius, Z_plane, offset)
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
            pos = [distance*cos(theta*(particle_iter-1)+offset);distance*sin(theta*(particle_iter-1)+offset);Z_plane];
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

    if isscalar(particles) % checks length is 1
        % disp("Using beam offset for force calc.")
        shape = ott.shapes.Shape.simple("sphere", particles{1}.radius);
        T_union = ott.TmatrixMie.simple(shape, ...
        'wavelength0', wavelength0, ...
        'index_medium', n_medium, ...
        'index_particle', n_particle ...
        );
        % ibeam = ibeam.translateXyz(-particles{1}.position);
        ibeam = translateXyz(ibeam, -particles{1}.position);
        beam_total = totalField(T_union*ibeam, ibeam);
        force = ForceCalc_quad([0;0;0], particles{focus_particle_index}.radius, beam_total);
        
    else
        shape_union = UnionAllParticles(particles);
        T_union = ott.TmatrixMie.simple(shape_union, ...
           'wavelength0', wavelength0, ...
           'index_medium', n_medium, ...
           'index_particle', n_particle ...
        );
        beam_total = totalField(T_union*ibeam, ibeam);
        force = ForceCalc_quad(particles{focus_particle_index}.position, particles{focus_particle_index}.radius, beam_total);
    end

    %Find scattered field
    

    %Find force on specific particle in radial direction
    %force = ForceCalc_modified(particles{focus_particle_index}.position,
    %particles{focus_particle_index}.radius, beam_total, force_samples); %%OLD
    
end



disp("Program end");