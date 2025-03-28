%{
This program has several options in a process you can follow in order to
help you find parameters suitable for trapping.

.(1) "View_Beam" = view different orders of laguerre beam, see where fringes
occur, swept through beam orders
.(2) "View_Particle_Dist" = Displays the particles in the beam, positions swept
through various radi
.(3) "View_Particle_Rad" = Displays particle at given position now, with
varying radius swept through
.(4) "Sweep_Particle_Dist" = Fixes all other values, sweeps pos only, gives
force on 0th particle in this configuration, manually check for trapping
.(5) "Sweep_Particle_Rad" = Fixes all other values, sweeps rad only, gives
force on 0th particle in this configuration, manually check for trapping
.(6) "Sweep_Particle_Dist_Rad" = Fixes all other values, sweeps pos AND rad only, gives
force on 0th particle in this configuration, manually check for trapping

%}
import java.awt.Robot;
import java.awt.event.*;

disp("Program start");

robot = Robot();

%%
%% Tmatrix = ott.TmatrixDda.simple('sphere', 0.1, 'index_relative', 1.2);
%%
%% - translation_method -- Method to use when calculating translations.
      %     Can either be 'Default' or 'NewBeamOffset', the latter calculates
      %     new beam shape coefficients for every new position.

%Startup
close all;
ott.warning('once');
ott.change_warnings('off');

%Queue up several tests in one run
mode_sequence = ["View_Beam"];%, "Sweep_Particle_Dist", "Sweep_Particle_Rad", "Sweep_Particle_Dist_Rad"];
%%MODES
for mode = mode_sequence
    if(mode == "View_Beam")
        disp("Performing View_Beam");
        view_range = 3.5e-6;
        View_Beam(view_range, 8, 8, true);
    
    elseif(mode == "Sweep_Particle_Dist")
        disp("Performing Sweep_Particle_Dist");
        view_range      = 3.0e-6;
        laguerreOrder   = 8;
        radius          = (1e-6)/2.0;
        target_dist     = 0.5e-6;
        Sweep_Particle_Dist(view_range, laguerreOrder, radius, target_dist, 0.4, 10);
    
    elseif(mode == "Sweep_Particle_Rad")
        disp("Performing Sweep_Particle_Rad");
        view_range      = 4.0e-6;
        laguerreOrder   = 8;
        distance        = 0.5e-6;
        target_rad      = (1.0e-6)/2.0;%2.4e-7;
        Sweep_Particle_Rad(view_range, laguerreOrder, distance, target_rad, 0.2, 10);
    
    elseif(mode == "Sweep_Particle_Dist_Rad")
        disp("Performing Sweep_Particle_Dist_Rad");
        laguerreOrder   = 8;
        target_dist     = 0.5e-6;
        target_rad      = (1.0e-6)/2.0;     %Diameter 1e-6
        %Sweep_Particle_Dist_Rad(laguerreOrder, target_dist, 0.5, 10, target_rad, 0.5, 10, "X");
        Sweep_Particle_Dist_Rad(laguerreOrder, target_dist, 0.1, 10, target_rad, 0.1, 10, "Y");
    
    elseif(mode == "All_Particles_R_Theta")
        disp("Performing All_Particles_R_Theta");
        laguerreOrder   = 8;
        particle_rad    = (1.0e-6)/2.0;    %Diameter 1e-6
        particle_dist   = 1.1e-6;
        num_particles   = 1;
        All_Particles_R_Theta(num_particles, laguerreOrder, particle_dist, particle_rad);
    elseif(mode == "All_Particles_R_Theta_Singular")
        disp("Performing All_Particles_R_Theta_Singular");
        laguerreOrder   = 8;
        particle_rad    = (1.0e-6)/2.0;    %Diameter 1e-6
        particle_dist   = 0.0;%1.1e-6;
        num_particles   = 1;
        %dist_iter_max, laguerreOrder, distance_target, distance_jump, radius
        All_Particles_R_Theta_Singular(20, laguerreOrder, particle_dist, 0.05e-6, particle_rad, false);
    elseif(mode == "All_Particles_R_Theta_Vary_Radius")
        disp("Performing All_Particles_R_Theta_Vary_Radius");
        laguerreOrder   = 8;
        particle_rad    = (1.0e-6)/2.0;    %Diameter 1e-6
        particle_dist   = 1.1e-6;
        %iter_max, laguerreOrder, target, jump, distance
        All_Particles_R_Theta_Vary_Radius(10, laguerreOrder, particle_rad, 0.05e-6, particle_dist);


    elseif(mode == "SingleParticle_Variations")
        disp("Performing SingleParticle_Variations");
        Nmax            = 10;
        n_medium_base   = 1.0;              %Low nMedium => converges with lower Nmax,  High nMedium => converges for larger Nmax
        n_particle      = 1.55;         %calcite~1.55
        wavelength0_base= 1064e-9;      %Vacuum wavelength
        laguerreOrder   = 8;
        radius_base     = (1.0e-6)/2.0; %Diameter 1e-6
        %particle_dist   = 1.1e-6;
        space_data = [ ...
            -2.2e-6, 2.2e-6, 5; ...     %2.2
            -2.2e-6, 2.2e-6, 5; ...
            -1.0e-6, 1.0e-6, 3 ...
        ];
        
        N = 0;
        %(1) Vary radius
        for var_iter = -N:N  %-N
            radius = radius_base +var_iter*0.1*radius_base;
            SingleParticle_Force_3dSpace("RAD="+var_iter, Nmax, laguerreOrder, n_medium_base, n_particle, wavelength0_base, radius, space_data, robot);
        end
        %{
        %(2) Vary n_medium
        for var_iter = -N:N 
            n_medium = n_medium_base +var_iter*0.1;
            SingleParticle_Force_3dSpace("N_MEDIUM="+var_iter, Nmax, laguerreOrder, n_medium, n_particle, wavelength0_base, radius_base, space_data, robot);
        end
        %(3) Vary wavelength
        for var_iter = -N:N
            wavelength0 = wavelength0_base +var_iter*0.1*wavelength0_base;
            SingleParticle_Force_3dSpace("LAMBDA="+var_iter, Nmax, laguerreOrder, n_medium_base, n_particle, wavelength0, radius_base, space_data, robot);
        end
        %}
    end
end


function output = WakeProgram(moveSwitch, robot)
    %{
    . Uncomment the lines below to move your mouse each MST calculation,
    which should prevent the computer sleeping (given that each MST calculation ~< 10 minutes)
    %}

    %{
    if(moveSwitch)
        robot.mouseMove(50,0);
    else
        robot.mouseMove(100,0);
    end
    %}
    output = ~moveSwitch;
end


function output = SingleParticle_Force_3dSpace(plotID, Nmax, laguerreOrder, n_medium, n_particle, wavelength0, radius, space_data, robot)
    %{
    space_data = [
        x_start, x_end, x_samplePoints;
        y_start, y_end, y_samplePoints;
        z_start, z_end, z_samplePoints;
    ]

    . Want to try to achieve trapping with a single particle
    . Parameters we know;
        - Beam NA, order, polarisation
        - Rough particle radius, rough particle position, refractive index
        of particle
    . To do this, we should test;
        - Particle inside a beam of varying diameter    | Done through loop over this function
            - Vary particle radius and distance from origin  | DONE in this function
            - Vary Z height in the beam                      |
        - Vary refractive index of medium               | Done through loop over this function
    %}

    %Wake computer
    moveSwitch = true;

    %Parameters
    beam_polarisation = [ 1 1i ];
    NA = 1.3;

    %Get beam
    beam = Get_Beam(laguerreOrder, n_medium, wavelength0, NA, beam_polarisation, Nmax);

    view_range = 1.2*space_data(1,2);    %1.2*space_data(1,2);
    %beam_local_check_figure = figure;
    fig = figure('Visible','on');
    subplot(1,2,1);
    beam.visualise('axis', 'z', ...
        'mask', @(xyz) vecnorm(xyz)<0.05*view_range, ...
        'range', [-1,1]*view_range );
    title("laguerre order="+laguerreOrder+",  Nmax="+Nmax);
    subplot(1,2,2);
    beam.visualise('axis', 'x', ...
        'mask', @(xyz) vecnorm(xyz)<0.05*view_range, ...
        'range', [-1,1]*view_range );
    %Save plot
    filename = "[ID: "+plotID+"] beamLocalCheck.fig";
    fullFilePath = fullfile("/home/james/Desktop/Light-Driven-Deformation/TMatrix_programs/OTT_programs/Trapping_Sweeper/SingleParticle_TrappingSweep_Figures", filename);
    saveas(fig, fullFilePath);
    close(fig);

    %Sample space
    disp("Looking over sample space...");
    x_jump = ( space_data(1,2)-space_data(1,1) )/(space_data(1,3) -1);
    y_jump = ( space_data(2,2)-space_data(2,1) )/(space_data(2,3) -1);
    z_jump = ( space_data(3,2)-space_data(3,1) )/(space_data(3,3) -1);
    forces = {};    %[ [Z Layer 0], [Z Layer 1], [Z Layer 2], ... ]
    for k = 1:space_data(3,3)
        disp("k= "+k);
        z_pos = space_data(3,1) +(k-1)*z_jump;
        forces_2d = {};
        for j = 1:space_data(2,3)
            disp("  j= "+j);
            y_pos = space_data(2,1) +(j-1)*y_jump;
            for i = 1:space_data(1,3)
                disp("    i= "+i);
                x_pos = space_data(1,1) +(i-1)*x_jump;
                particle = ott.shapes.Sphere(radius, [x_pos;y_pos;z_pos]);
                force = SingleParticle_Perform_MST_Calc(particle, beam, wavelength0, n_medium, n_particle);
                forces_2d{j, i} = force;
                moveSwitch = WakeProgram(moveSwitch, robot);
            end
        end
        forces{1, k} = forces_2d;
        forces{2, k} = z_pos;
    end

    %Plot as series of 2D surfaces
    disp("Plotting figures...");
    x_set = linspace(space_data(1,1), space_data(1,2), space_data(1,3));
    y_set = linspace(space_data(2,1), space_data(2,2), space_data(2,3));
    [X,Y] = meshgrid(x_set, y_set);
    %layered_XY_force_figure = figure;
    for z_layer = 1:size(forces,2)
        forces_2d = forces{1, z_layer};
        XY_FMag = cellfun( @(cell) sqrt( (cell(1))^2 + (cell(2))^2 ), forces_2d);  %XY magnitude of force, NOT including any Z comp.
        U = cellfun( @(cell) cell(1), forces_2d );
        V = cellfun( @(cell) cell(2), forces_2d );
        
        fig = figure('Visible','on');
        subplot(1,2,1);
        quiver(X,Y,U,V);
        %surf(X,Y,Z);
        xlabel("X");
        ylabel("Y");
        title("[ID: "+plotID+"] Forces in z="+forces{2, z_layer}+" (shifted 1 particle)");

        subplot(1,2,2);
        %quiver(X,Y,U,V);
        surf(X,Y,XY_FMag);
        xlabel("X");
        ylabel("Y");
        zlabel("XY Force Magnitude");
        title("[ID: "+plotID+"] Forces in z="+forces{2, z_layer}+" (shifted 1 particle)");

        %Save plot
        filename = plotID+"_XYLayerForce_"+z_layer+".fig";
        fullFilePath = fullfile("/home/james/Desktop/Light-Driven-Deformation/TMatrix_programs/OTT_programs/Trapping_Sweeper/SingleParticle_TrappingSweep_Figures", filename);
        saveas(fig, fullFilePath);
        close(fig);
    end

    output=1;
end

function force = SingleParticle_Perform_MST_Calc(particle, ibeam, wavelength0, n_medium, n_particle)
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

    %This T-Matrix gets recentered at the origin
    T_union = ott.TmatrixMie.simple(particle, ...
       'wavelength0', wavelength0, ...
       'index_medium', n_medium, ...
       'index_particle', n_particle ...
    );
    %This moves the beam to acount for offset instead
    beam_offset = translateXyz(ibeam, -particle.position);

    %Find scattered field
    beam_total = totalField(T_union*beam_offset, beam_offset);

    %Find force on specific particle in radial direction
    %Particle at origin now => consider about (0,0,0) NOT particle pos
    force = ForceCalc_quad([0;0;0], particle.radius, beam_total);
end








%####
%#CAN CAUSE CONFUSION, BETTER TO VISUALISE BEAM LOCALLY IN XY AND XZ PLANE
%####
function output = View_Beam(view_range, laguerreOrder_target_lower, laguerreOrder_target_upper, fix_view)
    Nmax =200;

    n_medium = 1.0;
    n_particle = 1.55;                  %calcite~1.55
    wavelength0_base = 1064e-9;         %Vacuum wavelength
    %wavelength_medium_base = wavelength0_base / n_medium;
    wavelength0 = wavelength0_base;
    beam_polarisation = [ 1 1i ];
    NA = 1.3;

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
        if(fix_view)
            beam.visualise('axis', 'z', ...
                'mask', @(xyz) vecnorm(xyz)<0.05*view_range, ...
                'range', [-1,1]*view_range );
            title("Laguerre beam order "+laguerreOrder_iter);
            subplot(subplot_width, subplot_height, subplot_index+1);
            beam.visualise('axis', 'x', ...
                'mask', @(xyz) vecnorm(xyz)<0.05*view_range, ...
                'range', [-1,1]*view_range );
            title("Laguerre beam order "+laguerreOrder_iter);
        else
            beam.visualise('axis', 'z');
            title("Laguerre beam order "+laguerreOrder_iter);
            subplot(subplot_width, subplot_height, subplot_index+1);
            beam.visualise('axis', 'x');
            title("Laguerre beam order "+laguerreOrder_iter);
        end
    end
    output=1;
end

function output = All_Particles_R_Theta(particle_iter_max, laguerreOrder, distance, radius)
    % OLD ARGS: target_dist, variation_dist, dist_iter_max, target_rad, variation_rad, rad_iter_max


    % Check if particles will overlap and make radii smaller if needed.
    particle_min_gap_factor = 0.9;
    r_max = distance * sin(pi/particle_iter_max) * (1-particle_min_gap_factor);
    if radius < r_max
        radius = r_max;
        disp("All_Particles_R_Theta: changed radius to "+radius);
    end

    Nmax = 15;  %20

    n_medium = 1.0;
    n_particle = 1.55;
    wavelength0_base = 1064e-9;         %Vacuum wavelength
    wavelength_medium_base = wavelength0_base / n_medium;
    wavelength0 = wavelength0_base;
    beam_polarisation = [ 1 1i ];
    NA = 1.3;

    % sweepData = {};     %Holds objects that are vectors [distance, radius, force on ith]
    focusses = linspace(1, particle_iter_max, particle_iter_max);
    rs = linspace(0, 0, particle_iter_max);
    thetas = linspace(0, 0, particle_iter_max);
    beam = Get_Beam(laguerreOrder, n_medium, wavelength0, NA, beam_polarisation, Nmax);

    beam_localCheck_figure = figure;
    focus_r_theta_figure = figure;

    particles = Get_Particles("circle", particle_iter_max, distance, radius, 0.0); %Get a set of particles
    
    figure(beam_localCheck_figure);
    beam.visualise( ...
        'axis', 'z', ...
        'mask', @(xyz) vecnorm(xyz-particles{1}.position)<particles{1}.radius ...
    );
    title("Local Beam Check");

    figure(focus_r_theta_figure);
    
    for focus_particle_iter = 1:particle_iter_max
        disp("Focussing particle " + focus_particle_iter);                   
        force = Perform_MST_Calc(particles, focus_particle_iter, beam, wavelength0, n_medium, n_particle);   %Calculate forces on ith particle in set
        force_rhat = Get_Specific_ForceScalar("r_hat", force, particles{focus_particle_iter});
        force_scalar = Get_Specific_ForceScalar("theta_hat", force, particles{focus_particle_iter});
        
        focusses(focus_particle_iter) = focus_particle_iter;
        rs(focus_particle_iter)       = force_rhat;
        thetas(focus_particle_iter)   = force_scalar;
        

    end
    
    hold on;
    if(size(particles, 2) > 1)
        plot(focusses, rs)
        plot(focusses, thetas)
        legend('rs','thetas')
    else
        disp("=== RESULTS ===");
        disp("    r_hat force= "+rs(1));
        disp("theta_hat force= "+thetas(1));
        disp("");
        disp("theta / r comp= "+thetas(1)/rs(1));
        disp("===============");
    end
    % Display_3D_Sweep(sweepDistRad_force_figure, target_dist, variation_dist, dist_iter_max, target_rad, variation_rad, rad_iter_max, "distance (m)", "radius (m)", "force on "+focus_particle_iter, sweepData);     %Display force with DIST and RAD sweeping
end

function output = All_Particles_R_Theta_Singular(dist_iter_max, laguerreOrder, distance_target, distance_jump, radius, z_vary)
    % OLD ARGS: target_dist, variation_dist, dist_iter_max, target_rad, variation_rad, rad_iter_max

    Nmax = 10;  %20

    n_medium = 1.0;
    n_particle = 1.55;
    wavelength0_base = 1064e-9;         %Vacuum wavelength
    wavelength_medium_base = wavelength0_base / n_medium;
    wavelength0 = wavelength0_base;
    beam_polarisation = [ 1 1i ];
    NA = 1.3;

    beam_localCheck_figure = figure;
    force_ratio_figure = figure;
    force_comps_figure = figure;
    dists  = [];
    r_comp = [];
    ratios = [];
    theta_comp = [];
    for iter = -dist_iter_max:dist_iter_max
        disp("  iter= "+iter);
        distance = distance_target +iter*distance_jump;
        if(z_vary)
            distance = distance_target;
        end
        beam = Get_Beam(laguerreOrder, n_medium, wavelength0, NA, beam_polarisation, Nmax);
        particles = Get_Particles("circle", 1, distance, radius, 0.0); %Get a set of particles
        if(z_vary)
            particles{1}.position(3) = iter*distance_jump;
        end
        
        figure(beam_localCheck_figure);
        beam.visualise( ...
            'axis', 'z', ...
            'mask', @(xyz) vecnorm(xyz-particles{1}.position)<particles{1}.radius ...
        );
        title("Local Beam Check");

        %Get forces
        force = Perform_MST_Calc(particles, 1, beam, wavelength0, n_medium, n_particle);   %Calculate forces on ith particle in set
        force_rhat = Get_Specific_ForceScalar("r_hat", force, particles{1});
        force_thetahat = Get_Specific_ForceScalar("theta_hat", force, particles{1});

        if(z_vary)
            dists(iter +dist_iter_max +1)       = particles{1}.position(3);
        else
            dists(iter +dist_iter_max +1)       = distance;
        end
        r_comp(iter +dist_iter_max +1)      = force_rhat;
        theta_comp(iter +dist_iter_max +1)  = force_thetahat;
        ratios(iter +dist_iter_max +1)      = abs(force_thetahat/force_rhat);
    end

    assignin("base", "r_comp", r_comp);
    assignin("base", "theta_comp", theta_comp);
    assignin("base", "ratios", ratios);

    figure(force_ratio_figure);
    plot(dists, ratios);
    if(z_vary)
        xlabel("dist in Z");
    else
        xlabel("dist from origin (XY)");
    end
    ylabel("ratio");
    title("Plot of theta/r component ratios of force");

    figure(force_comps_figure);
    hold on;
    plot(dists, theta_comp);
    plot(dists, r_comp);
    if(z_vary)
        xlabel("dist in Z");
    else
        xlabel("dist from origin (XY)");
    end
    ylabel("comp");
    title("Plot of theta and r components of force");
    legend("theta", "r");
    hold off;

    output=1;
end

function output = All_Particles_R_Theta_Vary_Radius(iter_max, laguerreOrder, target, jump, distance)
    % OLD ARGS: target_dist, variation_dist, dist_iter_max, target_rad, variation_rad, rad_iter_max

    Nmax = 10;  %20

    n_medium = 1.0;
    n_particle = 1.55;
    wavelength0_base = 1064e-9;         %Vacuum wavelength
    wavelength_medium_base = wavelength0_base / n_medium;
    wavelength0 = wavelength0_base;
    beam_polarisation = [ 1 1i ];
    NA = 1.3;

    beam_localCheck_figure = figure;
    force_ratio_figure = figure;
    force_comps_figure = figure;
    variation = [];
    r_comp = [];
    ratios = [];
    theta_comp = [];
    for iter = -iter_max:iter_max
        disp("  iter= "+iter);
        radius = target +iter*jump;
        beam = Get_Beam(laguerreOrder, n_medium, wavelength0, NA, beam_polarisation, Nmax);
        particles = Get_Particles("circle", 1, distance, radius, 0.0); %Get a set of particles
        
        figure(beam_localCheck_figure);
        beam.visualise( ...
            'axis', 'z', ...
            'mask', @(xyz) vecnorm(xyz-particles{1}.position)<particles{1}.radius ...
        );
        title("Local Beam Check");

        %Get forces
        force = Perform_MST_Calc(particles, 1, beam, wavelength0, n_medium, n_particle);   %Calculate forces on ith particle in set
        force_rhat = Get_Specific_ForceScalar("r_hat", force, particles{1});
        force_thetahat = Get_Specific_ForceScalar("theta_hat", force, particles{1});

        variation(iter +iter_max +1)   = radius;
        r_comp(iter +iter_max +1)      = force_rhat;
        theta_comp(iter +iter_max +1)  = force_thetahat;
        ratios(iter +iter_max +1)      = abs(force_thetahat/force_rhat);
    end

    figure(force_ratio_figure);
    plot(variation, ratios);
    xlabel("radius (at fixed dist)");
    ylabel("ratio");
    title("Plot of theta/r component ratios of force");

    figure(force_comps_figure);
    hold on;
    plot(variation, theta_comp);
    plot(variation, r_comp);
    xlabel("radius (at fixed dist)");
    ylabel("comp");
    title("Plot of theta and r components of force");
    legend("theta", "r");
    hold off;

    output=1;
end

function output = Sweep_Particle_Dist(view_range, laguerreOrder, radius, target_dist, variation_dist, dist_iter_max)
    Nmax = 90;

    n_medium = 1.0;
    n_particle = 1.55;
    wavelength0_base = 1064e-9;
    %wavelength_medium_base = wavelength0_base / n_medium;
    wavelength0 = wavelength0_base;
    beam_polarisation = [ 1 1i ];
    NA = 1.3;

    beam = Get_Beam(laguerreOrder, n_medium, wavelength0, NA, beam_polarisation, Nmax);

    dist_sweep_figure = figure;
    figure(dist_sweep_figure);
    subplot_width  = ceil(sqrt(dist_iter_max));
    subplot_height = ceil(sqrt(dist_iter_max));
    for dist_iter = 1:dist_iter_max
        distance = target_dist +( dist_iter-(dist_iter_max/2.0) )*(variation_dist*target_dist);
        particle_pos    = [distance; 0; 0];

        subplot(subplot_width, subplot_height, dist_iter);
        beam.visualise('axis', 'z', ...
            'mask', @(xyz) vecnorm(xyz-particle_pos)<radius, ...
            'range', [-1,1]*view_range );
        title("Distance= "+distance);
    end
    output=1;
end

function output = Sweep_Particle_Rad(view_range, laguerreOrder, distance, target_rad, variation_rad, rad_iter_max)
    Nmax = 90;

    n_medium = 1.3;
    n_particle = 1.6;
    wavelength0_base = 1064e-9;
    wavelength_medium_base = wavelength0_base / n_medium;
    wavelength0 = wavelength0_base;
    beam_polarisation = [ 1 1i ];
    NA = 1.3;

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
    Nmax = 90;

    n_medium = 1.3;
    n_particle = 1.6;
    particle_iter_max = 1;              %Max number of particles in the ring (1->N)
    focus_particle_iter = 1;            %Which particle index to observe in measurements (1->N)
    wavelength0_base = 1064e-9;         %Vacuum wavelength
    wavelength_medium_base = wavelength0_base / n_medium;
    wavelength0 = wavelength0_base;
    beam_polarisation = [ 1 1i ];
    NA = 1.3;

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
    - 'lg' -- Laguerre-Gauss  [ radial azimuthal ]
    .laguereeOrder = order of beam

    . Add option to vary beams
    %}
    beam = ott.BscPmGauss('lg', [ 0 laguerreOrder ], ...
        'polarisation', polarisation, ...
        'NA', NA, ...
        'index_medium', n_medium, ...
        'wavelength0', wavelength0, ...
        'Nmax', Nmax ...
    );
    beam.basis = 'regular';
end

function force_scalar = Get_Specific_ForceScalar(force_format, force, particle)
    %{
    .force_format = the specific type of component you want to be extracted
    from this net force
    .force = vector force that you want to pull a scalar from
    %}
    force_scalar = 0.0;
    if(force_format == "theta_hat")
        r_mag = vecnorm(particle.position);
        r_hat_vec = particle.position ./ r_mag;                %( x,y,z(=0))
        theta_hat_vec = [-r_hat_vec(2), r_hat_vec(1), 0.0];    %(-y,x,0)
        force_scalar = force(1)*theta_hat_vec(1) + force(2)*theta_hat_vec(2) + force(3)*theta_hat_vec(3);
        %{
        force_angle = atan2(force(2), force(1));
        pos_angle = atan2(particle.position(2), particle.position(1));
        force_scalar = sin(force_angle - pos_angle);
        %}
    elseif(force_format == "r_hat")
        r_mag = vecnorm(particle.position);
        r_hat_vec = particle.position ./ r_mag;                %(x,y,z(=0))
        force_scalar = force(1)*r_hat_vec(1) + force(2)*r_hat_vec(2) + force(3)*r_hat_vec(3);
        %{
        force_angle = atan2(force(2), force(1));
        pos_angle = atan2(particle.position(2), particle.position(1));
        force_scalar = cos(force_angle - pos_angle);
        %}
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
    %BUG FIXING
    %assignin("base", "particles", particles);
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
    %%
    %% TEST THIS WITH PARTICLE AND BEAM OFFSET INSTEAD
    %%

    %Find scattered field
    beam_total = totalField(T_union*ibeam, ibeam);

    %Find force on specific particle in radial direction
    %force = ForceCalc_modified(particles{focus_particle_index}.position,
    %particles{focus_particle_index}.radius, beam_total, force_samples); %%OLD
    force = ForceCalc_quad(particles{focus_particle_index}.position, particles{focus_particle_index}.radius, beam_total);
end


disp("Program end");