%{
This program has several options in a process you can follow in order to
help you find parameters suitable for trapping.
%}
disp("Program start");


%Startup
close all;
ott.warning('once');
ott.change_warnings('off');

tic
%%MODES % View_Beam, Time_Simulation, Time_Simulation_Multiple, DeltaXGrid, Scale_Beam
mode_sequence = "View_Beam";
for mode = mode_sequence
    if(mode == "View_Beam")
        disp("Performing View_Beam");
        view_range = 2.5e-6;
        minLg = 8;
        maxLg = 8;

        Nmaxs = [10,15];
        for Nmax = Nmaxs
            View_Beam(view_range, minLg, maxLg, Nmax);
        end
    
    elseif(mode == "Time_Simulation")
        disp("Performing Time_Simulation");
        d = 1e-6;
        angle = pi/6;
        initial_position = [d*cos(angle);d*sin(angle);0]; % [d;0;0]

        dt = 0.7e-4;
        t_max = 2*dt; % 150e-4
        radius = 0.5e-6;
        Nmax = 20;
        shouldLoad = false;
        shouldStartNew = true;
        laguerreOrder = 8;
        shouldUseDDA = false;
        figTitle = "radius = "+radius+ ", init pos = ("+ initial_position(1)+", "+initial_position(2)+", "+initial_position(3)+"), dt = "+dt+ ", Nmax = "+Nmax;

        Time_Simulation(initial_position, t_max, dt, radius, laguerreOrder, figTitle, Nmax, shouldLoad, shouldStartNew, shouldUseDDA)

    elseif(mode == "Time_Simulation_Multiple")
        disp("Performing Time_Simulation_Multiple");
        
        dt = 1e-4; % 1e-3
        t_max = 2*dt; 
        repeats = 1;
        
        laguerreOrder = 8;
        pos_angle = pi/6;
        Nmax = 15;

        shouldStartNew = true;
        shouldLoad = true;

        shouldUseDDA = false;

        distMin = 1.1e-6; % 0.9e-6
        distMax = 1.5e-6;
        distNum = 1;

        zMin = 0; % 0e-6 0.3e-6
        zMax = 1e-6;
        zNum = 1;

        radMin = 0.5e-6; % 0.35e-6
        radMax = 0.45e-6;
        radNum = 1;
    
        dispStr = "";
        zArr = InclusiveIterate_MinMax(zMin, zMax, zNum);
        distArr = InclusiveIterate_MinMax(distMin, distMax, distNum);
        radArr = InclusiveIterate_MinMax(radMin, radMax, radNum);
        
        for z_iter = 1:zNum
            z = zArr(z_iter);

            for dist_iter = 1:distNum
                distance = distArr(dist_iter);
                initial_position = [distance * cos(pos_angle); distance * sin(pos_angle); z];
    
                for rad_iter = 1:radNum
                    radius = radArr(rad_iter);
                    figNum = (z_iter-1) * radNum * distNum + (dist_iter-1) * radNum + rad_iter;
                    
                    titleStr = "Figure "+ figNum + " : radius = "+radius+ ", dt = "+dt+ ", Nmax = "+Nmax;
                    disp(titleStr)
                    dispStr = dispStr +"Figure "+ figNum + " : radius = "+radius+ ", distance = "+ distance + newline;
                    
                    disp(newline + "!!! --> shouldLoad is "+ shouldLoad)
                    disp("!!! --> shouldStartNew is "+ shouldStartNew)

                    for repeat = 1:repeats
                        Time_Simulation(initial_position, t_max, dt, radius, laguerreOrder, titleStr, Nmax, shouldLoad, shouldStartNew, shouldUseDDA)
                    end
                end
            end
        end
        disp(dispStr)

    elseif(mode == "DeltaXGrid")
        disp("Performing DeltaXGrid");
        
        xMax = 1.05e-6;
        xNum = 31;
        yMax = 1.05e-6;
        yNum = 31;
        zMin = 0.3e-6;
        zMax = 0.4e-6;
        zNum = 1;
        minDist = 0.7e-6;

        dt = 5e-3;
        radius = 0.35e-6;
        Nmax = 15;
        laguerreOrder = 8;
        shouldUseDDA = false;
        figTitle = "radius = "+radius+ ", dt = "+dt+ ", Nmax = "+Nmax;

        DeltaX_Grid(xMax, xNum, yMax, yNum, zMin, zMax, zNum, dt, radius, laguerreOrder, figTitle, Nmax, minDist, shouldUseDDA)
    
    elseif(mode == "Scale_Beam")
        disp("Performing Scale_Beam");

        dt = 5e-3;
        radius = 1e-8; % 0.35e-6
        figTitle = "radius = "+radius+ ", dt = "+dt+ ", Nmax = "+Nmax;

        plotSize = 4e-6;
        laguerreOrder = 8;

        Nmaxs = [10, 15];
        dispStr = "";
        shouldUseDDA = false;
        for Nmax = Nmaxs
            output = Scale_Beam(plotSize, laguerreOrder, Nmax, shouldUseDDA);
            dispStr = dispStr + output +newline;
        end
        disp(dispStr)


    end
end
toc

function output = Scale_Beam(plotSize, laguerreOrder, Nmax, shouldUseDDA)
    
    medium = "water"; % true:water, false: air
    if medium == "water"
        n_medium = 1.33;
        eta = 1e-3; % Pa s
        % figTitle = " In water, " + figTitle;
    else
        n_medium = 1.0;
        eta = 1.81e-5; % Pa s
        % figTitle = " In air, " + figTitle;
    end
    % disp("Calculating for: "+figTitle)

    n_medium = 1.0; 
    n_particle = 1.55;
    focus_particle_iter = 1;
    particle_iter_max = 1;
    wavelength0_base = 1064e-9;         %Vacuum wavelength
    wavelength_medium_base = wavelength0_base / n_medium;
    wavelength0 = wavelength0_base;
    beam_polarisation = [ 1 1i ];
    NA = 1.3;

    kb = 1.38e-23;
    temperature = 293;
    rng("default") % should keep the seed the same.
    scale_rand = 0e-1;
    
    xMaxScale = 2.5e-6;
    zMaxScale = 2.5e-6;
    
    beam = Get_Beam(laguerreOrder, n_medium, wavelength0, NA, beam_polarisation, Nmax);

    plotN = 100;
    reference_beam = Get_Beam(laguerreOrder, n_medium, wavelength0, NA, beam_polarisation, 10);
    global_scale_factor = Get_Beam_Scale(reference_beam, xMaxScale, zMaxScale);
    beam_scale = Get_Beam_Scale(beam, xMaxScale, zMaxScale);
    scale_factor = global_scale_factor / beam_scale;
    % disp("beam scale factor is "+scale_factor)
    Plot_Fields(beam, ["Z","Y"], [0,0], plotSize, plotN, scale_factor);
    % disp("Nmax is "+Nmax+" with scale "+beam_scale)

    particles = {};
    position = [1.0e-6;0;0];
    radius = 0.5e-6;
    particles{1} = ott.shapes.Sphere(radius, position);
    force = Perform_MST_Calc(particles, focus_particle_iter, beam, wavelength0, n_medium, n_particle, shouldUseDDA);   %Calculate forces on ith particle in set
    sforce = force * scale_factor;
    
    output = "Nmax is "+Nmax+" with force ("+ force(1)+", "+force(2)+", "+force(3)+") and scale "+scale_factor +newline+" so final ("+ sforce(1)+", "+sforce(2)+", "+sforce(3)+")";
    
    target_scale = 1e-13;
    force_scale = Get_Scale_From_Force(Nmax, target_scale, shouldUseDDA);
    fforce = force_scale * force;
    output = output + newline + "forceScale is " + force_scale + " gives (" + fforce(1)+", "+fforce(2)+", "+fforce(3)+")" + newline;
    disp(output)

end

function force_scale = Get_Scale_From_Force(Nmax, target_scale, shouldUseDDA)
    n_medium = 1.33;
    n_particle = 1.55;
    wavelength0 = 1064e-9;
    beam_polarisation = [ 1 1i ];
    NA = 1.3;
    laguerreOrder = 8;
    particles = {};
    position = [1.2e-6;0;0];
    radius = 0.5e-6;
    particles{1} = ott.shapes.Sphere(radius, position);
    beam = Get_Beam(laguerreOrder, n_medium, wavelength0, NA, beam_polarisation, Nmax);
    force = Perform_MST_Calc(particles, 1, beam, wavelength0, n_medium, n_particle, shouldUseDDA);   %Calculate forces on ith particle in set
    force_scale = target_scale / vecnorm(force);
end

function beam_scale = Get_Beam_Scale(field, xMaxScale, zMaxScale)
    f = setParamsBeamInteg(field);
    [result, ~] = quad2d(f,0, xMaxScale, 0, zMaxScale);
    beam_scale = result;
end

function f = setParamsBeamInteg(field)
    % return beamInteg with the parameters already set so f = f(Q1, Q2)
    f = @(X, Z) beamInteg(X, Z, field);
end

function integ = beamInteg(X, Z, field)
    sz = size(X);
    % disp(sz)
    integ = zeros(sz);
    for i = 1:sz(1)
        for j = sz(2)
            x = X(i,j);
            z = Z(i,j);
            [E, ~] = field.emFieldXyz([x;0;z]);
            integ(i,j) = vecnorm(E)^2; % square for intensity
        end
    end
end


function [X, Y, Zcol] = Plot_Fields(field, planeAxes, planeValues, plotSize, plotN, scale_factor)
    % Plots the magnitude of the field in the z=z_plane plane and returns
    % the X,Y,Z matrices used to make the surface. plotSize sets the x,y
    % limits and plotN gives the number of sampled points on each axis.
    % plane axis is "X", "Y" or "Z"
    numPlots = length(planeAxes);
    fig = figure();
    for plot_iter = 1:numPlots
        subplot(1, numPlots, plot_iter)
        planeAxis = planeAxes(plot_iter);
        planeValue = planeValues(plot_iter);

        linX = linspace(-plotSize, plotSize, plotN);
        linY = linspace(-plotSize, plotSize, plotN);
        Zcol = zeros(plotN, plotN);
        [X,Y] = meshgrid(linX, linY);
        
        positions = zeros(plotN^2,3);
        i=1;
        for x = linX
            for y = linY
                positions(i,:) = [x,y,planeValue];
                i = i+1;
            end
        end
    
        if planeAxis == "X"
             positions(:,[1,3]) = positions(:,[3,1]);
             % xlabel("y")
             % ylabel("z")
        elseif planeAxis == "Y"
            positions(:,[2,3]) = positions(:,[3,2]);
            % xlabel("x")
            % ylabel("y")
        elseif planeAxis == "Z"
            % xlabel("x")
            % ylabel("z")
        else
            disp("ERROR: planeAxis must be 'X', 'Y' or 'Z'.")
        end
        
        EMfield = field.emFieldXyz(positions.').';
        Zcol(:,:) = vecnorm(reshape(EMfield, [plotN, plotN,3]),2,3);
        % 
        % xlim([-plotSize, plotSize])
        % ylim([-plotSize, plotSize])
        % zlim([-plotSize, plotSize])
        surf(X, Y, Zcol * scale_factor)
        shading interp
        title("Field in "+planeAxis+" = "+planeValue+" plane")
        % View can set it to look down from above for a 2D plot.
        view(0,90)
    end
end



function output = Time_Simulation(initial_position, t_max, dt, radius, laguerreOrder, figTitle, Nmax, shouldLoad, shouldStartNew, shouldUseDDA)

    mediumIsWaterNotAir = true; % true:water, false: air
    if mediumIsWaterNotAir
        n_medium = 1.33;
        eta = 1e-3; % Pa s
        figTitle = " In water, " + figTitle;
    else
        n_medium = 1.0;
        eta = 1.81e-5; % Pa s
        figTitle = " In air, " + figTitle;
    end
    disp("Calculating for: "+figTitle)

    n_particle = 1.55;
    focus_particle_iter = 1;
    particle_iter_max = 1;
    wavelength0_base = 1064e-9;         %Vacuum wavelength
    wavelength_medium_base = wavelength0_base / n_medium;
    wavelength0 = wavelength0_base;
    beam_polarisation = [ 1 1i ];
    NA = 1.3;
    beam = Get_Beam(laguerreOrder, n_medium, wavelength0, NA, beam_polarisation, Nmax);
    mass = 2.7e3 * 4/3*pi*radius^3; % calcite density * volume (SI units)
    kb = 1.38e-23;
    temperature = 293;
    rng("default") % should keep the seed the same.
    scale_rand = 0e-1;    
    
    particles = {};
    num_times = floor(t_max/dt);
    times = linspace(0, t_max, num_times);
    coords = zeros(num_times, 3);
    forces = zeros(num_times, 3);
    
    % possible to start from the end of the last.
    if shouldLoad && ~shouldStartNew
        S = load('totCoordMatrix.mat');
        coordsFromFile = S.totCoords;
        initial_position = coordsFromFile(end,:)';
    end
    figTitle = figTitle + ", init pos = ("+ initial_position(1)+", "+initial_position(2)+", "+initial_position(3)+"), ";

    % Force scaling
    target_scale = 1e-13;
    force_scale = Get_Scale_From_Force(Nmax, target_scale, shouldUseDDA);
    disp("Force scale is "+force_scale)


    position = initial_position;

    for iter = 1:num_times
        % option to contrain the particle to z=0.
        % position(3) = initial_position(3);
        
        disp("Iteration "+iter+"/"+num_times)

        % Record values
        t = times(iter);
        particles{1} = ott.shapes.Sphere(radius, position);
        coords(iter, :) = particles{focus_particle_iter}.position;

        % Lorentz force calculation
        force = Perform_MST_Calc(particles, focus_particle_iter, beam, wavelength0, n_medium, n_particle, shouldUseDDA);   %Calculate forces on ith particle in set
        forcexyz = force_scale * force;

        % 0 = drag + Lorentz + Brownian
        % sqrt(2*zeta*kb*temperature) = 1.66e-15
        random_value = (2*rand(1)-1) * scale_rand;
        zeta = 6*pi*radius*eta;
        sum_forces = forcexyz + sqrt(2*zeta*kb*temperature)*random_value;

        forces(iter, :) = sum_forces;
        dx = dt/zeta * sum_forces;
        position = position + dx;

    end
    
    
    % option to load everything to plot.
    if shouldStartNew
        totCoords = coords;
        totForces = forces;
        save('totCoordMatrix.mat', 'totCoords');
        save('totForceMatrix.mat', 'totForces');
        

    elseif shouldLoad
        totCoords = [coordsFromFile; coords];
        save('totCoordMatrix.mat', 'totCoords');

        SF = load('totForceMatrix.mat');
        forcesFromFile = SF.totForces;
        totForces = [forcesFromFile; forces];
        save('totForceMatrix.mat', 'totForces');
    else
        totCoords = coords;
        totForces = forces;
    end


    figTitle = "steps " + length(totCoords) + figTitle;  
    save('coordMatrix.mat', 'coords'); % save after reading from old file.
    save('forceMatrix.mat', 'forces');
    force_scaling_type = "norm"; % Options: scale, norm, xy
    Plot_Force_Arrows(totCoords, totForces, force_scaling_type, figTitle)
end

function output = DeltaX_Grid(xMax, xNum, yMax, yNum, zMin, zMax, zNum, dt, radius, laguerreOrder, figTitle, Nmax, minDist, shouldUseDDA)

    mediumIsWaterNotAir = true; % true:water, false: air
    if mediumIsWaterNotAir
        n_medium = 1.33;
        eta = 1e-3; % Pa s
        figTitle = " In water, " + figTitle;
    else
        n_medium = 1.0;
        eta = 1.81e-5; % Pa s
        figTitle = " In air, " + figTitle;
    end

    n_particle = 1.55;
    focus_particle_iter = 1;
    particle_iter_max = 1;
    wavelength0_base = 1064e-9;         %Vacuum wavelength
    wavelength_medium_base = wavelength0_base / n_medium;
    wavelength0 = wavelength0_base;
    beam_polarisation = [ 1 1i ];
    NA = 1.3;
    beam = Get_Beam(laguerreOrder, n_medium, wavelength0, NA, beam_polarisation, Nmax);
    kb = 1.38e-23;
    temperature = 293;
    rng("default") % keep the seed the same.
    scale_rand = 0e-1; 
    
    particles = {};
    coords = zeros(xNum*yNum*zNum, 3);
    dxs = zeros(xNum*yNum*zNum, 3);

    xRange = InclusiveIterate_MinMax(-xMax, xMax, xNum);
    yRange = InclusiveIterate_MinMax(-yMax, yMax,  yNum);
    zRange = InclusiveIterate_MinMax(zMin, zMax, zNum);
    
    % check how many points are not near the origin.
    effNum = 0;
    minDistSquared = minDist^2;
    for z_iter = 1:zNum
        for x_iter = 1:xNum
            x = xRange(x_iter);

            % if x < 0
            %     continue
            % end

            for y_iter = 1:yNum
                y = yRange(y_iter);
                if x^2 + y^2 > minDistSquared
                    effNum = effNum + 1;
                end
            end
        end
    end
    
    iter = 0;
    for z_iter = 1:zNum
        z = zRange(z_iter);
        for x_iter = 1:xNum
            x = xRange(x_iter);

            % if x < 0
            %     continue
            % end

            for y_iter = 1:yNum
                y = yRange(y_iter);
                
                % Skip if too near the origin.
                if x^2 + y^2 < minDistSquared
                    continue
                end
                disp("Iteration "+iter+"/"+effNum+" :xyz = "+x+", "+y+", "+z)
                iter = iter + 1;

                position = [x;y;z];
                particles{1} = ott.shapes.Sphere(radius, position);
                coords(iter, :) = particles{focus_particle_iter}.position;
        
                % Lorentz force calculation
                force = Perform_MST_Calc(particles, focus_particle_iter, beam, wavelength0, n_medium, n_particle, shouldUseDDA);   %Calculate forces on ith particle in set
                forcexyz = Get_Specific_ForceScalar("XYZ", force, particles{focus_particle_iter});
                 
        
                % 0 = drag + Lorentz + Brownian
                % sqrt(2*zeta*kb*temperature) = 1.66e-15
                random_value = (2*rand(1)-1) * scale_rand;
                zeta = 6*pi*radius*eta;
                sum_forces = forcexyz + sqrt(2*zeta*kb*temperature)*random_value;

                dx = dt/zeta * sum_forces;
                dxs(iter, :) = dx;
            end
        end
    end

    force_scaling_type = "norm"; % Options: scale, norm, xy
    Plot_Force_Arrows(coords, dxs, force_scaling_type, figTitle) % note the print will name the dxs arrows as forces in the disp.
end



function output = View_Beam(view_range, laguerreOrder_target_lower, laguerreOrder_target_upper, Nmax)
    n_medium = 1.0; 
    n_particle = 1.55;
    wavelength0_base = 1064e-9;         %Vacuum wavelength
    wavelength_medium_base = wavelength0_base / n_medium;
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
        % beam = beam * 2;
        
        disp("---");
        disp("W; "+ subplot_width);
        disp("H; "+ subplot_height);
        disp("I: "+ subplot_index);
        subplot(subplot_width, subplot_height, subplot_index);
        beam.visualise('axis', 'z', ...
            'mask', @(xyz) vecnorm(xyz)<0.00*view_range, ...
            'range', [-1,1]*view_range );
        title("Laguerre beam order "+laguerreOrder_iter);
        subplot(subplot_width, subplot_height, subplot_index+1);
        beam.visualise('axis', 'x', ...
            'mask', @(xyz) vecnorm(xyz)<0.00*view_range, ...
            'range', [-1,1]*view_range );
        title("Laguerre beam order "+laguerreOrder_iter);
    end
    output=1;
end


function output = Plot_Force_Arrows(coords, forces, force_scaling_type, figTitle)
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
            disp("Force on particle "+i+" is ("+forces(i,1)+", "+forces(i,2)+", "+forces(i,3)+")")
        end


        for i = 1:particle_iter_max
            if force_scaling_type == "scale"
                scaled_forces = forces(i,:) .* scale * 0.2;
            elseif force_scaling_type == "norm"
                scaled_forces = forces(i,:) .* (0.06 * 1/sqrt(forces(i,1)^2+forces(i,2)^2+forces(i,3)^2) * sqrt(sum(coords(1,:).^2)) );
            elseif force_scaling_type == "xy"
                forces(i,3) = 0;
                scaled_forces = forces(i,:) .* (0.2 * 1/sqrt(forces(i,1)^2+forces(i,2)^2+forces(i,3)^2) * sqrt(sum(coords(1,:).^2)) );
            
            end
            
            disp("Position of particle "+i+" is ("+coords(i,1)+", "+coords(i,2)+", "+coords(i,3)+")")
           

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




function beam = Get_Beam(laguerreOrder, n_medium, wavelength0, NA, polarisation, Nmax)
    %{
    .laguerreOrder = order of beam

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


function force = Perform_MST_Calc(particles, focus_particle_index, ibeam, wavelength0, n_medium, n_particle, shouldUseDDA)
    if shouldUseDDA
        shape_union = UnionAllParticles(particles);
        spacing = wavelength0 / 10;
        voxels = shape_union.voxels(spacing, 'even_range', true);
        T_union = ott.TmatrixDda(voxels, ...
         'index_relative', n_particle/n_medium, ...
         'index_medium', n_medium, ...
         'spacing', spacing, ...    
         'z_rotational_symmetry', 1, ...
         'wavelength0', wavelength0, 'low_memory', true);
    
        beam_total = totalField(T_union*ibeam, ibeam);
        force = ForceCalc_quad(particles{focus_particle_index}.position, particles{focus_particle_index}.radius, beam_total);

    else
        if isscalar(particles) % checks length is 1
            
            shape = ott.shapes.Shape.simple("sphere", particles{1}.radius);
            T_union = ott.TmatrixMie.simple(shape, ...
            'wavelength0', wavelength0, ...
            'index_medium', n_medium, ...
            'index_particle', n_particle ...
            );
            ibeam = translateXyz(ibeam, -particles{1}.position);
            beam_total = totalField(T_union*ibeam, ibeam);
            force = ForceCalc_quad([0;0;0], particles{focus_particle_index}.radius, beam_total);
    
            % NSamples = 10;
            % force = ForceCalc_simpson([0;0;0], particles{focus_particle_index}.radius, beam_total, NSamples);
            
            
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
    end
    
end
    

function output = InclusiveIterate_MinMax(min, max, N)
    if N == 1
        output = min;
    else
        output = linspace(0,0,N);
        for iter = 1:N
            output(iter) = min + (iter-1)/(N-1) * (max - min);
        end
    end
end



disp("Program end");