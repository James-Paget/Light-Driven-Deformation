n_medium = 1.33;        % Water
n_particle = 1.59;      % Polystyrene
wavelength0 = 1064e-9;  % Vacuum wavelength
wavelength_medium = wavelength0 / n_medium;
NA = 1.02;              % Numerical aperture of beam
nrel = n_particle/n_medium;
beamPos = [0; 0; 0];
plotSize = 4e-6;

radius1 = 0.4e-6;
radius2 = 0.4e-6;
radius3 = 0.4e-6;
radius4 = 0.4e-6;
pos1 = [1.5e-6; 0; 0];
pos2 = [-1.5e-6; 0; 0];
pos3 = [0; 1.5e-6; 0];
pos4 = [0; -1.5e-6; 0];

radii = [radius1, radius2];
positions = [pos1, pos2];

% radii = [radius1, radius2, radius3, radius4];
% positions = [pos1, pos2, pos3, pos4];

union = ott.shapes.Union([ott.shapes.Shape.simple("sphere", radii(1))]);
if length(radii) > 1
    for i = 2:length(radii)
        shapei = ott.shapes.Shape.simple("sphere", radii(i));
        shapei.position = positions(:,i);
        union = ott.shapes.Union([union, shapei]);
    end
end


T = ott.TmatrixMie.simple(union, 'index_relative', nrel, ...
    'index_medium', 1.0, 'wavelength0', wavelength0);

lgmode = 3;
beam = ott.BscPmGauss('lg', [ 0 lgmode ], ...
     'polarisation', [ 1 1i ], 'NA', NA, ...
     'index_medium', n_medium, 'wavelength0', wavelength0);
beam.basis = 'regular';


% beam = beam.translateXyz(beamPos);
beam = beam + beam2.rotateXyz(3.1415926,0,0);

sbeam = T * beam;

% ### PLOTTING STARTS ###

plotField = false;
if plotField
    figure();
    subplot(1, 2, 1);
    sbeam.basis = 'outgoing';
    sbeam.visualise('axis', 'y', ...
       'mask', @(xyz) (vecnorm(xyz - pos1) < radius1) | (vecnorm(xyz - pos2) < radius2), ...
       'range', [1,1]*plotSize)
    title('Scattered field');
    
    subplot(1, 2, 2);
    tbeam = sbeam.totalField(beam);
    tbeam.basis = 'regular';
    tbeam.visualise('axis', 'y', ...
        'mask', @(xyz) (vecnorm(xyz - pos1) < radius1) | (vecnorm(xyz - pos2) < radius2), ...
        'range', [1,1]*plotSize)
    title('Total field');
end

plotBeam = true;
if plotBeam
    beam.basis = 'regular';
    figure();
    subplot(1, 2, 1);
    beam.visualise('axis', 'y');
    subplot(1, 2, 2);
    beam.visualise('axis', 'z');
end

plotForces = false;
if plotForces
    xyz = [0;0;1] .* linspace(-4, 4, 100).*1e-6;
    fxyz = ott.forcetorque(beam, T, 'position', xyz);
    figure();
    nPc = n_medium / 3e8;  % n / vacuum_speed
    plot(xyz(3, :), fxyz .* nPc);
    xlabel('Z position [m]');
    ylabel('Force [N/W]');
    legend({'Fx', 'Fy', 'Fz'});
end

surfForces = true;
z_intercept = 0e-6;
if surfForces
    gridNum = 30;
    linX = linspace(-plotSize, plotSize, gridNum);
    linY = linspace(-plotSize, plotSize, gridNum);
    linZ = linspace(0, 0, gridNum);
    [X, Y] = meshgrid(linX, linY);
    Zx = meshgrid(linZ, linZ);
    Zy = meshgrid(linZ, linZ);
    Zz = meshgrid(linZ, linZ);
    Zt = meshgrid(linZ, linZ);
    for x = 1:gridNum
        xyz = [0;1;0] .* linY + [linX(x);0;z_intercept];
        fxyz = ott.forcetorque(beam, T, 'position', xyz);
        Zx(x,:) = fxyz(1,:);
        Zy(x,:) = fxyz(2,:);
        Zz(x,:) = fxyz(3,:);
        Zt(x,:) = vecnorm(fxyz(:,:));
    end
       
    figure();
    subplot(2, 2, 1);
    s = surf(X,Y,Zx,'FaceAlpha',0.5);
    title('F_x');
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('F_x [N]');

    subplot(2, 2, 2);
    s = surf(X,Y,Zy,'FaceAlpha',0.5);
    title('F_y');xlabel('x [m]');
    ylabel('y [m]');
    zlabel('F_y [N]');

    subplot(2, 2, 3);
    s = surf(X,Y,Zz,'FaceAlpha',0.5);
    title('F_z');
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('F_z [N]');

    subplot(2, 2, 4);
    s = surf(X,Y,Zt,'FaceAlpha',0.5);
    title('F_{tot}');
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('F_{tot} [N]');

end


plotQuiver = true;
z_intercept = 0e-6;
if plotQuiver
    gridNum = 30;
    linX = linspace(-plotSize, plotSize, gridNum);
    linY = linspace(-plotSize, plotSize, gridNum);
    linZ = linspace(0, 0, gridNum);
    [X, Y] = meshgrid(linX, linY);
    Zx = meshgrid(linZ, linZ);
    Zy = meshgrid(linZ, linZ);
    Zz = meshgrid(linZ, linZ);
    Zt = meshgrid(linZ, linZ);
    Z0 = meshgrid(linZ, linZ);
    for x = 1:gridNum
        xyz = [0;1;0] .* linY + [linX(x);0;z_intercept];
        fxyz = ott.forcetorque(beam, T, 'position', xyz);
        Zx(x,:) = fxyz(1,:);
        Zy(x,:) = fxyz(2,:);
        Zz(x,:) = fxyz(3,:);
        Zt(x,:) = vecnorm(fxyz(:,:));
    end

    [U,V,W] = surfnorm(Zx,Zy,Zz);
       
    figure();
    quiver3(X,Y,Z0,U,V,W)
    title('F at each point in the plane');
    xlabel('x');
    ylabel('y');
    zlabel('F');
end