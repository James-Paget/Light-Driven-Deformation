% Set up constants
n_medium = 1.33;        % Water
n_particle = 1.59;      % Polystyrene
wavelength0 = 1064e-9;  % Vacuum wavelength
wavelength_medium = wavelength0;
radius = 1.0*wavelength_medium;
NA = 1.0;              % Numerical aperture of beam
nrel = n_particle/n_medium;
plotSize = 4e-6;

% Define system.
radius1 = 0.5e-6;
radius2 = 0.5e-6;
pos1 = [1e-6; 0; 0];
pos2 = [-1e-6; 0; 0];
beamPos = [0; 0; 0];

% Create spheres
shape1 = ott.shapes.Shape.simple("sphere", radius1);
shape2 = ott.shapes.Shape.simple("sphere", radius2);
shape1.position = pos1;
shape2.position = pos2;
union = ott.shapes.Union([shape1, shape2]);

T = ott.TmatrixMie.simple(union, 'index_relative', nrel, ...
    'index_medium', 1.0, 'wavelength0', wavelength0);


beam = ott.BscPmGauss('lg', [ 0 8 ], ...
     'polarisation', [ 1 1i ], 'NA', NA, ...
     'index_medium', n_medium, 'wavelength0', wavelength0);

% beam = ott.BscPmGauss('NA', NA, 'polarisation', [ 1 1i ], ...
%      'index_medium', n_medium, 'wavelength0', wavelength0);
% beam.basis = 'regular';
beam = beam.translateXyz(beamPos);

sbeam = T * beam;

% ### PLOTTING STARTS ###

plotField = true;
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
    % %{
    subplot(1, 2, 1);
    beam.visualise('axis', 'y');
    subplot(1, 2, 2);
    beam.visualise('axis', 'z');
    % %}
    %{
    subplot(1, 1, 1);
    beam.visualise('axis', 'z');
    %}
end

plotForces = true;
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
    for x = 1:gridNum
        xyz = [0;1;0] .* linY + [linX(x);0;z_intercept];
        fxyz = ott.forcetorque(beam, T, 'position', xyz);
        Zx(x,:) = fxyz(1,:);
        Zy(x,:) = fxyz(2,:);
        Zz(x,:) = fxyz(3,:);
        Zt(x,:) = vecnorm(fxyz(:,:));
    end

    [U,V,W] = surfnorm(X,Y,Zt);
       
    figure();
    s = surf(X,Y,Zx,'FaceAlpha',0.5);
    title('F_x');
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('F_x [N]');
end
