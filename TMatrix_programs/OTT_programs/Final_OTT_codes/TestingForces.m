% Testing the forces for different particle positions, beam positions,
% measured at different points.

close all;
Nmax = 15;
radius = 0.5e-6;
d = 1.9e-6;
% particle_poss = [[0;0;0], [d/2;0;0], [d;0;0], [0;0;0], [0;0;0] ];
% measure_poss = [[0;0;0], [d/2;0;0], [d;0;0], [d*3/4;0;0], [d;0;0]];
% beam_poss = [[0;0;0], [0;0;0], [0;0;0], [0;0;0],[0;0;0]];
particle_poss = [[1e-6;0;0], [0;0;0]];
measure_poss = [[1e-6;0;0], [0;0;0]];
beam_poss = [[0;0;0], [-1e-6;0;0]];
str = ["", "", "", "", ""];
disp("Nmax "+Nmax)
disp("For each listed position, tests what the force is there")

% disp("beam pos ("+ beam_pos(1)+", "+ beam_pos(2)+", "+ beam_pos(3)+")")

for i = 1:2
    particle_pos = particle_poss(:,i);
    disp("particle pos ("+ particle_pos(1)+", "+ particle_pos(2)+", "+ particle_pos(3)+")")
    measure_pos = measure_poss(:,i);
    beam_pos = beam_poss(:,i);
    beam = ott.BscPmGauss('lg', [ 0 8 ], ...
            'index_medium', 1.33, ...
            'wavelength0', 1064e-9, ...
            'NA', 1.3, ...
            'polarisation', [1 1i], ...
            'Nmax', Nmax ...
        );
    
    %{
    sumE = 0+0i;
    sumH = 0+0i;
    num=5;
    for angi = 1:num
        pos = measure_pos + [radius*cos(angi/num*2*pi);radius*sin(angi/num*2*pi);0];
        [E, H] = beam.emFieldXyz(pos);
        % disp("At ("+ pos(1)+", "+ pos(2)+", "+ pos(3)+"), E = ("+ E(1)+",   "+ E(2)+",   "+ E(3)+"), H= ("+ H(1)+",   "+ H(2)+",   "+ H(3)+")")
        % disp("At ("+ pos(1)+", "+ pos(2)+", "+ pos(3)+"), E = ("+ E(1)+"), H= ("+ H(1)+")")
        sumE = sumE + E;
        sumH = sumH + H;
    end
    % disp("avgs are "+sumE(1)/num+"   "+sumH(1)/num)
    %}

    shape = ott.shapes.Shape.simple("sphere", radius);
    shape.position = particle_pos;
    T = ott.TmatrixMie.simple(shape, ...
    'wavelength0', 1064e-9, ...
    'index_medium', 1.33, ...
    'index_particle', 1.55 ...
    );
    
    figure();
    beam.visualise('axis', 'z', 'range', [-1,1]*2.5e-6)
    scat_beam = T*beam;
    figure();
    scat_beam.visualise('axis', 'z', 'range', [-1,1]*2.5e-6)

    beam_total = totalField(T*beam, beam);
    figure();
    beam_total.visualise('axis', 'z', 'range', [-1,1]*2.5e-6)

    beam_total = translateXyz(beam_total, beam_pos);
    figure();
    beam_total.visualise('axis', 'z', 'range', [-1,1]*2.5e-6)

    force = ForceCalc_quad(measure_pos, radius, beam_total);
    disp(str(i)+"force is ("+ force(1)+", "+ force(2)+", "+ force(3)+")")

    if i == 2
        b2 = beam_total;
        % b2.visualise('axis', 'z', 'range', [-1,1]*2.5e-6, ...
            % 'mask', @(xyz) vecnorm(xyz)<0.0e-6)

        % beam.visualise('axis', 'z', 'range', [-1,1]*2.5e-6)
    end
end