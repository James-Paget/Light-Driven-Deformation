%{
. Considers a particle and a force field, sums the forces for a rigid total
force

. particle = shape forces being considered for
. beam = Incident beam acting on particle
. T = T-Matrix describing the ENTIRE setup (not just this shape)
. samples = Number of sections each axis bounding the particle is split
into in order to evaluate a grid of forces

%}
function output = FetchParticleNetForce(particle, beam, T, samples)
    bounding_width = particle.radius;
    x_points = particle.position(1) +[1;0;0] .* linspace(-bounding_width, bounding_width, samples);
    y_points = particle.position(2) +[0;1;0] .* linspace(-bounding_width, bounding_width, samples);
    %z_points = particle.position(3) +[0;0;1] .* linspace(-bounding_width, bounding_width, samples);
    z_plane = 0;
    force_net = [0;0;0];
    for i = 1:samples
        for j = 1:samples
            point = [x_points(1, i); y_points(2, j); z_plane];
            if( insideXyz(particle, point) )
                fxyz = ott.forcetorque(beam, T, 'position', point);
                force_net = force_net +fxyz;        %Vectorised comp. sum
            end
        end
    end

    %Divide to get average force (else would tend to infinity as samples tends to infinity)
    output = force_net/(samples^3);
end