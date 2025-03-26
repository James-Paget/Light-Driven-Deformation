%{
. Takes a list of particles, each of type ott.shapes.Shape.simple(...), and
returns a new shape that is the union of all these particles
%}
function output = UnionAllParticles(particles)
    particles_dim = size(particles);
    if(particles_dim(2) > 0)            %dim(2) is the number of columns (particles), dim(1) is the number of rows
        particle_union = particles{1};
        for index = 2:particles_dim(2)  %Note; starts from 2 as index 1 already accounted for
            particle_union = ott.shapes.Union([particle_union, particles{index}]);
        end
        output = particle_union;
    else
        disp("No particles to union; returning a simple sphere");
        output = ott.shapes.Shape.simple('sphere', 0.0);
    end
end