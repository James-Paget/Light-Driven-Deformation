%{
. Is parsed the positions on a display window and a set of particles, and
returns true for any position (xyz) within a particle's bounding area for
all particles in the set
. This allows a mask to be made for an arbitrary number of particles
%}
function output = CheckParticlesMaskCondition(xyz, particles)
    %Given explicitly to give 'withinAnyParticle' the correct structure to
    %handle xyz vectorised processing
    withinAnyParticle = vecnorm(xyz-particles{1}.position) < particles{1}.radius;
    particles_dim = size(particles);
    for index = 2:particles_dim(2)
        withinAnyParticle = withinAnyParticle | (vecnorm(xyz-particles{index}.position) < particles{index}.radius);
    end
    output = withinAnyParticle;
end