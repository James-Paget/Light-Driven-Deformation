%{
. Considers a resultant vector force and calculaets the enw position and
velocity a particle would have from this force applied to it in an eulerian
timestep jump
. This is inaccurate but can be a useful test to probe the action of
different forces
%}
function [pos_out, vel_out] = StepByForce_Euler(particle, original_velocity, force_net, timestep)
    acc = force_net/1.0;    %Assumed particle mass of 1kg -> Just a scale factor to the animation, true mass is not particularly important to resultant motion
    vel_out = original_velocity +acc*timestep;
    pos_out = particle.position +vel_out*timestep;
end