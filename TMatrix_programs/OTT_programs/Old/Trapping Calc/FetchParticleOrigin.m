%{
. Returns the position a particle should be placed at initially when
created, given you have N particles that need to be arranged in a given
shape, and this is the ith particle forming the shape

. layout_type = STRING: The name for the type of arrangement you want
. scale_length = FLOAT: Length that is relevant to the layout wanted (e.g. a width, radius, etc)
. particle_total = INT: Total number of particles in this arrangement
. particle_index = INT: Which particle number this is in your arrangement
%}
function output = FetchParticleOrigin(layout_type, scale_length, particle_total, particle_index, offset)
    if(layout_type == "line")
        %Spreads particle along the x-axis, from 0 to width
        width = scale_length;
        spacing = width/particle_total;
        output = [(particle_index-1)*spacing; 0; 0];    %Note;-1 as indices start at 1 in MATLAB
    elseif(layout_type == "circle")
        %Spread evenly along the circumference of a circle in the Z=0 plane
        theta = 2*pi/particle_total;
        output = [scale_length*cos(theta*(particle_index +offset));scale_length*sin(theta*(particle_index +offset));0];
    else
        %Default to coord origin if invalid layout_type given
        output = [0;0;0];
    end
end