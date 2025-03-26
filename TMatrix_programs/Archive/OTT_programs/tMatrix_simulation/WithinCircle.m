%{
. Checks if a 'point' is located within a circle located at position 'origin'
with radius 'radius'

**
NOT USED CURRENTLY
**
%}
function output = WithinCircle(point, origin, radius)
    separation = vecnorm(point-origin);
    output = separation <= radius;
end