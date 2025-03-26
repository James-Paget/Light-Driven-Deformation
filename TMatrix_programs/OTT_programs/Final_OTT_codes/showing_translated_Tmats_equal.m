wavelength0 = 1064e-9;  % Vacuum wavelength
nrel = 1.59/1.33; % polystyrene / water
radius = 0.5e-6;
pos1 = [-2*radius; 0; 0];
pos2 = [2*radius; 0; 0];

% make shapes and translate them to different positions
shape1 = ott.shapes.Shape.simple("sphere", radius);
shape1.position = pos1;
shape2 = ott.shapes.Shape.simple("sphere", radius);
shape2.position = pos2;

% Calculate T-matrix of each sphere
T1 = ott.TmatrixMie.simple(ott.shapes.Union(shape1), 'index_relative', nrel, ...
    'index_medium', 1.0, 'wavelength0', wavelength0);
T2 = ott.TmatrixMie.simple(ott.shapes.Union(shape2), 'index_relative', nrel, ...
    'index_medium', 1.0, 'wavelength0', wavelength0);

T1values = full(T1.data);
T2values = full(T2.data);

% Loop over matrix elements and print "not equal" if the element are different
fprintf("Test start\n")
for i = 1:length(T1.data)
    for j = 1:length(T1.data)
        val1 = T1values(i,j);
        val2 = T2values(i,j);
        if real(val1) ~= real(val2) || imag(val1) ~= imag(val2)
            fprintf("Not equal: %f + %fi and %f + %fi\n",real(val1),imag(val1),real(val2),imag(val2))
        end
    end
end
fprintf("Test end\n")



