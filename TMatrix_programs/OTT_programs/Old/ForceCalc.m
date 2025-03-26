function Ftot = ForceCalc(shape_position,radius, field, CubeN)


mu0 = 4e-7*pi;
e0 = 8.854e-12;

normals = [[1,0,0];[-1,0,0]; [0,1,0]; [0,-1,0]; [0,0,1]; [0,0,-1] ] ;
lin = linspace(-radius, radius, CubeN);
Ftot = [0;0;0];

% Loop over all faces
for cube_face = 1:6

    % Get the normal
    normal = normals(cube_face,:);
    
    
    % loop over some square
    for q1 = lin
        for q2 = lin
            % x face so y,z = q1,q2
            if normal(1) ~= 0
                if normal(1) > 0
                    x = radius;
                else
                    x = -radius;
                end
                y = q1;
                z = q2;
            
            % y face
            elseif normal(2) ~= 0
                if normal(2) > 0
                    y = radius;
                else
                    y = -radius;
                end
                x = q1;
                z = q2;
            % z face
            else
                if normal(3) > 0
                    z = radius;
                else
                    z = -radius;
                end
                x = q1;
                y = q2;
            end
            % shift coords to be around the sphere.
            x = x + shape_position(1);
            y = y + shape_position(2);
            z = z + shape_position(3);

            % Get EM field and then the tensor.
            [E, H] = field.emFieldXyz([x;y;z]);
            

            TM = 1/2 * real(e0*kron(E, conj(E).') + mu0*kron(H, conj(H).') ...
                - 1/2 * (e0 * dot(E,conj(E)) + mu0 * dot(H,conj(H))) * eye(3) );
            
            % TM multiplied by normal and element area.
            addition = TM * normal' * (2*radius/CubeN)^2;
            Ftot = Ftot + addition;
            
        end
    end
    
end