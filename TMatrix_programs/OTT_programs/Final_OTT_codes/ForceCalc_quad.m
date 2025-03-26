function resultant_force = ForceCalc_quad(shape_position, shape_radius, field)
    % calculate the force in a cube at shape_position with side length
    % 2*shape_radius due to the field.

    resultant_force = [0;0;0];
    normals = { [1;0;0],[-1;0;0], [0;1;0],[0;-1;0], [0;0;1],[0;0;-1] };
    
    % loop over each side of the cube
    for normal_index = 1:6
        normal = normals{normal_index};
        
        resultVec = [0;0;0];
        % quad2d is for scalar funcs so calculate each component separately
        for component = 1:3
            % f is the function forceTM but with the parameters put in.
            f = setParams(shape_radius, shape_position, normal, field, component);
            [result, ~] = quad2d(f,-shape_radius, shape_radius, -shape_radius, shape_radius);
            resultVec(component) = result;
        end
        
        % sum the contributions from each face.
        resultant_force = resultant_force + resultVec;

        % Print force on each face:
        % disp("Force on face with normal ("+ normal(1)+", "+ normal(2)+", "+ normal(3)+")")
        % disp(resultVec)
    end 
end




function f = setParams(shape_radius, shape_position, normal, field, component)
    % return forceTM with the parameters already set so f = f(Q1, Q2)
    f = @(Q1, Q2) forceTM(Q1, Q2, shape_radius, shape_position, normal, field, component);

end


function force = forceTM(Q1, Q2, shape_radius, shape_position, normal, field, component)
   
    shape = size(Q1);
    force = zeros(shape);

    for i = 1:shape(1)
        for j = 1:shape(2)
            q1 = Q1(i,j);
            q2 = Q2(i,j);
    
            mu0 = 4e-7*pi;
            e0 = 8.854e-12;
        
            % Work out which coords form the plane based on the given normal.
            % x face so y,z = q1,q2
            if normal(1) ~= 0
                if normal(1) > 0
                    x = shape_radius;
                else
                    x = -shape_radius;
                end
                y = q1;
                z = q2;
            
            % y face
            elseif normal(2) ~= 0
                if normal(2) > 0
                    y = shape_radius;
                else
                    y = -shape_radius;
                end
                x = q1;
                z = q2;
        
            % z face
            else
                if normal(3) > 0
                    z = shape_radius;
                else
                    z = -shape_radius;
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
            
            %{
            % Cutoff the field if the position is close to the origin.
            % Scale it to the magnitude of the field at the cutoff radius
            cutoff_rad = 0.6e-6;
            [E2, H2] = field.emFieldXyz([cutoff_rad;0;0]);
            rad = sqrt(x^2+y^2);
            if rad < cutoff_rad
                E = E./abs(E).*abs(E2);
                H = H./abs(H).*abs(H2);
            end
            %}
            
            TM = 1/2 * real(e0*kron(E, conj(E).') + mu0*kron(H, conj(H).') ...
                - 1/2 * (e0 * dot(E,conj(E)) + mu0 * dot(H,conj(H))) * eye(3) );
            
            % TM multiplied by normal.
            force(i,j) = dot(TM(component,:), normal);

        end
    end
end

