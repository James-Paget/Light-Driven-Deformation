function resultant_force = ForceCalc_modified(shape_position, shape_radius, field, samples)
    mu0 = 4e-7*pi;
    e0 = 8.854e-12;
    
    sample_width = (2.0*shape_radius)/samples;  %Width of sample segment in real units
    sample_area = (sample_width)^2;            %Small sample surface element area
    
    resultant_force = [0;0;0];
    normals = { [1;0;0],[-1;0;0], [0;1;0],[0;-1;0], [0;0;1],[0;0;-1] };
    for normal_index = 1:6
        normal = normals{normal_index};
        for j = 1:samples
            for i = 1:samples
                %Get coordinates of this point
                xyz = FetchBoundingCubeXyz(shape_position, shape_radius, sample_width, normal, j, i);
                %Calculate the stress matrix at each point
                [E, H] = field.emFieldXyz(xyz);
        
                %kron(A, B^T) => dyadic product for column vectors A,B
                TM = 1/2 * real(e0*kron(E, conj(E).') + mu0*kron(H, conj(H).') ...
                    - 1/2 * (e0 * dot(E,conj(E)) + mu0 * dot(H,conj(H))) * eye(3) );
        
                %Get integrand here (dot with normal)
                integrand = (TM * normal)*sample_area;
                resultant_force = resultant_force + integrand;
            end
        end
    end


    %{
    %% XY PLANES
    %Loop over surface of a face of bounding cube (TOP FACE)
    normal = [0;0;1];   %Rows not columns for correct matrix operations
    for j = 1:samples
        for i = 1:samples
            %Get coordinates of this point
            xyz = [shape_position(1)-shape_radius/2.0+i*sample_width; shape_position(2)-shape_radius/2.0+j*sample_width; shape_position(3)+shape_radius];
            %Calculate the stress matrix at each point
            [E, H] = field.emFieldXyz(xyz);
    
            %kron(A, B^T) => dyadic product for column vectors A,B
            TM = 1/2 * real(e0*kron(E, conj(E).') + mu0*kron(H, conj(H).') ...
                - 1/2 * (e0 * dot(E,conj(E)) + mu0 * dot(H,conj(H))) * eye(3) );
    
            %Get integrand here (dot with normal)
            integrand = (TM * normal)*sample_area;
            resultant_force = resultant_force + integrand;
        end
    end
    %Loop over surface of a face of bounding cube (BOTTOM FACE)
    normal = [0;0;-1];   %Rows not columns for correct matrix operations
    for j = 1:samples
        for i = 1:samples
            %Get coordinates of this point
            xyz = [shape_position(1)-shape_radius/2.0+i*sample_width; shape_position(2)-shape_radius/2.0+j*sample_width; shape_position(3)-shape_radius];
            %Calculate the stress matrix at each point
            [E, H] = field.emFieldXyz(xyz);
            TM = 1/2 * real(e0*kron(E, conj(E).') + mu0*kron(H, conj(H).') ...
                - 1/2 * (e0 * dot(E,conj(E)) + mu0 * dot(H,conj(H))) * eye(3) );
            %Get integrand here (dot with normal)
            integrand = (TM * normal)*sample_area;
            resultant_force = resultant_force + integrand;
        end
    end
    
    
    %% YZ PLANES
    %Loop over surface of a face of bounding cube (FRONT FACE)
    normal = [1;0;0];   %Rows not columns for correct matrix operations
    for j = 1:samples
        for k = 1:samples
            %Get coordinates of this point
            xyz = [shape_position(1)+shape_radius; shape_position(2)-shape_radius/2.0+j*sample_width; shape_position(3)-shape_radius/2.0+k*sample_width];
            %Calculate the stress matrix at each point
            [E, H] = field.emFieldXyz(xyz);
            TM = 1/2 * real(e0*kron(E, conj(E).') + mu0*kron(H, conj(H).') ...
                - 1/2 * (e0 * dot(E,conj(E)) + mu0 * dot(H,conj(H))) * eye(3) );
            %Get integrand here (dot with normal)
            integrand = (TM * normal)*sample_area;
            resultant_force = resultant_force + integrand;
        end
    end
    %Loop over surface of a face of bounding cube (BACK FACE)
    normal = [-1;0;0];   %Rows not columns for correct matrix operations
    for j = 1:samples
        for k = 1:samples
            %Get coordinates of this point
            xyz = [shape_position(1)-shape_radius; shape_position(2)-shape_radius/2.0+j*sample_width; shape_position(3)-shape_radius/2.0+k*sample_width];
            %Calculate the stress matrix at each point
            [E, H] = field.emFieldXyz(xyz);
            TM = 1/2 * real(e0*kron(E, conj(E).') + mu0*kron(H, conj(H).') ...
                - 1/2 * (e0 * dot(E,conj(E)) + mu0 * dot(H,conj(H))) * eye(3) );
            %Get integrand here (dot with normal)
            integrand = (TM * normal)*sample_area;
            resultant_force = resultant_force + integrand;
        end
    end
    
    
    %% XZ PLANES
    %Loop over surface of a face of bounding cube (RIGHT FACE)
    normal = [0;1;0];   %Rows not columns for correct matrix operations
    for i = 1:samples
        for k = 1:samples
            %Get coordinates of this point
            xyz = [shape_position(1)-shape_radius/2.0+i*sample_width; shape_position(2)+shape_radius; shape_position(3)-shape_radius/2.0+k*sample_width];
            %Calculate the stress matrix at each point
            [E, H] = field.emFieldXyz(xyz);
            TM = 1/2 * real(e0*kron(E, conj(E).') + mu0*kron(H, conj(H).') ...
                - 1/2 * (e0 * dot(E,conj(E)) + mu0 * dot(H,conj(H))) * eye(3) );
            %Get integrand here (dot with normal)
            integrand = (TM * normal)*sample_area;
            resultant_force = resultant_force + integrand;
        end
    end
    %Loop over surface of a face of bounding cube (BACK FACE)
    normal = [0;-1;0];   %Rows not columns for correct matrix operations
    for i = 1:samples
        for k = 1:samples
            %Get coordinates of this point
            xyz = [shape_position(1)-shape_radius/2.0+i*sample_width; shape_position(2)-shape_radius; shape_position(3)-shape_radius/2.0+k*sample_width];
            %Calculate the stress matrix at each point
            [E, H] = field.emFieldXyz(xyz);
            TM = 1/2 * real(e0*kron(E, conj(E).') + mu0*kron(H, conj(H).') ...
                - 1/2 * (e0 * dot(E,conj(E)) + mu0 * dot(H,conj(H))) * eye(3) );
            %Get integrand here (dot with normal)
            integrand = (TM * normal)*sample_area;
            resultant_force = resultant_force + integrand;
        end
    end
    %}
end

function xyz = FetchBoundingCubeXyz(shape_position, shape_radius, sample_width, normal, loop_a, loop_b)
    %{
    . loop_ and loop_b are the two indices used (2 from either i,j or k)
    . loop_a & b always ordered ascending (i->j->k)
    %}
    xyz = [NaN; NaN; NaN];
    if(abs(normal(1)) == 1)
        xyz(1) = shape_position(1)+normal(1)*shape_radius;
    elseif(abs(normal(2)) == 1)
        xyz(2) = shape_position(2)+normal(2)*shape_radius;
    elseif(abs(normal(3)) == 1)
        xyz(3) = shape_position(3)+normal(3)*shape_radius;
    end
    loop_a_found = false;
    for index = 1:3
        if( isnan(xyz(index)) )
            if(~loop_a_found)
                xyz(index) = shape_position(index)-shape_radius/2.0+loop_a*sample_width;
                loop_a_found = true;
            else
                xyz(index) = shape_position(index)-shape_radius/2.0+loop_b*sample_width;
            end
        end
    end
end