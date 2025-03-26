function resultant_force = ForceCalc_simpson(shape_position, shape_radius, field, samples)
    % round up to be even
    if mod(samples,2) == 1
        samples = samples + 1;
    end

    
    singleSamples = ForceCalc_modified(shape_position, shape_radius, field, samples/2);
    % double has 1/4 the error.
    doubleSamples = ForceCalc_modified(shape_position, shape_radius, field, samples);

    resultant_force = (doubleSamples*4 - singleSamples)/3;
    
end