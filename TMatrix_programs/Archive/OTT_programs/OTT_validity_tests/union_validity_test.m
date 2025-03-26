%{
. Performs tests to see if the ott.shapes.Union() function, when used with
a T-matrix, allows for correct inter-particle scattering
%}
disp("Program start");

close all;
ott.warning('once');
ott.change_warnings('off');

n_medium = 1.3;
n_particle = 1.6;
wavelength0 = 1064e-9;
wavelength_medium = wavelength0 / n_medium;

Test_1(n_medium, n_particle, wavelength0, wavelength_medium);

disp("Program end");


%% FUNCTION TO MANUALLY CALCULATE UNION
function output = Test_1(n_medium, n_particle, wavelength0, wavelength_medium)
    %{
    Want to find the;
    - Incident beam
    - Scattered beam from each particle
    - Scattered beam from cross-particle terms
    - Sum these for a total field
    - Compare this to the union result for total field
    - If they agree, then union works with multiple particles
    %}
    %Setup parameters
    view_range = 2e-6;
    radius = 0.5*wavelength_medium;
    pos_p1 = [-2.0*radius; 0; 0];
    pos_p2 = [+2.0*radius; 0; 0];
    %Specify particles
    shape_p1 = ott.shapes.Sphere(radius, pos_p1);
    shape_p2 = ott.shapes.Sphere(radius, pos_p2);
    %Specifies beam
    beam_inc = ott.BscPlane(0, 0, ...
        'polarisation', [ 1 0 ], ...
        'index_medium', n_medium, ...
        'wavelength0', wavelength0);
    beam_inc_p1 = translateXyz(beam_inc, pos_p1);   %Shifted as T-matrix centres the particle, so beam must be the one offset
    beam_inc_p2 = translateXyz(beam_inc, pos_p2);   %"" ""

    %Calculating union result -> Note; requires TmatrixMie, does NOT work
    %for Tmatrix standard
    shape_union = ott.shapes.Union([shape_p1, shape_p2]);
    T_union = ott.TmatrixMie.simple(shape_union, ...
       'wavelength0', wavelength0, ...
       'index_medium', n_medium, ...
       'index_particle', n_particle);
    beam_scat_union = T_union * beam_inc;
    beam_total_union = beam_inc +beam_scat_union;

    %Calculates T-Matrix for each particle
    T_p1 = ott.Tmatrix.simple(shape_p1, ...
       'wavelength0', wavelength0, ...
       'index_medium', n_medium, ...
       'index_particle', n_particle);
    T_p2 = ott.Tmatrix.simple(shape_p2, ...
       'wavelength0', wavelength0, ...
       'index_medium', n_medium, ...
       'index_particle', n_particle);;
    %Calculate scattering for each particle individually
    beam_scat_p1 = T_p1 * beam_inc_p1;
    beam_scat_p2 = T_p2 * beam_inc_p2;
    %Calculate cross-particle terms
    beam_scat_p1Fromp2 = T_p1 * beam_scat_p2;
    beam_scat_p2Fromp1 = T_p2 * beam_scat_p1;
    %Translate beams back to correct positions (particles at their positions, not origin)
    %% PASS

    %Plot each component
    figure();
    
    subplot(2,3,1);
    %% SHIFT MASK BY POS
    beam_inc_p1.visualise('axis', 'z', ...
        'mask', @(xyz) vecnorm(xyz) < radius, ...
        'range', [1,1]*view_range);
    title("P1 incident beam");
    subplot(2,3,2);
    %% SHIFT MASK BY POS
    beam_inc_p2.visualise('axis', 'z', ...
        'mask', @(xyz) vecnorm(xyz) < radius, ...
        'range', [1,1]*view_range);
    title("P2 incident beam");
    subplot(2,3,3);
    %% SHIFT MASK BY POS
    beam_scat_p1.visualise('axis', 'z', ...
        'mask', @(xyz) vecnorm(xyz) < radius, ...
        'range', [1,1]*view_range);
    title("P1 scattered beam");
    subplot(2,3,4);
    %% SHIFT MASK BY POS
    beam_scat_p2.visualise('axis', 'z', ...
        'mask', @(xyz) vecnorm(xyz) < radius, ...
        'range', [1,1]*view_range);
    title("P2 scattered beam");
    subplot(2,3,5);
    %% SHIFT MASK BY POS
    beam_scat_p1Fromp2.visualise('axis', 'z', ...
        'mask', @(xyz) vecnorm(xyz) < radius, ...
        'range', [1,1]*view_range);
    title("P1 From P2 scattered beam");
    subplot(2,3,6);
    %% SHIFT MASK BY POS
    beam_scat_p2Fromp1.visualise('axis', 'z', ...
        'mask', @(xyz) vecnorm(xyz) < radius, ...
        'range', [1,1]*view_range);
    title("P2 From P1 scattered beam");

    figure();
    subplot(1,2,1);
    beam_total_p1p2 = beam_inc + (beam_scat_p1 +beam_scat_p1Fromp2) + (beam_scat_p2 +beam_scat_p2Fromp1);
    beam_total_p1p2.visualise('axis', 'z', ...
        'mask', @(xyz) vecnorm(xyz) < radius, ...
        'range', [1,1]*view_range);
    title("Total field from MANUAL calc");
    subplot(1,2,2);
    beam_total_union.visualise('axis', 'z', ...
        'mask', @(xyz) vecnorm(xyz) < radius, ...
        'range', [1,1]*view_range);
    title("Total field from UNION calc");


    %%
    %% ATTEMPT TO SHIFT HERE
    %%

    %Get beam E values
    sample_number = 80;                                             %Number of samples to get for final SHIFTED grid
    shift = pos_p1;                                                 %Viewed in Z=0 plane => only allowed shifts in XY (EM only calculated for XY grid)
    distance_per_sample = 2.0*view_range / sample_number;
    shift_samples = floor(shift / distance_per_sample);
    max_shift = abs( shift(1) );    %%REPLACE WITH MAX
    additional_samples = abs(shift_samples(1));%floor((max_shift) / distance_per_sample);
    E_mag_sub   = zeros(sample_number, sample_number);              %The smaller set, of values to hold shifted
    E_mag_super = FetchBeam_EMField(beam_scat_p1, view_range +max_shift, sample_number +2.0*additional_samples); %The larger set, of values to pull to allow shifting

    %Shift beam E values
    %(1) Get larger set of E values in grid + empty grid to hold shifted
    %values
    %(2) Move over correct values to shifted
    for j =1:sample_number
        for i =1:sample_number
            %Go through each space in the shifted grid, pull correct value
            %from larger
            shift_j = additional_samples +j-shift_samples(2);
            shift_i = additional_samples +i-shift_samples(1);
            E_mag_sub(j,i) = E_mag_super(shift_j, shift_i);
        end
    end

    E_mag_super_extracted = E_mag_super( (additional_samples+1:sample_number+additional_samples) , (additional_samples+1:sample_number+additional_samples) );

    assignin("base", "E_sub", E_mag_sub);
    assignin("base", "E_super", E_mag_super);
    assignin("base", "E_super_extracted", E_mag_super_extracted);

    %Plot values (over exposed by large value at (0,0))
    figure();

    xrange = linspace(-1, 1, sample_number)*view_range;
    yrange = linspace(-1, 1, sample_number)*view_range;
    [xx, yy] = meshgrid(xrange, yrange);

    E_mag_sub = log(E_mag_sub);
    E_mag_super_extracted = log(E_mag_super_extracted);

    subplot(1,2,1);
    surf(xx, yy, E_mag_sub);
    view(0, 90);
    title("Shifted SUB plot");

    subplot(1,2,2);
    surf(xx, yy, E_mag_super_extracted);
    view(0, 90);
    title("Original SUPER plot");

    output = 0;
end

function output= FetchBeam_EMField(beam, range, sample_number)
    xrange = linspace(-1, 1, sample_number)*range;
    yrange = linspace(-1, 1, sample_number)*range;
    [xx, yy] = meshgrid(xrange, yrange);
    xyz = [xx(:) yy(:) zeros(size(xx(:)))].';
    [E, H] = beam.emFieldXyz(xyz);
    E_mag=reshape(sqrt(sum(real(E).^2,1)),[sample_number, sample_number]);
    output = E_mag;
end