import subprocess
import numpy as np
import random
import matplotlib.pyplot as plt
import math
import cmath
import scipy.integrate as integrate
import scipy.special as sp
from matplotlib import colormaps as cm
import matplotlib as mpl
from functools import partial
import sys
from make_animation_and_quiver import make_quiver_plot

"""
(1) Need to test T-avg force calc works like ADDA's
    -> Bring ADDA plane wave into here
    -> Make visualisation functions easier to use
    -> Compare forces here to ADDA
(2) Get Laguerre working fully
(3) Try trap particle in ADDA
(4) Project would be novel if trapped well, compared forces for assessement on close particle validty, and had extended ADDA software
"""

def PullAddaData_Parameters(filename):
    """
    . Pulls other miscellaneous parameters from the log file adda generates
    . Values are;
    [lambda, [refractive_index_real, refractive_index_imag], dipole_size, [beam_center_x, beam_center_y, beam_center_z], [beam_dir_x, beam_dir_y, beam_dir_z]]
    ####
    ## BEAM CENTERS NOT DONE YET BUT SAME SITUATION
    ####
    """
    #Load file
    file = open(filename, "r");
    #Check each line for wanted parameters
    parameter_set = {};
    for line in file:
        startIndex = 0;         #Used to separate substrings of values
        for char_index in range(0, len(line)):
            character = line[char_index:char_index+1];
            word = line[startIndex:char_index];
            #print("char= ",character);
            #print("word= ",word);
            match word:
                case "lambda: ":
                    value = -1.0;
                    try:
                        value = float(line[char_index:]); #Note; Remainder of line is all part of the float parameter, with no spaces => this approach works
                    except:
                        print("Could NOT cast parameter to float: lambda");
                    parameter_set.update({"lambda":value});
                    startIndex = char_index+1;
                case "refractive index: ":
                    value = [];
                    #Further refine parameter
                    unref_value = line[char_index:];
                    startIndex = 0      #Temporarily use startIndex for this unrelated task
                    for unref_char_index in range(0,len(unref_value)):
                        unref_char = unref_value[unref_char_index:unref_char_index+1];
                        match unref_char:
                            ####
                            ## COULD JUST DO AN IF AND COMBINE THE +- CASES --> MAYBE NICER TO SEE
                            ####
                            case "+":
                                #Implies the real part has been obtained
                                value.append(unref_value[startIndex:unref_char_index]); #Casted later
                                startIndex = unref_char_index+1;
                            case "-":
                                # "" ""
                                #Note; ONLY one of the + or - should trigger
                                value.append(unref_value[startIndex:unref_char_index]); #Casted later
                                startIndex = unref_char_index+1;
                            case "i":
                                #Implies the imaginary part has been obtained
                                value.append(unref_value[startIndex:unref_char_index]); #Casted later
                                startIndex = unref_char_index+1;
                    #Add this parameter
                    try:
                        value[0] = float(value[0]);
                        value[1] = float(value[1]);
                    except:
                        pass;
                        #print("Could NOT cast parameter to float: n");
                    parameter_set.update({"n":value});
                    startIndex = char_index+1;
                case "Dipole size: ":
                    value = -1.0;
                    #Further refine parameter
                    unref_value = line[char_index:];
                    startIndex = 0;
                    for unref_char_index in range(0, len(unref_value)):
                        unref_char = unref_value[unref_char_index:unref_char_index+1];
                        if(unref_char == " "):
                            #Implies end of value has been seen
                            value = unref_value[startIndex:unref_char_index];#Casted later
                            startIndex = unref_char_index+1;
                            break;
                    #Add this parameter
                    try:
                        value = float(value); #Note; Remainder of line is all part of the float parameter, with no spaces => this approach works
                    except:
                        print("Could NOT cast parameter to float: dipole_size");
                    parameter_set.update({"dipole_size":value});
                    startIndex = char_index+1;
                case "NEW PARAMETER CASE":
                    pass
    return parameter_set;

def PullAddaData_DipolePositions(filename):
    """
    ##
    ## NOTE; Somewhat obsolete as forces already give lattice positions in RAW as well
    ##
    . Performs this action for the dipoles of 1 particle
    . Pulls dipole positions from the .geom file generated by adda formatted as;
    . filename = name of file storing dipole position data, with the correct extension

    [ dip_pos, dip_pos, ... ]
    where dip_Pos = [x,y,z]
    
    .** Note; Assumes the file is in the same folder as this python script
    ####
    ## MULTIPLE EACH BY THE SPACING OF THE LATTICE, SO IT IS A RAW POSITION [NON-SCALED] --> DO THIS AT THE END ##
    ####
    """
    #Load file
    file = open(filename, "r");
    #Ignore description lines (first 4)
    for i in range(0,4):
        file.readline();
    #Begin reading and storing following dipole positions
    dipole_positions = [];
    for line in file:
        dipole_position = [];
        startIndex = 0;         #Used to separate substrings of values
        for char_index in range(0, len(line)):
            character = line[char_index:char_index+1];
            #print("char= ",character);
            if(character == " "):
                #Moving onto next coord => store current coord and reset startIndex
                value = -1;
                try:
                    value = int(line[startIndex:char_index]);
                except:
                    print(" -Could NOT cast dipole coordinate to int");
                dipole_position.append(value);   #Note; Disclude the space
                startIndex = char_index+1;
            elif(char_index == len(line)-1):
                #Handle last value separatly as it does NOT have a space after it
                value = -1;
                try:
                    value = int(line[startIndex:]);
                except:
                    print(" -Could NOT cast dipole coordinate to int");
                dipole_position.append(value);   #Note; Disclude the space
        dipole_positions.append(dipole_position);
    return dipole_positions;

##
## SHOULD JUST RENAME THIS A GENERALISED DATA PULL FOR [a,b,c,...] format, space sep.
##    
def PullAddaData_DipolePosForces(filename):
    """
    . Performs this action for the dipoles of 1 particle
    . Pulls the forces acting at each dipole lattice position generated by the 
    adda "RadForce-Y" file
    .filename = name of file storing dipole position data, with the correct extension

    Will format the output as such;
    [dip, force, dip_force, ...]
    where dip_force = [x,y,x, Fx,Fy,Fz] for its position and direction respecitively

    .** Note; Assumes the file is in the same folder as this python script
    .** Note; These are in RAW values (not relative)
    """
    #Load file
    file = open(filename, "r");
    #Ignore description lines (first 1)
    for i in range(0,1):
        file.readline();
    #Begin reading and storing following dipole positions
    dipole_posForce_set = [];
    for line in file:
        dipole_posForce = [];
        termIndex  = 0;     #Used to count which term is about to be added (some terms are discluded as unnecessary)
        startIndex = 0;     #Used to separate substrings of values
        for char_index in range(0, len(line)):
            character = line[char_index:char_index+1];
            #print("char= ",character);
            if(character == " "):
                #if(termIndex != 3): #add 0th, 1st, 2nd BUT NOT 3rd (magnitude, ...)
                #    pass
                #Moving onto next value => store current coord and reset startIndex
                value = -1.0;
                try:
                    value = float(line[startIndex:char_index]);
                except:
                    print(" -Could NOT cast dipole coordinate to int");
                dipole_posForce.append(value);   #Note; Disclude the space
                termIndex+=1;
                startIndex = char_index+1;
            elif(char_index == len(line)-1):
                #Handle last value separatly as it does NOT have a space after it
                value = -1.0;
                try:
                    value = float(line[startIndex:]);
                except:
                    print(" -Could NOT cast dipole coordinate to int");
                dipole_posForce.append(value);   #Note; Disclude the space
        dipole_posForce_set.append(dipole_posForce);
    return dipole_posForce_set;

def PullAddaData_DipolePosPolarisations(filename):
    """
    . Performs this action for the dipoles of 1 particle
    . Pulls the polarisations acting at each dipole lattice position generated by the 
    adda "DipPol-Y/X" file --> Only 1 required for systems with symmetry

    .filename = name of file storing dipole pos+polarisation data, with the correct extension

    Will format the output as such;
    [x,y,x, Px_re,Px_im,Py_re,Py_im,Pz_re,Pz_im] for its position and direction respectively

    .** Note; Assumes the file is in the same folder as this python script
    .** Note; These are in RAW values (not relative)
    """
    #Load file
    file = open(filename, "r");
    #Ignore description lines (first 1)
    for i in range(0,1):
        file.readline();
    #Begin reading and storing following dipole positions
    dipole_posPolarisation_set = [];
    for line in file:
        dipole_posPolarisation = [];
        termIndex  = 0;     #Used to count which term is about to be added (some terms are discluded as unnecessary)
        startIndex = 0;     #Used to separate substrings of values
        for char_index in range(0, len(line)):
            character = line[char_index:char_index+1];
            #print("char= ",character);
            if(character == " "):
                #if(termIndex != 3): #add 0th, 1st, 2nd BUT NOT 3rd (magnitude, ...)
                #    pass
                #Moving onto next value => store current coord and reset startIndex
                value = -1.0;
                try:
                    value = float(line[startIndex:char_index]);
                except:
                    print(" -Could NOT cast dipole coordinate to int");
                dipole_posPolarisation.append(value);   #Note; Disclude the space
                termIndex+=1;
                startIndex = char_index+1;
            elif(char_index == len(line)-1):
                #Handle last value separatly as it does NOT have a space after it
                value = -1.0;
                try:
                    value = float(line[startIndex:]);
                except:
                    print(" -Could NOT cast dipole coordinate to int");
                dipole_posPolarisation.append(value);   #Note; Disclude the space
        dipole_posPolarisation_set.append(dipole_posPolarisation);
    return dipole_posPolarisation_set;

def PullAddaData_DipolePosEInternal(filename):
    """
    . Performs this action for the dipoles of 1 particle
    . Pulls the E internal field from ADDA generated file

    .filename = name of file storing dipole pos+polarisation data, with the correct extension

    Will format the output as such;
    [x,y,x, Ex_re,Ex_im,Ey_re,Ey_im,Ez_re,Ez_im] for its position and direction respecitively

    .** Note; These are in RAW values (not relative, in micrometers and possibly CGS too)
    """
    #Load file
    file = open(filename, "r");
    #Ignore description lines (first 1)
    for i in range(0,1):
        file.readline();
    #Begin reading and storing following dipole positions
    dipole_posEInternal_set = [];
    for line in file:
        dipole_posEInternal = [];
        termIndex  = 0;     #Used to count which term is about to be added (some terms are discluded as unnecessary)
        startIndex = 0;     #Used to separate substrings of values
        for char_index in range(0, len(line)):
            character = line[char_index:char_index+1];
            #print("char= ",character);
            if(character == " "):
                #if(termIndex != 3): #add 0th, 1st, 2nd BUT NOT 3rd (magnitude, ...)
                #    pass
                #Moving onto next value => store current coord and reset startIndex
                value = -1.0;
                try:
                    value = float(line[startIndex:char_index]);
                except:
                    print(" -Could NOT cast dipole coordinate to int");
                dipole_posEInternal.append(value);   #Note; Disclude the space
                termIndex+=1;
                startIndex = char_index+1;
            elif(char_index == len(line)-1):
                #Handle last value separatly as it does NOT have a space after it
                value = -1.0;
                try:
                    value = float(line[startIndex:]);
                except:
                    print(" -Could NOT cast dipole coordinate to int");
                dipole_posEInternal.append(value);   #Note; Disclude the space
        dipole_posEInternal_set.append(dipole_posEInternal);
    return dipole_posEInternal_set;

def FetchResultantForce(posForce_set):
    """
    . Considers a particles, made of N dipoles, and sums the forces to get a total 
    force to act on the rigid body particle

    posForce_set = a list of the posForce objects, which are lists [x,y,z, |F|^2,Fx,Fy,Fz] for 
        the position and direction the force acts at
    """
    resultantForce = [0.0, 0.0, 0.0];
    for posForce in posForce_set:
        for coord in range(0,len(resultantForce)):
            resultantForce[coord] += posForce[coord+4];     #+3 to skip the x,y,z,|F|^2 components
    return resultantForce;

def RunAdda(input_shape, input_beam, dpl, n_real, n_imag, wavelength, output_folder):
    refractive_indices = str(n_real)+" "+str(n_imag)+" "+str(n_real)+" "+str(n_imag);
    adda_command = "./adda -dpl " +str(dpl)+ " -lambda " + str(wavelength)+ " -shape read " + input_shape + " -m "+ refractive_indices + " -store_int_field -store_dip_pol -beam "+input_beam+" -dir " + str(output_folder);
    adda_command = adda_command.split(" ");
    print("=== ADDA Log ===");
    result = subprocess.run(adda_command, stdout=subprocess.DEVNULL);
    print("=== ADDA Log ===");

def PullSingleParticleGeom(filename):
    #Pull out info from single shape
    singleShape_info = [];
    singleFile = open(filename,"r");
    singleFile_rawInfo = singleFile.readline();
    singleShape_info_width = 0;
    startIndex = 0;
    dipole_pos = [];
    while(len(singleFile_rawInfo) > 0):
        for char_index in range(0,len(singleFile_rawInfo)):
            character = singleFile_rawInfo[char_index: char_index+1];
            match(character):                                                                   #
                case(" "):                                                                      # Create a dipole position vector from each line of the .geom file
                    dipole_pos.append( int(singleFile_rawInfo[startIndex: char_index+1]) );     #
                    startIndex = char_index+1;                                                  #
                    if(len(dipole_pos) == 1):
                        #When you first get an X component for a dipole position, check if it states that the x-width should be any bigger
                        if(dipole_pos[0]+1 > singleShape_info_width):   #+1 in both because indices start at 0 => X position at 2 => has width 3
                            singleShape_info_width = dipole_pos[0]+1;   #
            if(char_index == len(singleFile_rawInfo)-1):
                dipole_pos.append( int(singleFile_rawInfo[startIndex: char_index+1]) );
                startIndex = 0;
        singleShape_info.append(dipole_pos);        #
        dipole_pos = [];                            # Reset parmeters for next line
        singleFile_rawInfo = singleFile.readline(); #
    singleFile.close();
    return singleShape_info, singleShape_info_width;


def Calculate_FarField_Grid(posPolarisations, exclusionRadius, space_data, k, permitivity_free):
    """
    . Finds the E far-field at the positions given by
    . Skips over positions within the particle (exclusionRadius)

    . posPolarisations = [[x, y, z, |P|^2, Px_re, Px_im, Py_re, Py_im, Pz_re, Px_im], ...]
    . exclusionRadius  = Float, skips calcualtions of grid within this area (used to ignore internal field of particle as would be inaccurate anyway => waste+misleading to calculate)
    . space_data = [[X_start, X_end, X_samples], [Y...], [Z...]]
    """
    E_far_field = [];

    x_set = np.linspace(space_data[0][0], space_data[0][1], space_data[0][2]);
    y_set = np.linspace(space_data[1][0], space_data[1][1], space_data[1][2]);
    z_set = np.linspace(space_data[2][0], space_data[2][1], space_data[2][2]);      #Often just 1 value, in the Z=a plane

    #Convert all to numpy arrays
    posPolarisations = np.array(posPolarisations);
    for posPolarisation in posPolarisations:
        np.array(posPolarisation);
    
    #Pull out positions and polarisations separately -> Pos in micrometers (+possible CGS units)
    dip_positions        = posPolarisations[:,0:3];         #[ [ x, y, z], ... ]
    dip_polarisations_Re = posPolarisations[:,4::2];        #[ [Px_re, Py_re, Pz_re], ... ]
    dip_polarisations_Im = posPolarisations[:,5::2];        #[ [Px_im, Py_im, Pz_im], ... ]
    dip_polarisations = dip_polarisations_Re +1j*dip_polarisations_Im;
    #print("dip_positions= ",dip_positions);
    #print("dip_polarisaions= ",dip_polarisations);

    for k_ind in range(0,len(z_set)):
        print("k= "+str(k_ind)+"+1/"+str(len(z_set)));
        E_far_field.append( [] );
        for j_ind in range(0,len(y_set)):
            E_far_field[k_ind].append( [] );
            for i_ind in range(0,len(x_set)):
                r_point  = np.array([x_set[i_ind], y_set[j_ind], z_set[k_ind]]);    #Vector disp. to point in space
                r_point_mag  = np.sqrt( np.sum( np.power(r_point, 2) ) );           #Distance to point in space
                if(r_point_mag > exclusionRadius):                                  #Exclude points near the origin / within the shape (not far-field scattering)
                    n_hat = r_point/r_point_mag;                                    #Normalised direction to target point in space
                    
                    F_n_polarisation = np.array([0.0 +0.0j, 0.0 +0.0j, 0.0 +0.0j]);         #Polarisation sum part of scattering amplitude
                    for dip_ind in range(0, len(dip_positions)):

                        r_dipole = np.array(dip_positions[dip_ind]);                        #Vector disp. to dipole
                        
                        #Sum over dipoles, get this scattering amplitude term --> Vectorised sum
                        F_n_polarisation += dip_polarisations[dip_ind]*np.exp(-1j*k*np.dot(r_dipole, n_hat) );      #Scattering amplitude
                    F_n = np.matmul((-1j*pow(k,3))*( np.identity(3) - np.outer(n_hat, n_hat) ), (F_n_polarisation));  #A matrix acting on a vector -> gives a vector
                
                    k_scat = k;     ############################### MAYBE 'k_scat' NEEDS TO BE * REFRACTIVE INDEX ##############################
                    E_far_field_gridPoint = ( (np.exp(1j*k_scat*r_point_mag))/(-1j*k_scat*r_point_mag) )*(F_n);    #'j' complex in python
                    E_far_field[k_ind][j_ind].append(E_far_field_gridPoint);
                else:
                    #Add 0 vector to points within exclusion radius
                    E_far_field[k_ind][j_ind].append( [0.0, 0.0, 0.0] );
    return E_far_field;

def Plot_2D_Array(slice_type, format, title_label, dataset, space_data, subplot_params=None):
    """
    . Data set formatted as 3D set of positions on grid points (grid points retrieved again from space_data, NOT explicitly stored in dataset)
    where each grid point holds a vector [e.g. (Ex, Ey, Ez) ]
    . This function will plot either the norm or components of that vector, specified through the 'format' variable
    . As this is a 2D plot, a surface plot will be done for each Z slice in the dataset, the color being based on the vector at that location (its magnitude)

    . slice_type = which plane is sliced through, e.g X,Y,Z, where Z=>looking at XY plane taken for each Z
    """
    print("Plotting data...");
    x_set = np.linspace(space_data[0][0], space_data[0][1], space_data[0][2]);
    y_set = np.linspace(space_data[1][0], space_data[1][1], space_data[1][2]);
    z_set = np.linspace(space_data[2][0], space_data[2][1], space_data[2][2]);      #Often just 1 value, in the Z=a plane

    if(slice_type == "Z"):
        #Do a separate plot for each Z plane
        for k_ind in range(0,len(z_set)):
            #For each data point
            Z_plane = np.format_float_scientific(z_set[k_ind],3);
            title_label_final = title_label+", Z="+str(Z_plane);
            data = np.zeros( (len(y_set), len(x_set)) );
            for j_ind in range(0, len(y_set)):
                for i_ind in range(0, len(x_set)):
                    #Pick parameter to plot
                    mod2_vector = np.absolute( np.array(dataset[k_ind][j_ind][i_ind]) );    #List of [[re,im], [re,im], [re,im]]
                    if(format == "norm"):
                        data[j_ind, i_ind] = np.sqrt(np.sum(np.power(mod2_vector, 2)));
                    elif(format == "X"):
                        data[j_ind, i_ind] = np.sqrt(np.sum(np.power(mod2_vector[0], 2)));
                    elif(format == "Y"):
                        data[j_ind, i_ind] = np.sqrt(np.sum(np.power(mod2_vector[1], 2)));
                    elif(format == "Z"):
                        data[j_ind, i_ind] = np.sqrt(np.sum(np.power(mod2_vector[2], 2)));
                    else:
                        print("Invalid format: ",format);

            #Plot
            #print(list(mpl.colormaps))
            X, Y = np.meshgrid(x_set, y_set);
            if(subplot_params!=None):
                plt.subplot(subplot_params[0], subplot_params[1], subplot_params[2]);
            plt.pcolor(X,Y,data, cmap=mpl.colormaps['viridis']);    #viridis, #binary_r, #binary
            plt.xlabel("X (μm)");
            plt.ylabel("Y (μm)");
            plt.colorbar();
            plt.title(title_label_final);
    if(slice_type == "Y"):
        #Do a separate plot for each Z plane
        #For each data point
        Y_plane = 0;#np.format_float_scientific(z_set[k_ind],3);
        title_label_final = title_label+", Yth plane="+str(Y_plane);  #Yth plane picked out
        data = np.zeros( (len(z_set), len(x_set)) );
        for k_ind in range(0,len(z_set)):
            for j_ind in range(0, len(y_set)):
                for i_ind in range(0, len(x_set)):
                    if(j_ind == Y_plane):
                        #Pick parameter to plot
                        mod2_vector = np.absolute( np.array(dataset[k_ind][j_ind][i_ind]) );    #List of [[re,im], [re,im], [re,im]]
                        if(format == "norm"):
                            data[i_ind, k_ind] = np.sqrt(np.sum(np.power(mod2_vector, 2)));
                        elif(format == "X"):
                            data[i_ind, k_ind] = np.sqrt(np.sum(np.power(mod2_vector[0], 2)));
                        elif(format == "Y"):
                            data[i_ind, k_ind] = np.sqrt(np.sum(np.power(mod2_vector[1], 2)));
                        elif(format == "Z"):
                            data[i_ind, k_ind] = np.sqrt(np.sum(np.power(mod2_vector[2], 2)));
                        else:
                            print("Invalid format: ",format);

        #Plot
        #print(list(mpl.colormaps))
        Z, X = np.meshgrid(z_set, x_set);
        if(subplot_params!=None):
            plt.subplot(subplot_params[0], subplot_params[1], subplot_params[2]);
        plt.pcolor(Z,X,data, cmap=mpl.colormaps['viridis']);    #viridis, #binary_r, #binary
        plt.xlabel("Z (μm)");
        plt.ylabel("X (μm)");
        plt.colorbar();
        plt.title(title_label_final);

def Format_PosVariable_List(posVariables):
    """
    . Takes a set of data in the format [x,y,z, var_x_re,var_x_im, var_y_re,var_y_im, var_z_re,var_z_im], and reformats it
    as a 3D grid of vectors
    . These 3D vector spaces can then be plotted in layers using "Plot_2D_Array"
    . Grid is a regular array, not numpy

    . grid_dim = [x,y,z] length in dipoles

    Data from ADDA is formatted in layers;
    Z=0, Sweep Xs for Y=0, Sweeps for Y=1, ..., Get plane
    """
    #Fully convert to np array
    posVariables = np.array(posVariables);
    for posVariable in posVariables:
        posVariable = np.array(posVariable);
    
    #Format nicer
    dip_positions   = posVariables[:,0:3];         #[ [ x, y, z], ... ]
    dip_variable_Re = posVariables[:,4::2];        #[ [Varx_re, Vary_re, Varz_re], ... ]
    dip_variable_Im = posVariables[:,5::2];        #[ [Varx_im, Vary_im, Varz_im], ... ]
    dip_variable    = dip_variable_Re +1j*dip_variable_Im;

    return dip_positions, dip_variable;     #[[x,y,z],...] and [[E_x,E_y,E_z],...] -> Complex E, pos and E and resp. to each other

def Generate_Beam(beam_spec, space_data):
    """
    . Generates the electric on the grid queried
    . Replacates the space data structure as you traverse through it, returns the field data

    . beam_spec = [beam_type, <params>]
    """
    field_data = [];
    space_jump = [ 
        (space_data[0][1] -space_data[0][0])/(space_data[0][2]-1), 
        (space_data[1][1] -space_data[1][0])/(space_data[1][2]-1), 
        0
    ];
    #REDO FOR NICE SOL FOR XYZ
    if(space_data[2][2] > 1):  #If want to measure for a single Z
        space_jump[2] = (space_data[2][1] -space_data[2][0])/(space_data[2][2]-1);
    else:
        space_jump[2] = 0.0;
    #REDO FOR NICE SOL FOR XYZ
    

    for k_ind in range(0, space_data[2][2]):
        field_data.append([]);
        for j_ind in range(0, space_data[1][2]):
            field_data[k_ind].append([]);
            for i_ind in range(0, space_data[0][2]):
                pos = [ 
                    space_data[0][0] +i_ind*space_jump[0], 
                    space_data[1][0] +j_ind*space_jump[1], 
                    space_data[2][0] +k_ind*space_jump[2]
                ];
                E_field = [0.0, 0.0, 0.0];
                if(beam_spec[0] == "test"):
                    E_field = [pos[0], pos[1], pos[2]];
                elif(beam_spec[0] == "laguerreGaussian"):
                    #beam_spec = [beam_type, wavelength, azi, radial, alpha, beta, w_0]
                    E_field = Get_E_LaguerreGaussian_Beam(pos, beam_spec[1], beam_spec[2], beam_spec[3], beam_spec[4], beam_spec[5], beam_spec[6]);
                elif(beam_spec[0] == "planeWave"):
                    #beam_spec = [beam_type, wavelength, azi, radial, alpha, beta, w_0]
                    E_field = Get_E_PlaneWave_Beam(pos, beam_spec[1], beam_spec[2], beam_spec[3]);
                else:
                    print("Invalid beam type");
                field_data[k_ind][j_ind].append(E_field);
    return field_data;

def Get_E_PlaneWave_Beam(pos, wavelength, beam_center, polarisation):
    """
    . Returns list for electric field vector components at a position, given by the plane wave beam equation (exact format copied from ADDA for comparison)
    . Polarisation either "X" or "Y" for now
    """
    #Raw ADDA; vSubtr(DipoleCoord+j,beam_center,r1);
    #          ctemp=imExp(WaveNum*DotProd(r1,prop));
    k = 2.0*math.pi/wavelength;
    r1 = np.array(pos)-np.array(beam_center);
    ctemp=cmath.exp(1j*k*r1[2]);    #We assume here that the 'prop' direction (for wavevector) is always in Z direction for simplicity
    if(polarisation=="X"):
        return [ctemp.real, 0.0, 0.0];
    elif(polarisation=="Y"):
        return [0.0, ctemp.real, 0.0];
    else:
        print("Invalid polarisation for plane wave: "+str(polarisation));
        return [0.0, 0.0, 0.0];

def Get_E_LaguerreGaussian_Beam(pos, wavelength, radial, azimuthal, alpha, beta, w_0, showVariables=True):
    """
    . pos = [x,y,z] point to get E at
    . radial/azimuthal = 2 parameters for laguerre beam, spin and orbital angular momentum order
        - Called m & n in paper "Gaussian, Hermite-Gaussian, and Laguerre-Gaussian beams: A primer"
    """
    z   = pos[2];
    phi = math.atan2(pos[1],pos[0]);
    rho = cmath.sqrt( pow(pos[0],2) + pow(pos[1],2) );

    l   = radial;       #azimuthal;#8.0;        #Orbital angular momentum number
    p   = azimuthal;    #radial;   #0.0;        #Spin angular momentum number
    k   = 2.0*math.pi/wavelength;
    z_R = k*pow(w_0, 2)/2.0;                    #Rayleigh factor

    if(showVariables):
        print("");
        print("=====");
        print("Variables");
        print("=====");
        print("l= "+str(l));
        print("p= "+str(p));
        print("k= "+str(k));
        print("w_0= "+str(w_0));
        print("z_R= "+str(z_R));
        print("alpha= "+str(alpha));
        print("beta= "+str(beta));
        print("");
    
    def w(zeta):
        return w_0*cmath.sqrt(1+pow(zeta,2));
    
    def L_l_p(value):
        """
        . Associated legendre polynomial
        """
        legendre_matrix = sp.clpmn(l,p,value);
        return legendre_matrix[0][l][p];        #Note; [0] here refers to generalisation of array value, hence for scalar value each index [n] is identical
   
    def u_l_p():
        # From same paper
        C = (math.factorial(p))*cmath.sqrt( (2)/(math.pi*math.factorial(p)*( math.factorial( abs(l) + p ) )) );   #From pg 8 of “Gaussian Beams In Optics of Course” paper
        return ( (C)/cmath.sqrt(1 +(pow(z,2)/pow(z_R,2))) )*pow( (rho*cmath.sqrt(2))/(w(z)), l)*(L_l_p( (2.0*pow(rho,2))/(pow(w(z),2)) ))*(cmath.exp( (-pow(rho,2))/(pow(w(z),2)) ))*(cmath.exp( (-1j*k*pow(rho,2)*z)/(2.0*(pow(z,2)+pow(z_R,2))) ))*(cmath.exp(1j*l*phi))*(cmath.exp( (1j*(2.0*p +l+1))*(cmath.atan( (z)/(z_R) )) ));
    
    def E(kappa):
        return u_l_p()*cmath.exp( (-k*pow(kappa,2)*z_R)/(2.0*(pow(k,2) - pow(kappa,2))) )*pow( (pow(kappa,2))/(pow(k,2) - pow(kappa,2)), (2.0*p+l+1.0)/(2) )*cmath.sqrt( (pow(k,2))/(pow(k,2) - pow(kappa,2)) );

    def E_vector_component(component, kappa):
        #Full vector version
        #E(kappa)*cmath.exp(1j*l*phi)*cmath.exp(1j*z*cmath.sqrt(pow(k,2)-pow(kappa,2)))*( np.array([alpha, beta, 0.0])*sp.jv(l, kappa*rho) + np.array([0.0, 0.0, 1.0])*( (kappa)/(2*cmath.sqrt(pow(k,2) - pow(kappa,2))) )*( (1j*alpha - beta)*(cmath.exp(-1j*phi))*sp.jv(l-1, kappa*rho) - (1j*alpha + beta)*(cmath.exp(1j*phi))*sp.jv(l+1, kappa*rho) ) );
        #sp.jv(l, x) = Bessel function of 1st kind, order l

        if(component == 0):     #X comp
            return E(kappa)*cmath.exp(1j*l*phi)*cmath.exp(1j*z*cmath.sqrt(pow(k,2)-pow(kappa,2)))*( alpha*sp.jv(l, kappa*rho) );
        elif(component == 1):   #Y comp
            return E(kappa)*cmath.exp(1j*l*phi)*cmath.exp(1j*z*cmath.sqrt(pow(k,2)-pow(kappa,2)))*( beta*sp.jv(l, kappa*rho) );
        elif(component == 2):   #Z comp
            return E(kappa)*cmath.exp(1j*l*phi)*cmath.exp(1j*z*cmath.sqrt(pow(k,2)-pow(kappa,2)))*(  (kappa)/(2*cmath.sqrt(pow(k,2) - pow(kappa,2))) )*( (1j*alpha - beta)*(cmath.exp(-1j*phi))*sp.jv(l-1, kappa*rho) - (1j*alpha + beta)*(cmath.exp(1j*phi))*sp.jv(l+1, kappa*rho) );
        else:
            return 0.0;

    E_X_Comp = integrate.quad(lambda x: E_vector_component(0, x), 0.0, k, complex_func=True, epsrel=1.49e-12);
    E_Y_Comp = integrate.quad(lambda x: E_vector_component(1, x), 0.0, k, complex_func=True, epsrel=1.49e-12);
    E_Z_Comp = integrate.quad(lambda x: E_vector_component(2, x), 0.0, k, complex_func=True, epsrel=1.49e-12);

    #integrate.quad returns [value, error] => want to pull just the value
    return [E_X_Comp[0].real, E_Y_Comp[0].real, E_Z_Comp[0].real];

def Sweep_Beam_Generation(space_data, wavelength):
    """
    . Sweeps through some alpha and beat values for a laguerre beam generated
    . Ideally they would be set by the intensity of the real bea generated, however here this sweep is intended to match the measured intensity to that of OTT
    """
    alpha_base = 1.0;
    beta_base  = 1.0;
    w_0_base = 0.6*wavelength;

    N_iter = 3;
    sub_width  = N_iter;
    sub_height = N_iter;

    w_0 = w_0_base;
    for beta_sweep in range(0, N_iter):
        beta = beta_base*(1+beta_sweep);
        for alpha_sweep in range(0, N_iter):
            alpha = alpha_base*(1+alpha_sweep);
            field_data = Generate_Beam(["laguerreGaussian",wavelength,0,8, alpha,beta, w_0], space_data);
            Plot_2D_Array(
                "Z",
                "norm", 
                "E norm, alpha="+str(alpha)+", beta="+str(beta)+", W_0="+str(w_0), 
                field_data, space_data, 
                subplot_params=[sub_width, sub_height, 1+beta_sweep +alpha_sweep*N_iter]
            );
    plt.show();

def Calculate_T_Averaged_Dipole_PosForce(dipole_size, beam_spec):
    """
    . Calculates a vector for the set of points given
    . Force is complex for generalisation

    . Assumes the internal fields calculated by ADDA were AS A RESULT OF the beam (specified through beam_spec) given
    """
    print("Calculating time-averaged forces for "+beam_spec[0]);
    #Pull internal E fields from ADDA
    posEInternal_set = PullAddaData_DipolePosEInternal("output_data/IntField-Y");     #Pulled directly from data
    
    #LDR is used in ADDA  y default, chosen with -pol <...>
    #LDR formulation for isotropic material (avoid matrices/tensors)
    #https://arxiv.org/pdf/astro-ph/0311304
    wavelength = beam_spec[1];
    n  = 1.0;   #Refractive index   ############ NEEDS TO BE PARSED IN ###############
    wave_number = 2.0*math.pi/wavelength;
    alpha_cm = ( (3.0*pow(dipole_size,3))/(4.0*math.pi) )*( (n**2 -1)/(n**2 +2) ); ##### CHECK UNITS --> SHOULD BE CGS TO MATCH ADDA I THINK #####
    b1 = -1.8915316;
    b2 =  0.1648469;
    b3 = -1.7700004;
    #########################
    ## CORRECT THESE VALUES ---> FOR LAGUERRE CIRUCLAR POALRISATION NEEDED
    #########################
    S  = 1.0; #X or Y linearly polarised
    polarisability = alpha_cm / (1.0 +(alpha_cm/pow(dipole_size,3))*( (b1 +pow(n,2)*b2 +pow(n,2)*b3*S)*(pow(wave_number*dipole_size,2)) -(2.0/3.0)*1j*(pow(wave_number*dipole_size,3)) ));

    #Find force for each dipole
    posForce_set = [];
    for dip_ind in range(0,len(posEInternal_set)):
        if(dip_ind % 5 == 0):
            print("Dipole "+str(dip_ind)+"/"+str(len(posEInternal_set)));
        posEInternal = posEInternal_set[dip_ind];
        posForce     = [posEInternal[0], posEInternal[1], posEInternal[2]];     #Start with positions in the posForce_set
        #Add Fx, Fy, Fz
        for i in range(0,3):    #For each force component
            #F   =       1/2 Re(polarisability*Grad( |E|^2 ))
            #F_i = SUM_j 1/2 Re(polarisability*E_j*del_i*E_j_conj)
            force_comp = 0.0 +0j;
            for j in range(0,3):
                ####
                ## MAKE REQUEST CORRECT BEAM => APPLICABLE TO ANY BEAM TYPE ##
                ####
                ##
                ## ---> This wants incident beam at dipole positions?
                ##      --> This is absoltely fine if true
                ##  ---> Problem is the gradient is needed too ---> could assume bad gradient between dipole positions, OR just assume beam is also replicated fully in here, the same as is in ADDA
                ##
                #pos, wavelength, radial, azimuthal, alpha, beta, w_0
                E_partial = partial( Get_E_LaguerreGaussian_Beam, wavelength=beam_spec[1], radial=beam_spec[2], azimuthal=beam_spec[3], alpha=beam_spec[4], beta=beam_spec[5], w_0=beam_spec[6] );
                E_comp_partial      = E_partial( posEInternal[0:3] )[j];                            #Complex valued
                E_grad_comp_partial = gradient(E_partial, posEInternal[0:3], i)[j];     # "" ""
                force_comp += (0.5)*( polarisability*E_comp_partial*(E_grad_comp_partial.conjugate()) ).real;
            posForce.append(np.absolute(force_comp));
        #Add |F|^2 back in afterwards
        mod_F_squared = np.sqrt(np.sum(np.array(posForce[3:6])**2));    #Indices 3,4,5 are Fx,Fy,Fz here
        #print("     posForce[3:6]= "+str(posForce[3:6]));
        #print("     modF^2= "+str(mod_F_squared));
        posForce.insert(3, np.absolute(mod_F_squared));
        #Add to full set of dipoles
        posForce_set.append(posForce);
    return posForce_set;

def gradient(function, pos, offset_component):
    pos = np.array(pos);
    pos_offset = np.array([0.0, 0.0, 0.0]);
    pos_offset[offset_component] = sys.float_info.epsilon;    #Accuracy between (0,1) -> closer to 0 => more accuracy
    
    #Default gradient approx
    #numerator_vector = np.array(function(list(pos))) - np.array(function(list(pos+pos_offset)));
    #gradient_vector = numerator_vector / sys.float_info.epsilon;

    #Central difference approx
    numerator_vector = np.array(function(list(pos+pos_offset))) - np.array(function(list(pos-pos_offset)));
    gradient_vector = numerator_vector / 2.0*sys.float_info.epsilon;

    return gradient_vector;  #Vector gradient

def Get_GreensTensor(r, r_dash, wavelength):
    """
    . Returns the Green's tensor, linking interaction between two points

    . r = position to view field at [x,y,z]
    . r_dash = position of the dipole contributing [x,y,z]
    """
    k = 2.0*math.pi/wavelength;
    R = np.array(r-r_dash);
    R_mag = np.sqrt(np.sum(R**2));
    if(R_mag != 0.0):
        #If valid location not on top of dipole, calculate its matrix
        factor = (cmath.exp(1j*k*R_mag))/(4.0*math.pi*R_mag);
        matrix_a = np.identity(3) - ((np.outer(R,R))/(pow(R_mag,2)));
        matrix_b = ( (1j*k*R_mag -1.0)/(pow(k,2)*pow(R_mag,2)) )*( np.identity(3) -3.0*((np.outer(R,R))/(pow(R_mag,2))) );
        return factor*(matrix_a +matrix_b);
    else:
        #If zero distance (considering on top of dipole) => ignore it
        return np.zeros( (3,3) );

def Calculate_NearField_Point(r, beam_spec, posPolarisations, permitivity_free):
    """
    . Calculates the electric field at some point within the particle at position r
    . This electric field is a complex vector

    . r = position where electric field is being calculated at, only valid within the particle (other use the far -field scattering)
    . posPolarisations = list of positions and the complex vector polarisation at that point
    """
    print("Calculating near-field E for "+beam_spec[0]);
    k = 2.0*math.pi/beam_spec[1];
    #Convert all to numpy arrays
    posPolarisations = np.array(posPolarisations);
    for posPolarisation in posPolarisations:
        np.array(posPolarisation);
    
    #Pull out positions and polarisations separately -> Pos in micrometers (+possible CGS units)
    #####
    ## MOVE THIS OUT AND PARSE IN FOR MORE EFFICIENCY
    #####
    dip_positions        = posPolarisations[:,0:3];         #[ [ x, y, z], ... ]
    dip_polarisations_Re = posPolarisations[:,4::2];        #[ [Px_re, Py_re, Pz_re], ... ]
    dip_polarisations_Im = posPolarisations[:,5::2];        #[ [Px_im, Py_im, Pz_im], ... ]
    dip_polarisations = dip_polarisations_Re +1j*dip_polarisations_Im;

    E_incident = Get_E_LaguerreGaussian_Beam(r, beam_spec[1], beam_spec[2], beam_spec[3], beam_spec[4], beam_spec[5], beam_spec[6]);
    E_interaction = np.array([0+0j, 0+0j, 0+0j]);
    for dip_ind in range(0,len(posPolarisations)):
        if(dip_ind % 5 == 0):
            print("Interaction with dipole "+str(dip_ind)+"/"+str(len(posPolarisations)));
        E_interaction += np.matmul(Get_GreensTensor(r, dip_positions[dip_ind], beam_spec[1]),dip_polarisations[dip_ind]);
    E_interaction = E_interaction*(pow(k,2)/(permitivity_free));
    #####
    ## INTERACTION TERM IS LOOKING VERY LARGE -> CHECK ACCURACY
    #####
    return E_incident + E_interaction;

def main():
    """
    . Performs an ADDA calculation for a given shape
    . Returns polarisations from this
    . Uses these polarisations to find the E far-field over a given discretised region

    ==========
    **Note; You will still retrieve plots for [invalid beams] specified if the previous output_data folder exists,
        as it will simply read that data
    ==========
    """
    print("Program Started");
    #Shape used in ADDA
    input_shape         = "cube5x5x5.geom";   #capsule10x10x15    #sphere5x5x5    #cube5x5x5
    #Name of beam used by ADDA, WITH ALL the required parameters
    input_beam          = "plane";            #lminus 1   #plane
    output_folder       = "output_data";    #Folder where results are placed
    wavelength          = 1.064;            #In μm
    dipoles_per_lambda  = 15;               #Number of dipoles to have within a certain length
    dipole_size = wavelength/dipoles_per_lambda;
    permitivity_free = 8.85*pow(10,-12);
    k = 2.0*math.pi/wavelength;

    space_data_XY = [  #Specifies where to sample far-field grid -> assumes particle centred at origin (bounding box centred)
        [-3.5, 3.5, 41],  # In μm, [X_start, X_end, X_samples]
        [-3.5, 3.5, 41],  # "" ""
        [-0.0, 0.0, 1]
    ];
    space_data_ZX = [  #Specifies where to sample far-field grid -> assumes particle centred at origin (bounding box centred)
        [-3.5, 3.5, 41],  # In μm, [X_start, X_end, X_samples]
        #[-3.5, 3.5, 41],  # "" ""
        [-0.0, 0.0,  2],
        [-3.5, 3.5, 41]  # "" ""
    ];

    ########################################################################
    ### SPLIT THIS INTO OTT STYLE OPTION SELECTOR, NOT ONE BIG MESSY LIST ##
    ########################################################################

    #Shapes start next to eachother, then move further apart
    print("Calculations started...");

    shape_info, shape_width = PullSingleParticleGeom(input_shape);   #Info on dipoles not needed here, just shape width for exclusion zone
    print("Shape width= ",shape_width);

    #Run adda
    RunAdda(input_shape, input_beam, dipoles_per_lambda, 1.5, 0, wavelength, output_folder);

    #Get polarisations out
    paramDict = PullAddaData_Parameters("log_laguerre08_realPart"); #output_data/log
    wavelength = paramDict["lambda"];   #Matches wavelength used by ADDA
    posPolarisations_X = PullAddaData_DipolePosForces("output_data/DipPol-X");
    #posPolarisations_Y = PullAddaData_DipolePosForces("DipPol-Y_laguerre08_realPart");   #output_data/DipPol-Y
    posPolarisations_Y = PullAddaData_DipolePosForces("DipPol-Y_laguerre08_imagPart");

    laguerre_beam_spec  = ["laguerreGaussian",wavelength,0,8, 0.4*pow(10,-8)*(1.0+1j), 0.4*pow(10,-8)*(1.0-1j), 0.6*wavelength];
    planeWave_beam_spec = ["planeWave",wavelength, [0,0,0], "X"];
    print("Wavelength used= ",wavelength);

    #-------------------
    # Force Calculation
    #-------------------
    #E_near_field = Calculate_NearField_Point(np.array(posPolarisations_X[0])[0:3], laguerre_beam_spec, posPolarisations_X, permitivity_free);
    #print("E_near_field= "+str(E_near_field));
    #beam_posForce_set = Calculate_T_Averaged_Dipole_PosForce(dipole_size, laguerre_beam_spec);
    #make_quiver_plot(beam_posForce_set, dipole_size, dipole_size, dipole_size, should_normalise_arrows=True);  #### MAY NEED TO SWAP AXES, MEHSGRID PROBLEM POSSIBLY ###

    #-----------------
    # E_Incident Plot
    #-----------------
    field_data_XY = Generate_Beam(planeWave_beam_spec, space_data_XY);   #laguerre_beam_spec
    field_data_ZX = Generate_Beam(planeWave_beam_spec, space_data_ZX);
    Plot_2D_Array("Z", "norm", "E norm", field_data_XY, space_data_XY);
    plt.show();
    Plot_2D_Array("Y", "norm", "E norm", field_data_ZX, space_data_ZX);
    plt.show();
    #Sweep_Beam_Generation(space_data, wavelength);
    #plt.show();


    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## Testing ADDA Plane-Wave Implementation ##
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    #(1) Calculate incident field at each dipole and compare
    print("---Performing ADDA Plane-Wave Test...");
    pythonCalc_posFields = [];
    for dipole in posPolarisations_X:
        pythonCalc_posField = [dipole[0], dipole[1], dipole[2]];    #Consider same position
        pythonCalc_incidentField = Get_E_PlaneWave_Beam(dipole[0:3], wavelength, [0.0, 0.0, 0.0], "X");
        pythonCalc_posField.append( np.sum(np.power(np.array(pythonCalc_incidentField), 2)) );    #|E|^2
        pythonCalc_posField.append(pythonCalc_incidentField[0]);    #Ex
        pythonCalc_posField.append(pythonCalc_incidentField[1]);    #Ey
        pythonCalc_posField.append(pythonCalc_incidentField[2]);    #Ez
        pythonCalc_posFields.append(pythonCalc_posField);           #Add to main set
    #(2) Calulate incident + scattered field at each dipole and compare
    #pass
    #(3) Calculate t-averaged force at each dipole and compare
    #pass
    print("     ADDA Pulled Values= ",posPolarisations_X);  #### NEEDS TO COMPARE TO INCIDENT FIELD ####
    print("");
    print("");
    print("     Python Calc Values= ",pythonCalc_posFields);
    print("---ADDA Plane-Wave Test Finished");

    
    #------------------------------
    # E_Far-Field Plot -- Laguerre
    #------------------------------
    exclusionRadius = 16*dipole_size;   #16x16x16 sphere used
    #print("X Polarisations");
    #E_far_field_X = Calculate_FarField_Grid(posPolarisations_Y, exclusionRadius, space_data, k, permitivity_free);    #Get 3D array of far-field E
    print("Y Polarisations");
    #E_far_field_Y = Calculate_FarField_Grid(posPolarisations_Y, exclusionRadius, space_data_XY, k, permitivity_free);    #Get 3D array of far-field E
    #Plot this grid of values
    #Plot_2D_Array("Z", "norm", "|E|, "+input_shape+", beam="+input_beam+", X polarisation", E_far_field_X, space_data);
    #Plot_2D_Array("Z", "norm", "|E|, "+input_shape+", beam="+input_beam+", Y polarisation", E_far_field_Y, space_data_XY);
    #plt.show();

    
    #------------------
    # E_Far-Field Plot
    #------------------
    #exclusionRadius = shape_width*dipole_size;
    #print("X Polarisations");
    #E_far_field_X = Calculate_FarField_Grid(posPolarisations_Y, exclusionRadius, space_data, k, permitivity_free);    #Get 3D array of far-field E
    #print("Y Polarisations");
    #E_far_field_Y = Calculate_FarField_Grid(posPolarisations_Y, exclusionRadius, space_data, k, permitivity_free);    #Get 3D array of far-field E
    #Plot this grid of values
    #Plot_2D_Array("Z", "norm", "|E|, "+input_shape+", beam="+input_beam+", X polarisation", E_far_field_X, space_data);
    #Plot_2D_Array("Z", "norm", "|E|, "+input_shape+", beam="+input_beam+", Y polarisation", E_far_field_Y, space_data);

    """
    #-----------------
    # E_Internal Plot
    #-----------------
    posEinternal = PullAddaData_DipolePosEInternal("output_data/IntField-Y");
    pos_int_field, E_int_field = Format_PosVariable_List(posEinternal);
    pos_int_field_layer1 = pos_int_field[0* 5*5 : 1* 5*5];    #For 5x5x5 shape
    E_int_field_layer1     = E_int_field[0* 5*5 : 1* 5*5];    #For 5x5x5 shape
    """

    print("Program ended");

main();