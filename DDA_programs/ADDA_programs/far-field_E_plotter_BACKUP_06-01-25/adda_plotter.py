import subprocess
import numpy as np
import random
import matplotlib.pyplot as plt
import math
import cmath
import scipy.integrate as integrate
import scipy.special as sp
from scipy.interpolate import griddata
from matplotlib import colormaps as cm
import matplotlib as mpl
from functools import partial
import sys
import re
from make_animation_and_quiver import make_quiver_plot

permitivity_free = 8.85*pow(10,-12);

def PullAddaData_Parameters(filename):
    """
    . Pulls other miscellaneous parameters from the log file adda generates
    . Values are;
    [lambda, [refractive_index_real, refractive_index_imag], dipole_size, [beam_center_x, beam_center_y, beam_center_z], [beam_dir_x, beam_dir_y, beam_dir_z]]
    ##
    ## Beam center pull not added yet but can be done easily
    ##
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
                case "CoupleConstant: ":
                    value = -1.0 +0.0j;
                    #Further refine parameter
                    unref_value = line[char_index:];
                    try:
                        # Extracts the real and imaginary components and their respect powers from the line, then converts to a complex value
                        value_terms = re.findall(r'[+-]?\d+(?:\.\d+)?', unref_value)
                        value = complex( float(value_terms[0]+"e"+value_terms[1]), float(value_terms[2]+"e"+value_terms[3]) )
                    except:
                        print("Could NOT cast parameter to float: polarisability");
                    parameter_set.update({"polarisability":value});
                    startIndex = char_index+1;
                case "NEW PARAMETER CASE":
                    pass
    return parameter_set;
  
def PullAddaData_General(filename):
    """
    . Pulls data from any ADDA generated data set, structured as separate particle information per line
    . Works for parameters such as forces (x, y, z, |F|^2, Fx, Fy, Fz), polarisations (x, y, z, |P|^2, Px, Py, Pz), etc

    .** Note; Assumes the file is in the same folder as this python script
    .** Note; These are in RAW values (not relative)
    """
    #Load file
    file = open(filename, "r");
    #Ignore description lines (first 1)
    for i in range(0,1):
        file.readline();
    #Begin reading and storing following dipole positions
    generalParam_set = [];
    for line in file:
        generalParams = [];
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
                    print(" -Could NOT cast general parameter to float");
                generalParams.append(value);   #Note; Disclude the space
                termIndex+=1;
                startIndex = char_index+1;
            elif(char_index == len(line)-1):
                #Handle last value separatly as it does NOT have a space after it
                value = -1.0;
                try:
                    value = float(line[startIndex:]);
                except:
                    print(" -Could NOT cast general parameter to float");
                generalParams.append(value);   #Note; Disclude the space
        generalParam_set.append(generalParams);
    return generalParam_set;

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

def RunAdda(input_shape, input_beam, input_beam_args, dpl, n_real, n_imag, wavelength, output_folder):
    beam_args = ""
    for arg in input_beam_args:
        beam_args += " "+str(arg)
    refractive_indices = str(n_real)+" "+str(n_imag)+" "+str(n_real)+" "+str(n_imag);
    adda_command = "./adda -dpl " +str(dpl)+ " -lambda " + str(wavelength)+ " -store_force -shape read " + input_shape + " -m "+ refractive_indices + " -store_beam -store_int_field -store_dip_pol -beam "+input_beam+beam_args+" -dir " + str(output_folder);
    adda_command = adda_command.split(" ");
    print("=== ADDA Log ===");
    result = subprocess.run(adda_command, stdout=subprocess.DEVNULL);

def PullSingleParticleGeom(filename):
    """
    . Fetch information about a .geom file that describes a single shape
    . Fetch parameters include;
        - Particle info (dipole positions)
        - Particle width (extremal)
    """
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

def Generate_Beam(beam_spec, pos_set):
    """
    . Generates the electric on the positions specified

    . beam_spec = [beam_type, <params>]
    . pos_set = [ [x,y,z], ... ] for each point of interest
    """
    field_data = [];
    for pos in pos_set:
        E_field = [0.0, 0.0, 0.0];
        match beam_spec[0]:
            case "test":
                E_field = [pos[0], pos[1], pos[2]]
            case "laguerre":
                #beam_spec = [beam_type, wavelength, azi, radial, alpha, beta, w_0]
                E_field = Get_E_LaguerreGaussian_Beam(pos, beam_spec[1], beam_spec[2], beam_spec[3], beam_spec[4], beam_spec[5], beam_spec[6])
            case "plane":
                #beam_spec = [beam_type, wavelength, azi, radial, alpha, beta, w_0]
                E_field = Get_E_PlaneWave_Beam(pos, beam_spec[1], beam_spec[2], beam_spec[3])
            case _:
                print("Invalid beam type; ",beam_spec[0])
        field_data.append( [pos[0], pos[1], pos[2], np.sqrt( pow(E_field[0],2) + pow(E_field[1],2) + pow(E_field[2],2) ), E_field[0], E_field[1], E_field[2]] )
    return field_data

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
    ##
    ## Implement with a Jones vector
    ##
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

def Fetch_Polarisability(wavelength, dipole_size, n):
    """
    . Calculates the polarisability of a dipole
    
    SH_DDA form used;
    polarisability = cm / (1 - (2 / 3) * 1j * k ** 3 * a0/(4*np.pi*8.85e-12))
    """
    # LDR formulation
    # https://arxiv.org/pdf/astro-ph/0311304
    permitivity_dielectric = (n.real**2)+(n.imag**2)
    volume = dipole_size**3
    cm = (permitivity_free)*(volume)*(permitivity_dielectric -1) / (permitivity_dielectric +2)
    b1 = -1.8915316
    b2 = 0.1648469
    b3 = -1.7700004
    ####
    ## THIS IS USING THE FACT THAT PLANE WAVE IS ONLY EVER X OR Y POLARISED, AND ALL BEAMS ARE Z PROPOGATING
    ##      => S always =0 for the plane wave cases
    ##          BUT will not equal 0 for the circular polarisation of laguerre for instance
    ##  **THIS WILL NEED CHANGING IN THE FUTURE**
    ####
    S  = (0.0*1.0)**2 # NOTE; Assumes only ever Z propogation
    k = 2.0*np.pi/wavelength
    alpha_LDR = cm / ( 1 + (cm/volume)*( (b1 +b2*permitivity_dielectric +b3*permitivity_dielectric*S)*((k*dipole_size)**2) -1j*(2.0/3.0)*((k*dipole_size)**3) ) )


    # SH_DDA formualtion
    ep1 = (n)**2
    ep2 = 1.0
    epm = 1.333 # water
    a0 = (4 * np.pi * 8.85e-12) * (volume) * ((ep1 - epm) / (ep1 + 2*epm))
    a = a0 / (1 - (2 / 3) * 1j * k ** 3 * a0/(4*np.pi*8.85e-12))


    print("cm = ",cm)
    print("alpha_LDR = ",alpha_LDR)
    print("a = ",a)
    return a

def Calculate_T_Averaged_Force(dip_positions, dip_polarisations, dipole_indices, n, dipole_size, polarisability, beam_spec, mode=0):
    """
    . Calculates the total force vector on each dipole
    . This uses the calculation for polarisation, internal E field and beam incident E field on each dipole is known (found in ADDA, or manually if wanted)
    . This force requires the E field gradient to be known, hence the beam used must be recreated here in Python (as well as be available in ADDA)
    . Force is complex for generalisation

    . dipole_set = list of dipoles to experience a force from the beam + scattering
    . dipole_indices = the indices of the dipoles from dipole_set that form the particle of interest; hence only these are summed over for the force (will be all indices for a single particle system)
    . n = complex refractive index of particle
    """
    print("Calculating time-averaged forces for "+beam_spec[0]);

    dipole_forces = np.zeros( (len(dip_positions), 7), dtype=complex)
    
    # LDR is used in ADDA  y default, chosen with -pol <...>
    #wavelength = beam_spec[1]
    #polarisability = Fetch_Polarisability(wavelength, dipole_size, n.real) #### NOTE; GIVEN BY LOG, NOT NEEDED

    for i in range(len(dipole_indices)):
        if(i % 20 == 0):
            print("Dipole "+str(i)+"/"+str(len(dipole_indices)))
        # For every dipole in particle of interest
        dip_ind = dipole_indices[i]
        pos = dip_positions[dip_ind]
        dip_force = np.zeros(7, dtype=complex)
        dip_force[:3] = pos
        E_func_totalNear = partial(Calculate_TotalField_NearField_Point, dip_positions=dip_positions, dip_polarisations=dip_polarisations, beam_spec=beam_spec, mode=0)     # Now is just a function of position to sample at
        for comp in range(3):
            # For each component
            # (1) Get (E.conj) gradient
            # (2) Get E force
            E_conj_gradient = conj_gradient(E_func_totalNear, pos, comp)   # Gradient of full E field at the postiion given
            dip_force[4+comp] += np.dot(dip_polarisations[dip_ind], (E_conj_gradient)) # Set [Fx, Fy, Fz]
            dip_force[4+comp] = ( (0.5)*dip_force[4+comp] ).real            # Set [Fx, Fy, Fz]
        dip_force[3] = np.dot(dip_force[4:].real, dip_force[4:].conjugate())     # Set |F|^2
        dipole_forces[i] = dip_force
    if(mode == 1):
        dipole_forces[:,:] = dipole_forces[:,:].real
        dipole_forces = np.array(dipole_forces, dtype=float)
    return dipole_forces

def gradient(function, pos, offset_component):
    """
    . Gives the gradient of the function result
    """
    pos = np.array(pos);
    pos_offset = np.array([0.0, 0.0, 0.0]);
    pos_offset[offset_component] = sys.float_info.epsilon;    # Accuracy between (0,1) -> closer to 0 => more accuracy

    # Central difference approx
    func_upper = ( np.array(function(list(pos+pos_offset)))[4:] )
    func_lower = ( np.array(function(list(pos-pos_offset)))[4:] )
    numerator_vector = func_upper-func_lower
    gradient_vector = numerator_vector / (2.0*sys.float_info.epsilon);
    return gradient_vector;  #Vector gradient

def conj_gradient(function, pos, offset_component):
    """
    . Gives the gradient of the conjugate of the function result
    """
    step = pow(10,-4) #sys.float_info.epsilon;
    pos = np.array(pos);
    pos_offset = np.array([0.0, 0.0, 0.0]);
    pos_offset[offset_component] = step;    # Accuracy between (0,1) -> closer to 0 => more accuracy

    # Central difference approx
    func_upper = ( np.array(function(list(pos+pos_offset)))[4:] ).conjugate()
    func_lower = ( np.array(function(list(pos-pos_offset)))[4:] ).conjugate()
    numerator_vector = func_upper-func_lower
    gradient_vector = numerator_vector / (2.0*step);

    return gradient_vector;  #Vector gradient

def Calculate_ScatteredField_NearField(dipole_set, request_indices, beam_spec, gridCalc=False):
    """
    . Calculates the scattered near field at the requested points from a full set of dipoles given

    . dipole_set = all the dipoles involved in the system to be included in the scattering in a list.
                    NOTE; This list should be [x,y,z,  |P|^2,  Px,r,Px,i,  Py,r,Py,i,  Pz,r,Pz,i], where P is the polarisation of each dipole
    . request_indicies = the index (within the dipole_set list) of all the positions to get the near field for

    ####
    ## MODIFY TO ALSO ALLOW ARBITRARY POINTS TO BE QUERIED FOR NEAR FIELD E FIELD --> WORRY ABOUT SINGULARITIES IF TOO CLOSE 
    ## TO OTHER DIPOLES / WITHIN ITS VOLUME PROBLEM?
    ##
    ## Pattern close but off in magnitude compared to ADDA -> Seems to be purely units problem -> should adjust for this
    ####
    """
    print("Calculating near-field E for "+beam_spec[0]);

    k = 2.0*math.pi/beam_spec[1]

    # Convert all to numpy arrays
    dipole_set = np.array(dipole_set);
    for dipole_set_value in dipole_set:
        np.array(dipole_set_value)

    # Separate important info
    dip_positions        = dipole_set[:,0:3]         #[ [ x, y, z], ... ]
    dip_polarisations_Re = dipole_set[:,4::2]        #[ [Px_re, Py_re, Pz_re], ... ]
    dip_polarisations_Im = dipole_set[:,5::2]        #[ [Px_im, Py_im, Pz_im], ... ]
    dip_polarisations = dip_polarisations_Re +1j*dip_polarisations_Im   # Combine into a single complex list

    # Calculate this scattered near-field at each requested point
    E_interactions = np.zeros((len(request_indices),3), dtype=complex)
    for j in range(len(request_indices)):
        r_ind = request_indices[j]
        ####
        ## Bug-fixing display
        if(j % 50 == 0):
            print("     dipole "+str(j)+"/"+str(len(request_indices)))
        ## Bug-fixing display
        ####
        E_interactions[j] = Calculate_ScatteredField_NearField_Point(dip_positions, dip_polarisations, dip_positions[r_ind], beam_spec)
    return E_interactions   # [E_x, E_y, E_z], Complex values

def Get_ScatteringAmplitude(wavelength, n, dip_positions, dip_polarisations):
    """
    n = unit vector to point considered
    """
    k = (2.0*np.pi)/(wavelength)
    term1 = np.identity(3) - np.outer(n,n)  # Matrix
    term2 = np.zeros(3, dtype=complex)      # Vector
    for i in range(len(dip_positions)):
        term2 += dip_polarisations[i]*np.exp(-1j*k*np.dot(dip_positions[i],n))
    return -1j*pow(k,3)*( np.dot(term1,term2) ) # Vector result

def Get_GreensTensor(r, r_dash, wavelength):
    """
    . Returns the Green's tensor, linking interaction between two points
    . r = position to view field at [x,y,z]
    . r_dash = position of the dipole contributing [x,y,z]
    """
    k = 2.0*np.pi/(wavelength);
    R = np.array(r-r_dash);
    R_mag = np.linalg.norm(R)   #np.sqrt(np.sum(R**2));
    if(abs(R_mag) > sys.float_info.epsilon):
        #If valid location not on top of dipole, calculate its matrix
        term1 = (cmath.exp(1j*k*R_mag))/(4.0*np.pi*R_mag);
        matrix_a = np.identity(3) - ((np.outer(R,R))/(pow(R_mag,2)));
        term2 = (1j*k*R_mag -1.0)/(pow(k,2)*pow(R_mag,2));
        matrix_b = term2*( np.identity(3) -3.0*((np.outer(R,R))/(pow(R_mag,2))) );
        return term1*(matrix_a +matrix_b);
    else:
        #If zero distance (considering on top of dipole) => ignore it
        return np.zeros( (3,3) );

def Calculate_ScatteredField_NearField_Point(dip_positions, dip_polarisations, pos, beam_spec):
    """
    . Near field scattering using Green's tensor form
    """
    k = 2.0*np.pi/(beam_spec[1])
    # Calculate this scattered near-field at each requested point
    E_interaction = np.zeros(3, dtype=complex)
    for i in range( len(dip_positions) ):
        E_interaction += np.matmul(Get_GreensTensor(pos, dip_positions[i], beam_spec[1]), dip_polarisations[i])
    E_interaction = E_interaction*(pow(k,2)/(permitivity_free))
    # NOTE; Unsure about this following line; ADDA states conversion from CGS->SI for |E|^2 is by multiplying by this value however 
    # the values given by ADDA stored files are supposed to be dimensionless, hence no CGS conversion required, but by another factor?
    E_interaction *= (4.0*np.pi*permitivity_free)
    return E_interaction   # [E_x, E_y, E_z], Complex values

def Calculate_ScatteredField_FarField_Point(dip_positions, dip_polarisations, pos, beam_spec):
    """
    . Far-field scattering with the scattering amplitude form
    """
    k = 2.0*np.pi/beam_spec[1]
    # Calculate this scattered near-field at each requested point
    r_mag = np.linalg.norm(pos)
    n = pos/r_mag   # Unit vector to point
    scat_amp = Get_ScatteringAmplitude(beam_spec[1], n, dip_positions, dip_polarisations)
    E_interaction = scat_amp*( (cmath.exp(1j*k*r_mag))/(-1j*k*r_mag) )
    return E_interaction   # [E_x, E_y, E_z], Complex values

def Calculate_TotalField_NearField_Point(pos, dip_positions, dip_polarisations, beam_spec, mode=0):
    """
    . Combines a scattered near field and incident field (from the beam) to give a total field within the particle
    """
    ####
    ## RENAME GEN BEAM TO CALACULATE_..._BEAM()
    ####
    E_beam = np.array(Generate_Beam(beam_spec, [pos]), dtype=complex)[0]    # Get 0th as only considering a single point
    E_scattered = Calculate_ScatteredField_NearField_Point(dip_positions, dip_polarisations, pos, beam_spec)
    E_total = np.zeros(len(E_beam), dtype=complex)
    E_total[:3] = pos                                   # Set position
    E_total[4:] += E_scattered[:]  +E_beam[4:]          # Set E components
    E_total[3]  = np.dot(E_total[4:], E_total[4:].conjugate()) # Set |E|^2, same as ADDA
    
    if(mode == 1):
        E_total[3:] = E_total[3:].real    # Take real parts of all components for plots if wanted
        E_total = np.array(E_total, dtype=float)
    return E_total

def Fetch_Separated_PolarisationSet(dipole_set):
    """
    . Considers a dipole_set of form [ [x,y,z,|P|^2,Px.r,Px.i,Py.r,Py.i,Pz.r,Pz.i], ... ]
    . Returns a dipoles formatted as two lists [ [x,y,z], ...] and [ [Px, Py, Pz], ...]
    """
    # Convert all to numpy arrays
    dipole_set = np.array(dipole_set);
    for dipole_set_value in dipole_set:
        np.array(dipole_set_value)

    # Separate important info
    dip_positions        = dipole_set[:,0:3]         #[ [ x, y, z], ... ]
    dip_polarisations_Re = dipole_set[:,4::2]        #[ [Px_re, Py_re, Pz_re], ... ]
    dip_polarisations_Im = dipole_set[:,5::2]        #[ [Px_im, Py_im, Pz_im], ... ]
    dip_polarisations = dip_polarisations_Re +1j*dip_polarisations_Im   # Combine into a single complex list
    return dip_positions, dip_polarisations

def Fetch_Plane_Subset(data_set, view_parameters, include_planePos=False):
    """
    . Considers a set of positions and extracts all the particle matching some conditions for being within a specified plane

    . view_parameters = [Coordinate_Plane_Index, Plane_Value, Plot_Type]
        Where;
            'Coordinate_Plane_Index' is either 0,1,2 for X,Y,Z plane to view in
            'Plane_Value' is the exact plane value to view at, e.g. [2, 0.0] => View XY plane (perp. to Z) at Z=0.0
    """
    # Get valid data points to this plot
    coord_set = []  # Holds 2 other coords (perp to plane) to plot against
    value_set = []  # Holds scalars, or vectors
    index_set = []  # Holds the index of found data
    for data_index in range(len(data_set)):
        # NOTE; Cutoff of 10^-4 of a micrometer is arbitrary here
        ##
        ## MAYBE TAKE MEASURE BASED ON LATTICE INDEX # RATHER THAN POSITION -> LESS ERROR PRONE + ALREADY ON GRID LATTICE => IS FINE
        ##
        data_point = data_set[data_index]
        if(abs(data_point[view_parameters[0]] - view_parameters[1]) < pow(10,-4)):   # Within some buffer of target plane
            # Add position to plot
            coords_temp = []
            for i in range(3):
                if(include_planePos):
                    coords_temp.append( data_point[i] )
                else:
                    if(i != view_parameters[0]):
                        coords_temp.append( data_point[i] )
                
            coord_set.append(coords_temp)
            value_set.append(data_point[3:])
            index_set.append(data_index)
    return coord_set, value_set, index_set

def Plot_3dLattice_CrossSection(data_set, lattice_spacing, view_parameters=[2, 0.0, "colour"], labels={"X":"", "Y":"", "Title":"", "args":[]}):
    """
    . Takes in a set of 3D coordinates on a cubic lattice
    . Plots a specified cross section of this lattice

    . view_parameters = [Coordinate_Plane_Index, Plane_Value, Plot_Type]
        Where;
            'Coordinate_Plane_Index' is either 0,1,2 for X,Y,Z plane to view in
            'Plane_Value' is the exact plane value to view at, e.g. [2, 0.0] => View XY plane (perp. to Z) at Z=0.0
            'Plot_Type' specifies how the data will be plotted, e.g. a 2D 'colour' plot, a 2D 'vector' plot, ...
    . labels = [Xlabel, Ylabel, Title, <extraArgs>]

    . Input data as; [ [x,y,z,<extraArgs>], ... ] for each dipole/position to consider
        Where <extraArgs> can be a scalar, vector, etc
    """
    # Get relevent set of data for the plane wanted
    #print("data_set = ",data_set)
    coord_set, value_set, index_set = Fetch_Plane_Subset(data_set, view_parameters)
    #print("coord_set = ",coord_set)
    #print("value_set = ",value_set)
    match view_parameters[2]:
        case "colour":
            value_set = np.array(value_set)[:, 0]
        case "vector":
            value_set = np.array(value_set)[:, 3:]
        case _:
            print("Plot type not recognised; ",view_parameters[2])

    # If any particles found do plot, if not leave it
    if( len(coord_set) > 0 ):
        # Arrange points found into a grid
        x_set = np.array(coord_set)[:,0]
        y_set = np.array(coord_set)[:,1]

        # Get bounding box for plot -> Number of points
        x_lattice_min = min(x_set) / lattice_spacing
        x_lattice_max = max(x_set) / lattice_spacing
        x_lattice_diff = int(np.floor(x_lattice_max-x_lattice_min))
        y_lattice_min = min(y_set) / lattice_spacing
        y_lattice_max = max(y_set) / lattice_spacing
        y_lattice_diff = int(np.floor(y_lattice_max-y_lattice_min))

        # Shift all lattice points to start at origin, and get in terms of lattice spacings not raw position
        x_set = np.floor( (x_set/lattice_spacing) -x_lattice_min )
        y_set = np.floor( (y_set/lattice_spacing) -y_lattice_min )

        # Create empty grid to populate
        grid_values = np.zeros( (x_lattice_diff, y_lattice_diff) )

        # Go through all grid points
        for j in range( grid_values.shape[0] ):
            for i in range( grid_values.shape[1] ):
                # Look for matches in relevent data
                matchFound = False
                for p in range( len(value_set) ):
                    xMatch = (int(x_set[p]) == i)
                    yMatch = (int(y_set[p]) == j)
                    if(xMatch and yMatch):
                        matchFound = True
                        grid_values[j,i] = value_set[p]
                        break
                if(not matchFound):             # Not needed, already zeroed
                    grid_values[j,i] = 0.0      #

        # Display plot
        fig, ax = plt.subplots()
        match view_parameters[2]:
            case "colour":
                im_plot=ax.imshow(
                    grid_values,
                    extent=(min(np.array(coord_set)[:,0]), max(np.array(coord_set)[:,0]), min(np.array(coord_set)[:,1]), max(np.array(coord_set)[:,1])),  # Set axes limits
                    origin='lower',   # Ensure correct orientation
                    cmap='viridis'    # Choose a colormap
                )
                plt.colorbar(im_plot, ax=ax)
                pass
            case "vector":
                ##
                ## Needs to be implemented
                ##
                pass
            case _:
                print("Plot type not recognised; ",view_parameters[2])

        # Label plot
        plt.xlabel( labels["X"] )
        plt.ylabel( labels["Y"] )
        plt.title(  labels["Title"] )
        plt.show()
    else:
        print("No particles to plot, perhaps the lattice is not in the plane specified?")


####
## MAKE BEAM SPEC A DICTIONARY
####
def main_plot_adda_E_incident(input_shape_name, input_beam_name, input_beam_args, output_folder, wavelength, dpl, material_refractive_index):
    """
    . Plots electric field at each dipole according to an ADDA calculation
    . Only includes electric field from the incident beam alone

    . dpl = dipoles per lambda

    NOTE; wavelength given in micrometers units
    """
    k = 2.0*math.pi/wavelength;
    dipole_size = wavelength/dpl;

    # Run ADDA simulation
    RunAdda(input_shape_name, input_beam_name, input_beam_args, dpl, material_refractive_index.real, material_refractive_index.imag, wavelength, output_folder);

    # Read output data from ADDA
    paramDict = PullAddaData_Parameters("output_data/log")
    lattice_spacing = paramDict["dipole_size"]
    
    #posEField_ADDA_Xpol = PullAddaData_General("output_data/IncBeam-X")     # Beam incident
    posEField_ADDA_Ypol = PullAddaData_General("output_data/IncBeam-Y")     #

    match input_beam_name:
        case "plane":
            beam_spec = [input_beam_name, wavelength, [0,0,0], "Y"]
        case "laguerreGaussian":
            beam_spec = [input_beam_name,wavelength,input_beam_args[0],input_beam_args[1], 0.4*pow(10,-8)*(1.0+1j), 0.4*pow(10,-8)*(1.0-1j), 0.6*wavelength];
        case _:
            print("No beam name match found, hence no beam_spec assigned; ",input_beam_name)
    pos_set = np.array(posEField_ADDA_Ypol)[:,:3]
    #posEField_Python_Xpol = Generate_Beam( [input_beam_name, wavelength, [0,0,0], "X"], pos_set)
    posEField_Python_Ypol = Generate_Beam( beam_spec, pos_set)

    # Plot data that has been read
    Plot_3dLattice_CrossSection(posEField_ADDA_Ypol, lattice_spacing, view_parameters=[2, lattice_spacing/2.0, "colour"], labels={"X":"X(10^-6m)", "Y":"Y(10^-6m)", "Title":"ADDA, Ypol, Plane Z="+str(lattice_spacing/2.0)+"(10^-6m)", "args":[]})
    Plot_3dLattice_CrossSection(posEField_Python_Ypol, lattice_spacing, view_parameters=[2, lattice_spacing/2.0, "colour"], labels={"X":"X(10^-6m)", "Y":"Y(10^-6m)", "Title":"Python, Ypol, Plane Z="+str(lattice_spacing/2.0)+"(10^-6m)", "args":[]})

def main_plot_adda_E_internal(input_shape_name, input_beam_name, input_beam_args, output_folder, wavelength, dpl, material_refractive_index):
    """
    . Plots electric field at each dipole according to an ADDA calculation
    . Includes contribution from incident beam and scattered parts from dipoles, within the particle at each dipole (internal)

    . dpl = dipoles per lambda

    NOTE; wavelength given in micrometers units
    """
    k = 2.0*math.pi/wavelength;
    dipole_size = wavelength/dpl;

    # Run ADDA simulation
    RunAdda(input_shape_name, input_beam_name, input_beam_args, dpl, material_refractive_index.real, material_refractive_index.imag, wavelength, output_folder);

    # Read output data from ADDA
    paramDict = PullAddaData_Parameters("output_data/log")
    lattice_spacing = paramDict["dipole_size"]
    
    posEField_ADDA_Xpol = PullAddaData_General("output_data/IntField-X")        # Internal field of each dipole
    posPolarisation_ADDA_Xpol = PullAddaData_General("output_data/DipPol-X")    # Polarisation of each dipole

    # Python calculations for scattered (near) field
    view_parameters = [2, lattice_spacing/2.0, "colour"]
    plane_beam_spec = [input_beam_name, wavelength, [0,0,0], "X"]
    pos_set, value_set, index_set = Fetch_Plane_Subset(posPolarisation_ADDA_Xpol, view_parameters, include_planePos=True)     # Fetch plane to be plotted from the data

    posEField_Python_Xpol = np.zeros( (len(pos_set),7), dtype=float )
    dip_positions, dip_polarisations = Fetch_Separated_PolarisationSet(posPolarisation_ADDA_Xpol)
    for i in range(len(pos_set)):
        posEField_Python_Xpol[i] = Calculate_TotalField_NearField_Point(pos_set[i], dip_positions, dip_polarisations, plane_beam_spec, mode=1)

    # Plot true field (ADDA) and Python calculation for comparison
    Plot_3dLattice_CrossSection(posEField_ADDA_Xpol, lattice_spacing, view_parameters, labels={"X":"X(10^-6m)", "Y":"Y(10^-6m)", "Title":"ADDA, Ypol, Plane Z="+str(lattice_spacing/2.0)+"(10^-6m)", "args":[]})
    Plot_3dLattice_CrossSection(posEField_Python_Xpol, lattice_spacing, view_parameters, labels={"X":"X(10^-6m)", "Y":"Y(10^-6m)", "Title":"Python, Ypol, Plane Z="+str(lattice_spacing/2.0)+"(10^-6m)", "args":[]})

def main_plot_adda_E_force(input_shape_name, input_beam_name, input_beam_args, output_folder, wavelength, dpl, material_refractive_index):
    """
    . Plots the force on each dipole from the ADDA calcualtion on a particle with a plane wave and from the Python 
        version found using ADDA values
    . This only works for plane waves as ADDA only supports force calcualtion for plane waves, and so can only be compared here

    NOTE; wavelength, lattice points, etc, given in micrometer units
    """
    k = 2.0*math.pi/wavelength;
    dipole_size = wavelength/dpl;

    # Run ADDA simulation
    RunAdda(input_shape_name, input_beam_name, input_beam_args, dpl, material_refractive_index.real, material_refractive_index.imag, wavelength, output_folder);

    # Read output data from ADDA
    paramDict = PullAddaData_Parameters("output_data/log")
    lattice_spacing = paramDict["dipole_size"]
    polarisability = paramDict["polarisability"]
    
    posForce_ADDA_Xpol = PullAddaData_General("output_data/RadForce-X")        # Internal field of each dipole
    posPolarisation_ADDA_Xpol = PullAddaData_General("output_data/DipPol-X")    # Polarisation of each dipole

    # Python calculations force on each
    view_parameters = [2, lattice_spacing/2.0, "colour"]
    plane_beam_spec = [input_beam_name, wavelength, [0,0,0], "X"]
    pos_set, value_set, index_set = Fetch_Plane_Subset(posPolarisation_ADDA_Xpol, view_parameters, include_planePos=True)     # Fetch plane to be plotted from the data
    dip_positions, dip_polarisations = Fetch_Separated_PolarisationSet(posPolarisation_ADDA_Xpol)
    posForce_Python_Xpol = Calculate_T_Averaged_Force(dip_positions, dip_polarisations, range(len(posPolarisation_ADDA_Xpol)), material_refractive_index, lattice_spacing, polarisability, plane_beam_spec, mode=1)

    # Plot forces for ADDA and Python calculations
    Plot_3dLattice_CrossSection(posForce_ADDA_Xpol, lattice_spacing, view_parameters, labels={"X":"X(10^-6m)", "Y":"Y(10^-6m)", "Title":"ADDA, Ypol, Plane Z="+str(lattice_spacing/2.0)+"(10^-6m)", "args":[]})
    Plot_3dLattice_CrossSection(posForce_Python_Xpol, lattice_spacing, view_parameters, labels={"X":"X(10^-6m)", "Y":"Y(10^-6m)", "Title":"Python, Ypol, Plane Z="+str(lattice_spacing/2.0)+"(10^-6m)", "args":[]})




#planeWave_beam_spec = [input_beam_name, wavelength, [0,0,0], "X"];
#laguerre_beam_spec  = ["laguerreGaussian",wavelength,0,8, 0.4*pow(10,-8)*(1.0+1j), 0.4*pow(10,-8)*(1.0-1j), 0.6*wavelength];

#-------------
# Program Run 
#-------------
print("Program Started");

if int(len(sys.argv)) != 2:
    print("Invalid arguments; format as python <COMMAND>");
    print("Valid <COMMAND> args are; 'plot_adda_Einc', ...");
    sys.exit()

match(sys.argv[1]):
    case 'plot_adda_Einc':
        #main_plot_adda_E_incident("sphere16x16x16.geom", "plane", "", "output_data", 1.064, 15, 1.5+0j);
        main_plot_adda_E_incident("sphere8x8x8.geom", "plane", "", "output_data", 1.064, 15, 1.5+0j);
    case 'plot_adda_Eint':
        main_plot_adda_E_internal("sphere8x8x8.geom", "plane", "", "output_data", 1.064, 15, 1.5+0j);
    case 'plot_adda_force':
        main_plot_adda_E_force("sphere8x8x8.geom", "plane", "", "output_data", 1.064, 15, 1.5+0j);
    case _:
        print("Invalid <COMMAND>; ",sys.argv[1])

print("Program ended");
