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
        generalParam_set.append(np.array(generalParams));
    return np.array(generalParam_set);

def RunAdda(input_shape, input_beam, input_beam_args, dpl, n_real, n_imag, wavelength, output_folder):
    beam_args = ""
    for arg in input_beam_args:
        beam_args += " "+str(arg)
    refractive_indices = str(n_real)+" "+str(n_imag)+" "+str(n_real)+" "+str(n_imag);
    adda_command = "./adda -dpl " +str(dpl)+ " -lambda " + str(wavelength)+ " -store_force -shape read " + input_shape + " -m "+ refractive_indices + " -store_beam -store_int_field -store_dip_pol -beam "+input_beam+beam_args+" -dir " + str(output_folder);
    adda_command = adda_command.split(" ");
    print("=== ADDA Log ===");
    result = subprocess.run(adda_command, stdout=subprocess.DEVNULL);

def Calculate_IncidentField(beam_spec, pos_set):
    field_data = []
    for i in range(len(pos_set)):
        Einc_field = Calculate_IncidentField_Point(beam_spec, pos_set[i])
        field_data.append(Einc_field)
    return np.array(field_data)

def Calculate_IncidentField_Point(beam_spec, pos):
    """
    . Generates the incident electric field from the beam used on the positions specified

    . beam_spec = {"type", ...}
    . pos_set = [ [x,y,z], ... ] for each point of interest
    """

    def Calculate_IncidentField_PlaneWave_Point(pos, wavelength, beam_center, polarisation):
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
            return np.array([ctemp.real, 0.0, 0.0]);
        elif(polarisation=="Y"):
            return np.array([0.0, ctemp.real, 0.0]);
        else:
            print("Invalid polarisation for plane wave: "+str(polarisation));
            return np.array([0.0, 0.0, 0.0]);

    def Calculate_IncidentField_LaguerreGaussian_Point(pos, wavelength, radial, azimuthal, alpha, beta, w_0):
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
        #return [E_X_Comp[0].real, E_Y_Comp[0].real, E_Z_Comp[0].real];
        return [E_X_Comp[0], E_Y_Comp[0], E_Z_Comp[0]];


    E_field = [0.0, 0.0, 0.0];
    match beam_spec["type"]:
        case "test":
            E_field = [pos[0], pos[1], pos[2]]
        case "LG":
            #beam_spec = [beam_type, wavelength, azi, radial, alpha, beta, w_0]
            E_field = Calculate_IncidentField_LaguerreGaussian_Point(pos, beam_spec["wavelength"], beam_spec["azi"], beam_spec["radial"], beam_spec["alpha"], beam_spec["beta"], beam_spec["w0"])
        case "plane":
            #beam_spec = [beam_type, wavelength, beam_center, polarisation]
            E_field = Calculate_IncidentField_PlaneWave_Point(pos, beam_spec["wavelength"], beam_spec["beam_center"], beam_spec["polarisation"])
        case _:
            print("Invalid beam type; ",beam_spec["type"])
    Einc_field = np.array([
        pos[0], 
        pos[1], 
        pos[2], 
        np.dot( E_field, np.conjugate(E_field) ),
        E_field[0], 
        E_field[1], 
        E_field[2]
    ])
    return Einc_field

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
    #### \/
    S  = (0.0*1.0)**2 # NOTE; Assumes only ever Z propogation
    k = 2.0*np.pi/wavelength
    alpha_LDR = cm / ( 1 + (cm/volume)*( (b1 +b2*permitivity_dielectric +b3*permitivity_dielectric*S)*((k*dipole_size)**2) -1j*(2.0/3.0)*((k*dipole_size)**3) ) )


    # SH_DDA formualtion
    # ep1 = (n)**2
    # ep2 = 1.0
    # epm = 1.333 # water
    # a0 = (4 * np.pi * 8.85e-12) * (volume) * ((ep1 - epm) / (ep1 + 2*epm))
    # a = a0 / (1 - (2 / 3) * 1j * k ** 3 * a0/(4*np.pi*8.85e-12))

    return alpha_LDR

def Calculate_T_Averaged_Force(dip_positions, dip_polarisations, dipole_indices, beam_spec, mode=0):
    """
    . Calculates the total force vector on each dipole
    . This uses the calculation for polarisation, internal E field and beam incident E field on each dipole is known (found in ADDA, or manually if wanted)
    . This force requires the E field gradient to be known, hence the beam used must be recreated here in Python (as well as be available in ADDA)
    . Force is complex for generalisation

    . dipole_set = list of dipoles to experience a force from the beam + scattering
    . dipole_indices = the indices of the dipoles from dipole_set that form the particle of interest; hence only these are summed over for the force (will be all indices for a single particle system)
    . n = complex refractive index of particle
    """
    print("Calculating time-averaged forces for "+beam_spec["type"]);

    dipole_forces = np.zeros( (len(dip_positions), 7), dtype=complex)

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
    """
    print("Calculating near-field E for "+beam_spec["type"]);

    k = 2.0*math.pi/beam_spec["wavelength"]

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
        if(j % 50 == 0):
            print("     dipole "+str(j)+"/"+str(len(request_indices)))
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
    k = 2.0*np.pi/(beam_spec["wavelength"])
    # Calculate this scattered near-field at each requested point
    E_interaction = np.zeros(3, dtype=complex)
    for i in range( len(dip_positions) ):
        E_interaction += np.matmul(Get_GreensTensor(pos, dip_positions[i], beam_spec["wavelength"]), dip_polarisations[i])
    E_interaction = E_interaction*(pow(k,2)/(permitivity_free))
    # NOTE; Unsure about this following line; ADDA states conversion from CGS->SI for |E|^2 is by multiplying by this value however 
    # the values given by ADDA stored files are supposed to be dimensionless, hence no CGS conversion required, but by another factor?
    E_interaction *= (4.0*np.pi*permitivity_free)
    return E_interaction   # [E_x, E_y, E_z], Complex values

def Calculate_ScatteredField_FarField_Point(dip_positions, dip_polarisations, pos, beam_spec):
    """
    . Far-field scattering with the scattering amplitude form
    """
    k = 2.0*np.pi/beam_spec["wavelength"]
    # Calculate this scattered near-field at each requested point
    r_mag = np.linalg.norm(pos)
    n = pos/r_mag   # Unit vector to point
    scat_amp = Get_ScatteringAmplitude(beam_spec["wavelength"], n, dip_positions, dip_polarisations)
    E_interaction = scat_amp*( (cmath.exp(1j*k*r_mag))/(-1j*k*r_mag) )
    return E_interaction   # [E_x, E_y, E_z], Complex values

def Calculate_TotalField_NearField_Point(pos, dip_positions, dip_polarisations, beam_spec, mode=0):
    """
    . Combines a scattered near field and incident field (from the beam) to give a total field within the particle
    """
    E_beam = np.array(Calculate_IncidentField(beam_spec, [pos]), dtype=complex)[0]    # Get 0th as only considering a single point
    E_scattered = Calculate_ScatteredField_NearField_Point(dip_positions, dip_polarisations, pos, beam_spec)
    E_total = np.zeros(len(E_beam), dtype=complex)
    E_total[:3] = pos                                   # Set position
    E_total[4:] += E_scattered[:]  +E_beam[4:]          # Set E components
    E_total[3]  = np.dot(E_total[4:], E_total[4:].conjugate()) # Set |E|^2, same as ADDA
    
    if(mode == 1):
        E_total[3:] = E_total[3:].real    # Take real parts of all components for plots if wanted
        E_total = np.array(E_total, dtype=float)
    return E_total








def Calculate_TAveragedForce_Point(point_index, dip_positions, dip_polarisations, beam_spec):
    """
    . point_index       = index of the point forces are being found for in dip_positions/dip_polarisations lists
    . dip_positions     = positions of all dipoles in full system
    . dip_polarisations = polarisations of all dipoles in full system
    . beam_spec         = {"type", ...}
        type = plane, LG, ...
    """
    dip_posForce = np.zeros(7, dtype=complex)
    dip_posForce[:3] = dip_positions[point_index]

    for i in range(len(dip_positions)): # Go through every other dipole, sum all interaction term forces
        if(i != point_index):
            E_func_totalNear = partial(Calculate_TotalField_NearField_Point, dip_positions=dip_positions, dip_polarisations=dip_polarisations, beam_spec=beam_spec)     # Now is just a function of position to sample at
            for j in range(3):  # For X,Y,Z component
                E_conj_gradient = conj_gradient(E_func_totalNear, dip_posForce[:3], j)   # Vector E gradient wrt each component X,Y,Z
                dip_posForce[4:] = ( (0.5)*dip_polarisations[i]*E_conj_gradient ).real           # Set [Fx, Fy, Fz]
    dip_posForce[3] = np.dot(dip_posForce[4:].real, np.conjugate(dip_posForce[4:]))     # Set |F|^2
    return np.array(dip_posForce, dtype=float)  # Forces, hence should be real values anyway (cast to remove plotting warning)

def Calculate_InternalEField_Point(dip_polarisation, polarisability):
    """
    . Internal field for a given dip polarisation
    """
    E_internal = dip_polarisation/polarisability
    return E_internal   # [E_x, E_y, E_z], Complex values


def get_beam_spec(input_beam_name, args={"wavelength":1.0}):
    #planeWave_beam_spec = ["plane", wavelength, [0,0,0], "X"];
    #laguerre_beam_spec  = ["LG",wavelength,0,8, 0.4*pow(10,-8)*(1.0+1j), 0.4*pow(10,-8)*(1.0-1j), 0.6*wavelength];
    match input_beam_name:
        case "plane":
            beam_spec = {"type":input_beam_name, "wavelength":args["wavelength"], "beam_center":[0,0,0], "polarisation":"Y"}
        case "LG":
            beam_spec = {"type":input_beam_name, "wavelength":args["wavelength"], "azi":0, "radial":8, "alpha":0.4*pow(10,-8)*(1.0+1j), "beta":0.4*pow(10,-8)*(1.0-1j), "w0":0.6*args["wavelength"]}
        case _:
            print("No beam name match found, hence no beam_spec assigned; ",input_beam_name)
            sys.exit()
    return beam_spec


def main_plot_adda_E_incident(input_shape_name, input_beam_name, input_beam_args, output_folder, wavelength, dpl, material_refractive_index):
    """
    . Plots electric field at each dipole according to an ADDA calculation
    . Only includes electric field from the incident beam alone

    . dpl = dipoles per lambda

    NOTE; wavelength given in micrometers units
    """
    # Run ADDA simulation
    #RunAdda(input_shape_name, input_beam_name, input_beam_args, dpl, material_refractive_index.real, material_refractive_index.imag, wavelength, output_folder);
    beam_spec = get_beam_spec(input_beam_name, args={"wavelength":wavelength})

    # Read output data from ADDA
    paramDict = PullAddaData_Parameters("output_data/log")
    lattice_spacing = paramDict["dipole_size"]
    
    posEField_ADDA_Ypol = PullAddaData_General("output_data/IncBeam-Y")      # Beam incident

    pos_set = np.array(posEField_ADDA_Ypol)[:,:3]
    posEField_Python_Ypol = Calculate_IncidentField( beam_spec, pos_set)

    # Plot data that has been read
    pos_raw = posEField_ADDA_Ypol[:,0:3]
    E_ADDA_Ypol_values_raw   = np.array(posEField_ADDA_Ypol)[:,3]
    E_Python_Ypol_values_raw = np.array(posEField_Python_Ypol)[:,3]
    for i in range(0,8,1):
        plot_params = {"axis":2, "lattice_value":i}
        plot_lattice_crosssection(lattice_spacing, pos_raw, E_ADDA_Ypol_values_raw  , plot_params)
        plot_lattice_crosssection(lattice_spacing, pos_raw, E_Python_Ypol_values_raw, plot_params)

def main_plot_adda_E_internal(input_shape_name, input_beam_name, input_beam_args, output_folder, wavelength, dpl, material_refractive_index):
    """
    . Plots electric field at each dipole according to an ADDA calculation
    . Includes contribution from incident beam and scattered parts from dipoles, within the particle at each dipole (internal)

    . dpl = dipoles per lambda
    . beam_spec = *IMPORTANT* ADDA run will assume this is just a plance wave here, as plane wave forces only can be calculated, therefore not generalised here (could be added if wanted here however, inside RunAdda())

    NOTE; wavelength given in micrometers units
    """
    # Run ADDA simulation
    RunAdda(input_shape_name, input_beam_name, input_beam_args, dpl, material_refractive_index.real, material_refractive_index.imag, wavelength, output_folder);

    # Read output data from ADDA
    paramDict = PullAddaData_Parameters("output_data/log")
    lattice_spacing = paramDict["dipole_size"]
    polarisability = paramDict["polarisability"]
    posEInt_ADDA_Xpol = PullAddaData_General("output_data/IntField-X")          # Internal field of each dipole; ADDA's result
    
    # Extract ADDA data to be used
    posPolarisation_ADDA_Xpol = PullAddaData_General("output_data/DipPol-X")    # Polarisation of each dipole; Used in Python functions here to try replicate the internal field ADDA found (will be required for force calculation later, which ADDA cannot perform in general)
    positions_raw = posPolarisation_ADDA_Xpol[:,0:3]
    dip_pol_real_raw = np.array(posPolarisation_ADDA_Xpol)[:,4::2]
    dip_pol_imag_raw = np.array(posPolarisation_ADDA_Xpol)[:,5::2]
    dip_pol_raw      = dip_pol_real_raw + 1j*dip_pol_imag_raw
    dipole_number = len(posPolarisation_ADDA_Xpol)

    # Generate Python output
    posEInt_Python_Xpol = np.zeros( (dipole_number,7), dtype=complex )
    for i in range(dipole_number):
        posEInt_Python_Xpol[i, :3] = positions_raw[i]
        posEInt_Python_Xpol[i, 4:] = Calculate_InternalEField_Point(dip_pol_raw[i], polarisability)#Calculate_TotalField_NearField_Point(positions_raw[i], positions_raw, dip_pol_raw, plane_beam_spec, mode=1)
        posEInt_Python_Xpol[i, 3]  = np.dot(posEInt_Python_Xpol[i, 4:], np.conjugate(posEInt_Python_Xpol[i, 4:]))

    # Plot true field (ADDA) and Python calculation for comparison
    for i in range(0,8,1):
        Eint_ADDA_Xpol_values_raw   = np.array(posEInt_ADDA_Xpol)[:,3]
        Eint_Python_Xpol_values_raw = np.array(posEInt_Python_Xpol)[:,3]
        plot_params = {"axis":2, "lattice_value":i}
        plot_lattice_crosssection(lattice_spacing, positions_raw, Eint_ADDA_Xpol_values_raw  , plot_params)
        plot_lattice_crosssection(lattice_spacing, positions_raw, Eint_Python_Xpol_values_raw, plot_params)

def main_plot_adda_E_force(input_shape_name, input_beam_name, input_beam_args, output_folder, wavelength, dpl, material_refractive_index):
    """
    . Plots the force on each dipole from the ADDA calcualtion on a particle with a plane wave and from the Python 
        version found using ADDA values
    . This only works for plane waves as ADDA only supports force calcualtion for plane waves, and so can only be compared here

    NOTE; wavelength, lattice points, etc, given in micrometer units
    """
    # Run ADDA simulation
    RunAdda(input_shape_name, input_beam_name, input_beam_args, dpl, material_refractive_index.real, material_refractive_index.imag, wavelength, output_folder);

    # Read output data from ADDA
    paramDict = PullAddaData_Parameters("output_data/log")
    lattice_spacing = paramDict["dipole_size"]
    polarisability = paramDict["polarisability"]
    wavelength = paramDict["lambda"]
    posForce_ADDA_Xpol = PullAddaData_General("output_data/RadForce-X")          # Force field at each dipole; ADDA's result
    
    # Extract ADDA data to be used
    posPolarisation_ADDA_Xpol = PullAddaData_General("output_data/DipPol-X")    # Polarisation of each dipole; Used in Python functions here to try replicate the force field ADDA found (as ADDA can only find forces for plane wave, therefore we verify validity on the plane wave then extend to general beams to reduce errors as much as possible)
    positions_raw = posPolarisation_ADDA_Xpol[:,0:3]
    dip_pol_real_raw = np.array(posPolarisation_ADDA_Xpol)[:,4::2]
    dip_pol_imag_raw = np.array(posPolarisation_ADDA_Xpol)[:,5::2]
    dip_pol_raw      = dip_pol_real_raw + 1j*dip_pol_imag_raw
    dipole_number = len(posPolarisation_ADDA_Xpol)

    # Generate Python output
    posForce_Python_Xpol = np.zeros( (dipole_number,7), dtype=complex )
    ##
    ## Can reduce this to just consider the indicies you want to plot -> For now plotting everything is required for full test
    ##
    plane_beam_spec = [input_beam_name, wavelength, [0,0,0], "X"]
    for i in range(dipole_number):
        if(i%10==0):
            print("Progress: "+str(i)+"/"+str(dipole_number))
        posForce_Python_Xpol[i] = Calculate_TAveragedForce_Point(i, positions_raw, dip_pol_raw, plane_beam_spec)

    # Plot true field (ADDA) and Python calculation for comparison
    for i in range(0,8,1):
        posForce_ADDA_Xpol_values_raw   = np.array(posForce_ADDA_Xpol)[:,3]
        posForce_Python_Xpol_values_raw = np.array(posForce_Python_Xpol)[:,3]
        plot_params = {"axis":2, "lattice_value":i}
        plot_lattice_crosssection(lattice_spacing, positions_raw, posForce_ADDA_Xpol_values_raw  , plot_params)
        plot_lattice_crosssection(lattice_spacing, positions_raw, posForce_Python_Xpol_values_raw, plot_params)

    # k = 2.0*math.pi/wavelength;
    # dipole_size = wavelength/dpl;

    # # Run ADDA simulation
    # RunAdda(input_shape_name, input_beam_name, input_beam_args, dpl, material_refractive_index.real, material_refractive_index.imag, wavelength, output_folder);

    # # Read output data from ADDA
    # paramDict = PullAddaData_Parameters("output_data/log")
    # lattice_spacing = paramDict["dipole_size"]
    # #polarisability = paramDict["polarisability"]
    
    # posForce_ADDA_Xpol = PullAddaData_General("output_data/RadForce-X")
    # posPolarisation_ADDA_Xpol = PullAddaData_General("output_data/DipPol-X")    # Polarisation of each dipole

    # # Python calculations force on each
    # view_parameters = [2, lattice_spacing/2.0, "colour"]
    # plane_beam_spec = [input_beam_name, wavelength, [0,0,0], "X"]
    # pos_set, value_set, index_set = Fetch_Plane_Subset(posPolarisation_ADDA_Xpol, view_parameters, include_planePos=True)     # Fetch plane to be plotted from the data
    # dip_positions, dip_polarisations = Fetch_Separated_PolarisationSet(posPolarisation_ADDA_Xpol)
    # posForce_Python_Xpol = Calculate_T_Averaged_Force(dip_positions, dip_polarisations, range(len(posPolarisation_ADDA_Xpol)), plane_beam_spec, mode=1)

    # # Plot forces for ADDA and Python calculations
    # pos_raw = posForce_ADDA_Xpol[:,0:3]
    # Eforce_ADDA_Xpol_values_raw   = np.array(posForce_ADDA_Xpol)[:,3]
    # Eforce_Python_Xpol_values_raw = np.array(posForce_ADDA_Xpol)[:,3]
    # for i in range(0,8,1):
    #     plot_params = {"axis":2, "lattice_value":i}
    #     plot_lattice_crosssection(lattice_spacing, pos_raw, Eforce_ADDA_Xpol_values_raw  , plot_params)
    #     plot_lattice_crosssection(lattice_spacing, pos_raw, Eforce_Python_Xpol_values_raw, plot_params)



def plot_lattice_crosssection(lattice_spacing, positions_raw, values_raw, plot_params):
    # (1) Bin the positions into indicies 0-N, lattice spacings
    positions_lattice_offset = np.round(2.0*np.array(positions_raw)/lattice_spacing)/2.0 # NOTE; np.array() to ensure deep copy made
    positions_lattice = positions_lattice_offset -np.array([min(positions_lattice_offset[:,0]), min(positions_lattice_offset[:,1]), min(positions_lattice_offset[:,2])])
    # (2) Pick an X,Y or Z layer to plot (between 0-N)
    bounds_lattice = [
        [int(min(positions_lattice[:,0])), int(max(positions_lattice[:,0]))], # X lower/upper
        [int(min(positions_lattice[:,1])), int(max(positions_lattice[:,1]))], # Y lower/upper
        [int(min(positions_lattice[:,2])), int(max(positions_lattice[:,2]))]  # Z lower/upper
    ]
    plot_indicies  = [] # Contains indices of relevent data (e.g. in the correct axis layer)
    plot_axis  = plot_params["axis"]
    plot_value = plot_params["lattice_value"]
    for i in range(len(positions_lattice)):
        if( positions_lattice[i,plot_axis] == plot_value ):
            plot_indicies.append(i)
    # (3) Setup plot space, fill it with non-zero values in the data set given
    plot_axis_X1 = int((plot_axis+1)%3)  # e.g. If plotting the Y=a plane, then get plot_axis=1, so X1=(1+1)%3=2=Z axis, and X2=(1+2)%3=0=X axis
    plot_axis_X2 = int((plot_axis+2)%3)  #
    X1, X2 = np.meshgrid( 
        (lattice_spacing)*(np.array(range(bounds_lattice[plot_axis_X1][0], bounds_lattice[plot_axis_X1][1]+1, 1)) +min(positions_lattice_offset[:,plot_axis_X1]) ),
        (lattice_spacing)*(np.array(range(bounds_lattice[plot_axis_X2][0], bounds_lattice[plot_axis_X2][1]+1, 1)) +min(positions_lattice_offset[:,plot_axis_X2]) )
    )
    X3 = np.zeros( (bounds_lattice[plot_axis_X1][1]-bounds_lattice[plot_axis_X1][0]+1, bounds_lattice[plot_axis_X2][1]-bounds_lattice[plot_axis_X2][0]+1) ) # Value magnitudes plotted
    for i in plot_indicies:
        X3[int(positions_lattice[i, plot_axis_X1]), int(positions_lattice[i, plot_axis_X2])] = values_raw[i]
    # (4) Plot data
    plt.pcolor(X1,X2,X3, cmap=mpl.colormaps['viridis'])    #viridis, #binary_r, #binary
    plt.xlabel("Axis Comp. "+str(plot_axis_X1)+" [m]")
    plt.ylabel("Axis Comp. "+str(plot_axis_X2)+" [m]")
    plt.colorbar()
    plt.title(f"Plot Parameters={plot_axis}, {plot_value}")
    plt.show()


#-------------
# Program Run 
#-------------
print("Program Started");

if int(len(sys.argv)) != 2:
    print("Invalid arguments; format as python <COMMAND>");
    print("Valid <COMMAND> args are; 'plot_adda_Einc', ...");
    sys.exit()

match(sys.argv[1]):
    #
    # Compares quantities calcualted in ADDA to Python, hence generalised forces calculated in Python should be accurate to ADDA's plane-wave only forces
    #
    case 'compare_adda_Einc':
        #main_plot_adda_E_incident("sphere16x16x16.geom", "plane", "", "output_data", 1.064, 15, 1.5+0j);
        main_plot_adda_E_incident("sphere16x16x16.geom", "plane", "", "output_data", 1.064, 15, 1.5+0j);
    case 'compare_adda_Eint':
        main_plot_adda_E_internal("sphere16x16x16.geom", "plane", "", "output_data", 1.064, 15, 1.5+0j);
    case 'compare_adda_force':
        main_plot_adda_E_force("sphere5x5x5.geom", "plane", "", "output_data", 1.064, 15, 1.5+0j);
    
    case "plot_test_validate":
        print("Running Test...")
        log_dir  = "output_data/log"
        data_dir = "output_data/DipPol-X"
        paramDict = PullAddaData_Parameters(log_dir)
        lattice_spacing = paramDict["dipole_size"]
        pulled_data = np.array(PullAddaData_General(data_dir))
        positions_raw = pulled_data[:,0:3]  # NOTE; Assumes format of [x, y, z, |V|^2, Re(V.x), Im(V.x), Re(V.y), Im(V.y), Re(V.z), Im(V.z)] here
        values_mod2_raw = pulled_data[:,3]  #
        for i in range(0,8):
            plot_params = {"axis":2, "lattice_value":i}
            plot_lattice_crosssection(lattice_spacing, positions_raw, values_mod2_raw, plot_params) # e.g. Plot the plane Z=4 using {2,4}, where this refers to the relative lattice position
    case _:
        print("Invalid <COMMAND>; ",sys.argv[1])

print("Program ended")