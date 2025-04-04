import subprocess
import numpy as np
import random
import matplotlib.pyplot as plt
import math

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
                case "NEW PARAMETER CASE":
                    pass
    return parameter_set;

def PullAddaData_DipolePositions(filename):
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

def SeparateDipolePosForces_2_X_Particles(posForces_joined):
    """
    . Separates a joined posForce set into 2 sets when you have JUST 2 particles separated in the X-axis ONLY
    . 0(n) operation, not the worst, not the best
    """
    split_set = [[], []];
    for dipole_index in range(0, len(posForces_joined)):
        dipole_posForce = posForces_joined[dipole_index];
        if(dipole_posForce[0] <= 0):    #If the X position is less than 0, it is one of the particles
            split_set[0].append(dipole_posForce);
        else:                           #If its X is greater than 0, it is the other particle
            split_set[1].append(dipole_posForce);
    return split_set;

def FetchResultantForce_NParticles(particles_posForce_set):
    particles_resultant_forces = [];
    for i in range(0, particles_posForce_set.size()):
        resultant_force = FetchResultantForce(particles_posForce_set[i]);
        particles_resultant_forces.append(resultant_force);
    return particles_resultant_forces;

def RunAdda(input_file, dpl, n_real, n_imag, wavelength, output_folder):
    adda_command = "./adda -dpl " +str(dpl)+ " -lambda " + str(wavelength)+ " -shape read " + input_file + " -m "+ str(n_real)+" "+str(n_imag)+" "+str(n_real)+" "+str(n_imag) + " -store_force -dir " + str(output_folder);
    adda_command = adda_command.split(" ");
    result = subprocess.run(adda_command, stdout=subprocess.DEVNULL);
def CreateSingleParticleGeom(type):
    """
    . Creates a shape.geom file for a single particle
    . This particle can be scaled with a system e.g. it stays a constant real size, even as dipole widths get smaller

    . lattice_size = size of a single dipole
    """
    singleFile = open("shape.geom", "w");
    match(type):
        case("sphere"):
            """
            lattice_number = math.floor(0/lattice_size);
            for j in range(0, lattice_number):
                for i in range(0, lattice_number):
                    if( math.sqrt( pow(i-lattice_number/2.0, 2) + pow(j-lattice_number/2.0, 2) + pow(k-lattice_number/2.0, 2) ) <= lattice_number/2.0 ):
                        dipole_pos = i+" "+j+" "+k+"\n";
                        singleFile.write(dipole_pos);
            """
        case("cube"):
            for k in range(0, 5):
                for j in range(0, 5):
                    for i in range(0, 5):
                        dipole_pos = str(i)+" "+str(j)+" "+str(k)+"\n";
                        singleFile.write(dipole_pos);
    singleFile.close();

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

def GenerateGeomFile(singleShape_info, singleShape_width, separation):
    """
    . Creates a shapeCombined.geom file that is the combination of two particles, whose relative dipole positions are pulled from 
    the shape.geom file
    . separation = Number of lattice spacings between the two particles, where 0=exactly touching

    . Note; The comments at the start of the shape.geom file should be removed before using this function to read this file
    """
    #Overwrite shapeCombined.geom file
    combinedFile = open("shapeCombined.geom","w");
    combinedFile.write("Nmat=2\n");
    for n in range(0,2):                                                #For the 0th and 1st shape (only ever deal with 2 shapes)
        for dipole_pos in singleShape_info:
            for pos_component_index in range(0,len(dipole_pos)):
                pos_component = dipole_pos[pos_component_index];
                component_string_value = str(pos_component);
                if(pos_component_index == 0):                           #Shift separation in X direction
                    component_string_value = str(pos_component +n*(singleShape_width +separation));
                combinedFile.write(component_string_value+" ");
            combinedFile.write(str(n+1));                               #To specify which particle this dipole belongs to
            combinedFile.write("\n");
        combinedFile.write("\n");                                       #Another line of spacing to make file easier to read
    combinedFile.close();

def GenerateGeomFile_circular(singleShape_info, singleShape_width, radial_distance, particle_number):
    """
    . Creates a combination of N particles in a ring
    . separation = Number of lattice spacings between the two particles, where 0=exactly touching

    . Note; The comments at the start of the shape.geom file should be removed before using this function to read this file
    """
    #Overwrite shapeCombined.geom file
    combinedFile = open("shapeCombined.geom","w");
    combinedFile.write("Nmat=2\n");
    for n in range(0,2):                                                #For each shape
        for dipole_pos in singleShape_info:
            for pos_component_index in range(0,len(dipole_pos)):
                pos_component = dipole_pos[pos_component_index];
                component_string_value = str(pos_component);
                if(pos_component_index == 0):                           #Shift separation in X direction
                    component_string_value = str(pos_component +n*(singleShape_width +separation));
                combinedFile.write(component_string_value+" ");
            combinedFile.write(str(n+1));                               #To specify which particle this dipole belongs to
            combinedFile.write("\n");
        combinedFile.write("\n");                                       #Another line of spacing to make file easier to read
    combinedFile.close();

def PlotData(plot_data, wavelength, dipole_size):
    """
    . Plots data about force experienced at different separations, for 2 different particles
    """
    plot_data = np.array(plot_data);
    particle_separation = plot_data[:,0]*dipole_size;
    particle_1_data = plot_data[:,1];
    particle_2_data = plot_data[:,2];

    plt.plot(particle_separation, particle_1_data, label ="LHP λ= "+str(wavelength)+"μm");
    #plt.plot(particle_separation, particle_2_data, label ="RHP λ= "+str(wavelength)+"μm");
    #plt.show();

def main():
    """
    . Calculates the force between two particles (spheres or cubes) as they are brought together in order judge

    . Note;
    -- Beam always sits at origin => always sits exactly halfway between the two particles
    """
    print("Program Started");
    input_file = "shapeCombined.geom";  #Shape used in DDA (for both particles)
    output_folder = "output_data";      #Folder where results are placed
    iterations = 100;                   #How many lattice spacings of separation to do
    separation_step = 5;                #How many dipoles the space increases by each cycle

    dipoles_per_lambda = 15;

    plot_data = [];     #Stores [res_force_xORy_particle1, res_force_xORy_particle2] for each separation, the first element being at 0 separation
    #CreateSingleParticleGeom("cube");
    singleShape_info, singleShape_info_width = PullSingleParticleGeom("shape.geom");    #Pull this data into a list once then use for all further operations

    #Shapes start next to eachother, then move further apart
    print("Calculations started...");
    print("single shape X-width= ",singleShape_info_width);
    plt.figure();
    for wavelength_factor in range(1, 2):
        wavelength    = 1.064;#*wavelength_factor*0.5;#1.0*wavelength_factor;   #In μm
        for iter in range(0, iterations):
            if(iter % 10 == 0):
                print("     Iter "+str(iter)+"/"+str(iterations));
            GenerateGeomFile(singleShape_info, singleShape_info_width, iter*separation_step);

            dipole_size = wavelength/dipoles_per_lambda;
            print("Dipole Size= ",dipole_size);
            RunAdda("shapeCombined.geom", dipoles_per_lambda, 1.5, 0, wavelength, "output_data");

            paramDict = PullAddaData_Parameters("output_data/log");
            particles_posForce_joinedSet = PullAddaData_DipolePosForces("output_data/RadForce-Y");

            particles_posForce_splitSet = SeparateDipolePosForces_2_X_Particles(particles_posForce_joinedSet);

            numpy_particles_posForce_splitSet_p1  = np.array(particles_posForce_splitSet[0]);
            #numpy_particles_posForce_splitSet_p1 /= paramDict["dipole_size"];
            #numpy_particles_posForce_splitSet_p1  = np.round(numpy_particles_posForce_splitSet_p1, 2);
            resultant_force_x_p1 = np.sum( numpy_particles_posForce_splitSet_p1[:,4] );

            numpy_particles_posForce_splitSet_p2  = np.array(particles_posForce_splitSet[1]);
            #numpy_particles_posForce_splitSet_p2 /= paramDict["dipole_size"];
            #numpy_particles_posForce_splitSet_p2  = np.round(numpy_particles_posForce_splitSet_p2, 2);
            resultant_force_x_p2 = np.sum( numpy_particles_posForce_splitSet_p2[:,4] );

            force_element = [iter*separation_step, resultant_force_x_p1, resultant_force_x_p2]; #Store force in X OR Y direction experienced by each particle at each separation
            #force_element = [iter*separation_step, random.random(), random.random()];
            plot_data.append(force_element);
        print("Plotting separation data...");
        PlotData(plot_data, wavelength, dipole_size);
        plot_data = [];
    plt.xlabel("Separation of particles in μm");
    plt.ylabel("X Force (10^-13N)");
    plt.title("X Force between 2 particles, dpl= "+str(dipoles_per_lambda));
    plt.legend();
    plt.show();
    print("Program ended");

main();