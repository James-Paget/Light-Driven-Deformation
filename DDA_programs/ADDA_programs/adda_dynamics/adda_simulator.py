import subprocess
import numpy as np
from make_animation_and_quiver import make_position_animation, make_quiver_plot

def PullAddaData():
    """
    . Pulls all (1 )position and (2) force data from an run instance of adda
    """
    dipolePositions = PullAddaData_DipolePositions("box.geom"); #### NEED TO MAKE THIS USE THE CUSTOM OBJECT WE WILL DEAL WITH -> ALWAYS INSTEAD ####
    dipoleForces = PullAddaData_DipoleForces("RadForce-Y");
    #pass
    return dipolePositions, dipoleForces;

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
                        print("Could NOT cast parameter to float");
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
                        print("Could NOT cast parameter to float");
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
                        print("Could NOT cast parameter to float");
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
    

def PullAddaData_DipoleForces(filename):
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

def FetchParticleResultantForce(posForce_set):
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

def FetchParticleOffset(timeStep, particleDetails, posForce_set):
    """
    . Calculates the position offset for the particle from force experienced
    . posForce_set = The posForce objects associated with each lattice point, formatted as 
        sets of the list [x,y,z, |F|^2,Fx,Fy,Fz]

    .**Note; ADDA only offsets the beam centre => it can be offset by -ve the particle's offset
    
    ########
    ## ASSUMED MASS=1, PERHAPS ASSUME LOWER (AND JUST USE AS A SCALING FACTOR) OR ACTUALLY CALCULATE
    ########
    """
    #Initialise values
    mass = 1.0;
    force = np.array(FetchParticleResultantForce(posForce_set));
    acc = force/mass;
    vel = np.array(particleDetails["velocity_history"][len(particleDetails["velocity_history"])-1]);  #Start from previous velocity
    pos = np.array([0.0, 0.0, 0.0]);  #Positional offset from just the force the particle experiences (not including current offset)

    #Calculate new values
    vel += acc*timeStep;
    pos += vel*timeStep

    #Store values
    particleDetails["velocity_history"].append(vel);
    particleDetails["offset_history"].append( 
        particleDetails["offset_history"][len(particleDetails["offset_history"]) -1] + pos 
    );
    #Return offset from JUST the force (particleDetails["offset_history"] has the new resultant offset stored)
    return list(pos);


def run_adda(input_file, output_folder, beam_offset, m1Re, path_to_adda):
    ####
    ## BEAM PROPOGATION DIRECTION DOES NOT APPEAR TO BE CHANGING RESULTS ??? Possibly just hard to see or some cancelling effect [SEEMS VERY UNLIKELY] ???
    ###
    adda_command = "./adda -shape read " + input_file + " -prop 1 0 0" + " -beam_center " + str(beam_offset[0]) + " " + str(beam_offset[1]) + " " + str(beam_offset[2]) + " -m " + str(m1Re) + " 0 -store_force -dir " + output_folder
    adda_command = adda_command.split(" ")
    result = subprocess.run(adda_command, cwd=path_to_adda+"/src/seq", stdout=subprocess.DEVNULL)
    return result.returncode


def main():

    particleDetails = {"offset_history":[np.array([0.0, 0.0, 0.0])], "velocity_history":[np.array([0.0, 0.0, 0.0])]}
    timeStep = 1e-3;      #Time step of the simulation, more accurate of middling values (~10^-4/10^-5)

    path_to_adda = "/home/james/Desktop/Abraham-Minkowski_Project/adda"  #/home/james/Desktop/Abraham-Minkowski_Project/ADDA_Sim
    input_file = "box.geom"
    output_folder = "output_folder"
    num_iterations = 100
    beam_offset = [0,0,0]
    m = 1.5

    for _ in range(num_iterations):

        return_code = run_adda(input_file, output_folder, beam_offset, m, path_to_adda)
        if return_code != 0:
            print("ADDA error")
            break

        result_folder_path = path_to_adda + "/src/seq/" + output_folder + "/"

        #Also assumes a [single] particle
        paramDict = PullAddaData_Parameters(result_folder_path + "log")        
        posForce_set = PullAddaData_DipoleForces(result_folder_path + "RadForce-Y")
        FetchParticleOffset(timeStep, particleDetails, posForce_set);         #How much the force moves this particle
        beam_offset = -particleDetails["offset_history"][ len(particleDetails["offset_history"])-1 ];   #Beam offset by opposite of what particle has been offset by in its frame of ref.
        
    x_spacing = paramDict["dipole_size"]
    y_spacing = paramDict["dipole_size"]
    z_spacing = paramDict["dipole_size"]

    data = np.swapaxes([particleDetails["offset_history"]],0,1)
    beam_wavelength = paramDict["lambda"]
    make_position_animation(data, timeStep, should_repeat=True, padding=beam_wavelength)

    posForce_set = PullAddaData_DipoleForces(result_folder_path + "RadForce-Y")
    make_quiver_plot(posForce_set, x_spacing, y_spacing, z_spacing)



main();