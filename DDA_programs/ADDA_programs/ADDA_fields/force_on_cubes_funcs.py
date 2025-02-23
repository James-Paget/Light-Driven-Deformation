import numpy as np
import subprocess
from matplotlib import pyplot as plt
from file_reading_funcs import PullAddaData_Parameters, PullAddaData_DipoleForces

def make_pos_force_grid(data, dipole_size):
    data = np.array(data)

    # Convert positions to an integer grid to avoid floats. Nee
    data[:,:3] = np.round(data[:,:3]/[dipole_size, dipole_size, dipole_size], 1)

    grid_x_min = np.min(data[:,0])
    grid_x_max = np.max(data[:,0])
    grid_y_min = np.min(data[:,1])
    grid_y_max = np.max(data[:,1])
    grid_z_min = np.min(data[:,2])
    grid_z_max = np.max(data[:,2])

    xs = np.arange(grid_x_min, grid_x_max+1, 1) * dipole_size
    ys = np.arange(grid_y_min, grid_y_max+1, 1) * dipole_size
    zs = np.arange(grid_z_min, grid_z_max+1, 1) * dipole_size
    x,y,z = np.meshgrid(xs, ys, zs)
    x,y,z = np.swapaxes(x, 0,1), np.swapaxes(y, 0,1), np.swapaxes(z, 0,1)

    grid_x_num = int(grid_x_max - grid_x_min + 1)
    grid_y_num = int(grid_y_max - grid_y_min + 1)
    grid_z_num = int(grid_z_max - grid_z_min + 1)

    # Make the direction data for the arrows
    shifted_data = data - [grid_x_min, grid_y_min, grid_z_min, 0, 0, 0, 0]
    forces = np.zeros((grid_x_num, grid_y_num, grid_z_num, 4))
    for dipole in shifted_data:
        forces[int(dipole[0]), int(dipole[1]), int(dipole[2])] = dipole[3:]

    Ftots = forces[:,:,:,0]
    Fxs = forces[:,:,:,1]
    Fys = forces[:,:,:,2]
    Fzs = forces[:,:,:,3]

    return np.array([x, y, z, Ftots, Fxs, Fys, Fzs])


def write_cubes_to_file(cube_positions, cube_lengths, filepath):
    """
    creates an input geometry of indices based on cube positions and sizes.
    Cube positions are integers specifying where the corner in the back, bottom left is (most towards negative in all axes)
    filepath is the path to save the txt file at.
    Raises error if cubes overlapping.
    """

    num_cubes = len(cube_positions)
    file_str = "Nmat="+str(num_cubes)
    record = []
    for cube_i in range(num_cubes):
        cube_pos = cube_positions[cube_i]
        cube_len = cube_lengths[cube_i]

        for i in range(cube_len):
            for j in range(cube_len):
                for k in range(cube_len):
                    x = cube_pos[0]+i
                    y = cube_pos[1]+j
                    z = cube_pos[2]+k

                    # if not (i-cube_len/2)**2 + (j-cube_len/2)**2 + (k-cube_len/2)**2 < (cube_len/2)**2:
                    #     continue


                    if [x,y,z] in record:
                        raise Exception("Cubes overlapping at ",x, y, z)
                    else:
                        record.append([x,y,z])
                        file_str += "\n" + str(x) + " " + str(y) + " " + str(z) + " " + str(cube_i+1)
    
        with open(filepath, 'w') as file:
            file.write(file_str)


def make_cube_indices(cube_position, cube_length):
    record = []
    for i in range(cube_length):
        for j in range(cube_length):
            for k in range(cube_length):
                x = cube_position[0]+i
                y = cube_position[1]+j
                z = cube_position[2]+k

                # if not (i-cube_length/2)**2 + (j-cube_length/2)**2 + (k-cube_length/2)**2 < (cube_length/2)**2:
                #         continue

                record.append([x,y,z])
    return record


def sum_forces(grid, index_list):
    F_net = np.zeros(3)
    for (i,j,k) in index_list:
        F_net += grid[4:,i,j,k]
    return F_net
 

def run_adda(input_file, output_folder, beam_offset, mRes, path_to_adda, Lambda, dpl):
    m_str = ""
    for m in mRes:
        m_str += str(m) + " 0 "

    beam_str = str(beam_offset[0]) + " " + str(beam_offset[1]) + " " + str(beam_offset[2])
        
    adda_command = "./adda -shape read " + input_file + " -lambda " + str(Lambda) + " -dpl " + str(dpl) + " -beam_center " + beam_str + " -m " + m_str + "-store_force -dir " + output_folder
    adda_command = adda_command.split(" ")
    result = subprocess.run(adda_command, cwd=path_to_adda+"/src/seq", stdout=subprocess.DEVNULL)
    # if result != 0:
    #     raise Exception("ADDA error")
    return result
    

def decay_analysis(plot_seps, plot_Fxs):
    # Attempt at peak finding for wavelength oscillation to get the attenuation rate.
    record = []
    finding_max = True
    lastSign = 1
    best = 0
    bestsep = -1
    for i in range(len(plot_seps)):
        if finding_max:
            if plot_Fxs[i] > best:
                best = plot_Fxs[i]
                bestsep = plot_seps[i]

        else:
            if plot_Fxs[i] < best:
                best = plot_Fxs[i]
                bestsep = plot_seps[i]

        if plot_Fxs[i] * lastSign < 0:
            # crosses y=0, so now want to find the other turning point. Record best and reset to 0.
            finding_max ^= True
            record.append((bestsep, best))
            best = 0
            bestsep = i
            lastSign *= -1

    record = np.array(record)
    record = record[record[:,1]!=0]
    c = 0
    plt.plot(record[:,0], record[:,1], label="Turning points")
    record[:,1] = np.abs(record[:,1]+c)
    plt.plot(record[:,0], record[:,1], label="Abs Turning points")
    plt.plot(plot_seps, plot_Fxs , label="Cube length 5 data")
    plt.legend()
    plt.show()



def GetFxsAndSeparations(cube_side_length, Lambda, dipole_size, sep_min, sep_max, path_to_adda, input_file, output_folder, beam_offset, ms, force_per_dipole):
    # Get [plot_seps, plot_Fxs]
    dpl = Lambda/dipole_size
    cube1_pos = [0,0,0]
    plot_seps = np.arange(sep_min, sep_max+1)
    plot_Fxs = []
    for sep in range(sep_min, sep_max+1):
        cube2_pos = [sep,0,0]

        cube_positions = [cube1_pos, cube2_pos]
        cube_lengths = [cube_side_length, cube_side_length]
        
        # Run ADDA.
        write_cubes_to_file(cube_positions, cube_lengths, path_to_adda+"/src/seq/"+input_file)
        run_adda(input_file, output_folder, beam_offset, ms, path_to_adda, Lambda, dpl)
        
        # Extract data from results folder.
        result_folder_path = path_to_adda + "/src/seq/" + output_folder + "/"
        paramDict = PullAddaData_Parameters(result_folder_path + "log")     # e.g. {'lambda': 6.283185307, 'n': [1.5, 0.0], 'dipole_size': 0.418879}    
        dipole_forces = PullAddaData_DipoleForces(result_folder_path + "RadForce-Y") # lists dipoles with 7: x,y,z, F^2, Fx,Fy,Fz
        # dipole_size = paramDict["dipole_size"]

        grid = make_pos_force_grid(dipole_forces, dipole_size)

        cube1_position = cube_positions[0]
        cube1_length = cube_lengths[0]
        cube1_indices = make_cube_indices(cube1_position, cube1_length)
        cube2_position = cube_positions[1]
        cube2_length = cube_lengths[1]
        cube2_indices = make_cube_indices(cube2_position, cube2_length)
        F_net = sum_forces(grid, cube2_indices)

        if force_per_dipole:
            plot_Fxs.append(F_net[0]/cube_side_length**3)
        else:
            plot_Fxs.append(F_net[0])

    return [plot_seps, plot_Fxs]