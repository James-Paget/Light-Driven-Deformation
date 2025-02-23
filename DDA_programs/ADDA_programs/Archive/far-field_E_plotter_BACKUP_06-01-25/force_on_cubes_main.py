import numpy as np
import subprocess
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from file_reading_funcs import PullAddaData_Parameters, PullAddaData_DipoleForces
from plot_funcs import make_quiver2_plot, make_quiver3_plot
from force_on_cubes_funcs import make_pos_force_grid, make_cube_indices, sum_forces, write_cubes_to_file, run_adda, decay_analysis, GetFxsAndSeparations





def force_with_sep_main(cube_side_lengths, sep_min, sep_max, Lambda, dipole_size, force_per_dipole=False):
    path_to_adda = "/Users/david/adda"
    input_file = "nmat_2_cubes.txt"
    output_folder = "output_folder"
    beam_offset = [0,0,0]
    m = 1.5
    ms = [m,m]

    dpl = Lambda/dipole_size

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    for cube_side_length in cube_side_lengths:
        if sep_min < cube_side_length:
            sep_min = cube_side_length
            print("force_with_sep_main: set sep min to", cube_side_length)


        plot_seps, plot_Fxs = GetFxsAndSeparations(cube_side_length, Lambda, dipole_size, sep_min, sep_max, path_to_adda, input_file, output_folder, beam_offset, ms, force_per_dipole)
        ax.plot(plot_seps, plot_Fxs, label=f"Cube number of dipoles: ${cube_side_length}^3$")

    
    if force_per_dipole:
        ax.set(xlabel="Separation in lattice spacings", ylabel="F_x", title="F_x per dipole against separation")
    else:
        ax.set(xlabel="Separation in lattice spacings", ylabel="F_x", title="F_x against separation")
    plt.legend()
    plt.show()


def force_with_sep_const_cube_size(cube_side_lengths, sep_min, sep_max, Lambda, cube_size, force_per_dipole=False):
    path_to_adda = "/Users/david/adda"
    input_file = "nmat_2_cubes.txt"
    output_folder = "output_folder"
    beam_offset = [0,0,0]
    m = 1.5
    ms = [m,m]

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    for cube_side_length in cube_side_lengths:
        # scale dipole size down so cube stays the same size.
        dipole_size = cube_size/cube_side_length
        dpl = Lambda/dipole_size
        if sep_min < cube_side_length:
            sep_min = cube_side_length
            print("force_with_sep_main: set sep min to", cube_side_length)

        plot_seps, plot_Fxs = GetFxsAndSeparations(cube_side_length, Lambda, dipole_size, sep_min, sep_max, path_to_adda, input_file, output_folder, beam_offset, ms, force_per_dipole)
        
        # cube_side_length increases, dipole size decreases, so lattice separations need to be scaled by dipole size.
        plot_seps = plot_seps * dipole_size

        # Option to truncate early
        trunc_num = sep_max/max(cube_side_lengths)
        plot_seps = np.array(plot_seps)
        plot_Fxs = np.array(plot_Fxs)
        filter = plot_seps <= trunc_num
        plot_seps = plot_seps[filter]
        plot_Fxs = plot_Fxs[filter]

        ax.plot(plot_seps, plot_Fxs, label=f"Cube number of dipoles: ${cube_side_length}^3$")

    if force_per_dipole:
        ax.set(xlabel="Separation /$10^{-6}$m", ylabel="$F_x$ /$10^{-13}$N", title="$F_x$ per dipole against cube separation")
    else:
        ax.set(xlabel="Separation /$10^{-6}$m", ylabel="$F_x$ /$10^{-13}$N", title="$F_x$ against cube separation")
    plt.legend()

    plt.savefig('F_vs_sep.eps', format='eps')
    plt.show()


def force_with_wavelength(cube_side_length, sep_min, sep_max, Lambdas, dipole_size, force_per_dipole=False):
    path_to_adda = "/Users/david/adda"
    input_file = "nmat_2_cubes.txt"
    output_folder = "output_folder"
    beam_offset = [0,0,0]
    m = 1.5
    ms = [m,m]

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    for Lambda in Lambdas:
        if sep_min < cube_side_length:
            sep_min = cube_side_length
            print("force_with_sep_main: set sep min to", cube_side_length)

        plot_seps, plot_Fxs = GetFxsAndSeparations(cube_side_length, Lambda, dipole_size, sep_min, sep_max, path_to_adda, input_file, output_folder, beam_offset, ms, force_per_dipole)
        ax.plot(plot_seps, plot_Fxs, label=f"Wavelength {Lambda}e-6m")
    
    if force_per_dipole:
        ax.set(xlabel="Separation in lattice spacings", ylabel="$F_x$ /$10^{-13}$N", title="$F_x$ per dipole against cube separation")
    else:
        ax.set(xlabel="Separation in lattice spacings", ylabel="$F_x$ /$10^{-13}$N", title="$F_x$ against cube separation")
    plt.legend()

    plt.savefig('F_vs_sep.eps', format='eps')
    plt.show()


    


        




def main():
    path_to_adda = "/Users/david/adda"
    input_file = "nmat_2_cubes.txt"
    output_folder = "output_folder"
    beam_offset = [0,0,0]
    m = 1.5
    ms = [m,m]
    Lambda = 6.28 # e-6 m
    dipole_size = Lambda/15#0.418
    dpl = Lambda/dipole_size
    cube_side_length = 5
    z_plane_i = 0
    cube_x_sep = 10

    # Cube positions are integers specifying where the corner in the back, bottom left is (most towards negative in all axes)
    cube1_pos = [0,0,0]
    cube2_pos = [cube_x_sep,0,0]

    cube_positions = [cube1_pos, cube2_pos]
    cube_lengths = [cube_side_length, cube_side_length]
    
    # Run ADDA.
    write_cubes_to_file(cube_positions, cube_lengths, path_to_adda+"/src/seq/"+input_file)
    run_adda(input_file, output_folder, beam_offset, ms, path_to_adda, Lambda, dpl)
    
        
    # Extract data from results folder.
    result_folder_path = path_to_adda + "/src/seq/" + output_folder + "/"
    paramDict = PullAddaData_Parameters(result_folder_path + "log")     # e.g. {'lambda': 6.283185307, 'n': [1.5, 0.0], 'dipole_size': 0.418879}    
    dipole_forces = PullAddaData_DipoleForces(result_folder_path + "RadForce-Y") # lists dipoles with 7: x,y,z, F^2, Fx,Fy,Fz
    dipole_size = paramDict["dipole_size"]
    grid = make_pos_force_grid(dipole_forces, dipole_size)

    for cube_i in range(len(cube_positions)):
        cube_position = cube_positions[cube_i]
        cube_length = cube_lengths[cube_i]
        cube_indices = make_cube_indices(cube_position, cube_length)
        F_net = sum_forces(grid, cube_indices)


    # make_quiver2_plot(dipole_forces, z_plane_i, dipole_size, cube_side_length)
    # make_quiver3_plot(grid)


    Lambda = 6.28 # e-6 m
    dipole_size = 0.418
    cube_side_lengths = [1,2,3,4] # in lattice spacings
    # force_with_sep_main(cube_side_lengths, 0, 100, Lambda, dipole_size, force_per_dipole=True)

    force_with_sep_const_cube_size(cube_side_lengths, 0, 200, Lambda, cube_size=1, force_per_dipole=False)

    Lambdas = [i for i in range(1,10,2)]
    cube_side_length = 5 # in lattice spacings
    dipole_size = 0.5
    # force_with_wavelength(cube_side_length, 0, 100, Lambdas, dipole_size, force_per_dipole=False)


main()
