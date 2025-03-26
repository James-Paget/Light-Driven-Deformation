import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation


def make_position_animation(data, dt, should_repeat=False, padding=0):
    """
    Animates the positions of the dipoles in time.
    
    data [[[float]]] should be [ [ [x,y,z,...] for each dipole ] for each time].
    dt float is the time step.
    should_repeat boolean determines if the animation should repeat.
    No return.
    """
    
    def update_graph(i):
        current_data = pos_data[i]
        graph._offsets3d = current_data[:,0],current_data[:,1], current_data[:,2]
        ax.set_title("Positions at time={:.3f} s".format(times[i]))

    # Set up data arrays.
    pos_data = np.array(data)[:,:,:3]
    num_times = len(data)
    times = np.arange(0,dt*num_times, dt)

    # Set graph bounds.
    x_min = np.min(pos_data[:,:,0])
    x_max = np.max(pos_data[:,:,0])
    y_min = np.min(pos_data[:,:,1])
    y_max = np.max(pos_data[:,:,1])
    z_min = np.min(pos_data[:,:,2])
    z_max = np.max(pos_data[:,:,2])

    # Initialise graph.
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set(xlabel="x", ylabel="y", zlabel="z", xlim=(x_min-padding,x_max+padding), ylim=(y_min-padding,y_max+padding), zlim=(z_min,z_max), title="Positions at time=0.000 s")

    initial_pos_data = pos_data[0]
    graph = ax.scatter(initial_pos_data[:,0], initial_pos_data[:,1], initial_pos_data[:,2])

    # Animate and plot.
    anim = matplotlib.animation.FuncAnimation(fig, update_graph, num_times, interval=40, blit=False, repeat=should_repeat)
    # anim.save("out.mp4", fps=40)
    plt.show()


## QUIVER

def make_quiver_plot(data, x_spacing, y_spacing, z_spacing, should_normalise_arrows=True):
    """
    Plots the positions and forces of a set of dipoles.
    
    data [[float]] should be [ [x,y,z,|F|^2, Fx, Fy, Fz] for each dipole ].
    No return.
    """
    ax = plt.figure().add_subplot(projection='3d')
    data = np.array(data)

    # Convert positions to an integer grid to avoid floats
    data[:,:3] = np.round(data[:,:3]/[x_spacing, y_spacing, z_spacing])

    grid_x_min = np.min(data[:,0])
    grid_x_max = np.max(data[:,0])
    grid_y_min = np.min(data[:,1])
    grid_y_max = np.max(data[:,1])
    grid_z_min = np.min(data[:,2])
    grid_z_max = np.max(data[:,2])

    xs = np.arange(grid_x_min, grid_x_max+1, 1) * x_spacing
    ys = np.arange(grid_y_min, grid_y_max+1, 1) * y_spacing
    zs = np.arange(grid_z_min, grid_z_max+1, 1) * z_spacing
    x,y,z = np.meshgrid(xs, ys, zs)

    grid_x_num = int(grid_x_max - grid_x_min + 1)
    grid_y_num = int(grid_y_max - grid_y_min + 1)
    grid_z_num = int(grid_z_max - grid_z_min + 1)

    # Make the direction data for the arrows
    shifted_data = data - [grid_x_min, grid_y_min, grid_z_min, 0, 0, 0, 0]
    forces = np.zeros((grid_x_num, grid_y_num, grid_z_num, 4))
    for dipole in shifted_data:
        forces[int(dipole[0]), int(dipole[1]), int(dipole[2])] = dipole[3:]

    Fxs = forces[:,:,:,1]
    Fys = forces[:,:,:,2]
    Fzs = forces[:,:,:,3]

    # Make colours
    F_max = np.max(forces[:,:,:,0])
    Fs = forces[:,:,:,0] / F_max
    cols = cm.viridis(Fs[:,:]).reshape(-1,4)
    
    ax.set(xlabel="x", ylabel="y", zlabel="z", xlim=(xs[0],xs[-1]), ylim=(ys[0],ys[-1]), zlim=(zs[0],zs[-1]), title="Quiver plot")
    ax.quiver(x, y, z, Fxs, Fys, Fzs, color=cols, length=0.1, normalize=should_normalise_arrows)

    plt.show()


# # ANIMATION SETUP
# my_num_times = 50
# my_num_particles = 500
# my_dt = 0.13
# my_data = np.random.rand(my_num_times, my_num_particles, 6)*10 - 3
# make_position_animation(my_data, my_dt, should_repeat=True)



# # QUIVER SETUP
# my_data_quiv = np.random.random((33,7)) # x,y,z, |F|^2, Fx, Fy, Fz
# my_data_quiv[:,3] = my_data_quiv[:,4]**2+my_data_quiv[:,5]**2+my_data_quiv[:,6]**2

# num = 0
# for i in range(1,4):
#     for j in range(1,4):
#         for k in range(1,4):
#             my_data_quiv[num,0] = i
#             my_data_quiv[num,1] = j
#             my_data_quiv[num,2] = k
#             num +=1

# my_data_quiv[num,0] = 0
# my_data_quiv[num,1] = 2
# my_data_quiv[num,2] = 2
# num +=1
# my_data_quiv[num,0] = 2
# my_data_quiv[num,1] = 2
# my_data_quiv[num,2] = 0
# num +=1
# my_data_quiv[num,0] = 2
# my_data_quiv[num,1] = 0
# my_data_quiv[num,2] = 2
# num +=1
# my_data_quiv[num,0] = 4
# my_data_quiv[num,1] = 2
# my_data_quiv[num,2] = 2
# num +=1
# my_data_quiv[num,0] = 2
# my_data_quiv[num,1] = 2
# my_data_quiv[num,2] = 4
# num +=1
# my_data_quiv[num,0] = 2
# my_data_quiv[num,1] = 4
# my_data_quiv[num,2] = 2
# num +=1

# # print(my_data_quiv[:,:3])

# my_spacing = 0.25
# my_data_quiv[:,:3] *= my_spacing

# make_quiver_plot(my_data_quiv, my_spacing, my_spacing, my_spacing, should_normalise_arrows=True)