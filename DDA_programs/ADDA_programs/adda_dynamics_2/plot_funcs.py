import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
import matplotlib.animation
from force_on_cubes_funcs import make_pos_force_grid


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


def make_quiver2_plot(dipole_forces, z_plane_i, dipole_size, cube_side_length):
     # require even side length and int+0.5 plane OR odd side length and int plane.
    if not (cube_side_length + 2*z_plane_i)%2 == 1:
        raise Exception("Invalide z_plane_i")
        
    # Filters dipole_forces to return data (np array) with the same z value. This is specified by z_plane_i which is the index in ADDA's geometry input e.g. 1, 2, 0.5...
    dipole_forces = np.array(dipole_forces)
    z_plane_pos = z_plane_i * dipole_size
    filter = (abs(dipole_forces[:,2] - z_plane_pos) < dipole_size/10)
    data = dipole_forces[filter]

    # pull out the x,y coords in integer multiples of the dipole size.
    data[:,:2] = np.round(data[:,:2]/[dipole_size, dipole_size])

    # get min and max to make x,y box
    grid_x_min = np.min(data[:,0])
    grid_x_max = np.max(data[:,0])
    grid_y_min = np.min(data[:,1])
    grid_y_max = np.max(data[:,1])
    xs = np.arange(grid_x_min, grid_x_max+1, 1) * dipole_size
    ys = np.arange(grid_y_min, grid_y_max+1, 1) * dipole_size

    x,y = np.meshgrid(xs, ys)

    grid_x_num = len(xs)
    grid_y_num = len(ys)

    # shift by x and y mins so indices start at 0. Then for those index pairs, assign the force.
    shifted_data = data - [grid_x_min, grid_y_min, 0, 0, 0, 0, 0]
    forces = np.zeros((grid_x_num, grid_y_num, 4))
    for dipole in shifted_data:
        forces[int(dipole[0]), int(dipole[1])] = dipole[3:]

    Fxs = forces[:,:,1]
    Fys = forces[:,:,2]

    # Make colours
    F_max = np.max(forces[:,:,0])
    Fs = forces[:,:,0] / F_max
    cols = cm.plasma(Fs[:,:]).reshape(-1,4)
    
    _, ax = plt.subplots()
    padding = 2
    ax.set(xlabel="x", ylabel="y", xlim=(xs[0]-padding,xs[-1]+padding), ylim=(ys[0]-padding,ys[-1]+padding))
    ax.quiver(x, y, np.transpose(Fxs), np.transpose(Fys), color=cols)
    ax.set_aspect("equal")

    plt.show()


## QUIVER



def make_quiver3_plot(grid, should_normalise_arrows=True):
    """
    Plots the positions and forces of a set of dipoles.
    
    grid is [x, y, z, Ftots, Fxs, Fys, Fzs] where each has a value for each point in the 3d box.
    No return.
    """
    Fscale = 1e2
    [x, y, z, Ftots, Fxs, Fys, Fzs] = grid
    ax = plt.figure().add_subplot(projection='3d')
    

    # Make colours
    Fcols = Ftots / np.max(Ftots)
    cols = cm.plasma(Fcols[:,:]).reshape(-1,4)
    
    ax.set(xlabel="x /μm", ylabel="y /μm", zlabel="z /μm", xlim=(x[0,0,0],x[-1,1,0]), ylim=(y[0,0,0],y[0,-1,0]), zlim=(z[0,0,0],z[0,0,-1]))
    # ax.quiver(x, y, z, Fxs, Fys, Fzs, color=cols, length=0.25, normalize=should_normalise_arrows)
    ax.quiver(x, y, z, Fxs*Fscale, Fys*Fscale, Fzs*Fscale, color=cols, normalize=False)
    
    ax.set_aspect("equal")
    plt.savefig('forces_quiver.eps', format='eps')
    plt.show()








def plot_voxelised_sphere(max_rad):
    # max_rad = sphere size

    def make_cube_surface(radius, centre):
        vertices = np.array([[0, 0, 0],  [1, 0, 0],[1, 1, 0],[0, 1, 0],[0, 0, 1],  [1, 0, 1],[1, 1, 1],[0, 1, 1]]) * radius + centre

        edges = [
        [vertices[0], vertices[1], vertices[5], vertices[4]], 
        [vertices[7], vertices[6], vertices[2], vertices[3]],  
        [vertices[0], vertices[1], vertices[2], vertices[3]], 
        [vertices[4], vertices[5], vertices[6], vertices[7]],
        [vertices[0], vertices[3], vertices[7], vertices[4]],
        [vertices[1], vertices[2], vertices[6], vertices[5]]  
        ]
        return edges


    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])

    for x in range(-max_rad, max_rad+1):
        for y in range(-max_rad, max_rad+1):
            for z in range(-max_rad, max_rad+1):
                r2 = x**2 + y**2 + z**2
                if r2 <= max_rad**2 and r2 > 0.5*max_rad**2:
                    edges = make_cube_surface(1, np.array([x,y,z]))
                    ax.add_collection3d(Poly3DCollection(edges, linewidths=0.5, edgecolors="gray", alpha=1))

    plt.savefig('voxels.eps', format='eps')
    plt.show()

plot_voxelised_sphere(6)