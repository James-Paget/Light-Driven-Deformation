"""
Particle creation routines
"""
import numpy as np

class ParticleCollection (object):
    num_particles = 0
    particle_spec = {'Silicon': [3.9, '#e0a105', '#ffe10a','S'],
                     'Sapphire': [2.5, '#b01240', '#e81e80','O'],
                     'Diamond': [2.417, '#b01240', '#e81e80','C'],
                     'SiN': [2.046, '#d8ff60', '#e8ff80','P',],
                     'Glass': [1.5, '#1870d0', '#2298f5','N'],
                     'LeadGlass': [1.6, '#1870d0', '#2298f5','N'],
                     'Low': [1.446, '#7e7e7e', '#aeaeae','H',],
                     'FusedSilica': [1.458+0.0j, '#7e7e7e', '#aeaeae','B',],
                     'FusedSilica0001': [1.458+0.0001j, '#7e7e7e', '#aeaeae','B'],
                     'FusedSilica0002': [1.458+0.0002j, '#7e7e7e', '#aeaeae','B'],
                     'FusedSilica0005': [1.458+0.0005j, '#7e7e7e', '#aeaeae','B'],
                     'FusedSilica001': [1.458+0.001j, '#7e7e7e', '#aeaeae','B'],
                     'FusedSilica002': [1.458+0.002j, '#7e7e7e', '#aeaeae','B'],
                     'FusedSilica005': [1.458+0.005j, '#7e7e7e', '#aeaeae','B'],
                     'FusedSilica01': [1.458+0.01j, '#7e7e7e', '#aeaeae','B'],
                     'Chromium': [3.41 + 3.57j, '#ffffff', '#eeeeee', 'Cr'],
                     'Air': [1.00, '#0f0f0f', '#1f1f1f', 'Ar'],
                     }
                     
    defaults = {"default_material":'FusedSilica',
                "default_radius":'200e-9',
                "default_density":'2200', # kg/m^3 typical glass value
                }
    default_material = defaults['default_material']
    default_radius = float(defaults['default_radius'])
    default_density = float(defaults['default_density'])

    
    def __init__(self,particleinfo):
        if particleinfo==None:
            # Set defaults
            ParticleCollection.num_particles = 1
            self.particle_type = np.asarray([ParticleCollection.default_material])
            self.particle_radius = np.asarray([ParticleCollection.default_radius])
            self.particle_indices = np.asarray([ParticleCollection.particle_spec[ParticleCollection.default_material][0]],dtype=complex)
            self.particle_colour = np.asarray([ParticleCollection.particle_spec[ParticleCollection.default_material][1]])
            self.particle_vtfcolour = np.asarray([ParticleCollection.particle_spec[ParticleCollection.default_material][3]])
            self.particle_density = np.asarray([ParticleCollection.default_density])
            self.particle_positions = np.zeros((1,3),dtype=float)
        else:
            # Read from file
            self.default_material = particleinfo.get('default_material',ParticleCollection.default_material)
            self.default_radius = float(particleinfo.get('default_radius',ParticleCollection.default_radius))
            self.particle_list = particleinfo.get('particle_list',None)
            if self.particle_list==None or self.particle_list==False:
                # Set defaults
                ParticleCollection.num_particles = 1
                self.particle_type = np.asarray([ParticleCollection.default_material])
                self.particle_radius = np.asarray([ParticleCollection.default_radius])
                self.particle_indices = np.asarray([ParticleCollection.particle_spec[ParticleCollection.default_material][0]],dtype=complex)
                self.particle_colour = np.asarray([ParticleCollection.particle_spec[ParticleCollection.default_material][1]])
                self.particle_vtfcolour = np.asarray([ParticleCollection.particle_spec[ParticleCollection.default_material][3]])
                self.particle_density = np.asarray([ParticleCollection.default_density])
                self.particle_positions = np.zeros((1,3),dtype=float)
            else:
                # Read individual particles
                i=0
                self.particle_type = []
                self.particle_radius = []
                self.particle_colour = []
                self.particle_vtfcolour = []
                self.particle_positions = []
                self.particle_indices = []
                self.particle_density = []
                for newparticle in self.particle_list:
                    particle = self.particle_list[newparticle]
                    print("Loading particle",particle)
                    if particle != None:
                        self.particle_type.append(particle.get('material',self.default_material))
                        self.particle_radius.append(float(particle.get('radius',self.default_radius)))
                        self.altcolour = bool(particle.get('altcolour',False))
                        if self.altcolour==False:
                            self.particle_colour.append(ParticleCollection.particle_spec[self.particle_type[i]][1])
                        else:
                            self.particle_colour.append(ParticleCollection.particle_spec[self.particle_type[i]][2])
                        self.particle_vtfcolour.append(ParticleCollection.particle_spec[self.particle_type[i]][3])
                        self.particle_indices.append(ParticleCollection.particle_spec[self.particle_type[i]][0])
                        self.coords = particle.get('coords',"0.0 0.0 0.0")
                        self.fields = self.coords.split(" ")
                        if self.fields[0]=="None":
                            self.particle_positions.append(np.array((0.0,0.0,0.0),dtype=np.float64))
                        else:
                            self.particle_positions.append(np.array((0.0,0.0,0.0),dtype=np.float64))
                            for j in range(min(len(self.fields),3)):
                                self.particle_positions[i][j] = float(self.fields[j])
                        self.particle_density.append(ParticleCollection.default_density)
                    else:
                        ParticleCollection.num_particles = 1
                        self.particle_type.append(ParticleCollection.default_material)
                        self.particle_radius.append(ParticleCollection.default_radius)
                        self.particle_indices.append(ParticleCollection.particle_spec[ParticleCollection.default_material][0])
                        self.particle_colour.append(ParticleCollection.particle_spec[ParticleCollection.default_material][1])
                        self.particle_vtfcolour.append(ParticleCollection.particle_spec[ParticleCollection.default_material][3])
                        self.particle_density.append(ParticleCollection.default_density)
                        self.particle_positions.append(np.array((0.0,0.0,0.0),dtype=np.float64))
                    i+=1
                ParticleCollection.num_particles = i

    def get_refractive_indices(self):
        return np.asarray(self.particle_indices,dtype=complex)
        
    def get_particle_types(self):
        return np.asarray(self.particle_type)
        
    def get_particle_colours(self):
        return np.asarray(self.particle_colour)
        
    def get_particle_vtfcolours(self):
        return np.asarray(self.particle_vtfcolour)
        
    def get_particle_radii(self):
        return np.asarray(self.particle_radius,dtype=float)

    def get_particle_density(self):
        return np.asarray(self.particle_density,dtype=float)

    def get_particle_positions(self):
        return np.asarray(self.particle_positions,dtype=float).reshape((ParticleCollection.num_particles,3))
    
