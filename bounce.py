from vpython import *
from csv import writer

# Create environment
def make_grid(xmax, dx, scale=15):
    for x in range(-xmax, xmax + dx, dx):  # Create vertical lines
        curve(pos=[vector(x / scale, 0, xmax / scale), vector(x / scale, 0, -xmax / scale)], radius=0.001)
    for z in range(-xmax, xmax + dx, dx):  # Create horizontal lines
        curve(pos=[vector(xmax / scale, 0, z / scale), vector(-xmax / scale, 0, z / scale)], radius=0.001)

scene.height = 720
scene.width = 1280
make_grid(xmax=4, dx=1)
lamp = local_light(pos=vector(1,1,1), color=color.white)
color = color.green

# Create mass object
class Mass:
    def __init__(self, m, p, v, a):
        self.m = m
        self.pos = p
        self.vel = v
        self.acc = a

# Create spring object
class Spring:
    def __init__(self, L0, m1, m2, i1, i2, k=100000):
        self.k = k
        self.L0 = L0
        self.axis = m2.pos - m1.pos
        self.L = sqrt(self.axis.x ** 2 + self.axis.y ** 2 + self.axis.z ** 2)
        self.m1 = m1
        self.m2 = m2
        self.pos = m1.pos
        # keep track of indices of associated masses
        self.i1 = i1
        self.i2 = i2

    def update_length(self, m1, m2):
        self.axis = m2.pos - m1.pos
        self.L = sqrt(self.axis.x ** 2 + self.axis.y ** 2 + self.axis.z ** 2)

    def potential(self):
        dL = spring.L0 - spring.L
        return 0.5 * spring.k * dL ** 2

    # Find force spring exerts on particular mass
    def get_force(self, mass):
        F_k = vector(0, 0, 0)
        dL = spring.L - spring.L0
        k_mag = spring.k * (dL)  # magnitude of spring force
        # Compression
        if dL < 0:
            if spring.m1 == mass:
                k_dir = spring.m2.pos - spring.m1.pos  # direction of spring force
            else:
                k_dir = spring.m1.pos - spring.m2.pos
        # Tension
        elif dL > 0:
            if spring.m2 == mass:
                k_dir = spring.m1.pos - spring.m2.pos  # direction of spring force
            else:
                k_dir = - spring.m1.pos + spring.m2.pos
        else:
            k_dir = spring.m2.pos - spring.m1.pos

        k_unit = k_dir / spring.L
        F_k += k_mag * k_unit
        return F_k


# Initialize static cube
off = 0.3 # offset
og_mass_positions = [vector(0, 0 + off, 0), vector(0, 0 + off, 0.1), vector(0, 0.1 + off, 0), vector(0, 0.1 + off, 0.1),
                     vector(0.1, 0 + off, 0), vector(0.1, 0 + off, 0.1), vector(0.1, 0.1 + off, 0.1), vector(0.1, 0.1 + off, 0)]

# Add masses into array
cube_masses = [] # basically an array of Mass objects
for pos in og_mass_positions:
    cube_mass = Mass(m=0.1, p=pos, v=vector(0, 0, 0), a=vector(0, 0, 0))
    cube_masses.append(cube_mass)

# Combine every mass with a spring
def connect_masses(masses):
    cube_springs = []
    for i in range(len(masses)):
        for j in range(i + 1, len(masses)):
            dist_vec = masses[i].pos - masses[j].pos
            dist = dist_vec.mag
            cube_spring = Spring(L0=dist, m1=masses[i], m2=masses[j], i1=i, i2=j)
            cube_springs.append(cube_spring)
    return cube_springs
cube_springs = connect_masses(cube_masses) # an array of Spring objects connecting elements in cube_massses

# Plot masses
def display_masses(mass_array):
    spheres = []
    for mass in mass_array:
        pt = sphere(pos=mass.pos, radius=0.01, color=color)
        spheres.append([pt, mass])
    return spheres # returns a list of Sphere objects that are plotted in vpython
spheres = display_masses(cube_masses)
print(spheres)

# Plot springs
def display_springs(cube_springs):
    lines = []
    for spring in cube_springs:
        cyl = cylinder(pos=spring.m1.pos, axis=spring.axis, color=color, radius=0.001)
        lines.append([cyl, spring])
    return lines # returns a list of Cylinder objects that are plotted in vpython
rods = display_springs(cube_springs)

# Global constants
dt = 0.0001
T = 0
g = vector(0, -9.81, 0)
k_c = 100000
counter = 0

# Clear file for storing energies
with open('energy.csv', 'r+') as file:
    file.truncate(0)

while 1:
    KE = 0
    PE = 0
    E = 0

    print(T)
    rate(1000)
    counter += 1
    mass_positions = [] # initialize array to store this iterations mass positions
    for mass in cube_masses:
        # Update mass parameters based on forces
        F_k = vector(0, 0, 0)
        # Find spring force
        for spring in cube_springs:
            if spring.m1 == mass or spring.m2 == mass:
                F_k += spring.get_force(mass)
                PE += spring.potential()
        F_c = vector(0, 0, 0)
        F_g = mass.m * g  # gravitational force
        F_e = vector(0, 0, 0)  # external forces
        if mass.pos.y < 0:
            F_c = vector(0, -k_c * mass.pos.y, 0)  # restoration force
            PE += abs(0.5 * k_c * mass.pos.y ** 2)
        F = F_g + F_e + F_c + F_k  # sum up forces
        mass.acc = F / mass.m
        mass.vel += mass.acc * dt
        mass.pos += mass.vel * dt
        mass_positions.append(mass.pos)
        T += dt
        KE += 0.5 * mass.m * mass.vel.mag ** 2
        PE += mass.pos.y * mass.m * 9.81

    for i in range(len(rods)):
        rod = rods[i][0]
        spring = rods[i][1]
        spring.update_length(cube_masses[spring.i1], cube_masses[spring.i2])

    if counter % 30 == 0:
        for i in range(len(mass_positions)):
            spheres[i][0].pos = mass_positions[i]

        for i in range(len(rods)):
            rod = rods[i][0]
            spring = rods[i][1]
            rod.pos = mass_positions[spring.i1]
            rod.axis = mass_positions[spring.i2] - mass_positions[spring.i1]

    if counter % 100 == 0:
        E = KE + PE
        with open('energy.csv', 'a', newline='') as file:
            write = writer(file)
            write.writerow([round(T,2), KE, PE, E])
            file.close()

