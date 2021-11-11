from vpython import *

# Create environment
def make_grid(xmax, dx, scale=15):
    for x in range(-xmax, xmax + dx, dx):  # Create vertical lines
        curve(pos=[vector(x / scale, 0, xmax / scale), vector(x / scale, 0, -xmax / scale)], radius=0.001)
    for z in range(-xmax, xmax + dx, dx):  # Create horizontal lines
        curve(pos=[vector(xmax / scale, 0, z / scale), vector(-xmax / scale, 0, z / scale)], radius=0.001)

scene.height = 720
scene.width = 1280
make_grid(xmax=4, dx=1)
color = color.green

# Create mass object
class Mass:
    def __init__(self, m, p, v, a):
        self.m = m
        self.pos = p
        self.vel = v
        self.acc = a

a = 0
b = 150
w = .1
c = 1
# Create spring object
class Spring:
    def __init__(self, L0, m1, m2, i1, i2, k=10000, t=0):
        self.k = a+b*sin(w*t+c)
        self.L0 = L0
        self.axis = m1.pos - m2.pos
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

    def update_k(self, T):
        self.k = a+b*sin(w*T+c)
    # Find force spring exerts on particular mass
    def get_force(self, mass):
        F_k = vector(0, 0, 0)
        if spring.m1 == mass or spring.m2 == mass:
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
off = 0 # offset
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
            cube_spring = Spring(k=10000, L0=dist, m1=masses[i], m2=masses[j], i1=i, i2=j)
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
k_c = 1000000
counter = 0

while 1:
    print(T)
    rate(10000)
    counter += 1
    mass_positions = [] # initialize array to store this iterations mass positions
    for mass in cube_masses:
        # Update mass parameters based on forces
        F_k = vector(0, 0, 0)
        # Find spring force
        for spring in cube_springs:
            F_k += spring.get_force(mass)
        # print("spring force")
        # print(F_k)
        F_c = vector(0, 0, 0)
        F_g = mass.m * g  # gravitational force
        F_e = vector(0, 0, 0)  # external forces
        if mass.pos.y < 0:
            F_c = vector(0, -k_c * mass.pos.y, 0)  # restoration force
        F = F_g + F_e + F_c + F_k  # sum up forces
        print(F)
        mass.acc = F / mass.m
        mass.vel += mass.acc * dt
        mass.pos += mass.vel * dt
        mass_positions.append(mass.pos)
        T += dt

    for i in range(len(mass_positions)):
        spheres[i][0].pos = mass_positions[i]

    for i in range(len(rods)):
        rod = rods[i][0]
        spring = rods[i][1]
        spring.update_length(cube_masses[spring.i1], cube_masses[spring.i2])
        rod.pos = mass_positions[spring.i1]
        print("rod pos")
        print(rod.pos)
        rod.axis = mass_positions[spring.i2] - mass_positions[spring.i1]
        spring.update_k(T)
