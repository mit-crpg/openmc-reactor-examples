#Importing the Necessary Libraries
import matplotlib.pyplot as plt
from math import pi, sin, cos, sqrt
import numpy as np
import openmc

"""Materials"""

#2.0% Enriched Fuel
fuel = openmc.Material(name='2.0% Fuel')
fuel.add_element('U', 1.0, enrichment = 2.0)
fuel.add_nuclide('O16', 2.0)
fuel.set_density('g/cm3', 10.400)

#Light water
water = openmc.Material(name='Water')
water.add_nuclide('H1', 2.0)
water.add_nuclide('O16', 1.0)
water.add_s_alpha_beta('c_H_in_H2O')
water.set_density('g/cm3', 1.0)

#Zircaloy
zircaloy = openmc.Material(name='Zircaloy')
zircaloy.add_element('Zr', 0.99)
zircaloy.add_element('Nb', 0.1)
zircaloy.set_density('g/cm3', 8.59)

#Carrier Rod material
wall = openmc.Material(name='Carrier Rod Material')
wall.add_element('Zr', 0.975)
wall.add_element('Nb', 0.025)
wall.set_density('g/cm3', 8.59)

#Helium
helium = openmc.Material(name = 'Helium')
helium.add_element('He', 1)
helium.set_density('g/cm3', 0.178)

#Instantiate a Materials collection
materials_file = openmc.Materials([fuel, water, zircaloy, helium, wall])

#Export to "materials.xml"
materials_file.export_to_xml()


"""Creating the Fuel Rod"""

#Making the boundary planes for the sleeves
sleeve_min_z = openmc.ZPlane(z0=-12, boundary_type='transmission')
sleeve_max_z = openmc.ZPlane(z0=+12, boundary_type='transmission')

# Create cylinders for the fuel and clad
fuel_inner_radius = openmc.ZCylinder(x0=0, y0=0, r=0.55) #0.3 mm clearing between fuel pellets and tube
fuel_outer_radius = openmc.ZCylinder(x0=0, y0=0, r=0.566)    #68) #13.6 mm diameter
clad_inner_radius = openmc.ZCylinder(x0=0, y0=0, r=0.596)       #0.71) #0.3 mm clearing between fuel an cladding
clad_outer_radius = openmc.ZCylinder(x0=0, y0=0, r=0.68)    #796) #0.86 mm wall thickness

# Create a Universe to encapsulate a fuel pin
pin_cell_universe = openmc.Universe(name='2.0% Fuel Pin') 

# Create fuel Cell
fuel_cell = openmc.Cell(name='2.0% Fuel')
fuel_cell.fill = fuel
fuel_cell.region = -fuel_inner_radius
pin_cell_universe.add_cell(fuel_cell)

# Create Void Space
void_space = openmc.Cell(name='empty_space')
void_space.fill = helium
void_space.region = +fuel_inner_radius & -fuel_outer_radius
pin_cell_universe.add_cell(void_space)

# Create a clad Cell
clad_cell = openmc.Cell(name='2.0% Clad')
clad_cell.fill = zircaloy
clad_cell.region = +clad_inner_radius & -clad_outer_radius
pin_cell_universe.add_cell(clad_cell)

# Create a 'sleeve' Cell that slices through the middle of the elongated element
sleeve_cell = openmc.Cell(name='Sleeve')
sleeve_cell.fill = helium
sleeve_cell.region = -clad_outer_radius & +sleeve_min_z & -sleeve_max_z
pin_cell_universe.add_cell(sleeve_cell)

# Create a moderator Cell
moderator_cell = openmc.Cell(name='Moderator')
moderator_cell.fill = water
moderator_cell.region = +clad_outer_radius
pin_cell_universe.add_cell(moderator_cell)


"""Creating the Carrier Rod"""

#Create Cylinders
rod_inner_radius = openmc.ZCylinder(x0=0, y0=0, r=0.75) 
rod_outer_radius = openmc.ZCylinder(x0=0, y0=0, r=0.7625)

#Create Universe
rod_cell_universe = openmc.Universe(name='Carrier Rod')

# Inner Rod
inner_rod = openmc.Cell(name='Inner Rod')
inner_rod.fill = helium
inner_rod.region = -rod_inner_radius
rod_cell_universe.add_cell(inner_rod)

# Outer Rod
outer_rod = openmc.Cell(name='Outer Rod')
outer_rod.fill = wall
outer_rod.region = +rod_inner_radius & -rod_outer_radius
rod_cell_universe.add_cell(outer_rod)


"""CircLattice of the Bundle"""

#Defining boundary planes
element_min_z = openmc.ZPlane(z0=-364, boundary_type='reflective')
element_max_z = openmc.ZPlane(z0=+364, boundary_type='reflective')
carrier_min_z = openmc.ZPlane(z0=-600, boundary_type='reflective')
carrier_max_z = openmc.ZPlane(z0=+4000, boundary_type='reflective')

#Creating circular lattice
circlat = openmc.Universe(name='Circular Lattice')

#Making the channel and cell
channel = openmc.ZCylinder(r=4.0, boundary_type='reflective')

channel_cell = openmc.Cell(name='Channel')
channel_cell.fill = water
channel_cell.region = -channel & +rod_outer_radius & +element_min_z & -element_max_z

#Calculations before finding pin placement
r_channel = 4.0
r_element = 0.796
r_outer = r_channel-r_element

padding = (pi/12)*r_outer-r_element
r_inner = r_outer*cos(pi/12)-sqrt((r_element*2+padding)**2-(r_element+padding)**2)

#Adding the fuel pins into the channel
num_in_rings = [6, 12]
radii = [r_inner, r_outer]
angles = [0, pi/12]

for index in range(len(num_in_rings)):
	elem_num = num_in_rings[index]
	ring_radius = radii[index]
	theta_0 = angles[index]
	theta = 2*pi/elem_num
	for element in range(elem_num):
		x = ring_radius*cos(element*theta + theta_0)
		y = ring_radius*sin(element*theta + theta_0)

		pin_boundary = openmc.ZCylinder(x0=x, y0=y, r=r_element)
		pin = openmc.Cell(fill=pin_cell_universe, region=-pin_boundary & +element_min_z & -element_max_z)
		pin.translation = (x,y,0)
		pin.id = (index+1)*100 + element
		channel_cell.region &= ~pin.region
		circlat.add_cell(pin)

#Adding the 3 Components of the Carrier Rod
rod_bound = openmc.ZCylinder(x0=0, y0=0, r=0.7625, boundary_type='reflective')

rod_above = openmc.Cell(fill=rod_cell_universe, region=-rod_bound & -carrier_max_z & +element_max_z)
channel_cell.region &= ~rod_above.region
rod_mid = openmc.Cell(fill=rod_cell_universe, region=-rod_outer_radius & +element_min_z & -element_max_z)
channel_cell.region &= ~rod_mid.region
rod_below = openmc.Cell(fill=rod_cell_universe, region=-rod_bound & -element_min_z & +carrier_min_z)
channel_cell.region &= ~rod_below.region

circlat.add_cell(rod_above)
circlat.add_cell(rod_mid)
circlat.add_cell(rod_below)
circlat.add_cell(channel_cell)

# Create root Universe
geometry = openmc.Geometry(circlat)
geometry.export_to_xml()

#plt.figure(figsize=(12,12))
#circlat.plot(basis='yz', origin=(0,0,0), width=(10,1000), color_by='material', pixels=[1000,1000])
#plt.savefig('Vertical')
#plt.savefig('RBMK')

# OpenMC simulation parameters
batches = 100
inactive = 10
particles = 5000

# Instantiate a Settings object
settings_file = openmc.Settings()
settings_file.batches = batches
settings_file.inactive = inactive
settings_file.particles = particles

# Create an initial uniform spatial source distribution over fissionable zones
bounds = [-4, -4, -3.4, 4, 4, 3.4]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings_file.source = openmc.Source(space=uniform_dist)

# Export to "settings.xml"
settings_file.export_to_xml()

openmc.run()
