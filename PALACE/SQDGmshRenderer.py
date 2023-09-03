from qiskit_metal.qgeometries.qgeometries_handler import QGeometryTables
from qiskit_metal.renderers.renderer_mpl.mpl_renderer import QMplRenderer
import gmsh
import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import shapely
import qiskit_metal as metal
from qiskit_metal import designs, draw
from qiskit_metal import MetalGUI, Dict, open_docs
from qiskit_metal import qgeometries
from qiskit_metal.toolbox_metal import math_and_overrides
from qiskit_metal.qlibrary.core import QComponent
from collections import OrderedDict


class Palace_Gmsh_Renderer:

    #lists to store entites which comprise the dielectric gaps and metals in the design
    dielectric_gap_surface = []
    metal_surfaces = []

    #class variable to store ID (int) representing the chip base
    dielectric_box = None
    
    #fragmented dielectric vol
    frag_dielectric_vol = None

    #list to store launch pads to create ports on
    launch_pads_list = []

    #lumped element ports
    ports_list = []

    #meshing parameter
    lc = None

    #constructor takes qiskit-metal design
    def __init__(self, design):
        self.design = design

        #get dimensions of chip base and convert to design units in 'mm'
        self.center_x = self.design.parse_value(self.design.chips['main'].size.center_x)
        self.center_y = self.design.parse_value(self.design.chips['main'].size.center_y)
        self.center_z = self.design.parse_value(self.design.chips['main'].size.center_z)
        self.size_x = self.design.parse_value(self.design.chips['main'].size.size_x)
        self.size_y = self.design.parse_value(self.design.chips['main'].size.size_y)
        self.size_z = self.design.parse_value(self.design.chips['main'].size.size_z)



    def convert_design_to_gmsh(self):
        '''convert qiskit metal design to geometry in Gmsh'''

        #Start gmsh and add model
        gmsh.initialize()
        gmsh.model.add('qiskit_to_gmsh')

        #create list with component names
        component_list = self.design.all_component_names_id()
        component_names = []
        for i in range(len(component_list)):
            component_names.append(component_list[i][0])

        #draw all components into Gmsh
        for component_name in component_names:
            if(not self.design.qgeometry.get_component(component_name)['path'].empty):
                path_to_draw = self.design.qgeometry.get_component(component_name)['path']
                print('drawing path:', component_name)
                for index in path_to_draw.index:
                    self.draw_path(path_to_draw.loc[index], component_name, index, flag = 'path')
            if(not self.design.qgeometry.get_component(component_name)['poly'].empty):
                poly_to_draw = self.design.qgeometry.get_component(component_name)['poly']
                print('drawing polygon:', component_name)
                for index in poly_to_draw.index:
                    self.draw_polygon(poly_to_draw.loc[index], component_name, index, flag = 'poly')
            if(not self.design.qgeometry.get_component(component_name)['junction'].empty):
                junc_to_draw = self.design.qgeometry.get_component(component_name)['junction']
                print('drawing junction:', component_name)
                for index in junc_to_draw.index:
                    self.draw_path(junc_to_draw.loc[index], component_name, index, flag = 'junction')

        #add in physical groups for design in gmsh
        dielectric_gaps = []
        for i,value in enumerate(Palace_Gmsh_Renderer.dielectric_gap_surface):
            dielectric_gaps.append(value[1])

        #create lumped element ports
        for i,value in enumerate(Palace_Gmsh_Renderer.launch_pads_list):
            self.create_ports_on_launchpad(value)

        self.add_ground_plane() 
        self.draw_air_box()

        #update geometries
        gmsh.model.geo.synchronize()
        gmsh.model.occ.synchronize()

    

    def draw_chip_base(self):
        '''This method draws the chip base given the dimensions defined by the user'''

        #half values of the sizes
        half_size_x = self.size_x/2
        half_size_y = self.size_y/2

        #store coordinates for the top surface of the chip
        surface_1 = {'point1': [self.center_x + half_size_x, self.center_y + half_size_y, self.center_z],
                     'point2': [self.center_x + half_size_x, self.center_y - half_size_y, self.center_z],
                     'point3': [self.center_x - half_size_x, self.center_y - half_size_y, self.center_z],
                     'point4': [self.center_x - half_size_x, self.center_y + half_size_y, self.center_z]}
        
        #define lists to store points, lines and surfaces
        points = []
        lines = []
        surfaces = []

        #add points for chip base
        for i,value in enumerate(surface_1):
            point = gmsh.model.occ.addPoint(surface_1[value][0], surface_1[value][1], surface_1[value][2])
            points.append(point)

        #draw lines for chip base
        for j,value in enumerate(points):
            if(j<len(points)-1):
                line = gmsh.model.occ.add_line(points[j], points[j+1])
                lines.append(line)        
        line = gmsh.model.occ.add_line(points[len(points)-1], points[0])
        lines.append(line)        

        #create curved loop
        curve_loop = gmsh.model.occ.add_curve_loop(lines)
        
        #create_surface
        base_surface = gmsh.model.occ.add_plane_surface([curve_loop])

        #create volume from top surface using extrude function in gmsh
        Palace_Gmsh_Renderer.dielectric_box = gmsh.model.occ.extrude([(2, base_surface)],0,0,self.size_z)

        #update model with the created geometry items
        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()

        

    def draw_polygon(self, polygon, component_name, index, flag):
        '''takes a shapely polygon object or pandas series 
            in as an argument and then draws it in Gmsh''' 

        #gets polygon coordinates depending on the type of argument fed to the function
        if (type(polygon) == shapely.geometry.polygon.Polygon):
            coords = polygon.exterior.coords.xy
        elif (type(polygon) == pd.Series):
            coords = polygon.geometry.exterior.coords.xy
        else:
            print('correct object not input') #need to change to throw exception

        #define lists to store points and lines
        points = []
        lines = []

        #num of points to draw, subtract 1 because the starting point is repeated in the coords list
        num_of_points = len(coords[0]) - 1
    
        #create 2D points in gmsh
        for i in range(num_of_points):
            point = gmsh.model.occ.addPoint(coords[0][i], coords[1][i], 
                                            self.design.parse_value(self.design.chips['main'].size.center_z))
            points.append(point)

        #index to access points list
        index_points = num_of_points-1
    
        #create lines from points
        for j in range(index_points):
            line = gmsh.model.occ.add_line(points[j], points[j+1])
            lines.append(line)        
        #add line connecting final point to beginning point
        line = gmsh.model.occ.add_line(points[num_of_points-1], points[0])
        lines.append(line) 
        
        #create curved loop
        curve_loop = gmsh.model.occ.add_curve_loop(lines)
        
        #create_surface
        surface = gmsh.model.occ.add_plane_surface([curve_loop])
        
        #check which metal elements have been created
        metal_list = []
        for i,value in enumerate(Palace_Gmsh_Renderer.metal_surfaces):
            metal_list.append(value[0])

        #add gmsh physical groups
        if self.design.qgeometry.get_component(component_name)[flag].loc[index].at['subtract'] == True:
            Palace_Gmsh_Renderer.dielectric_gap_surface.append((component_name, surface))
        elif component_name in metal_list:
            gmsh.model.addPhysicalGroup(2, [surface], name = component_name + '_' + str(surface))
            Palace_Gmsh_Renderer.metal_surfaces.append((component_name, surface))
        else:
            gmsh.model.addPhysicalGroup(2, [surface], name = component_name)
            Palace_Gmsh_Renderer.metal_surfaces.append((component_name, surface))

        #update model with the created geometry items
        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()



    def draw_path(self, path: pd.Series, component_name, index, flag):
        '''takes a pandas series in as an argument and then draws it in Gmsh'''

        #get width and buffer amount for path
        width = path.width
        buffer_amt = width/2

        #fillet path using QMplRenderer
        qmpl = QMplRenderer(None, self.design, None)
        path_filleted = qmpl.fillet_path(path)

        #buffer the path by the user defined width
        poly_path = shapely.buffer(path_filleted, distance = buffer_amt, cap_style = 'flat')

        #draw the polygon in Gmsh
        self.draw_polygon(poly_path, component_name, index, flag)



    def draw_air_box(self):

        #for air box choose increase in dimensions in 'mm'
        air_box_delta_x = (1/3) * self.size_x
        air_box_delta_y = (1/3) * self.size_y
        air_box_delta_z = 3 * self.size_z

        #half values of the sizes plus add increase for air box
        half_size_x = self.size_x/2 + air_box_delta_x/2
        half_size_y = self.size_y/2 + air_box_delta_y/2

        #store coordinates for surfaces of chip
        air_box_surface = {'point1': [self.center_x + half_size_x, self.center_y + half_size_y, self.center_z + air_box_delta_z],
                            'point2': [self.center_x + half_size_x, self.center_y - half_size_y, self.center_z + air_box_delta_z],
                            'point3': [self.center_x - half_size_x, self.center_y - half_size_y, self.center_z + air_box_delta_z],
                            'point4': [self.center_x - half_size_x, self.center_y + half_size_y, self.center_z + air_box_delta_z]}
        
        #define lists to store points, lines and surfaces
        points = []
        lines = []
        surfaces = []

        #add points for air_box
        for i,value in enumerate(air_box_surface):
            point = gmsh.model.occ.addPoint(air_box_surface[value][0], air_box_surface[value][1], 
                                    air_box_surface[value][2])
            points.append(point)

        #draw lines for airbox
        for j,value in enumerate(points):
            if(j<len(points)-1):
                line = gmsh.model.occ.add_line(points[j], points[j+1])
                lines.append(line)        
        line = gmsh.model.occ.add_line(points[len(points)-1], points[0])
        lines.append(line)        

        #create curved loop
        curve_loop = gmsh.model.occ.add_curve_loop(lines)
        
        #create_surface
        surface = gmsh.model.occ.add_plane_surface([curve_loop])

        #create volume from top surface of airbox using extrude function in gmsh
        air_box = gmsh.model.occ.extrude([(2, surface)],0,0,-2*air_box_delta_z)
        
        print('air_box:', air_box)

        #add physical group for dielectric chip base
        #gmsh.model.addPhysicalGroup(2, [air_box[1][1]], name = 'air_box')

        #cut out chip base from airbox
        air_box_cutout, air_box_cutout_map = gmsh.model.occ.fragment([air_box[1]], [(3, Palace_Gmsh_Renderer.frag_dielectric_vol)], 
                                                 removeObject=True, removeTool=False)

        #update model with the created geometry items
        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()

        #add physical group for dielectric chip base
        gmsh.model.addPhysicalGroup(3, [air_box_cutout[1][1]], name = 'air_box')

        #add physical group for the surfaces of the air box which will represent the far field bondary conditions
        gmsh.model.addPhysicalGroup(2, [surface, air_box[0][1], air_box[2][1], air_box[3][1], 
                                        air_box[4][1], air_box[5][1]], name = 'far_field')



    def add_ground_plane(self):
        '''add in ground plane on top of the dielectric box'''

        #list to store all dielectric gaps
        dielectric_gap_list = []
        metal_list = []
        cutout_list = []
        joint_list = []

        #create list of gaps to be cut from the base of the ground plane
        for i,value in enumerate(Palace_Gmsh_Renderer.dielectric_gap_surface):
            dielectric_gap_list.append((2,value[1]))
            joint_list.append((2,value[1]))

        #create list of gaps to be cut from the base of the ground plane
        for i,value in enumerate(Palace_Gmsh_Renderer.metal_surfaces):
            metal_list.append((2,value[1]))
            joint_list.append((2,value[1]))

        #update model with the created geometry items
        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()

        #create surface of ground plane
        rec2 = gmsh.model.occ.add_rectangle(self.center_x - self.size_x/2, 
                                            self.center_y - self.size_y/2, self.center_z,
                                            self.size_x, self.size_y)

        #update model with the created geometry items
        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()

        #cut out dielectric gaps from metal surface
        cutout_elements, cutout_map = gmsh.model.occ.cut([(2,rec2)], 
                                                                  dielectric_gap_list,
                                                                 removeObject=True, removeTool=True)

        #add physical group for ground plane
        ground_plane = []
        for i, value in enumerate(cutout_elements):
            ground_plane.append(value[1])
        gmsh.model.addPhysicalGroup(2, ground_plane, name = 'ground_plane')
        

        #update model with the created geometry items
        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()

        #add metal elements (eg: resonator, feed line, etc) to ground plane
        metal, metal_map = gmsh.model.occ.fragment(cutout_elements, 
                                                                 metal_list,
                                                                removeObject=True, removeTool=True)
        
        #update model with the created geometry items
        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()

        #uncomment for debugging pruposes
        print('metal elements in design:', metal_list)
        print('new metal elements:', metal)
        print('dielectric elements in design:', dielectric_gap_list)
        print('elements in ground plane:',cutout_elements)

        #create volume of dielectric base
        self.draw_chip_base()
        
        #set up elements to fragment with dielectric volume
        elements_to_fragment = metal_list
        elements_to_fragment.extend(cutout_elements)

        #append ports to list to fragment
        for i,value in enumerate(Palace_Gmsh_Renderer.ports_list):
            elements_to_fragment.append((2,value))

        #fragment the newly created elements with the dielectric volume
        chip_base, chip_base_map = gmsh.model.occ.fragment([Palace_Gmsh_Renderer.dielectric_box[1]], 
                                                                 elements_to_fragment,
                                                                removeObject=True, removeTool=True)
        
        Palace_Gmsh_Renderer.frag_dielectric_vol = chip_base[0][1]

        #add dielectric volume as chip base
        gmsh.model.addPhysicalGroup(3, [chip_base[0][1]], name = 'dielectric_base')



    def add_ports_on_launchpad(self, launch_pad):
        '''store launch pad objects to have ports created on'''
        Palace_Gmsh_Renderer.launch_pads_list.append(launch_pad)




    def create_ports_on_launchpad(self, launch_pad):
        '''create lumped port on launch pad for RF simulation'''

        #get indices for lp and pocket
        lp_index = self.design.qgeometry.get_component(launch_pad.name)['poly'].index[0]
        pocket_index = self.design.qgeometry.get_component(launch_pad.name)['poly'].index[1]

        #get pandas geodataframe for the launch pad
        lp_coords = self.design.qgeometry.get_component(launch_pad.name)['poly'].loc[lp_index].geometry.exterior.coords.xy
        pocket_coords = self.design.qgeometry.get_component(launch_pad.name)['poly'].loc[pocket_index].geometry.exterior.coords.xy

        #coordinates for first port
        start_point_lp_x_1 = lp_coords[0][2]
        start_point_lp_y_1 = lp_coords[1][2]
        start_point_poc_x = pocket_coords[0][2]
        start_point_poc_y = pocket_coords[1][2]
        x_diff = start_point_poc_x - start_point_lp_x_1
        y_diff = start_point_poc_y - start_point_lp_y_1

        #coordinates for second port
        start_point_lp_x_2 = lp_coords[0][3]
        start_point_lp_y_2 = lp_coords[1][3]

        #coordinates for second port
        scale = 0.5

        #Check orientation of the launchpad to define direction of lumped element ports
        if(launch_pad.options['orientation'] == '0' or launch_pad.options['orientation'] == '360'):
            dx1 = scale*-x_diff
            dy1 = y_diff
            dx2 = scale*-x_diff
            dy2 = -y_diff
        elif (launch_pad.options['orientation'] == '180' or launch_pad.options['orientation'] == '-180'):
            dx1 = scale*-x_diff
            dy1 = y_diff
            dx2 = scale*-x_diff
            dy2 = -y_diff
        elif (launch_pad.options['orientation'] == '270' or launch_pad.options['orientation'] == '-90'):
            dx1 = x_diff
            dy1 = scale*-y_diff
            dx2 = -x_diff
            dy2 = scale*-y_diff
        elif (launch_pad.options['orientation'] == '90' or launch_pad.options['orientation'] == '-270'):
            dx1 = x_diff
            dy1 = -scale*y_diff
            dx2 = -x_diff
            dy2 = -scale*-y_diff
        
        #draw the ports
        port1 = gmsh.model.occ.add_rectangle(start_point_lp_x_1, start_point_lp_y_1, self.center_z, dx1, dy1)
        port2 = gmsh.model.occ.add_rectangle(start_point_lp_x_2, start_point_lp_y_2, self.center_z, dx2, dy2)

        Palace_Gmsh_Renderer.ports_list.append(port1)
        Palace_Gmsh_Renderer.ports_list.append(port2)

        gmsh.model.addPhysicalGroup(2, [port1], name = 'port_'+launch_pad.name+'a')
        gmsh.model.addPhysicalGroup(2, [port2], name = 'port_'+launch_pad.name+'b')

        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()



    def view_design_components(self):
        '''view the wireframe with the physical groups of the design displayed'''

        gmsh.option.setNumber('Geometry.CurveWidth',0.1)
        gmsh.option.setNumber('Geometry.LabelType', 4)
        gmsh.option.setNumber('Geometry.Points', 0)
        gmsh.option.setNumber('Geometry.SurfaceLabels', 1)
        gmsh.option.setNumber('General.GraphicsFontSize', 14)
        gmsh.option.setColor('Geometry.Color.Surfaces',178,37,23, 190)
        gmsh.fltk.run()
    


    def print_physical_groups(self):
        '''print out the physical groups which are the boundary conditions for the simulation'''

        physical_groups = gmsh.model.get_physical_groups()

        for i, value in enumerate(physical_groups):
            phys_group_name = gmsh.model.get_physical_name(value[0], value[1])
            print('name:', phys_group_name + '\t\t', 'identifier:', value[1])

        
        

        