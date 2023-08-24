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

    #list to store entites which comprise the dielectric gaps and metals in the design
    dielectric_gap_surface = []
    metal_surfaces = []

    #class variable to store ID (int) representing the chip base
    dielectric_box = None
    
    #examine
    cut_ground_plane = None 
    cut_map = None

    metal_elements = None
    metal_map = None

    dielec_cut_list_elements = []

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
                print(component_name)
                for index in path_to_draw.index:
                    self.draw_path(path_to_draw.loc[index], component_name, index, flag = 'path')
            elif(not self.design.qgeometry.get_component(component_name)['poly'].empty):
                poly_to_draw = self.design.qgeometry.get_component(component_name)['poly']
                print(component_name)
                for index in poly_to_draw.index:
                    self.draw_polygon(poly_to_draw.loc[index], component_name, index, flag = 'poly')
            elif(not self.design.qgeometry.get_component(component_name)['junction'].empty):
                junc_to_draw = self.design.qgeometry.get_component(component_name)['junction']
                print(component_name)
                for index in junc_to_draw.index:
                    pass
        
        #add in physical groups for design in gmsh
        dielectric_gaps = []
        for i,value in enumerate(Palace_Gmsh_Renderer.dielectric_gap_surface):
            dielectric_gaps.append(value[1])
    
        
        self.add_ground_plane()
        #self.add_dielectric_cutout() 
        #self.draw_chip_base()
        #self.draw_air_box()

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

        #add physical group for dielectric chip base
        #gmsh.model.addPhysicalGroup(3, Palace_Gmsh_Renderer.dielectric_box[1], name = 'dielectric_base')

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
        
        #add gmsh physical groups
        if self.design.qgeometry.get_component(component_name)[flag].loc[index].at['subtract'] == True:
            Palace_Gmsh_Renderer.dielectric_gap_surface.append((component_name, surface))
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
        
        #cut out chip base from airbox
        air_box_cutout, air_box_cutout_map = gmsh.model.occ.cut([air_box[1]], [Palace_Gmsh_Renderer.dielectric_box[1]], 
                                                 removeObject=True, removeTool=False)

        #add physical group for dielectric chip base
        gmsh.model.addPhysicalGroup(3, [air_box_cutout[0][1]], name = 'air_box')

        #update model with the created geometry items
        gmsh.model.occ.synchronize()


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


        # dielectric_cutout_list = self.add_dielectric_cutout()
        # #create list of gaps to be cut from the base of the ground plane
        # for i,value in enumerate(dielectric_cutout_list):
        #     cutout_list.append((2,value[1]))

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

        cutout_elements, cutout_map = gmsh.model.occ.cut([(2,rec2)], 
                                                                  dielectric_gap_list,
                                                                 removeObject=True, removeTool=True)

        ground_plane = cutout_elements[-1]
        #add physical group for dielectric chip base
        gmsh.model.addPhysicalGroup(2, [ground_plane[1]], name = 'ground_plane')

        self.cut_ground_plane, self.cut_map = gmsh.model.occ.fragment([ground_plane], 
                                                                 metal_list,
                                                                removeObject=True, removeTool=True)
        

        # self.metal_elements, self.metal_map = gmsh.model.occ.fragment([(2,rec2)], 
        #                                                           metal_list,
        #                                                          removeObject=True, removeTool=False)

        #remove dielectric elements now that ground plane is created
        #gmsh.model.remove_entities(dielectric_gap_list, recursive=True)

        print(rec2)

        #update model with the created geometry items
        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()
        print('metal tags:', metal_list)
        print('new metal tags:', self.cut_ground_plane)
        print('cutout tags:', dielectric_gap_list)
        print(' new cut out:',cutout_elements)
        #print(cutout_elements)

        #add physical group for dielectric chip base
        #gmsh.model.addPhysicalGroup(2, [self.cut_ground_plane[0][1]], name = 'ground_plane')

        self.draw_chip_base()
        
        elements_to_fragment = metal_list
        elements_to_fragment.append(ground_plane)

        chip_base, chip_base_map = gmsh.model.occ.fragment([Palace_Gmsh_Renderer.dielectric_box[1]], 
                                                                 elements_to_fragment,
                                                                removeObject=True, removeTool=True)
        
        print(chip_base)
        gmsh.model.addPhysicalGroup(3, [chip_base[0][1]], name = 'dielectric_base')



    def add_dielectric_cutout(self):
        '''Method to add in the dielectric cutout around the metal elements'''
        
        dielec_cut_list = []
        dielec_cut_list_elements = []
        for i, value_dielec in enumerate(self.dielectric_gap_surface):
            for j, value_metal in enumerate(self.metal_surfaces):
                #cut out metal from dielectric
                if value_dielec[0] == value_metal[0]:
                    dielec_cut, dielec_cut_map = gmsh.model.occ.cut([(2,value_dielec[1])], [(2,value_metal[1])], removeObject=False, removeTool=False)
                    dielec_cut_list.append(dielec_cut)
                #cut out open to ground regions
                elif (not self.design.qgeometry.get_component(value_dielec[0])['poly'].empty and 
                      self.design.qgeometry.get_component(value_dielec[0])['poly'].iloc[0][1] == 'open_to_ground'):
                    Palace_Gmsh_Renderer.dielec_cut_list_elements.append(value_dielec[1])

        gmsh.model.geo.synchronize()
        gmsh.model.occ.synchronize()

        #extract gmsh element tags from dielec_cut_list
        for i,value_1 in enumerate(dielec_cut_list):
            for j,value_2 in enumerate(value_1):
                Palace_Gmsh_Renderer.dielec_cut_list_elements.append(value_2[1])

        #create physical group
        #gmsh.model.addPhysicalGroup(2, dielec_cut_list_elements, name = "dielectric_cutout")

        return dielec_cut_list_elements

    def fuse_elements(self, element_1, element_2):
        '''fuse elements defined in the physical group'''

        

        