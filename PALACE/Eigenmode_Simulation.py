from SQDMetal.PALACE.Model import PALACE_Simulation_Base
import matplotlib.pyplot as plt
import numpy as np


class PALACE_Eigenmode_Simulation(PALACE_Simulation_Base):
    
    default_options = {'mesh_refinement': 1,
                            }                   

    def __init__(self, name, simulation_options):
        self.name = name