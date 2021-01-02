import os
import sys
from natsort import natsorted

##############################################################
#####################- FUNCTIONS - ###########################
##############################################################


def sim_file_sorting(file_extension, work_path):
    # the aim of this function is to sort the various simulation files into the order that Lumerical executes them, i.e. sim_1_var_1; sim_1_var_2; sim_2_var_1; etc
    # this will be useful as it will already sort and store the results in the same order as the simulations were ran.
    list_of_files = []
    for file in os.listdir(os.path.abspath(work_path)):  # list all files in the directory that was reached through the absolute path of the supplied directory
        if file.endswith(file_extension):
            list_of_files.append(file)
    return natsorted(list_of_files)

def mesh_density():
    #this function serves to calculate the number of yee cells per volume unit.
    #this is relevant because the yee cell mesh itself absorbs some power.
    #the purpose of this metric is to get a rough idea if a certain set of simulations can be directly compared between each other or not.
    #if the mesh density is too far apart (several orders of magnitude), the simulations cannot be directly compared without proper study of the impact of the mesh
    x_mesh =len()