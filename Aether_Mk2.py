import os
import sys
from natsort import natsorted


def sim_file_sorting(file_extension, work_path):

    list_of_files = [] 
    for file in os.listdir(os.path.abspath(work_path)):
        if file.endswith(file_extension):
            list_of_files.append(file)

    return natsorted(list_of_files)
