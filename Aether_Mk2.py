from PyQt5 import QtWidgets#, QtGui

from PyQt5.QtWidgets import QFileDialog,QTableWidgetItem# QAction,
from Aether_GUI import Ui_MainWindow  # importing our generated file

import sys
import importlib.util
import numpy as np
import pandas
from scipy.constants import nu2lambda
import scipy.constants as sc
from scipy import interpolate

import scipy.integrate as integrate
import os
from natsort import natsorted
import configparser
import time
from pathlib import Path
import openpyxl

config = configparser.ConfigParser()

def mesh_per_volume():
    mesh_x=len(fdtd.getresult("FDTD","x"))
    mesh_y = len(fdtd.getresult("FDTD", "y"))
    mesh_z = len(fdtd.getresult("FDTD", "z"))
    mesh_volume=mesh_x*mesh_y*mesh_z

    fdtd.select("FDTD")
    simulation_x = fdtd.get("x span")
    simulation_y = fdtd.get("y span")
    simulation_z = fdtd.get("z span")
    simulation_volume=simulation_x*simulation_y*simulation_z

    mesh_volume_ratio=mesh_volume/simulation_volume
    return mesh_volume_ratio

def File_Sorting(file_extension,work_path):
    os.chdir(work_path)
    list_of_simulations = []
    for file in os.listdir(os.path.abspath(os.getcwd())):
        if file.endswith(file_extension):
            list_of_simulations.append(file)
    return natsorted(list_of_simulations)

def CIGS_EQE(monitor_CIGS, material_CIGS, generation_toggle,file_name,tolerance):
    global gen_dataset
    if generation_toggle:
        gen_dataset = fdtd.getresult(monitor_CIGS, "G_export")
    pabs_dataset = fdtd.getresult(monitor_CIGS, "Pabs")
    simulated_index = fdtd.getresult(monitor_CIGS + "::index", "index_detail")
    frequency = np.squeeze(fdtd.getresult(monitor_CIGS, "f"))
    length_f = len(pabs_dataset['f'])
    fourD_material_map = np.empty(shape=(len(pabs_dataset['x']), len(pabs_dataset['y']), len(pabs_dataset['z']), length_f))
    truncated_frequency = frequency[0:1]
    # ------------------------------------------------#
    material_index_data = fdtd.getfdtdindex(material_CIGS, truncated_frequency, min(truncated_frequency), max(truncated_frequency))
    material_real_index = np.real(material_index_data)
    material_imag_index = np.imag(material_index_data)
    # ------------------------------------------------#
    index_z=simulated_index['index_z'] #replaces the last Z value by 0 in order to eliminate the surface errros caused by interpolation
    index_z[:,:,-1,:]=0
    reshaped_simulated_index = index_z[:, :, :, 0]
    simulated_real_index = np.real(reshaped_simulated_index)
    simulated_imag_index = np.imag(reshaped_simulated_index)
    # ------------------------------------------------#

    material_map_real = np.isclose(material_real_index[0], simulated_real_index, atol=tolerance)
    material_map_imag = np.isclose(material_imag_index[0], simulated_imag_index, atol=tolerance)
    material_map_real = material_map_real.astype(int)
    material_map_imag = material_map_imag.astype(int)
    material_map = np.multiply(material_map_real, material_map_imag)
    # ------------------------------------------------#
    for i in range(length_f):
        fourD_material_map[:, :, :, i] = material_map #repeates the same material map throughout all frequencies since the structure does not depend on the λ
        #to obtain the EQE, we need the fourD_material_map, in order to get the material map for each frequency
    fdtd.select(monitor_CIGS)
    periods = int(fdtd.get("periods"))

    dimensions_integrate = np.arange(1, 4) #a numpy range from 1 to 4 (x,y,z,f)
    CIGS_Pabs_map = pabs_dataset['Pabs'] * fourD_material_map
    CIGS_EQE = fdtd.integrate2(CIGS_Pabs_map, dimensions_integrate, pabs_dataset['x'], pabs_dataset['y'], pabs_dataset['z'])

    if generation_toggle:
        if fdtd.get("generation_frequency"):
            fourD_material_map = np.tile(fourD_material_map[:, :, :,:], (periods,1, 1, 1))
            generation_map = gen_dataset['G'][0:, :, :,:] * fourD_material_map
        else:
            material_map = np.tile(material_map[0:, :, :], (periods, 1, 1))
            generation_map = gen_dataset['G'][0:, :, :] * material_map
    # ------------------------------------------------#


    # ------------------------------------------------#
    if generation_toggle:
        # now we need to send the generation dataset back to FDTD and save it
        gen_dataset.update({'G': generation_map})  # here I replace the "G" on the original dataset with the new Gen Map. #This allows me to keep the dataset structure and allows for a easier import into Charge
        try:
            fdtd.clear()
            #fdtd.clearexcept("gen_data_"+str(file_name[:-4]));  # remove all variables on the workspace
        except:
            pass #if the FDTD workspace is already clear, it will return an error. Which will crash if python script.
            # So far, it seems to be safe to handle it by just "passing"

        fdtd.putv("gen_data_"+str(file_name[:-4]), gen_dataset)  # we send the dataset back to FDTD #[:-4] to remove the ".fsp" part of the string and prevent an error

        fdtd.matlabsave("gen_data_"+str(file_name[:-4]))  # FDTD will export all the variables in the workspace into a .mat file. But since the fdtd.clear erased all vars, it's just the desired generation dataset that is exported
        gen_dataset=None

    return CIGS_EQE

def parasitic_eqe(monitor_substrate, material_CIGS, tolerance):
    pabs_dataset = fdtd.getresult(monitor_substrate, "Pabs")
    simulated_index = fdtd.getresult(monitor_substrate + "::index", "index_x")
    frequency = np.squeeze(fdtd.getresult(monitor_substrate + "::index", "f"))
    length_f = len(pabs_dataset['f'])
    fourD_material_map = np.empty(shape=(len(pabs_dataset['x']), len(pabs_dataset['y']), len(pabs_dataset['z']), length_f))
    truncated_frequency = frequency[0:1]
    # ------------------------------------------------#
    material_index_data = fdtd.getfdtdindex(material_CIGS, truncated_frequency, min(truncated_frequency), max(truncated_frequency))
    material_real_index = np.real(material_index_data)
    material_imag_index = np.imag(material_index_data)
    # ------------------------------------------------#
    reshaped_simulated_index = simulated_index[:, :, :, 0]
    simulated_real_index = np.real(reshaped_simulated_index)
    simulated_imag_index = np.imag(reshaped_simulated_index)
    # ------------------------------------------------#
    material_map_real = np.invert(np.isclose(material_real_index[0], simulated_real_index, atol=tolerance))
    material_map_imag = np.invert(np.isclose(material_imag_index[0], simulated_imag_index, atol=tolerance))
    material_map_real = material_map_real.astype(int)
    material_map_imag = material_map_imag.astype(int)
    material_map = np.multiply(material_map_real, material_map_imag)

    # ------------------------------------------------#
    for i in range(length_f):
        fourD_material_map[:, :, :, i] = material_map
    # ------------------------------------------------#
    dimensions_integrate = np.arange(1, 4)
    parasitic_Pabs_map = pabs_dataset['Pabs'] * fourD_material_map
    parasitic_EQE = fdtd.integrate2(parasitic_Pabs_map, dimensions_integrate, pabs_dataset['x'], pabs_dataset['y'], pabs_dataset['z'])

    return parasitic_EQE

def reflection_results(reflection_monitor):
    reflection_data=fdtd.getresult(reflection_monitor,'T')

    reflection_eqe=1+reflection_data['T']
    return reflection_eqe

def electrical_results(contact_name):

   electrical_result=charge.getresult('CHARGE',contact_name)
   current=electrical_result['I']
   voltage=electrical_result['V_'+contact_name]
   length_x=charge.getnamed('simulation region','x span')
   length_y=charge.getnamed('simulation region','y span')
   density=-(current/(10*length_x*length_y))
   #If the simulation just ran 1 value, that means that the power, voc, etc cannot be calculated
   #But the function needs to returns something anyway. So it will return a numpy array with a 0
    #Why a numpy array? Because the dataframes code is expecting a numpy array and will crash otherwise.
   if len(electrical_result['I']) !=1:
       interpolator=interpolate.UnivariateSpline(voltage,density,s=0)
       Voc=interpolator.roots()
       Jsc=[interpolator(0)]
       #smallest_density_index = np.argmin(abs(density))
       Power=voltage*interpolator(voltage)
       eta=[np.amax(Power)]
       FF=np.amax(Power)/(Voc*Jsc)
       Pmax=[np.amax(Power)]
   else:
       Power=np.zeros(1)
       Voc=np.zeros(1)
       Jsc=density
       eta=np.zeros(1)
       FF=np.zeros(1)
       Pmax=np.zeros(1)
   return voltage,density,Power,Voc,Jsc,eta,FF,Pmax

class mywindow(QtWidgets.QMainWindow):


    def __init__(self):
        super(mywindow, self).__init__()
        self.table = QtWidgets.QTableView()
        self.setFixedSize(824,924)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        self.ui.pushButton_browse.clicked.connect(self.lumapiClicked)

        self.ui.data_browse.clicked.connect(self.dataBrowseClicked)

        self.ui.actionExit.triggered.connect(self.exitTriggered)

        self.ui.actionSave.triggered.connect(self.saveTriggered)

        self.ui.actionOpen.triggered.connect(self.openTriggered)

        self.ui.pushButton_start.clicked.connect(self.startClicked)

        self.ui.pushButton_start_CHARGE.clicked.connect(self.startChargeClicked)

        self.ui.data_browse_CSV.clicked.connect(self.dataCSVClicked)

        self.ui.pushButton_process_charge.clicked.connect(self.process_CHARGE)

        self.ui.data_browse_CHARGE_Process.clicked.connect(self.ChargeBrowseClicked)
        #self.ui.data_browse_CHARGE.clicked.connect(self.dataCHARGEClicked)
        #self.ui.pushButton_TESTING.clicked.connect(self.TESTING)

        self.ui.data_browse_CHARGE_BLUEPRINT.clicked.connect(self.dataCHARGE_BLUEPRINT_Clicked)
       # self. ui.pushButton_Force_Gen.clicked.connect(self.forceGenLoad)
    def lumapiClicked(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        global lumapiFile
        lumapiFile, _ = QFileDialog.getOpenFileName(self, "Lumapi API Location?", "", "Python Files (*.py)", options=options) # , _ is needed or if the user does not input a file the program will crash
        if lumapiFile:
            self.ui.lineEdit_lumapi.setText(lumapiFile)

    def dataBrowseClicked(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        data_folder= QFileDialog.getExistingDirectory(self,"Select a folder with FDTD files:",'',options=QFileDialog.ShowDirsOnly)
        if data_folder:

            # list_of_simulations = []
            # for file in os.listdir(data_folder):
            #     if file.endswith(".fsp"):
            #         list_of_simulations.append(file)
            global sorted_FDTD
            sorted_FDTD = File_Sorting(".fsp",data_folder)
            self.ui.lineEdit_data_folder.setText(data_folder)
            self.ui.Files_List.setRowCount(len(sorted_FDTD))
            self.ui.Files_List.setColumnCount(2)

            self.ui.Results_List.setRowCount(len(sorted_FDTD))
            self.ui.Results_List.setColumnCount(3)

            self.ui.Files_List.setHorizontalHeaderLabels(('Simulation Files', 'Status'))
            for i,sims_names in enumerate(sorted_FDTD):
                self.ui.Files_List.setItem(i,0,QTableWidgetItem(sims_names))
                self.ui.Files_List.setItem(i,1,QTableWidgetItem("Not Processed"))

    def dataCSVClicked(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        CSVFile, _ = QFileDialog.getOpenFileName(self, "CSV file with model data", "", "CSV Files (*.csv)", options=options)
        if CSVFile:
            self.ui.lineEdit_data_CSV.setText(CSVFile)
            model_dataframe = pandas.read_csv(CSVFile)
            self.ui.Files_List_CSV.setRowCount(len(model_dataframe.index))
            self.ui.Files_List_CSV.setColumnCount(len(model_dataframe.columns))
            self.ui.Files_List_CSV.setHorizontalHeaderLabels(model_dataframe.columns.values.tolist())
            for i in range(len(model_dataframe.index)):
                for j in range(len(model_dataframe.columns)):
                    self.ui.Files_List_CSV.setItem(i,j,QTableWidgetItem(str(model_dataframe.iat[i,j])))
                    #self.ui.Files_List_CSV.setItem(i,j,QTableWidgetItem(str(round(model_dataframe.iat[i,j],3) if np.char.isnumeric(str(model_dataframe.iat[i,j])) else str(model_dataframe.iat[i,j])  )))
    def dataCHARGE_BLUEPRINT_Clicked(self):
        options = QFileDialog.Options()
        #options |= QFileDialog.DontUseNativeDialog
        global blueprint_file
        blueprint_file, _ = QFileDialog.getOpenFileName(self, "Select the CHARGE blueprint file:",'', 'LDEV Files (*.ldev)', options=options)
        if blueprint_file:
            self.ui.lineEdit_sims_CHARGE_BLUEPRINT.setText(blueprint_file)

    def exitTriggered(self):
       sys.exit(app.exec())

    def openTriggered(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        configFile, _ = QFileDialog.getOpenFileName(self, "QFileDialog.getOpenFileName()", "", "INI Files (*.ini)", options=options)
        if configFile:
            config.read(configFile)
            self.ui.lineEdit_lumapi.setText(config['lumapi']['path'])
            self.ui.lineEdit_Mat_CIGS.setText(config['Monitors']['CIGS_Material']) #reads the value on the config file and sets it to the lineedit object
            self.ui.lineEdit_Substr_M.setText(config['Monitors']['Substrate_Monitors'])
            self.ui.lineEdit_CIGS_M.setText(config['Monitors']['CIGS_Monitors'])
            self.ui.lineEdit_Reflection_M.setText(config['Monitors']['Reflection_Monitors'])

            self.ui.checkBox_EQE.setChecked(config.getboolean('Settings FDTD','Calculate_CIGS_EQE'))
            self.ui.checkBox_Gen.setChecked(config.getboolean('Settings FDTD','Export_Generation'))
            self.ui.checkBox_Parasit.setChecked(config.getboolean('Settings FDTD','Calculate_Substrate_EQE'))
            self.ui.checkBox_export_charge.setChecked(config.getboolean('Settings FDTD', 'Export_to_Charge'))
            self.ui.checkBox_Mesh.setChecked(config.getboolean('Settings FDTD','Export_Mesh'))
            self.ui.doubleSpinBox.setValue(config.getfloat('Settings FDTD','Tolerance'))

            self.ui.checkBox_Auto_Charge.setChecked(config.getboolean('Settings CHARGE','Auto_Charge'))
            self.ui.checkBox_Force_Gen.setChecked(config.getboolean('Settings CHARGE', 'Force_Charge'))

    def saveTriggered(self):
        config['lumapi']={}
        if self.ui.lineEdit_lumapi.text():
            config['lumapi']['path']=self.ui.lineEdit_lumapi.text()

        config['Monitors'] = {}
        if self.ui.lineEdit_CIGS_M.text():
            config['Monitors']['CIGS_Monitors'] = self.ui.lineEdit_CIGS_M.text()
        else:
            config['Monitors']['CIGS_Monitors']=''

        if self.ui.lineEdit_Reflection_M.text():
            config['Monitors']['Reflection_Monitors'] = self.ui.lineEdit_Reflection_M.text()
        else:
            config['Monitors']['Reflection_Monitors'] = ''

        if self.ui.lineEdit_Substr_M.text():
            config['Monitors']['Substrate_Monitors'] = self.ui.lineEdit_Substr_M.text()
        else:
            config['Monitors']['Substrate_Monitors'] =''

        if self.ui.lineEdit_Mat_CIGS.text():
            config['Monitors']['CIGS_Material'] = self.ui.lineEdit_Mat_CIGS.text()
        else:
            config['Monitors']['CIGS_Material']=''


        config['Settings FDTD']={}
        config['Settings FDTD']['Export_Generation']=str(self.ui.checkBox_Gen.isChecked())
        config['Settings FDTD']['Calculate_CIGS_EQE']=str(self.ui.checkBox_EQE.isChecked())
        config['Settings FDTD']['Calculate_Substrate_EQE'] = str(self.ui.checkBox_Parasit.isChecked())
        config['Settings FDTD']['Export_to_Charge'] = str(self.ui.checkBox_export_charge.isChecked())
        config['Settings FDTD']['Export_Mesh']=str(self.ui.checkBox_Mesh.isChecked())
        config['Settings FDTD']['Tolerance']=self.ui.doubleSpinBox.text()

        config['Settings CHARGE']={}
        config['Settings CHARGE']['Auto_Charge'] = str(self.ui.checkBox_Auto_Charge.isChecked())
        config['Settings CHARGE']['Force_Charge'] = str(self.ui.checkBox_Force_Gen.isChecked())
        with open('Aether_Config.ini', 'w') as Configfile:
            config.write(Configfile)

    def startClicked (self):
        # place here the location of the Lumerical "lumapi" library
        spec_win = importlib.util.spec_from_file_location("lumapi", self.ui.lineEdit_lumapi.text())
        # Functions that perform the actual loading
        lumapi = importlib.util.module_from_spec(spec_win)
        spec_win.loader.exec_module(lumapi)
        global fdtd, CIGS_Reflection_dataframe

        #fdtd = lumapi.FDTD(hide=True)

        CIGS_Monitors=self.ui.lineEdit_CIGS_M.text()
        CIGS_Material=self.ui.lineEdit_Mat_CIGS.text()
        Substrate_Monitors=self.ui.lineEdit_Substr_M.text()
        Reflection_Monitors=self.ui.lineEdit_Reflection_M.text()

        for i,files in enumerate(sorted_FDTD):


            self.ui.Results_List.setItem(i,0,QTableWidgetItem(str(files)))
            self.ui.Results_List.setItem(i,1, QTableWidgetItem("0"))
            self.ui.Results_List.setItem(i,2, QTableWidgetItem("0"))

            if i == 0:
                self.ui.Results_List.setHorizontalHeaderLabels(('Simulation Files', 'CIGS Jsc', 'parasitic Jsc')) #this only needs to run once to set the name of the headers

                fdtd = lumapi.FDTD(hide=False, filename=files)
                if not fdtd.havedata(CIGS_Monitors):
                    fdtd.runanalysis
                    fdtd.save

                if self.ui.checkBox_EQE.isChecked():

                    CIGS_EQE_dataframe = pandas.DataFrame(nu2lambda(np.squeeze(fdtd.getresult(CIGS_Monitors, "f"))), columns=['λ'])
                    CIGS_Jsc_list = []
                    CIGS_Reflection_dataframe = pandas.DataFrame(nu2lambda(np.squeeze(fdtd.getresult(Reflection_Monitors, "f"))), columns=['λ'])
                if self.ui.checkBox_Parasit.isChecked():
                    parasitic_EQE_dataframe = pandas.DataFrame(nu2lambda(np.squeeze(fdtd.getresult(CIGS_Monitors,"f"))),columns=['λ'])
                    parasitic_Jsc_list = []

                if self.ui.checkBox_export_charge.isChecked():
                    model_dataframe=pandas.DataFrame(columns=['FDTD_name','Metal_thickness','Dielectric_thickness','Encapsulation_gap','contact_dimension'])
                    fdtd.select("::model")
                    model_dataframe.loc[i]=[files,fdtd.get('Metal_thickness'),fdtd.get('Dielectric_thickness'),fdtd.get('Encapsulation_gap'),fdtd.get('contact_dimension')]

                if self.ui.checkBox_Mesh.isChecked():
                    CIGS_Mesh_list=[]

                # calculate stuff for the jsc
                AM15G = fdtd.solar(1)
                AM15G = np.asarray(AM15G).squeeze()
                wavelength = fdtd.solar(0)
                wavelength = np.asarray(wavelength).squeeze()
                qhc = (sc.e / (sc.h * sc.c))

                frequency = nu2lambda(fdtd.getresult(CIGS_Monitors, "f"))
                frequency = frequency.flatten()
                interpolator = interpolate.interp1d(wavelength, AM15G)
                interpolated_AM15G = interpolator(frequency)


            if i != 0:
                fdtd.load(files)
                if not fdtd.havedata(CIGS_Monitors):
                    fdtd.runanalysis
                    fdtd.save

            if self.ui.checkBox_export_charge.isChecked():
                fdtd.select("::model")
                model_dataframe.loc[i] = [files, fdtd.get('Metal_thickness'), fdtd.get('Dielectric_thickness'), fdtd.get('Encapsulation_gap'), fdtd.get('contact_dimension')]

            if self.ui.checkBox_EQE.isChecked():

                tolerance=float(self.ui.doubleSpinBox.text()) #takes the tolerance to the isclose function from the Spinbox UI element
                CIGS_EQE_data=CIGS_EQE(CIGS_Monitors,CIGS_Material,self.ui.checkBox_Gen.isChecked(),files,tolerance)
                CIGS_EQE_data=CIGS_EQE_data.flatten()
                CIGS_EQE_dataframe.insert(i+1,"EQE "+files,CIGS_EQE_data)
                Jsc_integrated=integrate.simps(frequency*qhc*interpolated_AM15G*CIGS_EQE_data,frequency)*0.1 #mA/cm^2
                CIGS_Jsc_list.append(Jsc_integrated)
                self.ui.Results_List.setItem(i, 1, QTableWidgetItem(str(Jsc_integrated)))
                CIGS_Reflection_data=reflection_results(Reflection_Monitors)

                CIGS_Reflection_dataframe.insert(i+1,"Reflection EQE "+files,CIGS_Reflection_data)

            if self.ui.checkBox_Parasit.isChecked():
                tolerance = float(self.ui.doubleSpinBox.text())
                parasitic_EQE_data = parasitic_eqe(Substrate_Monitors, CIGS_Material, tolerance)
                parasitic_EQE_data =parasitic_EQE_data.flatten()
                parasitic_EQE_dataframe.insert(i+1,"parasitic EQE "+files,parasitic_EQE_data)
                parasitic_Jsc_integrated= integrate.simps(frequency*qhc*interpolated_AM15G*parasitic_EQE_data,frequency)*0.1
                parasitic_Jsc_list.append(parasitic_Jsc_integrated)
                self.ui.Results_List.setItem(i, 2, QTableWidgetItem(str(parasitic_Jsc_integrated)))

            if self.ui.checkBox_Mesh.isChecked():
                Mesh_result= mesh_per_volume()
                CIGS_Mesh_list.append(Mesh_result)

            self.ui.Files_List.setItem(i, 1, QTableWidgetItem('Processed'))



        if self.ui.checkBox_EQE.isChecked():

            CIGS_EQE_dataframe.to_excel('CIGS_EQE.xlsx', index=False)
            CIGS_Reflection_dataframe.to_excel('Reflection.xlsx',index=False)
            CIGS_Jsc_dict = dict(zip(sorted_FDTD, CIGS_Jsc_list))
            CIGS_Jsc_dataframe = pandas.DataFrame([CIGS_Jsc_dict])
            CIGS_Jsc_dataframe.to_excel('CIGS_Jsc.xlsx', index=False)

        if self.ui.checkBox_Parasit.isChecked():
            parasitic_EQE_dataframe.to_excel('parasitic_eqe.xlsx', index=False)
            parasitic_Jsc_dict = dict(zip(sorted_FDTD, parasitic_Jsc_list))
            parasitic_Jsc_dataframe = pandas.DataFrame([parasitic_Jsc_dict])
            parasitic_Jsc_dataframe.to_excel('parasitic_Jsc.xlsx', index=False)

        if self.ui.checkBox_export_charge.isChecked():
            model_dataframe.to_csv("model_dataframe.csv", index=False)

        if self.ui.checkBox_Mesh.isChecked():
            CIGS_Mesh_dict=dict(zip(sorted_FDTD, CIGS_Mesh_list))
            CIGS_Mesh_dataframe=pandas.DataFrame([CIGS_Mesh_dict])
            CIGS_Mesh_dataframe.to_excel('Mesh_Volume.xlsx', index=False)

        fdtd.close()


    def startChargeClicked (self):
        os.chdir(self.ui.lineEdit_sims_CHARGE_Process.text()) #change the directory to the working folder
        # place here the location of the Lumerical "lumapi" library
        spec_win = importlib.util.spec_from_file_location("lumapi", self.ui.lineEdit_lumapi.text())
        # Functions that perform the actual loading
        lumapi = importlib.util.module_from_spec(spec_win)
        spec_win.loader.exec_module(lumapi)
        global charge
        charge= lumapi.DEVICE(hide=False,filename=blueprint_file)
        try:
            charge.deletesweep('Auto_Sweep')
        except:
            pass

        if self.ui.checkBox_Auto_Charge.isChecked():

            model_dataframe=pandas.read_csv(self.ui.lineEdit_data_CSV.text())
            charge.addsweep()
            time.sleep(2) #the sleep timers are for the script to wait 2 seconds, to prevent CHARGE from doing weird stuff
            charge.setsweep("sweep", "name", "Auto_Sweep");
            time.sleep(2)
            charge.setsweep("Auto_Sweep",'type','Values')
            time.sleep(2)
            charge.setsweep("Auto_Sweep",'number of points',len(model_dataframe.index))
            time.sleep(2)
            for (columnName, columnData) in model_dataframe.iteritems():

                if columnName == 'Metal_thickness':
                    sweep_parameter = {}
                    sweep_parameter['Parameter'] = '::model::Metal_thickness'
                    sweep_parameter['Type'] = 'Length'
                    for i, name in enumerate(model_dataframe['Metal_thickness']):
                        sweep_parameter['Value_' + str(i + 1)] =model_dataframe['Metal_thickness'][i]
                    charge.addsweepparameter('Auto_Sweep',sweep_parameter)

                ####################    Dielectric_thickness ##############
                if columnName == 'Dielectric_thickness':
                    sweep_parameter = {}
                    sweep_parameter['Parameter'] = '::model::Dielectric_thickness'
                    sweep_parameter['Type'] = 'Length'
                    for i, name in enumerate(model_dataframe['Dielectric_thickness']):
                        sweep_parameter['Value_' + str(i + 1)] = model_dataframe['Dielectric_thickness'][i]
                    charge.addsweepparameter('Auto_Sweep', sweep_parameter)

                ##################   Encapsulation_gap ################
                if columnName == 'Encapsulation_gap':
                    sweep_parameter = {}
                    sweep_parameter['Parameter'] = '::model::Encapsulation_gap'
                    sweep_parameter['Type'] = 'Length'
                    for i, name in enumerate(model_dataframe['Encapsulation_gap']):
                        sweep_parameter['Value_' + str(i + 1)] = model_dataframe['Encapsulation_gap'][i]
                    charge.addsweepparameter('Auto_Sweep', sweep_parameter)


                ##################   contact_dimension ################
                if columnName == 'contact_dimension':
                    sweep_parameter = {}
                    sweep_parameter['Parameter'] = '::model::contact_dimension'
                    sweep_parameter['Type'] = 'Length'
                    for i, name in enumerate(model_dataframe['contact_dimension']):
                        sweep_parameter['Value_' + str(i + 1)] = model_dataframe['contact_dimension'][i]
                    charge.addsweepparameter('Auto_Sweep', sweep_parameter)

            charge.savesweep()

            ####### GEN DATASETS ###############
            # list_of_gen_datasets = []
            # for file in os.listdir(os.path.abspath(os.getcwd())):
            #     if file.endswith(".mat"):
            #         list_of_gen_datasets.append(file)

            sorted_GEN_DATASET = File_Sorting(".mat",os.getcwd())
            if self.ui.checkBox_Force_Gen.isChecked():

                os.chdir(os.path.abspath(os.getcwd())+"\\"+Path(self.ui.lineEdit_sims_CHARGE_BLUEPRINT.text()).stem + '_Auto_Sweep')
                ####### CHARGE SIMS ###############
                sorted_CHARGE = File_Sorting(".ldev",os.getcwd())

                for i, sim_name in enumerate(sorted_CHARGE):

                    charge.cd(os.path.abspath(os.getcwd()))#jump into the sweep folder to get the sweep files
                    charge.load(sim_name)
                    charge.addimportgen()
                    charge.select('CHARGE::generation')
                    os.chdir('..')#jump OUT of the sweep folder to grab the generation datasets
                    charge.cd(os.path.abspath(os.getcwd()))
                    charge.importdataset(sorted_GEN_DATASET[i])
                    os.chdir(os.path.abspath(os.getcwd()) + "\\" + Path(self.ui.lineEdit_sims_CHARGE_BLUEPRINT.text()).stem + '_Auto_Sweep') #jump back into the sweep folder to repeat the process
                    charge.addjob(sim_name, "CHARGE")
                    charge.save()

        charge.runjobs('CHARGE')
        charge.save()

    def ChargeBrowseClicked(self):
        options = QFileDialog.Options()
        # options |= QFileDialog.DontUseNativeDialog
        global data_folder_charge
        data_folder_charge = QFileDialog.getExistingDirectory(self, "Select a folder with CHARGE files:", '', options=QFileDialog.ShowDirsOnly)
        if data_folder_charge:
            self.ui.lineEdit_sims_CHARGE_Process.setText(data_folder_charge)

        os.chdir(self.ui.lineEdit_sims_CHARGE_Process.text())
        global  sorted_CHARGE
        sorted_CHARGE = File_Sorting(".ldev",data_folder_charge)

        self.ui.Files_List_CHARGE.setRowCount(len(sorted_CHARGE))
        self.ui.Files_List_CHARGE.setColumnCount(2)
        self.ui.Files_List_CHARGE.setHorizontalHeaderLabels(('Simulation Files', 'Status'))

        for i, sims_names in enumerate(sorted_CHARGE):
            self.ui.Files_List_CHARGE.setItem(i, 0, QTableWidgetItem(sims_names))
            self.ui.Files_List_CHARGE.setItem(i, 1, QTableWidgetItem("Not Processed"))

    def process_CHARGE(self):
        os.chdir(self.ui.lineEdit_sims_CHARGE_Process.text())
        spec_win = importlib.util.spec_from_file_location("lumapi", self.ui.lineEdit_lumapi.text())
        # Functions that perform the actual loading
        lumapi = importlib.util.module_from_spec(spec_win)
        spec_win.loader.exec_module(lumapi)
        global charge
        charge = lumapi.DEVICE(hide=False)


        self.ui.Results_List_CHARGE.setRowCount(len(sorted_CHARGE))
        self.ui.Results_List_CHARGE.setColumnCount(6)
        self.ui.Results_List_CHARGE.setHorizontalHeaderLabels(('Simulation Files', 'Voc', 'Jsc', 'Eff', 'FF', 'Pmax'))


        for i, sim_name in enumerate(sorted_CHARGE):
            charge.load(sim_name)
            results=electrical_results('Mo_contact')
            #voltage,density,Power,Voc,Jsc,eta,FF,Pmax

            if i == 0:
                density_dataframe = pandas.DataFrame(results[0],columns=['V (Volts)'])
                density_dataframe.insert(i + 1, sim_name + ' Jsc (mA/cm^2)', results[1])

                power_dataframe=pandas.DataFrame(results[0],columns=[' V (Volts)'])
                power_dataframe = pandas.concat([power_dataframe, pandas.DataFrame(results[2])], ignore_index=True, axis=1)
                #getting weird behaviour if I just "insert" the data without 1st creating the dataframe WITH data
                Voc_dataframe=pandas.DataFrame(results[3],columns=[sim_name+' Voc (Volts)'])
                Jsc_dataframe=pandas.DataFrame(results[4], columns=[sim_name+' Jsc (mA)'])
                eta_dataframe=pandas.DataFrame(results[5], columns=[sim_name + ' Efficiency'])
                FF_dataframe = pandas.DataFrame(results[6],columns=[sim_name + ' Fill-Factor'])
                Pmax_dataframe = pandas.DataFrame(results[7],columns=[sim_name + ' Max Power(mW/cm^2)'])

            if i != 0:

                density_dataframe.insert(i+1,sim_name+' Jsc (mA/cm^2)',pandas.Series(list(results[1]))) #Converts the tuple from results to a list, so that the pd.series can add
                #collumns with diferent lenghts to the dataframe, in case the simulations have diferent voltages
                power_dataframe=pandas.concat([power_dataframe,pandas.DataFrame(results[2])],ignore_index=True, axis=1)
                Voc_dataframe.insert(i,sim_name+' Voc (Volts)',results[3])
                Jsc_dataframe.insert(i,sim_name+' Jsc (mA)',results[4])
                eta_dataframe.insert(i, sim_name + ' Efficiency', results[5])
                FF_dataframe.insert(i, sim_name + ' Fill-Factor', results[6])
                Pmax_dataframe.insert(i, sim_name + ' Max Power(mW/cm^2)', results[7])

            self.ui.Results_List_CHARGE.setItem(i, 0, QTableWidgetItem(str(sorted_CHARGE[i])))
            self.ui.Results_List_CHARGE.setItem(i, 1, QTableWidgetItem(str(results[3][0])))
            self.ui.Results_List_CHARGE.setItem(i, 2, QTableWidgetItem(str(results[4][0])))
            self.ui.Results_List_CHARGE.setItem(i, 3, QTableWidgetItem(str(results[5][0])))
            self.ui.Results_List_CHARGE.setItem(i, 4, QTableWidgetItem(str(results[6][0])))
            self.ui.Results_List_CHARGE.setItem(i, 5, QTableWidgetItem(str(results[7][0])))

            self.ui.Files_List_CHARGE.setItem(i, 1, QTableWidgetItem('Processed'))

        density_dataframe.to_excel('Density_Charge.xlsx', index=False)
        power_dataframe.to_excel('Power_Charge.xlsx', index=False)
        Voc_dataframe.to_excel('Voc_Charge.xlsx', index=False)
        Jsc_dataframe.to_excel('Jsc_Charge.xlsx', index=False)
        eta_dataframe.to_excel('Eff_Charge.xlsx', index=False)
        FF_dataframe.to_excel('FF_Charge.xlsx', index=False)
        Pmax_dataframe.to_excel('Pmax_Charge.xlsx', index=False)

        charge.close()


app = QtWidgets.QApplication([])
application = mywindow()
application.show()

sys.exit(app.exec())
