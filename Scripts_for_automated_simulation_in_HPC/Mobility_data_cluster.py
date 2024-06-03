############# Import required Libraries ########################
import pandas as pd
import numpy as np
from ovito.io import *
from ovito.modifiers import *
from ovito.data import *
from ovito.pipeline import *
import matplotlib.pyplot as plt
import tqdm
import pyvista as pv
import os  # for makedirs
import homcloud.interface as hc  # HomCloud 
from tqdm.notebook import tqdm  # For progressbar
import sklearn.linear_model as lm  # Machine learning
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
from sklearn.decomposition import PCA ,KernelPCA # for PCA
from sklearn.decomposition import FastICA  # for PCA
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler, StandardScaler

from sklearn.neighbors import KernelDensity
from scipy.signal import find_peaks
from scipy.stats import linregress
import pickle


######### Creat classes ###################

class GB_info:
    
    '''
    Given the GB type and simulation condition, this class identifies from a series of Lammps dump files:
    1: The loaction of GB at a every simulation step
    2: The loacation of atoms within a distant d  of the average GB loaction
    3: The relative displacement of top grain with respect to bottom grains
    
    Input:
        1:subfolder -> directory for the location of the simulation data
        2:name -> Name of the simulated grain boundary
        3:temp -> Temperature at which the simulation is performed
        4:p -> The applied driving force during the simulation
        5:p_type -> The type of applied driving force 
         ("constant_stress" for applied shear stress or  "constant_potential for jump in chemical potential")
        6:save_gb_atoms -> A flag to indicate if location of atoms within +- 15 Ang should be saved as numpy array
         (1 for yes or 0 for no)
        7:compute_pd -> A flag to indicate if peristant diagram is computed
         (1 for yes or 0 for no)
  
        
    '''
    def __init__(self, subfolder,name, temp,p,p_type,save_gb_atoms=1,compute_pd=1):
        self.subfolder=subfolder
        self.name = name
        self.temp= temp
        self.p=p
        self.p_mod=p
        self.p_type=p_type
        self.compute_pd=compute_pd
        self.save_gb_atoms=save_gb_atoms
    
        
        if p_type=="constant_stress":
            
            self.gb_name=f'{subfolder}{temp}/{p}/{p_type}/{name}/dump.{name}_shear_{Driving_force}_{temp}.wrap'
            self.save_dir=f'{subfolder}{temp}/{p}/{p_type}/{name}/PD'
        if p_type=="constant_potential":
            self.gb_name=f'{subfolder}{temp}/{p}/{p_type}/{name}/dump.{name}_potential_{p}_{temp}.wrap'
            self.save_dir=f'{subfolder}{temp}/{p}/{p_type}/{name}/PD'
        
    def GB_property_evaluation(self):
        step=[]
        N_atom=[]
        GB_location=[]
        Shear_x=[]
        Shear_y=[]
        
        pipeline = import_file(self.gb_name)
        N_frame=pipeline.source.num_frames
        os.makedirs(self.save_dir, exist_ok=True) 
        pipeline.modifiers.append(CentroSymmetryModifier())
        # Polyhedral template matching:
        pipeline.modifiers.append(PolyhedralTemplateMatchingModifier(
            rmsd_cutoff = 0.6, 
            output_rmsd = True, 
            output_interatomic_distance = True, 
            output_orientation = True))
        pipeline.modifiers.append(GrainSegmentationModifier())
        for frame in range(N_frame):
            data = pipeline.compute(frame)
            Sim_cell=np.array(data.cell)[0:3,0:3]
            Cell_origin=np.array(data.cell)[0:3,3:4]
            coords = data.particles.positions
            flag=np.array(data.particles['Periodic Image'])
            if self.p_type=="constant_potential":
                ui1=data.particles['f_synt[2]'][data.particles['Grain']==1]
                phi_u1=np.mean(0.5*(1-np.cos(np.pi*ui1)))
                ui2=data.particles['f_synt[2]'][data.particles['Grain']==2]
                phi_u2=np.mean(0.5*(1-np.cos(np.pi*ui2)))
                u_free=np.abs(phi_u1-phi_u2)
                self.p_mod=u_free
            #RMSD=data.particles['RMSD']
            CSP=data.particles['Centrosymmetry']
            #pipeline.modifiers[0].enabled = False
            #pipeline.modifiers[1].enabled = False
            
            #pipeline.modifiers.append(ExpressionSelectionModifier(expression = 'Position.Z > CellSize.Z-12.5 && Position.Z < CellSize.Z'))
            
            z=coords[:,2]
            imagex=flag[:,0]
            imagey=flag[:,1]
            box_top_location= Sim_cell[2,2]+Cell_origin[2]
            #print(box_top_location)
            box_bottom_location=Cell_origin[2]
            #print(box_bottom_location)
            is_top=np.where(np.logical_and(coords[:,2]>=box_top_location-12.5, coords[:,2]<=box_top_location))
            
            is_bottom=np.where(np.logical_and(coords[:,2]>=box_bottom_location, coords[:,2]<=box_bottom_location+12.5))
            x_top_unwrap=coords[:,0][is_top]+imagex[is_top]*Sim_cell[0,0]
            x_bottom_unwrap=coords[:,0][is_bottom]+imagex[is_bottom]*Sim_cell[0,0]
            
            y_top_unwrap=coords[:,1][is_top]+imagey[is_top]*Sim_cell[1,1]
            y_bottom_unwrap=coords[:,1][is_bottom]+imagey[is_bottom]*Sim_cell[1,1]
            
            Relative_displacement_x=np.mean(x_top_unwrap)-np.mean(x_bottom_unwrap)
            Relative_displacement_y=np.mean(y_top_unwrap)-np.mean(y_bottom_unwrap)
            
            rmsd=CSP[:]
            #rmsd1=rmsd.reshape(-1,1 )
            ban_diam=4.0699
            kde = KernelDensity(kernel='gaussian', bandwidth=ban_diam).fit(z.reshape(-1, 1),sample_weight=1+rmsd)
            s = np.linspace(min(z)+3*ban_diam,max(z)-3*ban_diam,200)
            e = np.exp(kde.score_samples(s.reshape(-1, 1)))
          
            peak_indices, peak_dict = find_peaks(e,height=np.mean(e))
            peak_heights = peak_dict['peak_heights']
            highest_peak_index = peak_indices[np.argmax(peak_heights)]
            s_max = s[highest_peak_index]
            
            GB_width_parameter=15
            GB_center=coords[np.where(np.logical_and(coords[:,2]>=s_max-GB_width_parameter, coords[:,2]<=s_max+GB_width_parameter))]
            N_atom.append(len(GB_center))
            GB_location.append(s_max)
            Shear_x.append(Relative_displacement_x)
            Shear_y.append(Relative_displacement_y)
            step.append(frame)
            
            if self.save_gb_atoms==1:
                np.save(f'{self.save_dir}/{frame}.npy', GB_center)
                
            if self.compute_pd==1:
                hc.PDList.from_alpha_filtration(GB_center, save_boundary_map=True, save_to=self.save_dir+"/{0}.pdgm".format(frame))
                
            print(f"completed frame {frame} for GB {self.name}")
            
            if np.abs(box_top_location-s_max) < ban_diam*5:
                print("Grain boundary reached top slab! Moving to another GB")
                break
            
            elif np.abs(box_bottom_location-s_max) < ban_diam*5:
                print("Grain boundary reached bottom slab! Moving to another GB")
                break
            
            
            
        self.N_atom=N_atom
        self.GB_location=GB_location
        self.Shear_x=Shear_x
        self.Shear_y=Shear_y
        self.step=step
        
        def find_normal_velocity(self,dt=5):
            x=np.array(self.step)*dt
            y=np.array(self.GB_location)
            slope_normal_velocity, intercept, r, p, se = linregress(x, y)
            self.normal_velocity=slope_normal_velocity
            self.mobility=slope_normal_velocity/self.p
            
        def find_shear_velocity_x(self,dt=5):
            x=np.array(self.step)*dt
            y=np.array(self.Shear_x)
            slope_shear_x, intercept, r, p, se = linregress(x, y)
            self.shear_velocity_x=slope_shear_x
            
        def shear_velocity_y(self,dt=5):
            x=np.array(self.step)*dt
            y=np.array(self.Shear_y)
            slope_shear_y, intercept, r, p, se = linregress(x, y)
            self.shear_velocity_y=slope_shear_y
        def find_coupling_factor_x(self):
            if self.shear_velocity_x==0:
                self.coupling_x="NAN"
            else:
                self.coupling_x=self.shear_velocity_x/self.normal_velocity
                
        def find_coupling_factor_y(self):
            if self.shear_velocity_y==0:
                self.coupling_y="NAN"
            else:
                self.coupling_x=self.shear_velocity_y/self.normal_velocity
            
            
            #hc.PDList.from_alpha_filtration(GB_center, save_boundary_map=True, save_to=save_dir+"/{0}.pdgm".format(frame))


############  Specify directory and GB type to construct persistant diagram ####################

Folder_name="D:/CGBS_and_E/Computed_GBE_and_Structure_Al/"
Subfolder=Folder_name+"Gen_GB_min_Final_109_all_simulations/"
Temp=500
Driving_force=500
Driving_force_type="constant_stress"
GB_index=0
compute_pd=1
save_gb_atoms=1
with open(Folder_name+'GB_names.txt') as gn:
    GB_names = [line.rstrip()[:-4] for line in gn][GB_index:GB_index+98]
    
    
Sim_results_potential = dict.fromkeys(GB_names, None)
Sim_results_stress = dict.fromkeys(GB_names, None)
### Directory GB migration simulation data for constant tau and Psi for given GB under specified condition  and save dir for PD###
for GB_type in  GB_names:    
    gb_name_stress=f'{Subfolder}{Temp}/constant_stress/{Driving_force}/{GB_type}/{GB_type}_shear_{Driving_force}_{Temp}.wrap'
    gb_name_potential=f'{Subfolder}{Temp}/constant_potential/{Driving_force}/{GB_type}/{GB_type}_potential_{Driving_force}_{Temp}.wrap'

    #save_dir_stress=f'{Subfolder}{Temp}/constant_stress/{Driving_force}/{GB_type}/PD'
    #save_dir_potential=f'{Subfolder}{Temp}/constant_potential/{Driving_force}/{GB_type}/PD'
    try:
        Sim=GB_info(Subfolder,GB_type,Temp,Driving_force,Driving_force_type,save_gb_atoms,compute_pd)
        Sim.GB_property_evaluation()
        if Driving_force_type=="constant_potential":
            Sim_results_potential[GB_type]=Sim
        if Driving_force_type=="constant_stress":
            Sim_results_stress[GB_type]=Sim
    except:
        continue
if Driving_force_type=="constant_potential":    
    with open(f'{Subfolder}{Temp}/constant_stress/{Driving_force}/Sim_results_potential.pkl','wb+') as f:
        pickle.dump(Sim_results_potential,f)
        
if Driving_force_type=="constant_stress":    
    with open(f'{Subfolder}{Temp}/constant_stress/{Driving_force}/Sim_results_stress.pkl','wb+') as f:
        pickle.dump(Sim_results_stress,f)

# for n,s in Sim_results_stress.items():
    
#      plt.scatter(s.step,s.GB_location)



### Import GB migration simulation data for constant tau and Psi for given GB under specified condition ###
