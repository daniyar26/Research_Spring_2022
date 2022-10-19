# -*- coding: utf-8 -*-
"""
This library is used to analyze the mechanical properties and energy absorption
characteristics of the 3D printed lattice structures. 

Previously similar code was used to analyze plastic lattices. This is updated
version that uses OOP for the analysis and is more efficient

The code is open and can be edited. No restrictions are imposed on its
distribution and/or use. 

All rights are reserved

For any query please contact: Daniyar Syrlybayev
Email: daniyar.syrlybayev@nu.edu.kz
Phone/WhatsApp: +77054955819
Telegram: @sdrfem

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def read_xlsx(path):
    data = pd.read_excel(path, 
                     skiprows = [0], usecols=lambda x: x if 
                     not x.startswith('Unnamed') else None)
    return data

class lattice_structure:
    def __init__(self, path, dimen, rel_dens, height, mass, label):
        self.data = read_xlsx(path)
        self.area = dimen[1] * dimen[0]
        self.rel_dens = rel_dens
        self.height = height
        self.props = {'Lattice': label,
                      'Mass (g)': mass,
                      'Relative Density': rel_dens}
    
    def __str__(self):
        return str(self.props)
    
    def stress_strain(self):
        self.stress, self.strain = (np.array(self.data["Load"])/self.area, 
            np.array(self.data["Displ."])/self.height)
    
    def find_strengh(self, limit = 0.6):
        k = list(abs(self.strain - limit)).index(min(abs(self.strain - limit)))
        list_stress = list(self.stress)
        self.props['Strengh'] = max(list_stress[0 : k])
        ind = list_stress.index(self.props['Strengh'])
        self.props['Disp'] = self.strain[ind]
    
    
    def find_elastic_modulus(self, limit1, limit2):
        ind1 = list(abs(self.strain - limit1)).index(min(
            abs(self.strain - limit1)))
        ind2 = list(abs(self.strain - limit2)).index(min(abs(
            self.strain - limit2)))
        
        y = self.stress[ind1:ind2]
        x = self.strain[ind1:ind2]
        
        self.f = np.polyfit(x, y, 1)
        self.expression = np.poly1d(self.f)
        self.props['Elastic Modulus'] = self.f[0]
        
    def energy_absorption_efficiency(self, frequency = 50):
        self.energy_absorption_eff = []
        self.new_strain1 = []

        for i in range(1, self.strain.shape[0]):
            if i % frequency != 1:
                continue
            
            en_abs = np.trapz(self.stress[0:i], self.strain[0:i])/1000
            efficiency = en_abs/(max(self.stress[0:i])) * 100
            self.energy_absorption_eff.append(efficiency)
            self.new_strain1.append(self.strain[i])
            
        
        self.props['Maximum Energy Absorption Efficiency'] = max(
            self.energy_absorption_eff)
        ind = self.energy_absorption_eff.index(self.props[
            'Maximum Energy Absorption Efficiency'])
        self.props['Densification Strain'] = self.new_strain1[ind]
        self.props['Maximum Energy Absorption Efficiency'] /= (
            self.props['Densification Strain'])
        self.energy_absorption_eff /= self.props['Densification Strain']
    
    def energy_aborption(self, frequency = 50):
        self.energy_absorption = []
        self.new_strain = []
        ind = list(self.strain - self.props['Densification Strain']
                   ).index(min(abs(self.strain - self.props[
                       'Densification Strain'])))
        
        for i in range(1, ind):
            if i % frequency != 1:
                continue
            
            en_abs = np.trapz(self.stress[0:i], self.strain[0:i])/1000
            self.energy_absorption.append(en_abs)
            self.new_strain.append(self.strain[i])
        
        self.props['Maximum Energy Absorption'] = max(self.energy_absorption)
    
    def plateau_stress(self):
        dens_index = list(self.strain - self.props['Densification Strain']
                   ).index(min(abs(self.strain - self.props[
                       'Densification Strain'])))
        displ_index = list(self.strain - self.props['Disp']
                   ).index(min(abs(self.strain - self.props[
                       'Disp'])))
        plat_stress = np.trapz(self.stress[displ_index:dens_index], 
                               self.strain[displ_index:dens_index])/(
                                   self.props['Densification Strain'] -
                                   self.props['Disp'])
        self.props['Plateau Stress'] = plat_stress
        
    def specific_props(self):
        not_in = ['Disp', 'Maximum Energy Absorption Efficiency', 
                  'Densification Strain', 'Mass (g)', 'Relative Density', 
                  'Lattice']
        params = [i for i in self.props.keys() if i not in not_in][:]
        for i in params:
            self.props['Specific ' + i] = (self.props[i]/
                                          self.props['Mass (g)'])


def process_lattice(path, dimen, rel_dens, height, mass, label, limit1, 
                    limit2):
    lattice = lattice_structure(path, dimen, rel_dens, height, mass, label)
    lattice.stress_strain()
    lattice.find_strengh(limit = 0.3)
    lattice.find_elastic_modulus(limit1, limit2)
    lattice.energy_absorption_efficiency()
    lattice.energy_aborption()
    lattice.plateau_stress()
    return lattice


def plot_stress_strain_curve(*lattice:[lattice_structure], 
                             savepath, save_name, labels):
    fig, ax = plt.subplots()
    ax.set_xlabel('Strain $\epsilon$')
    ax.set_ylabel('Stress $\sigma$ $MPa$')
 
    for i in range(len(lattice)):
        ax.plot(lattice[i].strain, lattice[i].stress, label = labels[i])
        ax.plot(lattice[i].props['Disp'], lattice[i].props['Strengh'], 
                    marker = 'x')
    plt.legend(ncol = 2, prop={'size': 8})
    fig.savefig(savepath + '//' + save_name + '.pdf')
    plt.show()
    plt.close('all')
    return fig, ax

def plot_energy_absorption_curve(*lattice, savepath, save_name, labels):
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('Strain $\epsilon$')
    ax1.set_ylabel('Energy Absorption $\psi$')
    for i in range(len(lattice)):
        ax1.plot(lattice[i].new_strain, lattice[i].energy_absorption, 
                label = labels[i])
    plt.legend(ncol = 2, prop={'size': 8})
    fig.savefig(savepath + '//' + save_name + '.pdf')
    plt.show()
    plt.close('all')
    return fig, ax1

def plot_energy_absorption_eff_curve(*lattice, savepath, save_name, labels):
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('Strain $\epsilon$')
    ax1.set_ylabel('Energy Absorption Efficiency $\eta$, %')
    for i in range(len(lattice)):
        ax1.plot(lattice[i].new_strain1, lattice[i].energy_absorption_eff, 
                label = labels[i])
    plt.legend(ncol = 2, prop={'size': 8})
    fig.savefig(savepath + '//' + save_name + '.pdf')
    plt.show()
    plt.close('all')
    return fig, ax1
        
     
# #-------------------------------------Test------------------------------------

# A = process_lattice(r"C:\Users\Daniyar Syrlybayev\Desktop\Masters Thesis\Gyroid\Gyroid_0_25\Samp1_Gyoid_0_25.xlsx", 
#             (30, 30), 0.25, 30, 19.26, 'Gyroid', limit1 = 0.024, limit2 = 0.025)
# plot_stress_strain_curve(A,
#     savepath = r'C:\Users\Daniyar Syrlybayev\Desktop\Masters Thesis\Gyroid\Gyroid_0_25', 
#     save_name = 'Gyroid', labels = ['Gyroid'])

# plot_energy_absorption_curve(A,
#     savepath = r'C:\Users\Daniyar Syrlybayev\Desktop\Masters Thesis\Gyroid\Gyroid_0_25', 
#     save_name = 'Gyroid_en_abs', labels = ['Gyroid'])

# plot_energy_absorption_eff_curve(A,
#     savepath = r'C:\Users\Daniyar Syrlybayev\Desktop\Masters Thesis\Gyroid\Gyroid_0_25', 
#     save_name = 'Gyroid_en_abs_eff', labels = ['Gyroid'])


# print(A)


    