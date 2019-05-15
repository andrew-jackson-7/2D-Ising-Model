import numpy as np
import scipy
import matplotlib
import sys
import math
import time as tm

filename = 'ParametersFile.txt'

file = open(filename, 'r')
parameters = file.read()
parasplit = parameters.split('\n')

print(parasplit)

NumParaVals = [float(col) for col in parasplit]

print(NumParaVals)

lattice_length = int(NumParaVals[0])

J = float(NumParaVals[1])

H = float(NumParaVals[2])

delta_H = float(NumParaVals[3])

kBT = 0.5

class Lattice(object):
       
    def __init__(self, lattice_length):

        self.lattice_length = lattice_length
    
        self.lattice = np.zeros((self.lattice_length, self.lattice_length))
    
        self.Bonds = {}

    def lattice_initialize(self):

        for i in range(self.lattice_length):

            for j in range(self.lattice_length):

                if round(np.random.uniform(0.0, 1.0), 2) < 0.5:

                    self.lattice[i][j] = -1.0

                else:

                    self.lattice[i][j] = 1.0

    def bonds_initialize(self):

        k = 0
        
        for o in range(2):

            for i in range(1, self.lattice_length + 1):

                for j in range(1, self.lattice_length + 1):
    
                    if k <= (self.lattice_length ** 2) - 1:
    
                        if self.lattice_length % j == 0:

                            self.Bonds["bond " + str(k)] = np.array([self.lattice[i - 1][j - 1], self.lattice[i - 1][self.lattice_length - j]])
            
                            k = k + 1
            
                        else:
            
                            self.Bonds["bond " + str(k)] = np.array([self.lattice[i - 1][j - 1], self.lattice[i - 1][j]])
            
                            k = k + 1
            
                    else:

                        if self.lattice_length % i == 0:

                            self.Bonds["bond " + str(k)] = np.array([self.lattice[i - 1][j - 1], self.lattice[self.lattice_length - i][j - 1]])

                            k = k + 1

                        else:
            
                            self.Bonds["bond " + str(k)] = np.array([self.lattice[i - 1][j - 1], self.lattice[i][j - 1]])

                            k = k + 1
    
    def get_lattice_site(self, i, j):

        return self.lattice[i - 1][j - 1]
    
    def set_lattice_site(self, i, j, value):
    
        self.lattice[i - 1][j - 1] = value
        
        if (i == 1 and j != 1):
        
            self.Bonds["bond " + str((i - 1) + (j - 1))][0] = value
        
            self.Bonds["bond " + str((i - 1) + (j - 2))][1] = value
        
            self.Bonds["bond " + str((self.lattice_length ** 2) + (j - 1))][0] = value
        
            self.Bonds["bond " + str(((2 * self.lattice_length - 1) * self.lattice_length) + (j - 1))][1] = value
        
        elif (i != 1 and j == 1):
        
            self.Bonds["bond " + str((i - 1) * self.lattice_length + (j - 1))][0] = value
        
            self.Bonds["bond " + str((i - 1) * self.lattice_length + (self.lattice_length - 1))][1] = value
        
            self.Bonds["bond " + str((i + self.lattice_length - 1) * self.lattice_length + (j - 1))][0] = value
        
            self.Bonds["bond " + str((i + self.lattice_length - 2) * self.lattice_length + (j - 1))][1] = value
        
        elif (i == 1 and j == 1):
        
            self.Bonds["bond " + str((i - 1) + (j - 1))][0] = value
        
            self.Bonds["bond " + str((i - 1) + (self.lattice_length - 1))][1] = value
        
            self.Bonds["bond " + str((i + self.lattice_length - 1) * self.lattice_length + (j - 1))][0] = value
        
            self.Bonds["bond " + str(((2 * self.lattice_length - 1) * self.lattice_length) + (j - 1))][1] = value
            
        else:
    
            self.Bonds["bond " + str((i - 1) * self.lattice_length + (j - 1))][0] = value
    
            self.Bonds["bond " + str((i - 1) * self.lattice_length + (j - 2))][1] = value
    
            self.Bonds["bond " + str((i + self.lattice_length - 1) * self.lattice_length + (j - 1))][0] = value
    
            self.Bonds["bond " + str((i + self.lattice_length - 2) * self.lattice_length + (j - 1))][1] = value

    def get_bond(self, i):

        return self.Bonds["bond " + str(i - 1)]

    def sum_of_sites(self):
    
        Site_Term_Total = 0.0

        for i in range(self.lattice_length):

            for j in range(self.lattice_length):

                Site_Term_Total = Site_Term_Total + self.lattice[i][j]

        return Site_Term_Total

    def sum_of_NN(self):

        NN_Term_Total = 0.0

        for i in range(0, len(self.Bonds)):

            NN_Term_Total = NN_Term_Total + (self.Bonds.get("bond " + str(i))[0] * self.Bonds.get("bond " + str(i))[1])

        return NN_Term_Total

    def return_lattice(self):

        return self.lattice

    def return_bonds(self):

        return self.Bonds
        
    @staticmethod
    def sum_total_NN_bonds(array_of_sites, latticelength, i, j):
    
        if (i == 1 and j != 1):
        
            return (array_of_sites["bond " + str((i - 1) + (j - 1))][0] * array_of_sites["bond " + str((i - 1) + (j - 1))][1]) + (array_of_sites["bond " + str((i - 1) + (j - 2))][0] * array_of_sites["bond " + str((i - 1) + (j - 2))][1]) + (array_of_sites["bond " + str((latticelength ** 2) + (j - 1))][0] * array_of_sites["bond " + str((latticelength ** 2) + (j - 1))][1]) + (array_of_sites["bond " + str(((2 * latticelength - 1) * latticelength) + (j - 1))][0] * array_of_sites["bond " + str(((2 * latticelength - 1) * latticelength) + (j - 1))][1])
        
        elif (i != 1 and j == 1):
        
            return (array_of_sites["bond " + str((i - 1) * latticelength + (j - 1))][0] * array_of_sites["bond " + str((i - 1) * latticelength + (j - 1))][1]) + (array_of_sites["bond " + str((i - 1) * latticelength + (latticelength - 1))][0] * array_of_sites["bond " + str((i - 1) * latticelength + (latticelength - 1))][1]) + (array_of_sites["bond " + str((i + latticelength - 1) * latticelength + (j - 1))][0] * array_of_sites["bond " + str((i + latticelength - 1) * latticelength + (j - 1))][1]) + (array_of_sites["bond " + str((i + latticelength - 2) * latticelength + (j - 1))][0] * array_of_sites["bond " + str((i + latticelength - 2) * latticelength + (j - 1))][1])
        
        elif (i == 1 and j == 1):
        
            return (array_of_sites["bond " + str((i - 1) + (j - 1))][0] * array_of_sites["bond " + str((i - 1) + (j - 1))][1]) + (array_of_sites["bond " + str((i - 1) + (latticelength - 1))][0] * array_of_sites["bond " + str((i - 1) + (latticelength - 1))][1]) + (array_of_sites["bond " + str((i + latticelength - 1) * latticelength + (j - 1))][0] * array_of_sites["bond " + str((i + latticelength - 1) * latticelength + (j - 1))][1]) + (array_of_sites["bond " + str(((2 * latticelength - 1) * latticelength) + (j - 1))][0] * array_of_sites["bond " + str(((2 * latticelength - 1) * latticelength) + (j - 1))][1])
            
        else:
        
            return (array_of_sites["bond " + str((i - 1) * latticelength + (j - 1))][0] * array_of_sites["bond " + str((i - 1) * latticelength + (j - 1))][1]) + (array_of_sites["bond " + str((i - 1) * latticelength + (j - 2))][0] * array_of_sites["bond " + str((i - 1) * latticelength + (j - 2))][1]) + (array_of_sites["bond " + str((i + latticelength - 1) * latticelength + (j - 1))][0] * array_of_sites["bond " + str((i + latticelength - 1) * latticelength + (j - 1))][1]) + (array_of_sites["bond " + str((i + latticelength - 2) * latticelength + (j - 1))][0] * array_of_sites["bond " + str((i + latticelength - 2) * latticelength + (j - 1))][1])


lattice_1 = Lattice(lattice_length)

lattice_1.lattice_initialize()

lattice_1.bonds_initialize()

Energy_Config = {}

NN_Term_Total = lattice_1.sum_of_NN()

Site_Term_Total = lattice_1.sum_of_sites()

Energy_Config["EnergyOfConfig_0"] = -J * (NN_Term_Total) - H * (Site_Term_Total)

q = 1

for i in range(1, 500000):

    random_site_i = np.random.randint(1, lattice_length + 1)
    
    random_site_j = np.random.randint(1, lattice_length + 1)
    
    site_value = lattice_1.get_lattice_site(random_site_i, random_site_j)
    
    bonds_array = lattice_1.return_bonds()
    
    Deducted_NN_Term_Total = NN_Term_Total - Lattice.sum_total_NN_bonds(bonds_array, lattice_length, random_site_i, random_site_j)
    
    Deducted_Site_Term_Total = Site_Term_Total - site_value
    
    if (site_value < 0):

        lattice_1.set_lattice_site(random_site_i, random_site_j, 1.0)

        Proposed_NN_Term_Total = Deducted_NN_Term_Total + Lattice.sum_total_NN_bonds(bonds_array, lattice_length, random_site_i, random_site_j)

        Proposed_Site_Term_Total = Deducted_Site_Term_Total + lattice_1.get_lattice_site(random_site_i, random_site_j)

        Proposed_Energy_Config = -J * (Proposed_NN_Term_Total) - H * (Proposed_Site_Term_Total)

    else:

        lattice_1.set_lattice_site(random_site_i, random_site_j, -1.0)

        Proposed_NN_Term_Total = Deducted_NN_Term_Total + Lattice.sum_total_NN_bonds(bonds_array, lattice_length, random_site_i, random_site_j)

        Proposed_Site_Term_Total = Deducted_Site_Term_Total + lattice_1.get_lattice_site(random_site_i, random_site_j)

        Proposed_Energy_Config = -J * (Proposed_NN_Term_Total) - H * (Proposed_Site_Term_Total)

    if (Proposed_Energy_Config - Energy_Config["EnergyOfConfig_" + str(q - 1)] < 0):

        Energy_Config["EnergyOfConfig_" + str(q)] = Proposed_Energy_Config

        NN_Term_Total = Proposed_NN_Term_Total

        Site_Term_Total = Proposed_Site_Term_Total

        q = q + 1

    else:

        temp_random_value = np.random.uniform(low = 0.0, high = 1.0 + sys.float_info.epsilon)

        if (temp_random_value <= math.exp(-(1.0 / kBT) * (Proposed_Energy_Config - Energy_Config["EnergyOfConfig_" + str(q - 1)]))):

            Energy_Config["EnergyOfConfig_" + str(q)] = Proposed_Energy_Config
        
            NN_Term_Total = Proposed_NN_Term_Total
        
            Site_Term_Total = Proposed_Site_Term_Total
        
            q = q + 1

        else:

            if (lattice_1.get_lattice_site(random_site_i, random_site_j) == 1.0):
            
                lattice_1.set_lattice_site(random_site_i, random_site_j, -1.0)
                
            else:
            
                lattice_1.set_lattice_site(random_site_i, random_site_j, 1.0)

filename = 'Energy_Config_Array.txt'

file = open(filename, 'w')

p = 0

for item in Energy_Config:

    if ("EnergyOfConfig_" + str(p)) in Energy_Config:
    
        file.write('%s, %s\n' % ("EnergyOfConfig_" + str(p), Energy_Config["EnergyOfConfig_" + str(p)]))

    p = p + 1

file.close()
