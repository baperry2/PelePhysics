import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# Plot cp and dcp/dT as a function of T for a specified mechanism from PelePhysics using Cantera
#
# Requires: Cantera, Numpy, Matplotlib

#### Inputs ####

mech = "LiDryer"
npts = 1001
Tmin = 300
Tmax = 2300
crossT = 1000
epsT = 1e-10
pres = 1.0 * ct.one_atm
specs = ["OH","H2O","H2","O2","N2","CO2","CO","CH4","HO2","H","O"]

#### Start of Script ####

plt.figure('normalized dcp/dT')
plt.clf()
plt.figure('normalized cp')
plt.clf()

gas = ct.Solution(mech+ "/mechanism.yaml")
specs = [spec for spec in specs if spec in gas.species_names]
for spec in specs:
    comp = spec + ":1.0"
    data = ct.SolutionArray(gas,npts)
    Temps = np.linspace(Tmin, Tmax, npts)
    data.TPY = Temps, pres, comp
    cp = data.cp
    dcpdT = np.diff(cp)/np.diff(Temps)
    plt.figure('normalized cp')
    plt.plot(Temps, cp/np.max(cp), label=spec)
    plt.figure('normalized dcp/dT')
    plt.plot(Temps[1:], dcpdT/np.max(np.abs(dcpdT)), label=spec)
    plt.legend(frameon=False)

for spec in gas.species_names:
    comp = spec + ":1.0"
    gas.TPY = crossT - epsT, pres, comp
    cp_low = gas.cp
    h_low = gas.h
    gas.TPY = crossT + epsT, pres, comp
    cp_high = gas.cp
    h_high = gas.h
    print(f"{spec} has C0 relative discontinuity of {np.abs(cp_high-cp_low)/cp_low}")
    print(f"{spec} has C0 relative discontinuity of {np.abs(h_high-h_low)/h_low}")

plt.figure('normalized cp')
plt.xlabel('T (K)')
plt.ylabel('normalized cp')
plt.legend(frameon=False)
plt.tight_layout()

plt.figure('normalized dcp/dT')
plt.xlabel('T (K)')
plt.ylabel('normalized dcp/dT')
plt.legend(frameon=False)
plt.tight_layout()

plt.show()
