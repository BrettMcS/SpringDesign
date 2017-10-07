# -*- coding: utf-8 -*-
"""
PR20 Primary Spring Nest Design.
"""
import numpy as np
import helicoil as hc
import CoilSpringDesign


def ODerror(solid_stress_reserve, ODtarget, nest_data):
    """
    """
    nest_data["solid_stress_reserve"] = solid_stress_reserve

    nest = CoilSpringDesign.TwoCoilSetLength(**nest_data)

    solved,OC_data,IC_data = nest.get_solution()

    if solved:
        return OC_data.OD - ODtarget




# Specify the spring parameters.See also CLN.1615

nest_data = dict(
    material = hc.prEN10089,
    end_condition =    0.7,  # the spring end condition.
    axial_rate =     280.0,  # N/mm Axial rate
    design_load =  34760.0,  # N at maximum services
    design_length =  367.8,  # mm. Spring modelled at maximum service load
    radial_coil_gap =  8.0,  # mm. radial gap, inner to outer coils
    max_compression = 45.0,  # to check clearance from the solid length
    lo_cycle_defln_amplitude = 35.0,
    hi_cycle_defln_amplitude = 25.0,
    compression_defln_reserve = 20.0,  # mm
    solid_stress_reserve = 70.0,  # MPa
    estimated_solution = [200.0, 30.0, 8.0, 500.0, 20.0, 11.0, 500.0])
#                          OD     do   No    L0o    di    Ni    L0i

nest = CoilSpringDesign.TwoCoilSetLength(**nest_data)

solved,OC_data,IC_data = nest.get_solution()

if solved:
    open("PR20_Outer-coil.txt", 'w').write(hc.coil_data_csv(OC_data))
    open("PR20_Inner-coil.txt", 'w').write(hc.coil_data_csv(IC_data))
    print("\nsolution written")
else:
    print("No solution found. Solver data returned:")
    print(OC_data)

print("Solid stress reserve, Nest OD,  Free length OC-IC")
for solid_stress_reserve in np.arange(10.0, 120.0, 10.0):
    nest_data["solid_stress_reserve"] = solid_stress_reserve
    nest = CoilSpringDesign.TwoCoilSetLength(**nest_data)
    solved,OC_data,IC_data = nest.get_solution()
    freeDel = OC_data.L0 - IC_data.L0
    print(f"{solid_stress_reserve:7.1f}, {OC_data.OD:7.2f}, {freeDel:7.2f}")