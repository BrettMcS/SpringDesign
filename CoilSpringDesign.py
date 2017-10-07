#CoilSpringDesign.py
"""
Design of Springs to EN 13906-1:2014.

"""
import numpy as np
from scipy.optimize import root
import helicoil as hc


class TwoCoilSetLength:
    """
    A two-coil spring nest with a defined design lengh.
    """
    def __init__(self, material, axial_rate, design_load, design_length,
                 radial_coil_gap, max_compression,
                 lo_cycle_defln_amplitude, hi_cycle_defln_amplitude, 
                 end_condition, compression_defln_reserve, solid_stress_reserve,
                 estimated_solution):
        """
        """
        self.mat = material
        self.end_con = end_condition
        self.R = axial_rate
        self.F = design_load
        self.L = design_length
        self.coil_gap = radial_coil_gap
        self.max_comp = max_compression
        self.lo_cycle_amp = lo_cycle_defln_amplitude
        self.hi_cycle_amp = hi_cycle_defln_amplitude
        self.comp_res = compression_defln_reserve
        self.solid_str_res = solid_stress_reserve
        self.x0 = estimated_solution


    def error(self, x):
        """
        Return the vector of residuals of the target solution of a 2-coil design

        :Outline:
            A two-coil design with a defined axial rate and design length.
            The routine targets a solution where the nest has the target
            minimum length, the solid stress reserves of the two coils are equal,
            and the inner coil has a target stability value.
            The free lengths may be unequal.

        :Parameters:
            x: np array, independent variables [OD,do,No,L0o,di,Ni,L0i]
            
        :Returns:
            v: np array, residual from target solution
        """
        OD,do,No,L0o,di,Ni,L0i = x[:]

        min_service_len = self.L - self.max_comp  # min length in service
        min_length = min_service_len - self.comp_res # the target min length   
        Do = OD - do
        Di = Do - do - 2*self.coil_gap - di

        Lco = hc.solid_length(No, do)
        Lci = hc.solid_length(Ni, di)

        Sao = hc.Sa_min_reserve_length(Do, do, No)
        Sai = hc.Sa_min_reserve_length(Di, di, Ni)
        
        minLength_o = Lco + Sao
        minLength_i = Lci + Sai

        Ro = hc.axial_rate(self.mat.G, do, Do, No)
        Ri = hc.axial_rate(self.mat.G, di, Di, Ni)

        if L0o > L0i:
            Fo1 = Ro*(L0o - L0i)
            Fi1 = 0.0
            L1 = L0i
        else:
            Fo1 = 0.0
            Fi1 = Ri*(L0i - L0o)
            L1 = L0o
            
        F2 = self.F - Fo1 - Fi1
        L_at_F = L1 - F2/(Ro + Ri)

        solidStr_o = hc.axial_stress_static(Do, do, Ro*(L0o-Lco))
        solidStr_i = hc.axial_stress_static(Di, di, Ri*(L0i-Lci))

        solidStrRes_o = self.mat.solid_stress_limit(do) - solidStr_o
        solidStrRes_i = self.mat.solid_stress_limit(di) - solidStr_i

        buck_defln = hc.buckling_deflection(self.mat.G, self.mat.E, 
                                            Di, L0i, self.end_con)
        buckling_len = L0i - buck_defln

        # Calculate the errors in the solution:
        v = np.zeros(len(x))

        v[0] = Ro + Ri - self.R  # total rate to be R
        v[1] = minLength_o - min_length
        v[2] = minLength_i - min_length
        v[3] = L_at_F - self.L
        v[4] = solidStrRes_o - self.solid_str_res
        v[5] = solidStrRes_i - self.solid_str_res
        v[6] = buckling_len - min_service_len
        
        return v


    def get_solution(self):
        """
        """
        sol = root(self.error, self.x0)

        if not sol.success:
            return False, sol, None

        OD,do,No,L0o,di,Ni,L0i = sol.x

        Do = OD - do
        Di = Do - do - 2*self.coil_gap - di

        Ro = hc.axial_rate(self.mat.G, do, Do, No)
        Ri = hc.axial_rate(self.mat.G, di, Di, Ni)

        if L0o > L0i:
            Fo1 = Ro*(L0o - L0i)
            Fi1 = 0.0
            L1 = L0i
        else:
            Fo1 = 0.0
            Fi1 = Ri*(L0i - L0o)
            L1 = L0o

        Fo = Fo1 + Ro*(L1 - self.L)
        Fi = Fi1 + Ri*(L1 - self.L)

        L_min = self.L - self.max_comp

        OC_data = hc.coil_data("Outer Coil", do, Do, No, Fo, self.L, L_min, 
                            self.lo_cycle_amp, self.hi_cycle_amp, self.mat)

        IC_data = hc.coil_data("Inner Coil", di, Di, Ni, Fi, self.L, L_min, 
                            self.lo_cycle_amp, self.hi_cycle_amp, self.mat)

        return True, OC_data, IC_data


if __name__ == "__main__":

    # Specify the spring parameters.See also CLN.1615

    nest_data = dict(
        material = hc.prEN10089,
        axial_rate = 280.0,  # N/mm Axial rate
        design_load = 34760.0,  # N at maximum services
        design_length = 367.8,  # mm. Spring modelled at maximum service load
        radial_coil_gap = 8.0,  # mm. radial gap, inner to outer coils
        max_compression = 45.0,  # to check clearance from the solid length
        lo_cycle_defln_amplitude = 35.0,
        hi_cycle_defln_amplitude = 25.0,
        end_condition = 0.7,  # the spring end condition.
        compression_defln_reserve = 20.0,  # mm
        solid_stress_reserve = 70.0,  # MPa
        estimated_solution = [200.0, 30.0, 8.0, 500.0, 20.0, 11.0, 500.0],
    )    #                     OD     do   No    L0o    di    Ni    L0i

    nest = TwoCoilSetLength(**nest_data)

    solved,OC_data,IC_data = nest.get_solution()

    if solved:
        open("test_Outer-coil.txt", 'w').write(hc.coil_data_csv(OC_data))
        open("test_Inner-coil.txt", 'w').write(hc.coil_data_csv(IC_data))
        print("\nsolution written")
    else:
        print("No solution found. Solver data returned:")
        print(OC_data)

    print("Solid stress reserve, IC coil gap/do,  OC coil gap/di")
    for solid_stress_reserve in np.arange(10.0, 100.0, 10.0):
        nest_data["solid_stress_reserve"] = solid_stress_reserve
        nest = TwoCoilSetLength(**nest_data)
        solved,OC_data,IC_data = nest.get_solution()
        OC_gap_ratio = OC_data.coil_gap / OC_data.d
        IC_gap_ratio = IC_data.coil_gap / IC_data.d
        print(f"{solid_stress_reserve:7.1f}, {OC_gap_ratio:7.2f}, {IC_gap_ratio:7.2f}")