import meep as mp
import numpy as np
import math
import os
from grating_antenna_gds import waveguide_antenna


def simulation(sim, phis, thetas,fcen, npts, sx, sy, sz, dpml, dair,r, case):
    nearfield = sim.add_near2far(fcen, 0, 1, mp.Near2FarRegion(mp.Vector3(0,0.5*sy-dpml-2.5*dair,0), size=mp.Vector3(sx-2*dpml, 0, sz-2*dpml)))

    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ex, mp.Vector3(), 1e-3))
    
    angles = phis

    E = np.zeros((npts,3),dtype=np.complex128)
    H = np.zeros((npts,3),dtype=np.complex128)
    for n in range(npts):
        ff = sim.get_farfield(nearfield,
                          mp.Vector3(r*math.cos(angles[n]),
                                     r*math.sin(angles[n]),0))
        E[n,:] = [ff[j] for j in range(3)]
        H[n,:] = [ff[j+3] for j in range(3)]
    
    Px = np.real(np.conj(E[:, 1]) * H[:, 2] - np.conj(E[:, 2]) * H[:, 1])
    Py = np.real(np.conj(E[:, 2]) * H[:, 0] - np.conj(E[:, 0]) * H[:, 2])
    Pz = np.real(np.conj(E[:, 0]) * H[:, 1] - np.conj(E[:, 1]) * H[:, 0])
    Pr = np.sqrt(np.square(Px) + np.square(Pz))
    Pr = Pr/np.max(Pr)
    idx = np.where((phis > np.deg2rad(45)) & (phis < np.deg2rad(70)))
    sPr = Pr[idx]
    val = (np.sum(sPr))/(np.sum(Pr))
    
    return val


def setup(wg, h, case):
    grid_step = 0.002  # example resolution unit
    wg = round(wg / grid_step) * grid_step
    h = round(h / grid_step) * grid_step
    # p = round(p / grid_step) * grid_step
    # d = round(d / grid_step) * grid_step

    
    dpml = 1.0
    dair = 1

    t_si3n4 = 1
    t_si = 0.45
    t_sio2 = 3*(t_si3n4)
    t_sub = t_sio2
    
    resolution = 10


    si = mp.Medium(index=3.6)
    sio2 = mp.Medium(index=1.44)
    si3n4 = mp.Medium(index=2.0)

    wvl = 1.55
    fcen = 1 / wvl
    width = 0.1
    fwidth = width * fcen

    src = mp.GaussianSource(frequency=fcen, fwidth=fwidth)

    r = 1000/fcen  # 1000 wavelengths out from the source
    npts = 100  # number of points in [0,2*pi) range of angles
    # define field projection angle ranges
    thetas = np.linspace(0, np.pi/3, npts)  # polar angle
    phis = np.linspace(0, 2*np.pi, npts)  # azimuthal angle
  
    l_g, w, gdsII_file = waveguide_antenna(wg, h)

    geometry = (
        mp.get_GDSII_prisms(si3n4, gdsII_file, 4, -t_si3n4/2, t_si3n4/2)+
        mp.get_GDSII_prisms(si, gdsII_file, 3, -t_si/2, t_si/2) +
        mp.get_GDSII_prisms(sio2, gdsII_file, 2, -t_sio2/2, t_sio2/2) +
        mp.get_GDSII_prisms(si, gdsII_file, 1, -t_sub/2, t_sub/2)
    
    )
    sx, sy, sz = l_g + 2 * dpml, w + 2 * dpml + 3.5 * dair, t_sub + 2 * dpml
    cell = mp.Vector3(sx, sy, sz)

    for g in geometry:
        #w0+w1 = 4

        g.center -= mp.Vector3(l_g/2, 4, 0)
        
    #define the source
    #size = (0, waveguide width, 0)
    sources = [mp.EigenModeSource(src=mp.GaussianSource(fcen, fwidth=0.2*fcen),
                                size=mp.Vector3(0,wg,t_sub),
                                center=mp.Vector3(-0.5*sx + dpml,wg/2,0),
                                eig_match_freq=True,
                                eig_parity=mp.ODD_Y)]

    # Use PML boundary layers on all sides
    pml_layers = [mp.PML(dpml)]
    
    # Create simulation object
    sim = mp.Simulation(
        cell_size=cell,
        boundary_layers=pml_layers,
        geometry=geometry,
        sources=sources,
        resolution=resolution,
        dimensions=3
    )
    #to visualize the setup
    val = simulation(sim, phis, thetas,fcen, npts, sx, sy, sz, dpml, dair, r, case)
    os.remove(gdsII_file)
    return val


def main(input, grad):
    mp.verbosity(0)
    try:
        wg = input[0]
        h = input[1]
        # p = input[2]
        # d = input[3]
        val = setup(wg, h, "noplot")
        print(f"Evaluated at wg={wg}, h={h}, val={val}")
        
        return val
    except Exception as e:
        print(f"Error during evaluation: {e}")
        # Return a very low value so NLopt avoids this area
        return 1e-6

    return val
    
        

