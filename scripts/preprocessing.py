import ehtim as eh
import numpy as np
import matplotlib.pyplot as plt
import ehtim.scattering.stochastic_optics as so
import sys
import os

#uvfits_path = './sgra_standardized_uvfits_hops/uvfits/3599/hops/hops_3599_SGRA_LO_netcal_LMTcal_10s.uvfits'
#uvfits_path = './ER6/3.netcal/hops_3601_SGRA_LO_netcal_LMTcal_10s.uvfits'
uvfits_path = sys.argv[1]

def generate_scans(obs, averaging_mode, timeres=60.):

    if averaging_mode == 'scan-average':
        obs.add_scans(dt=2*timeres/3600.)
        obs = obs.avg_coherent(timeres)
        obs = obs.avg_coherent(0., scan_avg=True)
    elif averaging_mode == 'scan-subaverage':
        obs.add_scans(dt=2*timeres/3600.)
        obs = obs.avg_coherent(timeres)
    elif averaging_mode == 'snapshot':
        obs = obs.avg_coherent(timeres)
        obs.add_scans(dt=0.)
    else:
        raise ValueError('Averaging mode not recognized.')

    # Only include data points with at least 80% of averaging time  

    todel = []
    for i in range(len(obs.data)):
        if obs.data['tint'][i] < 0.8*timeres:
            todel.append(i)
    obs.data = np.delete(obs.data, todel)

    
    # Plot to check scans
    #nscans = len(obs.tlist(scan_gather=True))
    #print(nscans)
    #for i in range(nscans):
    #    obs_scan = obs.copy()
    #    obs_scan.data = obs_scan.tlist(scan_gather=True)[i]
    #    plt.scatter(obs_scan.data['time'],np.zeros(len(obs_scan.data)))
    #plt.show()
    #input('...')

    return obs

def add_noise_deblur(obs, scan_number, noise_frac=0.05, refr_floor_frac=0.005, zbl=2.):

    obs_scan = obs.copy()
    obs_scan.data = obs_scan.tlist(scan_gather=True)[int(scan_number)]

    # Add noise
    # Fractional systematic noise
    obs_scan = obs_scan.add_fractional_noise(noise_frac)

    # Refractive noise floor
    # Should we estimate zbl from the data?
    obs_scan.data['sigma'] = (obs_scan.data['sigma']**2 + (refr_floor_frac*zbl)**2)**0.5

    # Deblur
    sm = so.ScatteringModel()
    obs_scan = sm.Deblur_obs(obs_scan.switch_polrep("stokes"))

    return obs_scan

averaging_modes = ['scan-average', 'scan-subaverage', 'snapshot']
noise_fracs = [0.02, 0.05]

for averaging_mode in averaging_modes:
    for noise_frac in noise_fracs:
        # Set up paths
        obs_filename = uvfits_path.split('/')[-1].split('.uvf')[0]
        obs_dir = '/'.join(uvfits_path.split('/')[:-1])
        out_dir = obs_dir + '/preprocessed_test/%s/noisefrac%s/'%(averaging_mode, noise_frac)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        # Load data
        obs = eh.obsdata.load_uvfits(uvfits_path)

        # Divide into scans
        obs_scans = generate_scans(obs, averaging_mode)

        # Add noise, collect scan info, and export scans
        nscans = len(obs_scans.tlist(scan_gather=True))
        tstart = np.zeros(nscans)
        tstop = np.zeros(nscans)
        sites = []
        for scan_number in range(nscans):
            obs_scan_processed = add_noise_deblur(obs_scans, scan_number, noise_frac=noise_frac)

            tstart[scan_number] = obs_scan_processed.data['time'][0] 
            tstop[scan_number] = obs_scan_processed.data['time'][-1] 

            scan_sites = []
            for site in obs_scan_processed.tarr['site']:
                if site in obs_scan_processed.data['t1'] or site in obs_scan_processed.data['t2']:
                    scan_sites.append(site)
            sites.append(scan_sites)

            out_file = out_dir + obs_filename + '_preprocessed_%s_noisefrac%s_scan%s.uvfits'%(averaging_mode, noise_frac, scan_number)
            obs_scan_processed.save_uvfits(out_file)

        # Write out scan info
        if noise_frac == noise_fracs[0]:
            with open(obs_dir + '/preprocessed_test/%s/scaninfo_%s.txt'%(averaging_mode, averaging_mode), 'w') as f:
                print('# Scan number, tstart, tstop, sites', end='\n', file=f)
                for i in range(nscans):
                    print('%03d %.5f %.5f %s'%(i, tstart[i], tstop[i], ','.join(sites[i])), end='\n', file=f)
            f.close()



