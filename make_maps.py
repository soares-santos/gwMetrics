import numpy as np
import healpy as hp
import matplotlib.pyplot as plt


def pix2radec(ipix,nside=512):
    theta, phi = hp.pix2ang(nside,ipix,nest=False)
    dec = np.rad2deg(0.5 * np.pi - theta) 
    ra = np.rad2deg(phi)
    return ra,dec


def gaussian2d_from_sample_map(ra_coord=0.,dec_coord=0.,sigma_ra=5.,sigma_dec=5.,nside=128,sample_map_file_name="datafiles/bayestar.fits.gz"):

    m = hp.read_map(sample_map_file_name,verbose=False)

    if isinstance(ra_coord, (list, tuple, np.ndarray)) and len(ra_coord) == len(dec_coord):

        maps = []

        if isinstance(sigma_ra, (list, tuple, np.ndarray)) and len(sigma_ra) == len(sigma_dec):
            for i in range(len(ra_coord)):
                maps.append(one_gaussian2d_from_sample_map(m,ra_coord[i],dec_coord[i],sigma_ra[i],sigma_dec[i],nside))

        if isinstance(sigma_ra, (float,int)) and isinstance(sigma_dec, (float,int)):
            for i in range(len(ra_coord)):
                maps.append(one_gaussian2d_from_sample_map(m,ra_coord[i],dec_coord[i],sigma_ra,sigma_dec,nside)) 

        return maps

    if isinstance(ra_coord, (float,int)) and isinstance(dec_coord, (float,int)) and isinstance(sigma_ra, (float,int)) and isinstance(sigma_dec, (float,int)):
        return one_gaussian2d_from_sample_map(m,ra_coord,dec_coord,sigma_ra,sigma_dec,nside_p)


def one_gaussian2d_from_sample_map(hpx,ra,dec,sigma_ra,sigma_dec,nside_p):

    #hpx = hp.read_map("bayestar.fits.gz",verbose=False)

    #ra = 0.
    #dec = 0.

    #sigma_ra = 5.
    #sigma_dec = 5.

    #nside_p = 128

    s = np.random.normal(ra,sigma_ra,len(hpx))
    count, bins, ignored = plt.hist(s, nside_p*10, normed=False)
    ra_prob = count/np.sum(count)
    ra_bin_edges = bins[:-1]

    s = np.random.normal(dec,sigma_dec,len(hpx))
    count, bins, ignored = plt.hist(s, nside_p*10, normed=False)
    dec_prob = count/np.sum(count)
    dec_bin_edges = bins[:-1]

    hpx_p=hpx

    nside = hp.npix2nside(len(hpx))

    if nside > nside_p : hpx_p=hp.ud_grade(hpx_p, nside_p)

    for pix in range(len(hpx_p)):
        ra_pix,dec_pix=pix2radec(pix,nside_p)
        if ra_pix > 180. : ra_pix = 360. - ra_pix 
        ra_idx=(np.abs(ra_bin_edges-ra_pix)).argmin()
        dec_idx=(np.abs(dec_bin_edges-dec_pix)).argmin()
        ra_p=ra_prob[ra_idx]
        dec_p=dec_prob[dec_idx]
        hpx_p[pix]=hpx_p[pix] * 0.0 + (ra_p * dec_p) 

    if nside > nside_p : hpx_p=hp.ud_grade(hpx_p, nside, power=-2)

    #hp.mollview(hpx_p)

    return hpx_p

