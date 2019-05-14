import copy
import pickle
import numpy as np
import glob

# matplotlib
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

# astropy
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.utils.data import download_file
from astropy.coordinates import SkyCoord
from astropy import wcs
import astropy.units as u



def download_image_save_cutout(position, size, cutout_filename = 'example_cutout.fits'):
    
    # Download the image
    # commented out since I already have file, could be useful just change filename argument to url
    # filename = download_file(url)
    field = cutout_filename[:2].lower()
    if field == 'co':
        f125w = r"C:\Users\Kyle\Pictures\COSMOS\cos_2epoch_wfc3_f125w_060mas_v1.0_drz.fits"
        f160w = r"C:\Users\Kyle\Pictures\COSMOS\cos_2epoch_wfc3_f160w_060mas_v1.0_drz.fits"
    elif field == 'eg':
        f125w = r"C:\Users\Kyle\Pictures\AEGIS...EGS\aegis_3dhst.v4.0.F125W_orig_sci.fits.gz"
        f160w = r"C:\Users\Kyle\Pictures\AEGIS...EGS\aegis_3dhst.v4.0.F160W_orig_sci.fits.gz"
    elif field == 'ud':
        f125w = r"C:\Users\Kyle\Pictures\UDS\uds_3dhst.v4.0.F125W_orig_sci.fits.gz"
        f160w = r"C:\Users\Kyle\Pictures\UDS\uds_3dhst.v4.0.F160W_orig_sci.fits.gz"
    elif field == 'gs':
        f125w = r"C:\Users\Kyle\Pictures\GOODS-S\goodss_3dhst.v4.0.F125W_orig_sci.fits.gz"
        f160w = r"C:\Users\Kyle\Pictures\GOODS-S\goodss_3dhst.v4.0.F160W_orig_sci.fits.gz"
    elif field == 'gn':
        f125w = r"C:\Users\Kyle\Pictures\GOODS-N\goodsn_3dhst.v4.0.F125W_orig_sci.fits.gz"
        f160w = r"C:\Users\Kyle\Pictures\GOODS-N\goodsn_3dhst.v4.0.F160W_orig_sci.fits.gz"
    filename = f160w

    # Load the image and the WCS
    hdu = fits.open(filename)[0]
    w = wcs.WCS(hdu.header)
    
    # the image is already drizzled (distortion corrected) 
    # we shouldn't be using sip distortion coefficients even though they are still in the header
    # w.sip = None temporarily removes these sip coefficients from a wcs obj
    # note maybe future should have an if drizzled = True: to decide on applying this 
    w.sip = None
    # still get warning about sip coefficients but think this is okay
    # the wcs object is created from hdu.header which has the sip coefficients 
    # might be  better to remove the sip coefficients from the hdu header which suspect would get rid of warning 
    # however the sip coefficients should now be removed from the wcs object  

    # Make the cutout, including the WCS
    cutout = Cutout2D(hdu.data, position=position, size=size, wcs=w)

    # Put the cutout image in the FITS HDU 
    hdu.data = cutout.data

    # Update FITS header with the cutout WCS
    hdu.header.update(cutout.wcs.to_header())

    # Write the cutout to a new FITS file
    #cutout_filename = 'example_cutout.fits'
    hdu.writeto('cutouts/'+cutout_filename, overwrite=True)

def ellipse(file,possible_hosts,sn_position,title='',save = False,show=False):
    
    # access the file
    fits_image = fits.open('cutouts/'+file)
    hdu = fits_image[0]
    
    # the deg/pixel for the fits image along ra and dec
    # the ellipse (add_patch) needs to know pixel dimensions of a,b to plot on the image 
    ra_conv = hdu.header.get('CD1_1')
    dec_conv = hdu.header.get('CD1_1')
    
    # the wcs obj from header of the image;
    # the image is drizzled ~ distortion corrected, and should be  w.sip = None already
    # just doing it again for consistency/peace of mind
    w = wcs.WCS(hdu.header)
    w.sip = None
    
    # initialize the figure
    f= plt.figure()
    ax=f.gca(projection=w)
    #ax.imshow(hdu.data,vmin=0,vmax=0.015)
    # try with norm
    #ax.imshow(hdu.data,norm=matplotlib.colors.Normalize(vmin=-.015,vmax=.015))
    print(np.amin(hdu.data ))
    norms = [matplotlib.colors.LogNorm(vmin=np.abs(np.amin(hdu.data)),vmax=np.amax(hdu.data)+2*np.abs(np.amin(hdu.data))),matplotlib.colors.PowerNorm(gamma=2.)] 
    # diverging color maps...'handwave switches color'
    cmaps_d = ['PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
            'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']
    # sequential color maps...'handwave doesn't switch color'
    cmaps_s = ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
            'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
            'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']
    ax.imshow(hdu.data,norm=norms[0],cmap=cmaps_d[5])

    ax.set_xlabel(r'$\alpha$')
    ax.set_ylabel(r'$\delta$')
    ax.set_title(title)
    
    sn_pix_coords = wcs.utils.skycoord_to_pixel(sn_position, w,origin=1,mode='all')
    ax.scatter([int(sn_pix_coords[0])],[int(sn_pix_coords[1])],marker='X',c='white',edgecolor='white')
    pix_coords = []
    
    # now into the possible hosts to get major,minor,position angle and plot ellipse at loc 
    for i in range(len(possible_hosts)):
        gal_position = possible_hosts[i][0] # SkyCoord ra,dec for the host  
        # get the pixel coords for the center of the host using wcs
        gal_pixel_coords = wcs.utils.skycoord_to_pixel(gal_position, w, origin=1, mode='all')
        x = int(gal_pixel_coords[0])
        y = int(gal_pixel_coords[1])
        # the major,minor are given in pixels...
        # gal_a, gal_b already in right pixels for putting on the image since I just cropped same image
        gal_a = np.float(possible_hosts[i][1]) # pixel
        gal_b = np.float(possible_hosts[i][2]) # pixel
        gal_theta = np.float(possible_hosts[i][3]) # deg ccw/x
        ax.add_patch(Ellipse((x,y),width = 5*gal_a, height = 5*gal_b, angle = gal_theta,color='black',facecolor=None,fill=False))
    if save == True:
        plt.savefig('nearby_host_ellipses/'+file[:8])
    if show == True:
        plt.show()
    return 










with open (r'C:/Users/Kyle/Documents/school/research/summer18/Catalog/supernovae/sourcefiles/candels_sn', 'rb') as fp:
    candels = pickle.load(fp)
"""
cutout_filenames = [candels[i]['name'] for i in range(len(candels))]
coords = [SkyCoord(candels[i]['RA'],candels[i]['DEC'], unit=(u.hourangle,u.deg)) for i in range(len(candels))]
print(cutout_filenames)
print(coords)
for i in range(len(candels)):
    sized = u.Quantity((20, 20), u.arcsec)
    download_image_save_cutout(coords[i],sized,cutout_filename=cutout_filenames[i]+'.fits')
"""

with open (r'C:/Users/Kyle/Documents/school/research/summer18/Catalog/candels_possible_hosts', 'rb') as fp:
    candels_possible_hosts = pickle.load(fp)
with open('candels_possible_hosts.pkl','rb') as fp:
    candels_possible_hosts_full = pickle.load(fp)










nicks = [candels_possible_hosts[i]['snid'][0] for i in range(len(candels_possible_hosts))]
for nick in nicks:
    #print('Looking at SN:',nick.lower())
    ellip_params = []

    shapes = []
    #print('Number of hosts within 5 arcsec are:',len(candels_possible_hosts_full[nicks.index(nick)]))
    for i in candels_possible_hosts_full[nicks.index(nick)]:
        if nick[:2].lower() == 'co':
            shapes.append(i[-1][-36:-27])
        elif nick[:2].lower() == 'eg':
            shapes.append(i[-1][457:466])
        elif nick[:2].lower() == 'gs':
            shapes.append(i[-1][393:401])
        elif nick[:2].lower() == 'gn':
            shapes.append(i[-1][83:86])
        elif nick[:2].lower() == 'ud':
            shapes.append(i[-1][91:99])
    # cxx,cyy,cxy,a,b,da,db,theta,dtheta...a,b,theta what we need
    # cos ''_full[-1][-36:-27] ...3,4,7
    # egs ''_full[-1][457:466]... 3,4,7
    # goods-s ''_full[-1][393:401]...3,4,7
    # uds ''_full[-1][91:99]...3,5,7
    # goods-n ''_full[-1][83:86]...1,2,3 

    #print('They are at coords and have a,b,theta:')
    #for i in range(len(candels_possible_hosts[nicks.index(nick)])):
    for i in range(3):
        if nick[:2].lower() == 'co' or nick[:2].lower() == 'eg' or nick[:2].lower() == 'gs':
            a = 3
            b = 4
            theta = 7
        elif nick[:2].lower() == 'gn':
            a = 0
            b = 1
            theta  = 2
        elif nick[:2].lower() == 'ud':
            a = 3
            b = 5
            theta = 7
        ellip_params.append([candels_possible_hosts_full[nicks.index(nick)][i][0][0],
                     shapes[i][a],shapes[i][b],shapes[i][theta]])
    #print(ellip_params)
    #print(candels_possible_hosts[nicks.index(nick)]['snid'][0])
    nicksy = [candels[i]['name'] for i in range(len(candels))]
    #print(candels[0])
    sn_coords = SkyCoord(candels[nicksy.index(nick)]['RA'],candels[nicksy.index(nick)]['DEC'], unit=(u.hourangle,u.deg))
    #print('SN is at:',sn_coords)
    #print(candels[nicksy.index(nick)]['name'])

    ellipse(nick.lower()+'.fits',ellip_params,sn_position=sn_coords,title='temp_scale',save=False,show=True)