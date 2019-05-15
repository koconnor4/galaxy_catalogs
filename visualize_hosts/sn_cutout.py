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

def ellipse(file,possible_hosts,sn_position,title='', save = False,show=False):
    
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
    #print(np.amin(hdu.data ))
    norms = [matplotlib.colors.LogNorm(vmin=np.abs(np.amin(hdu.data)),vmax=np.amax(hdu.data)+2*np.abs(np.amin(hdu.data))),matplotlib.colors.PowerNorm(gamma=2.)] 
    # diverging color maps...'handwave switches color'
    cmaps_d = ['PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
            'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']
    # sequential color maps...'handwave doesn't switch color'
    cmaps_s = ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
            'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
            'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']
    

    ax.imshow(hdu.data,norm=norms[0],cmap=cmaps_s[0])

    ax.set_xlabel(r'$\alpha$')
    ax.set_ylabel(r'$\delta$')
    ax.set_title(title)
    
    sn_pix_coords = wcs.utils.skycoord_to_pixel(sn_position, w,origin=1,mode='all')
    ax.scatter([int(sn_pix_coords[0])],[int(sn_pix_coords[1])],marker='X',c='black',edgecolor='white',label='SN')
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
        galid = np.float(possible_hosts[i][4])

        colors = ['red','blue','green','cyan','violet','turquoise','purple','yellow','orange','lime','navy','pink','saddlebrown','silver','gold','navajowhite']
        ax.add_patch(Ellipse((x,y),width = 5*gal_a, height = 5*gal_b, angle = gal_theta,color=colors[i],label=galid,facecolor=None,fill=False))

    plt.legend(loc='upper right', bbox_to_anchor=(1.4, 1.0))
    if save == True:
        plt.savefig('nearby_host_ellipses/'+file[:8])
    if show == True:
        plt.show()
    return 







