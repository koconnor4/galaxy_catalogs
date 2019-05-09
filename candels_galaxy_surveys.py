import numpy as np
import csv
from astropy.io import ascii
from astropy import cosmology as cosmo
from astropy.coordinates import SkyCoord,Angle  
from astropy.table import Table,vstack
import astropy.units as u
import pickle
import os
import matplotlib.pyplot as plt 
from operator import itemgetter

__current_dir__=os.path.abspath(os.getcwd())
__data_dir__=os.path.join(__current_dir__,"..","galaxy_data")

__all__=['galaxy_survey','galaxy_catalog']
class galaxy_survey(object):
	def __init__(self,source,physpar,mass,photom):
		self.source = source
		self.physpar = physpar
		self.mass = mass
		self.photom = photom
		print(self.source,source)
		if source == 'Santini2015; goods-s_candels; https://archive.stsci.edu/prepds/candels/':
			self.redshift_names='zbest'
			self.ra_names='RA (F160W coordinate, J2000, degree)'
			self.ra_units=u.deg 
			self.dec_names='DEC (F160W coordinate, J2000, degree)'
			self.dec_units=u.deg 
			self.kron_radius_names = 'KRON_RADIUS Kron apertures in units of A or B'
			self.cxx='CXX_IMAGE Cxx object ellipse parameter [pixel**(-2)]' 
			self.cyy='CYY_IMAGE Cyy object ellipse parameter [pixel**(-2)]' 
			self.cxy='CXY_IMAGE Cxy object ellipse parameter [pixel**(-2)]'
			self.a_ellip='A_IMAGE Profile RMS along major axis [pixel]' 
			self.b_ellip='B_IMAGE Profile RMS along major axis [pixel]'
			self.da_ellip='ERRA_IMAGE RMS position error along major axis [pixel]'
			self.db_ellip='ERRB_IMAGE RMS position error along major axis [pixel]' 
			self.theta_ellip='THETA_IMAGE Position angle (CCW/x) [deg]' 
			self.dtheta_ellip='ERRTHETA_IMAGE Error ellipse position angle (CCW/x) [deg]' 
			# DAFUQ is tau, method 2a tau, other methods? 
			self.age_2a_tau='age_2a_tau'
			self.sfr_2a_tau='SFR_2a_tau'
			self.tau_2a_tau='tau_2a_tau'
		elif source == 'Stefanon2017; egs_candels; https://archive.stsci.edu/prepds/candels/':
			self.redshift_names='zbest'
			self.ra_names='  3 RA                            (deg) '
			self.ra_units=u.deg 
			self.dec_names='  4 DEC                           (deg) '
			self.dec_units=u.deg 
			self.cxx='458 CXX_IMAGE                     (pixel**(-2)) '
			self.cyy='459 CYY_IMAGE                     (pixel**(-2)) ' 
			self.cxy='460 CXY_IMAGE                     (pixel**(-2)) '
			self.a_ellip='461 A_IMAGE                       (pixel) ' 
			self.b_ellip='462 B_IMAGE                       (pixel) '
			self.da_ellip='463 ERRA_IMAGE                    (pixel) '
			self.db_ellip='464 ERRB_IMAGE                    (pixel) ' 
			self.theta_ellip='465 THETA_IMAGE                   (deg) ' 
			self.dtheta_ellip='466 ERRTHETA_IMAGE                (deg) ' 
			self.kron_radius_names = '522 KRON_RADIUS                   (pixel) '
			self.age_2a_tau='  2 age_2a_tau        (dex(yr)) '
			self.sfr_2a_tau='  5 SFR_2a_tau        (solMass/yr) '
			self.tau_2a_tau='  3 tau_2a_tau        (Gyr) '
		elif source == 'Nayyeri2017; cos_candels; https://archive.stsci.edu/prepds/candels/':
			self.redshift_names='zbest'
			self.ra_names='RA (deg) '
			self.ra_units=u.deg 
			self.dec_names='DEC (deg) '
			self.dec_units=u.deg 
			self.kron_radius_names = 'KRON_RADIUS (pixel) '
			self.cxx='CXX_IMAGE (pixel**(-2)) ' 
			self.cyy='CYY_IMAGE (pixel**(-2)) ' 
			self.cxy='CXY_IMAGE (pixel**(-2)) '
			self.a_ellip='A_IMAGE (pixel) ' 
			self.b_ellip='B_IMAGE (pixel) '
			self.da_ellip='ERRA_IMAGE (pixel) '
			self.db_ellip='ERRB_IMAGE (pixel) ' 
			self.theta_ellip='THETA_IMAGE (deg) ' 
			self.dtheta_ellip='ERRTHETA_IMAGE (deg) '
			self.age_2a_tau='age_2a_tau (dex(yr)) '
			self.sfr_2a_tau='SFR_2a_tau (solMass/yr) '
			self.tau_2a_tau='tau_2a_tau (Gyr) ' 
		elif source == 'Santini2015; uds_candels; https://archive.stsci.edu/prepds/candels/':
			self.redshift_names='zbest'
			self.ra_names='RA'
			self.ra_units=u.deg 
			self.dec_names='Dec'
			self.dec_units=u.deg 
			self.cxx='cxx_image'
			self.cyy='cyy_image' 
			self.cxy='cxy_image'
			self.a_ellip='a_image' 
			self.b_ellip='b_image'
			self.da_ellip='erra_image'
			self.db_ellip='errb_image' 
			self.theta_ellip='theta_image' 
			self.dtheta_ellip='errtheta_image' 
			self.kron_radius_names = 'kron_radius'
			self.age_2a_tau='age_2a_tau'
			self.sfr_2a_tau='SFR_2a_tau'
			self.tau_2a_tau='tau_2a_tau'
		elif source == 'Skelton14; goods-n_3dhst; https://3dhst.research.yale.edu/Data.php':
			self.redshift_names='z'
			self.ra_names='ra'
			self.ra_units=u.deg 
			self.dec_names='dec'
			self.dec_units=u.deg
			self.kron_radius_names = 'kron_radius'
			self.a_ellip='a_image'
			self.b_ellip='b_image'
			# the following are columns in the mass table; others have similar cols in physpar table don't know method but isnt necessarily same as 2a tau as all others. 
			self.lage='lage'
			self.ltau='ltau'
			self.lsfr='lsfr'
			self.lmass='lmass'
		else: print('whoa now, slow down, you didnt go into any of the sources')
		
	def __str__(self):
		print('howdy partner I am a galaxy survey')
		return('ok good talk, take care now')

	def gal_skycoords(self):
		locs = []
		for i in range(len(self.mass)):
			loc = SkyCoord(self.photom[i][self.ra_names],self.photom[i][self.dec_names],unit=(self.ra_units,self.dec_units))
			locs.append(loc)		
		
		return locs


	def nearest(self,loc,conversion = 3600.0/0.06):
		# x,y in deg so 3600 '' / deg then f160w image has .06 ''/pixel, we need cancellation with cij in pixel^-2; returns R ~ units of A or B semi-axes the elliptical radii away
		locs = self.gal_skycoords()
		Rs = []
		for i in range(len(self.photom)):
			dra, ddec = loc.spherical_offsets_to(locs[i])
			x = dra.value
			y = ddec.value
			x*=conversion
			y*=conversion
			R = (np.float(self.photom[i][self.cxx])*x**2 + np.float(self.photom[i][self.cyy])*y**2 + np.float(self.photom[i][self.cxy])*x*y)**(1/2)
			Rs.append(R)
		nearest = min(enumerate(Rs), key=itemgetter(1))[0] # gets the index of the minimum R
		return [galaxy_survey(self.source,self.physpar[nearest],self.mass[nearest],self.photom[nearest]),Rs[nearest],Rs]




class galaxy_catalog(object):
	def __init__(self,surveys):
		if not isinstance(surveys,(list,tuple)):
			surveys = [surveys]
		self.survey_dict=dict([])
		self.redshift_names=dict([])
		self.ra_names=dict([])
		self.ra_units=dict([])
		self.dec_names=dict([])
		self.dec_units=dict([])
		self.kron_radius_names=dict([])
		self.cxx=dict([])
		self.cyy=dict([])
		self.cxy=dict([])
		self.cyy=dict([])
		self.a_ellip=dict([])
		self.b_ellip=dict([])
		self.da_ellip=dict([])
		self.db_ellip=dict([])
		self.theta_ellip=dict([])
		self.dtheta_ellip=dict([])

		for survey in surveys:
			if survey == 'goodss':
				with open(os.path.join(__data_dir__,survey+'.pkl'),'rb') as fp:
					self.goodss = pickle.load(fp)
					self.survey_dict[survey]=self.goodss
					self.redshift_names[survey]='zbest'
					self.ra_names[survey]='RA (F160W coordinate, J2000, degree)'
					self.ra_units[survey]=u.deg 
					self.dec_names[survey]='DEC (F160W coordinate, J2000, degree)'
					self.dec_units[survey]=u.deg 
					self.kron_radius_names[survey] = 'KRON_RADIUS Kron apertures in units of A or B'
					self.cxx[survey]='CXX_IMAGE Cxx object ellipse parameter [pixel**(-2)]' 
					self.cyy[survey]='CYY_IMAGE Cyy object ellipse parameter [pixel**(-2)]' 
					self.cxy[survey]='CXY_IMAGE Cxy object ellipse parameter [pixel**(-2)]'
					self.a_ellip[survey]='A_IMAGE Profile RMS along major axis [pixel]' 
					self.b_ellip[survey]='B_IMAGE Profile RMS along major axis [pixel]'
					self.da_ellip[survey]='ERRA_IMAGE RMS position error along major axis [pixel]'
					self.db_ellip[survey]='ERRB_IMAGE RMS position error along major axis [pixel]' 
					self.theta_ellip[survey]='THETA_IMAGE Position angle (CCW/x) [deg]' 
					self.dtheta_ellip[survey]='ERRTHETA_IMAGE Error ellipse position angle (CCW/x) [deg]' 
			elif survey == 'egs':
				with open(os.path.join(__data_dir__,survey+'.pkl'),'rb') as fp:
					self.egs = pickle.load(fp)
					self.survey_dict[survey]=self.egs
					self.redshift_names[survey]='zbest'
					self.ra_names[survey]='  3 RA                            (deg) '
					self.ra_units[survey]=u.deg 
					self.dec_names[survey]='  4 DEC                           (deg) '
					self.dec_units[survey]=u.deg 
					self.cxx[survey]='458 CXX_IMAGE                     (pixel**(-2)) '
					self.cyy[survey]='459 CYY_IMAGE                     (pixel**(-2)) ' 
					self.cxy[survey]='460 CXY_IMAGE                     (pixel**(-2)) '
					self.a_ellip[survey]='461 A_IMAGE                       (pixel) ' 
					self.b_ellip[survey]='462 B_IMAGE                       (pixel) '
					self.da_ellip[survey]='463 ERRA_IMAGE                    (pixel) '
					self.db_ellip[survey]='464 ERRB_IMAGE                    (pixel) ' 
					self.theta_ellip[survey]='465 THETA_IMAGE                   (deg) ' 
					self.dtheta_ellip[survey]='466 ERRTHETA_IMAGE                (deg) ' 
					self.kron_radius_names[survey] = '522 KRON_RADIUS                   (pixel) '
			elif survey == 'cos':
				with open(os.path.join(__data_dir__,survey+'.pkl'),'rb') as fp:
					self.cos = pickle.load(fp)
					self.survey_dict[survey]=self.cos
					self.redshift_names[survey]='zbest'
					self.ra_names[survey]='RA (deg) '
					self.ra_units[survey]=u.deg 
					self.dec_names[survey]='DEC (deg) '
					self.dec_units[survey]=u.deg 
					self.kron_radius_names[survey] = 'KRON_RADIUS (pixel) '
					self.cxx[survey]='CXX_IMAGE (pixel**(-2)) ' 
					self.cyy[survey]='CYY_IMAGE (pixel**(-2)) ' 
					self.cxy[survey]='CXY_IMAGE (pixel**(-2)) '
					self.a_ellip[survey]='A_IMAGE (pixel) ' 
					self.b_ellip[survey]='B_IMAGE (pixel) '
					self.da_ellip[survey]='ERRA_IMAGE (pixel) '
					self.db_ellip[survey]='ERRB_IMAGE (pixel) ' 
					self.theta_ellip[survey]='THETA_IMAGE (deg) ' 
					self.dtheta_ellip[survey]='ERRTHETA_IMAGE (deg) ' 
			elif survey == 'uds':
				with open(os.path.join(__data_dir__,survey+'.pkl'),'rb') as fp:
					self.uds = pickle.load(fp)
					self.survey_dict[survey]=self.uds
					self.redshift_names[survey]='zbest'
					self.ra_names[survey]='RA'
					self.ra_units[survey]=u.deg 
					self.dec_names[survey]='Dec'
					self.dec_units[survey]=u.deg 
					self.cxx[survey]='cxx_image'
					self.cyy[survey]='cyy_image' 
					self.cxy[survey]='cxy_image'
					self.a_ellip[survey]='a_image' 
					self.b_ellip[survey]='b_image'
					self.da_ellip[survey]='erra_image'
					self.db_ellip[survey]='errb_image' 
					self.theta_ellip[survey]='theta_image' 
					self.dtheta_ellip[survey]='errtheta_image' 
					self.kron_radius_names[survey] = 'kron_radius'
			elif survey == 'goodsn':
				with open(os.path.join(__data_dir__,survey+'.pkl'),'rb') as fp:
					self.goodsn = pickle.load(fp)
					self.survey_dict[survey]=self.goodsn
					self.redshift_names[survey]='z'
					self.ra_names[survey]='ra'
					self.ra_units[survey]=u.deg 
					self.dec_names[survey]='dec'
					self.dec_units[survey]=u.deg
					self.kron_radius_names[survey] = 'kron_radius'
					self.a_ellip[survey]='a_image'
					self.b_ellip[survey]='b_image'
			else: print('I didnt go into any of these!')

	def __str__(self):
		print('hello, I am a galaxy catalog')
		return('bye')
	
	def redshift_distribution(self,surveys,minz=0,maxz=3,plot=False):
		if not isinstance(surveys,(list,tuple)):
			surveys = [surveys]
		for survey in surveys:
			#print(len(self.survey_dict[survey].mass))
			#for i in range(len(self.survey_dict[survey].mass)):
			masses_in_z = []
			photoms_in_z = []
			physpars_in_z = []
			for i in range(len(self.survey_dict[survey].mass)):
				if np.float(self.survey_dict[survey].mass[i][self.redshift_names[survey]]) > minz and np.float(self.survey_dict[survey].mass[i][self.redshift_names[survey]]) < maxz:
					masses_in_z.append(self.survey_dict[survey].mass[i])
					photoms_in_z.append(self.survey_dict[survey].photom[i])
					physpars_in_z.append(self.survey_dict[survey].physpar[i])
			a = vstack([i for i in physpars_in_z])
			b = vstack([i for i in masses_in_z])
			c = vstack([i for i in photoms_in_z])
			source = self.survey_dict[survey].source
			in_z_survey = galaxy_survey(source,a,b,c)
			if plot == True:
				plt.hist([in_z_survey.mass[self.redshift_names[survey]]])
				plt.xlabel('Redshift')
				plt.show()

			return in_z_survey
	
	def gal_skycoords(self,surveys):
		if not isinstance(surveys,(list,tuple)):
			surveys = [surveys]
		locs = []
		for survey in surveys:
			for i in range(len(self.survey_dict[survey].mass)):
				loc = SkyCoord(self.survey_dict[survey].photom[i][self.ra_names[survey]],self.survey_dict[survey].photom[i][self.dec_names[survey]],unit=(self.ra_units[survey],self.dec_units[survey]))
				locs.append(loc)		
		
		return locs
	
	def nearby(self,surveys,loc,within=Angle(5/3600,u.deg)):
		if not isinstance(surveys,(list,tuple)):
			surveys = [surveys]
		near_masses = []
		near_photoms = []
		near_physpars = []
		for survey in surveys:
			locs = [i for i in self.gal_skycoords(survey)]
			for i in range(len(locs)):
				if loc.separation(locs[i]) < within:
					near_masses.append(self.survey_dict[survey].mass[i])
					near_photoms.append(self.survey_dict[survey].photom[i])
					near_physpars.append(self.survey_dict[survey].physpar[i])
			a = vstack([i for i in near_physpars])
			b = vstack([i for i in near_masses])
			c = vstack([i for i in near_photoms])
			source = self.survey_dict[survey].source
			nearby_gals = galaxy_survey(source,a,b,c)
		
		return nearby_gals

	# this not useful here... faster to get nearby using this catalog object; it will output survey of those nearby. Then only do calc to find nearest using nearby survey. 
	def nearest(self,surveys,loc,conversion = 1):
		if not isinstance(surveys,(list,tuple)):
			surveys = [surveys]
		Rs = []
		for survey in surveys:
			locs = [i for i in self.gal_skycoords(survey)]
			for i in range(len(locs)):
				dra, ddec = loc.spherical_offsets_to(locs[i])
				x = np.abs(dra.value)
				y = np.abs(ddec.value)
				R = (np.float(self.survey_dict[survey].photom[i][self.cxx[survey]])*x**2 + np.float(self.survey_dict[survey].photom[i][self.cyy[survey]])*y**2 + np.float(self.survey_dict[survey].photom[i][self.cxy[survey]])*x*y)*conversion
				Rs.append(R)
		nearest = min(enumerate(Rs), key=itemgetter(1))[0] # gets the index of the minimum R
		return galaxy_survey(self.survey_dict[survey].source,self.survey_dict[survey].physpar[nearest],self.survey_dict[survey].mass[nearest],self.survey_dict[survey].photom[nearest])


def output_survey(survey):
	with open(r'C:/Users/Kyle/Documents/school/research/summer18/Catalog/galaxies/'+survey+'0to3','rb') as fp:
		surveyData = pickle.load(fp)
	
	if survey == 'goodss':
		source = 'Santini2015; goods-s_candels; https://archive.stsci.edu/prepds/candels/'
		goodssphysparheader = np.loadtxt('goodssphysparheader.txt', dtype='str',delimiter='\n')
		goodsstotmultiheader = np.loadtxt('goodsstotmultiheader.txt',dtype='str',delimiter='\n')
		goodssmasscatheader = np.loadtxt('goodssmasscatheader.txt',dtype='str',delimiter='\n')
		phys_row = np.array([surveyData[i][-2] for i in range(len(surveyData))])
		photom_row = np.array([surveyData[i][-1] for i in range(len(surveyData))])
		mass_row = np.array([surveyData[i][-3] for i in range(len(surveyData))])
		phys_test = Table(phys_row,names=goodssphysparheader)
		photom_test = Table(photom_row,names=goodsstotmultiheader)
		mass_test = Table(mass_row,names=goodssmasscatheader)
		my_survey = galaxy_survey(source,phys_test,mass_test,photom_test)
	elif survey == 'egs':
		source = 'Stefanon2017; egs_candels; https://archive.stsci.edu/prepds/candels/'
		egsphysparheader = np.loadtxt('egsphysparheader.txt', dtype='str',delimiter='\n')
		egstotmultiheader = np.loadtxt('egstotmultiheader.txt',dtype='str',delimiter='\n')
		egsmasscatheader = np.loadtxt('egsmasscatheader.txt',dtype='str',delimiter='\n')
		phys_row = np.array([surveyData[i][-2] for i in range(len(surveyData))])
		photom_row = np.array([surveyData[i][-1] for i in range(len(surveyData))])
		mass_row = np.array([surveyData[i][-3] for i in range(len(surveyData))])
		phys_test = Table(phys_row,names=egsphysparheader)
		photom_test = Table(photom_row,names=egstotmultiheader)
		mass_test = Table(mass_row,names=egsmasscatheader)
		my_survey = galaxy_survey(source,phys_test,mass_test,photom_test)
	elif survey == 'cos':
		source = 'Nayyeri2017; cos_candels; https://archive.stsci.edu/prepds/candels/'
		cosphysparheader = np.loadtxt('cosphysparheader.txt', dtype='str',delimiter='\n')
		costotmultiheader = np.loadtxt('costotmultiheader.txt',dtype='str',delimiter='\n')
		cosmasscatheader = np.loadtxt('cosmasscatheader.txt',dtype='str',delimiter='\n')
		phys_row = np.array([surveyData[i][-2] for i in range(len(surveyData))])
		photom_row = np.array([surveyData[i][-1] for i in range(len(surveyData))])
		mass_row = np.array([surveyData[i][-3] for i in range(len(surveyData))])
		phys_test = Table(phys_row,names=cosphysparheader)
		photom_test = Table(photom_row,names=costotmultiheader)
		mass_test = Table(mass_row,names=cosmasscatheader)
		my_survey = galaxy_survey(source,phys_test,mass_test,photom_test)
	elif survey == 'uds':
		source = 'Santini2015; uds_candels; https://archive.stsci.edu/prepds/candels/'
		udsphysparheader = np.loadtxt('udsphysparheader.txt', dtype='str',delimiter='\n')
		udstotmultiheader = np.loadtxt('udstotmultiheader.txt',dtype='str',delimiter='\n')
		udsmasscatheader = np.loadtxt('udsmasscatheader.txt',dtype='str',delimiter='\n')
		phys_row = np.array([surveyData[i][-2] for i in range(len(surveyData))])
		photom_row = np.array([surveyData[i][-1] for i in range(len(surveyData))])
		mass_row = np.array([surveyData[i][-3] for i in range(len(surveyData))])
		print(mass_row[0])
		print(len(mass_row),mass_row.shape)
		print(len(phys_row),phys_row.shape)
		print(len(photom_row),photom_row.shape)
		phys_test = Table(phys_row,names=udsphysparheader)
		photom_test = Table(photom_row,names=udstotmultiheader)
		mass_test = Table(mass_row,names=udsmasscatheader)
		my_survey = galaxy_survey(source,phys_test,mass_test,photom_test)
	elif survey == 'goodsn':
		source = 'Skelton14; goods-n_3dhst; https://3dhst.research.yale.edu/Data.php'
		goodsnphysparheader = np.loadtxt('goodsnphysparheader.txt', dtype='str',delimiter='\n')
		goodsntotmultiheader = np.loadtxt('goodsntotmultiheader.txt',dtype='str',delimiter='\n')
		goodsnmasscatheader = np.loadtxt('goodsnmasscatheader.txt',dtype='str',delimiter='\n')
		phys_row = np.array([surveyData[i][-2] for i in range(len(surveyData))])
		photom_row = np.array([surveyData[i][-1] for i in range(len(surveyData))])
		mass_row = np.array([surveyData[i][-3] for i in range(len(surveyData))])
		phys_test = Table(phys_row,names=goodsnphysparheader)
		photom_test = Table(photom_row,names=goodsntotmultiheader)
		mass_test = Table(mass_row,names=goodsnmasscatheader)
		my_survey = galaxy_survey(source,phys_test,mass_test,photom_test)

	with open(survey+".pkl",'wb') as fp:
		pickle.dump(my_survey,fp)


def update_survey(galaxy_survey):
	return galaxy_survey(galaxy_survey.source,galaxy_survey.physpar,galaxy_survey.mass,galaxy_survey.photom)
