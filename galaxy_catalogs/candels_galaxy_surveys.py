import numpy as np
import csv
from astropy.io import ascii
from astropy import cosmology as cosmo
from astropy.coordinates import SkyCoord  
from astropy.table import Table
import astropy.units as u
import pickle


class galaxy_survey(object):
	def __init__(self,source,physpar,mass,photom):
		self.source = source
		self.physpar = physpar
		self.mass = mass
		self.photom = photom
		#self.coords = coords
		#self.z = z
		#self.name = name
class galaxy_catalog(object):
	def __init__(self,surveys):
		if not isinstance(surveys,(list,tuple)):
			surveys = [surveys]
		for survey in surveys:
			if survey == 'goodss':
				with open(survey+'.pkl','rb') as fp:
					self.goodss = pickle.load(fp)
			elif survey == 'egs':
				with open(survey+'.pkl','rb') as fp:
					self.egs = pickle.load(fp)
			elif survey == 'cos':
				with open(survey+'.pkl','rb') as fp:
					self.cos = pickle.load(fp)
			elif survey == 'uds':
				with open(survey+'.pkl','rb') as fp:
					self.uds = pickle.load(fp)
			elif survey == 'goodsn':
				with open(survey+'.pkl','rb') as fp:
					self.goodsn = pickle.load(fp)





def output_survey(survey):
	with open(r'C:/Users/Kyle/Documents/school/research/summer18/Catalog/galaxies/'+survey+'0to3','rb') as fp:
		surveyData = pickle.load(fp)
	#print(len(surveyData[0]))
	#print(surveyData[0])
	#print(len(surveyData[0][-1]))
	#print(len(surveyData[0][-2]))
	#print(len(surveyData[0][-3]))
	
	if survey == 'goodss':
		source = 'Santini2015; goods-s_candels; https://archive.stsci.edu/prepds/candels/'
		goodssphysparheader = np.loadtxt('goodssphysparheader.txt', dtype='str',delimiter='\n')
		goodsstotmultiheader = np.loadtxt('goodsstotmultiheader.txt',dtype='str',delimiter='\n')
		goodssmasscatheader = np.loadtxt('goodssmasscatheader.txt',dtype='str',delimiter='\n')
		phys_row = np.array([surveyData[i][-2] for i in range(len(surveyData))])
		photom_row = np.array([surveyData[i][-1] for i in range(len(surveyData))])
		mass_row = np.array([surveyData[i][-3] for i in range(len(surveyData))])
		#print(phys_row.shape)
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
		#print(phys_row.shape)
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
		#print(phys_row.shape)
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
		#print(phys_row.shape)
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


#survey = 'uds'
#output_survey(survey)
uds_catalog = galaxy_catalog('uds')
#egs_catalog = galaxy_catalog('egs')
#goodss_catalog = galaxy_catalog('goodss')
#goodsn_catalog = galaxy_catalog('goodsn')

galaxy_catalogs = galaxy_catalog(['uds','egs','goodss','goodsn'])
print(galaxy_catalogs.uds.mass)