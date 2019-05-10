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
from candels_galaxy_surveys import galaxy_survey,galaxy_catalog



def update_survey(survey):
	return galaxy_survey(survey.source,survey.physpar,survey.mass,survey.photom)