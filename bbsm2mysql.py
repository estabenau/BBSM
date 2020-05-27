#!/usr/bin/python3
######################
# Read through BBSM salinity output binary (*.16) file
# selecting data by node number and naming output by station
# and create output formatted for easy loading to modeland via generic load
#
# fort.16 format: binary, has 13075 elements and 6857 nodes in the model domain
# so data is stored as simple series of values with subsequent
# time steps written after the previous step.
# Model outputs on 20 min. interval 
# example: the measurement_value stored as the 6858th item in the list is for element #1, time step #2
###############################

import os, sys
import numpy as np
import datetime as dt

debug = 0 # set debug flag

if len(sys.argv) != 5:
	sys.exit("\nUsage: %s model version alt test_YN\n ie %s BBSM 3 CEPP_alternative N\n " %(sys.argv[0],sys.argv[0]))

model=sys.argv[1]
version=sys.argv[2]
alt=sys.argv[3]
testyn=sys.argv[4]
mt=dt.datetime(1996,1,1,00,20,00) 		#CHECK must match start date in model


#DATA directory, data source, working directory, data requested

data_dir = "/data/models_out/%s/v%s/%s/" % (model,version,alt)
infilename = "fort_%s.16" % (alt)  
working_dir = "/data/models/BBSM/util/" 	# input location for the station list
#station_node = "serafy_stations.txt"

#########################################
# Determine element number of sites to be extracted
# from list built from modeland. Could do this at
# time of extract to make more flexible
#########################################
stations = np.genfromtxt(working_dir + "bbsm2modeland_sites.csv",delimiter=',',dtype=[('station_ID','|S10'),('station_name','|S10'),('long','<f8'),('lat','<f8'),('element','<i8'),('node1','<i8'),('node2','<i8'),('node3','<i8')])
#stations = np.genfromtxt(working_dir + station_node,delimiter=',',dtype=[('station_ID','|S10'),('node1','<i8')])
#print stations['station_ID'] # saved this as handy way to get handle on data

nnodes = 6857    # number of nodes in the model
nodes = list(range(nnodes))	# make a list
nodes = np.asarray(nodes)	# convert the list to array using numpy
nodes = nodes + 1		# add one to make array from 1 to number of nodes	

##################################
# GET THE DATA
##################################

infile = open(data_dir + infilename, mode='rb')

#steps=2 # number of 20 min. time steps that will be extracted (420480 in 16 years)
#steps=420480 # number of 20 min. time steps that will be extracted (420480 in 16 years)
if (testyn == 'Y'):
	steps=2
	print("test = %s, steps= %s" % (testyn,steps))
elif (testyn == 'N'):
	steps=420480
	#print "test = %s, steps= %s" % (testyn,steps)
else:
	print("test_YN must be either Y or N")
	exit()

for j in range(steps): # loop through model time steps
	if debug: print("mytime : ",j)
	for i in range(nnodes): # loop through all elements in each model time step
		salinity_value = float(np.fromstring(infile.read(4),dtype='float32'))
		if (i + 1) in stations['node1']:
			print("%s,%s,%s,BBSM_%s,%s-%02d-%02d,%02d%02d,%8.5f" % (model,version,alt,str(i+1),mt.year,mt.month,mt.day,mt.hour,mt.minute,salinity_value))
	mt = mt + dt.timedelta(0,1200) # add 1200 seconds (20 minutes) per time step

infile.close()
exit()
