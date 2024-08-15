# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 13:43:18 2024

@author: ziqinding94
"""

#%%
#Define NLDN txt and ENTLN npy path
#All NLDN and ENTLN timing is compensated with propagation time delay to LOG
#Find all parameters and store them in Data_ENTLN.pickle,Data_NLDN.pickle 
#%%
#import Libaries,clear variables,close figs    
from plot_distance_func import *
from IPython import get_ipython
get_ipython().magic('reset -sf') #clear variables

from matplotlib import pyplot as plt
plt.close('all') #close figs

import numpy as np
from datetime import datetime
import time
TimeStamp_0 = time.time()


#%%
# NLDN file reading and prepreprocessing 
# Find all parameters and store them in Data_NLDN        

# from haversine import haversine
from datetime import timedelta
Para_Lightning_Speed = 299792458#m/s
Para_LOG_LatLon = (29.64229939988503, -82.347214315579) # (lat, lon)
print('NLDN data reading ...')
FileName_NLDN = "NLDN.txt"
FilePath_NLDN = "E:/streamer zone project/NLDN/"
Data_NLDN = {}
Temp_TriggerTime_Array = []
Temp_Latitude_Array = []
Temp_Longitude_Array = []
Temp_PeakCurrent_Array = []
Temp_Semi_major = []
Temp_Semi_minor = []
Temp_angle = []
Temp_StrokeType_Array = []
Temp_Distance_Array = []
f = open(FilePath_NLDN+FileName_NLDN, "r")
temp =f.read()
NLDN_txt = temp.split()
f.close()
#%%
for i in range(int(len(NLDN_txt)/10)):
    # Temp_Distance_Array.append(haversine(Para_LOG_LatLon, (float(NLDN_txt[i*6+2]),float(NLDN_txt[i*6+3]))))
    # Temp_TriggerTime_Array.append(datetime.strptime(((NLDN_txt[i*6]+' '+NLDN_txt[i*6+1][0:15])),'%Y-%m-%d %H:%M:%S.%f')+timedelta(seconds = (haversine(Para_LOG_LatLon, (float(NLDN_txt[i*6+2]),float(NLDN_txt[i*6+3])))/Para_Lightning_Speed*1000)))
    Temp_Latitude_Array.append(float(NLDN_txt[i*10+2]))
    Temp_Longitude_Array.append(float(NLDN_txt[i*10+3]))
    Temp_PeakCurrent_Array.append(float(NLDN_txt[i*10+4]))
    
    Temp_Semi_major.append(float(NLDN_txt[i*10+5]))
    Temp_Semi_minor.append(float(NLDN_txt[i*10+6]))
    Temp_angle.append(float(NLDN_txt[i*10+7]))
        
    Temp_StrokeType_Array.append(NLDN_txt[i*6+9])
    
Data_NLDN['TriggerTime'] = Temp_TriggerTime_Array
Data_NLDN['Latitude'] = Temp_Latitude_Array
Data_NLDN['Longitude'] = Temp_Longitude_Array
Data_NLDN['PeakCurrent'] = Temp_PeakCurrent_Array

Data_NLDN['SMA'] = Temp_Semi_major
Data_NLDN['SMI'] = Temp_Semi_minor
Data_NLDN['Angle'] = Temp_angle

Data_NLDN['StrokeType'] = Temp_StrokeType_Array
# Data_NLDN['Distance'] = Temp_Distance_Array

print('Done')
print('A total of ',len(Temp_Semi_major),'Pulses are found in NLDN Database')
# del temp ,NLDN_txt,Temp_TriggerTime_Array,Temp_Latitude_Array,Temp_Longitude_Array,Temp_Semi_major,Temp_Semi_minor,Temp_angle,Temp_PeakCurrent_Array,Temp_StrokeType_Array,f,i,Para_Lightning_Speed,Para_LOG_LatLon#Clear local variables 
# del FileName_NLDN,FilePath_NLDN
import pickle
with open('Data_NLDN.pickle', 'wb') as handle:
    pickle.dump(Data_NLDN, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
# del Data_NLDN,handle
#%%
#plot the ecllipse
x1 = np. zeros((1,len(Temp_Semi_major)))
y1 = np. zeros((1,len(Temp_Semi_major)))

garage = [29.642334926867328, -82.35093637620214]
garage_x1 = my_function(Para_LOG_LatLon[0],Para_LOG_LatLon[1],Para_LOG_LatLon[0],garage[1])*(garage[1]-Para_LOG_LatLon[1])/abs(Para_LOG_LatLon[1] - garage[1])
garage_y1 = my_function(Para_LOG_LatLon[0],Para_LOG_LatLon[1],garage[0],Para_LOG_LatLon[1])*(garage[0]-Para_LOG_LatLon[0])/abs(Para_LOG_LatLon[0] - garage[0])
water_facility = [29.642533671646277, -82.34994742699345]
water_facility_x1 = my_function(Para_LOG_LatLon[0],Para_LOG_LatLon[1],Para_LOG_LatLon[0],water_facility[1])*(water_facility[1]-Para_LOG_LatLon[1])/abs(Para_LOG_LatLon[1] - water_facility[1])
water_facility_y1 = my_function(Para_LOG_LatLon[0],Para_LOG_LatLon[1],water_facility[0],Para_LOG_LatLon[1])*(water_facility[0]-Para_LOG_LatLon[0])/abs(Para_LOG_LatLon[0] - water_facility[0])

for i in range(0, len(Temp_Semi_major),1):

    x1[0,i] = my_function(Para_LOG_LatLon[0],Para_LOG_LatLon[1],Para_LOG_LatLon[0],Temp_Longitude_Array[i])*(Temp_Longitude_Array[i]-Para_LOG_LatLon[1])/abs(Para_LOG_LatLon[1] - Temp_Longitude_Array[i])
    y1[0,i] = my_function(Para_LOG_LatLon[0],Para_LOG_LatLon[1],Temp_Latitude_Array[i],Para_LOG_LatLon[1])*(Temp_Latitude_Array[i]-Para_LOG_LatLon[0])/abs(Para_LOG_LatLon[0] - Temp_Latitude_Array[i])


 # Create a figure and axes

fig, ax = plt.subplots()
plt.plot(x1*1000, y1*1000,'ro')
plt.plot(garage_x1*1000,garage_y1*1000,'ks')
plt.plot(water_facility_x1*1000,water_facility_y1*1000,'ks')
plt.plot(0,0,'ks')
plt.xlabel('Distance relative to LOG (m)')
plt.ylabel('Distance relative to LOG (m)')

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse



# Create an ellipse
for i in range(0, len(Temp_Semi_major),1):
    semi_major_axis = Temp_Semi_major[i]*1000
    semi_minor_axis = Temp_Semi_minor[i]*1000
    angle = Temp_angle[i]  # in degrees
    angle = 90-angle#relative to north, convert it respect to x axis
    cx = x1[0,i]*1000
    cy = y1[0,i]*1000
# Convert angle to radians
    theta = np.deg2rad(angle)

# Parametric equation of the ellipse
    t = np.linspace(0, 2*np.pi, 100)
    x = semi_major_axis * np.cos(t)
    y = semi_minor_axis * np.sin(t)

# Rotation matrix
    x_rot = x * np.cos(theta) - y * np.sin(theta)
    y_rot = x * np.sin(theta) + y * np.cos(theta)
    x_rot += cx
    y_rot += cy
# Plot the ellipse
    plt.plot(x_rot, y_rot)
    


plt.axis('equal')
plt.show()
# plt.savefig("E:/Yanan third project/figures/LMA_check/"+name+"_"+str(t)+'.tif')




