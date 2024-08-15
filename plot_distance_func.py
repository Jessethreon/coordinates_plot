# -*- coding: utf-8 -*-
"""
Created on Sat Aug 27 15:56:02 2022

@author: ziqinding94
"""
import math

def my_function(la1,lo1,la2,lo2):
    
    R = 6373;
    
    point1la = la1;
    point1lo = lo1;
    point2la = la2;
    point2lo = lo2;
    deltala = math.radians(point1la-point2la);
    deltalo = math.radians(point1lo - point2lo);
    a = (math.sin(deltala/2))**2 + math.cos(math.radians(point1la)) * math.cos(math.radians(point2la))*(math.sin(deltalo/2))**2;
    c = 2*math.atan2( math.sqrt(a), math.sqrt(1-a) );
    d = R * c;

    return d