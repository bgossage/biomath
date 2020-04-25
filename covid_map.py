#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 14:24:45 2020

@author: bgossage
"""

import numpy
import pandas
import plotly
##import plotly.graph_objects

data_path = ("../COVID-19/csse_covid_19_data/csse_covid_19_time_series/"
             "time_series_covid19_deaths_US.csv")

us_data = pandas.read_csv( data_path, delimiter=',', header=0 )

print( us_data[['Lat', 'Long_']] )

print( us_data[['4/18/20']])
