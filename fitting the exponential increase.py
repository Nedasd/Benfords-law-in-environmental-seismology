#!/usr/bin/python
# -*- coding: UTF-8 -*-txt')
from scipy.stats import iqr, ks_2samp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from obspy import read, Stream, UTCDateTime, signal

#plt.rc('font',family="Arial")
plt.rcParams.update( {'font.size':7, 'font.family': "Arial"} )#, 'font.weight':'bold'
dir = "/Users/qizhou/#file/2_projects/BL_paper/fitting R2/"

# all sac is store in the driver,/Volumes/Section_4.7/1seismic_data/
sac_dir = "/Volumes/Section_4.7/1seismic_data/GFZ/ILLgraben/sac/"
sac_dir1 = "/Volumes/Section_4.7/1seismic_data/WSL/all components/"
# BL data, all features is stored here, maybe the id_clou is different
bl_dir = "/Users/qizhou/Desktop/plot/2013-2014/all_2013-2014.txt"
# this dir is for output
out_dir = "/Users/qizhou/Desktop/plot/"


def func_in(x, a, b, c):
    return a * np.exp(b * x) + c

def func_de(x, a, b, c):
    return a * np.exp(-b * x) + c


#st = read(sac_dir + "/2014/193/IGB02.14.193.14.00.00.BHZ.SAC")
#st1.plot(outfile='/Users/qizhou/#file/2_projects/BL_paper/fitting R2/1.png', size=(1200, 800))

# 2013-2014 events
df = pd.read_table("/Users/qizhou/#file/2_projects/BL_paper/fig/fitting R2/max_goodness_2013-2014.txt", sep=",", header=0)
def fit_exponential_curve(date_max_goodness, event):
    date = UTCDateTime(date_max_goodness)
    sac_name0 = "/{}/{}/IGB02.{}.{}.{}.00.00.BHZ.SAC".format(date.year, str(date.julday).zfill(3), str(date.year)[2:4],
                                                            str(date.julday).zfill(3), str(date.hour-1).zfill(2) )
    sac_name1 = "/{}/{}/IGB02.{}.{}.{}.00.00.BHZ.SAC".format(date.year, str(date.julday).zfill(3), str(date.year)[2:4],
                                                            str(date.julday).zfill(3), str(date.hour).zfill(2) )
    # how many minutes before the max goodness of fit
    r2_list = []
    a_list = []
    b_list = []
    c_list = []
    for n in range (1, 6):
        amp_iq = []
        ## get the 1s-interval interquartile range
        for step in range(0, 60 * (n + 1)):
            st = Stream()
            st += read(sac_dir + sac_name0)
            st += read(sac_dir + sac_name1)

            d1 = UTCDateTime(date_max_goodness) - 60 * n + step  # 1 min before max goodness
            d2 = UTCDateTime(date_max_goodness) - 60 * n + 1 + step
            st1 = st.trim(d1, d2, nearest_sample=False)
            m = abs(st1[0].data)
            m = iqr(m)
            # print(d1, m)
            amp_iq.append(m)
        ## pick the data
        y = amp_iq[0: amp_iq.index(max(amp_iq)) + 1]
        ## fit the data
        try:
            p0 = [5000, 0.1, -4000]
            popt, pcov = curve_fit(func_in, np.arange(0, len(y), 1), y, p0=p0, maxfev=5000)
            a = popt[0]
            b = popt[1]
            c = popt[2]
        except:
            a = 0
            b = 0
            c = 0
        fit_y = func_in(np.arange(0, len(y), 1), a, b, c)
        r2 = r2_score(y, fit_y)
        print(r2)
        r2_list.append(r2)
        a_list.append(a)
        b_list.append(b)
        c_list.append(c)

        record_data = "{},{},{},{},{},{},{}".format(event, date_max_goodness, n, r2, a, b, c)
        f = open(dir + "fit.txt", 'a')
        f.write(record_data + "\n")
        f.close()

    # find the best fit
    id = r2_list.index(max(r2_list))
    record_data = "{},{},{},{},{},{},{}".format(event, date_max_goodness, id+1, r2_list[id], a_list[id], b_list[id], c_list[id])
    f = open(dir + "fit_best.txt", 'a')
    f.write(record_data + "\n")
    f.close()

for event in range (0, 24):
    # for event
    date_max_goodness = df.iloc[event, 1]
    # for noise
    #date_max_goodness = df.iloc[event, 0]
    #fit_exponential_curve(date_max_goodness=date_max_goodness, event = event)
    print(event)


# 2017-2019 events
df1 = pd.read_table("/Users/qizhou/#file/2_projects/BL_paper/fitting R2/max_goodness_2017-2019.txt", sep=",", header=0)
def fit_exponential_curve1(date_max_goodness, event):
    date = UTCDateTime(date_max_goodness)
    if event<=3:
        sac_name = "{}/ILL02/EHZ.D_new/XP.ILL02..EHZ.D.{}.{}".format(date.year, date.year, str(date.julday).zfill(3) )
    else:
        sac_name = "{}/ILL12/EHZ.D_new/XP.ILL12..EHZ.D.{}.{}".format(date.year, date.year, str(date.julday).zfill(3))
    # how many minutes before the max goodness of fit
    r2_list = []
    a_list = []
    b_list = []
    c_list = []
    for n in range (1, 6):
        amp_iq = []
        ## get the 1s-interval interquartile range
        for step in range(0, 60 * (n + 1)):
            d1 = UTCDateTime(date_max_goodness) - 60 * n + step  # 1 min before max goodness
            d2 = UTCDateTime(date_max_goodness) - 60 * n + 1 + step
            st = read(sac_dir1 + sac_name)
            st1 = st.trim(d1, d2, nearest_sample=False)
            m = abs(st1[0].data)
            m = iqr(m)
            # print(m)
            amp_iq.append(m)
        ## pick the data
        y = amp_iq[0: amp_iq.index(max(amp_iq)) + 1]
        ## fit the data
        try:
            p0 = [100, 0.1, -1000]
            popt, pcov = curve_fit(func_in, np.arange(0, len(y), 1), y, p0=p0, maxfev=5000)
            a = popt[0]
            b = popt[1]
            c = popt[2]
        except:
            a = 0
            b = 0
            c = 0
        fit_y = func_in(np.arange(0, len(y), 1), a, b, c)
        r2 = r2_score(y, fit_y)
        print(r2)
        r2_list.append(r2)
        a_list.append(a)
        b_list.append(b)
        c_list.append(c)

        record_data = "{},{},{},{},{},{},{}".format(event,date_max_goodness, n, r2, a, b, c)
        f = open(dir + "fit.txt", 'a')
        f.write(record_data + "\n")
        f.close()

    # find the best fit
    id = r2_list.index(max(r2_list))
    record_data = "{},{},{},{},{},{},{}".format(event, date_max_goodness, id+1, r2_list[id], a_list[id], b_list[id], c_list[id])
    f = open(dir + "fit_best.txt", 'a')
    f.write(record_data + "\n")
    f.close()

for event in range (0, 21):
    # for event
    #date_max_goodness = df1.iloc[event, 1]
    # for noise
    date_max_goodness = df1.iloc[event, 0]
    fit_exponential_curve1(date_max_goodness=date_max_goodness, event = event)
    print(event)



# Visualize the data/curve 2013-2014
df = pd.read_table("/Users/qizhou/#file/2_projects/BL_paper/fitting R2/max_goodness_2013-2014.txt", sep=",", header=0)
df_event = pd.read_table("/Users/qizhou/#file/2_projects/BL_paper/fitting R2/R3/fit_best_2013-2014event.txt", sep=",", header=None)
df_noise = pd.read_table("/Users/qizhou/#file/2_projects/BL_paper/fitting R2/R3/fit_best_2013-2014noise.txt", sep=",", header=None)
def fit_exponential_visualize(date_max_goodness, event, n):
    amp_iq = []
    date = UTCDateTime(date_max_goodness)
    sac_name0 = "/{}/{}/IGB02.{}.{}.{}.00.00.BHZ.SAC".format(date.year, str(date.julday).zfill(3), str(date.year)[2:4],
                                                            str(date.julday).zfill(3), str(date.hour-1).zfill(2) )
    sac_name1 = "/{}/{}/IGB02.{}.{}.{}.00.00.BHZ.SAC".format(date.year, str(date.julday).zfill(3), str(date.year)[2:4],
                                                            str(date.julday).zfill(3), str(date.hour).zfill(2) )
    # how many minutes before the max goodness of fit
    for step in range(0, 60*(n+1) ):
        st = Stream()
        st += read(sac_dir + sac_name0)
        st += read(sac_dir + sac_name1)

        d1 = UTCDateTime(date_max_goodness) -60*n + step # 1 min before max goodness
        d2 = UTCDateTime(date_max_goodness) -60*n + 1 + step
        st1 = st.trim(d1, d2, nearest_sample=False)
        m = abs(st1[0].data)
        m = iqr(m)
        #print(d1, m)
        amp_iq.append(m)

    # fit the exponential curve by the data before max iq
    plt.scatter(np.arange(0, len(amp_iq), 1), amp_iq, label = "observed")
    y = amp_iq[0: amp_iq.index(max(amp_iq)) + 1]
    plt.scatter(np.arange(0, len(y), 1), y, label = "used data")
    #plt.show()

    try:
        p0 = [5000, 0.1, -4000]
        popt, pcov = curve_fit(func_in, np.arange(0, len(y), 1), y, p0=p0, maxfev=5000)
        a = popt[0]
        b = popt[1]
        c = popt[2]
    except:
        a = 0
        b = 0
        c = 0
    fit_y = func_in(np.arange(0, len(y), 1), a, b, c)
    r2 = r2_score(y, fit_y)
    print(r2)

    #plt.scatter(np.arange(0, len(fit_y), 1), y, color="black", label="observed", s=15)
    plt.scatter(np.arange(0, len(fit_y), 1), fit_y, color="green", alpha=0.8, label="fitting")
    plt.text(x=0, y=max(fit_y), s="y=" + str(a) + " * e^" + "(" + str(b) + " * x)+" + str(c))
    plt.text(x=0, y=min(fit_y), s="R2: " + str(r2) + "  phi: " + str(df.iloc[event, 2]) )
    plt.xlabel( "From: "+str(UTCDateTime(date_max_goodness) -60) + " unit(second)" )
    plt.ylabel( "Interquartile range of 1 second")
    plt.legend()

    plt.savefig(dir + str(event) + "_" +str(date_max_goodness) + "_raw.png", dpi=600)
    plt.show()

for event in range (0, 24):
    # for event
    #date_max_goodness = df.iloc[event, 1]
    #n = df_event.iloc[event, 2]
    # for noise
    date_max_goodness = df.iloc[event, 0]
    n = df_noise.iloc[event, 2]
    fit_exponential_visualize(date_max_goodness=date_max_goodness, event = event, n = n)
    print(event)


# Visualize the data/curve 2017-2019
df1 = pd.read_table("/Users/qizhou/#file/2_projects/BL_paper/fitting R2/max_goodness_2017-2019.txt", sep=",", header=0)
df1_event = pd.read_table("/Users/qizhou/#file/2_projects/BL_paper/fitting R2/R3/fit_best_2017-2019event.txt", sep=",", header=None)
df1_noise = pd.read_table("/Users/qizhou/#file/2_projects/BL_paper/fitting R2/R3/fit_best_2017-2019noise.txt", sep=",", header=None)
def fit_exponential_visualize1(date_max_goodness, event, n):
    amp_iq = []
    date = UTCDateTime(date_max_goodness)
    if event<=3:
        sac_name = "{}/ILL02/EHZ.D_new/XP.ILL02..EHZ.D.{}.{}".format(date.year, date.year, str(date.julday).zfill(3) )
    else:
        sac_name = "{}/ILL12/EHZ.D_new/XP.ILL12..EHZ.D.{}.{}".format(date.year, date.year, str(date.julday).zfill(3))
    # how many minutes before the max goodness of fit, n
    for step in range(0, 60*(n+1) ):
        d1 = UTCDateTime(date_max_goodness) -60*n + step # 1 min before max goodness
        d2 = UTCDateTime(date_max_goodness) -60*n + 1 + step
        st = read(sac_dir1 + sac_name)
        st1 = st.trim(d1, d2, nearest_sample=False)
        m = abs(st1[0].data)
        m = iqr(m)
        #print(m)
        amp_iq.append(m)

    # fit the exponential curve by the data before max iq
    plt.scatter(np.arange(0, len(amp_iq), 1), amp_iq, label = "observed")
    y = amp_iq[0: amp_iq.index(max(amp_iq)) + 1]
    plt.scatter(np.arange(0, len(y), 1), y, label = "used data")
    #plt.show()

    try:
        p0 = [100, 0.1, -1000]
        popt, pcov = curve_fit(func_in, np.arange(0, len(y), 1), y, p0=p0, maxfev=5000)
        a = popt[0]
        b = popt[1]
        c = popt[2]
    except:
        a = 0
        b = 0
        c = 0
    fit_y = func_in(np.arange(0, len(y), 1), a, b, c)
    r2 = r2_score(y, fit_y)
    print(r2)

    #plt.scatter(np.arange(0, len(fit_y), 1), y, color="black", label="observed", s=15)
    plt.scatter(np.arange(0, len(fit_y), 1), fit_y, color="green", alpha=0.8, label="fitting")
    plt.text(x=0, y=max(fit_y), s="y=" + str(a) + " * e^" + "(" + str(b) + " * x)+" + str(c))
    plt.text(x=0, y=min(fit_y), s="R2: " + str(r2) + "  phi: " + str(df.iloc[event, 2]) )
    plt.xlabel( "From: "+str(UTCDateTime(date_max_goodness) -60) + " unit(second)" )
    plt.ylabel( "Interquartile range of 1 second")
    plt.legend()

    plt.savefig(dir + str(event) + "_" +str(date_max_goodness) + ".png", dpi=600)
    plt.show()

for event in range (0, 21):
    # for event
    #date_max_goodness = df1.iloc[event, 1]
    #n = df1_event.iloc[event, 2]
    # for noise1
    date_max_goodness = df1.iloc[event, 0]
    n = df1_noise.iloc[event, 2]
    fit_exponential_visualize1(date_max_goodness=date_max_goodness, event = event, n=n)
    print(event)
