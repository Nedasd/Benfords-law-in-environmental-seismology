#!/usr/bin/python
# -*- coding: UTF-8 -*-txt')
#import numba as nb
import os
import numpy as np
import pandas as pd
from scipy.stats import iqr, ks_2samp, chi2_contingency
from obspy import read, Stream, read_inventory, signal
from obspy.core import UTCDateTime


import datetime
s1 = datetime.datetime.now()

sac_dir = "/home/qizhou/#data/2017/ILL02/EHZ.D/"
output_dir = "/home/qizhou/1projects/1cal_allBL/2017/"

cal_year = 2017
cal_window = 60 # unit second
station = "ILL02"
component = "EHZ"
ruler = 1e3

BL_f = [0.301, 0.176, 0.125, 0.097, 0.079, 0.067, 0.058, 0.051, 0.046]#nb.typed.List( [0.301, 0.176, 0.125, 0.097, 0.079, 0.067, 0.058, 0.051, 0.046] )
null_data = ["Na","Na","Na","Na","Na","Na","Na","Na","Na","Na","Na","Na","Na","Na",
             "Na","Na","Na","Na","Na","Na","Na","Na","Na","Na","Na","Na","Na","Na"]# 28
file_dir = os.listdir(sac_dir)

def process_seismic_signal (sac_dir, j_day):
    j_day = str(j_day).zfill(3)
    file_name = sac_dir + "XP."+ str(station) + ".." + str(component) + ".D." + str(cal_year) +"." + str(j_day)
    stream = read(file_name)
    stream = stream.detrend("spline", order=2, dspline=3600)
    return stream

def process_seismic_signal2 (sac_dir, j_day):
    file_name = sac_dir + file_dir[j_day]
    stream = read(file_name)
    stream = stream.detrend("spline", order=2, dspline=3600)
    return stream


def cal_BL_index(data):
    ## raw iq, max, min
    # <editor-fold desc="Hide the following">
    raw_iq = float("{:.2f}".format( iqr(data) ))
    raw_max = float("{:.2f}".format( np.max(data) ))
    raw_min = float("{:.2f}".format( np.min(data) ))
    # </editor-fold>


    ## get the 1-9 observed probability, and iq, max, min
    # <editor-fold desc="Hide the following">
    if raw_min >= ruler: # move the data to the same base line
        data1 = data - (raw_min - ruler)
    else:
        data1 = data + (ruler - raw_min)
    iq = float("{:.2f}".format( iqr(data1) ))
    max = float("{:.2f}".format( np.max(data1) ))
    min = float("{:.2f}".format( np.min(data1) ))

    amp_data = pd.DataFrame( data1 )
    amp_data = amp_data.astype(str)

    d = (amp_data.iloc[:, 0]).str[0: 1]
    d = list(d)

    digit_l = []
    for digit in range(1, 10):
        first_digit = d.count(str(digit))
        digit_l.append(first_digit)
    digit_f = digit_l / np.sum(digit_l)
    digit_f = [float('{:.3f}'.format(i)) for i in digit_f]
    # </editor-fold>

    # get goodness, ks, chi-squared, alpha
    # <editor-fold desc="Hide the following">
    frequency = []
    for a in range(0, 9):
        first_digit_f = pow((digit_f[a] - BL_f[a]), 2) / BL_f[a]
        frequency.append(first_digit_f)
    goodness = ( 1 - pow(sum(frequency), 0.5) ) * 100
    goodness = float( "{:.3f}".format(goodness) )
    #digit_f_nb = nb.typed.List(digit_f)
    #goodness = cal_BL_index_phi (digit_f_nb)
    #goodness = float("{:.3f}".format(goodness))

    ks = ks_2samp(BL_f, digit_f, alternative='two-sided', method='exact')
    ks = float( "{:.3f}".format(ks[1]) )# pvalue
    chi = chi2_contingency( np.array([BL_f, digit_f]) )
    chi = float( "{:.3f}".format(chi[1]) )

    sum_d = []
    y_min = np.min( data1 )
    for s in range (0, len(data1) ) :
        i = np.log(data1[s] / y_min)
        sum_d.append( i )
    alpha = 1 + len(data1) / np.sum( sum_d )
    alpha = float("{:.2f}".format(alpha))
    #data1_nb = nb.typed.List(data1)
    #alpha = cal_BL_index_alpha (data1_nb)
    #alpha = float("{:.2f}".format(alpha))
    # </editor-fold>

    trans = digit_l + digit_f + [raw_max, raw_min, raw_iq, max, min, iq, goodness, ks, chi, alpha]
    return trans



for event_number in range ( 0, len(file_dir) ):
    print(j_day)
    try: # there are have seismic data
        st = process_seismic_signal2(sac_dir=sac_dir, j_day=event_number)
        df = pd.read_table("/home/qizhou/#data/title.txt", sep=",", header=None)
        j_day = ( st[0].stats.starttime ).julday
        for step in range(0, 1440):
            d1 = UTCDateTime(year=cal_year, julday=j_day) + step * cal_window
            d2 = UTCDateTime(year=cal_year, julday=j_day) + cal_window + step * cal_window
            tr = st.copy()
            try: # hava the data except this minute
                tr1 = tr.trim(starttime=d1, endtime=d2, nearest_sample=False)
                data = abs(tr1[0].data)
                out_BL_data = cal_BL_index2(data)
            except:
                out_BL_data = null_data
            all_data = [str(step) + str(station) + str(component), str(d1.strftime('%Y-%m-%d %H:%M:%S'))] + out_BL_data
            df1 = pd.DataFrame(all_data).T
            df = pd.concat([df, df1])
        # df.reset_index()
        df.to_csv(output_dir + str(UTCDateTime(year=cal_year, julday=j_day).strftime('%Y-%m-%d')) + ".txt", index=False,header=False)
        df.to_csv(output_dir + "2017" + "_all.txt", index=False, header=False, mode='a')
    except: # there are do have data
        df = pd.read_table("/home/qizhou/#data/title.txt", sep=",", header=None)
        for step in range(0, 1440):
            d1 = UTCDateTime(year=cal_year, julday=j_day) + step * cal_window
            d2 = UTCDateTime(year=cal_year, julday=j_day) + cal_window + step * cal_window
            out_BL_data = null_data
            all_data = [str(step) + str(station) + str(component), str(d1.strftime('%Y-%m-%d %H:%M:%S'))] + out_BL_data
            df1 = pd.DataFrame(all_data).T
            df = pd.concat([df, df1])
        df.to_csv(output_dir + str(cal_year) + "_all.txt", index=False, header=False,mode='a')
    print( "done: " + str(UTCDateTime(year=cal_year, julday=j_day).strftime('%Y-%m-%d')) )

