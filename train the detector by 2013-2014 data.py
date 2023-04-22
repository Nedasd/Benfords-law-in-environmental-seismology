#!/usr/bin/python
# -*- coding: UTF-8 -*-txt')
import numpy as np
import pandas as pd
import platform
import numba as nb


mean_number = 20 # not used

true_positive = []
false_positive = []
true_positive_wsl = []

if platform.system() == 'Darwin':
    dirs = "/Users/qizhou/Downloads/confusion matrix/test/"
elif platform.system() == 'Linux':
    dirs = "/home/qizhou/1projects/2confusion_matrix/2_new_method/training/"### please change here
else:
    print("set your file dir")

@nb.jit(nopython=True)
def pick_event_sub1(diff_slope, iq, date, len_date, threshold_d, threshold_iq_ratio, threshold_slope):

    conforming_date = nb.typed.List( ["Nothing"] )
    conforming_slope = [0]

    step = 20#mean_number
    while step < (len_date - threshold_d):
        condition1 = diff_slope[step]
        # (iq[step] / np.mean(iq[step - mean_number:step]))
        condition2 = iq[step] / ( sum(iq[step - mean_number:step]) / len(iq[step - mean_number:step]) )
        # np.mean(diff_slope[step:step + threshold_d])
        condition3 = sum(diff_slope[step:step + threshold_d]) / len(diff_slope[step:step + threshold_d])

        if condition1 <= threshold_slope and condition2 >= threshold_iq_ratio and  condition3 <= threshold_slope:
            print(condition3)
            print(date[step:step + threshold_d])
            conforming_slope.append( condition3 )
            conforming_date = date[step:step + threshold_d]
            step = step + threshold_d
        else:
            step = step + 1
    return (conforming_date[1:], conforming_slope[1:])

# find the event by BL
def pick_event1(output_space, name, start_id, end_id, threshold_slope, threshold_d, threshold_iq_ratio, test_number):

    ## load the data base, all processed data is store here, both 2013 and 2014
    diff_slope_nb = diff_slope_nb_all[ int(start_id) : int(end_id+1) ]
    iq_nb = iq_nb_all[ int(start_id) : int(end_id+1) ]
    date_nb = date_nb_all[ int(start_id) : int(end_id+1) ]
    goodness_nb = goodness_nb_all[ int(start_id) : int(end_id+1) ]
    len_date = end_id-start_id+1

    output = pick_event_sub1(diff_slope=diff_slope_nb, iq=iq_nb, date=date_nb, len_date=len_date,
                             threshold_d=threshold_d, threshold_iq_ratio=threshold_iq_ratio, threshold_slope=threshold_slope)
    conforming_date = output[0]
    conforming_slope = output[1]
    if len(conforming_date) != 0:
        # make the date short
        conforming_date2 = []
        for add_date in conforming_date:
            conforming_date2.append(add_date[11:16])

        record_data = "threshold_slope:{},{},{},{},{},{},{}".format( name, conforming_date[0], max(goodness_nb), conforming_date[-1], max(iq_nb), str(conforming_date2), str(conforming_slope) )
        # <editor-fold desc="record data">
        f = open(output_space, 'a')
        f.write(record_data + "\n")
        f.close()
        # </editor-fold>
        # <editor-fold desc="calculate confusion matrix">
        if test_number <= 9:  # the first 24 number is TP
            true_positive.append(1)
            true_positive_wsl.append(1)
        elif test_number <= 23:
            true_positive.append(1)
        else:
            false_positive.append(1)
        # </editor-fold>


### calculate the True positive number and True nagative number,
def calculate_TP_TN (con_p, n_t, name):
    TP = np.sum(true_positive)
    FP = np.sum(false_positive)

    FN = con_p - TP
    TN = con_p*n_t - FP
    TPR = format(TP/con_p, ".3f")
    FPR = format(FP/(con_p*n_t), ".3f")
    F1_score = format( 2*TP/(2*TP+FP+FN), ".3f" )
    TS = format( TP/(TP+FP+FN), ".3f" )
    confusion_matrix = "con_p,{},con_n,{},WSL_TP,{},TP,{},TN,{},FP,{},FN,{},TPR,{},FPR,{},TS,{},F1_score,{}".format(con_p,con_p*n_t,sum(true_positive_wsl),TP,TN,FP,FN,TPR,FPR,TS,F1_score)

    f = open(output_dir + "#confusion_matrix_2013_2014.txt", 'a')
    f.write(str(name)+","+confusion_matrix + "\n")
    f.close()



output_dir = dirs + "round/"
con_positive_negative_space = dirs +"#P24_N1200_2013-2014.txt"
data_base = dirs + "2013_2014.txt"
# the testing parameter is stored here
df_p = pd.read_table(dirs + "training_parameters.txt", delimiter=",", header=0)

# the human marked events were stored here, con_positive_negative_space
df = pd.read_table(con_positive_negative_space, delimiter=",", header=0)
df1 = pd.read_table(data_base, delimiter=",", header=0, low_memory=False)  # all data is stored here
event_date_all = np.array(df1.iloc[:, 1])## all the date-time is stored here

### process the data
diff_slope_all = df1.iloc[:, 29]
diff_slope_all = diff_slope_all.replace("Na", 10)
diff_slope_all = np.array(diff_slope_all.replace("inf", 10), dtype=float)

iq_all = df1.iloc[:, 25]
iq_all = np.array(iq_all.replace("Na", 100), dtype=float)

goodness_all = df1.iloc[:, 26]
goodness_all = np.array(goodness_all.replace("Na", -52.39), dtype=float)

diff_slope_nb_all = nb.typed.List(diff_slope_all)
iq_nb_all = nb.typed.List(iq_all)
goodness_nb_all = nb.typed.List(goodness_all)
date_nb_all = nb.typed.List(event_date_all)


for p in range(0, len(df_p) ):#len(df_p)
    del true_positive[:]
    del false_positive[:]
    del true_positive_wsl[:]

    threshold_slope = df_p.iloc[p, 0]
    threshold_d = df_p.iloc[p, 1]
    threshold_iq_ratio = df_p.iloc[p, 2]

    name = "threshold_slope_{},threshold_d_{},threshold_iq_ratio_{}".format(threshold_slope, threshold_d,threshold_iq_ratio)
    store_space_output = output_dir + "detail1/" + str(name) + ".txt"

    for test_number in range(0, len(df) ):
        start = df.iloc[test_number, 1]
        end = df.iloc[test_number, 2]
        start_id = np.where(event_date_all==start)[0]   - mean_number
        end_id = np.where(event_date_all==end)[0]

        pick_event1(output_space=store_space_output, name=name, start_id=start_id, end_id=end_id,
                   threshold_slope=threshold_slope, test_number=test_number, threshold_d=threshold_d,
                   threshold_iq_ratio=threshold_iq_ratio)
        #print(str(test_number) + str(name))
    # </editor-fold>

    ## 2. calculate confusion matrix
    calculate_TP_TN(con_p=24, n_t=50, name=name)
    print(str(name))
