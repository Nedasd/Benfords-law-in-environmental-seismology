#!/usr/bin/python
# -*- coding: UTF-8 -*-txt')

import numpy as np
import pandas as pd
import time
import random
import datetime

P2N = 50# Positive to Negative ratio

# get positive data
df1 = pd.read_table("/Users/qizhou/Downloads/confusion matrix/0.data/2017-2019_22.txt", delimiter=",", header=0)
positive_event_list = []
for step in range(0, 22):
    positive_event_start = df1.iloc[step, 1].strip()#+":00"
    positive_event_start = datetime.datetime.strptime(positive_event_start, '%Y-%m-%d %H:%M:%S')

    event_duration = int(df1.iloc[step, 3])+1
    for a in range(0, event_duration):
        time1 =  positive_event_start + datetime.timedelta(minutes = a)
        positive_event_list.append(str(time1))
    print(step)

# get all available data from IGB02 or ILL02/ILL12
df2 = pd.read_table("/Users/qizhou/Downloads/confusion matrix/0.data/all_2017-2019(2017_ILL02).txt", sep=",",  header=0, usecols=[1])
all_event_list = []
for a in range(0, len(df2)):
    time2 = df2.iloc[a, 0]
    time2 = time2.strip()

    time2 = datetime.datetime.strptime(time2, '%Y-%m-%d %H:%M:%S')
    all_event_list.append(str(time2))

# select the negative event time/duration randomly
min_duration=20 # 1/3 hour
max_duration=10*60 # 6 hours
n_list = []
for a in range(0, P2N*22+1):
    i = 1
    while i !=0 :
        random_time = random.randint(0, len(df2) - max_duration)
        random_duration = random.randint(min_duration, max_duration)

        start_id = random_time
        end_id = random_time + random_duration

        negative_list = all_event_list[start_id:end_id]
        inter1 = list(set(negative_list) & set(positive_event_list))# do not take positive data
        inter2 = list(set(negative_list) & set(n_list))# do not duplicate data / take same data twice
        #print(len(inter1), len(inter2))
        if len(inter1) == 0 and len(inter2) == 0:
            f = open("/Users/qizhou/Downloads/confusion matrix/P1_N50_BHZ.txt", "a")
            record_data = "{}, {}, {}, {}, {}".format(a, all_event_list[start_id], all_event_list[end_id], random_duration, "none")
            f.write(record_data + "\n")
            for data in negative_list:
                n_list.append(data)
            i=0
        else:
            i=1
