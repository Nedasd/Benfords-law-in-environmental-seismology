library(eseis)
Sys.setenv(TZ='UTC')
# notice change the file directory
working_directory="/data/STORAGE/Seismic/Europe/Switzerland/Illgraben/data/seismic/sac/"


### for statistics
d0<-c(0.301,0.176,0.125,0.097,0.079,0.067,0.058,0.051,0.046)
d0<-as.numeric(unlist(d0))

### please change this
station_name<-"IGB02"
dir<-"/data/DATA/qzhou/#BL_fig/3_BL events/5.event_power_law/"
cal_step <<- 24*60
cal_windows <<- 60


cal_first_digit<-function(cal_start_time, station_name, start_date) {# cal_first_digit
  
  step<--1
  while (step<cal_step) {# while (step<cal_step)
    step<-step+1
    cal_start_time<-as.POSIXct(x =cal_start_time,TZ = "UTC")
    if (step==0){
      cal_start_time<-cal_start_time # cal_start_time, 2 hours before the event start/event_date
    } #1@
    else{
      cal_start_time<-cal_start_time+cal_windows# 4 min windows, 4*60=60s
    }#2@
    
    
    got_the_event <- aux_getevent(start = cal_start_time, 
                                  duration = cal_windows, # unit second
                                  component = "BHZ",
                                  station = station_name, 
                                  dir = working_directory)
    #### attention this value is calculated by abs
    y <- abs(got_the_event$signal)
    
    # section 1 y_max, y_min, iq ----
    iq<-round( IQR(y), 0)
    y_max<-max(y)
    y_min<-min(y)
    
    # section 2 1-9 and frequency ----
    sub_digit<- substr(y, 1, 1)
    digit<-c()
    for (a in 1:9) {
      first_digit <- sum(sub_digit == a)
      digit<-append(digit,first_digit)
    }
    
    frequency<-c()
    for (a in 1:9) {
      first_digit_f<- { (digit[a]/sum(digit) )- d0[a] }^2/ d0[a]
      frequency<-append(frequency, first_digit_f)
    }
    
    # section 3 fit of goodness, sh ----
    goodness<-{1-sqrt(sum(frequency)) }*100
    goodness<-round(goodness,3)
    
    
    # calculate the shannon entropy
    digit_frequency<-digit/sum(digit)
    shanoon_entropy<-c()
    for (a in 1:9) {
      if (digit[a]==0) {
        entropy<-0
      } else {
        entropy<-digit_frequency[a]*log2(digit_frequency[a])*(-1)
      }
      shanoon_entropy<-append(shanoon_entropy,entropy)
    }
    #print(shanoon_entropy)
    sh<-sum(shanoon_entropy)
    sh<-round(sh,3)
    
    # section 4 ks,chi ----
    d1<-as.numeric(unlist(digit/sum(digit)))
    ks<-ks.test(d0, d1, alternative = "two.sided")
    if (ks[[1]]<ks[2]) {
      ks<-paste( round(ks[[1]],3),  round(ks[[2]],3), "1", sep=",")
    } else {
      ks<-paste(round(ks[[1]],3),  round(ks[[2]],3), "0", sep=",")
    }
    try(
      chi<-chisq.test(d1, p =d0, rescale.p=TRUE, simulate.p.value = TRUE)
    )
    if (chi[1]<=15.5 & chi[3]==1 ) {
      chi_t<-paste( round(chi[[1]],3),  round(chi[[3]], 3), "1", sep=",")
    } else {
      chi_t<-paste( round(chi[[1]],3),  round(chi[[3]], 3), "0", sep=",")
    }
    
    # section 5 power slope ----
    y <- sort(y)
    y <- y[y != 0]
    x<-1:length(y)
    a<-log(x)
    b<-log(y)
    d <-data.frame(x=log(x), y = log(y))
    fit<-lm(b~a, data=d)
    
    sum_d = 0
    for (s in 1:length(y) ) {
      i = log(y[s]/y_min)
      sum_d =  sum_d +i
    }
    alpha = 1 + 12000/sum_d
    
    # section 6 write data ----
    digit<-toString(digit)
    digit_frequency<-round(digit_frequency, 2)
    digit_frequency<-toString(digit_frequency)
    
    
    cal_start_time<-as.character.Date(cal_start_time, format="%Y-%m-%d %H:%M:%S")
    step1<-paste(step, station_name, "BHZ", sep = "")
    data<-paste(step1,
                cal_start_time,
                digit,
                y_max,
                y_min,
                iq,
                goodness,
                ks, 
                chi_t,
                coefficients(fit)[1],
                coefficients(fit)[2],
                alpha, sep=",")
    txt_name<-paste(dir, "/2013/", start_date, " ", station_name,".txt", sep="")
    write.table(data,file = txt_name,sep = '\n',row.names=F, col.names=F, append=T,quote = FALSE)
    
    if (step !=1){
      txt_name<-paste(dir, "all_2013.txt",  sep="")
      write.table(data,file = txt_name,sep = '\n',row.names=F, col.names=F, append=T,quote = FALSE)
    }
    print(paste(cal_start_time, "done: ",step))
    
    
  }# while (step<cal_step) 
}#cal_first_digit


year <<- "2013"
jday1 = seq(from = 1, to = 108, by = 1)#056 incomplete data
jday2 = seq(from = 184, to = 219, by = 1)#219 incomplete data


jday = c(jday1,jday2)
jday = as.character(jday)



for (kp in 1:length(jday)) {
  
  if (nchar(jday[kp]) == 1) {
    jday[kp] = paste("00", jday[kp], sep = "")
  }
  if (nchar(jday[kp]) == 2) {
    jday[kp] = paste("0", jday[kp], sep = "")
  }
  ori_date = paste(year, "-01-01", sep = "")
  start_date = as.Date(strtoi(jday[kp], base = 10L) - 1, origin = as.Date(ori_date))
  print(start_date)
  
  cal_start_time =as.POSIXct(x =start_date,TZ = "UTC")
  
  tryCatch ( # tryCatch 
    
    expr = { 
      cal_first_digit(cal_start_time=cal_start_time, station_name=station_name, start_date=start_date)
    },
    
    error = function(e) {
      data<-paste (cal_start_time, jday[kp], "error", sep = ",") 
      print( data )
      txt_name<-paste(dir, "error_2013.txt", sep="")
      write.table(data,file = txt_name,sep = '\n',row.names=F, col.names=F, append=T,quote = FALSE)
    }
  )
  
  gc()
}