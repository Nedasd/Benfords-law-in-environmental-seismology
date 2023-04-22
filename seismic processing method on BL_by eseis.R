library(eseis)

# get the amplitude plot
Sys.setenv(TZ='UTC')
working_directory="/data/STORAGE/Seismic/Europe/Switzerland/Illgraben/data/seismic/sac/"
d0<-c(0.301,0.176,0.125,0.097,0.079,0.067,0.058,0.051,0.046)# Theoretical values by Benford's law
d0<-as.numeric(unlist(d0))

got_the_event <- aux_getevent(start =as.POSIXct(x ="2014-07-12 13:00:00", TZ = "UTC"), 
                              duration = 6*3600, 
                              component = "BHZ",
                              station = "IGB02",
                              dir = working_directory)
#got_the_event <- signal_deconvolve(data = got_the_event, sensor = "TC120s", logger = "Cube3extBOB")
#got_the_event <- signal_demean(data = got_the_event)
#got_the_event <- signal_detrend(data = got_the_event)
#got_the_event <- signal_filter(data = got_the_event, f = c(1, 50))
#got_the_event <- signal_taper(data = got_the_event, p = 0.02)
#got_the_event <- signal_envelope(data=got_the_event, p = 0.05)


g<- c()

for (step in 0:361) {
  s1<- as.POSIXct( "2014-07-12 13:00:00", tz = "UTC")+step*60
  s2<- as.POSIXct( "2014-07-12 13:01:00", tz = "UTC")+step*60
  li<- as.POSIXct( c(s1, s2), tz = "UTC")
  sac<-signal_clip(data = got_the_event, limits = li)
  
  y <- abs(sac$signal) # when detrend, please take "y/1e10"
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
  
  goodness<-{1-sqrt(sum(frequency)) }*100
  goodness<-round(goodness,3)
  g<-append(g, goodness)
  
}

limits<- as.POSIXct( c("2014-07-12 14:40", "2014-07-12 18:00"), tz = "UTC")

png(file = paste("/data/DATA/qzhou/your_figure.png", sep=""), bg="white", width = 800, height = 400)
par(mfrow=c(1,2))

plot(got_the_event, xlab = "Time", ylab = "Amplitude(m/s)")
abline(v = limits, col="red", lwd=2)


ss = as.POSIXct(x ="2014-07-12 13:00:00",TZ = "UTC")
ee = as.POSIXct(x ="2014-07-12 19:00:00",TZ = "UTC")

x_label<-seq.POSIXt(from = as.POSIXct(ss), to = as.POSIXct(ee), by = '1 hour')
x_label<-format(x_label,"%H:%M")#"%Y-%m-%d %H:%M:%S"
plot(g, type = "l", xlab = "Time", ylab = "Fit of goodness", xaxt="n")
axis(1,at=seq(from=1, to=length(g), by=60),labels=x_label)

abline(v = c(100, 300), col="red", lwd=2)

dev.off()