# for aim3,
# divide samples into 2 groups
# based on time to R.

# preliminary example
# 10 samples; 1-6: se, 7:10, pe

### y vector generation:

## divide samples into 2 groups by 
## the median of time to R
timeToR=c(3.8,4.3, 3.1,2.6,3.2,1.5,2.1,1.0,3.6,0.8)
names(timeToR)=1:10
cutoff=median(timeToR)

slow.idx=which(timeToR>cutoff)
fast.idx=which(timeToR<cutoff)

status=rep(1,10); status[fast.idx]=-1

### x feature matrix
## cli feature matrix part
gender=c(1,1,-1,1,-1,-1,1,1,-1,1)
race=c(1,1,1,1,-1,0,1,0,1,1)
age=c(16,15.8,14.3,6,17,7.3,1.9,18,13,16)

x.c=cbind(gender,race,age)