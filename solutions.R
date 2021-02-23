

# /////////////////////////////   Question 1   ////////////////////////////////////



#define the chi squared function curve
curve(dchisq(x, df = 8), from = 0, to = 30,
      main = expression("X~χ"[v]^"2"),
      cex.main = 1.5, bty="l",xaxt="n",yaxt="n",ylab="", xlab="x",xaxs="i", yaxs="i", lwd=1.2) 

axis(1, at=0, "0") # zero to x axis
axis(1, at=13, "a") # a to x axis
axis(1,at=20, "b") # b to x axis


#create vector of x values
x_vector <- seq(13,20)

#create vector of chi-square density values
p_vector <- dchisq(x_vector, df = 8)

#fill in portion of the density plot from 13 to 20
polygon(c(x_vector, rev(x_vector)), c(p_vector, rep(0, length(p_vector))),
        col = adjustcolor('blue', alpha=1.0), border = 1)

#adding text to the plot
text(20, 0.09, expression(f(x) == frac(1, 2^frac(v,2)*gamma*(frac(v,2)))
                          *x^(frac(v,2)*-1)*e^-frac(x,2)*' ,   '*x*'>'*0),cex = 1.2)
text(25, 0.02, expression(P('a'*'<'*X*'<'*'b') == integral(f(x)*d(x),'b','a')),cex = 1.2)

#Adding arrow to the plot
arrows(19.5, 0.02, 16, 0.0145,col = 'red',lwd = 2) 




# /////////////////////////////   Question 2   ////////////////////////////////////

library(bmp)

#a) Calculate the average of the three faces, and plotting an image of it.

#Read bmp images
p1.bmp <- read.bmp("/Volumes/data/Maths & Stats/img1.bmp")
p2.bmp <- read.bmp("/Volumes/data/Maths & Stats/img2.bmp")
p3.bmp <- read.bmp("/Volumes/data/Maths & Stats/img3.bmp")
avg <- (p1.bmp + p2.bmp + p3.bmp)/3

#Rotate & Plot the average face
image(t(apply(avg,2,rev)),col = gray((0:32)/32),axes=F)

#b) plot images of the difference of each face from the average face on one figure.

par(mfrow=c(2,2))

#Rotate the images and plot the difference faces
differ1<- p1.bmp-avg
image(t(apply(differ1,2,rev)),col = gray((0:32)/32),axes=F)
differ2 <- p2.bmp-avg
image(t(apply(differ2,2,rev)),col = gray((0:32)/32),axes=F)
differ3 <- p3.bmp-avg
image(t(apply(differ3,2,rev)),col = gray((0:32)/32),axes=F)

#c) Based on the covariance matrix of the differences, calculate eigenfaces (equal to eigenvec- tors) 
#and plot images of the first 3 eigenfaces.

#bind all the difference faces
difference<- rbind(as.vector(differ1),as.vector(differ2),as.vector(differ3))

#Calculate covariance matrix and eigen vectors
covmatrix<- cov(difference)
eig<- eigen(covmatrix)
eigfaces<- eig$vectors

#Retreive First 3 Eigen Faces
eigface1<- matrix(eigfaces[,1], nrow=51, byrow=TRUE)
eigface2<- matrix(eigfaces[,2], nrow=51, byrow=TRUE)
eigface3<- matrix(eigfaces[,3], nrow=51, byrow=TRUE)

# Plot First 3 Eigen Faces 
par(mfrow=c(2,2))
image(t(apply(t(eigface1),2,rev)), col=gray((0:32)/32), axes=F)
image(t(apply(t(eigface2),2,rev)), col=gray((0:32)/32), axes=F)
image(t(apply(t(eigface3),2,rev)), col=gray((0:32)/32), axes=F)





# /////////////////////////////   Question 3   ////////////////////////////////////



A<-cbind(c(1,1,1,1),c(2,0,0,0)) 
rslt<- svd(A)
rslt



# /////////////////////////////   Question 4   ////////////////////////////////////


par(mfrow=c(1,1))
library(deSolve)
#Differential Equation System
Ovini <-c(y=0,y1=1, y2=(-5/2))
derivs.T <- function (t,Ov,parms){
  with(as.list(c(Ov,parms)), {
    dy<-y1
    dy2 <-y2
    dy3 <- (exp(-t)-dy2+dy+y)
    list(c(dy,dy2, dy3))
  })
}

#Set the time interval
times <- seq(from = 0, to = 5, by = 0.01)

#Solve and Plot the numerical Solution
out.T <- ode(y = Ovini, times = times, func = derivs.T, parms = NULL)
plot(out.T[,"time"],out.T[,"y"], type = "l",xlab="time",ylab="y",col="green",lwd=2,
     main = "Theoretical and numerical solutions for t (0, 5).")

#Plot of Theoretical Solution
Tmp<-function(t) -(1/4*exp(-t)*(t^2))+((1)*exp(-t)*t)
cuv<- curve(Tmp,0,12,add=TRUE,lty=2,col="red")
legend("topright", c("numerical", "theoretical"), lty = c(1,2), col = c("green","red"), box.lwd = 0)





# /////////////////////////////   Question 5   ////////////////////////////////////


library(deSolve)
#(i) Write R code to solve this system of equations.

#Creating the System of Equation
LotVmod <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx<- (a*x)-(g*x^2)-(b*x*y)
    dy <- (-c*y)+(d*x*y)
    return(list(c(dx, dy)))
  })
}

#Run you code to obtain the numerical solution for t in the interval (0,5) 
Pars <- c(a = 5, b = 0.01, c=100, d=0.01, g=0.0001)
State <- c(x = 10000, y = 60)
Time <- seq(0, 5, by = 0.1)
out <- as.data.frame(ode(func = LotVmod, y = State, parms = Pars, times = Time))

#Plot showing how the human population change with the fish population
matplot(out[,"x"],out[,"y"],main="Human population change w.r.t Fish population",
        type = "l",xlab="Fish population",ylab="Human population",col="black",lwd=1)

#Plot showing how both population change over time
matplot(out[,1],out[,2:3],main="Both population change over Time",
        type = "l", xlab = "time", ylab = "population")
legend("topright", c("Fish", "Human"), lty = c(1,2), col = c(1,2), box.lwd = 0)


#(iii) Find the values of x and y for which the system achieves equilibrium.
Pars <- c(a = 5, b = 0.01, c=100, d=0.01, g=0.0001)
State <- c(x = 10000, y = 60)
Time <- seq(0, 14, by = 0.1)
out <- as.data.frame(ode(func = LotVmod, y = State, parms = Pars, times = Time))
#Plot showing how both population change over time
matplot(out[,1],out[,2:3],main="Both population change over Time", 
        type = "l", xlab = "time", ylab = "population")
legend("topright", c("Fish", "Human"), lty = c(1,2), col = c(1,2), box.lwd = 0)
#Looking at the plot it seems that system achieves equilbrium at t=12

#Finding x at which system achieves equilibrium
Pars <- c(a = 5, b = 0.01, c=100, d=0.01, g=0.0001)
State <- c(x = 10000, y = 60)
Time <- seq(0, 12, by = 0.01)
out <- as.data.frame(ode(func = LotVmod, y = State, parms = Pars, times = Time))

#Plot showing how both population change over time
matplot(out[,1],out[,2], main="Fish population change over Time",
        type = "l", xlab = "time", ylab = "population")
legend("topright", c("Fish","x"), lty = 1, col = c("black","green"), box.lwd = 0)
abline(h=10000, col="green", lwd = 3)

#Finding y at which system achieves equilibrium
Pars <- c(a = 5, b = 0.01, c=100, d=0.01, g=0.0001)
State <- c(x = 10000, y = 60)
Time <- seq(0, 12, by = 0.01)
out <- as.data.frame(ode(func = LotVmod, y = State, parms = Pars, times = Time))

#Plot showing how both population change over time
matplot(out[,1],out[,3],main="Human population change over Time",
        type = "l", xlab = "time", ylab = "population")
legend("topright", c("Human","y"), lty = 1, col = c("black","green"), box.lwd = 0)
abline(h=400, col="green", lwd = 3)

#So value of x=10000 and y=400





# /////////////////////////////   Question 6   ////////////////////////////////////


#Creating and Producing the plot of the function
func <- function(x,y) (x^2 + y - 11)^2 + (x + y^2 - 7)^2
x<-seq(-4.5,4.5,length.out = 35) 
y<-seq(-4.5,4.5,length.out =35)
z <- outer(x,y,func)

#plotting the function
persp(x,y,z,theta =-45,phi=45,expand=0.5,col="cyan",shade =0.70,ticktype ="detailed",xlab ="x",ylab ="y",zlab="z")
image(x,y,z)
contour(x,y,z)


func2<-function(x) (x[1]^2 + x[2] - 11)^2 + (x[1] + x[2]^2 - 7)^2
# All Minimums
min1<- optim(c(-4,-4),func2,hessian = T)$par
min2<- optim(c(2,-2),func2,hessian = T)$par
min3<- optim(c(2,2),func2,hessian = T)$par
min4<- optim(c(4,4),func2,hessian = T)$par

minval<- rbind(min1,min2,min3,min4)
minval

#All Maximums
max<- optim(c(0,-0),func2,hessian = T,control=list(fnscale=-1))$par
max




# /////////////////////////////   Question 7  ////////////////////////////////////



#creating Probability function
prob <- function(n){
  loop<-seq(1,100000,1) # Define repitition
  object<-seq(1,n,1) # Actors sequence
  cnt<-0 
  for(k in loop) # for all iterations
  {
    child<-sample(object) # random samples for actors babies photographs 
    for(i in 1:length(object)) # for all the n values
    {
      if(object[i] == child[i]) # if the actor matches the baby at the same index
      {
        cnt<-cnt+1 # increment the count
        break
      }
    }
  }
  return(cnt / length(loop)) # average number of matches
}

run <-seq(2,15,1) # create n sequence vector
results<-vector() 
for(val in run) # for all the test values
{
  results<-c(results, prob(val)) # run the experiment
}

plot(run, results, col = "blue", 
     main = "Probabilities for N 2 to 15", xlab = "n", 
     ylab = "Probability of atleast one", type = "b", lwd = 3) # plot the results





# /////////////////////////////   Question 8   ////////////////////////////////////


# function for generating m datasets with size nn from a binomial distributions
creator<-function(nn,m)
{
  p<-0.01 # set p value
  t<-10 # set the trial size
  bidist<-matrix(rbinom(m*nn, t, p), nrow=nn, ncol=m) # binomial distrbution function to generate matrix
  gen<-apply(bidist, 2, mean) # apply method generate means
  return(gen)
}

# 10000 means from two modes
mean_sam<-10000 # initialize number of means
modes<-c(20,100) # initialize group sizes
means<-lapply(modes,creator,mean_sam) # generate the means

par(mfrow=c(1,length(means)))
# Plot both histograms
for(val in seq(1, length(modes),1)) 
{
  hist(means[[val]], main=paste("Sample Size", modes[val]), xlab="means", col = c("blue","cyan"))
}
# unique means percentage in samples, 
for(val in seq(1,length(modes),1))
{
  out<-(length(unique(means[[val]])) / modes[val] * 100)
  print(out)
}






# /////////////////////////////   Question 9   ////////////////////////////////////

par(mfrow=c(1,1))
#initializing x with random sample of size 20 mean 100 and SD 4
x <- rnorm(20,100,4)

#initializing e with random sample of size 20 mean 0 and SD 1
e <- rnorm(20,0,1)

#initializing y of 20 size
y<- vector()

for(i in 1:length(x))
{
  y[i]<- 2 + x[i] + e[i]
}

#plot x and y on the same graph
plot(x,y,type="p",col="black", xlab = "samples", ylab = "observations"
     , main = "Plot of points x & y")

#plot the lines showing means of x and y
abline(v = mean(x), col="blue", lwd=3, lty=1)
abline(h = mean(y), col="black", lwd=3, lty=1)
legend("topright", c("y-mean","x-mean"), lty = 1, 
       col = c("black","blue"), box.lwd = 1)

#performing 2 sample t test
tsam<- t.test(x,y, conf.level = 0.95)
tsam
#performing matched pairs t test
ptsam <- t.test(x,y, paired = TRUE)
ptsam






# /////////////////////////////   Question 10   ////////////////////////////////////

# read the data
cheese<-read.table("/Volumes/data/Maths & Stats/cheese.txt", header=TRUE, sep="\t")
colnames(cheese)<-abbreviate(colnames(cheese)) # Abbreviate Column names
head(cheese)

# Part i
processed<-cheese[-1] # Remove name data
cheese.pca<-prcomp(processed,scale=TRUE) # Peform PCA
cheese.pca
plot(cheese.pca, xlim=rev(c(1,12)), col = "#ECEF95", 
     main = "Singular Values") # Scree plot lowest to highest

# Part ii
biplot(cheese.pca, main = "Affect of nutrients on PCA") # Biplot how nutrients affect PCA

# Part iii
summary(cheese.pca) # Get the summary
# 10 components

# Part iv
component<-cheese.pca$rotation[,1] # Get First Principal Component
component
# Cheddar on index 6, Edam on index 15
calcScore<-function(A, B) # Calculate score function
{
  total<-0 
  for(i in 1:length(A)) # for all the elements in the vector
  {
    total<-total + (A[i] * B[i]) # times each element by the same element in the principle component score vector
  }
  return(total) # return the score
}
cheddar<-cheese[6,-1]
edam<-cheese[15, -1]
cheddarScore<-calcScore(component, cheddar)
edamScore<-calcScore(component, edam)
cheddarScore
edamScore
