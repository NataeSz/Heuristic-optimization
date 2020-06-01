# FIREFLY ALGORITHM
library(R6)
library('ggplot2')
library(RColorBrewer)
remove(list=ls())
citation()

# 
show_plots=TRUE # show_plots==TRUE to show iteration subplotsjesli chcemy wykresy z ustawieniem swietlikow z iteracji
escape=FALSE # escape==TRUE if fireflies can cross boudaries
max_iter<-10
lower_boundary<--10
upper_boundary<-10
liczba_swietlikow<-25 #or round((upper_boundary-lower_boundary)*2)
plot_step<-5 #round(max_iter/5) # iteration step of subplots
#
set.seed(1)
#

# func<-function(x){
#   return(x[1] + ((sin(x[1]^2-x[2]))))
# }

#Ackley function
func<-function(x){
 return(-20*exp(-0.2*sqrt((x[1]**2+x[2]**2)/2)) - exp(0.5*(cos(2*pi*x[1]) + cos(2*pi*x[2]))) + 20 + exp(1))
}

  

##########
# Class: firefly
Firefly<-R6Class(classname="Firefly",
                 public=list(
                   #pola
                   beta0=NULL,
                   beta=NULL,
                   gamma=NULL,
                   intensity =NULL,
                   function_dim=2,
                   position = NULL,
                   
                   #inicjalizacja
                   initialize=function(beta0=0.5, gamma=0.2,function_dim=2, position=NA) {
                     self$beta0=beta0 #initial and max attractiveness
                     self$beta=beta0/2 # attractiveness
                     self$gamma=gamma # absorption coef.
                     self$function_dim=function_dim 
                     self$position = runif(self$function_dim, lower_boundary, upper_boundary)
                     self$intensity=self$update_intensity()
                   },
                   update_intensity=function(){
                     if(func(self$position)==0){
                       self$intensity=1/(self$position+0.01)
                     }else{
                       self$intensity=1/func(self$position)
                     }
                   },
                     # if firefliess are not allowed to cross boundaries
                   check_boundaries=function(){
                     if(self$function_dim==1){
                       if (self$position>upper_boundary){
                         self$position=upper_boundary
                       }
                       if (self$position<lower_boundary){
                         self$position=lower_boundary
                       }
                     }else{
                       for (i in 1:self$function_dim){
                         if (self$position[i]>upper_boundary){
                           self$position[i]=upper_boundary
                         }
                         if (self$position[i]<lower_boundary){
                           self$position[i]=lower_boundary
                         }
                       }
                      }
                   }
                   
                   ))


df<-data.frame()
t<-0

################
# preprocessing 
fun_seq<-seq(lower_boundary,upper_boundary,(upper_boundary-lower_boundary)/200)
f <- expand.grid(x=fun_seq,y=fun_seq)
f$z<-apply(f[,c('x','y')],1,func)

###############
# fireflies positions plot
wykres<-function(df){ ## f, lower_boundary i upper_boundary are global variables
  print(ggplot() +
          geom_tile(data = f, aes(x = x, y = y, fill = z)) +
          stat_contour(data = f, aes(x = x, y = y,z=z),color="white", size=0.25)+
          scale_fill_gradientn(colors=(brewer.pal(6,"Greens")))+
          labs(color = "Firefly\nattractiveness", fill='Function\nvalue',size=F) +  #title=paste0("Fireflies - iteration ",max(df[,'t']))
          geom_point(df, mapping=aes(x=df[,'x'],y=df[,'y'],color=(df[,'beta']),size=(df[,'beta'])), alpha=.7)+
          xlim(lower_boundary,upper_boundary)+
          ylim(lower_boundary,upper_boundary)+
          theme_classic()+
          guides(size=F)+
          ggtitle(paste0("Fireflies - iteration ",df[1,'t'])))
}

# create firefly swarm
x=list()
for (i in 1:liczba_swietlikow){
  x[[i]]<-Firefly$new() #mozna podac wlasne parametry
  df<-rbind(df,data.frame(i=i,intensity=x[[i]]$intensity,beta=x[[i]]$beta, x=x[[i]]$position[1], y=x[[i]]$position[2],t=t))
  #print(x[[i]])
}


# initial fireflies positioning plot
wykres(df)


iteracja<- df[df[,'beta']==max(df[,'beta']),] # iteration '0'


# FireFly Algorithm

for (t in 1:max_iter){
  for (i in 1:liczba_swietlikow){
    for (j in 1:liczba_swietlikow){
      if(i!=j){
        if (func(x[[j]]$position)<func(x[[i]]$position)){ # move i firefly to the direction of j
          # count distance between them
          distance<-sqrt((x[[i]]$position[1]-x[[j]]$position[1])^2 + (x[[i]]$position[2]-x[[j]]$position[2])^2) 
          # modify attractiveness 
          x[[i]]$beta<-x[[i]]$beta0*exp(-x[[i]]$gamma*distance^2) 
          x[[i]]$update_intensity()
          
          # random vector
          u<-runif(2, lower_boundary, upper_boundary)
          # new position
          x[[i]]$position= x[[i]]$position + x[[i]]$intensity*(x[[j]]$position-x[[i]]$position) +runif(1, -0.5,0.5)
          if(escape==FALSE){
            x[[i]]$check_boundaries()
          }
        }
      }
    }
    df<-rbind(df,data.frame(i=i,intensity=x[[i]]$intensity, beta=x[[i]]$beta, x=x[[i]]$position[1], y=x[[i]]$position[2],t=t))
  }
  najlepszy=min(df[df[,'intensity']==max(df[,'intensity']),'i'])
  iteracja<-df[df[,'t']==t,]
  cat("Najlepszy swietlik w iteracji",t,":\n")
  print(iteracja[iteracja[,'intensity']==max(iteracja[,'intensity']),])
  cat("\n")
  
  # plot
  
  if(show_plots==T || t==max_iter){
    if(t%%plot_step==0){
      wykres(df[df[,'t']==t,])
    }
  }
  
  # the best firefly is moved randomly
  x[[najlepszy]]$position<- x[[najlepszy]]$position +runif(1, -0.5,0.5)
  x[[najlepszy]]$update_intensity()
}


if (max_iter<=50){
  ggplot(df, aes(t,beta)) +
    theme_classic()+
    geom_line(aes(group = t), col="snow4") +
    geom_point(aes(color = x), size=df$beta*7, alpha=.3)+
    xlab("iteration")+
    ggtitle("Attractiveness value of iterations")
}


# cost function plot
df$t <- as.factor(df$t)
temp<-do.call(rbind, lapply(split(df,df$t), function(x) {return(x[which.max(x$beta),])}))
temp$f<-apply(temp[,c('x','y')],1,func)

ggplot(temp, aes(x=as.numeric(t),y=f)) +
  geom_line(size=1, color='snow4')+
  labs(title="Cost function for the best fireflies of each iteration")+
  xlab("Iteration")+
  ylab("Cost Function Value")

cat("Najlepszy swietlik w ostatniej iteracji:\n")
print(iteracja[iteracja[,'intensity']==max(iteracja[,'intensity']),])
cat("\n")

cat("Najlepszy swietlik ze wszystkich iteracji:\n")
print(df[df[,'intensity']==max(df[,'intensity']),])
cat("\n")

iteracja$f<-apply(iteracja[,c('x','y')],1,func) #
cat(noquote("Podsumowanie ostatniej iteracji:\n"))
summary(iteracja[c('intensity', 'beta', 'x','y','f')])