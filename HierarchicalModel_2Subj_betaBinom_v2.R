#################
### 1 factory, 2 coins
#Example for a hierarchical model where the theta parameters for binomial distributions (coin biases)
#are themselves drawn from a beta distibution (factory)
#I am trying to clearly visualize the joint distributions here

library('patchwork')
library('data.table')
library('plotly')
rm(list=ls())
#dev.off()

######################################
plot4D= function(gridDF, n_sampPoints= 1000, n_droplines= 100,info=NULL) {
  #Function to visualize three parameters and their joint probability in 3 dimensions
    #I use a visualization where the density of points in space and their color represents the
     #probability density in the region
  #input: 3D-array with probability densities as values; dimension names need to state the
   #value of the variable on the dimension
   #n_sampPoints: points sampled from distribution
  
  #input dataframe is on a x-y-z-grid => melt
  DF_3D_melt= as.data.frame(ftable(gridDF))
   #https://stackoverflow.com/questions/63311405/seeking-r-function-to-melt-5-dimensional-array-like-pivot-longer
  colnames(DF_3D_melt)= c('theta1','theta2','omega','prob')
  DF_3D_melt$theta1= as.numeric(as.character(DF_3D_melt$theta1))
  DF_3D_melt$theta2= as.numeric(as.character(DF_3D_melt$theta2))
  DF_3D_melt$omega= as.numeric(as.character(DF_3D_melt$omega))
  head(DF_3D_melt)
  
  #resample according to probability
  set.seed(123)
  sampIdx= sample(1:dim(DF_3D_melt)[1],size=n_sampPoints,replace=TRUE,prob=DF_3D_melt$prob) #sample row indices
  DF_3D_melt_resamp= DF_3D_melt[sampIdx,]
  
  #add jitter so that points do not overlap completely
  DF_3D_melt_resampNois= DF_3D_melt_resamp
  DF_3D_melt_resampNois[,1:3]= DF_3D_melt_resampNois[,1:3]+rnorm(1,0,0.1)
  DF_3D_melt_resampNois= DF_3D_melt_resampNois[
    apply(DF_3D_melt_resampNois[,1:3],1,function(x){all(x>0) & all(x<1)}),] #remove entries where parameters are out of [0,1] range due to scatter
  
  #plot with density representing probability

  #Generate data for droplines
  #https://community.plotly.com/t/droplines-from-points-in-3d-scatterplot/4113/4
  #dataframe: repeat x twice, then NA; repeat y twice, then NA; z: start, end, NA
  #NA causes paths to terminate and start new one
  set.seed(333)
  dropIdx= sample(1:dim(DF_3D_melt_resampNois)[1],n_droplines) #only subset of droplines so that points do not get hidden too much
  DF_3D_melt_dropSubset= DF_3D_melt_resampNois[dropIdx,]
  droplineDF= cbind.data.frame(
    x= as.vector(t( cbind.data.frame(DF_3D_melt_dropSubset$theta1,DF_3D_melt_dropSubset$theta1,
                                     rep(NA,length(DF_3D_melt_dropSubset$theta1))))),
    y= as.vector(t( cbind.data.frame(DF_3D_melt_dropSubset$theta2,DF_3D_melt_dropSubset$theta2,
                                     rep(NA,length(DF_3D_melt_dropSubset$theta2))))),
    z= as.vector(t( cbind.data.frame(DF_3D_melt_dropSubset$omega,
                                     rep(0,length(DF_3D_melt_dropSubset$omega)),
                                     rep(NA,length(DF_3D_melt_dropSubset$omega)))))
  )
 
  #scene = list(camera = list(eye = list(x = 1, y = .2, z = 1.25)))
  scene = list(camera = list(eye = list(x = 2.1, y = .42, z = 1.5)),
               aspectratio = list(x = 1, y = 1, z = 1),
               xaxis=list(title="\u03981",range=list(0,1),tickvals = list(0,.2,.4,.6,.8,1)),   #
               yaxis=list(title="\u03982",range=list(0,1),tickvals = list(0,.2,.4,.6,.8,1)),
               zaxis=list(title="\u03C9",range=list(0,1)),tickvals = list(0,.2,.4,.6,.8,1))
  plotly::plot_ly() %>%
    add_paths(data=droplineDF, x= ~x, y= ~y, z= ~z) %>%  #alpha=0.3, line=list(dash='dot')
    add_markers(data= DF_3D_melt_resampNois,
                x= ~theta1,y= ~theta2,
                z= ~omega, type='scatter3d',
                color=DF_3D_melt_resampNois$prob,mode='markers',
                size=0.1) %>%
    layout(showlegend= FALSE) %>%
    hide_colorbar() %>%
    layout(scene= scene, title= paste(info))
    #add_segments(x= ~theta1,y= ~theta2,z= ~omega,xend= ~theta1,yend=~theta2, zend=0)  #apparently not supported in 3D
}
######################################

#data
k1=15 #3 
n1=75 #15 
k2= 20 #4 #
n2= 25 #5 #

#omega is the mode of the beta distribution that generates the thetas
omega_vec= seq(0,1,0.02)
A_omega= 2  #parameter 'A' for beta-distribution for omega
B_omega= 2  #parameter 'B' for beta-distribution for omega
omega_prior= dbeta(omega_vec,A_omega,B_omega)
plot(omega_prior,omega_vec,xlab=expression(paste('P(',omega,')')),ylab=expression(omega),main='Prior')
plot(omega_vec,omega_prior,xlab=expression(omega),ylab=expression(paste('P(',omega,')')),main='Prior',type='l')

Kappa= 75
theta1_vec= seq(0,1,0.01)
theta2_vec= seq(0,1,0.01)

#thetas given omega
#P(theta1, theta2 | omega): P(theta1, theta2, omega)= P(theta1 | omega) * P(theta2 | omega)
thetaGivOmega_3D= array(data=NA, dim= c(length(theta1_vec),length(theta2_vec),length(omega_vec)),
                        dimnames= list(theta1_vec,theta2_vec,omega_vec))
for(i in 1:length(theta1_vec)){
  for(j in 1:length(theta2_vec)){
    for(k in 1:length(omega_vec)){
      a= omega_vec[k]* (Kappa-2)+1
      b= (1-omega_vec[k])*(Kappa-2)+1
      thetaGivOmega_3D[i,j,k]= dbeta(theta1_vec[i],a,b) * dbeta(theta2_vec[j],a,b)
    }
  }
}
thetaGivOmega_3D[1:5,1:5,1:5]
plot4D(thetaGivOmega_3D, n_sampPoints= 1000, n_droplines= 100)


#prior: P(theta1, theta2, omega)= P(theta1 | omega) * P(theta2 | omega) * P(omega)
prior_3D= array(data=NA, dim= c(length(theta1_vec),length(theta2_vec),length(omega_vec)),
                dimnames= list(theta1_vec,theta2_vec,omega_vec))
for(i in 1:length(theta1_vec)){
  for(j in 1:length(theta2_vec)){
    for(k in 1:length(omega_vec)){
      a= omega_vec[k]* (Kappa-2)+1
      b= (1-omega_vec[k])*(Kappa-2)+1
      prior_3D[i,j,k]= dbeta(theta1_vec[i],a,b) * dbeta(theta2_vec[j],a,b) * omega_prior[k]
    }
  }
}
prior_3D= prior_3D/sum(prior_3D)

prior_3D[1:5,1:5,1:5]
plot4D(prior_3D, n_sampPoints= 1000, n_droplines=100,info='PRIOR')

png(filename= paste('Priors_indiv.png'),height=12,width=4,unit='in',res=350)
par(mfrow=c(3,1))
margPrior_theta1= apply(prior_3D,1,sum)  #marginalize over theta1
print( plot(theta1_vec,margPrior_theta1,xlab='\u03981',ylab='P(\u03981)',type='l',main='Prior') )
margPrior_theta2= apply(prior_3D,2,sum)  #marginalize over theta2
print( plot(theta2_vec,margPrior_theta2,xlab='\u03982',ylab='P(\u03982)',type='l',main='Prior') )
margPrior_omega= apply(prior_3D,3,sum)  #marginalize over omega
print(  plot(omega_vec,margPrior_omega,xlab='\u03C9',ylab='P(\u03C9)',type='l',main='Prior')  )  #same as above
dev.off()


# likelihood
lik_3D= array(data=NA, dim= c(length(theta1_vec),length(theta2_vec),length(omega_vec)),
              dimnames= list(theta1_vec,theta2_vec,omega_vec))
for (i in 1:length(theta1_vec)){
  for (j in 1:length(theta2_vec)){
    for(k in 1:length(omega_vec)){
      lik= dbinom(k1,n1,theta1_vec[i]) * dbinom(k2,n2,theta2_vec[j])
      lik_3D[i,j,k]= lik
    }
  }
}
lik_3D_norm= lik_3D/sum(lik_3D)
plot4D(lik_3D_norm, n_sampPoints= 300, n_droplines=50,info='LIKELIHOOD')

apply(lik_3D,c(1,2),sum)


# posterior
post_3D= prior_3D * lik_3D
post_3D= post_3D/sum(post_3D)
post_3D[1:5,1:5,1:5]
plot4D(post_3D, n_sampPoints= 300, n_droplines=100,info='POSTERIOR')

par(mfrow=c(1,1))
margPost_theta1_omega= apply(post_3D,c(1,3),sum)  #marginalize over theta2
persp(as.numeric(rownames(margPost_theta1_omega)),as.numeric(colnames(margPost_theta1_omega)),
      margPost_theta1_omega, 
      xlab= paste('\u0398',1,sep=''), ylab='\u03C9',zlab='P(\u03981 , \u03C9)',
      main= expression(paste('P (',theta,'1 , ',omega,')')),
      theta=30,phi=50,expand=0.5,col='lightblue')

margPost_theta2_omega= apply(post_3D,c(2,3),sum)  #marginalize over theta1
persp(as.numeric(rownames(margPost_theta2_omega)),as.numeric(colnames(margPost_theta2_omega)),
      margPost_theta2_omega, 
      xlab= paste('\u0398',2,sep=''), ylab='\u03C9',zlab='P(\u03982 , \u03C9)',
      main= expression(paste('P (',theta,'2 , ',omega,')')),
      theta=30,phi=50,expand=0.5,col='lightblue')

margPost_omega= apply(post_3D,3,sum)  #marginalize over theta1 and theta2
plot(omega_vec,margPost_omega,xlab='\u03C9',ylab='P(\u03C9)',type='l',main='Posterior')
text(x=0.1,y=0.05,paste('k1 =',k1,'\nn1 =',n1,'\nk2 =',k2,'\nn2 =',n2))

margPost_theta1_theta2= apply(post_3D,c(1,2),sum)  #marginalize over omega
persp(as.numeric(rownames(margPost_theta1_theta2)),as.numeric(colnames(margPost_theta1_theta2)),
      margPost_theta1_theta2, 
      xlab= paste('\u0398',1,sep=''), ylab=paste('\u0398',2,sep=''),zlab='P(\u03981 , \u03982)',
      main= expression(paste('P (',theta,'1 , ',theta,'2)')),
      theta=30,phi=50,expand=0.5,col='lightblue')

png(filename= paste('Posteriors_indiv.png'),height=12,width=4,unit='in',res=350)
par(mfrow=c(3,1))
margPost_theta1= apply(post_3D,1,sum)  #marginalize over theta2 and omega
plot(theta1_vec,margPost_theta1,xlab='\u03981',ylab='P(\u03981 | D)',type='l',main='Posterior')
abline(v=k1/n1,col='blue',lty=2)
text(x=0.1,y=0.035,paste('k1 =',k1,'\nn1 =',n1))

margPost_theta2= apply(post_3D,2,sum)  #marginalize over theta1 and omega
plot(theta2_vec,margPost_theta2,xlab='\u03982',ylab='P(\u03982 | D)',type='l',main='Posterior')
abline(v=k2/n2,col='blue',lty=2)
text(x=0.1,y=0.03,paste('k2 =',k2,'\nn2 =',n2))

margPost_omega= apply(post_3D,3,sum)  #marginalize over theta1 and theta2
plot(omega_vec,margPost_omega,xlab='\u03C9',ylab='P(\u03C9 | D)',type='l',main='Posterior')
text(x=0.1,y=0.03,paste('k1 =',k1,'\nn1 =',n1,'\nk2 =',k2,'\nn2 =',n2))
dev.off()

#what happens if we repeatedly get inconsistent likelihoods between the two coins ?,
 #i.e. coin1 always low and coin2 always high proportion of heads
  #is the omega-distribution going to split into bimodal distribution?

data= matrix(rep(c(3,15,4,5),times=100),ncol=4,byrow=TRUE)
colnames(data)= c('k1','n1','k2','n2')
prior_3D_Loop= prior_3D  #helps avoiding errors when running script not from top
postModes_DF= matrix(NA,ncol=2,nrow=dim(data)[1]) #records marginal posterior modes of theta1 and theta2
colnames(postModes_DF)= c('theta1','theta2')
for (roundNo in 1:dim(data)[1]){
  print(paste('i=',roundNo))
  #prior
  
  #likelihood
  k1= data[roundNo,1]
  n1= data[roundNo,2]
  k2= data[roundNo,3]
  n2= data[roundNo,4]
  lik_3D= array(data=NA, dim= c(length(theta1_vec),length(theta2_vec),length(omega_vec)),
                dimnames= list(theta1_vec,theta2_vec,omega_vec))
  for (i in 1:length(theta1_vec)){
    for (j in 1:length(theta2_vec)){
      for(k in 1:length(omega_vec)){
        lik= dbinom(k1,n1,theta1_vec[i]) * dbinom(k2,n2,theta2_vec[j])
        lik_3D[i,j,k]= lik
      }
    }
  }
  lik_3D_norm= lik_3D/sum(lik_3D)
  
  #posterior
  post_3D= prior_3D_Loop * lik_3D
  post_3D= post_3D/sum(post_3D)
  post_3D[1:5,1:5,1:5]
  if (roundNo %% 5 == 0) {  #only every 5th round
    png(filename= paste('Posterior3D_round',roundNo,'.png'))
    print( plot4D(post_3D, n_sampPoints= 300, n_droplines=100, info=paste('round:',as.character(roundNo))) )
    dev.off()
  }
  
  prior_3D_Loop= post_3D #for next round
  
  png(filename= paste('Posteriors_indiv_round',roundNo,'.png'),height=12,width=4,unit='in',res=350)  
  par(mfrow=c(3,1))
  margPost_theta1= apply(post_3D,1,sum)  #marginalize over theta1
  print( plot(theta1_vec,margPost_theta1,xlab='\u03981',ylab='P(\u03981 | D)',type='l',main='Posterior') )
  abline(v=k1/n1,col='blue',lty=2)
  text(x=0.1,y=0.035,paste('round:',roundNo,'\nk1 =',k1,'\nn1 =',n1))

  margPost_theta2= apply(post_3D,2,sum)  #marginalize over theta2
  print( plot(theta2_vec,margPost_theta2,xlab='\u03982',ylab='P(\u03982 | D)',type='l',main='Posterior') )
  abline(v=k2/n2,col='blue',lty=2)
  text(x=0.1,y=0.03,paste('round:',roundNo,'\nk2 =',k2,'\nn2 =',n2))
  
  margPost_omega= apply(post_3D,3,sum)  #marginalize over omega
  print(  plot(omega_vec,margPost_omega,xlab='\u03C9',ylab='P(\u03C9 | D)',type='l',main='Posterior')  )
  text(x=0.1,y=0.05,paste('round:',roundNo,'\nk1 =',k1,'\nn1 =',n1,'\nk2 =',k2,'\nn2 =',n2))
  dev.off()
  
  postModes_DF[roundNo,1]= (theta1_vec[margPost_theta1==max(margPost_theta1)])[1]
  postModes_DF[roundNo,2]= (theta2_vec[margPost_theta2==max(margPost_theta2)])[1]
}

theta1_priorMode= theta1_vec[which(margPrior_theta1==max(margPrior_theta1))[1]]
theta2_priorMode= theta2_vec[which(margPrior_theta2==max(margPrior_theta2))[1]]
priorPostModes_DF= rbind.data.frame( c(theta1_priorMode,theta2_priorMode), postModes_DF)
print(
  plot(x=0,y=0,type='n',xlim=c(0,1),ylim= c(0,dim(priorPostModes_DF)[1]-1),xlab=expression(theta),ylab='round' ) )
lines(priorPostModes_DF[,1], 0:((dim(priorPostModes_DF)[1])-1),col='blue',type='b',pch=19 )
lines(priorPostModes_DF[,2], 0:((dim(priorPostModes_DF)[1])-1),col='red',type='b',pch=19 )  
legend(x=0.8,y=20,legend=c(expression(paste(theta,'1',sep='')),
                          expression(paste(theta,'2',sep=''))),col=c('blue','red'),lty=1)
