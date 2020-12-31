library(R.matlab)
library(gplots)
library(network)
library(magic)

setwd("...") #replace "..." by the directory of the BSBM toolbox

source('BSBM_Network_Source.R')


## load the outputs from BSBM_ContPnl.m

avg_nwk_Gam <- readMat("Simulation1_Gam_m.mat")

avg_nwk1 <- avg_nwk_Gam$Gam.result #the clustering probabilities
avg_nwk2 <- avg_nwk_Gam$m.result   #the network edge probabilities


conn <- avg_nwk1
conn2 <- avg_nwk2

d <- nrow(conn)


## Set thresholds for the clustering probabilities and the network edge probabilities

#Choice of thresholds is described in Section 2.4 of our paper in detail;
#The values we use in this example are calculated following the procedure we described in Section 2.4.

null_Gam = 0.2704  #threshold for the clustering probabilities
null_m = 0.3475    #threshold for the network edge probabilities


#########

## Clustering

ClusterList <- CreateClusterList(d,conn,conn2,null_m)
DifCluster <- CreateDifCluster(d,ClusterList) 
length(DifCluster[[1]])	
sort( DifCluster[[1]])
L_DifCluster <- length(DifCluster)

clusterlabel<-numeric(d)
for(k in 1:L_DifCluster){
  clusterlabel[DifCluster[[k]]]<- k
}

## Set edges for within cluster

connW <- avg_nwk1

if (L_DifCluster>1){
  for(l1 in 1:(L_DifCluster-1)){
    for(l2 in (l1+1):L_DifCluster){
      connW[DifCluster[[l1]],DifCluster[[l2]]] <-0
      connW[DifCluster[[l2]],DifCluster[[l1]]] <-0
    }
  }
}


## Set edges for between cluster

connB <- avg_nwk1
connB[connW != 0] = 0
diag(connB) <- diag(avg_nwk1)


# Set 2D coordinates of the network nodes

X<-c(1:5-3,1:5-3,1:5-3,1:5+5,1:5+5,1:5+5,1:5+5,1:5+5,1:5+5,1:5+5)
X = X+2

Y<-c(numeric(5)+28, numeric(5)+25,numeric(5)+22,numeric(5)+37, numeric(5)+34,numeric(5)+31,
     numeric(5)+19,numeric(5)+16,numeric(5)+13, numeric(5)+10)

X <- as.data.frame(X)
Y <- as.data.frame(Y)


X$label <- c(rep(0,50))
X$focus <- c(rep(0,50))


Y$label <- c(rep(0,50))
Y$focus <- c(rep(0,50))


p_trh <- null_Gam

## Within-cluster edges

x0<-NULL
y0<-NULL

x1<-NULL
y1<-NULL


num_clu <- 10

for (s1 in 1:num_clu){
  for (s2 in 1:num_clu){
    
    if (s1 == s2){  
      
      for(i in 1:d) {
        for(j in 1:d) {
          if(i!=j){
            if((conn[i,j] > p_trh) && (connW[i,j] > 0)  && (i %in% seq(from = 1, to = 50, by = 1)[clusterlabel==s1]) && (j %in% seq(from = 1, to = 50, by = 1)[clusterlabel==s2])){
              x0<-c(x0,X[j,1])
              y0<-c(y0,Y[j,1])
              x1<-c(x1,X[i,1])
              y1<-c(y1,Y[i,1])
            }
          }
        }
      }
    }
  }
}


## Between-cluster edges

z0<-NULL
w0<-NULL

z1<-NULL
w1<-NULL

for(i in 1:d) {
  for(j in 1:d) {
    if(i!=j){
      if((conn[i,j] > p_trh)  && (connB[i,j] > 0) ){
        z0<-c(z0,X[j,1])
        w0<-c(w0,Y[j,1])
        z1<-c(z1,X[i,1])
        w1<-c(w1,Y[i,1])
      }
    }
  }
}


## Assign colors for different clusters

X_cluster1 <- X[clusterlabel==2, ]
Y_cluster1 <- Y[clusterlabel==2, ]
X_cluster2 <- X[clusterlabel==3, ]
Y_cluster2 <- Y[clusterlabel==3, ]
X_cluster3 <- X[clusterlabel==1, ]
Y_cluster3 <- Y[clusterlabel==1, ]
X_cluster4 <- X[clusterlabel==4, ]
Y_cluster4 <- Y[clusterlabel==4, ]
X_cluster5 <- X[clusterlabel==5, ]
Y_cluster5 <- Y[clusterlabel==5, ]

X_cluster6 <- X[clusterlabel==6, ]
Y_cluster6 <- Y[clusterlabel==6, ]
X_cluster7 <- X[clusterlabel==7, ]
Y_cluster7 <- Y[clusterlabel==7, ]
X_cluster8 <- X[clusterlabel==8, ]
Y_cluster8 <- Y[clusterlabel==8, ]
X_cluster9 <- X[clusterlabel==9, ]
Y_cluster9 <- Y[clusterlabel==9, ]
X_cluster10 <- X[clusterlabel==10, ]
Y_cluster10 <- Y[clusterlabel==10, ]


## Make plot

plot(X$X,Y$Y,pch=16,cex=6,lwd=20,ylim=c(7,40),xlim=c(0,12),col=5,
     xlab="",ylab="",cex.main=1.5,cex.axis=1.5,xaxt="n",yaxt="n", frame.plot=FALSE)


network.arrow(x0,y0,x1,y1, curve = F, width = 0.048,
              length = 0.4, angle = 30, offset.head = 0, offset.tail = 0, col=8, border=8)
network.arrow(z0,w0,z1,w1, curve = F, width = 0.048,
              length = 0.4, angle = 30, offset.head = 0, offset.tail = 0, col=6, border=6)



points(X_cluster1[X_cluster1$focus != 1, ]$X, Y_cluster1[Y_cluster1$focus != 1, ]$Y,
       pch=16,lwd=20,col='red',cex=6)
points(X_cluster1[X_cluster1$focus == 1, ]$X, Y_cluster1[Y_cluster1$focus == 1, ]$Y,
       pch=18,lwd=20,col='red',cex=7)

if(dim(X_cluster2)[1] > 1){
  
  points(X_cluster2[X_cluster2$focus != 1, ]$X, Y_cluster2[Y_cluster2$focus != 1, ]$Y,
         pch=16,lwd=20,col='blue',cex=6)
  points(X_cluster2[X_cluster2$focus == 1, ]$X, Y_cluster2[Y_cluster2$focus == 1, ]$Y,
         pch=18,lwd=20,col='blue',cex=7)
}

if(dim(X_cluster3)[1] > 1){
  points(X_cluster3[X_cluster3$focus != 1, ]$X, Y_cluster3[Y_cluster3$focus != 1, ]$Y,
         pch=16,lwd=20,col='green',cex=6)
  points(X_cluster3[X_cluster3$focus == 1, ]$X, Y_cluster3[Y_cluster3$focus == 1, ]$Y,
         pch=18,lwd=20,col='green',cex=7)
}

if(dim(X_cluster4)[1] > 1){
  points(X_cluster4[X_cluster4$focus != 1, ]$X, Y_cluster4[Y_cluster4$focus != 1, ]$Y,
         pch=16,lwd=20,col='yellow',cex=6)
  points(X_cluster4[X_cluster4$focus == 1, ]$X, Y_cluster4[Y_cluster4$focus == 1, ]$Y,
         pch=18,lwd=20,col='yellow',cex=7)
}

if(dim(X_cluster5)[1] > 1){
  points(X_cluster5[X_cluster5$focus != 1, ]$X, Y_cluster5[Y_cluster5$focus != 1, ]$Y,
         pch=16,lwd=20,col='pink',cex=6)
  points(X_cluster5[X_cluster5$focus == 1, ]$X, Y_cluster5[Y_cluster5$focus == 1, ]$Y,
         pch=18,lwd=20,col='pink',cex=7)
}

if(dim(X_cluster6)[1] > 1){
  points(X_cluster6[X_cluster6$focus != 1, ]$X, Y_cluster6[Y_cluster6$focus != 1, ]$Y,
         pch=16,lwd=20,col='orange',cex=6)
  points(X_cluster6[X_cluster6$focus == 1, ]$X, Y_cluster6[Y_cluster6$focus == 1, ]$Y,
         pch=18,lwd=20,col='orange',cex=7)
}

if(dim(X_cluster7)[1] > 1){
  points(X_cluster7[X_cluster7$focus != 1, ]$X, Y_cluster7[Y_cluster7$focus != 1, ]$Y,
         pch=16,lwd=20,col='purple',cex=6)
  points(X_cluster7[X_cluster7$focus == 1, ]$X, Y_cluster7[Y_cluster7$focus == 1, ]$Y,
         pch=18,lwd=20,col='purple',cex=7)
}

if(dim(X_cluster8)[1] > 1){
  points(X_cluster8[X_cluster8$focus != 1, ]$X, Y_cluster8[Y_cluster8$focus != 1, ]$Y,
         pch=16,lwd=20,col='aquamarine4',cex=6)
  points(X_cluster8[X_cluster8$focus == 1, ]$X, Y_cluster8[Y_cluster8$focus == 1, ]$Y,
         pch=18,lwd=20,col='aquamarine4',cex=7)
}

if(dim(X_cluster9)[1] > 1){
  points(X_cluster9[X_cluster9$focus != 1, ]$X, Y_cluster9[Y_cluster9$focus != 1, ]$Y,
         pch=16,lwd=20,col='darkmagenta',cex=6)
  points(X_cluster9[X_cluster9$focus == 1, ]$X, Y_cluster9[Y_cluster9$focus == 1, ]$Y,
         pch=18,lwd=20,col='darkmagenta',cex=7)
}

if(dim(X_cluster10)[1] > 1){
  points(X_cluster10[X_cluster10$focus != 1, ]$X, Y_cluster10[Y_cluster10$focus != 1, ]$Y,
         pch=16,lwd=20,col='thistle',cex=6)
  points(X_cluster10[X_cluster10$focus == 1, ]$X, Y_cluster10[Y_cluster10$focus == 1, ]$Y,
         pch=18,lwd=20,col='thistle',cex=7)
}



text(X$X, Y$Y-0.6, labels=c(1:50), cex=2, col='black')


## Save the network plot

dev.copy(png, 'Simulation1_Network.png',width=10, height=10, units="in",res=500)
dev.off()


