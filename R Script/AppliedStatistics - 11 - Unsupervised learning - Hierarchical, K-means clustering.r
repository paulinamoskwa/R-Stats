
# Unsupervised: Hierarchical, K-means clustering

# ------------------------------------------------------------------------------------
# 0. Settings
# ------------------------------------------------------------------------------------
library(mvtnorm)
library(rgl)
library(car)
load("mcshapiro.test.RData")

### ----------------------------------------------------------------------------------
### 1. HIERARCHICAL CLUSTERING
### (generic p,g, in this case: p=4, g=3)
### ----------------------------------------------------------------------------------
# Data
data_grezzi = iris

# Dimensions
n = dim(data_grezzi)[1]
p = dim(data_grezzi)[2]-1

# Data without labels
data = data_grezzi[,1:4]

# Data'slabels
species.name = as.factor(data_grezzi[,5])
g = length(levels(species.name))

# Simple EDA
pairs(data)

## If there is no label:
## n = dim(data)[1]
## p = dim(data)[2]

# ------------------------------------------------------------------------------------
# 1.1. Compute dissimilarities
# methods: "euclidean", "manhattan", "canberra"
# ------------------------------------------------------------------------------------
iris.e = dist(data, method="euclidean")
iris.m = dist(data, method="manhattan")
iris.c = dist(data, method="canberra")

par(mfrow=c(1,3))
image(1:n,1:n,as.matrix(iris.e), main='metrics: Euclidean', asp=1, xlab='i', ylab='j' )
image(1:n,1:n,as.matrix(iris.c), main='metrics: Canberra', asp=1, xlab='i', ylab='j' )
image(1:n,1:n,as.matrix(iris.m), main='metrics: Manhattan', asp=1, xlab='i', ylab='j' )

# Comment
#   light colors = small values
#   dark colors  = large values
# ------------------------------------------------------------------------------------
# 1.2. Hierarchical clustering
# distance: "euclidean"
# linkages: "single", "average", "complete", "ward"
# ------------------------------------------------------------------------------------
iris.es = hclust(iris.e, method='single')
iris.ea = hclust(iris.e, method='average')
iris.ec = hclust(iris.e, method='complete')
iris.wa = hclust(iris.e, method="ward.D2")

# Order of aggregation
iris.es$merge
# Comment

# [128,] -119 80 (e.g.)
# It means that at the step 128 the unity 119 has been aggregated to the cluster
# produced at step 80. If both values are positive, e.g. [140,] 90 95, then we are
# aggregating the clusters produced at the step 90 and 95. If both are negative then
# we are creating a new cluster.

# Distance at which we have aggregations
iris.es$height

# Ordering that allows to avoid intersection in the dendrogram
iris.es$order

# ------------------------------------------------------------------------------------
# 1.3. Plot of the dendrograms
# ------------------------------------------------------------------------------------
par(mfrow=c(2,2))
plot(iris.es, main='euclidean-single', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
plot(iris.ec, main='euclidean-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
plot(iris.ea, main='euclidean-average', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
plot(iris.wa, main='euclidean-ward', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')

# ------------------------------------------------------------------------------------
# 1.4. Cutting the dendrogram (k=2 clusters)
# ------------------------------------------------------------------------------------
par(mfrow=c(2,2))
plot(iris.es, main='euclidean-single', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(iris.es, k=2)
plot(iris.ec, main='euclidean-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(iris.ec, k=2)
plot(iris.ea, main='euclidean-average', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(iris.ea, k=2)
plot(iris.wa, main='euclidean-ward', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(iris.wa, k=2)

# How to cut a dendrogram? (k=2)
# We generate vectors of labels through the command cutree()
cluster.ec = cutree(iris.ec, k=2)
cluster.es = cutree(iris.es, k=2)
cluster.ea = cutree(iris.ea, k=2)
cluster.wa = cutree(iris.wa, k=2)

# ------------------------------------------------------------------------------------
# 1.5. How good is the clustering?
# ------------------------------------------------------------------------------------
# Did it aggregate coherently with the dissimilarity matrix or not?

# Cophenetic Matrices
coph.es <- cophenetic(iris.es)
coph.ec <- cophenetic(iris.ec)
coph.ea <- cophenetic(iris.ea)
coph.wa <- cophenetic(iris.wa)

# Compare with dissimilarity matrix (Euclidean distance)
#layout(rbind(c(0,1,0),c(2,3,4,5)))
#image(as.matrix(iris.e), main='Euclidean', asp=1 )
#image(as.matrix(coph.es), main='Single', asp=1 )
#image(as.matrix(coph.ec), main='Complete', asp=1 )
#image(as.matrix(coph.ea), main='Average', asp=1 )
#image(as.matrix(coph.wa), main='Ward', asp=1)

# Cophenetic Coefficients
es = cor(iris.e, coph.es)
ec = cor(iris.e, coph.ec)
ea = cor(iris.e, coph.ea)
ew = cor(iris.e, coph.wa)
c("Eucl-Single"=es,"Eucl-Compl."=ec,"Eucl-Ave."=ea, "Eucl-Ward"=ew)

# Interpret the clusters (ONLY if we have the true labels)
table(label.true = species.name, label.cluster = cluster.es)
table(label.true = species.name, label.cluster = cluster.ec)
table(label.true = species.name, label.cluster = cluster.ea)
table(label.true = species.name, label.cluster = cluster.wa)

# Plot
plot(data, col=ifelse(cluster.es==1,'red','blue'), pch=19)
plot(data, col=ifelse(cluster.ec==1,'red','blue'), pch=19)
plot(data, col=ifelse(cluster.ea==1,'red','blue'), pch=19)
plot(data, col=ifelse(cluster.wa==1,'red','blue'), pch=19)

# ------------------------------------------------------------------------------------
# 1.6. If p=3 -> plot 3d
# ------------------------------------------------------------------------------------
# Data
plot3d(data, size=3, col='orange', aspect=F)
# Single linkage
plot3d(data, size=3, col=cluster.es+1, aspect=F)
# Average linkage
plot3d(data, size=3, col=cluster.ea+1, aspect=F)
# Complete linkage
plot3d(data, size=3, col=cluster.ec+1, aspect=F)
# Ward linkage
plot3d(data, size=3, col=cluster.wa+1, aspect=F)

### ----------------------------------------------------------------------------------
### 2. K-MEAN CLUSTERING (p=4, g=3)
### ----------------------------------------------------------------------------------
# Data
data_grezzi = iris

# Dimensions
n = dim(data_grezzi)[1]
p = dim(data_grezzi)[2]-1

# Data without labels
data = data_grezzi[,1:4]

# Data'slabels
species.name = as.factor(data_grezzi[,5])
g = length(levels(species.name))
pairs(data)

## If there is no label:
## n = dim(data)[1]
## p = dim(data)[2]

# ------------------------------------------------------------------------------------
# 2.1. K-means
# ------------------------------------------------------------------------------------
# We fix the number of centers
result.k = kmeans(data, centers=2)

# Labels of cluster
result.k$cluster

# Centers of clusters
result.k$centers

# Tot sum of squares
result.k$totss

# Sum of squares within clusters
result.k$withinss

# Sum(sum of squares in clusters)
result.k$tot.withinss

# Sum of squares between clusters
result.k$betweenss

# Dimension of the clusters
result.k$size

# Plot
plot(data, col = result.k$cluster+1)

# Plot 3d (if p=3)
plot3d(data, size=3, col=result.k$cluster+1, aspect = F)
points3d(result.k$centers,size=10)

# ------------------------------------------------------------------------------------
# 2.2 How to choose k
# ------------------------------------------------------------------------------------
# Evaluate the variability between the groups with respect to the variability
# withing the groups
b = NULL
w = NULL

for(k in 1:10){
  result.k = kmeans(data, k)
  w = c(w, sum(result.k$wit))
  b = c(b, result.k$bet)
}
matplot(1:10, w/(w+b), pch='', xlab='clusters', ylab='within/tot',
        main='Choice of k', ylim=c(0,1))
lines(1:10, w/(w+b), type='b', lwd=2)