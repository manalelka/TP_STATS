library(sp)
library(maps)
library(microbenchmark)
library(TSP)
library(TSPpackage)
set.seed(1998)

#1.1 Longueur des chemins
#
repetitive<-rep(0,50)
nearest_insertion<-rep(0,50)
two_opt<-rep(0,50)
nearest<-rep(0,50)
branch<-rep(0,50)

for (i in 1:50) {
  sommets <- data.frame(x = runif(10), y = runif(10))
  couts <- distance(sommets)
  repetitive[i]<-TSPsolve(couts,'repetitive_nn')
  nearest_insertion[i]<-TSPsolve(couts,'nearest_insertion')
  two_opt[i]<-TSPsolve(couts,'two_opt')
  nearest[i]<-TSPsolve(couts,'nearest')
  branch[i]<-TSPsolve(couts,'branch')
}
mat <- cbind(repetitive,nearest_insertion,two_opt,nearest,branch)
colnames(mat)<-c("repetitive","near_insertion","two_opt","nearest","branch")
par(mfrow=c(1,1))
boxplot(mat,notch=TRUE)

#test
t.test(nearest-branch,NULL,"greater")

#test deux à deux
results<-as.vector(mat)
methods<-rep(colnames(mat),each=50)
pairwise.t.test(results,methods,p.adjust.method='bonferroni')

#1.2 Temps de calcul
microbenchmark<-microbenchmark(TSPsolve(couts,'repetitive_nn'),TSPsolve(couts,'branch'),TSPsolve(couts,'nearest')
               ,TSPsolve(couts,'nearest_insertion'),TSPsolve(couts,'two_opt')
               ,times=20, setup={  sommets <- data.frame(x = runif(10), y = runif(10))
               couts <- distance(sommets)})
summary(microbenchmark)


#2.1 Comportement par rapport au nombre de sommets : premier modèle
seqn <- seq(4,13,1)
temps<-matrix(0,length(seqn),10)
for(i in 1:length(seqn))
temps[i,]<-microbenchmark(TSPsolve(couts, method = 'branch'),
                      times = 10,
                      setup = { n <- seqn[i]
                      couts <- distance(cbind(x = runif(n), y = runif(n)))
                      }
)$time
par(mfrow=c(1,2)) # 2 graphiques sur 1 ligne
par(mar=c(1,1,1,1))
matplot(seqn, temps, xlab='seqn', ylab='temps')
matplot(seqn, log(temps)^2, xlab='seqn', ylab='expression(log(temps)^2)')
vect_temps <- log(as.vector(temps))^2
vect_dim <- rep(seqn,times=10)
temps.lm <- lm(vect_temps ~ vect_dim)
summary(temps.lm)
#résidu
par(mfrow=c(2,2)) # 4 graphiques, sur 2 lignes et 2 colonnes
par(mar=c(1,1,1,1))
plot(temps.lm)

shapiro.test(residuals(temps.lm))

#2.2 Comportement par rapport au nombre de sommets : étude du comportement moyen
temps.moy <- rowMeans(temps)
par(mfrow=c(1,2)) # 2 graphiques sur 1 ligne
par(mar=c(1,1,1,1))
matplot(seqn, temps.moy, xlab='seqn', ylab='temps_moy')
matplot(seqn, log(temps.moy)^2, xlab='seqn', ylab='expression(log(temps_moy)^2)')
vect_temps_moy <- log(as.vector(temps.moy))^2
temps.lm.moy <- lm(vect_temps_moy ~ seqn)
summary(temps.lm.moy)

#résidu
par(mfrow=c(2,2)) # 4 graphiques, sur 2 lignes et 2 colonnes
par(mar=c(1,1,1,1))
plot(temps.lm.moy)

shapiro.test(residuals(temps.lm.moy))

#2.3 Comportement par rapport à la structure du graphe 

data.graph <- data.frame(read.csv('DonneesTSP.csv'))

data_temps <- log(data.graph$tps)
data.graph$dim <- sqrt(data.graph$dim)
data.graph$tps<-NULL
data_temps.lm <-lm(data_temps~.,data = data.graph)
summary(data_temps.lm)

new_lm <-step(data_temps.lm)
summary(new_lm)

#résidu
par(mfrow=c(2,2)) # 4 graphiques, sur 2 lignes et 2 colonnes
par(mar=c(1,1,1,1))
plot(new_lm)

shapiro.test(residuals(new_lm))



















