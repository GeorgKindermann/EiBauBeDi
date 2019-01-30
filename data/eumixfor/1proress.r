#https://doi.org/10.5061/dryad.8v04m
#https://doi.org/10.1007/s13595-017-0660-z

temp <- tempfile()
download.file("https://datadryad.org/bitstream/handle/10255/dryad.146072/EuMIXFOR%20Scots%20pine%20-%20European%20beech%20data..zip?sequence=1",temp)
datT <- read.csv(unz(temp, "Trees.txt"))
datR <- read.csv(unz(temp, "Cores.txt"))
unlink(temp)

datT <- datT[datT$Triplet == 1042 & datT$Plot == "pibe",]
datR <- datR[datR$Triplet == 1042 & datR$Plot == "pibe",]

t1 <- with(datR[datR$year > 2008,], aggregate(list(rw=rw), by=list(tree=Nr), FUN=mean))
t2 <- with(datT, data.frame(tree=Nr,x,y,z=0,species,d=dbh,h))
me <- merge(t2, t1, all.x=T)
me$species <- as.character(me$species)
me$species[me$species == "Fagus sylvatica"] <- "FASY"
me$species[me$species == "Pinus sylvestris"] <- "PISY"
me$species <- as.factor(me$species)

#Estimate height increment using yield table curves
me$ih5 <- NA
t <- 40
c1 <- -10.6; c2 <- 2.82; c3 <- 0.598
a <- nls(h ~ c0*(log(1 + exp(c1)*t^c2))^c3, data=me[me$species=="FASY",], start=list(c0=20))
me$ih5[me$species=="FASY"] <- (coef(a)[1]*(log(1 + exp(c1)*(t-5)^c2))^c3) / (coef(a)[1]*(log(1 + exp(c1)*(t)^c2))^c3)
c1 <- -13.3; c2 <- 3.31; c3 <- 0.316
a <- nls(h ~ c0*(log(1 + exp(c1)*t^c2))^c3, data=me[me$species=="PISY",], start=list(c0=20))
me$ih5[me$species=="PISY"] <- (coef(a)[1]*(log(1 + exp(c1)*(t-5)^c2))^c3) / (coef(a)[1]*(log(1 + exp(c1)*(t)^c2))^c3)

#Estimate missing tree heights
with(me, plot(d, h, col=species, xlim=range(me$d)))
c3 <- 3
a1 <- nls(h ~ 1.3 + c0*log(1 + exp(c2)*(d)^c3)^c1, data=me[me$species=="FASY",], start=list(c0=25,c1=0.2,c2=-10), trace=T, algo="port", control= nls.control(maxiter = 99, warnOnly = T))
c3 <- 3
a2 <- nls(h ~ 1.3 + c0*log(1 + exp(c2)*(d)^c3)^c1, data=me[me$species=="PISY",], start=list(c0=25,c1=0.2,c2=-10), trace=T, algo="port", control= nls.control(maxiter = 999, warnOnly = T))
fun <- function(x) {predict(a1, list(d=x))}
curve(fun, from=0, to=max(me$d), add=T)
fun <- function(x) {predict(a2, list(d=x))}
curve(fun, from=0, to=max(me$d), add=T, col=2)
me$hm <- NA
me$hm[me$species=="FASY"] <- round(predict(a1, newdata=me[me$species=="FASY",]),1)
me$hm[me$species=="PISY"] <- round(predict(a1, newdata=me[me$species=="PISY",]),1)
me$h[is.na(me$h)] <- me$hm[is.na(me$h)]
me$hm <- NULL

#Estimate id5 for all not measured trees
me$rwm <- NA
library(randomForest)
set.seed(0)
a <- randomForest(rw ~ d + h, data=me[me$species=="FASY",], na.action=na.omit)
me$rwm[me$species=="FASY"] <- round(predict(a, newdata=me[me$species=="FASY",]),3)
a <- randomForest(rw ~ d + h, data=me[me$species=="PISY",], na.action=na.omit)
me$rwm[me$species=="PISY"] <- round(predict(a, newdata=me[me$species=="PISY",]),3)
me$rw[is.na(me$rw)] <- me$rwm[is.na(me$rw)]
me$rwm <- NULL

me$d2013 <- me$d
me$h2013 <- me$h
me$d2008 <- round(pmax(0, me$d - me$rw*10), 1)
me$h2008 <- round(me$h * me$ih5, 1)
me$d <- NULL
me$h <- NULL
me$rw <- NULL
me$ih5 <- NULL


#There are no plot corners, so they need to be estimated
library(alphahull)
tt <- with(me, unique(data.frame(x,y)))
plot(tt)
tt <- ahull(tt, alpha = 30)
tt <- with(data.frame(tt$ashape.obj$edges), unique(data.frame(x=c(x1, x2), y=c(y1,y2), ind1=c(ind1,ind2), ind2=c(ind2,ind1))))
chain <- c(tt$ind1[1], tt$ind2[1])
while(chain[1] != chain[length(chain)]) {
  for(k in 1:NROW(tt)) {
    if(chain[length(chain)] == tt$ind1[k] && chain[length(chain)-1] != tt$ind2[k]) {
      chain <- c(chain, tt$ind2[k]); break;
    }
  }
}
chain <- data.frame(ind1 = chain[-1])
chain$cnr <- 1:NROW(chain)
tt <- merge(unique(tt[,c("x","y","ind1")]), chain)[,c("x","y","cnr")]
tt <- tt[order(tt$cnr),]
polygon(tt, lwd=2, lty=2)
me2 <- merge(me, tt, all.x=T)
me2$gew <- 1
t1 <- sum(!is.na(me2$cnr))
me2$gew[!is.na(me2$cnr)] <- (t1/2 - 1) / t1

names(me2)[names(me2)=="tree"] <- "#tree"
me2 <- me2[order(me2$"#tree"),]
write.table(me2[,c("#tree", "x", "y", "z", "gew", "species", "d2008", "h2008", "d2013", "h2013")], "/tmp/treesEmf.txt", row.names = FALSE, quote = FALSE)

names(tt)[names(tt)=="cnr"] <- "#corner"
tt <- tt[order(tt$"#corner"),]
write.table(tt[,c("#corner","x","y")], "/tmp/cornersEmf.txt", row.names = FALSE, quote = FALSE)
area <- 0
previous <- NROW(tt)
for(current in 1:NROW(tt)) {
  area <- area + tt$x[previous] * tt$y[current];
  area <- area - tt$y[previous] * tt$x[current];
  previous <- current;
}
(abs(area/2)) #552.4176
