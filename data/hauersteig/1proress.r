#https://doi.org/10.1038/sdata.2018.77
#https://doi.org/10.6084/m9.figshare.c.3971019.v1

dat <- read.csv("https://s3-eu-west-1.amazonaws.com/pstorage-npg-968563215/10596865/corners.csv")
t1 <- with(dat[dat$plot==2 & dat$cor==TRUE,], data.frame(corner, x, y))
names(t1)[1] <- paste0("#", names(t1)[1])
t1 <- t1[order(t1$"#corner"),]
write.table(t1[,c("#corner","x","y")], "/tmp/cornersHst.txt", row.names = FALSE, quote = FALSE)
area <- 0
previous <- NROW(t1)
for(current in 1:NROW(t1)) {
  area <- area + t1$x[previous] * t1$y[current];
  area <- area - t1$y[previous] * t1$x[current];
  previous <- current;
}
(abs(area/2)) #2474.808
rm(dat)

dat <- read.csv("https://s3-eu-west-1.amazonaws.com/pstorage-npg-968563215/10597051/dhcComplSmooth.csv")
datP <- read.csv("https://s3-eu-west-1.amazonaws.com/pstorage-npg-968563215/10596871/pos.csv")

t1 <- with(dat[dat$plot==2 & dat$year %in% c(1989, 1993),], data.frame(tree,year,d,h))
t2 <- reshape(t1, timevar = "year", idvar = "tree", direction="wide", sep="")
me <- merge(datP[datP$plot==2 & datP$core==TRUE, c("tree","x","y","z","species")], t2)
me$gew <- 1
names(me)[1] <- paste0("#", names(me)[1])
me <- me[order(me$"#tree"),]
write.table(me[,c("#tree", "x", "y", "z", "gew", "species", "d1989", "h1989", "d1993", "h1993")], "/tmp/treesHst.txt", row.names = FALSE, quote = FALSE)
