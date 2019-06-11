dxy <- .1    #Grid resolution
xMax <- 40   #Plot size
yMax <- 40

n <- xMax * yMax / 5 #Number of trees
set.seed(0)
td <- data.frame(x = runif(n)*xMax, y = runif(n)*yMax, d = rexp(n, .1), sp=rbinom(n, 1, .4)) #Tree data: x,y..Position, d..diameter, sp..species
td$w <- td$d/4  #Weight

fun <- function(x,y) { #Basic "Kreisbogen" equation
  ff <- function(x,y) {which.min(sqrt((td$x-x)^2 + (td$y-y)^2)/td$w)}
  mapply(ff, x, y)
}

x <- seq(0, xMax, by=dxy)
y <- seq(0, yMax, by=dxy)

#precalculate standarea and match tree size to stand area (as data is random)
res <- outer(x, y, fun)
tt <- table(res)*dxy^2
td$area[as.integer(names(tt))] <- tt
td$area[is.na(td$area)] <- 0
td[, c("x", "y")] <- td[order(td$area), c("x", "y")]
td$area <- NULL

res <- outer(x, y, fun)
png("/tmp/img.png", width = dim(res)[1]*2, height = dim(res)[2]*2, antialias="none", pointsize = 1)
image(x, y, rbind(FALSE, diff(res)!=0) | t(rbind(FALSE, diff(t(res))!=0)), asp=1, axes=F, xlab="", ylab="", col=c(0,1)) #Borders between trees are black
tt <- rbind(FALSE, diff(array(td$sp[res],dim(res)))!=0) | t(rbind(FALSE, diff(t(array(td$sp[res],dim(res))))!=0))
tt[!tt] <- NA
image(x, y, tt, col=c(3), add=T) #Borders between different species are green
with(td, symbols(x=x, y=y, circles=d/100, inches=F, add=T, fg=sp+1, bg=sp+1))
dev.off()

#Calculate stad area
tt <- table(res)*dxy^2
td$area[as.integer(names(tt))] <- tt
td$area[is.na(td$area)] <- 0

