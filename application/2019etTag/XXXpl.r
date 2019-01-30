x <- read.table("../../data/treesEmf.txt", comment.char = "", header =T)
y <- read.table("../../data/cornersEmf.txt", comment.char = "", header =T)

with(x, plot(x, y, asp=1))
with(y, lines(x,y))

symbols(x=-3.76, y=10.5, circles=4, inches=F, add=T)
symbols(x=.04, y=12.4, circles=7, inches=F, add=T)
symbols(x=1.94, y=13.5, circles=4, inches=F, add=T)

with(x[x$X.tree==68,], points(x,y,col=2))
points(x=.04, y=12.4, col=3)
