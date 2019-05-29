x <- read.table("/tmp/resHst.txt")
y <- read.table("../../data/treesHst.txt")
#x <- read.table("/tmp/resEmf.txt")
#y <- read.table("../../data/treesEmf.txt")

names(x) <- c("method", "tree", "dens", "area")
x$area <- x$area/10000
y$V7 <- y$V7 + 1.3
y$V9 <- y$V9 + 1.3
ig <- data.frame(tree=y$V1, g=y$V7^2*pi/40000, ig=(y$V9^2 - y$V7^2)*pi/40000/5, d=y$V7)

me <- merge(x[,c("method", "tree", "area")], ig)
me <- me[me$method!="hba",]
me$method <- factor(me$method)
#coplot(ig/area ~ g/area | method, data=me, panel=panel.smooth, xlab="Grundfläche [m²/ha]", ylab="Grundflächenzuwachs [m²/ha/Jahr]")
pdf("/tmp/igg.pdf", width=7, height=6)
coplot(ig/area ~ g/area | method, data=me, panel=panel.smooth, xlab=c("Grundfläche [m²/ha]", NA), ylab="Grundflächenzuwachs [m²/ha/Jahr]", xlim=c(20,100), ylim=c(0,2), show.given = F)
dev.off()
pdf("/tmp/igd.pdf", width=7, height=6)
coplot(ig/area ~ d | method, data=me, panel=panel.smooth, xlab=c("BHD [cm]", NA), ylab="Grundflächenzuwachs [m²/ha/Jahr]", ylim=c(0,2), show.given = F, main="")
dev.off()

iv <- data.frame(tree=y$V1, g=y$V7^2*pi/40000, v=y$V7^2*pi/40000*y$V8*.5, iv=(y$V9^2*y$V10 - y$V7^2*y$V8)*pi/40000/5*.5, d=y$V7)
me <- merge(x[,c("method", "tree", "area")], iv)
me <- me[me$method!="hba",]
me$method <- factor(me$method)
pdf("/tmp/ivg.pdf", width=7, height=6)
coplot(iv/area ~ g/area | method, data=me, panel=panel.smooth, xlab=c("Grundfläche [m²/ha]", NA), ylab="Volumszuwachs [m³/ha/Jahr]", xlim=c(20,100), ylim=c(0,45), show.given = F)
dev.off()
pdf("/tmp/ivd.pdf", width=7, height=6)
coplot(iv/area ~ d | method, data=me, panel=panel.smooth, xlab=c("BHD [cm]", NA), ylab="Volumszuwachs [m³/ha/Jahr]", ylim=c(0,45), show.given = F, main="")
dev.off()

find /tmp/ -name "i[gv][gd].pdf" -exec pdfcrop {} \;
cp /tmp/i[gv][gd]-crop.pdf ~/praesentationen/19/ertragskunde/pic
