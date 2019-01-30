x <- read.table("/tmp/resHst.txt")
summary(x)
aggregate(V4 ~ V1, FUN=sum)
#Summen ergeben nicht automatisch die Bestandesflaeche -> dort hin gewichten
#Freiflaechen noch nicht beruecksichtigt

