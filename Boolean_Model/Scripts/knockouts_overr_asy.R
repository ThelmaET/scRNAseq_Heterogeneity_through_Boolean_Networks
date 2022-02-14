####################
#                  #
#    Knockouts     #
#                  #
####################

getwd()
setwd("C:/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th")
getwd()

library(BoolNet)
library(BoolNetPerturb)
network_mcts <- loadNetwork("network_mcts_c1.txt")


KOAURKA<-fixGenes(network_mcts,"AURKA", 0)
KOAURKA.attr<-getAttractors(KOAURKA, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(KOAURKA.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_KOAURKA.csv")

originalNet <- fixGenes(KOAURKA, "AURKA", -1)



KOBRCA1<-fixGenes(network_mcts,"BRCA1", 0)
KOBRCA1.attr<-getAttractors(KOBRCA1, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(KOBRCA1.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_KOBRCA1.csv")

originalNet <- fixGenes(KOBRCA1, "BRCA1", -1)



KOBIRC5<-fixGenes(network_mcts,"BIRC5", 0)
KOBIRC5.attr<-getAttractors(KOBIRC5, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(KOBIRC5.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_KOBIRC5.csv")

originalNet <- fixGenes(KOBIRC5, "BIRC5", -1)



KOBTG2<-fixGenes(network_mcts,"BTG2", 0)
KOBTG2.attr<-getAttractors(KOBTG2, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(KOBTG2.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_KOBTG2.csv")

originalNet <- fixGenes(KOBTG2, "BTG2", -1)



KOCCNB1<-fixGenes(network_mcts,"CCNB1", 0)
KOCCNB1.attr<-getAttractors(KOCCNB1, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(KOCCNB1.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_KOCCNB1.csv")

originalNet <- fixGenes(KOCCNB1, "CCNB1", -1)



KOCCND1<-fixGenes(network_mcts,"CCND1", 0)
KOCCND1.attr<-getAttractors(KOCCND1, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(KOCCND1.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_KOCCND1.csv")

originalNet <- fixGenes(KOCCND1, "CCND1", -1)



KOCREB1<-fixGenes(network_mcts,"CREB1", 0)
KOCREB1.attr<-getAttractors(KOCREB1, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(KOCREB1.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_KOCREB1.csv")

originalNet <- fixGenes(KOCREB1, "CREB1", -1)



KOESR1<-fixGenes(network_mcts,"ESR1", 0)
KOESR1.attr<-getAttractors(KOESR1, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(KOESR1.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_KOESR1.csv")

originalNet <- fixGenes(KOESR1, "ESR1", -1)



KOFOXA1<-fixGenes(network_mcts,"FOXA1", 0)
KOFOXA1.attr<-getAttractors(KOFOXA1, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(KOFOXA1.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_KOFOXA1.csv")

originalNet <- fixGenes(KOFOXA1, "FOXA1", -1)



KOFOXM1<-fixGenes(network_mcts,"FOXM1", 0)
KOFOXM1.attr<-getAttractors(KOFOXM1, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(KOFOXM1.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_KOFOXM1.csv")

originalNet <- fixGenes(KOFOXM1, "FOXM1", -1)



KOGATA3<-fixGenes(network_mcts,"GATA3", 0)
KOGATA3.attr<-getAttractors(KOGATA3, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(KOGATA3.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_KOGATA3.csv")

originalNet <- fixGenes(KOGATA3, "GATA3", -1)



KOIL20<-fixGenes(network_mcts,"IL20", 0)
KOIL20.attr<-getAttractors(KOIL20, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(KOIL20.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_KOIL20.csv")

originalNet <- fixGenes(KOIL20, "IL20", -1)



KOMKI67<-fixGenes(network_mcts,"MKI67", 0)
KOMKI67.attr<-getAttractors(KOMKI67, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(KOMKI67.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_KOMKI67.csv")

originalNet <- fixGenes(KOMKI67, "MKI67", -1)



KOMYC<-fixGenes(network_mcts,"MYC", 0)
KOMYC.attr<-getAttractors(KOMYC, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(KOMYC.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_KOMYC.csv")

originalNet <- fixGenes(KOMYC, "MYC", -1)



KONFE2L2<-fixGenes(network_mcts,"NFE2L2", 0)
KONFE2L2.attr<-getAttractors(KONFE2L2, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(KONFE2L2.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_KONFE2L2.csv")

originalNet <- fixGenes(KONFE2L2, "NFE2L2", -1)



KOPTTG1<-fixGenes(network_mcts,"PTTG1", 0)
KOPTTG1.attr<-getAttractors(KOPTTG1, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(KOPTTG1.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_KOPTTG1.csv")

originalNet <- fixGenes(KOPTTG1, "PTTG1", -1)



KOSTAT1<-fixGenes(network_mcts,"STAT1", 0)
KOSTAT1.attr<-getAttractors(KOSTAT1, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(KOSTAT1.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_KOSTAT1.csv")

originalNet <- fixGenes(KOSTAT1, "STAT1", -1)



KOSTAT3<-fixGenes(network_mcts,"STAT3", 0)
KOSTAT3.attr<-getAttractors(KOSTAT3, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(KOSTAT3.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_KOSTAT3.csv")

originalNet <- fixGenes(KOSTAT3, "STAT3", -1)



KOSNAIL<-fixGenes(network_mcts,"SNAIL", 0)
KOSNAIL.attr<-getAttractors(KOSNAIL, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(KOSNAIL.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_KOSNAIL.csv")

originalNet <- fixGenes(KOSNAIL, "SNAIL", -1)



KOSOD2<-fixGenes(network_mcts,"SOD2", 0)
KOSOD2.attr<-getAttractors(KOSOD2, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(KOSOD2.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_KOSOD2.csv")

originalNet <- fixGenes(KOSOD2, "SOD2", -1)



KOS1007<-fixGenes(network_mcts,"S1007", 0)
KOS1007.attr<-getAttractors(KOS1007, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(KOS1007.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_KOS1007.csv")

originalNet <- fixGenes(KOS1007, "S1007", -1)



KOTP53<-fixGenes(network_mcts,"TP53", 0)
KOTP53.attr<-getAttractors(KOTP53, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(KOTP53.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_KOTP53.csv")

originalNet <- fixGenes(KOTP53, "TP53", -1)



KOTWIST1<-fixGenes(network_mcts,"TWIST1", 0)
KOTWIST1.attr<-getAttractors(KOTWIST1, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(KOTWIST1.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_KOTWIST1.csv")

originalNet <- fixGenes(KOTWIST1, "TWIST1", -1)






#####################
#                   #
#  Overexpression   #
#                   #
#####################

getwd()
setwd("C:/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th")
getwd()

library(BoolNet)
library(BoolNetPerturb)
network_mcts <- loadNetwork("network_mcts_c1.txt")


OVEAURKA<-fixGenes(network_mcts,"AURKA", 1)
OVEAURKA.attr<-getAttractors(OVEAURKA, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVEAURKA.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVEAURKA.csv")

originalNet <- fixGenes(OVEAURKA, "AURKA", -1)



OVEBRCA1<-fixGenes(network_mcts,"BRCA1", 1)
OVEBRCA1.attr<-getAttractors(OVEBRCA1, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVEBRCA1.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVEBRCA1.csv")

originalNet <- fixGenes(OVEBRCA1, "BRCA1", -1)



OVEBIRC5<-fixGenes(network_mcts,"BIRC5", 1)
OVEBIRC5.attr<-getAttractors(OVEBIRC5, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVEBIRC5.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVEBIRC5.csv")

originalNet <- fixGenes(OVEBIRC5, "BIRC5", -1)



OVEBTG2<-fixGenes(network_mcts,"BTG2", 1)
OVEBTG2.attr<-getAttractors(OVEBTG2, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVEBTG2.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVEBTG2.csv")

originalNet <- fixGenes(OVEBTG2, "BTG2", -1)



OVECCNB1<-fixGenes(network_mcts,"CCNB1", 1)
OVECCNB1.attr<-getAttractors(OVECCNB1, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVECCNB1.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVECCNB1.csv")

originalNet <- fixGenes(OVECCNB1, "CCNB1", -1)



OVECCND1<-fixGenes(network_mcts,"CCND1", 1)
OVECCND1.attr<-getAttractors(OVECCND1, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVECCND1.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVECCND1.csv")

originalNet <- fixGenes(OVECCND1, "CCND1", -1)



OVECREB1<-fixGenes(network_mcts,"CREB1", 1)
OVECREB1.attr<-getAttractors(OVECREB1, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVECREB1.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVECREB1.csv")

originalNet <- fixGenes(OVECREB1, "CREB1", -1)



OVEESR1<-fixGenes(network_mcts,"ESR1", 1)
OVEESR1.attr<-getAttractors(OVEESR1, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVEESR1.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVEESR1.csv")

originalNet <- fixGenes(OVEESR1, "ESR1", -1)



OVEFOXA1<-fixGenes(network_mcts,"FOXA1", 1)
OVEFOXA1.attr<-getAttractors(OVEFOXA1, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVEFOXA1.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVEFOXA1.csv")

originalNet <- fixGenes(OVEFOXA1, "FOXA1", -1)



OVEFOXM1<-fixGenes(network_mcts,"FOXM1", 1)
OVEFOXM1.attr<-getAttractors(OVEFOXM1, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVEFOXM1.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVEFOXM1.csv")

originalNet <- fixGenes(OVEFOXM1, "FOXM1", -1)



OVEGATA3<-fixGenes(network_mcts,"GATA3", 1)
OVEGATA3.attr<-getAttractors(OVEGATA3, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVEGATA3.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVEGATA3.csv")

originalNet <- fixGenes(OVEGATA3, "GATA3", -1)



OVEIL20<-fixGenes(network_mcts,"IL20", 1)
OVEIL20.attr<-getAttractors(OVEIL20, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVEIL20.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVEIL20.csv")

originalNet <- fixGenes(OVEIL20, "IL20", -1)



OVEMKI67<-fixGenes(network_mcts,"MKI67", 1)
OVEMKI67.attr<-getAttractors(OVEMKI67, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVEMKI67.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVEMKI67.csv")

originalNet <- fixGenes(OVEMKI67, "MKI67", -1)



OVEMYC<-fixGenes(network_mcts,"MYC", 1)
OVEMYC.attr<-getAttractors(OVEMYC, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVEMYC.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVEMYC.csv")

originalNet <- fixGenes(OVEMYC, "MYC", -1)



OVENFE2L2<-fixGenes(network_mcts,"NFE2L2", 1)
OVENFE2L2.attr<-getAttractors(OVENFE2L2, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVENFE2L2.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVENFE2L2.csv")

originalNet <- fixGenes(OVENFE2L2, "NFE2L2", -1)



OVEPTTG1<-fixGenes(network_mcts,"PTTG1", 1)
OVEPTTG1.attr<-getAttractors(OVEPTTG1, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVEPTTG1.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVEPTTG1.csv")

originalNet <- fixGenes(OVEPTTG1, "PTTG1", -1)



OVESTAT1<-fixGenes(network_mcts,"STAT1", 1)
OVESTAT1.attr<-getAttractors(OVESTAT1, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVESTAT1.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVESTAT1.csv")

originalNet <- fixGenes(OVESTAT1, "STAT1", -1)



OVESTAT3<-fixGenes(network_mcts,"STAT3", 1)
OVESTAT3.attr<-getAttractors(OVESTAT3, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVESTAT3.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVESTAT3.csv")

originalNet <- fixGenes(OVESTAT3, "STAT3", -1)



OVESNAIL<-fixGenes(network_mcts,"SNAIL", 1)
OVESNAIL.attr<-getAttractors(OVESNAIL, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVESNAIL.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVESNAIL.csv")

originalNet <- fixGenes(OVESNAIL, "SNAIL", -1)



OVESOD2<-fixGenes(network_mcts,"SOD2", 1)
OVESOD2.attr<-getAttractors(OVESOD2, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVESOD2.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVESOD2.csv")

originalNet <- fixGenes(OVESOD2, "SOD2", -1)



OVES1007<-fixGenes(network_mcts,"S1007", 1)
OVES1007.attr<-getAttractors(OVES1007, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVES1007.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVES1007.csv")

originalNet <- fixGenes(OVES1007, "S1007", -1)



OVETP53<-fixGenes(network_mcts,"TP53", 1)
OVETP53.attr<-getAttractors(OVETP53, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVETP53.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVETP53.csv")

originalNet <- fixGenes(OVETP53, "TP53", -1)



OVETWIST1<-fixGenes(network_mcts,"TWIST1", 1)
OVETWIST1.attr<-getAttractors(OVETWIST1, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVETWIST1.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVETWIST1.csv")

originalNet <- fixGenes(OVETWIST1, "TWIST1", -1)





####################################
#                                  #
#  Simultaneous perturbations      #
#                                  #
####################################

OVEBRCA1_TP53 <- fixGenes(network_mcts, c("BRCA1","TP53"), c(1,1))
OVEBRCA1_TP53.attr<-getAttractors(OVEBRCA1_TP53, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVEBRCA1_TP53.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVEBRCA1_TP53.csv")

originalNet <- fixGenes(OVEBRCA1_TP53, c("BRCA1","TP53"), -1)



OVEBRCA1_CCNB1 <- fixGenes(network_mcts, c("BRCA1","CCNB1"), c(1,1))
OVEBRCA1_CCNB1.attr<-getAttractors(OVEBRCA1_CCNB1, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVEBRCA1_CCNB1.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVEBRCA1_CCNB1.csv")

originalNet <- fixGenes(OVEBRCA1_CCNB1, c("BRCA1","CCNB1"), -1)



OVEBRCA1_CCNB1_TP53 <- fixGenes(network_mcts, c("BRCA1","CCNB1", "TP53"), c(1,1,1))
OVEBRCA1_CCNB1_TP53.attr<-getAttractors(OVEBRCA1_CCNB1_TP53, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVEBRCA1_CCNB1_TP53.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVEBRCA1_CCNB1_TP53.csv")

originalNet <- fixGenes(OVEBRCA1_CCNB1_TP53, c("BRCA1","CCNB1", "TP53"), -1)



OVETP53_KOTWIST1 <- fixGenes(network_mcts, c("TP53","TWIST1"), c(1,0))
OVETP53_KOTWIST1.attr<-getAttractors(OVETP53_KOTWIST1, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVETP53_KOTWIST1.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVETP53_KOTWIST1.csv")

originalNet <- fixGenes(OVETP53_KOTWIST1, c("TP53","TWIST1"), -1)



OVETP53_KOS1007 <- fixGenes(network_mcts, c("TP53","S1007"), c(1,0))
OVETP53_KOS1007.attr<-getAttractors(OVETP53_KOS1007, method = "random", type="asynchronous", startStates = 10000)
attr.df <- attractorToDataframe(OVETP53_KOS1007.attr, Boolean = TRUE)
head(attr.df)
write.csv(attr.df, "/Users/thelm/OneDrive/Documentos/maestria/proyecto/boolnet/red_th/perturbaciones_asy/atractores_asy_c1_OVETP53_KOS1007.csv")

originalNet <- fixGenes(OVETP53_KOS1007, c("TP53","S1007"), -1)
