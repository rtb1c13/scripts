# R script to run anovas and print Tukey HSD results
# Used in conjunction with reorder_results.py
# Reads in data frames in format required by aov

for(i in 1:47)
{
  ligand = read.table(paste("lig",i,".txt",sep=""))
  results=aov(ligand$V1 ~ ligand$V2,data=ligand)
  Tukeyresults = TukeyHSD(results,conf.level=0.95)
  sink(paste("TukeyHSD_",i,".txt",sep=""),append=FALSE,split=FALSE)
  print(Tukeyresults)
  sink()
}
