# Adapted from Zohra Ouaray, 2015
#Read in structure data - this should be list of coordinates (one line)

data<-scan("ALL_6G08_MD.txt")

#number of coordinates per snapshot (ie. number of atoms *3)
coords=1317
#number of snapshots
n=1000
# Transform to correct array and do PCA
data_matrix<-array(data, dim=c(coords,n))
matrix_trans<-t(data_matrix)
library(bio3d)
pcadata<-pca.xyz(matrix_trans)

# Print animations of principal components
mktrj.pca(pcadata, file="PDB-files/AllStruct_PC1_MD.pdb", pc=1, mag=1,step=0.050)
mktrj.pca(pcadata, file="PDB-files/AllStruct_PC2_MD.pdb", pc=2, mag=1,step=0.050)
mktrj.pca(pcadata, file="PDB-files/AllStruct_PC3_MD.pdb", pc=3, mag=1,step=0.050)

# First plot histograms of full structure sets, 1 -> n
PC1<-pcadata$z[1:n]
PC2<-pcadata$z[(n+1):(n*2)]
PC3<-pcadata$z[((n*2)+1):(n*3)]
pdf("PCA_histograms.pdf")
par(mfrow=c(2,2))
hist(PC1,breaks=seq(-200,200,10),ylim=c(0,200))
hist(PC2,breaks=seq(-200,200,10),ylim=c(0,200))
hist(PC3,breaks=seq(-200,200,10),ylim=c(0,200))

#Now plot histogram of total overlaid with smoothed density plots of each of the colours
pdf("PC1_densities.pdf")
hist(PC1,breaks=seq(-200,200,10),probability=TRUE,ylim=c(0,0.025),col="grey")
saxsfits<-scan("saxs_fits.txt") # Reads a column of saxs fit values
saxsfits_bool <- saxsfits<3.00  #Black
PC1_black <- PC1[saxsfits_bool]
lines(density(PC1_black,adjust=0.5),col="black",lwd=2)
saxsfits_bool <- saxsfits>=3.00 & saxsfits<6.00  #Blue
PC1_blue <- PC1[saxsfits_bool]
lines(density(PC1_blue,adjust=0.5),col="blue",lwd=2)
saxsfits_bool <- saxsfits>=6.00 & saxsfits<9.00  #Green
PC1_green <- PC1[saxsfits_bool]
lines(density(PC1_green,adjust=0.5),col="green",lwd=2)
saxsfits_bool <- saxsfits>=9.00  #Red
PC1_red <- PC1[saxsfits_bool]
lines(density(PC1_red,adjust=0.5),col="red",lwd=2)
q() # QUITS HERE CURRENTLY!

#plot a graph of the pcadata
pdf("ALL_6G08_MD_ONLY_10.pdf",width=14,height=14)
plot.pca(pcadata,col=terrain.colors(n*1.5))

PC1<-pcadata$z[1:n]
PC2<-pcadata$z[(n+1):(n*2)]
PC3<-pcadata$z[((n*2)+1):(n*3)]

#barplot(PC1,xlim=c(0,191),ylim=c(-100,100), main = "Residues contribution to PC1",xlab = "residue number", ylab = "PC1")
#barplot(PC2,xlim=c(0,191),ylim=c(-100,100), main = "Residues contribution to PC2",xlab = "residue number", ylab = "PC2")
#barplot(PC3,xlim=c(0,191),ylim=c(-100,100), main = "Residues contribution to PC3",xlab = "residue number", ylab = "PC3")

#pc1<-scan("PC1_eigvalue.out")
#h1<-hist(PC1)
#pc2<-scan("PC2_eigvalue.out")
#h2<-hist(PC2)
#plot PC1 vs PC2 for trajectory of model2
#plot(PC1,PC2,xlim=c(-10,25),ylim=c(-30,10), add=TRUE)
#plot(h1, col="black")
#plot(h2, col="black")


#allow to add projection on the same plot:
#par(new=TRUE)
plot(PC1,PC2,xlim=c(-150,200),ylim=c(-40,50))

xCA_6G08_fab_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/CA_6G08_fab.txt")
xCA_6G08_fab_proj<-project.pca(xCA_6G08_fab_raw,pcadata)
points(xCA_6G08_fab_proj[1],xCA_6G08_fab_proj[2],pch=20)
text(xCA_6G08_fab_proj[1],xCA_6G08_fab_proj[2],pos=4,label="6G08",col="purple")

x6G08_MD_1_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_1.txt")
x6G08_MD_1_proj<-project.pca(x6G08_MD_1_raw,pcadata)
points(x6G08_MD_1_proj[1],x6G08_MD_1_proj[2],pch=20)
text(x6G08_MD_1_proj[1],x6G08_MD_1_proj[2],pos=4,label="1",col="blue")

x6G08_MD_2_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_2.txt")
x6G08_MD_2_proj<-project.pca(x6G08_MD_2_raw,pcadata)
points(x6G08_MD_2_proj[1],x6G08_MD_2_proj[2],pch=20)
text(x6G08_MD_2_proj[1],x6G08_MD_2_proj[2],pos=4,label="2",col="blue")

x6G08_MD_3_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_3.txt")
x6G08_MD_3_proj<-project.pca(x6G08_MD_3_raw,pcadata)
points(x6G08_MD_3_proj[1],x6G08_MD_3_proj[2],pch=20)
text(x6G08_MD_3_proj[1],x6G08_MD_3_proj[2],pos=4,label="3",col="blue")

x6G08_MD_4_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_4.txt")
x6G08_MD_4_proj<-project.pca(x6G08_MD_4_raw,pcadata)
points(x6G08_MD_4_proj[1],x6G08_MD_4_proj[2],pch=20)
text(x6G08_MD_4_proj[1],x6G08_MD_4_proj[2],pos=4,label="4",col="blue")

x6G08_MD_5_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_5.txt")
x6G08_MD_5_proj<-project.pca(x6G08_MD_5_raw,pcadata)
points(x6G08_MD_5_proj[1],x6G08_MD_5_proj[2],pch=20)
text(x6G08_MD_5_proj[1],x6G08_MD_5_proj[2],pos=4,label="5",col="blue")

x6G08_MD_6_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_6.txt")
x6G08_MD_6_proj<-project.pca(x6G08_MD_6_raw,pcadata)
points(x6G08_MD_6_proj[1],x6G08_MD_6_proj[2],pch=20)
text(x6G08_MD_6_proj[1],x6G08_MD_6_proj[2],pos=4,label="6",col="black")

x6G08_MD_7_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_7.txt")
x6G08_MD_7_proj<-project.pca(x6G08_MD_7_raw,pcadata)
points(x6G08_MD_7_proj[1],x6G08_MD_7_proj[2],pch=20)
text(x6G08_MD_7_proj[1],x6G08_MD_7_proj[2],pos=4,label="7",col="blue")

x6G08_MD_8_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_8.txt")
x6G08_MD_8_proj<-project.pca(x6G08_MD_8_raw,pcadata)
points(x6G08_MD_8_proj[1],x6G08_MD_8_proj[2],pch=20)
text(x6G08_MD_8_proj[1],x6G08_MD_8_proj[2],pos=4,label="8",col="blue")

x6G08_MD_9_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_9.txt")
x6G08_MD_9_proj<-project.pca(x6G08_MD_9_raw,pcadata)
points(x6G08_MD_9_proj[1],x6G08_MD_9_proj[2],pch=20)
text(x6G08_MD_9_proj[1],x6G08_MD_9_proj[2],pos=4,label="9",col="green")

x6G08_MD_10_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_10.txt")
x6G08_MD_10_proj<-project.pca(x6G08_MD_10_raw,pcadata)
points(x6G08_MD_10_proj[1],x6G08_MD_10_proj[2],pch=20)
text(x6G08_MD_10_proj[1],x6G08_MD_10_proj[2],pos=4,label="10",col="red")

x6G08_MD_11_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_11.txt")
x6G08_MD_11_proj<-project.pca(x6G08_MD_11_raw,pcadata)
points(x6G08_MD_11_proj[1],x6G08_MD_11_proj[2],pch=20)
text(x6G08_MD_11_proj[1],x6G08_MD_11_proj[2],pos=4,label="11",col="blue")

x6G08_MD_12_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_12.txt")
x6G08_MD_12_proj<-project.pca(x6G08_MD_12_raw,pcadata)
points(x6G08_MD_12_proj[1],x6G08_MD_12_proj[2],pch=20)
text(x6G08_MD_12_proj[1],x6G08_MD_12_proj[2],pos=4,label="12",col="black")

x6G08_MD_13_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_13.txt")
x6G08_MD_13_proj<-project.pca(x6G08_MD_13_raw,pcadata)
points(x6G08_MD_13_proj[1],x6G08_MD_13_proj[2],pch=20)
text(x6G08_MD_13_proj[1],x6G08_MD_13_proj[2],pos=4,label="13",col="black")

x6G08_MD_14_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_14.txt")
x6G08_MD_14_proj<-project.pca(x6G08_MD_14_raw,pcadata)
points(x6G08_MD_14_proj[1],x6G08_MD_14_proj[2],pch=20)
text(x6G08_MD_14_proj[1],x6G08_MD_14_proj[2],pos=4,label="14",col="blue")

x6G08_MD_15_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_15.txt")
x6G08_MD_15_proj<-project.pca(x6G08_MD_15_raw,pcadata)
points(x6G08_MD_15_proj[1],x6G08_MD_15_proj[2],pch=20)
text(x6G08_MD_15_proj[1],x6G08_MD_15_proj[2],pos=4,label="15",col="black")

x6G08_MD_16_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_16.txt")
x6G08_MD_16_proj<-project.pca(x6G08_MD_16_raw,pcadata)
points(x6G08_MD_16_proj[1],x6G08_MD_16_proj[2],pch=20)
text(x6G08_MD_16_proj[1],x6G08_MD_16_proj[2],pos=4,label="16",col="black")

x6G08_MD_17_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_17.txt")
x6G08_MD_17_proj<-project.pca(x6G08_MD_17_raw,pcadata)
points(x6G08_MD_17_proj[1],x6G08_MD_17_proj[2],pch=20)
text(x6G08_MD_17_proj[1],x6G08_MD_17_proj[2],pos=4,label="17",col="blue")

x6G08_MD_18_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_18.txt")
x6G08_MD_18_proj<-project.pca(x6G08_MD_18_raw,pcadata)
points(x6G08_MD_18_proj[1],x6G08_MD_18_proj[2],pch=20)
text(x6G08_MD_18_proj[1],x6G08_MD_18_proj[2],pos=4,label="18",col="black")

x6G08_MD_19_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_19.txt")
x6G08_MD_19_proj<-project.pca(x6G08_MD_19_raw,pcadata)
points(x6G08_MD_19_proj[1],x6G08_MD_19_proj[2],pch=20)
text(x6G08_MD_19_proj[1],x6G08_MD_19_proj[2],pos=4,label="19",col="black")

x6G08_MD_20_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_20.txt")
x6G08_MD_20_proj<-project.pca(x6G08_MD_20_raw,pcadata)
points(x6G08_MD_20_proj[1],x6G08_MD_20_proj[2],pch=20)
text(x6G08_MD_20_proj[1],x6G08_MD_20_proj[2],pos=4,label="20",col="black")

x6G08_MD_21_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_21.txt")
x6G08_MD_21_proj<-project.pca(x6G08_MD_21_raw,pcadata)
points(x6G08_MD_21_proj[1],x6G08_MD_21_proj[2],pch=20)
text(x6G08_MD_21_proj[1],x6G08_MD_21_proj[2],pos=4,label="21",col="blue")

x6G08_MD_22_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_22.txt")
x6G08_MD_22_proj<-project.pca(x6G08_MD_22_raw,pcadata)
points(x6G08_MD_22_proj[1],x6G08_MD_22_proj[2],pch=20)
text(x6G08_MD_22_proj[1],x6G08_MD_22_proj[2],pos=4,label="22",col="blue")

x6G08_MD_23_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_23.txt")
x6G08_MD_23_proj<-project.pca(x6G08_MD_23_raw,pcadata)
points(x6G08_MD_23_proj[1],x6G08_MD_23_proj[2],pch=20)
text(x6G08_MD_23_proj[1],x6G08_MD_23_proj[2],pos=4,label="23",col="blue")

x6G08_MD_24_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_24.txt")
x6G08_MD_24_proj<-project.pca(x6G08_MD_24_raw,pcadata)
points(x6G08_MD_24_proj[1],x6G08_MD_24_proj[2],pch=20)
text(x6G08_MD_24_proj[1],x6G08_MD_24_proj[2],pos=4,label="24",col="blue")

x6G08_MD_25_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_25.txt")
x6G08_MD_25_proj<-project.pca(x6G08_MD_25_raw,pcadata)
points(x6G08_MD_25_proj[1],x6G08_MD_25_proj[2],pch=20)
text(x6G08_MD_25_proj[1],x6G08_MD_25_proj[2],pos=4,label="25",col="blue")

x6G08_MD_26_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_26.txt")
x6G08_MD_26_proj<-project.pca(x6G08_MD_26_raw,pcadata)
points(x6G08_MD_26_proj[1],x6G08_MD_26_proj[2],pch=20)
text(x6G08_MD_26_proj[1],x6G08_MD_26_proj[2],pos=4,label="26",col="blue")

x6G08_MD_27_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_27.txt")
x6G08_MD_27_proj<-project.pca(x6G08_MD_27_raw,pcadata)
points(x6G08_MD_27_proj[1],x6G08_MD_27_proj[2],pch=20)
text(x6G08_MD_27_proj[1],x6G08_MD_27_proj[2],pos=4,label="27",col="blue")

x6G08_MD_28_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_28.txt")
x6G08_MD_28_proj<-project.pca(x6G08_MD_28_raw,pcadata)
points(x6G08_MD_28_proj[1],x6G08_MD_28_proj[2],pch=20)
text(x6G08_MD_28_proj[1],x6G08_MD_28_proj[2],pos=4,label="28",col="blue")

x6G08_MD_29_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_29.txt")
x6G08_MD_29_proj<-project.pca(x6G08_MD_29_raw,pcadata)
points(x6G08_MD_29_proj[1],x6G08_MD_29_proj[2],pch=20)
text(x6G08_MD_29_proj[1],x6G08_MD_29_proj[2],pos=4,label="29",col="blue")

x6G08_MD_30_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_30.txt")
x6G08_MD_30_proj<-project.pca(x6G08_MD_30_raw,pcadata)
points(x6G08_MD_30_proj[1],x6G08_MD_30_proj[2],pch=20)
text(x6G08_MD_30_proj[1],x6G08_MD_30_proj[2],pos=4,label="30",col="black")

x6G08_MD_31_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_31.txt")
x6G08_MD_31_proj<-project.pca(x6G08_MD_31_raw,pcadata)
points(x6G08_MD_31_proj[1],x6G08_MD_31_proj[2],pch=20)
text(x6G08_MD_31_proj[1],x6G08_MD_31_proj[2],pos=4,label="31",col="black")

x6G08_MD_32_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_32.txt")
x6G08_MD_32_proj<-project.pca(x6G08_MD_32_raw,pcadata)
points(x6G08_MD_32_proj[1],x6G08_MD_32_proj[2],pch=20)
text(x6G08_MD_32_proj[1],x6G08_MD_32_proj[2],pos=4,label="32",col="black")

x6G08_MD_33_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_33.txt")
x6G08_MD_33_proj<-project.pca(x6G08_MD_33_raw,pcadata)
points(x6G08_MD_33_proj[1],x6G08_MD_33_proj[2],pch=20)
text(x6G08_MD_33_proj[1],x6G08_MD_33_proj[2],pos=4,label="33",col="black")

x6G08_MD_34_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_34.txt")
x6G08_MD_34_proj<-project.pca(x6G08_MD_34_raw,pcadata)
points(x6G08_MD_34_proj[1],x6G08_MD_34_proj[2],pch=20)
text(x6G08_MD_34_proj[1],x6G08_MD_34_proj[2],pos=4,label="34",col="blue")

x6G08_MD_35_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_35.txt")
x6G08_MD_35_proj<-project.pca(x6G08_MD_35_raw,pcadata)
points(x6G08_MD_35_proj[1],x6G08_MD_35_proj[2],pch=20)
text(x6G08_MD_35_proj[1],x6G08_MD_35_proj[2],pos=4,label="35",col="black")

x6G08_MD_36_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_36.txt")
x6G08_MD_36_proj<-project.pca(x6G08_MD_36_raw,pcadata)
points(x6G08_MD_36_proj[1],x6G08_MD_36_proj[2],pch=20)
text(x6G08_MD_36_proj[1],x6G08_MD_36_proj[2],pos=4,label="36",col="black")

x6G08_MD_37_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_37.txt")
x6G08_MD_37_proj<-project.pca(x6G08_MD_37_raw,pcadata)
points(x6G08_MD_37_proj[1],x6G08_MD_37_proj[2],pch=20)
text(x6G08_MD_37_proj[1],x6G08_MD_37_proj[2],pos=4,label="37",col="blue")

x6G08_MD_38_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_38.txt")
x6G08_MD_38_proj<-project.pca(x6G08_MD_38_raw,pcadata)
points(x6G08_MD_38_proj[1],x6G08_MD_38_proj[2],pch=20)
text(x6G08_MD_38_proj[1],x6G08_MD_38_proj[2],pos=4,label="38",col="blue")

x6G08_MD_39_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_39.txt")
x6G08_MD_39_proj<-project.pca(x6G08_MD_39_raw,pcadata)
points(x6G08_MD_39_proj[1],x6G08_MD_39_proj[2],pch=20)
text(x6G08_MD_39_proj[1],x6G08_MD_39_proj[2],pos=4,label="39",col="blue")

x6G08_MD_40_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_40.txt")
x6G08_MD_40_proj<-project.pca(x6G08_MD_40_raw,pcadata)
points(x6G08_MD_40_proj[1],x6G08_MD_40_proj[2],pch=20)
text(x6G08_MD_40_proj[1],x6G08_MD_40_proj[2],pos=4,label="40",col="black")

x6G08_MD_41_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_41.txt")
x6G08_MD_41_proj<-project.pca(x6G08_MD_41_raw,pcadata)
points(x6G08_MD_41_proj[1],x6G08_MD_41_proj[2],pch=20)
text(x6G08_MD_41_proj[1],x6G08_MD_41_proj[2],pos=4,label="41",col="black")

x6G08_MD_42_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_42.txt")
x6G08_MD_42_proj<-project.pca(x6G08_MD_42_raw,pcadata)
points(x6G08_MD_42_proj[1],x6G08_MD_42_proj[2],pch=20)
text(x6G08_MD_42_proj[1],x6G08_MD_42_proj[2],pos=4,label="42",col="black")

x6G08_MD_43_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_43.txt")
x6G08_MD_43_proj<-project.pca(x6G08_MD_43_raw,pcadata)
points(x6G08_MD_43_proj[1],x6G08_MD_43_proj[2],pch=20)
text(x6G08_MD_43_proj[1],x6G08_MD_43_proj[2],pos=4,label="43",col="black")

x6G08_MD_44_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_44.txt")
x6G08_MD_44_proj<-project.pca(x6G08_MD_44_raw,pcadata)
points(x6G08_MD_44_proj[1],x6G08_MD_44_proj[2],pch=20)
text(x6G08_MD_44_proj[1],x6G08_MD_44_proj[2],pos=4,label="44",col="black")

x6G08_MD_45_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_45.txt")
x6G08_MD_45_proj<-project.pca(x6G08_MD_45_raw,pcadata)
points(x6G08_MD_45_proj[1],x6G08_MD_45_proj[2],pch=20)
text(x6G08_MD_45_proj[1],x6G08_MD_45_proj[2],pos=4,label="45",col="black")

x6G08_MD_46_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_46.txt")
x6G08_MD_46_proj<-project.pca(x6G08_MD_46_raw,pcadata)
points(x6G08_MD_46_proj[1],x6G08_MD_46_proj[2],pch=20)
text(x6G08_MD_46_proj[1],x6G08_MD_46_proj[2],pos=4,label="46",col="black")

x6G08_MD_47_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_47.txt")
x6G08_MD_47_proj<-project.pca(x6G08_MD_47_raw,pcadata)
points(x6G08_MD_47_proj[1],x6G08_MD_47_proj[2],pch=20)
text(x6G08_MD_47_proj[1],x6G08_MD_47_proj[2],pos=4,label="47",col="black")

x6G08_MD_48_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_48.txt")
x6G08_MD_48_proj<-project.pca(x6G08_MD_48_raw,pcadata)
points(x6G08_MD_48_proj[1],x6G08_MD_48_proj[2],pch=20)
text(x6G08_MD_48_proj[1],x6G08_MD_48_proj[2],pos=4,label="48",col="black")

x6G08_MD_49_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_49.txt")
x6G08_MD_49_proj<-project.pca(x6G08_MD_49_raw,pcadata)
points(x6G08_MD_49_proj[1],x6G08_MD_49_proj[2],pch=20)
text(x6G08_MD_49_proj[1],x6G08_MD_49_proj[2],pos=4,label="49",col="black")

x6G08_MD_50_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_50.txt")
x6G08_MD_50_proj<-project.pca(x6G08_MD_50_raw,pcadata)
points(x6G08_MD_50_proj[1],x6G08_MD_50_proj[2],pch=20)
text(x6G08_MD_50_proj[1],x6G08_MD_50_proj[2],pos=4,label="50",col="black")

x6G08_MD_51_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_51.txt")
x6G08_MD_51_proj<-project.pca(x6G08_MD_51_raw,pcadata)
points(x6G08_MD_51_proj[1],x6G08_MD_51_proj[2],pch=20)
text(x6G08_MD_51_proj[1],x6G08_MD_51_proj[2],pos=4,label="51",col="black")

x6G08_MD_52_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_52.txt")
x6G08_MD_52_proj<-project.pca(x6G08_MD_52_raw,pcadata)
points(x6G08_MD_52_proj[1],x6G08_MD_52_proj[2],pch=20)
text(x6G08_MD_52_proj[1],x6G08_MD_52_proj[2],pos=4,label="52",col="black")

x6G08_MD_53_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_53.txt")
x6G08_MD_53_proj<-project.pca(x6G08_MD_53_raw,pcadata)
points(x6G08_MD_53_proj[1],x6G08_MD_53_proj[2],pch=20)
text(x6G08_MD_53_proj[1],x6G08_MD_53_proj[2],pos=4,label="53",col="black")

x6G08_MD_54_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_54.txt")
x6G08_MD_54_proj<-project.pca(x6G08_MD_54_raw,pcadata)
points(x6G08_MD_54_proj[1],x6G08_MD_54_proj[2],pch=20)
text(x6G08_MD_54_proj[1],x6G08_MD_54_proj[2],pos=4,label="54",col="black")

x6G08_MD_55_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_55.txt")
x6G08_MD_55_proj<-project.pca(x6G08_MD_55_raw,pcadata)
points(x6G08_MD_55_proj[1],x6G08_MD_55_proj[2],pch=20)
text(x6G08_MD_55_proj[1],x6G08_MD_55_proj[2],pos=4,label="55",col="black")

x6G08_MD_56_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_56.txt")
x6G08_MD_56_proj<-project.pca(x6G08_MD_56_raw,pcadata)
points(x6G08_MD_56_proj[1],x6G08_MD_56_proj[2],pch=20)
text(x6G08_MD_56_proj[1],x6G08_MD_56_proj[2],pos=4,label="56",col="black")

x6G08_MD_57_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_57.txt")
x6G08_MD_57_proj<-project.pca(x6G08_MD_57_raw,pcadata)
points(x6G08_MD_57_proj[1],x6G08_MD_57_proj[2],pch=20)
text(x6G08_MD_57_proj[1],x6G08_MD_57_proj[2],pos=4,label="57",col="black")

x6G08_MD_58_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_58.txt")
x6G08_MD_58_proj<-project.pca(x6G08_MD_58_raw,pcadata)
points(x6G08_MD_58_proj[1],x6G08_MD_58_proj[2],pch=20)
text(x6G08_MD_58_proj[1],x6G08_MD_58_proj[2],pos=4,label="58",col="blue")

x6G08_MD_59_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_59.txt")
x6G08_MD_59_proj<-project.pca(x6G08_MD_59_raw,pcadata)
points(x6G08_MD_59_proj[1],x6G08_MD_59_proj[2],pch=20)
text(x6G08_MD_59_proj[1],x6G08_MD_59_proj[2],pos=4,label="59",col="black")

x6G08_MD_60_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_60.txt")
x6G08_MD_60_proj<-project.pca(x6G08_MD_60_raw,pcadata)
points(x6G08_MD_60_proj[1],x6G08_MD_60_proj[2],pch=20)
text(x6G08_MD_60_proj[1],x6G08_MD_60_proj[2],pos=4,label="60",col="black")

x6G08_MD_61_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_61.txt")
x6G08_MD_61_proj<-project.pca(x6G08_MD_61_raw,pcadata)
points(x6G08_MD_61_proj[1],x6G08_MD_61_proj[2],pch=20)
text(x6G08_MD_61_proj[1],x6G08_MD_61_proj[2],pos=4,label="61",col="black")

x6G08_MD_62_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_62.txt")
x6G08_MD_62_proj<-project.pca(x6G08_MD_62_raw,pcadata)
points(x6G08_MD_62_proj[1],x6G08_MD_62_proj[2],pch=20)
text(x6G08_MD_62_proj[1],x6G08_MD_62_proj[2],pos=4,label="62",col="black")

x6G08_MD_63_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_63.txt")
x6G08_MD_63_proj<-project.pca(x6G08_MD_63_raw,pcadata)
points(x6G08_MD_63_proj[1],x6G08_MD_63_proj[2],pch=20)
text(x6G08_MD_63_proj[1],x6G08_MD_63_proj[2],pos=4,label="63",col="black")

x6G08_MD_64_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_64.txt")
x6G08_MD_64_proj<-project.pca(x6G08_MD_64_raw,pcadata)
points(x6G08_MD_64_proj[1],x6G08_MD_64_proj[2],pch=20)
text(x6G08_MD_64_proj[1],x6G08_MD_64_proj[2],pos=4,label="64",col="black")

x6G08_MD_65_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_65.txt")
x6G08_MD_65_proj<-project.pca(x6G08_MD_65_raw,pcadata)
points(x6G08_MD_65_proj[1],x6G08_MD_65_proj[2],pch=20)
text(x6G08_MD_65_proj[1],x6G08_MD_65_proj[2],pos=4,label="65",col="black")

x6G08_MD_66_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_66.txt")
x6G08_MD_66_proj<-project.pca(x6G08_MD_66_raw,pcadata)
points(x6G08_MD_66_proj[1],x6G08_MD_66_proj[2],pch=20)
text(x6G08_MD_66_proj[1],x6G08_MD_66_proj[2],pos=4,label="66",col="black")

x6G08_MD_67_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_67.txt")
x6G08_MD_67_proj<-project.pca(x6G08_MD_67_raw,pcadata)
points(x6G08_MD_67_proj[1],x6G08_MD_67_proj[2],pch=20)
text(x6G08_MD_67_proj[1],x6G08_MD_67_proj[2],pos=4,label="67",col="green")

x6G08_MD_68_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_68.txt")
x6G08_MD_68_proj<-project.pca(x6G08_MD_68_raw,pcadata)
points(x6G08_MD_68_proj[1],x6G08_MD_68_proj[2],pch=20)
text(x6G08_MD_68_proj[1],x6G08_MD_68_proj[2],pos=4,label="68",col="blue")

x6G08_MD_69_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_69.txt")
x6G08_MD_69_proj<-project.pca(x6G08_MD_69_raw,pcadata)
points(x6G08_MD_69_proj[1],x6G08_MD_69_proj[2],pch=20)
text(x6G08_MD_69_proj[1],x6G08_MD_69_proj[2],pos=4,label="69",col="blue")

x6G08_MD_70_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_70.txt")
x6G08_MD_70_proj<-project.pca(x6G08_MD_70_raw,pcadata)
points(x6G08_MD_70_proj[1],x6G08_MD_70_proj[2],pch=20)
text(x6G08_MD_70_proj[1],x6G08_MD_70_proj[2],pos=4,label="70",col="black")

x6G08_MD_71_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_71.txt")
x6G08_MD_71_proj<-project.pca(x6G08_MD_71_raw,pcadata)
points(x6G08_MD_71_proj[1],x6G08_MD_71_proj[2],pch=20)
text(x6G08_MD_71_proj[1],x6G08_MD_71_proj[2],pos=4,label="71",col="black")

x6G08_MD_72_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_72.txt")
x6G08_MD_72_proj<-project.pca(x6G08_MD_72_raw,pcadata)
points(x6G08_MD_72_proj[1],x6G08_MD_72_proj[2],pch=20)
text(x6G08_MD_72_proj[1],x6G08_MD_72_proj[2],pos=4,label="72",col="black")

x6G08_MD_73_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_73.txt")
x6G08_MD_73_proj<-project.pca(x6G08_MD_73_raw,pcadata)
points(x6G08_MD_73_proj[1],x6G08_MD_73_proj[2],pch=20)
text(x6G08_MD_73_proj[1],x6G08_MD_73_proj[2],pos=4,label="73",col="blue")

x6G08_MD_74_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_74.txt")
x6G08_MD_74_proj<-project.pca(x6G08_MD_74_raw,pcadata)
points(x6G08_MD_74_proj[1],x6G08_MD_74_proj[2],pch=20)
text(x6G08_MD_74_proj[1],x6G08_MD_74_proj[2],pos=4,label="74",col="blue")

x6G08_MD_75_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_75.txt")
x6G08_MD_75_proj<-project.pca(x6G08_MD_75_raw,pcadata)
points(x6G08_MD_75_proj[1],x6G08_MD_75_proj[2],pch=20)
text(x6G08_MD_75_proj[1],x6G08_MD_75_proj[2],pos=4,label="75",col="black")

x6G08_MD_76_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_76.txt")
x6G08_MD_76_proj<-project.pca(x6G08_MD_76_raw,pcadata)
points(x6G08_MD_76_proj[1],x6G08_MD_76_proj[2],pch=20)
text(x6G08_MD_76_proj[1],x6G08_MD_76_proj[2],pos=4,label="76",col="black")

x6G08_MD_77_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_77.txt")
x6G08_MD_77_proj<-project.pca(x6G08_MD_77_raw,pcadata)
points(x6G08_MD_77_proj[1],x6G08_MD_77_proj[2],pch=20)
text(x6G08_MD_77_proj[1],x6G08_MD_77_proj[2],pos=4,label="77",col="black")

x6G08_MD_78_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_78.txt")
x6G08_MD_78_proj<-project.pca(x6G08_MD_78_raw,pcadata)
points(x6G08_MD_78_proj[1],x6G08_MD_78_proj[2],pch=20)
text(x6G08_MD_78_proj[1],x6G08_MD_78_proj[2],pos=4,label="78",col="blue")

x6G08_MD_79_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_79.txt")
x6G08_MD_79_proj<-project.pca(x6G08_MD_79_raw,pcadata)
points(x6G08_MD_79_proj[1],x6G08_MD_79_proj[2],pch=20)
text(x6G08_MD_79_proj[1],x6G08_MD_79_proj[2],pos=4,label="79",col="green")

x6G08_MD_80_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_80.txt")
x6G08_MD_80_proj<-project.pca(x6G08_MD_80_raw,pcadata)
points(x6G08_MD_80_proj[1],x6G08_MD_80_proj[2],pch=20)
text(x6G08_MD_80_proj[1],x6G08_MD_80_proj[2],pos=4,label="80",col="blue")

x6G08_MD_81_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_81.txt")
x6G08_MD_81_proj<-project.pca(x6G08_MD_81_raw,pcadata)
points(x6G08_MD_81_proj[1],x6G08_MD_81_proj[2],pch=20)
text(x6G08_MD_81_proj[1],x6G08_MD_81_proj[2],pos=4,label="81",col="black")

x6G08_MD_82_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_82.txt")
x6G08_MD_82_proj<-project.pca(x6G08_MD_82_raw,pcadata)
points(x6G08_MD_82_proj[1],x6G08_MD_82_proj[2],pch=20)
text(x6G08_MD_82_proj[1],x6G08_MD_82_proj[2],pos=4,label="82",col="black")

x6G08_MD_83_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_83.txt")
x6G08_MD_83_proj<-project.pca(x6G08_MD_83_raw,pcadata)
points(x6G08_MD_83_proj[1],x6G08_MD_83_proj[2],pch=20)
text(x6G08_MD_83_proj[1],x6G08_MD_83_proj[2],pos=4,label="83",col="green")

x6G08_MD_84_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_84.txt")
x6G08_MD_84_proj<-project.pca(x6G08_MD_84_raw,pcadata)
points(x6G08_MD_84_proj[1],x6G08_MD_84_proj[2],pch=20)
text(x6G08_MD_84_proj[1],x6G08_MD_84_proj[2],pos=4,label="84",col="blue")

x6G08_MD_85_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_85.txt")
x6G08_MD_85_proj<-project.pca(x6G08_MD_85_raw,pcadata)
points(x6G08_MD_85_proj[1],x6G08_MD_85_proj[2],pch=20)
text(x6G08_MD_85_proj[1],x6G08_MD_85_proj[2],pos=4,label="85",col="black")

x6G08_MD_86_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_86.txt")
x6G08_MD_86_proj<-project.pca(x6G08_MD_86_raw,pcadata)
points(x6G08_MD_86_proj[1],x6G08_MD_86_proj[2],pch=20)
text(x6G08_MD_86_proj[1],x6G08_MD_86_proj[2],pos=4,label="86",col="blue")

x6G08_MD_87_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_87.txt")
x6G08_MD_87_proj<-project.pca(x6G08_MD_87_raw,pcadata)
points(x6G08_MD_87_proj[1],x6G08_MD_87_proj[2],pch=20)
text(x6G08_MD_87_proj[1],x6G08_MD_87_proj[2],pos=4,label="87",col="blue")

x6G08_MD_88_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_88.txt")
x6G08_MD_88_proj<-project.pca(x6G08_MD_88_raw,pcadata)
points(x6G08_MD_88_proj[1],x6G08_MD_88_proj[2],pch=20)
text(x6G08_MD_88_proj[1],x6G08_MD_88_proj[2],pos=4,label="88",col="green")

x6G08_MD_89_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_89.txt")
x6G08_MD_89_proj<-project.pca(x6G08_MD_89_raw,pcadata)
points(x6G08_MD_89_proj[1],x6G08_MD_89_proj[2],pch=20)
text(x6G08_MD_89_proj[1],x6G08_MD_89_proj[2],pos=4,label="89",col="green")

x6G08_MD_90_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_90.txt")
x6G08_MD_90_proj<-project.pca(x6G08_MD_90_raw,pcadata)
points(x6G08_MD_90_proj[1],x6G08_MD_90_proj[2],pch=20)
text(x6G08_MD_90_proj[1],x6G08_MD_90_proj[2],pos=4,label="90",col="green")

x6G08_MD_91_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_91.txt")
x6G08_MD_91_proj<-project.pca(x6G08_MD_91_raw,pcadata)
points(x6G08_MD_91_proj[1],x6G08_MD_91_proj[2],pch=20)
text(x6G08_MD_91_proj[1],x6G08_MD_91_proj[2],pos=4,label="91",col="green")

x6G08_MD_92_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_92.txt")
x6G08_MD_92_proj<-project.pca(x6G08_MD_92_raw,pcadata)
points(x6G08_MD_92_proj[1],x6G08_MD_92_proj[2],pch=20)
text(x6G08_MD_92_proj[1],x6G08_MD_92_proj[2],pos=4,label="92",col="green")

x6G08_MD_93_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_93.txt")
x6G08_MD_93_proj<-project.pca(x6G08_MD_93_raw,pcadata)
points(x6G08_MD_93_proj[1],x6G08_MD_93_proj[2],pch=20)
text(x6G08_MD_93_proj[1],x6G08_MD_93_proj[2],pos=4,label="93",col="blue")

x6G08_MD_94_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_94.txt")
x6G08_MD_94_proj<-project.pca(x6G08_MD_94_raw,pcadata)
points(x6G08_MD_94_proj[1],x6G08_MD_94_proj[2],pch=20)
text(x6G08_MD_94_proj[1],x6G08_MD_94_proj[2],pos=4,label="94",col="black")

x6G08_MD_95_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_95.txt")
x6G08_MD_95_proj<-project.pca(x6G08_MD_95_raw,pcadata)
points(x6G08_MD_95_proj[1],x6G08_MD_95_proj[2],pch=20)
text(x6G08_MD_95_proj[1],x6G08_MD_95_proj[2],pos=4,label="95",col="green")

x6G08_MD_96_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_96.txt")
x6G08_MD_96_proj<-project.pca(x6G08_MD_96_raw,pcadata)
points(x6G08_MD_96_proj[1],x6G08_MD_96_proj[2],pch=20)
text(x6G08_MD_96_proj[1],x6G08_MD_96_proj[2],pos=4,label="96",col="green")

x6G08_MD_97_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_97.txt")
x6G08_MD_97_proj<-project.pca(x6G08_MD_97_raw,pcadata)
points(x6G08_MD_97_proj[1],x6G08_MD_97_proj[2],pch=20)
text(x6G08_MD_97_proj[1],x6G08_MD_97_proj[2],pos=4,label="97",col="green")

x6G08_MD_98_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_98.txt")
x6G08_MD_98_proj<-project.pca(x6G08_MD_98_raw,pcadata)
points(x6G08_MD_98_proj[1],x6G08_MD_98_proj[2],pch=20)
text(x6G08_MD_98_proj[1],x6G08_MD_98_proj[2],pos=4,label="98",col="black")

x6G08_MD_99_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_99.txt")
x6G08_MD_99_proj<-project.pca(x6G08_MD_99_raw,pcadata)
points(x6G08_MD_99_proj[1],x6G08_MD_99_proj[2],pch=20)
text(x6G08_MD_99_proj[1],x6G08_MD_99_proj[2],pos=4,label="99",col="blue")

x6G08_MD_100_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_100.txt")
x6G08_MD_100_proj<-project.pca(x6G08_MD_100_raw,pcadata)
points(x6G08_MD_100_proj[1],x6G08_MD_100_proj[2],pch=20)
text(x6G08_MD_100_proj[1],x6G08_MD_100_proj[2],pos=4,label="100",col="black")

x6G08_MD_101_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_101.txt")
x6G08_MD_101_proj<-project.pca(x6G08_MD_101_raw,pcadata)
points(x6G08_MD_101_proj[1],x6G08_MD_101_proj[2],pch=20)
text(x6G08_MD_101_proj[1],x6G08_MD_101_proj[2],pos=4,label="101",col="black")

x6G08_MD_102_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_102.txt")
x6G08_MD_102_proj<-project.pca(x6G08_MD_102_raw,pcadata)
points(x6G08_MD_102_proj[1],x6G08_MD_102_proj[2],pch=20)
text(x6G08_MD_102_proj[1],x6G08_MD_102_proj[2],pos=4,label="102",col="green")

x6G08_MD_103_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_103.txt")
x6G08_MD_103_proj<-project.pca(x6G08_MD_103_raw,pcadata)
points(x6G08_MD_103_proj[1],x6G08_MD_103_proj[2],pch=20)
text(x6G08_MD_103_proj[1],x6G08_MD_103_proj[2],pos=4,label="103",col="blue")

x6G08_MD_104_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_104.txt")
x6G08_MD_104_proj<-project.pca(x6G08_MD_104_raw,pcadata)
points(x6G08_MD_104_proj[1],x6G08_MD_104_proj[2],pch=20)
text(x6G08_MD_104_proj[1],x6G08_MD_104_proj[2],pos=4,label="104",col="blue")

x6G08_MD_105_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_105.txt")
x6G08_MD_105_proj<-project.pca(x6G08_MD_105_raw,pcadata)
points(x6G08_MD_105_proj[1],x6G08_MD_105_proj[2],pch=20)
text(x6G08_MD_105_proj[1],x6G08_MD_105_proj[2],pos=4,label="105",col="blue")

x6G08_MD_106_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_106.txt")
x6G08_MD_106_proj<-project.pca(x6G08_MD_106_raw,pcadata)
points(x6G08_MD_106_proj[1],x6G08_MD_106_proj[2],pch=20)
text(x6G08_MD_106_proj[1],x6G08_MD_106_proj[2],pos=4,label="106",col="black")

x6G08_MD_107_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_107.txt")
x6G08_MD_107_proj<-project.pca(x6G08_MD_107_raw,pcadata)
points(x6G08_MD_107_proj[1],x6G08_MD_107_proj[2],pch=20)
text(x6G08_MD_107_proj[1],x6G08_MD_107_proj[2],pos=4,label="107",col="green")

x6G08_MD_108_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_108.txt")
x6G08_MD_108_proj<-project.pca(x6G08_MD_108_raw,pcadata)
points(x6G08_MD_108_proj[1],x6G08_MD_108_proj[2],pch=20)
text(x6G08_MD_108_proj[1],x6G08_MD_108_proj[2],pos=4,label="108",col="black")

x6G08_MD_109_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_109.txt")
x6G08_MD_109_proj<-project.pca(x6G08_MD_109_raw,pcadata)
points(x6G08_MD_109_proj[1],x6G08_MD_109_proj[2],pch=20)
text(x6G08_MD_109_proj[1],x6G08_MD_109_proj[2],pos=4,label="109",col="black")

x6G08_MD_110_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_110.txt")
x6G08_MD_110_proj<-project.pca(x6G08_MD_110_raw,pcadata)
points(x6G08_MD_110_proj[1],x6G08_MD_110_proj[2],pch=20)
text(x6G08_MD_110_proj[1],x6G08_MD_110_proj[2],pos=4,label="110",col="black")

x6G08_MD_111_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_111.txt")
x6G08_MD_111_proj<-project.pca(x6G08_MD_111_raw,pcadata)
points(x6G08_MD_111_proj[1],x6G08_MD_111_proj[2],pch=20)
text(x6G08_MD_111_proj[1],x6G08_MD_111_proj[2],pos=4,label="111",col="black")

x6G08_MD_112_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_112.txt")
x6G08_MD_112_proj<-project.pca(x6G08_MD_112_raw,pcadata)
points(x6G08_MD_112_proj[1],x6G08_MD_112_proj[2],pch=20)
text(x6G08_MD_112_proj[1],x6G08_MD_112_proj[2],pos=4,label="112",col="blue")

x6G08_MD_113_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_113.txt")
x6G08_MD_113_proj<-project.pca(x6G08_MD_113_raw,pcadata)
points(x6G08_MD_113_proj[1],x6G08_MD_113_proj[2],pch=20)
text(x6G08_MD_113_proj[1],x6G08_MD_113_proj[2],pos=4,label="113",col="blue")

x6G08_MD_114_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_114.txt")
x6G08_MD_114_proj<-project.pca(x6G08_MD_114_raw,pcadata)
points(x6G08_MD_114_proj[1],x6G08_MD_114_proj[2],pch=20)
text(x6G08_MD_114_proj[1],x6G08_MD_114_proj[2],pos=4,label="114",col="blue")

x6G08_MD_115_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_115.txt")
x6G08_MD_115_proj<-project.pca(x6G08_MD_115_raw,pcadata)
points(x6G08_MD_115_proj[1],x6G08_MD_115_proj[2],pch=20)
text(x6G08_MD_115_proj[1],x6G08_MD_115_proj[2],pos=4,label="115",col="black")

x6G08_MD_116_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_116.txt")
x6G08_MD_116_proj<-project.pca(x6G08_MD_116_raw,pcadata)
points(x6G08_MD_116_proj[1],x6G08_MD_116_proj[2],pch=20)
text(x6G08_MD_116_proj[1],x6G08_MD_116_proj[2],pos=4,label="116",col="black")

x6G08_MD_117_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_117.txt")
x6G08_MD_117_proj<-project.pca(x6G08_MD_117_raw,pcadata)
points(x6G08_MD_117_proj[1],x6G08_MD_117_proj[2],pch=20)
text(x6G08_MD_117_proj[1],x6G08_MD_117_proj[2],pos=4,label="117",col="blue")

x6G08_MD_118_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_118.txt")
x6G08_MD_118_proj<-project.pca(x6G08_MD_118_raw,pcadata)
points(x6G08_MD_118_proj[1],x6G08_MD_118_proj[2],pch=20)
text(x6G08_MD_118_proj[1],x6G08_MD_118_proj[2],pos=4,label="118",col="green")

x6G08_MD_119_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_119.txt")
x6G08_MD_119_proj<-project.pca(x6G08_MD_119_raw,pcadata)
points(x6G08_MD_119_proj[1],x6G08_MD_119_proj[2],pch=20)
text(x6G08_MD_119_proj[1],x6G08_MD_119_proj[2],pos=4,label="119",col="blue")

x6G08_MD_120_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_120.txt")
x6G08_MD_120_proj<-project.pca(x6G08_MD_120_raw,pcadata)
points(x6G08_MD_120_proj[1],x6G08_MD_120_proj[2],pch=20)
text(x6G08_MD_120_proj[1],x6G08_MD_120_proj[2],pos=4,label="120",col="black")

x6G08_MD_121_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_121.txt")
x6G08_MD_121_proj<-project.pca(x6G08_MD_121_raw,pcadata)
points(x6G08_MD_121_proj[1],x6G08_MD_121_proj[2],pch=20)
text(x6G08_MD_121_proj[1],x6G08_MD_121_proj[2],pos=4,label="121",col="black")

x6G08_MD_122_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_122.txt")
x6G08_MD_122_proj<-project.pca(x6G08_MD_122_raw,pcadata)
points(x6G08_MD_122_proj[1],x6G08_MD_122_proj[2],pch=20)
text(x6G08_MD_122_proj[1],x6G08_MD_122_proj[2],pos=4,label="122",col="black")

x6G08_MD_123_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_123.txt")
x6G08_MD_123_proj<-project.pca(x6G08_MD_123_raw,pcadata)
points(x6G08_MD_123_proj[1],x6G08_MD_123_proj[2],pch=20)
text(x6G08_MD_123_proj[1],x6G08_MD_123_proj[2],pos=4,label="123",col="black")

x6G08_MD_124_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_124.txt")
x6G08_MD_124_proj<-project.pca(x6G08_MD_124_raw,pcadata)
points(x6G08_MD_124_proj[1],x6G08_MD_124_proj[2],pch=20)
text(x6G08_MD_124_proj[1],x6G08_MD_124_proj[2],pos=4,label="124",col="black")

x6G08_MD_125_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_125.txt")
x6G08_MD_125_proj<-project.pca(x6G08_MD_125_raw,pcadata)
points(x6G08_MD_125_proj[1],x6G08_MD_125_proj[2],pch=20)
text(x6G08_MD_125_proj[1],x6G08_MD_125_proj[2],pos=4,label="125",col="blue")

x6G08_MD_126_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_126.txt")
x6G08_MD_126_proj<-project.pca(x6G08_MD_126_raw,pcadata)
points(x6G08_MD_126_proj[1],x6G08_MD_126_proj[2],pch=20)
text(x6G08_MD_126_proj[1],x6G08_MD_126_proj[2],pos=4,label="126",col="black")

x6G08_MD_127_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_127.txt")
x6G08_MD_127_proj<-project.pca(x6G08_MD_127_raw,pcadata)
points(x6G08_MD_127_proj[1],x6G08_MD_127_proj[2],pch=20)
text(x6G08_MD_127_proj[1],x6G08_MD_127_proj[2],pos=4,label="127",col="black")

x6G08_MD_128_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_128.txt")
x6G08_MD_128_proj<-project.pca(x6G08_MD_128_raw,pcadata)
points(x6G08_MD_128_proj[1],x6G08_MD_128_proj[2],pch=20)
text(x6G08_MD_128_proj[1],x6G08_MD_128_proj[2],pos=4,label="128",col="black")

x6G08_MD_129_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_129.txt")
x6G08_MD_129_proj<-project.pca(x6G08_MD_129_raw,pcadata)
points(x6G08_MD_129_proj[1],x6G08_MD_129_proj[2],pch=20)
text(x6G08_MD_129_proj[1],x6G08_MD_129_proj[2],pos=4,label="129",col="black")

x6G08_MD_130_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_130.txt")
x6G08_MD_130_proj<-project.pca(x6G08_MD_130_raw,pcadata)
points(x6G08_MD_130_proj[1],x6G08_MD_130_proj[2],pch=20)
text(x6G08_MD_130_proj[1],x6G08_MD_130_proj[2],pos=4,label="130",col="black")

x6G08_MD_131_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_131.txt")
x6G08_MD_131_proj<-project.pca(x6G08_MD_131_raw,pcadata)
points(x6G08_MD_131_proj[1],x6G08_MD_131_proj[2],pch=20)
text(x6G08_MD_131_proj[1],x6G08_MD_131_proj[2],pos=4,label="131",col="black")

x6G08_MD_132_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_132.txt")
x6G08_MD_132_proj<-project.pca(x6G08_MD_132_raw,pcadata)
points(x6G08_MD_132_proj[1],x6G08_MD_132_proj[2],pch=20)
text(x6G08_MD_132_proj[1],x6G08_MD_132_proj[2],pos=4,label="132",col="black")

x6G08_MD_133_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_133.txt")
x6G08_MD_133_proj<-project.pca(x6G08_MD_133_raw,pcadata)
points(x6G08_MD_133_proj[1],x6G08_MD_133_proj[2],pch=20)
text(x6G08_MD_133_proj[1],x6G08_MD_133_proj[2],pos=4,label="133",col="black")

x6G08_MD_134_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_134.txt")
x6G08_MD_134_proj<-project.pca(x6G08_MD_134_raw,pcadata)
points(x6G08_MD_134_proj[1],x6G08_MD_134_proj[2],pch=20)
text(x6G08_MD_134_proj[1],x6G08_MD_134_proj[2],pos=4,label="134",col="black")

x6G08_MD_135_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_135.txt")
x6G08_MD_135_proj<-project.pca(x6G08_MD_135_raw,pcadata)
points(x6G08_MD_135_proj[1],x6G08_MD_135_proj[2],pch=20)
text(x6G08_MD_135_proj[1],x6G08_MD_135_proj[2],pos=4,label="135",col="black")

x6G08_MD_136_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_136.txt")
x6G08_MD_136_proj<-project.pca(x6G08_MD_136_raw,pcadata)
points(x6G08_MD_136_proj[1],x6G08_MD_136_proj[2],pch=20)
text(x6G08_MD_136_proj[1],x6G08_MD_136_proj[2],pos=4,label="136",col="black")

x6G08_MD_137_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_137.txt")
x6G08_MD_137_proj<-project.pca(x6G08_MD_137_raw,pcadata)
points(x6G08_MD_137_proj[1],x6G08_MD_137_proj[2],pch=20)
text(x6G08_MD_137_proj[1],x6G08_MD_137_proj[2],pos=4,label="137",col="black")

x6G08_MD_138_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_138.txt")
x6G08_MD_138_proj<-project.pca(x6G08_MD_138_raw,pcadata)
points(x6G08_MD_138_proj[1],x6G08_MD_138_proj[2],pch=20)
text(x6G08_MD_138_proj[1],x6G08_MD_138_proj[2],pos=4,label="138",col="black")

x6G08_MD_139_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_139.txt")
x6G08_MD_139_proj<-project.pca(x6G08_MD_139_raw,pcadata)
points(x6G08_MD_139_proj[1],x6G08_MD_139_proj[2],pch=20)
text(x6G08_MD_139_proj[1],x6G08_MD_139_proj[2],pos=4,label="139",col="black")

x6G08_MD_140_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_140.txt")
x6G08_MD_140_proj<-project.pca(x6G08_MD_140_raw,pcadata)
points(x6G08_MD_140_proj[1],x6G08_MD_140_proj[2],pch=20)
text(x6G08_MD_140_proj[1],x6G08_MD_140_proj[2],pos=4,label="140",col="black")

x6G08_MD_141_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_141.txt")
x6G08_MD_141_proj<-project.pca(x6G08_MD_141_raw,pcadata)
points(x6G08_MD_141_proj[1],x6G08_MD_141_proj[2],pch=20)
text(x6G08_MD_141_proj[1],x6G08_MD_141_proj[2],pos=4,label="141",col="green")

x6G08_MD_142_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_142.txt")
x6G08_MD_142_proj<-project.pca(x6G08_MD_142_raw,pcadata)
points(x6G08_MD_142_proj[1],x6G08_MD_142_proj[2],pch=20)
text(x6G08_MD_142_proj[1],x6G08_MD_142_proj[2],pos=4,label="142",col="black")

x6G08_MD_143_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_143.txt")
x6G08_MD_143_proj<-project.pca(x6G08_MD_143_raw,pcadata)
points(x6G08_MD_143_proj[1],x6G08_MD_143_proj[2],pch=20)
text(x6G08_MD_143_proj[1],x6G08_MD_143_proj[2],pos=4,label="143",col="black")

x6G08_MD_144_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_144.txt")
x6G08_MD_144_proj<-project.pca(x6G08_MD_144_raw,pcadata)
points(x6G08_MD_144_proj[1],x6G08_MD_144_proj[2],pch=20)
text(x6G08_MD_144_proj[1],x6G08_MD_144_proj[2],pos=4,label="144",col="blue")

x6G08_MD_145_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_145.txt")
x6G08_MD_145_proj<-project.pca(x6G08_MD_145_raw,pcadata)
points(x6G08_MD_145_proj[1],x6G08_MD_145_proj[2],pch=20)
text(x6G08_MD_145_proj[1],x6G08_MD_145_proj[2],pos=4,label="145",col="black")

x6G08_MD_146_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_146.txt")
x6G08_MD_146_proj<-project.pca(x6G08_MD_146_raw,pcadata)
points(x6G08_MD_146_proj[1],x6G08_MD_146_proj[2],pch=20)
text(x6G08_MD_146_proj[1],x6G08_MD_146_proj[2],pos=4,label="146",col="black")

x6G08_MD_147_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_147.txt")
x6G08_MD_147_proj<-project.pca(x6G08_MD_147_raw,pcadata)
points(x6G08_MD_147_proj[1],x6G08_MD_147_proj[2],pch=20)
text(x6G08_MD_147_proj[1],x6G08_MD_147_proj[2],pos=4,label="147",col="green")

x6G08_MD_148_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_148.txt")
x6G08_MD_148_proj<-project.pca(x6G08_MD_148_raw,pcadata)
points(x6G08_MD_148_proj[1],x6G08_MD_148_proj[2],pch=20)
text(x6G08_MD_148_proj[1],x6G08_MD_148_proj[2],pos=4,label="148",col="green")

x6G08_MD_149_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_149.txt")
x6G08_MD_149_proj<-project.pca(x6G08_MD_149_raw,pcadata)
points(x6G08_MD_149_proj[1],x6G08_MD_149_proj[2],pch=20)
text(x6G08_MD_149_proj[1],x6G08_MD_149_proj[2],pos=4,label="149",col="blue")

x6G08_MD_150_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_150.txt")
x6G08_MD_150_proj<-project.pca(x6G08_MD_150_raw,pcadata)
points(x6G08_MD_150_proj[1],x6G08_MD_150_proj[2],pch=20)
text(x6G08_MD_150_proj[1],x6G08_MD_150_proj[2],pos=4,label="150",col="red")

x6G08_MD_151_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_151.txt")
x6G08_MD_151_proj<-project.pca(x6G08_MD_151_raw,pcadata)
points(x6G08_MD_151_proj[1],x6G08_MD_151_proj[2],pch=20)
text(x6G08_MD_151_proj[1],x6G08_MD_151_proj[2],pos=4,label="151",col="black")

x6G08_MD_152_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_152.txt")
x6G08_MD_152_proj<-project.pca(x6G08_MD_152_raw,pcadata)
points(x6G08_MD_152_proj[1],x6G08_MD_152_proj[2],pch=20)
text(x6G08_MD_152_proj[1],x6G08_MD_152_proj[2],pos=4,label="152",col="black")

x6G08_MD_153_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_153.txt")
x6G08_MD_153_proj<-project.pca(x6G08_MD_153_raw,pcadata)
points(x6G08_MD_153_proj[1],x6G08_MD_153_proj[2],pch=20)
text(x6G08_MD_153_proj[1],x6G08_MD_153_proj[2],pos=4,label="153",col="black")

x6G08_MD_154_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_154.txt")
x6G08_MD_154_proj<-project.pca(x6G08_MD_154_raw,pcadata)
points(x6G08_MD_154_proj[1],x6G08_MD_154_proj[2],pch=20)
text(x6G08_MD_154_proj[1],x6G08_MD_154_proj[2],pos=4,label="154",col="black")

x6G08_MD_155_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_155.txt")
x6G08_MD_155_proj<-project.pca(x6G08_MD_155_raw,pcadata)
points(x6G08_MD_155_proj[1],x6G08_MD_155_proj[2],pch=20)
text(x6G08_MD_155_proj[1],x6G08_MD_155_proj[2],pos=4,label="155",col="black")

x6G08_MD_156_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_156.txt")
x6G08_MD_156_proj<-project.pca(x6G08_MD_156_raw,pcadata)
points(x6G08_MD_156_proj[1],x6G08_MD_156_proj[2],pch=20)
text(x6G08_MD_156_proj[1],x6G08_MD_156_proj[2],pos=4,label="156",col="black")

x6G08_MD_157_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_157.txt")
x6G08_MD_157_proj<-project.pca(x6G08_MD_157_raw,pcadata)
points(x6G08_MD_157_proj[1],x6G08_MD_157_proj[2],pch=20)
text(x6G08_MD_157_proj[1],x6G08_MD_157_proj[2],pos=4,label="157",col="black")

x6G08_MD_158_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_158.txt")
x6G08_MD_158_proj<-project.pca(x6G08_MD_158_raw,pcadata)
points(x6G08_MD_158_proj[1],x6G08_MD_158_proj[2],pch=20)
text(x6G08_MD_158_proj[1],x6G08_MD_158_proj[2],pos=4,label="158",col="black")

x6G08_MD_159_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_159.txt")
x6G08_MD_159_proj<-project.pca(x6G08_MD_159_raw,pcadata)
points(x6G08_MD_159_proj[1],x6G08_MD_159_proj[2],pch=20)
text(x6G08_MD_159_proj[1],x6G08_MD_159_proj[2],pos=4,label="159",col="black")

x6G08_MD_160_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_160.txt")
x6G08_MD_160_proj<-project.pca(x6G08_MD_160_raw,pcadata)
points(x6G08_MD_160_proj[1],x6G08_MD_160_proj[2],pch=20)
text(x6G08_MD_160_proj[1],x6G08_MD_160_proj[2],pos=4,label="160",col="black")

x6G08_MD_161_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_161.txt")
x6G08_MD_161_proj<-project.pca(x6G08_MD_161_raw,pcadata)
points(x6G08_MD_161_proj[1],x6G08_MD_161_proj[2],pch=20)
text(x6G08_MD_161_proj[1],x6G08_MD_161_proj[2],pos=4,label="161",col="black")

x6G08_MD_162_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_162.txt")
x6G08_MD_162_proj<-project.pca(x6G08_MD_162_raw,pcadata)
points(x6G08_MD_162_proj[1],x6G08_MD_162_proj[2],pch=20)
text(x6G08_MD_162_proj[1],x6G08_MD_162_proj[2],pos=4,label="162",col="black")

x6G08_MD_163_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_163.txt")
x6G08_MD_163_proj<-project.pca(x6G08_MD_163_raw,pcadata)
points(x6G08_MD_163_proj[1],x6G08_MD_163_proj[2],pch=20)
text(x6G08_MD_163_proj[1],x6G08_MD_163_proj[2],pos=4,label="163",col="black")

x6G08_MD_164_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_164.txt")
x6G08_MD_164_proj<-project.pca(x6G08_MD_164_raw,pcadata)
points(x6G08_MD_164_proj[1],x6G08_MD_164_proj[2],pch=20)
text(x6G08_MD_164_proj[1],x6G08_MD_164_proj[2],pos=4,label="164",col="black")

x6G08_MD_165_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_165.txt")
x6G08_MD_165_proj<-project.pca(x6G08_MD_165_raw,pcadata)
points(x6G08_MD_165_proj[1],x6G08_MD_165_proj[2],pch=20)
text(x6G08_MD_165_proj[1],x6G08_MD_165_proj[2],pos=4,label="165",col="green")

x6G08_MD_166_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_166.txt")
x6G08_MD_166_proj<-project.pca(x6G08_MD_166_raw,pcadata)
points(x6G08_MD_166_proj[1],x6G08_MD_166_proj[2],pch=20)
text(x6G08_MD_166_proj[1],x6G08_MD_166_proj[2],pos=4,label="166",col="black")

x6G08_MD_167_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_167.txt")
x6G08_MD_167_proj<-project.pca(x6G08_MD_167_raw,pcadata)
points(x6G08_MD_167_proj[1],x6G08_MD_167_proj[2],pch=20)
text(x6G08_MD_167_proj[1],x6G08_MD_167_proj[2],pos=4,label="167",col="black")

x6G08_MD_168_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_168.txt")
x6G08_MD_168_proj<-project.pca(x6G08_MD_168_raw,pcadata)
points(x6G08_MD_168_proj[1],x6G08_MD_168_proj[2],pch=20)
text(x6G08_MD_168_proj[1],x6G08_MD_168_proj[2],pos=4,label="168",col="black")

x6G08_MD_169_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_169.txt")
x6G08_MD_169_proj<-project.pca(x6G08_MD_169_raw,pcadata)
points(x6G08_MD_169_proj[1],x6G08_MD_169_proj[2],pch=20)
text(x6G08_MD_169_proj[1],x6G08_MD_169_proj[2],pos=4,label="169",col="black")

x6G08_MD_170_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_170.txt")
x6G08_MD_170_proj<-project.pca(x6G08_MD_170_raw,pcadata)
points(x6G08_MD_170_proj[1],x6G08_MD_170_proj[2],pch=20)
text(x6G08_MD_170_proj[1],x6G08_MD_170_proj[2],pos=4,label="170",col="red")

x6G08_MD_171_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_171.txt")
x6G08_MD_171_proj<-project.pca(x6G08_MD_171_raw,pcadata)
points(x6G08_MD_171_proj[1],x6G08_MD_171_proj[2],pch=20)
text(x6G08_MD_171_proj[1],x6G08_MD_171_proj[2],pos=4,label="171",col="red")

x6G08_MD_172_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_172.txt")
x6G08_MD_172_proj<-project.pca(x6G08_MD_172_raw,pcadata)
points(x6G08_MD_172_proj[1],x6G08_MD_172_proj[2],pch=20)
text(x6G08_MD_172_proj[1],x6G08_MD_172_proj[2],pos=4,label="172",col="blue")

x6G08_MD_173_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_173.txt")
x6G08_MD_173_proj<-project.pca(x6G08_MD_173_raw,pcadata)
points(x6G08_MD_173_proj[1],x6G08_MD_173_proj[2],pch=20)
text(x6G08_MD_173_proj[1],x6G08_MD_173_proj[2],pos=4,label="173",col="black")

x6G08_MD_174_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_174.txt")
x6G08_MD_174_proj<-project.pca(x6G08_MD_174_raw,pcadata)
points(x6G08_MD_174_proj[1],x6G08_MD_174_proj[2],pch=20)
text(x6G08_MD_174_proj[1],x6G08_MD_174_proj[2],pos=4,label="174",col="blue")

x6G08_MD_175_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_175.txt")
x6G08_MD_175_proj<-project.pca(x6G08_MD_175_raw,pcadata)
points(x6G08_MD_175_proj[1],x6G08_MD_175_proj[2],pch=20)
text(x6G08_MD_175_proj[1],x6G08_MD_175_proj[2],pos=4,label="175",col="blue")

x6G08_MD_176_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_176.txt")
x6G08_MD_176_proj<-project.pca(x6G08_MD_176_raw,pcadata)
points(x6G08_MD_176_proj[1],x6G08_MD_176_proj[2],pch=20)
text(x6G08_MD_176_proj[1],x6G08_MD_176_proj[2],pos=4,label="176",col="blue")

x6G08_MD_177_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_177.txt")
x6G08_MD_177_proj<-project.pca(x6G08_MD_177_raw,pcadata)
points(x6G08_MD_177_proj[1],x6G08_MD_177_proj[2],pch=20)
text(x6G08_MD_177_proj[1],x6G08_MD_177_proj[2],pos=4,label="177",col="blue")

x6G08_MD_178_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_178.txt")
x6G08_MD_178_proj<-project.pca(x6G08_MD_178_raw,pcadata)
points(x6G08_MD_178_proj[1],x6G08_MD_178_proj[2],pch=20)
text(x6G08_MD_178_proj[1],x6G08_MD_178_proj[2],pos=4,label="178",col="blue")

x6G08_MD_179_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_179.txt")
x6G08_MD_179_proj<-project.pca(x6G08_MD_179_raw,pcadata)
points(x6G08_MD_179_proj[1],x6G08_MD_179_proj[2],pch=20)
text(x6G08_MD_179_proj[1],x6G08_MD_179_proj[2],pos=4,label="179",col="blue")

x6G08_MD_180_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_180.txt")
x6G08_MD_180_proj<-project.pca(x6G08_MD_180_raw,pcadata)
points(x6G08_MD_180_proj[1],x6G08_MD_180_proj[2],pch=20)
text(x6G08_MD_180_proj[1],x6G08_MD_180_proj[2],pos=4,label="180",col="black")

x6G08_MD_181_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_181.txt")
x6G08_MD_181_proj<-project.pca(x6G08_MD_181_raw,pcadata)
points(x6G08_MD_181_proj[1],x6G08_MD_181_proj[2],pch=20)
text(x6G08_MD_181_proj[1],x6G08_MD_181_proj[2],pos=4,label="181",col="black")

x6G08_MD_182_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_182.txt")
x6G08_MD_182_proj<-project.pca(x6G08_MD_182_raw,pcadata)
points(x6G08_MD_182_proj[1],x6G08_MD_182_proj[2],pch=20)
text(x6G08_MD_182_proj[1],x6G08_MD_182_proj[2],pos=4,label="182",col="black")

x6G08_MD_183_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_183.txt")
x6G08_MD_183_proj<-project.pca(x6G08_MD_183_raw,pcadata)
points(x6G08_MD_183_proj[1],x6G08_MD_183_proj[2],pch=20)
text(x6G08_MD_183_proj[1],x6G08_MD_183_proj[2],pos=4,label="183",col="black")

x6G08_MD_184_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_184.txt")
x6G08_MD_184_proj<-project.pca(x6G08_MD_184_raw,pcadata)
points(x6G08_MD_184_proj[1],x6G08_MD_184_proj[2],pch=20)
text(x6G08_MD_184_proj[1],x6G08_MD_184_proj[2],pos=4,label="184",col="green")

x6G08_MD_185_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_185.txt")
x6G08_MD_185_proj<-project.pca(x6G08_MD_185_raw,pcadata)
points(x6G08_MD_185_proj[1],x6G08_MD_185_proj[2],pch=20)
text(x6G08_MD_185_proj[1],x6G08_MD_185_proj[2],pos=4,label="185",col="blue")

x6G08_MD_186_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_186.txt")
x6G08_MD_186_proj<-project.pca(x6G08_MD_186_raw,pcadata)
points(x6G08_MD_186_proj[1],x6G08_MD_186_proj[2],pch=20)
text(x6G08_MD_186_proj[1],x6G08_MD_186_proj[2],pos=4,label="186",col="red")

x6G08_MD_187_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_187.txt")
x6G08_MD_187_proj<-project.pca(x6G08_MD_187_raw,pcadata)
points(x6G08_MD_187_proj[1],x6G08_MD_187_proj[2],pch=20)
text(x6G08_MD_187_proj[1],x6G08_MD_187_proj[2],pos=4,label="187",col="blue")

x6G08_MD_188_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_188.txt")
x6G08_MD_188_proj<-project.pca(x6G08_MD_188_raw,pcadata)
points(x6G08_MD_188_proj[1],x6G08_MD_188_proj[2],pch=20)
text(x6G08_MD_188_proj[1],x6G08_MD_188_proj[2],pos=4,label="188",col="blue")

x6G08_MD_189_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_189.txt")
x6G08_MD_189_proj<-project.pca(x6G08_MD_189_raw,pcadata)
points(x6G08_MD_189_proj[1],x6G08_MD_189_proj[2],pch=20)
text(x6G08_MD_189_proj[1],x6G08_MD_189_proj[2],pos=4,label="189",col="blue")

x6G08_MD_190_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_190.txt")
x6G08_MD_190_proj<-project.pca(x6G08_MD_190_raw,pcadata)
points(x6G08_MD_190_proj[1],x6G08_MD_190_proj[2],pch=20)
text(x6G08_MD_190_proj[1],x6G08_MD_190_proj[2],pos=4,label="190",col="blue")

x6G08_MD_191_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_191.txt")
x6G08_MD_191_proj<-project.pca(x6G08_MD_191_raw,pcadata)
points(x6G08_MD_191_proj[1],x6G08_MD_191_proj[2],pch=20)
text(x6G08_MD_191_proj[1],x6G08_MD_191_proj[2],pos=4,label="191",col="black")

x6G08_MD_192_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_192.txt")
x6G08_MD_192_proj<-project.pca(x6G08_MD_192_raw,pcadata)
points(x6G08_MD_192_proj[1],x6G08_MD_192_proj[2],pch=20)
text(x6G08_MD_192_proj[1],x6G08_MD_192_proj[2],pos=4,label="192",col="black")

x6G08_MD_193_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_193.txt")
x6G08_MD_193_proj<-project.pca(x6G08_MD_193_raw,pcadata)
points(x6G08_MD_193_proj[1],x6G08_MD_193_proj[2],pch=20)
text(x6G08_MD_193_proj[1],x6G08_MD_193_proj[2],pos=4,label="193",col="blue")

x6G08_MD_194_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_194.txt")
x6G08_MD_194_proj<-project.pca(x6G08_MD_194_raw,pcadata)
points(x6G08_MD_194_proj[1],x6G08_MD_194_proj[2],pch=20)
text(x6G08_MD_194_proj[1],x6G08_MD_194_proj[2],pos=4,label="194",col="blue")

x6G08_MD_195_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_195.txt")
x6G08_MD_195_proj<-project.pca(x6G08_MD_195_raw,pcadata)
points(x6G08_MD_195_proj[1],x6G08_MD_195_proj[2],pch=20)
text(x6G08_MD_195_proj[1],x6G08_MD_195_proj[2],pos=4,label="195",col="blue")

x6G08_MD_196_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_196.txt")
x6G08_MD_196_proj<-project.pca(x6G08_MD_196_raw,pcadata)
points(x6G08_MD_196_proj[1],x6G08_MD_196_proj[2],pch=20)
text(x6G08_MD_196_proj[1],x6G08_MD_196_proj[2],pos=4,label="196",col="black")

x6G08_MD_197_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_197.txt")
x6G08_MD_197_proj<-project.pca(x6G08_MD_197_raw,pcadata)
points(x6G08_MD_197_proj[1],x6G08_MD_197_proj[2],pch=20)
text(x6G08_MD_197_proj[1],x6G08_MD_197_proj[2],pos=4,label="197",col="black")

x6G08_MD_198_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_198.txt")
x6G08_MD_198_proj<-project.pca(x6G08_MD_198_raw,pcadata)
points(x6G08_MD_198_proj[1],x6G08_MD_198_proj[2],pch=20)
text(x6G08_MD_198_proj[1],x6G08_MD_198_proj[2],pos=4,label="198",col="black")

x6G08_MD_199_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_199.txt")
x6G08_MD_199_proj<-project.pca(x6G08_MD_199_raw,pcadata)
points(x6G08_MD_199_proj[1],x6G08_MD_199_proj[2],pch=20)
text(x6G08_MD_199_proj[1],x6G08_MD_199_proj[2],pos=4,label="199",col="black")

x6G08_MD_200_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_200.txt")
x6G08_MD_200_proj<-project.pca(x6G08_MD_200_raw,pcadata)
points(x6G08_MD_200_proj[1],x6G08_MD_200_proj[2],pch=20)
text(x6G08_MD_200_proj[1],x6G08_MD_200_proj[2],pos=4,label="200",col="black")

x6G08_MD_201_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_201.txt")
x6G08_MD_201_proj<-project.pca(x6G08_MD_201_raw,pcadata)
points(x6G08_MD_201_proj[1],x6G08_MD_201_proj[2],pch=20)
text(x6G08_MD_201_proj[1],x6G08_MD_201_proj[2],pos=4,label="201",col="blue")

x6G08_MD_202_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_202.txt")
x6G08_MD_202_proj<-project.pca(x6G08_MD_202_raw,pcadata)
points(x6G08_MD_202_proj[1],x6G08_MD_202_proj[2],pch=20)
text(x6G08_MD_202_proj[1],x6G08_MD_202_proj[2],pos=4,label="202",col="black")

x6G08_MD_203_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_203.txt")
x6G08_MD_203_proj<-project.pca(x6G08_MD_203_raw,pcadata)
points(x6G08_MD_203_proj[1],x6G08_MD_203_proj[2],pch=20)
text(x6G08_MD_203_proj[1],x6G08_MD_203_proj[2],pos=4,label="203",col="blue")

x6G08_MD_204_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_204.txt")
x6G08_MD_204_proj<-project.pca(x6G08_MD_204_raw,pcadata)
points(x6G08_MD_204_proj[1],x6G08_MD_204_proj[2],pch=20)
text(x6G08_MD_204_proj[1],x6G08_MD_204_proj[2],pos=4,label="204",col="blue")

x6G08_MD_205_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_205.txt")
x6G08_MD_205_proj<-project.pca(x6G08_MD_205_raw,pcadata)
points(x6G08_MD_205_proj[1],x6G08_MD_205_proj[2],pch=20)
text(x6G08_MD_205_proj[1],x6G08_MD_205_proj[2],pos=4,label="205",col="blue")

x6G08_MD_206_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_206.txt")
x6G08_MD_206_proj<-project.pca(x6G08_MD_206_raw,pcadata)
points(x6G08_MD_206_proj[1],x6G08_MD_206_proj[2],pch=20)
text(x6G08_MD_206_proj[1],x6G08_MD_206_proj[2],pos=4,label="206",col="red")

x6G08_MD_207_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_207.txt")
x6G08_MD_207_proj<-project.pca(x6G08_MD_207_raw,pcadata)
points(x6G08_MD_207_proj[1],x6G08_MD_207_proj[2],pch=20)
text(x6G08_MD_207_proj[1],x6G08_MD_207_proj[2],pos=4,label="207",col="blue")

x6G08_MD_208_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_208.txt")
x6G08_MD_208_proj<-project.pca(x6G08_MD_208_raw,pcadata)
points(x6G08_MD_208_proj[1],x6G08_MD_208_proj[2],pch=20)
text(x6G08_MD_208_proj[1],x6G08_MD_208_proj[2],pos=4,label="208",col="blue")

x6G08_MD_209_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_209.txt")
x6G08_MD_209_proj<-project.pca(x6G08_MD_209_raw,pcadata)
points(x6G08_MD_209_proj[1],x6G08_MD_209_proj[2],pch=20)
text(x6G08_MD_209_proj[1],x6G08_MD_209_proj[2],pos=4,label="209",col="black")

x6G08_MD_210_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_210.txt")
x6G08_MD_210_proj<-project.pca(x6G08_MD_210_raw,pcadata)
points(x6G08_MD_210_proj[1],x6G08_MD_210_proj[2],pch=20)
text(x6G08_MD_210_proj[1],x6G08_MD_210_proj[2],pos=4,label="210",col="blue")

x6G08_MD_211_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_211.txt")
x6G08_MD_211_proj<-project.pca(x6G08_MD_211_raw,pcadata)
points(x6G08_MD_211_proj[1],x6G08_MD_211_proj[2],pch=20)
text(x6G08_MD_211_proj[1],x6G08_MD_211_proj[2],pos=4,label="211",col="blue")

x6G08_MD_212_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_212.txt")
x6G08_MD_212_proj<-project.pca(x6G08_MD_212_raw,pcadata)
points(x6G08_MD_212_proj[1],x6G08_MD_212_proj[2],pch=20)
text(x6G08_MD_212_proj[1],x6G08_MD_212_proj[2],pos=4,label="212",col="blue")

x6G08_MD_213_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_213.txt")
x6G08_MD_213_proj<-project.pca(x6G08_MD_213_raw,pcadata)
points(x6G08_MD_213_proj[1],x6G08_MD_213_proj[2],pch=20)
text(x6G08_MD_213_proj[1],x6G08_MD_213_proj[2],pos=4,label="213",col="blue")

x6G08_MD_214_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_214.txt")
x6G08_MD_214_proj<-project.pca(x6G08_MD_214_raw,pcadata)
points(x6G08_MD_214_proj[1],x6G08_MD_214_proj[2],pch=20)
text(x6G08_MD_214_proj[1],x6G08_MD_214_proj[2],pos=4,label="214",col="blue")

x6G08_MD_215_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_215.txt")
x6G08_MD_215_proj<-project.pca(x6G08_MD_215_raw,pcadata)
points(x6G08_MD_215_proj[1],x6G08_MD_215_proj[2],pch=20)
text(x6G08_MD_215_proj[1],x6G08_MD_215_proj[2],pos=4,label="215",col="blue")

x6G08_MD_216_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_216.txt")
x6G08_MD_216_proj<-project.pca(x6G08_MD_216_raw,pcadata)
points(x6G08_MD_216_proj[1],x6G08_MD_216_proj[2],pch=20)
text(x6G08_MD_216_proj[1],x6G08_MD_216_proj[2],pos=4,label="216",col="blue")

x6G08_MD_217_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_217.txt")
x6G08_MD_217_proj<-project.pca(x6G08_MD_217_raw,pcadata)
points(x6G08_MD_217_proj[1],x6G08_MD_217_proj[2],pch=20)
text(x6G08_MD_217_proj[1],x6G08_MD_217_proj[2],pos=4,label="217",col="blue")

x6G08_MD_218_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_218.txt")
x6G08_MD_218_proj<-project.pca(x6G08_MD_218_raw,pcadata)
points(x6G08_MD_218_proj[1],x6G08_MD_218_proj[2],pch=20)
text(x6G08_MD_218_proj[1],x6G08_MD_218_proj[2],pos=4,label="218",col="green")

x6G08_MD_219_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_219.txt")
x6G08_MD_219_proj<-project.pca(x6G08_MD_219_raw,pcadata)
points(x6G08_MD_219_proj[1],x6G08_MD_219_proj[2],pch=20)
text(x6G08_MD_219_proj[1],x6G08_MD_219_proj[2],pos=4,label="219",col="red")

x6G08_MD_220_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_220.txt")
x6G08_MD_220_proj<-project.pca(x6G08_MD_220_raw,pcadata)
points(x6G08_MD_220_proj[1],x6G08_MD_220_proj[2],pch=20)
text(x6G08_MD_220_proj[1],x6G08_MD_220_proj[2],pos=4,label="220",col="red")

x6G08_MD_221_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_221.txt")
x6G08_MD_221_proj<-project.pca(x6G08_MD_221_raw,pcadata)
points(x6G08_MD_221_proj[1],x6G08_MD_221_proj[2],pch=20)
text(x6G08_MD_221_proj[1],x6G08_MD_221_proj[2],pos=4,label="221",col="blue")

x6G08_MD_222_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_222.txt")
x6G08_MD_222_proj<-project.pca(x6G08_MD_222_raw,pcadata)
points(x6G08_MD_222_proj[1],x6G08_MD_222_proj[2],pch=20)
text(x6G08_MD_222_proj[1],x6G08_MD_222_proj[2],pos=4,label="222",col="green")

x6G08_MD_223_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_223.txt")
x6G08_MD_223_proj<-project.pca(x6G08_MD_223_raw,pcadata)
points(x6G08_MD_223_proj[1],x6G08_MD_223_proj[2],pch=20)
text(x6G08_MD_223_proj[1],x6G08_MD_223_proj[2],pos=4,label="223",col="red")

x6G08_MD_224_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_224.txt")
x6G08_MD_224_proj<-project.pca(x6G08_MD_224_raw,pcadata)
points(x6G08_MD_224_proj[1],x6G08_MD_224_proj[2],pch=20)
text(x6G08_MD_224_proj[1],x6G08_MD_224_proj[2],pos=4,label="224",col="red")

x6G08_MD_225_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_225.txt")
x6G08_MD_225_proj<-project.pca(x6G08_MD_225_raw,pcadata)
points(x6G08_MD_225_proj[1],x6G08_MD_225_proj[2],pch=20)
text(x6G08_MD_225_proj[1],x6G08_MD_225_proj[2],pos=4,label="225",col="red")

x6G08_MD_226_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_226.txt")
x6G08_MD_226_proj<-project.pca(x6G08_MD_226_raw,pcadata)
points(x6G08_MD_226_proj[1],x6G08_MD_226_proj[2],pch=20)
text(x6G08_MD_226_proj[1],x6G08_MD_226_proj[2],pos=4,label="226",col="blue")

x6G08_MD_227_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_227.txt")
x6G08_MD_227_proj<-project.pca(x6G08_MD_227_raw,pcadata)
points(x6G08_MD_227_proj[1],x6G08_MD_227_proj[2],pch=20)
text(x6G08_MD_227_proj[1],x6G08_MD_227_proj[2],pos=4,label="227",col="red")

x6G08_MD_228_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_228.txt")
x6G08_MD_228_proj<-project.pca(x6G08_MD_228_raw,pcadata)
points(x6G08_MD_228_proj[1],x6G08_MD_228_proj[2],pch=20)
text(x6G08_MD_228_proj[1],x6G08_MD_228_proj[2],pos=4,label="228",col="green")

x6G08_MD_229_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_229.txt")
x6G08_MD_229_proj<-project.pca(x6G08_MD_229_raw,pcadata)
points(x6G08_MD_229_proj[1],x6G08_MD_229_proj[2],pch=20)
text(x6G08_MD_229_proj[1],x6G08_MD_229_proj[2],pos=4,label="229",col="red")

x6G08_MD_230_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_230.txt")
x6G08_MD_230_proj<-project.pca(x6G08_MD_230_raw,pcadata)
points(x6G08_MD_230_proj[1],x6G08_MD_230_proj[2],pch=20)
text(x6G08_MD_230_proj[1],x6G08_MD_230_proj[2],pos=4,label="230",col="green")

x6G08_MD_231_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_231.txt")
x6G08_MD_231_proj<-project.pca(x6G08_MD_231_raw,pcadata)
points(x6G08_MD_231_proj[1],x6G08_MD_231_proj[2],pch=20)
text(x6G08_MD_231_proj[1],x6G08_MD_231_proj[2],pos=4,label="231",col="green")

x6G08_MD_232_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_232.txt")
x6G08_MD_232_proj<-project.pca(x6G08_MD_232_raw,pcadata)
points(x6G08_MD_232_proj[1],x6G08_MD_232_proj[2],pch=20)
text(x6G08_MD_232_proj[1],x6G08_MD_232_proj[2],pos=4,label="232",col="red")

x6G08_MD_233_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_233.txt")
x6G08_MD_233_proj<-project.pca(x6G08_MD_233_raw,pcadata)
points(x6G08_MD_233_proj[1],x6G08_MD_233_proj[2],pch=20)
text(x6G08_MD_233_proj[1],x6G08_MD_233_proj[2],pos=4,label="233",col="blue")

x6G08_MD_234_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_234.txt")
x6G08_MD_234_proj<-project.pca(x6G08_MD_234_raw,pcadata)
points(x6G08_MD_234_proj[1],x6G08_MD_234_proj[2],pch=20)
text(x6G08_MD_234_proj[1],x6G08_MD_234_proj[2],pos=4,label="234",col="green")

x6G08_MD_235_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_235.txt")
x6G08_MD_235_proj<-project.pca(x6G08_MD_235_raw,pcadata)
points(x6G08_MD_235_proj[1],x6G08_MD_235_proj[2],pch=20)
text(x6G08_MD_235_proj[1],x6G08_MD_235_proj[2],pos=4,label="235",col="green")

x6G08_MD_236_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_236.txt")
x6G08_MD_236_proj<-project.pca(x6G08_MD_236_raw,pcadata)
points(x6G08_MD_236_proj[1],x6G08_MD_236_proj[2],pch=20)
text(x6G08_MD_236_proj[1],x6G08_MD_236_proj[2],pos=4,label="236",col="blue")

x6G08_MD_237_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_237.txt")
x6G08_MD_237_proj<-project.pca(x6G08_MD_237_raw,pcadata)
points(x6G08_MD_237_proj[1],x6G08_MD_237_proj[2],pch=20)
text(x6G08_MD_237_proj[1],x6G08_MD_237_proj[2],pos=4,label="237",col="green")

x6G08_MD_238_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_238.txt")
x6G08_MD_238_proj<-project.pca(x6G08_MD_238_raw,pcadata)
points(x6G08_MD_238_proj[1],x6G08_MD_238_proj[2],pch=20)
text(x6G08_MD_238_proj[1],x6G08_MD_238_proj[2],pos=4,label="238",col="green")

x6G08_MD_239_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_239.txt")
x6G08_MD_239_proj<-project.pca(x6G08_MD_239_raw,pcadata)
points(x6G08_MD_239_proj[1],x6G08_MD_239_proj[2],pch=20)
text(x6G08_MD_239_proj[1],x6G08_MD_239_proj[2],pos=4,label="239",col="red")

x6G08_MD_240_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_240.txt")
x6G08_MD_240_proj<-project.pca(x6G08_MD_240_raw,pcadata)
points(x6G08_MD_240_proj[1],x6G08_MD_240_proj[2],pch=20)
text(x6G08_MD_240_proj[1],x6G08_MD_240_proj[2],pos=4,label="240",col="green")

x6G08_MD_241_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_241.txt")
x6G08_MD_241_proj<-project.pca(x6G08_MD_241_raw,pcadata)
points(x6G08_MD_241_proj[1],x6G08_MD_241_proj[2],pch=20)
text(x6G08_MD_241_proj[1],x6G08_MD_241_proj[2],pos=4,label="241",col="green")

x6G08_MD_242_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_242.txt")
x6G08_MD_242_proj<-project.pca(x6G08_MD_242_raw,pcadata)
points(x6G08_MD_242_proj[1],x6G08_MD_242_proj[2],pch=20)
text(x6G08_MD_242_proj[1],x6G08_MD_242_proj[2],pos=4,label="242",col="green")

x6G08_MD_243_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_243.txt")
x6G08_MD_243_proj<-project.pca(x6G08_MD_243_raw,pcadata)
points(x6G08_MD_243_proj[1],x6G08_MD_243_proj[2],pch=20)
text(x6G08_MD_243_proj[1],x6G08_MD_243_proj[2],pos=4,label="243",col="green")

x6G08_MD_244_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_244.txt")
x6G08_MD_244_proj<-project.pca(x6G08_MD_244_raw,pcadata)
points(x6G08_MD_244_proj[1],x6G08_MD_244_proj[2],pch=20)
text(x6G08_MD_244_proj[1],x6G08_MD_244_proj[2],pos=4,label="244",col="blue")

x6G08_MD_245_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_245.txt")
x6G08_MD_245_proj<-project.pca(x6G08_MD_245_raw,pcadata)
points(x6G08_MD_245_proj[1],x6G08_MD_245_proj[2],pch=20)
text(x6G08_MD_245_proj[1],x6G08_MD_245_proj[2],pos=4,label="245",col="green")

x6G08_MD_246_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_246.txt")
x6G08_MD_246_proj<-project.pca(x6G08_MD_246_raw,pcadata)
points(x6G08_MD_246_proj[1],x6G08_MD_246_proj[2],pch=20)
text(x6G08_MD_246_proj[1],x6G08_MD_246_proj[2],pos=4,label="246",col="green")

x6G08_MD_247_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_247.txt")
x6G08_MD_247_proj<-project.pca(x6G08_MD_247_raw,pcadata)
points(x6G08_MD_247_proj[1],x6G08_MD_247_proj[2],pch=20)
text(x6G08_MD_247_proj[1],x6G08_MD_247_proj[2],pos=4,label="247",col="blue")

x6G08_MD_248_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_248.txt")
x6G08_MD_248_proj<-project.pca(x6G08_MD_248_raw,pcadata)
points(x6G08_MD_248_proj[1],x6G08_MD_248_proj[2],pch=20)
text(x6G08_MD_248_proj[1],x6G08_MD_248_proj[2],pos=4,label="248",col="black")

x6G08_MD_249_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_249.txt")
x6G08_MD_249_proj<-project.pca(x6G08_MD_249_raw,pcadata)
points(x6G08_MD_249_proj[1],x6G08_MD_249_proj[2],pch=20)
text(x6G08_MD_249_proj[1],x6G08_MD_249_proj[2],pos=4,label="249",col="blue")

x6G08_MD_250_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_250.txt")
x6G08_MD_250_proj<-project.pca(x6G08_MD_250_raw,pcadata)
points(x6G08_MD_250_proj[1],x6G08_MD_250_proj[2],pch=20)
text(x6G08_MD_250_proj[1],x6G08_MD_250_proj[2],pos=4,label="250",col="black")

x6G08_MD_251_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_251.txt")
x6G08_MD_251_proj<-project.pca(x6G08_MD_251_raw,pcadata)
points(x6G08_MD_251_proj[1],x6G08_MD_251_proj[2],pch=20)
text(x6G08_MD_251_proj[1],x6G08_MD_251_proj[2],pos=4,label="251",col="black")

x6G08_MD_252_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_252.txt")
x6G08_MD_252_proj<-project.pca(x6G08_MD_252_raw,pcadata)
points(x6G08_MD_252_proj[1],x6G08_MD_252_proj[2],pch=20)
text(x6G08_MD_252_proj[1],x6G08_MD_252_proj[2],pos=4,label="252",col="blue")

x6G08_MD_253_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_253.txt")
x6G08_MD_253_proj<-project.pca(x6G08_MD_253_raw,pcadata)
points(x6G08_MD_253_proj[1],x6G08_MD_253_proj[2],pch=20)
text(x6G08_MD_253_proj[1],x6G08_MD_253_proj[2],pos=4,label="253",col="blue")

x6G08_MD_254_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_254.txt")
x6G08_MD_254_proj<-project.pca(x6G08_MD_254_raw,pcadata)
points(x6G08_MD_254_proj[1],x6G08_MD_254_proj[2],pch=20)
text(x6G08_MD_254_proj[1],x6G08_MD_254_proj[2],pos=4,label="254",col="black")

x6G08_MD_255_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_255.txt")
x6G08_MD_255_proj<-project.pca(x6G08_MD_255_raw,pcadata)
points(x6G08_MD_255_proj[1],x6G08_MD_255_proj[2],pch=20)
text(x6G08_MD_255_proj[1],x6G08_MD_255_proj[2],pos=4,label="255",col="blue")

x6G08_MD_256_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_256.txt")
x6G08_MD_256_proj<-project.pca(x6G08_MD_256_raw,pcadata)
points(x6G08_MD_256_proj[1],x6G08_MD_256_proj[2],pch=20)
text(x6G08_MD_256_proj[1],x6G08_MD_256_proj[2],pos=4,label="256",col="blue")

x6G08_MD_257_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_257.txt")
x6G08_MD_257_proj<-project.pca(x6G08_MD_257_raw,pcadata)
points(x6G08_MD_257_proj[1],x6G08_MD_257_proj[2],pch=20)
text(x6G08_MD_257_proj[1],x6G08_MD_257_proj[2],pos=4,label="257",col="green")

x6G08_MD_258_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_258.txt")
x6G08_MD_258_proj<-project.pca(x6G08_MD_258_raw,pcadata)
points(x6G08_MD_258_proj[1],x6G08_MD_258_proj[2],pch=20)
text(x6G08_MD_258_proj[1],x6G08_MD_258_proj[2],pos=4,label="258",col="green")

x6G08_MD_259_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_259.txt")
x6G08_MD_259_proj<-project.pca(x6G08_MD_259_raw,pcadata)
points(x6G08_MD_259_proj[1],x6G08_MD_259_proj[2],pch=20)
text(x6G08_MD_259_proj[1],x6G08_MD_259_proj[2],pos=4,label="259",col="green")

x6G08_MD_260_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_260.txt")
x6G08_MD_260_proj<-project.pca(x6G08_MD_260_raw,pcadata)
points(x6G08_MD_260_proj[1],x6G08_MD_260_proj[2],pch=20)
text(x6G08_MD_260_proj[1],x6G08_MD_260_proj[2],pos=4,label="260",col="blue")

x6G08_MD_261_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_261.txt")
x6G08_MD_261_proj<-project.pca(x6G08_MD_261_raw,pcadata)
points(x6G08_MD_261_proj[1],x6G08_MD_261_proj[2],pch=20)
text(x6G08_MD_261_proj[1],x6G08_MD_261_proj[2],pos=4,label="261",col="blue")

x6G08_MD_262_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_262.txt")
x6G08_MD_262_proj<-project.pca(x6G08_MD_262_raw,pcadata)
points(x6G08_MD_262_proj[1],x6G08_MD_262_proj[2],pch=20)
text(x6G08_MD_262_proj[1],x6G08_MD_262_proj[2],pos=4,label="262",col="blue")

x6G08_MD_263_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_263.txt")
x6G08_MD_263_proj<-project.pca(x6G08_MD_263_raw,pcadata)
points(x6G08_MD_263_proj[1],x6G08_MD_263_proj[2],pch=20)
text(x6G08_MD_263_proj[1],x6G08_MD_263_proj[2],pos=4,label="263",col="blue")

x6G08_MD_264_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_264.txt")
x6G08_MD_264_proj<-project.pca(x6G08_MD_264_raw,pcadata)
points(x6G08_MD_264_proj[1],x6G08_MD_264_proj[2],pch=20)
text(x6G08_MD_264_proj[1],x6G08_MD_264_proj[2],pos=4,label="264",col="blue")

x6G08_MD_265_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_265.txt")
x6G08_MD_265_proj<-project.pca(x6G08_MD_265_raw,pcadata)
points(x6G08_MD_265_proj[1],x6G08_MD_265_proj[2],pch=20)
text(x6G08_MD_265_proj[1],x6G08_MD_265_proj[2],pos=4,label="265",col="blue")

x6G08_MD_266_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_266.txt")
x6G08_MD_266_proj<-project.pca(x6G08_MD_266_raw,pcadata)
points(x6G08_MD_266_proj[1],x6G08_MD_266_proj[2],pch=20)
text(x6G08_MD_266_proj[1],x6G08_MD_266_proj[2],pos=4,label="266",col="blue")

x6G08_MD_267_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_267.txt")
x6G08_MD_267_proj<-project.pca(x6G08_MD_267_raw,pcadata)
points(x6G08_MD_267_proj[1],x6G08_MD_267_proj[2],pch=20)
text(x6G08_MD_267_proj[1],x6G08_MD_267_proj[2],pos=4,label="267",col="blue")

x6G08_MD_268_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_268.txt")
x6G08_MD_268_proj<-project.pca(x6G08_MD_268_raw,pcadata)
points(x6G08_MD_268_proj[1],x6G08_MD_268_proj[2],pch=20)
text(x6G08_MD_268_proj[1],x6G08_MD_268_proj[2],pos=4,label="268",col="blue")

x6G08_MD_269_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_269.txt")
x6G08_MD_269_proj<-project.pca(x6G08_MD_269_raw,pcadata)
points(x6G08_MD_269_proj[1],x6G08_MD_269_proj[2],pch=20)
text(x6G08_MD_269_proj[1],x6G08_MD_269_proj[2],pos=4,label="269",col="blue")

x6G08_MD_270_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_270.txt")
x6G08_MD_270_proj<-project.pca(x6G08_MD_270_raw,pcadata)
points(x6G08_MD_270_proj[1],x6G08_MD_270_proj[2],pch=20)
text(x6G08_MD_270_proj[1],x6G08_MD_270_proj[2],pos=4,label="270",col="black")

x6G08_MD_271_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_271.txt")
x6G08_MD_271_proj<-project.pca(x6G08_MD_271_raw,pcadata)
points(x6G08_MD_271_proj[1],x6G08_MD_271_proj[2],pch=20)
text(x6G08_MD_271_proj[1],x6G08_MD_271_proj[2],pos=4,label="271",col="green")

x6G08_MD_272_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_272.txt")
x6G08_MD_272_proj<-project.pca(x6G08_MD_272_raw,pcadata)
points(x6G08_MD_272_proj[1],x6G08_MD_272_proj[2],pch=20)
text(x6G08_MD_272_proj[1],x6G08_MD_272_proj[2],pos=4,label="272",col="blue")

x6G08_MD_273_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_273.txt")
x6G08_MD_273_proj<-project.pca(x6G08_MD_273_raw,pcadata)
points(x6G08_MD_273_proj[1],x6G08_MD_273_proj[2],pch=20)
text(x6G08_MD_273_proj[1],x6G08_MD_273_proj[2],pos=4,label="273",col="blue")

x6G08_MD_274_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_274.txt")
x6G08_MD_274_proj<-project.pca(x6G08_MD_274_raw,pcadata)
points(x6G08_MD_274_proj[1],x6G08_MD_274_proj[2],pch=20)
text(x6G08_MD_274_proj[1],x6G08_MD_274_proj[2],pos=4,label="274",col="blue")

x6G08_MD_275_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_275.txt")
x6G08_MD_275_proj<-project.pca(x6G08_MD_275_raw,pcadata)
points(x6G08_MD_275_proj[1],x6G08_MD_275_proj[2],pch=20)
text(x6G08_MD_275_proj[1],x6G08_MD_275_proj[2],pos=4,label="275",col="green")

x6G08_MD_276_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_276.txt")
x6G08_MD_276_proj<-project.pca(x6G08_MD_276_raw,pcadata)
points(x6G08_MD_276_proj[1],x6G08_MD_276_proj[2],pch=20)
text(x6G08_MD_276_proj[1],x6G08_MD_276_proj[2],pos=4,label="276",col="blue")

x6G08_MD_277_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_277.txt")
x6G08_MD_277_proj<-project.pca(x6G08_MD_277_raw,pcadata)
points(x6G08_MD_277_proj[1],x6G08_MD_277_proj[2],pch=20)
text(x6G08_MD_277_proj[1],x6G08_MD_277_proj[2],pos=4,label="277",col="green")

x6G08_MD_278_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_278.txt")
x6G08_MD_278_proj<-project.pca(x6G08_MD_278_raw,pcadata)
points(x6G08_MD_278_proj[1],x6G08_MD_278_proj[2],pch=20)
text(x6G08_MD_278_proj[1],x6G08_MD_278_proj[2],pos=4,label="278",col="red")

x6G08_MD_279_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_279.txt")
x6G08_MD_279_proj<-project.pca(x6G08_MD_279_raw,pcadata)
points(x6G08_MD_279_proj[1],x6G08_MD_279_proj[2],pch=20)
text(x6G08_MD_279_proj[1],x6G08_MD_279_proj[2],pos=4,label="279",col="green")

x6G08_MD_280_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_280.txt")
x6G08_MD_280_proj<-project.pca(x6G08_MD_280_raw,pcadata)
points(x6G08_MD_280_proj[1],x6G08_MD_280_proj[2],pch=20)
text(x6G08_MD_280_proj[1],x6G08_MD_280_proj[2],pos=4,label="280",col="blue")

x6G08_MD_281_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_281.txt")
x6G08_MD_281_proj<-project.pca(x6G08_MD_281_raw,pcadata)
points(x6G08_MD_281_proj[1],x6G08_MD_281_proj[2],pch=20)
text(x6G08_MD_281_proj[1],x6G08_MD_281_proj[2],pos=4,label="281",col="green")

x6G08_MD_282_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_282.txt")
x6G08_MD_282_proj<-project.pca(x6G08_MD_282_raw,pcadata)
points(x6G08_MD_282_proj[1],x6G08_MD_282_proj[2],pch=20)
text(x6G08_MD_282_proj[1],x6G08_MD_282_proj[2],pos=4,label="282",col="green")

x6G08_MD_283_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_283.txt")
x6G08_MD_283_proj<-project.pca(x6G08_MD_283_raw,pcadata)
points(x6G08_MD_283_proj[1],x6G08_MD_283_proj[2],pch=20)
text(x6G08_MD_283_proj[1],x6G08_MD_283_proj[2],pos=4,label="283",col="green")

x6G08_MD_284_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_284.txt")
x6G08_MD_284_proj<-project.pca(x6G08_MD_284_raw,pcadata)
points(x6G08_MD_284_proj[1],x6G08_MD_284_proj[2],pch=20)
text(x6G08_MD_284_proj[1],x6G08_MD_284_proj[2],pos=4,label="284",col="green")

x6G08_MD_285_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_285.txt")
x6G08_MD_285_proj<-project.pca(x6G08_MD_285_raw,pcadata)
points(x6G08_MD_285_proj[1],x6G08_MD_285_proj[2],pch=20)
text(x6G08_MD_285_proj[1],x6G08_MD_285_proj[2],pos=4,label="285",col="red")

x6G08_MD_286_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_286.txt")
x6G08_MD_286_proj<-project.pca(x6G08_MD_286_raw,pcadata)
points(x6G08_MD_286_proj[1],x6G08_MD_286_proj[2],pch=20)
text(x6G08_MD_286_proj[1],x6G08_MD_286_proj[2],pos=4,label="286",col="red")

x6G08_MD_287_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_287.txt")
x6G08_MD_287_proj<-project.pca(x6G08_MD_287_raw,pcadata)
points(x6G08_MD_287_proj[1],x6G08_MD_287_proj[2],pch=20)
text(x6G08_MD_287_proj[1],x6G08_MD_287_proj[2],pos=4,label="287",col="red")

x6G08_MD_288_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_288.txt")
x6G08_MD_288_proj<-project.pca(x6G08_MD_288_raw,pcadata)
points(x6G08_MD_288_proj[1],x6G08_MD_288_proj[2],pch=20)
text(x6G08_MD_288_proj[1],x6G08_MD_288_proj[2],pos=4,label="288",col="green")

x6G08_MD_289_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_289.txt")
x6G08_MD_289_proj<-project.pca(x6G08_MD_289_raw,pcadata)
points(x6G08_MD_289_proj[1],x6G08_MD_289_proj[2],pch=20)
text(x6G08_MD_289_proj[1],x6G08_MD_289_proj[2],pos=4,label="289",col="blue")

x6G08_MD_290_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_290.txt")
x6G08_MD_290_proj<-project.pca(x6G08_MD_290_raw,pcadata)
points(x6G08_MD_290_proj[1],x6G08_MD_290_proj[2],pch=20)
text(x6G08_MD_290_proj[1],x6G08_MD_290_proj[2],pos=4,label="290",col="black")

x6G08_MD_291_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_291.txt")
x6G08_MD_291_proj<-project.pca(x6G08_MD_291_raw,pcadata)
points(x6G08_MD_291_proj[1],x6G08_MD_291_proj[2],pch=20)
text(x6G08_MD_291_proj[1],x6G08_MD_291_proj[2],pos=4,label="291",col="blue")

x6G08_MD_292_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_292.txt")
x6G08_MD_292_proj<-project.pca(x6G08_MD_292_raw,pcadata)
points(x6G08_MD_292_proj[1],x6G08_MD_292_proj[2],pch=20)
text(x6G08_MD_292_proj[1],x6G08_MD_292_proj[2],pos=4,label="292",col="black")

x6G08_MD_293_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_293.txt")
x6G08_MD_293_proj<-project.pca(x6G08_MD_293_raw,pcadata)
points(x6G08_MD_293_proj[1],x6G08_MD_293_proj[2],pch=20)
text(x6G08_MD_293_proj[1],x6G08_MD_293_proj[2],pos=4,label="293",col="blue")

x6G08_MD_294_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_294.txt")
x6G08_MD_294_proj<-project.pca(x6G08_MD_294_raw,pcadata)
points(x6G08_MD_294_proj[1],x6G08_MD_294_proj[2],pch=20)
text(x6G08_MD_294_proj[1],x6G08_MD_294_proj[2],pos=4,label="294",col="red")

x6G08_MD_295_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_295.txt")
x6G08_MD_295_proj<-project.pca(x6G08_MD_295_raw,pcadata)
points(x6G08_MD_295_proj[1],x6G08_MD_295_proj[2],pch=20)
text(x6G08_MD_295_proj[1],x6G08_MD_295_proj[2],pos=4,label="295",col="green")

x6G08_MD_296_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_296.txt")
x6G08_MD_296_proj<-project.pca(x6G08_MD_296_raw,pcadata)
points(x6G08_MD_296_proj[1],x6G08_MD_296_proj[2],pch=20)
text(x6G08_MD_296_proj[1],x6G08_MD_296_proj[2],pos=4,label="296",col="blue")

x6G08_MD_297_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_297.txt")
x6G08_MD_297_proj<-project.pca(x6G08_MD_297_raw,pcadata)
points(x6G08_MD_297_proj[1],x6G08_MD_297_proj[2],pch=20)
text(x6G08_MD_297_proj[1],x6G08_MD_297_proj[2],pos=4,label="297",col="red")

x6G08_MD_298_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_298.txt")
x6G08_MD_298_proj<-project.pca(x6G08_MD_298_raw,pcadata)
points(x6G08_MD_298_proj[1],x6G08_MD_298_proj[2],pch=20)
text(x6G08_MD_298_proj[1],x6G08_MD_298_proj[2],pos=4,label="298",col="red")

x6G08_MD_299_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_299.txt")
x6G08_MD_299_proj<-project.pca(x6G08_MD_299_raw,pcadata)
points(x6G08_MD_299_proj[1],x6G08_MD_299_proj[2],pch=20)
text(x6G08_MD_299_proj[1],x6G08_MD_299_proj[2],pos=4,label="299",col="green")

x6G08_MD_300_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_300.txt")
x6G08_MD_300_proj<-project.pca(x6G08_MD_300_raw,pcadata)
points(x6G08_MD_300_proj[1],x6G08_MD_300_proj[2],pch=20)
text(x6G08_MD_300_proj[1],x6G08_MD_300_proj[2],pos=4,label="300",col="green")

x6G08_MD_301_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_301.txt")
x6G08_MD_301_proj<-project.pca(x6G08_MD_301_raw,pcadata)
points(x6G08_MD_301_proj[1],x6G08_MD_301_proj[2],pch=20)
text(x6G08_MD_301_proj[1],x6G08_MD_301_proj[2],pos=4,label="301",col="blue")

x6G08_MD_302_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_302.txt")
x6G08_MD_302_proj<-project.pca(x6G08_MD_302_raw,pcadata)
points(x6G08_MD_302_proj[1],x6G08_MD_302_proj[2],pch=20)
text(x6G08_MD_302_proj[1],x6G08_MD_302_proj[2],pos=4,label="302",col="green")

x6G08_MD_303_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_303.txt")
x6G08_MD_303_proj<-project.pca(x6G08_MD_303_raw,pcadata)
points(x6G08_MD_303_proj[1],x6G08_MD_303_proj[2],pch=20)
text(x6G08_MD_303_proj[1],x6G08_MD_303_proj[2],pos=4,label="303",col="green")

x6G08_MD_304_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_304.txt")
x6G08_MD_304_proj<-project.pca(x6G08_MD_304_raw,pcadata)
points(x6G08_MD_304_proj[1],x6G08_MD_304_proj[2],pch=20)
text(x6G08_MD_304_proj[1],x6G08_MD_304_proj[2],pos=4,label="304",col="green")

x6G08_MD_305_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_305.txt")
x6G08_MD_305_proj<-project.pca(x6G08_MD_305_raw,pcadata)
points(x6G08_MD_305_proj[1],x6G08_MD_305_proj[2],pch=20)
text(x6G08_MD_305_proj[1],x6G08_MD_305_proj[2],pos=4,label="305",col="red")

x6G08_MD_306_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_306.txt")
x6G08_MD_306_proj<-project.pca(x6G08_MD_306_raw,pcadata)
points(x6G08_MD_306_proj[1],x6G08_MD_306_proj[2],pch=20)
text(x6G08_MD_306_proj[1],x6G08_MD_306_proj[2],pos=4,label="306",col="green")

x6G08_MD_307_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_307.txt")
x6G08_MD_307_proj<-project.pca(x6G08_MD_307_raw,pcadata)
points(x6G08_MD_307_proj[1],x6G08_MD_307_proj[2],pch=20)
text(x6G08_MD_307_proj[1],x6G08_MD_307_proj[2],pos=4,label="307",col="green")

x6G08_MD_308_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_308.txt")
x6G08_MD_308_proj<-project.pca(x6G08_MD_308_raw,pcadata)
points(x6G08_MD_308_proj[1],x6G08_MD_308_proj[2],pch=20)
text(x6G08_MD_308_proj[1],x6G08_MD_308_proj[2],pos=4,label="308",col="red")

x6G08_MD_309_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_309.txt")
x6G08_MD_309_proj<-project.pca(x6G08_MD_309_raw,pcadata)
points(x6G08_MD_309_proj[1],x6G08_MD_309_proj[2],pch=20)
text(x6G08_MD_309_proj[1],x6G08_MD_309_proj[2],pos=4,label="309",col="green")

x6G08_MD_310_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_310.txt")
x6G08_MD_310_proj<-project.pca(x6G08_MD_310_raw,pcadata)
points(x6G08_MD_310_proj[1],x6G08_MD_310_proj[2],pch=20)
text(x6G08_MD_310_proj[1],x6G08_MD_310_proj[2],pos=4,label="310",col="green")

x6G08_MD_311_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_311.txt")
x6G08_MD_311_proj<-project.pca(x6G08_MD_311_raw,pcadata)
points(x6G08_MD_311_proj[1],x6G08_MD_311_proj[2],pch=20)
text(x6G08_MD_311_proj[1],x6G08_MD_311_proj[2],pos=4,label="311",col="blue")

x6G08_MD_312_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_312.txt")
x6G08_MD_312_proj<-project.pca(x6G08_MD_312_raw,pcadata)
points(x6G08_MD_312_proj[1],x6G08_MD_312_proj[2],pch=20)
text(x6G08_MD_312_proj[1],x6G08_MD_312_proj[2],pos=4,label="312",col="blue")

x6G08_MD_313_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_313.txt")
x6G08_MD_313_proj<-project.pca(x6G08_MD_313_raw,pcadata)
points(x6G08_MD_313_proj[1],x6G08_MD_313_proj[2],pch=20)
text(x6G08_MD_313_proj[1],x6G08_MD_313_proj[2],pos=4,label="313",col="blue")

x6G08_MD_314_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_314.txt")
x6G08_MD_314_proj<-project.pca(x6G08_MD_314_raw,pcadata)
points(x6G08_MD_314_proj[1],x6G08_MD_314_proj[2],pch=20)
text(x6G08_MD_314_proj[1],x6G08_MD_314_proj[2],pos=4,label="314",col="green")

x6G08_MD_315_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_315.txt")
x6G08_MD_315_proj<-project.pca(x6G08_MD_315_raw,pcadata)
points(x6G08_MD_315_proj[1],x6G08_MD_315_proj[2],pch=20)
text(x6G08_MD_315_proj[1],x6G08_MD_315_proj[2],pos=4,label="315",col="blue")

x6G08_MD_316_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_316.txt")
x6G08_MD_316_proj<-project.pca(x6G08_MD_316_raw,pcadata)
points(x6G08_MD_316_proj[1],x6G08_MD_316_proj[2],pch=20)
text(x6G08_MD_316_proj[1],x6G08_MD_316_proj[2],pos=4,label="316",col="blue")

x6G08_MD_317_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_317.txt")
x6G08_MD_317_proj<-project.pca(x6G08_MD_317_raw,pcadata)
points(x6G08_MD_317_proj[1],x6G08_MD_317_proj[2],pch=20)
text(x6G08_MD_317_proj[1],x6G08_MD_317_proj[2],pos=4,label="317",col="green")

x6G08_MD_318_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_318.txt")
x6G08_MD_318_proj<-project.pca(x6G08_MD_318_raw,pcadata)
points(x6G08_MD_318_proj[1],x6G08_MD_318_proj[2],pch=20)
text(x6G08_MD_318_proj[1],x6G08_MD_318_proj[2],pos=4,label="318",col="green")

x6G08_MD_319_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_319.txt")
x6G08_MD_319_proj<-project.pca(x6G08_MD_319_raw,pcadata)
points(x6G08_MD_319_proj[1],x6G08_MD_319_proj[2],pch=20)
text(x6G08_MD_319_proj[1],x6G08_MD_319_proj[2],pos=4,label="319",col="green")

x6G08_MD_320_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_320.txt")
x6G08_MD_320_proj<-project.pca(x6G08_MD_320_raw,pcadata)
points(x6G08_MD_320_proj[1],x6G08_MD_320_proj[2],pch=20)
text(x6G08_MD_320_proj[1],x6G08_MD_320_proj[2],pos=4,label="320",col="green")

x6G08_MD_321_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_321.txt")
x6G08_MD_321_proj<-project.pca(x6G08_MD_321_raw,pcadata)
points(x6G08_MD_321_proj[1],x6G08_MD_321_proj[2],pch=20)
text(x6G08_MD_321_proj[1],x6G08_MD_321_proj[2],pos=4,label="321",col="blue")

x6G08_MD_322_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_322.txt")
x6G08_MD_322_proj<-project.pca(x6G08_MD_322_raw,pcadata)
points(x6G08_MD_322_proj[1],x6G08_MD_322_proj[2],pch=20)
text(x6G08_MD_322_proj[1],x6G08_MD_322_proj[2],pos=4,label="322",col="green")

x6G08_MD_323_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_323.txt")
x6G08_MD_323_proj<-project.pca(x6G08_MD_323_raw,pcadata)
points(x6G08_MD_323_proj[1],x6G08_MD_323_proj[2],pch=20)
text(x6G08_MD_323_proj[1],x6G08_MD_323_proj[2],pos=4,label="323",col="red")

x6G08_MD_324_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_324.txt")
x6G08_MD_324_proj<-project.pca(x6G08_MD_324_raw,pcadata)
points(x6G08_MD_324_proj[1],x6G08_MD_324_proj[2],pch=20)
text(x6G08_MD_324_proj[1],x6G08_MD_324_proj[2],pos=4,label="324",col="green")

x6G08_MD_325_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_325.txt")
x6G08_MD_325_proj<-project.pca(x6G08_MD_325_raw,pcadata)
points(x6G08_MD_325_proj[1],x6G08_MD_325_proj[2],pch=20)
text(x6G08_MD_325_proj[1],x6G08_MD_325_proj[2],pos=4,label="325",col="red")

x6G08_MD_326_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_326.txt")
x6G08_MD_326_proj<-project.pca(x6G08_MD_326_raw,pcadata)
points(x6G08_MD_326_proj[1],x6G08_MD_326_proj[2],pch=20)
text(x6G08_MD_326_proj[1],x6G08_MD_326_proj[2],pos=4,label="326",col="red")

x6G08_MD_327_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_327.txt")
x6G08_MD_327_proj<-project.pca(x6G08_MD_327_raw,pcadata)
points(x6G08_MD_327_proj[1],x6G08_MD_327_proj[2],pch=20)
text(x6G08_MD_327_proj[1],x6G08_MD_327_proj[2],pos=4,label="327",col="blue")

x6G08_MD_328_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_328.txt")
x6G08_MD_328_proj<-project.pca(x6G08_MD_328_raw,pcadata)
points(x6G08_MD_328_proj[1],x6G08_MD_328_proj[2],pch=20)
text(x6G08_MD_328_proj[1],x6G08_MD_328_proj[2],pos=4,label="328",col="red")

x6G08_MD_329_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_329.txt")
x6G08_MD_329_proj<-project.pca(x6G08_MD_329_raw,pcadata)
points(x6G08_MD_329_proj[1],x6G08_MD_329_proj[2],pch=20)
text(x6G08_MD_329_proj[1],x6G08_MD_329_proj[2],pos=4,label="329",col="red")

x6G08_MD_330_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_330.txt")
x6G08_MD_330_proj<-project.pca(x6G08_MD_330_raw,pcadata)
points(x6G08_MD_330_proj[1],x6G08_MD_330_proj[2],pch=20)
text(x6G08_MD_330_proj[1],x6G08_MD_330_proj[2],pos=4,label="330",col="green")

x6G08_MD_331_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_331.txt")
x6G08_MD_331_proj<-project.pca(x6G08_MD_331_raw,pcadata)
points(x6G08_MD_331_proj[1],x6G08_MD_331_proj[2],pch=20)
text(x6G08_MD_331_proj[1],x6G08_MD_331_proj[2],pos=4,label="331",col="green")

x6G08_MD_332_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_332.txt")
x6G08_MD_332_proj<-project.pca(x6G08_MD_332_raw,pcadata)
points(x6G08_MD_332_proj[1],x6G08_MD_332_proj[2],pch=20)
text(x6G08_MD_332_proj[1],x6G08_MD_332_proj[2],pos=4,label="332",col="red")

x6G08_MD_333_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_333.txt")
x6G08_MD_333_proj<-project.pca(x6G08_MD_333_raw,pcadata)
points(x6G08_MD_333_proj[1],x6G08_MD_333_proj[2],pch=20)
text(x6G08_MD_333_proj[1],x6G08_MD_333_proj[2],pos=4,label="333",col="green")

x6G08_MD_334_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_334.txt")
x6G08_MD_334_proj<-project.pca(x6G08_MD_334_raw,pcadata)
points(x6G08_MD_334_proj[1],x6G08_MD_334_proj[2],pch=20)
text(x6G08_MD_334_proj[1],x6G08_MD_334_proj[2],pos=4,label="334",col="green")

x6G08_MD_335_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_335.txt")
x6G08_MD_335_proj<-project.pca(x6G08_MD_335_raw,pcadata)
points(x6G08_MD_335_proj[1],x6G08_MD_335_proj[2],pch=20)
text(x6G08_MD_335_proj[1],x6G08_MD_335_proj[2],pos=4,label="335",col="green")

x6G08_MD_336_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_336.txt")
x6G08_MD_336_proj<-project.pca(x6G08_MD_336_raw,pcadata)
points(x6G08_MD_336_proj[1],x6G08_MD_336_proj[2],pch=20)
text(x6G08_MD_336_proj[1],x6G08_MD_336_proj[2],pos=4,label="336",col="green")

x6G08_MD_337_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_337.txt")
x6G08_MD_337_proj<-project.pca(x6G08_MD_337_raw,pcadata)
points(x6G08_MD_337_proj[1],x6G08_MD_337_proj[2],pch=20)
text(x6G08_MD_337_proj[1],x6G08_MD_337_proj[2],pos=4,label="337",col="green")

x6G08_MD_338_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_338.txt")
x6G08_MD_338_proj<-project.pca(x6G08_MD_338_raw,pcadata)
points(x6G08_MD_338_proj[1],x6G08_MD_338_proj[2],pch=20)
text(x6G08_MD_338_proj[1],x6G08_MD_338_proj[2],pos=4,label="338",col="blue")

x6G08_MD_339_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_339.txt")
x6G08_MD_339_proj<-project.pca(x6G08_MD_339_raw,pcadata)
points(x6G08_MD_339_proj[1],x6G08_MD_339_proj[2],pch=20)
text(x6G08_MD_339_proj[1],x6G08_MD_339_proj[2],pos=4,label="339",col="black")

x6G08_MD_340_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_340.txt")
x6G08_MD_340_proj<-project.pca(x6G08_MD_340_raw,pcadata)
points(x6G08_MD_340_proj[1],x6G08_MD_340_proj[2],pch=20)
text(x6G08_MD_340_proj[1],x6G08_MD_340_proj[2],pos=4,label="340",col="blue")

x6G08_MD_341_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_341.txt")
x6G08_MD_341_proj<-project.pca(x6G08_MD_341_raw,pcadata)
points(x6G08_MD_341_proj[1],x6G08_MD_341_proj[2],pch=20)
text(x6G08_MD_341_proj[1],x6G08_MD_341_proj[2],pos=4,label="341",col="blue")

x6G08_MD_342_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_342.txt")
x6G08_MD_342_proj<-project.pca(x6G08_MD_342_raw,pcadata)
points(x6G08_MD_342_proj[1],x6G08_MD_342_proj[2],pch=20)
text(x6G08_MD_342_proj[1],x6G08_MD_342_proj[2],pos=4,label="342",col="green")

x6G08_MD_343_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_343.txt")
x6G08_MD_343_proj<-project.pca(x6G08_MD_343_raw,pcadata)
points(x6G08_MD_343_proj[1],x6G08_MD_343_proj[2],pch=20)
text(x6G08_MD_343_proj[1],x6G08_MD_343_proj[2],pos=4,label="343",col="black")

x6G08_MD_344_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_344.txt")
x6G08_MD_344_proj<-project.pca(x6G08_MD_344_raw,pcadata)
points(x6G08_MD_344_proj[1],x6G08_MD_344_proj[2],pch=20)
text(x6G08_MD_344_proj[1],x6G08_MD_344_proj[2],pos=4,label="344",col="blue")

x6G08_MD_345_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_345.txt")
x6G08_MD_345_proj<-project.pca(x6G08_MD_345_raw,pcadata)
points(x6G08_MD_345_proj[1],x6G08_MD_345_proj[2],pch=20)
text(x6G08_MD_345_proj[1],x6G08_MD_345_proj[2],pos=4,label="345",col="blue")

x6G08_MD_346_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_346.txt")
x6G08_MD_346_proj<-project.pca(x6G08_MD_346_raw,pcadata)
points(x6G08_MD_346_proj[1],x6G08_MD_346_proj[2],pch=20)
text(x6G08_MD_346_proj[1],x6G08_MD_346_proj[2],pos=4,label="346",col="black")

x6G08_MD_347_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_347.txt")
x6G08_MD_347_proj<-project.pca(x6G08_MD_347_raw,pcadata)
points(x6G08_MD_347_proj[1],x6G08_MD_347_proj[2],pch=20)
text(x6G08_MD_347_proj[1],x6G08_MD_347_proj[2],pos=4,label="347",col="blue")

x6G08_MD_348_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_348.txt")
x6G08_MD_348_proj<-project.pca(x6G08_MD_348_raw,pcadata)
points(x6G08_MD_348_proj[1],x6G08_MD_348_proj[2],pch=20)
text(x6G08_MD_348_proj[1],x6G08_MD_348_proj[2],pos=4,label="348",col="black")

x6G08_MD_349_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_349.txt")
x6G08_MD_349_proj<-project.pca(x6G08_MD_349_raw,pcadata)
points(x6G08_MD_349_proj[1],x6G08_MD_349_proj[2],pch=20)
text(x6G08_MD_349_proj[1],x6G08_MD_349_proj[2],pos=4,label="349",col="black")

x6G08_MD_350_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_350.txt")
x6G08_MD_350_proj<-project.pca(x6G08_MD_350_raw,pcadata)
points(x6G08_MD_350_proj[1],x6G08_MD_350_proj[2],pch=20)
text(x6G08_MD_350_proj[1],x6G08_MD_350_proj[2],pos=4,label="350",col="black")

x6G08_MD_351_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_351.txt")
x6G08_MD_351_proj<-project.pca(x6G08_MD_351_raw,pcadata)
points(x6G08_MD_351_proj[1],x6G08_MD_351_proj[2],pch=20)
text(x6G08_MD_351_proj[1],x6G08_MD_351_proj[2],pos=4,label="351",col="blue")

x6G08_MD_352_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_352.txt")
x6G08_MD_352_proj<-project.pca(x6G08_MD_352_raw,pcadata)
points(x6G08_MD_352_proj[1],x6G08_MD_352_proj[2],pch=20)
text(x6G08_MD_352_proj[1],x6G08_MD_352_proj[2],pos=4,label="352",col="black")

x6G08_MD_353_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_353.txt")
x6G08_MD_353_proj<-project.pca(x6G08_MD_353_raw,pcadata)
points(x6G08_MD_353_proj[1],x6G08_MD_353_proj[2],pch=20)
text(x6G08_MD_353_proj[1],x6G08_MD_353_proj[2],pos=4,label="353",col="blue")

x6G08_MD_354_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_354.txt")
x6G08_MD_354_proj<-project.pca(x6G08_MD_354_raw,pcadata)
points(x6G08_MD_354_proj[1],x6G08_MD_354_proj[2],pch=20)
text(x6G08_MD_354_proj[1],x6G08_MD_354_proj[2],pos=4,label="354",col="black")

x6G08_MD_355_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_355.txt")
x6G08_MD_355_proj<-project.pca(x6G08_MD_355_raw,pcadata)
points(x6G08_MD_355_proj[1],x6G08_MD_355_proj[2],pch=20)
text(x6G08_MD_355_proj[1],x6G08_MD_355_proj[2],pos=4,label="355",col="black")

x6G08_MD_356_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_356.txt")
x6G08_MD_356_proj<-project.pca(x6G08_MD_356_raw,pcadata)
points(x6G08_MD_356_proj[1],x6G08_MD_356_proj[2],pch=20)
text(x6G08_MD_356_proj[1],x6G08_MD_356_proj[2],pos=4,label="356",col="green")

x6G08_MD_357_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_357.txt")
x6G08_MD_357_proj<-project.pca(x6G08_MD_357_raw,pcadata)
points(x6G08_MD_357_proj[1],x6G08_MD_357_proj[2],pch=20)
text(x6G08_MD_357_proj[1],x6G08_MD_357_proj[2],pos=4,label="357",col="blue")

x6G08_MD_358_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_358.txt")
x6G08_MD_358_proj<-project.pca(x6G08_MD_358_raw,pcadata)
points(x6G08_MD_358_proj[1],x6G08_MD_358_proj[2],pch=20)
text(x6G08_MD_358_proj[1],x6G08_MD_358_proj[2],pos=4,label="358",col="black")

x6G08_MD_359_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_359.txt")
x6G08_MD_359_proj<-project.pca(x6G08_MD_359_raw,pcadata)
points(x6G08_MD_359_proj[1],x6G08_MD_359_proj[2],pch=20)
text(x6G08_MD_359_proj[1],x6G08_MD_359_proj[2],pos=4,label="359",col="black")

x6G08_MD_360_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_360.txt")
x6G08_MD_360_proj<-project.pca(x6G08_MD_360_raw,pcadata)
points(x6G08_MD_360_proj[1],x6G08_MD_360_proj[2],pch=20)
text(x6G08_MD_360_proj[1],x6G08_MD_360_proj[2],pos=4,label="360",col="blue")

x6G08_MD_361_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_361.txt")
x6G08_MD_361_proj<-project.pca(x6G08_MD_361_raw,pcadata)
points(x6G08_MD_361_proj[1],x6G08_MD_361_proj[2],pch=20)
text(x6G08_MD_361_proj[1],x6G08_MD_361_proj[2],pos=4,label="361",col="red")

x6G08_MD_362_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_362.txt")
x6G08_MD_362_proj<-project.pca(x6G08_MD_362_raw,pcadata)
points(x6G08_MD_362_proj[1],x6G08_MD_362_proj[2],pch=20)
text(x6G08_MD_362_proj[1],x6G08_MD_362_proj[2],pos=4,label="362",col="red")

x6G08_MD_363_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_363.txt")
x6G08_MD_363_proj<-project.pca(x6G08_MD_363_raw,pcadata)
points(x6G08_MD_363_proj[1],x6G08_MD_363_proj[2],pch=20)
text(x6G08_MD_363_proj[1],x6G08_MD_363_proj[2],pos=4,label="363",col="red")

x6G08_MD_364_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_364.txt")
x6G08_MD_364_proj<-project.pca(x6G08_MD_364_raw,pcadata)
points(x6G08_MD_364_proj[1],x6G08_MD_364_proj[2],pch=20)
text(x6G08_MD_364_proj[1],x6G08_MD_364_proj[2],pos=4,label="364",col="red")

x6G08_MD_365_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_365.txt")
x6G08_MD_365_proj<-project.pca(x6G08_MD_365_raw,pcadata)
points(x6G08_MD_365_proj[1],x6G08_MD_365_proj[2],pch=20)
text(x6G08_MD_365_proj[1],x6G08_MD_365_proj[2],pos=4,label="365",col="red")

x6G08_MD_366_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_366.txt")
x6G08_MD_366_proj<-project.pca(x6G08_MD_366_raw,pcadata)
points(x6G08_MD_366_proj[1],x6G08_MD_366_proj[2],pch=20)
text(x6G08_MD_366_proj[1],x6G08_MD_366_proj[2],pos=4,label="366",col="red")

x6G08_MD_367_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_367.txt")
x6G08_MD_367_proj<-project.pca(x6G08_MD_367_raw,pcadata)
points(x6G08_MD_367_proj[1],x6G08_MD_367_proj[2],pch=20)
text(x6G08_MD_367_proj[1],x6G08_MD_367_proj[2],pos=4,label="367",col="green")

x6G08_MD_368_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_368.txt")
x6G08_MD_368_proj<-project.pca(x6G08_MD_368_raw,pcadata)
points(x6G08_MD_368_proj[1],x6G08_MD_368_proj[2],pch=20)
text(x6G08_MD_368_proj[1],x6G08_MD_368_proj[2],pos=4,label="368",col="green")

x6G08_MD_369_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_369.txt")
x6G08_MD_369_proj<-project.pca(x6G08_MD_369_raw,pcadata)
points(x6G08_MD_369_proj[1],x6G08_MD_369_proj[2],pch=20)
text(x6G08_MD_369_proj[1],x6G08_MD_369_proj[2],pos=4,label="369",col="red")

x6G08_MD_370_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_370.txt")
x6G08_MD_370_proj<-project.pca(x6G08_MD_370_raw,pcadata)
points(x6G08_MD_370_proj[1],x6G08_MD_370_proj[2],pch=20)
text(x6G08_MD_370_proj[1],x6G08_MD_370_proj[2],pos=4,label="370",col="green")

x6G08_MD_371_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_371.txt")
x6G08_MD_371_proj<-project.pca(x6G08_MD_371_raw,pcadata)
points(x6G08_MD_371_proj[1],x6G08_MD_371_proj[2],pch=20)
text(x6G08_MD_371_proj[1],x6G08_MD_371_proj[2],pos=4,label="371",col="blue")

x6G08_MD_372_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_372.txt")
x6G08_MD_372_proj<-project.pca(x6G08_MD_372_raw,pcadata)
points(x6G08_MD_372_proj[1],x6G08_MD_372_proj[2],pch=20)
text(x6G08_MD_372_proj[1],x6G08_MD_372_proj[2],pos=4,label="372",col="blue")

x6G08_MD_373_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_373.txt")
x6G08_MD_373_proj<-project.pca(x6G08_MD_373_raw,pcadata)
points(x6G08_MD_373_proj[1],x6G08_MD_373_proj[2],pch=20)
text(x6G08_MD_373_proj[1],x6G08_MD_373_proj[2],pos=4,label="373",col="blue")

x6G08_MD_374_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_374.txt")
x6G08_MD_374_proj<-project.pca(x6G08_MD_374_raw,pcadata)
points(x6G08_MD_374_proj[1],x6G08_MD_374_proj[2],pch=20)
text(x6G08_MD_374_proj[1],x6G08_MD_374_proj[2],pos=4,label="374",col="blue")

x6G08_MD_375_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_375.txt")
x6G08_MD_375_proj<-project.pca(x6G08_MD_375_raw,pcadata)
points(x6G08_MD_375_proj[1],x6G08_MD_375_proj[2],pch=20)
text(x6G08_MD_375_proj[1],x6G08_MD_375_proj[2],pos=4,label="375",col="red")

x6G08_MD_376_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_376.txt")
x6G08_MD_376_proj<-project.pca(x6G08_MD_376_raw,pcadata)
points(x6G08_MD_376_proj[1],x6G08_MD_376_proj[2],pch=20)
text(x6G08_MD_376_proj[1],x6G08_MD_376_proj[2],pos=4,label="376",col="blue")

x6G08_MD_377_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_377.txt")
x6G08_MD_377_proj<-project.pca(x6G08_MD_377_raw,pcadata)
points(x6G08_MD_377_proj[1],x6G08_MD_377_proj[2],pch=20)
text(x6G08_MD_377_proj[1],x6G08_MD_377_proj[2],pos=4,label="377",col="blue")

x6G08_MD_378_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_378.txt")
x6G08_MD_378_proj<-project.pca(x6G08_MD_378_raw,pcadata)
points(x6G08_MD_378_proj[1],x6G08_MD_378_proj[2],pch=20)
text(x6G08_MD_378_proj[1],x6G08_MD_378_proj[2],pos=4,label="378",col="black")

x6G08_MD_379_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_379.txt")
x6G08_MD_379_proj<-project.pca(x6G08_MD_379_raw,pcadata)
points(x6G08_MD_379_proj[1],x6G08_MD_379_proj[2],pch=20)
text(x6G08_MD_379_proj[1],x6G08_MD_379_proj[2],pos=4,label="379",col="blue")

x6G08_MD_380_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_380.txt")
x6G08_MD_380_proj<-project.pca(x6G08_MD_380_raw,pcadata)
points(x6G08_MD_380_proj[1],x6G08_MD_380_proj[2],pch=20)
text(x6G08_MD_380_proj[1],x6G08_MD_380_proj[2],pos=4,label="380",col="blue")

x6G08_MD_381_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_381.txt")
x6G08_MD_381_proj<-project.pca(x6G08_MD_381_raw,pcadata)
points(x6G08_MD_381_proj[1],x6G08_MD_381_proj[2],pch=20)
text(x6G08_MD_381_proj[1],x6G08_MD_381_proj[2],pos=4,label="381",col="black")

x6G08_MD_382_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_382.txt")
x6G08_MD_382_proj<-project.pca(x6G08_MD_382_raw,pcadata)
points(x6G08_MD_382_proj[1],x6G08_MD_382_proj[2],pch=20)
text(x6G08_MD_382_proj[1],x6G08_MD_382_proj[2],pos=4,label="382",col="blue")

x6G08_MD_383_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_383.txt")
x6G08_MD_383_proj<-project.pca(x6G08_MD_383_raw,pcadata)
points(x6G08_MD_383_proj[1],x6G08_MD_383_proj[2],pch=20)
text(x6G08_MD_383_proj[1],x6G08_MD_383_proj[2],pos=4,label="383",col="blue")

x6G08_MD_384_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_384.txt")
x6G08_MD_384_proj<-project.pca(x6G08_MD_384_raw,pcadata)
points(x6G08_MD_384_proj[1],x6G08_MD_384_proj[2],pch=20)
text(x6G08_MD_384_proj[1],x6G08_MD_384_proj[2],pos=4,label="384",col="black")

x6G08_MD_385_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_385.txt")
x6G08_MD_385_proj<-project.pca(x6G08_MD_385_raw,pcadata)
points(x6G08_MD_385_proj[1],x6G08_MD_385_proj[2],pch=20)
text(x6G08_MD_385_proj[1],x6G08_MD_385_proj[2],pos=4,label="385",col="blue")

x6G08_MD_386_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_386.txt")
x6G08_MD_386_proj<-project.pca(x6G08_MD_386_raw,pcadata)
points(x6G08_MD_386_proj[1],x6G08_MD_386_proj[2],pch=20)
text(x6G08_MD_386_proj[1],x6G08_MD_386_proj[2],pos=4,label="386",col="blue")

x6G08_MD_387_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_387.txt")
x6G08_MD_387_proj<-project.pca(x6G08_MD_387_raw,pcadata)
points(x6G08_MD_387_proj[1],x6G08_MD_387_proj[2],pch=20)
text(x6G08_MD_387_proj[1],x6G08_MD_387_proj[2],pos=4,label="387",col="blue")

x6G08_MD_388_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_388.txt")
x6G08_MD_388_proj<-project.pca(x6G08_MD_388_raw,pcadata)
points(x6G08_MD_388_proj[1],x6G08_MD_388_proj[2],pch=20)
text(x6G08_MD_388_proj[1],x6G08_MD_388_proj[2],pos=4,label="388",col="red")

x6G08_MD_389_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_389.txt")
x6G08_MD_389_proj<-project.pca(x6G08_MD_389_raw,pcadata)
points(x6G08_MD_389_proj[1],x6G08_MD_389_proj[2],pch=20)
text(x6G08_MD_389_proj[1],x6G08_MD_389_proj[2],pos=4,label="389",col="green")

x6G08_MD_390_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_390.txt")
x6G08_MD_390_proj<-project.pca(x6G08_MD_390_raw,pcadata)
points(x6G08_MD_390_proj[1],x6G08_MD_390_proj[2],pch=20)
text(x6G08_MD_390_proj[1],x6G08_MD_390_proj[2],pos=4,label="390",col="blue")

x6G08_MD_391_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_391.txt")
x6G08_MD_391_proj<-project.pca(x6G08_MD_391_raw,pcadata)
points(x6G08_MD_391_proj[1],x6G08_MD_391_proj[2],pch=20)
text(x6G08_MD_391_proj[1],x6G08_MD_391_proj[2],pos=4,label="391",col="green")

x6G08_MD_392_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_392.txt")
x6G08_MD_392_proj<-project.pca(x6G08_MD_392_raw,pcadata)
points(x6G08_MD_392_proj[1],x6G08_MD_392_proj[2],pch=20)
text(x6G08_MD_392_proj[1],x6G08_MD_392_proj[2],pos=4,label="392",col="green")

x6G08_MD_393_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_393.txt")
x6G08_MD_393_proj<-project.pca(x6G08_MD_393_raw,pcadata)
points(x6G08_MD_393_proj[1],x6G08_MD_393_proj[2],pch=20)
text(x6G08_MD_393_proj[1],x6G08_MD_393_proj[2],pos=4,label="393",col="green")

x6G08_MD_394_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_394.txt")
x6G08_MD_394_proj<-project.pca(x6G08_MD_394_raw,pcadata)
points(x6G08_MD_394_proj[1],x6G08_MD_394_proj[2],pch=20)
text(x6G08_MD_394_proj[1],x6G08_MD_394_proj[2],pos=4,label="394",col="green")

x6G08_MD_395_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_395.txt")
x6G08_MD_395_proj<-project.pca(x6G08_MD_395_raw,pcadata)
points(x6G08_MD_395_proj[1],x6G08_MD_395_proj[2],pch=20)
text(x6G08_MD_395_proj[1],x6G08_MD_395_proj[2],pos=4,label="395",col="green")

x6G08_MD_396_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_396.txt")
x6G08_MD_396_proj<-project.pca(x6G08_MD_396_raw,pcadata)
points(x6G08_MD_396_proj[1],x6G08_MD_396_proj[2],pch=20)
text(x6G08_MD_396_proj[1],x6G08_MD_396_proj[2],pos=4,label="396",col="green")

x6G08_MD_397_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_397.txt")
x6G08_MD_397_proj<-project.pca(x6G08_MD_397_raw,pcadata)
points(x6G08_MD_397_proj[1],x6G08_MD_397_proj[2],pch=20)
text(x6G08_MD_397_proj[1],x6G08_MD_397_proj[2],pos=4,label="397",col="blue")

x6G08_MD_398_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_398.txt")
x6G08_MD_398_proj<-project.pca(x6G08_MD_398_raw,pcadata)
points(x6G08_MD_398_proj[1],x6G08_MD_398_proj[2],pch=20)
text(x6G08_MD_398_proj[1],x6G08_MD_398_proj[2],pos=4,label="398",col="blue")

x6G08_MD_399_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_399.txt")
x6G08_MD_399_proj<-project.pca(x6G08_MD_399_raw,pcadata)
points(x6G08_MD_399_proj[1],x6G08_MD_399_proj[2],pch=20)
text(x6G08_MD_399_proj[1],x6G08_MD_399_proj[2],pos=4,label="399",col="green")

x6G08_MD_400_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_400.txt")
x6G08_MD_400_proj<-project.pca(x6G08_MD_400_raw,pcadata)
points(x6G08_MD_400_proj[1],x6G08_MD_400_proj[2],pch=20)
text(x6G08_MD_400_proj[1],x6G08_MD_400_proj[2],pos=4,label="400",col="blue")

x6G08_MD_401_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_401.txt")
x6G08_MD_401_proj<-project.pca(x6G08_MD_401_raw,pcadata)
points(x6G08_MD_401_proj[1],x6G08_MD_401_proj[2],pch=20)
text(x6G08_MD_401_proj[1],x6G08_MD_401_proj[2],pos=4,label="401",col="green")

x6G08_MD_402_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_402.txt")
x6G08_MD_402_proj<-project.pca(x6G08_MD_402_raw,pcadata)
points(x6G08_MD_402_proj[1],x6G08_MD_402_proj[2],pch=20)
text(x6G08_MD_402_proj[1],x6G08_MD_402_proj[2],pos=4,label="402",col="green")

x6G08_MD_403_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_403.txt")
x6G08_MD_403_proj<-project.pca(x6G08_MD_403_raw,pcadata)
points(x6G08_MD_403_proj[1],x6G08_MD_403_proj[2],pch=20)
text(x6G08_MD_403_proj[1],x6G08_MD_403_proj[2],pos=4,label="403",col="red")

x6G08_MD_404_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_404.txt")
x6G08_MD_404_proj<-project.pca(x6G08_MD_404_raw,pcadata)
points(x6G08_MD_404_proj[1],x6G08_MD_404_proj[2],pch=20)
text(x6G08_MD_404_proj[1],x6G08_MD_404_proj[2],pos=4,label="404",col="blue")

x6G08_MD_405_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_405.txt")
x6G08_MD_405_proj<-project.pca(x6G08_MD_405_raw,pcadata)
points(x6G08_MD_405_proj[1],x6G08_MD_405_proj[2],pch=20)
text(x6G08_MD_405_proj[1],x6G08_MD_405_proj[2],pos=4,label="405",col="red")

x6G08_MD_406_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_406.txt")
x6G08_MD_406_proj<-project.pca(x6G08_MD_406_raw,pcadata)
points(x6G08_MD_406_proj[1],x6G08_MD_406_proj[2],pch=20)
text(x6G08_MD_406_proj[1],x6G08_MD_406_proj[2],pos=4,label="406",col="red")

x6G08_MD_407_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_407.txt")
x6G08_MD_407_proj<-project.pca(x6G08_MD_407_raw,pcadata)
points(x6G08_MD_407_proj[1],x6G08_MD_407_proj[2],pch=20)
text(x6G08_MD_407_proj[1],x6G08_MD_407_proj[2],pos=4,label="407",col="red")

x6G08_MD_408_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_408.txt")
x6G08_MD_408_proj<-project.pca(x6G08_MD_408_raw,pcadata)
points(x6G08_MD_408_proj[1],x6G08_MD_408_proj[2],pch=20)
text(x6G08_MD_408_proj[1],x6G08_MD_408_proj[2],pos=4,label="408",col="red")

x6G08_MD_409_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_409.txt")
x6G08_MD_409_proj<-project.pca(x6G08_MD_409_raw,pcadata)
points(x6G08_MD_409_proj[1],x6G08_MD_409_proj[2],pch=20)
text(x6G08_MD_409_proj[1],x6G08_MD_409_proj[2],pos=4,label="409",col="red")

x6G08_MD_410_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_410.txt")
x6G08_MD_410_proj<-project.pca(x6G08_MD_410_raw,pcadata)
points(x6G08_MD_410_proj[1],x6G08_MD_410_proj[2],pch=20)
text(x6G08_MD_410_proj[1],x6G08_MD_410_proj[2],pos=4,label="410",col="green")

x6G08_MD_411_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_411.txt")
x6G08_MD_411_proj<-project.pca(x6G08_MD_411_raw,pcadata)
points(x6G08_MD_411_proj[1],x6G08_MD_411_proj[2],pch=20)
text(x6G08_MD_411_proj[1],x6G08_MD_411_proj[2],pos=4,label="411",col="red")

x6G08_MD_412_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_412.txt")
x6G08_MD_412_proj<-project.pca(x6G08_MD_412_raw,pcadata)
points(x6G08_MD_412_proj[1],x6G08_MD_412_proj[2],pch=20)
text(x6G08_MD_412_proj[1],x6G08_MD_412_proj[2],pos=4,label="412",col="red")

x6G08_MD_413_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_413.txt")
x6G08_MD_413_proj<-project.pca(x6G08_MD_413_raw,pcadata)
points(x6G08_MD_413_proj[1],x6G08_MD_413_proj[2],pch=20)
text(x6G08_MD_413_proj[1],x6G08_MD_413_proj[2],pos=4,label="413",col="green")

x6G08_MD_414_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_414.txt")
x6G08_MD_414_proj<-project.pca(x6G08_MD_414_raw,pcadata)
points(x6G08_MD_414_proj[1],x6G08_MD_414_proj[2],pch=20)
text(x6G08_MD_414_proj[1],x6G08_MD_414_proj[2],pos=4,label="414",col="green")

x6G08_MD_415_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_415.txt")
x6G08_MD_415_proj<-project.pca(x6G08_MD_415_raw,pcadata)
points(x6G08_MD_415_proj[1],x6G08_MD_415_proj[2],pch=20)
text(x6G08_MD_415_proj[1],x6G08_MD_415_proj[2],pos=4,label="415",col="green")

x6G08_MD_416_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_416.txt")
x6G08_MD_416_proj<-project.pca(x6G08_MD_416_raw,pcadata)
points(x6G08_MD_416_proj[1],x6G08_MD_416_proj[2],pch=20)
text(x6G08_MD_416_proj[1],x6G08_MD_416_proj[2],pos=4,label="416",col="red")

x6G08_MD_417_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_417.txt")
x6G08_MD_417_proj<-project.pca(x6G08_MD_417_raw,pcadata)
points(x6G08_MD_417_proj[1],x6G08_MD_417_proj[2],pch=20)
text(x6G08_MD_417_proj[1],x6G08_MD_417_proj[2],pos=4,label="417",col="red")

x6G08_MD_418_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_418.txt")
x6G08_MD_418_proj<-project.pca(x6G08_MD_418_raw,pcadata)
points(x6G08_MD_418_proj[1],x6G08_MD_418_proj[2],pch=20)
text(x6G08_MD_418_proj[1],x6G08_MD_418_proj[2],pos=4,label="418",col="red")

x6G08_MD_419_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_419.txt")
x6G08_MD_419_proj<-project.pca(x6G08_MD_419_raw,pcadata)
points(x6G08_MD_419_proj[1],x6G08_MD_419_proj[2],pch=20)
text(x6G08_MD_419_proj[1],x6G08_MD_419_proj[2],pos=4,label="419",col="green")

x6G08_MD_420_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_420.txt")
x6G08_MD_420_proj<-project.pca(x6G08_MD_420_raw,pcadata)
points(x6G08_MD_420_proj[1],x6G08_MD_420_proj[2],pch=20)
text(x6G08_MD_420_proj[1],x6G08_MD_420_proj[2],pos=4,label="420",col="red")

x6G08_MD_421_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_421.txt")
x6G08_MD_421_proj<-project.pca(x6G08_MD_421_raw,pcadata)
points(x6G08_MD_421_proj[1],x6G08_MD_421_proj[2],pch=20)
text(x6G08_MD_421_proj[1],x6G08_MD_421_proj[2],pos=4,label="421",col="red")

x6G08_MD_422_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_422.txt")
x6G08_MD_422_proj<-project.pca(x6G08_MD_422_raw,pcadata)
points(x6G08_MD_422_proj[1],x6G08_MD_422_proj[2],pch=20)
text(x6G08_MD_422_proj[1],x6G08_MD_422_proj[2],pos=4,label="422",col="red")

x6G08_MD_423_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_423.txt")
x6G08_MD_423_proj<-project.pca(x6G08_MD_423_raw,pcadata)
points(x6G08_MD_423_proj[1],x6G08_MD_423_proj[2],pch=20)
text(x6G08_MD_423_proj[1],x6G08_MD_423_proj[2],pos=4,label="423",col="blue")

x6G08_MD_424_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_424.txt")
x6G08_MD_424_proj<-project.pca(x6G08_MD_424_raw,pcadata)
points(x6G08_MD_424_proj[1],x6G08_MD_424_proj[2],pch=20)
text(x6G08_MD_424_proj[1],x6G08_MD_424_proj[2],pos=4,label="424",col="red")

x6G08_MD_425_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_425.txt")
x6G08_MD_425_proj<-project.pca(x6G08_MD_425_raw,pcadata)
points(x6G08_MD_425_proj[1],x6G08_MD_425_proj[2],pch=20)
text(x6G08_MD_425_proj[1],x6G08_MD_425_proj[2],pos=4,label="425",col="red")

x6G08_MD_426_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_426.txt")
x6G08_MD_426_proj<-project.pca(x6G08_MD_426_raw,pcadata)
points(x6G08_MD_426_proj[1],x6G08_MD_426_proj[2],pch=20)
text(x6G08_MD_426_proj[1],x6G08_MD_426_proj[2],pos=4,label="426",col="red")

x6G08_MD_427_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_427.txt")
x6G08_MD_427_proj<-project.pca(x6G08_MD_427_raw,pcadata)
points(x6G08_MD_427_proj[1],x6G08_MD_427_proj[2],pch=20)
text(x6G08_MD_427_proj[1],x6G08_MD_427_proj[2],pos=4,label="427",col="green")

x6G08_MD_428_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_428.txt")
x6G08_MD_428_proj<-project.pca(x6G08_MD_428_raw,pcadata)
points(x6G08_MD_428_proj[1],x6G08_MD_428_proj[2],pch=20)
text(x6G08_MD_428_proj[1],x6G08_MD_428_proj[2],pos=4,label="428",col="red")

x6G08_MD_429_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_429.txt")
x6G08_MD_429_proj<-project.pca(x6G08_MD_429_raw,pcadata)
points(x6G08_MD_429_proj[1],x6G08_MD_429_proj[2],pch=20)
text(x6G08_MD_429_proj[1],x6G08_MD_429_proj[2],pos=4,label="429",col="green")

x6G08_MD_430_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_430.txt")
x6G08_MD_430_proj<-project.pca(x6G08_MD_430_raw,pcadata)
points(x6G08_MD_430_proj[1],x6G08_MD_430_proj[2],pch=20)
text(x6G08_MD_430_proj[1],x6G08_MD_430_proj[2],pos=4,label="430",col="red")

x6G08_MD_431_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_431.txt")
x6G08_MD_431_proj<-project.pca(x6G08_MD_431_raw,pcadata)
points(x6G08_MD_431_proj[1],x6G08_MD_431_proj[2],pch=20)
text(x6G08_MD_431_proj[1],x6G08_MD_431_proj[2],pos=4,label="431",col="green")

x6G08_MD_432_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_432.txt")
x6G08_MD_432_proj<-project.pca(x6G08_MD_432_raw,pcadata)
points(x6G08_MD_432_proj[1],x6G08_MD_432_proj[2],pch=20)
text(x6G08_MD_432_proj[1],x6G08_MD_432_proj[2],pos=4,label="432",col="green")

x6G08_MD_433_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_433.txt")
x6G08_MD_433_proj<-project.pca(x6G08_MD_433_raw,pcadata)
points(x6G08_MD_433_proj[1],x6G08_MD_433_proj[2],pch=20)
text(x6G08_MD_433_proj[1],x6G08_MD_433_proj[2],pos=4,label="433",col="green")

x6G08_MD_434_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_434.txt")
x6G08_MD_434_proj<-project.pca(x6G08_MD_434_raw,pcadata)
points(x6G08_MD_434_proj[1],x6G08_MD_434_proj[2],pch=20)
text(x6G08_MD_434_proj[1],x6G08_MD_434_proj[2],pos=4,label="434",col="red")

x6G08_MD_435_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_435.txt")
x6G08_MD_435_proj<-project.pca(x6G08_MD_435_raw,pcadata)
points(x6G08_MD_435_proj[1],x6G08_MD_435_proj[2],pch=20)
text(x6G08_MD_435_proj[1],x6G08_MD_435_proj[2],pos=4,label="435",col="green")

x6G08_MD_436_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_436.txt")
x6G08_MD_436_proj<-project.pca(x6G08_MD_436_raw,pcadata)
points(x6G08_MD_436_proj[1],x6G08_MD_436_proj[2],pch=20)
text(x6G08_MD_436_proj[1],x6G08_MD_436_proj[2],pos=4,label="436",col="red")

x6G08_MD_437_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_437.txt")
x6G08_MD_437_proj<-project.pca(x6G08_MD_437_raw,pcadata)
points(x6G08_MD_437_proj[1],x6G08_MD_437_proj[2],pch=20)
text(x6G08_MD_437_proj[1],x6G08_MD_437_proj[2],pos=4,label="437",col="green")

x6G08_MD_438_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_438.txt")
x6G08_MD_438_proj<-project.pca(x6G08_MD_438_raw,pcadata)
points(x6G08_MD_438_proj[1],x6G08_MD_438_proj[2],pch=20)
text(x6G08_MD_438_proj[1],x6G08_MD_438_proj[2],pos=4,label="438",col="red")

x6G08_MD_439_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_439.txt")
x6G08_MD_439_proj<-project.pca(x6G08_MD_439_raw,pcadata)
points(x6G08_MD_439_proj[1],x6G08_MD_439_proj[2],pch=20)
text(x6G08_MD_439_proj[1],x6G08_MD_439_proj[2],pos=4,label="439",col="red")

x6G08_MD_440_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_440.txt")
x6G08_MD_440_proj<-project.pca(x6G08_MD_440_raw,pcadata)
points(x6G08_MD_440_proj[1],x6G08_MD_440_proj[2],pch=20)
text(x6G08_MD_440_proj[1],x6G08_MD_440_proj[2],pos=4,label="440",col="red")

x6G08_MD_441_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_441.txt")
x6G08_MD_441_proj<-project.pca(x6G08_MD_441_raw,pcadata)
points(x6G08_MD_441_proj[1],x6G08_MD_441_proj[2],pch=20)
text(x6G08_MD_441_proj[1],x6G08_MD_441_proj[2],pos=4,label="441",col="red")

x6G08_MD_442_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_442.txt")
x6G08_MD_442_proj<-project.pca(x6G08_MD_442_raw,pcadata)
points(x6G08_MD_442_proj[1],x6G08_MD_442_proj[2],pch=20)
text(x6G08_MD_442_proj[1],x6G08_MD_442_proj[2],pos=4,label="442",col="red")

x6G08_MD_443_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_443.txt")
x6G08_MD_443_proj<-project.pca(x6G08_MD_443_raw,pcadata)
points(x6G08_MD_443_proj[1],x6G08_MD_443_proj[2],pch=20)
text(x6G08_MD_443_proj[1],x6G08_MD_443_proj[2],pos=4,label="443",col="red")

x6G08_MD_444_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_444.txt")
x6G08_MD_444_proj<-project.pca(x6G08_MD_444_raw,pcadata)
points(x6G08_MD_444_proj[1],x6G08_MD_444_proj[2],pch=20)
text(x6G08_MD_444_proj[1],x6G08_MD_444_proj[2],pos=4,label="444",col="red")

x6G08_MD_445_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_445.txt")
x6G08_MD_445_proj<-project.pca(x6G08_MD_445_raw,pcadata)
points(x6G08_MD_445_proj[1],x6G08_MD_445_proj[2],pch=20)
text(x6G08_MD_445_proj[1],x6G08_MD_445_proj[2],pos=4,label="445",col="red")

x6G08_MD_446_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_446.txt")
x6G08_MD_446_proj<-project.pca(x6G08_MD_446_raw,pcadata)
points(x6G08_MD_446_proj[1],x6G08_MD_446_proj[2],pch=20)
text(x6G08_MD_446_proj[1],x6G08_MD_446_proj[2],pos=4,label="446",col="red")

x6G08_MD_447_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_447.txt")
x6G08_MD_447_proj<-project.pca(x6G08_MD_447_raw,pcadata)
points(x6G08_MD_447_proj[1],x6G08_MD_447_proj[2],pch=20)
text(x6G08_MD_447_proj[1],x6G08_MD_447_proj[2],pos=4,label="447",col="red")

x6G08_MD_448_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_448.txt")
x6G08_MD_448_proj<-project.pca(x6G08_MD_448_raw,pcadata)
points(x6G08_MD_448_proj[1],x6G08_MD_448_proj[2],pch=20)
text(x6G08_MD_448_proj[1],x6G08_MD_448_proj[2],pos=4,label="448",col="red")

x6G08_MD_449_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_449.txt")
x6G08_MD_449_proj<-project.pca(x6G08_MD_449_raw,pcadata)
points(x6G08_MD_449_proj[1],x6G08_MD_449_proj[2],pch=20)
text(x6G08_MD_449_proj[1],x6G08_MD_449_proj[2],pos=4,label="449",col="red")

x6G08_MD_450_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_450.txt")
x6G08_MD_450_proj<-project.pca(x6G08_MD_450_raw,pcadata)
points(x6G08_MD_450_proj[1],x6G08_MD_450_proj[2],pch=20)
text(x6G08_MD_450_proj[1],x6G08_MD_450_proj[2],pos=4,label="450",col="red")

x6G08_MD_451_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_451.txt")
x6G08_MD_451_proj<-project.pca(x6G08_MD_451_raw,pcadata)
points(x6G08_MD_451_proj[1],x6G08_MD_451_proj[2],pch=20)
text(x6G08_MD_451_proj[1],x6G08_MD_451_proj[2],pos=4,label="451",col="red")

x6G08_MD_452_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_452.txt")
x6G08_MD_452_proj<-project.pca(x6G08_MD_452_raw,pcadata)
points(x6G08_MD_452_proj[1],x6G08_MD_452_proj[2],pch=20)
text(x6G08_MD_452_proj[1],x6G08_MD_452_proj[2],pos=4,label="452",col="red")

x6G08_MD_453_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_453.txt")
x6G08_MD_453_proj<-project.pca(x6G08_MD_453_raw,pcadata)
points(x6G08_MD_453_proj[1],x6G08_MD_453_proj[2],pch=20)
text(x6G08_MD_453_proj[1],x6G08_MD_453_proj[2],pos=4,label="453",col="red")

x6G08_MD_454_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_454.txt")
x6G08_MD_454_proj<-project.pca(x6G08_MD_454_raw,pcadata)
points(x6G08_MD_454_proj[1],x6G08_MD_454_proj[2],pch=20)
text(x6G08_MD_454_proj[1],x6G08_MD_454_proj[2],pos=4,label="454",col="green")

x6G08_MD_455_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_455.txt")
x6G08_MD_455_proj<-project.pca(x6G08_MD_455_raw,pcadata)
points(x6G08_MD_455_proj[1],x6G08_MD_455_proj[2],pch=20)
text(x6G08_MD_455_proj[1],x6G08_MD_455_proj[2],pos=4,label="455",col="green")

x6G08_MD_456_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_456.txt")
x6G08_MD_456_proj<-project.pca(x6G08_MD_456_raw,pcadata)
points(x6G08_MD_456_proj[1],x6G08_MD_456_proj[2],pch=20)
text(x6G08_MD_456_proj[1],x6G08_MD_456_proj[2],pos=4,label="456",col="red")

x6G08_MD_457_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_457.txt")
x6G08_MD_457_proj<-project.pca(x6G08_MD_457_raw,pcadata)
points(x6G08_MD_457_proj[1],x6G08_MD_457_proj[2],pch=20)
text(x6G08_MD_457_proj[1],x6G08_MD_457_proj[2],pos=4,label="457",col="red")

x6G08_MD_458_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_458.txt")
x6G08_MD_458_proj<-project.pca(x6G08_MD_458_raw,pcadata)
points(x6G08_MD_458_proj[1],x6G08_MD_458_proj[2],pch=20)
text(x6G08_MD_458_proj[1],x6G08_MD_458_proj[2],pos=4,label="458",col="red")

x6G08_MD_459_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_459.txt")
x6G08_MD_459_proj<-project.pca(x6G08_MD_459_raw,pcadata)
points(x6G08_MD_459_proj[1],x6G08_MD_459_proj[2],pch=20)
text(x6G08_MD_459_proj[1],x6G08_MD_459_proj[2],pos=4,label="459",col="red")

x6G08_MD_460_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_460.txt")
x6G08_MD_460_proj<-project.pca(x6G08_MD_460_raw,pcadata)
points(x6G08_MD_460_proj[1],x6G08_MD_460_proj[2],pch=20)
text(x6G08_MD_460_proj[1],x6G08_MD_460_proj[2],pos=4,label="460",col="blue")

x6G08_MD_461_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_461.txt")
x6G08_MD_461_proj<-project.pca(x6G08_MD_461_raw,pcadata)
points(x6G08_MD_461_proj[1],x6G08_MD_461_proj[2],pch=20)
text(x6G08_MD_461_proj[1],x6G08_MD_461_proj[2],pos=4,label="461",col="blue")

x6G08_MD_462_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_462.txt")
x6G08_MD_462_proj<-project.pca(x6G08_MD_462_raw,pcadata)
points(x6G08_MD_462_proj[1],x6G08_MD_462_proj[2],pch=20)
text(x6G08_MD_462_proj[1],x6G08_MD_462_proj[2],pos=4,label="462",col="green")

x6G08_MD_463_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_463.txt")
x6G08_MD_463_proj<-project.pca(x6G08_MD_463_raw,pcadata)
points(x6G08_MD_463_proj[1],x6G08_MD_463_proj[2],pch=20)
text(x6G08_MD_463_proj[1],x6G08_MD_463_proj[2],pos=4,label="463",col="blue")

x6G08_MD_464_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_464.txt")
x6G08_MD_464_proj<-project.pca(x6G08_MD_464_raw,pcadata)
points(x6G08_MD_464_proj[1],x6G08_MD_464_proj[2],pch=20)
text(x6G08_MD_464_proj[1],x6G08_MD_464_proj[2],pos=4,label="464",col="green")

x6G08_MD_465_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_465.txt")
x6G08_MD_465_proj<-project.pca(x6G08_MD_465_raw,pcadata)
points(x6G08_MD_465_proj[1],x6G08_MD_465_proj[2],pch=20)
text(x6G08_MD_465_proj[1],x6G08_MD_465_proj[2],pos=4,label="465",col="green")

x6G08_MD_466_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_466.txt")
x6G08_MD_466_proj<-project.pca(x6G08_MD_466_raw,pcadata)
points(x6G08_MD_466_proj[1],x6G08_MD_466_proj[2],pch=20)
text(x6G08_MD_466_proj[1],x6G08_MD_466_proj[2],pos=4,label="466",col="blue")

x6G08_MD_467_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_467.txt")
x6G08_MD_467_proj<-project.pca(x6G08_MD_467_raw,pcadata)
points(x6G08_MD_467_proj[1],x6G08_MD_467_proj[2],pch=20)
text(x6G08_MD_467_proj[1],x6G08_MD_467_proj[2],pos=4,label="467",col="blue")

x6G08_MD_468_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_468.txt")
x6G08_MD_468_proj<-project.pca(x6G08_MD_468_raw,pcadata)
points(x6G08_MD_468_proj[1],x6G08_MD_468_proj[2],pch=20)
text(x6G08_MD_468_proj[1],x6G08_MD_468_proj[2],pos=4,label="468",col="blue")

x6G08_MD_469_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_469.txt")
x6G08_MD_469_proj<-project.pca(x6G08_MD_469_raw,pcadata)
points(x6G08_MD_469_proj[1],x6G08_MD_469_proj[2],pch=20)
text(x6G08_MD_469_proj[1],x6G08_MD_469_proj[2],pos=4,label="469",col="black")

x6G08_MD_470_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_470.txt")
x6G08_MD_470_proj<-project.pca(x6G08_MD_470_raw,pcadata)
points(x6G08_MD_470_proj[1],x6G08_MD_470_proj[2],pch=20)
text(x6G08_MD_470_proj[1],x6G08_MD_470_proj[2],pos=4,label="470",col="blue")

x6G08_MD_471_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_471.txt")
x6G08_MD_471_proj<-project.pca(x6G08_MD_471_raw,pcadata)
points(x6G08_MD_471_proj[1],x6G08_MD_471_proj[2],pch=20)
text(x6G08_MD_471_proj[1],x6G08_MD_471_proj[2],pos=4,label="471",col="green")

x6G08_MD_472_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_472.txt")
x6G08_MD_472_proj<-project.pca(x6G08_MD_472_raw,pcadata)
points(x6G08_MD_472_proj[1],x6G08_MD_472_proj[2],pch=20)
text(x6G08_MD_472_proj[1],x6G08_MD_472_proj[2],pos=4,label="472",col="blue")

x6G08_MD_473_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_473.txt")
x6G08_MD_473_proj<-project.pca(x6G08_MD_473_raw,pcadata)
points(x6G08_MD_473_proj[1],x6G08_MD_473_proj[2],pch=20)
text(x6G08_MD_473_proj[1],x6G08_MD_473_proj[2],pos=4,label="473",col="red")

x6G08_MD_474_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_474.txt")
x6G08_MD_474_proj<-project.pca(x6G08_MD_474_raw,pcadata)
points(x6G08_MD_474_proj[1],x6G08_MD_474_proj[2],pch=20)
text(x6G08_MD_474_proj[1],x6G08_MD_474_proj[2],pos=4,label="474",col="green")

x6G08_MD_475_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_475.txt")
x6G08_MD_475_proj<-project.pca(x6G08_MD_475_raw,pcadata)
points(x6G08_MD_475_proj[1],x6G08_MD_475_proj[2],pch=20)
text(x6G08_MD_475_proj[1],x6G08_MD_475_proj[2],pos=4,label="475",col="green")

x6G08_MD_476_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_476.txt")
x6G08_MD_476_proj<-project.pca(x6G08_MD_476_raw,pcadata)
points(x6G08_MD_476_proj[1],x6G08_MD_476_proj[2],pch=20)
text(x6G08_MD_476_proj[1],x6G08_MD_476_proj[2],pos=4,label="476",col="red")

x6G08_MD_477_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_477.txt")
x6G08_MD_477_proj<-project.pca(x6G08_MD_477_raw,pcadata)
points(x6G08_MD_477_proj[1],x6G08_MD_477_proj[2],pch=20)
text(x6G08_MD_477_proj[1],x6G08_MD_477_proj[2],pos=4,label="477",col="green")

x6G08_MD_478_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_478.txt")
x6G08_MD_478_proj<-project.pca(x6G08_MD_478_raw,pcadata)
points(x6G08_MD_478_proj[1],x6G08_MD_478_proj[2],pch=20)
text(x6G08_MD_478_proj[1],x6G08_MD_478_proj[2],pos=4,label="478",col="blue")

x6G08_MD_479_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_479.txt")
x6G08_MD_479_proj<-project.pca(x6G08_MD_479_raw,pcadata)
points(x6G08_MD_479_proj[1],x6G08_MD_479_proj[2],pch=20)
text(x6G08_MD_479_proj[1],x6G08_MD_479_proj[2],pos=4,label="479",col="black")

x6G08_MD_480_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_480.txt")
x6G08_MD_480_proj<-project.pca(x6G08_MD_480_raw,pcadata)
points(x6G08_MD_480_proj[1],x6G08_MD_480_proj[2],pch=20)
text(x6G08_MD_480_proj[1],x6G08_MD_480_proj[2],pos=4,label="480",col="blue")

x6G08_MD_481_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_481.txt")
x6G08_MD_481_proj<-project.pca(x6G08_MD_481_raw,pcadata)
points(x6G08_MD_481_proj[1],x6G08_MD_481_proj[2],pch=20)
text(x6G08_MD_481_proj[1],x6G08_MD_481_proj[2],pos=4,label="481",col="blue")

x6G08_MD_482_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_482.txt")
x6G08_MD_482_proj<-project.pca(x6G08_MD_482_raw,pcadata)
points(x6G08_MD_482_proj[1],x6G08_MD_482_proj[2],pch=20)
text(x6G08_MD_482_proj[1],x6G08_MD_482_proj[2],pos=4,label="482",col="black")

x6G08_MD_483_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_483.txt")
x6G08_MD_483_proj<-project.pca(x6G08_MD_483_raw,pcadata)
points(x6G08_MD_483_proj[1],x6G08_MD_483_proj[2],pch=20)
text(x6G08_MD_483_proj[1],x6G08_MD_483_proj[2],pos=4,label="483",col="blue")

x6G08_MD_484_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_484.txt")
x6G08_MD_484_proj<-project.pca(x6G08_MD_484_raw,pcadata)
points(x6G08_MD_484_proj[1],x6G08_MD_484_proj[2],pch=20)
text(x6G08_MD_484_proj[1],x6G08_MD_484_proj[2],pos=4,label="484",col="black")

x6G08_MD_485_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_485.txt")
x6G08_MD_485_proj<-project.pca(x6G08_MD_485_raw,pcadata)
points(x6G08_MD_485_proj[1],x6G08_MD_485_proj[2],pch=20)
text(x6G08_MD_485_proj[1],x6G08_MD_485_proj[2],pos=4,label="485",col="blue")

x6G08_MD_486_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_486.txt")
x6G08_MD_486_proj<-project.pca(x6G08_MD_486_raw,pcadata)
points(x6G08_MD_486_proj[1],x6G08_MD_486_proj[2],pch=20)
text(x6G08_MD_486_proj[1],x6G08_MD_486_proj[2],pos=4,label="486",col="black")

x6G08_MD_487_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_487.txt")
x6G08_MD_487_proj<-project.pca(x6G08_MD_487_raw,pcadata)
points(x6G08_MD_487_proj[1],x6G08_MD_487_proj[2],pch=20)
text(x6G08_MD_487_proj[1],x6G08_MD_487_proj[2],pos=4,label="487",col="blue")

x6G08_MD_488_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_488.txt")
x6G08_MD_488_proj<-project.pca(x6G08_MD_488_raw,pcadata)
points(x6G08_MD_488_proj[1],x6G08_MD_488_proj[2],pch=20)
text(x6G08_MD_488_proj[1],x6G08_MD_488_proj[2],pos=4,label="488",col="blue")

x6G08_MD_489_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_489.txt")
x6G08_MD_489_proj<-project.pca(x6G08_MD_489_raw,pcadata)
points(x6G08_MD_489_proj[1],x6G08_MD_489_proj[2],pch=20)
text(x6G08_MD_489_proj[1],x6G08_MD_489_proj[2],pos=4,label="489",col="blue")

x6G08_MD_490_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_490.txt")
x6G08_MD_490_proj<-project.pca(x6G08_MD_490_raw,pcadata)
points(x6G08_MD_490_proj[1],x6G08_MD_490_proj[2],pch=20)
text(x6G08_MD_490_proj[1],x6G08_MD_490_proj[2],pos=4,label="490",col="blue")

x6G08_MD_491_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_491.txt")
x6G08_MD_491_proj<-project.pca(x6G08_MD_491_raw,pcadata)
points(x6G08_MD_491_proj[1],x6G08_MD_491_proj[2],pch=20)
text(x6G08_MD_491_proj[1],x6G08_MD_491_proj[2],pos=4,label="491",col="blue")

x6G08_MD_492_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_492.txt")
x6G08_MD_492_proj<-project.pca(x6G08_MD_492_raw,pcadata)
points(x6G08_MD_492_proj[1],x6G08_MD_492_proj[2],pch=20)
text(x6G08_MD_492_proj[1],x6G08_MD_492_proj[2],pos=4,label="492",col="blue")

x6G08_MD_493_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_493.txt")
x6G08_MD_493_proj<-project.pca(x6G08_MD_493_raw,pcadata)
points(x6G08_MD_493_proj[1],x6G08_MD_493_proj[2],pch=20)
text(x6G08_MD_493_proj[1],x6G08_MD_493_proj[2],pos=4,label="493",col="green")

x6G08_MD_494_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_494.txt")
x6G08_MD_494_proj<-project.pca(x6G08_MD_494_raw,pcadata)
points(x6G08_MD_494_proj[1],x6G08_MD_494_proj[2],pch=20)
text(x6G08_MD_494_proj[1],x6G08_MD_494_proj[2],pos=4,label="494",col="green")

x6G08_MD_495_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_495.txt")
x6G08_MD_495_proj<-project.pca(x6G08_MD_495_raw,pcadata)
points(x6G08_MD_495_proj[1],x6G08_MD_495_proj[2],pch=20)
text(x6G08_MD_495_proj[1],x6G08_MD_495_proj[2],pos=4,label="495",col="red")

x6G08_MD_496_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_496.txt")
x6G08_MD_496_proj<-project.pca(x6G08_MD_496_raw,pcadata)
points(x6G08_MD_496_proj[1],x6G08_MD_496_proj[2],pch=20)
text(x6G08_MD_496_proj[1],x6G08_MD_496_proj[2],pos=4,label="496",col="green")

x6G08_MD_497_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_497.txt")
x6G08_MD_497_proj<-project.pca(x6G08_MD_497_raw,pcadata)
points(x6G08_MD_497_proj[1],x6G08_MD_497_proj[2],pch=20)
text(x6G08_MD_497_proj[1],x6G08_MD_497_proj[2],pos=4,label="497",col="green")

x6G08_MD_498_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_498.txt")
x6G08_MD_498_proj<-project.pca(x6G08_MD_498_raw,pcadata)
points(x6G08_MD_498_proj[1],x6G08_MD_498_proj[2],pch=20)
text(x6G08_MD_498_proj[1],x6G08_MD_498_proj[2],pos=4,label="498",col="black")

x6G08_MD_499_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_499.txt")
x6G08_MD_499_proj<-project.pca(x6G08_MD_499_raw,pcadata)
points(x6G08_MD_499_proj[1],x6G08_MD_499_proj[2],pch=20)
text(x6G08_MD_499_proj[1],x6G08_MD_499_proj[2],pos=4,label="499",col="green")

x6G08_MD_500_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_500.txt")
x6G08_MD_500_proj<-project.pca(x6G08_MD_500_raw,pcadata)
points(x6G08_MD_500_proj[1],x6G08_MD_500_proj[2],pch=20)
text(x6G08_MD_500_proj[1],x6G08_MD_500_proj[2],pos=4,label="500",col="red")

x6G08_MD_501_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_501.txt")
x6G08_MD_501_proj<-project.pca(x6G08_MD_501_raw,pcadata)
points(x6G08_MD_501_proj[1],x6G08_MD_501_proj[2],pch=20)
text(x6G08_MD_501_proj[1],x6G08_MD_501_proj[2],pos=4,label="501",col="red")

x6G08_MD_502_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_502.txt")
x6G08_MD_502_proj<-project.pca(x6G08_MD_502_raw,pcadata)
points(x6G08_MD_502_proj[1],x6G08_MD_502_proj[2],pch=20)
text(x6G08_MD_502_proj[1],x6G08_MD_502_proj[2],pos=4,label="502",col="blue")

x6G08_MD_503_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_503.txt")
x6G08_MD_503_proj<-project.pca(x6G08_MD_503_raw,pcadata)
points(x6G08_MD_503_proj[1],x6G08_MD_503_proj[2],pch=20)
text(x6G08_MD_503_proj[1],x6G08_MD_503_proj[2],pos=4,label="503",col="black")

x6G08_MD_504_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_504.txt")
x6G08_MD_504_proj<-project.pca(x6G08_MD_504_raw,pcadata)
points(x6G08_MD_504_proj[1],x6G08_MD_504_proj[2],pch=20)
text(x6G08_MD_504_proj[1],x6G08_MD_504_proj[2],pos=4,label="504",col="blue")

x6G08_MD_505_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_505.txt")
x6G08_MD_505_proj<-project.pca(x6G08_MD_505_raw,pcadata)
points(x6G08_MD_505_proj[1],x6G08_MD_505_proj[2],pch=20)
text(x6G08_MD_505_proj[1],x6G08_MD_505_proj[2],pos=4,label="505",col="green")

x6G08_MD_506_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_506.txt")
x6G08_MD_506_proj<-project.pca(x6G08_MD_506_raw,pcadata)
points(x6G08_MD_506_proj[1],x6G08_MD_506_proj[2],pch=20)
text(x6G08_MD_506_proj[1],x6G08_MD_506_proj[2],pos=4,label="506",col="blue")

x6G08_MD_507_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_507.txt")
x6G08_MD_507_proj<-project.pca(x6G08_MD_507_raw,pcadata)
points(x6G08_MD_507_proj[1],x6G08_MD_507_proj[2],pch=20)
text(x6G08_MD_507_proj[1],x6G08_MD_507_proj[2],pos=4,label="507",col="blue")

x6G08_MD_508_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_508.txt")
x6G08_MD_508_proj<-project.pca(x6G08_MD_508_raw,pcadata)
points(x6G08_MD_508_proj[1],x6G08_MD_508_proj[2],pch=20)
text(x6G08_MD_508_proj[1],x6G08_MD_508_proj[2],pos=4,label="508",col="blue")

x6G08_MD_509_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_509.txt")
x6G08_MD_509_proj<-project.pca(x6G08_MD_509_raw,pcadata)
points(x6G08_MD_509_proj[1],x6G08_MD_509_proj[2],pch=20)
text(x6G08_MD_509_proj[1],x6G08_MD_509_proj[2],pos=4,label="509",col="blue")

x6G08_MD_510_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_510.txt")
x6G08_MD_510_proj<-project.pca(x6G08_MD_510_raw,pcadata)
points(x6G08_MD_510_proj[1],x6G08_MD_510_proj[2],pch=20)
text(x6G08_MD_510_proj[1],x6G08_MD_510_proj[2],pos=4,label="510",col="blue")

x6G08_MD_511_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_511.txt")
x6G08_MD_511_proj<-project.pca(x6G08_MD_511_raw,pcadata)
points(x6G08_MD_511_proj[1],x6G08_MD_511_proj[2],pch=20)
text(x6G08_MD_511_proj[1],x6G08_MD_511_proj[2],pos=4,label="511",col="blue")

x6G08_MD_512_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_512.txt")
x6G08_MD_512_proj<-project.pca(x6G08_MD_512_raw,pcadata)
points(x6G08_MD_512_proj[1],x6G08_MD_512_proj[2],pch=20)
text(x6G08_MD_512_proj[1],x6G08_MD_512_proj[2],pos=4,label="512",col="blue")

x6G08_MD_513_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_513.txt")
x6G08_MD_513_proj<-project.pca(x6G08_MD_513_raw,pcadata)
points(x6G08_MD_513_proj[1],x6G08_MD_513_proj[2],pch=20)
text(x6G08_MD_513_proj[1],x6G08_MD_513_proj[2],pos=4,label="513",col="black")

x6G08_MD_514_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_514.txt")
x6G08_MD_514_proj<-project.pca(x6G08_MD_514_raw,pcadata)
points(x6G08_MD_514_proj[1],x6G08_MD_514_proj[2],pch=20)
text(x6G08_MD_514_proj[1],x6G08_MD_514_proj[2],pos=4,label="514",col="black")

x6G08_MD_515_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_515.txt")
x6G08_MD_515_proj<-project.pca(x6G08_MD_515_raw,pcadata)
points(x6G08_MD_515_proj[1],x6G08_MD_515_proj[2],pch=20)
text(x6G08_MD_515_proj[1],x6G08_MD_515_proj[2],pos=4,label="515",col="blue")

x6G08_MD_516_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_516.txt")
x6G08_MD_516_proj<-project.pca(x6G08_MD_516_raw,pcadata)
points(x6G08_MD_516_proj[1],x6G08_MD_516_proj[2],pch=20)
text(x6G08_MD_516_proj[1],x6G08_MD_516_proj[2],pos=4,label="516",col="black")

x6G08_MD_517_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_517.txt")
x6G08_MD_517_proj<-project.pca(x6G08_MD_517_raw,pcadata)
points(x6G08_MD_517_proj[1],x6G08_MD_517_proj[2],pch=20)
text(x6G08_MD_517_proj[1],x6G08_MD_517_proj[2],pos=4,label="517",col="black")

x6G08_MD_518_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_518.txt")
x6G08_MD_518_proj<-project.pca(x6G08_MD_518_raw,pcadata)
points(x6G08_MD_518_proj[1],x6G08_MD_518_proj[2],pch=20)
text(x6G08_MD_518_proj[1],x6G08_MD_518_proj[2],pos=4,label="518",col="black")

x6G08_MD_519_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_519.txt")
x6G08_MD_519_proj<-project.pca(x6G08_MD_519_raw,pcadata)
points(x6G08_MD_519_proj[1],x6G08_MD_519_proj[2],pch=20)
text(x6G08_MD_519_proj[1],x6G08_MD_519_proj[2],pos=4,label="519",col="black")

x6G08_MD_520_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_520.txt")
x6G08_MD_520_proj<-project.pca(x6G08_MD_520_raw,pcadata)
points(x6G08_MD_520_proj[1],x6G08_MD_520_proj[2],pch=20)
text(x6G08_MD_520_proj[1],x6G08_MD_520_proj[2],pos=4,label="520",col="blue")

x6G08_MD_521_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_521.txt")
x6G08_MD_521_proj<-project.pca(x6G08_MD_521_raw,pcadata)
points(x6G08_MD_521_proj[1],x6G08_MD_521_proj[2],pch=20)
text(x6G08_MD_521_proj[1],x6G08_MD_521_proj[2],pos=4,label="521",col="blue")

x6G08_MD_522_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_522.txt")
x6G08_MD_522_proj<-project.pca(x6G08_MD_522_raw,pcadata)
points(x6G08_MD_522_proj[1],x6G08_MD_522_proj[2],pch=20)
text(x6G08_MD_522_proj[1],x6G08_MD_522_proj[2],pos=4,label="522",col="black")

x6G08_MD_523_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_523.txt")
x6G08_MD_523_proj<-project.pca(x6G08_MD_523_raw,pcadata)
points(x6G08_MD_523_proj[1],x6G08_MD_523_proj[2],pch=20)
text(x6G08_MD_523_proj[1],x6G08_MD_523_proj[2],pos=4,label="523",col="black")

x6G08_MD_524_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_524.txt")
x6G08_MD_524_proj<-project.pca(x6G08_MD_524_raw,pcadata)
points(x6G08_MD_524_proj[1],x6G08_MD_524_proj[2],pch=20)
text(x6G08_MD_524_proj[1],x6G08_MD_524_proj[2],pos=4,label="524",col="black")

x6G08_MD_525_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_525.txt")
x6G08_MD_525_proj<-project.pca(x6G08_MD_525_raw,pcadata)
points(x6G08_MD_525_proj[1],x6G08_MD_525_proj[2],pch=20)
text(x6G08_MD_525_proj[1],x6G08_MD_525_proj[2],pos=4,label="525",col="blue")

x6G08_MD_526_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_526.txt")
x6G08_MD_526_proj<-project.pca(x6G08_MD_526_raw,pcadata)
points(x6G08_MD_526_proj[1],x6G08_MD_526_proj[2],pch=20)
text(x6G08_MD_526_proj[1],x6G08_MD_526_proj[2],pos=4,label="526",col="green")

x6G08_MD_527_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_527.txt")
x6G08_MD_527_proj<-project.pca(x6G08_MD_527_raw,pcadata)
points(x6G08_MD_527_proj[1],x6G08_MD_527_proj[2],pch=20)
text(x6G08_MD_527_proj[1],x6G08_MD_527_proj[2],pos=4,label="527",col="blue")

x6G08_MD_528_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_528.txt")
x6G08_MD_528_proj<-project.pca(x6G08_MD_528_raw,pcadata)
points(x6G08_MD_528_proj[1],x6G08_MD_528_proj[2],pch=20)
text(x6G08_MD_528_proj[1],x6G08_MD_528_proj[2],pos=4,label="528",col="green")

x6G08_MD_529_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_529.txt")
x6G08_MD_529_proj<-project.pca(x6G08_MD_529_raw,pcadata)
points(x6G08_MD_529_proj[1],x6G08_MD_529_proj[2],pch=20)
text(x6G08_MD_529_proj[1],x6G08_MD_529_proj[2],pos=4,label="529",col="blue")

x6G08_MD_530_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_530.txt")
x6G08_MD_530_proj<-project.pca(x6G08_MD_530_raw,pcadata)
points(x6G08_MD_530_proj[1],x6G08_MD_530_proj[2],pch=20)
text(x6G08_MD_530_proj[1],x6G08_MD_530_proj[2],pos=4,label="530",col="blue")

x6G08_MD_531_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_531.txt")
x6G08_MD_531_proj<-project.pca(x6G08_MD_531_raw,pcadata)
points(x6G08_MD_531_proj[1],x6G08_MD_531_proj[2],pch=20)
text(x6G08_MD_531_proj[1],x6G08_MD_531_proj[2],pos=4,label="531",col="blue")

x6G08_MD_532_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_532.txt")
x6G08_MD_532_proj<-project.pca(x6G08_MD_532_raw,pcadata)
points(x6G08_MD_532_proj[1],x6G08_MD_532_proj[2],pch=20)
text(x6G08_MD_532_proj[1],x6G08_MD_532_proj[2],pos=4,label="532",col="red")

x6G08_MD_533_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_533.txt")
x6G08_MD_533_proj<-project.pca(x6G08_MD_533_raw,pcadata)
points(x6G08_MD_533_proj[1],x6G08_MD_533_proj[2],pch=20)
text(x6G08_MD_533_proj[1],x6G08_MD_533_proj[2],pos=4,label="533",col="blue")

x6G08_MD_534_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_534.txt")
x6G08_MD_534_proj<-project.pca(x6G08_MD_534_raw,pcadata)
points(x6G08_MD_534_proj[1],x6G08_MD_534_proj[2],pch=20)
text(x6G08_MD_534_proj[1],x6G08_MD_534_proj[2],pos=4,label="534",col="blue")

x6G08_MD_535_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_535.txt")
x6G08_MD_535_proj<-project.pca(x6G08_MD_535_raw,pcadata)
points(x6G08_MD_535_proj[1],x6G08_MD_535_proj[2],pch=20)
text(x6G08_MD_535_proj[1],x6G08_MD_535_proj[2],pos=4,label="535",col="red")

x6G08_MD_536_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_536.txt")
x6G08_MD_536_proj<-project.pca(x6G08_MD_536_raw,pcadata)
points(x6G08_MD_536_proj[1],x6G08_MD_536_proj[2],pch=20)
text(x6G08_MD_536_proj[1],x6G08_MD_536_proj[2],pos=4,label="536",col="black")

x6G08_MD_537_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_537.txt")
x6G08_MD_537_proj<-project.pca(x6G08_MD_537_raw,pcadata)
points(x6G08_MD_537_proj[1],x6G08_MD_537_proj[2],pch=20)
text(x6G08_MD_537_proj[1],x6G08_MD_537_proj[2],pos=4,label="537",col="blue")

x6G08_MD_538_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_538.txt")
x6G08_MD_538_proj<-project.pca(x6G08_MD_538_raw,pcadata)
points(x6G08_MD_538_proj[1],x6G08_MD_538_proj[2],pch=20)
text(x6G08_MD_538_proj[1],x6G08_MD_538_proj[2],pos=4,label="538",col="red")

x6G08_MD_539_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_539.txt")
x6G08_MD_539_proj<-project.pca(x6G08_MD_539_raw,pcadata)
points(x6G08_MD_539_proj[1],x6G08_MD_539_proj[2],pch=20)
text(x6G08_MD_539_proj[1],x6G08_MD_539_proj[2],pos=4,label="539",col="green")

x6G08_MD_540_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_540.txt")
x6G08_MD_540_proj<-project.pca(x6G08_MD_540_raw,pcadata)
points(x6G08_MD_540_proj[1],x6G08_MD_540_proj[2],pch=20)
text(x6G08_MD_540_proj[1],x6G08_MD_540_proj[2],pos=4,label="540",col="red")

x6G08_MD_541_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_541.txt")
x6G08_MD_541_proj<-project.pca(x6G08_MD_541_raw,pcadata)
points(x6G08_MD_541_proj[1],x6G08_MD_541_proj[2],pch=20)
text(x6G08_MD_541_proj[1],x6G08_MD_541_proj[2],pos=4,label="541",col="red")

x6G08_MD_542_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_542.txt")
x6G08_MD_542_proj<-project.pca(x6G08_MD_542_raw,pcadata)
points(x6G08_MD_542_proj[1],x6G08_MD_542_proj[2],pch=20)
text(x6G08_MD_542_proj[1],x6G08_MD_542_proj[2],pos=4,label="542",col="green")

x6G08_MD_543_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_543.txt")
x6G08_MD_543_proj<-project.pca(x6G08_MD_543_raw,pcadata)
points(x6G08_MD_543_proj[1],x6G08_MD_543_proj[2],pch=20)
text(x6G08_MD_543_proj[1],x6G08_MD_543_proj[2],pos=4,label="543",col="green")

x6G08_MD_544_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_544.txt")
x6G08_MD_544_proj<-project.pca(x6G08_MD_544_raw,pcadata)
points(x6G08_MD_544_proj[1],x6G08_MD_544_proj[2],pch=20)
text(x6G08_MD_544_proj[1],x6G08_MD_544_proj[2],pos=4,label="544",col="green")

x6G08_MD_545_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_545.txt")
x6G08_MD_545_proj<-project.pca(x6G08_MD_545_raw,pcadata)
points(x6G08_MD_545_proj[1],x6G08_MD_545_proj[2],pch=20)
text(x6G08_MD_545_proj[1],x6G08_MD_545_proj[2],pos=4,label="545",col="blue")

x6G08_MD_546_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_546.txt")
x6G08_MD_546_proj<-project.pca(x6G08_MD_546_raw,pcadata)
points(x6G08_MD_546_proj[1],x6G08_MD_546_proj[2],pch=20)
text(x6G08_MD_546_proj[1],x6G08_MD_546_proj[2],pos=4,label="546",col="red")

x6G08_MD_547_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_547.txt")
x6G08_MD_547_proj<-project.pca(x6G08_MD_547_raw,pcadata)
points(x6G08_MD_547_proj[1],x6G08_MD_547_proj[2],pch=20)
text(x6G08_MD_547_proj[1],x6G08_MD_547_proj[2],pos=4,label="547",col="green")

x6G08_MD_548_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_548.txt")
x6G08_MD_548_proj<-project.pca(x6G08_MD_548_raw,pcadata)
points(x6G08_MD_548_proj[1],x6G08_MD_548_proj[2],pch=20)
text(x6G08_MD_548_proj[1],x6G08_MD_548_proj[2],pos=4,label="548",col="blue")

x6G08_MD_549_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_549.txt")
x6G08_MD_549_proj<-project.pca(x6G08_MD_549_raw,pcadata)
points(x6G08_MD_549_proj[1],x6G08_MD_549_proj[2],pch=20)
text(x6G08_MD_549_proj[1],x6G08_MD_549_proj[2],pos=4,label="549",col="red")

x6G08_MD_550_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_550.txt")
x6G08_MD_550_proj<-project.pca(x6G08_MD_550_raw,pcadata)
points(x6G08_MD_550_proj[1],x6G08_MD_550_proj[2],pch=20)
text(x6G08_MD_550_proj[1],x6G08_MD_550_proj[2],pos=4,label="550",col="green")

x6G08_MD_551_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_551.txt")
x6G08_MD_551_proj<-project.pca(x6G08_MD_551_raw,pcadata)
points(x6G08_MD_551_proj[1],x6G08_MD_551_proj[2],pch=20)
text(x6G08_MD_551_proj[1],x6G08_MD_551_proj[2],pos=4,label="551",col="red")

x6G08_MD_552_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_552.txt")
x6G08_MD_552_proj<-project.pca(x6G08_MD_552_raw,pcadata)
points(x6G08_MD_552_proj[1],x6G08_MD_552_proj[2],pch=20)
text(x6G08_MD_552_proj[1],x6G08_MD_552_proj[2],pos=4,label="552",col="red")

x6G08_MD_553_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_553.txt")
x6G08_MD_553_proj<-project.pca(x6G08_MD_553_raw,pcadata)
points(x6G08_MD_553_proj[1],x6G08_MD_553_proj[2],pch=20)
text(x6G08_MD_553_proj[1],x6G08_MD_553_proj[2],pos=4,label="553",col="green")

x6G08_MD_554_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_554.txt")
x6G08_MD_554_proj<-project.pca(x6G08_MD_554_raw,pcadata)
points(x6G08_MD_554_proj[1],x6G08_MD_554_proj[2],pch=20)
text(x6G08_MD_554_proj[1],x6G08_MD_554_proj[2],pos=4,label="554",col="green")

x6G08_MD_555_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_555.txt")
x6G08_MD_555_proj<-project.pca(x6G08_MD_555_raw,pcadata)
points(x6G08_MD_555_proj[1],x6G08_MD_555_proj[2],pch=20)
text(x6G08_MD_555_proj[1],x6G08_MD_555_proj[2],pos=4,label="555",col="blue")

x6G08_MD_556_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_556.txt")
x6G08_MD_556_proj<-project.pca(x6G08_MD_556_raw,pcadata)
points(x6G08_MD_556_proj[1],x6G08_MD_556_proj[2],pch=20)
text(x6G08_MD_556_proj[1],x6G08_MD_556_proj[2],pos=4,label="556",col="green")

x6G08_MD_557_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_557.txt")
x6G08_MD_557_proj<-project.pca(x6G08_MD_557_raw,pcadata)
points(x6G08_MD_557_proj[1],x6G08_MD_557_proj[2],pch=20)
text(x6G08_MD_557_proj[1],x6G08_MD_557_proj[2],pos=4,label="557",col="blue")

x6G08_MD_558_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_558.txt")
x6G08_MD_558_proj<-project.pca(x6G08_MD_558_raw,pcadata)
points(x6G08_MD_558_proj[1],x6G08_MD_558_proj[2],pch=20)
text(x6G08_MD_558_proj[1],x6G08_MD_558_proj[2],pos=4,label="558",col="blue")

x6G08_MD_559_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_559.txt")
x6G08_MD_559_proj<-project.pca(x6G08_MD_559_raw,pcadata)
points(x6G08_MD_559_proj[1],x6G08_MD_559_proj[2],pch=20)
text(x6G08_MD_559_proj[1],x6G08_MD_559_proj[2],pos=4,label="559",col="blue")

x6G08_MD_560_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_560.txt")
x6G08_MD_560_proj<-project.pca(x6G08_MD_560_raw,pcadata)
points(x6G08_MD_560_proj[1],x6G08_MD_560_proj[2],pch=20)
text(x6G08_MD_560_proj[1],x6G08_MD_560_proj[2],pos=4,label="560",col="blue")

x6G08_MD_561_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_561.txt")
x6G08_MD_561_proj<-project.pca(x6G08_MD_561_raw,pcadata)
points(x6G08_MD_561_proj[1],x6G08_MD_561_proj[2],pch=20)
text(x6G08_MD_561_proj[1],x6G08_MD_561_proj[2],pos=4,label="561",col="red")

x6G08_MD_562_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_562.txt")
x6G08_MD_562_proj<-project.pca(x6G08_MD_562_raw,pcadata)
points(x6G08_MD_562_proj[1],x6G08_MD_562_proj[2],pch=20)
text(x6G08_MD_562_proj[1],x6G08_MD_562_proj[2],pos=4,label="562",col="green")

x6G08_MD_563_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_563.txt")
x6G08_MD_563_proj<-project.pca(x6G08_MD_563_raw,pcadata)
points(x6G08_MD_563_proj[1],x6G08_MD_563_proj[2],pch=20)
text(x6G08_MD_563_proj[1],x6G08_MD_563_proj[2],pos=4,label="563",col="red")

x6G08_MD_564_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_564.txt")
x6G08_MD_564_proj<-project.pca(x6G08_MD_564_raw,pcadata)
points(x6G08_MD_564_proj[1],x6G08_MD_564_proj[2],pch=20)
text(x6G08_MD_564_proj[1],x6G08_MD_564_proj[2],pos=4,label="564",col="blue")

x6G08_MD_565_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_565.txt")
x6G08_MD_565_proj<-project.pca(x6G08_MD_565_raw,pcadata)
points(x6G08_MD_565_proj[1],x6G08_MD_565_proj[2],pch=20)
text(x6G08_MD_565_proj[1],x6G08_MD_565_proj[2],pos=4,label="565",col="blue")

x6G08_MD_566_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_566.txt")
x6G08_MD_566_proj<-project.pca(x6G08_MD_566_raw,pcadata)
points(x6G08_MD_566_proj[1],x6G08_MD_566_proj[2],pch=20)
text(x6G08_MD_566_proj[1],x6G08_MD_566_proj[2],pos=4,label="566",col="red")

x6G08_MD_567_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_567.txt")
x6G08_MD_567_proj<-project.pca(x6G08_MD_567_raw,pcadata)
points(x6G08_MD_567_proj[1],x6G08_MD_567_proj[2],pch=20)
text(x6G08_MD_567_proj[1],x6G08_MD_567_proj[2],pos=4,label="567",col="red")

x6G08_MD_568_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_568.txt")
x6G08_MD_568_proj<-project.pca(x6G08_MD_568_raw,pcadata)
points(x6G08_MD_568_proj[1],x6G08_MD_568_proj[2],pch=20)
text(x6G08_MD_568_proj[1],x6G08_MD_568_proj[2],pos=4,label="568",col="blue")

x6G08_MD_569_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_569.txt")
x6G08_MD_569_proj<-project.pca(x6G08_MD_569_raw,pcadata)
points(x6G08_MD_569_proj[1],x6G08_MD_569_proj[2],pch=20)
text(x6G08_MD_569_proj[1],x6G08_MD_569_proj[2],pos=4,label="569",col="green")

x6G08_MD_570_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_570.txt")
x6G08_MD_570_proj<-project.pca(x6G08_MD_570_raw,pcadata)
points(x6G08_MD_570_proj[1],x6G08_MD_570_proj[2],pch=20)
text(x6G08_MD_570_proj[1],x6G08_MD_570_proj[2],pos=4,label="570",col="green")

x6G08_MD_571_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_571.txt")
x6G08_MD_571_proj<-project.pca(x6G08_MD_571_raw,pcadata)
points(x6G08_MD_571_proj[1],x6G08_MD_571_proj[2],pch=20)
text(x6G08_MD_571_proj[1],x6G08_MD_571_proj[2],pos=4,label="571",col="red")

x6G08_MD_572_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_572.txt")
x6G08_MD_572_proj<-project.pca(x6G08_MD_572_raw,pcadata)
points(x6G08_MD_572_proj[1],x6G08_MD_572_proj[2],pch=20)
text(x6G08_MD_572_proj[1],x6G08_MD_572_proj[2],pos=4,label="572",col="blue")

x6G08_MD_573_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_573.txt")
x6G08_MD_573_proj<-project.pca(x6G08_MD_573_raw,pcadata)
points(x6G08_MD_573_proj[1],x6G08_MD_573_proj[2],pch=20)
text(x6G08_MD_573_proj[1],x6G08_MD_573_proj[2],pos=4,label="573",col="black")

x6G08_MD_574_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_574.txt")
x6G08_MD_574_proj<-project.pca(x6G08_MD_574_raw,pcadata)
points(x6G08_MD_574_proj[1],x6G08_MD_574_proj[2],pch=20)
text(x6G08_MD_574_proj[1],x6G08_MD_574_proj[2],pos=4,label="574",col="black")

x6G08_MD_575_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_575.txt")
x6G08_MD_575_proj<-project.pca(x6G08_MD_575_raw,pcadata)
points(x6G08_MD_575_proj[1],x6G08_MD_575_proj[2],pch=20)
text(x6G08_MD_575_proj[1],x6G08_MD_575_proj[2],pos=4,label="575",col="red")

x6G08_MD_576_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_576.txt")
x6G08_MD_576_proj<-project.pca(x6G08_MD_576_raw,pcadata)
points(x6G08_MD_576_proj[1],x6G08_MD_576_proj[2],pch=20)
text(x6G08_MD_576_proj[1],x6G08_MD_576_proj[2],pos=4,label="576",col="blue")

x6G08_MD_577_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_577.txt")
x6G08_MD_577_proj<-project.pca(x6G08_MD_577_raw,pcadata)
points(x6G08_MD_577_proj[1],x6G08_MD_577_proj[2],pch=20)
text(x6G08_MD_577_proj[1],x6G08_MD_577_proj[2],pos=4,label="577",col="black")

x6G08_MD_578_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_578.txt")
x6G08_MD_578_proj<-project.pca(x6G08_MD_578_raw,pcadata)
points(x6G08_MD_578_proj[1],x6G08_MD_578_proj[2],pch=20)
text(x6G08_MD_578_proj[1],x6G08_MD_578_proj[2],pos=4,label="578",col="blue")

x6G08_MD_579_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_579.txt")
x6G08_MD_579_proj<-project.pca(x6G08_MD_579_raw,pcadata)
points(x6G08_MD_579_proj[1],x6G08_MD_579_proj[2],pch=20)
text(x6G08_MD_579_proj[1],x6G08_MD_579_proj[2],pos=4,label="579",col="blue")

x6G08_MD_580_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_580.txt")
x6G08_MD_580_proj<-project.pca(x6G08_MD_580_raw,pcadata)
points(x6G08_MD_580_proj[1],x6G08_MD_580_proj[2],pch=20)
text(x6G08_MD_580_proj[1],x6G08_MD_580_proj[2],pos=4,label="580",col="green")

x6G08_MD_581_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_581.txt")
x6G08_MD_581_proj<-project.pca(x6G08_MD_581_raw,pcadata)
points(x6G08_MD_581_proj[1],x6G08_MD_581_proj[2],pch=20)
text(x6G08_MD_581_proj[1],x6G08_MD_581_proj[2],pos=4,label="581",col="green")

x6G08_MD_582_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_582.txt")
x6G08_MD_582_proj<-project.pca(x6G08_MD_582_raw,pcadata)
points(x6G08_MD_582_proj[1],x6G08_MD_582_proj[2],pch=20)
text(x6G08_MD_582_proj[1],x6G08_MD_582_proj[2],pos=4,label="582",col="green")

x6G08_MD_583_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_583.txt")
x6G08_MD_583_proj<-project.pca(x6G08_MD_583_raw,pcadata)
points(x6G08_MD_583_proj[1],x6G08_MD_583_proj[2],pch=20)
text(x6G08_MD_583_proj[1],x6G08_MD_583_proj[2],pos=4,label="583",col="blue")

x6G08_MD_584_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_584.txt")
x6G08_MD_584_proj<-project.pca(x6G08_MD_584_raw,pcadata)
points(x6G08_MD_584_proj[1],x6G08_MD_584_proj[2],pch=20)
text(x6G08_MD_584_proj[1],x6G08_MD_584_proj[2],pos=4,label="584",col="red")

x6G08_MD_585_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_585.txt")
x6G08_MD_585_proj<-project.pca(x6G08_MD_585_raw,pcadata)
points(x6G08_MD_585_proj[1],x6G08_MD_585_proj[2],pch=20)
text(x6G08_MD_585_proj[1],x6G08_MD_585_proj[2],pos=4,label="585",col="blue")

x6G08_MD_586_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_586.txt")
x6G08_MD_586_proj<-project.pca(x6G08_MD_586_raw,pcadata)
points(x6G08_MD_586_proj[1],x6G08_MD_586_proj[2],pch=20)
text(x6G08_MD_586_proj[1],x6G08_MD_586_proj[2],pos=4,label="586",col="blue")

x6G08_MD_587_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_587.txt")
x6G08_MD_587_proj<-project.pca(x6G08_MD_587_raw,pcadata)
points(x6G08_MD_587_proj[1],x6G08_MD_587_proj[2],pch=20)
text(x6G08_MD_587_proj[1],x6G08_MD_587_proj[2],pos=4,label="587",col="black")

x6G08_MD_588_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_588.txt")
x6G08_MD_588_proj<-project.pca(x6G08_MD_588_raw,pcadata)
points(x6G08_MD_588_proj[1],x6G08_MD_588_proj[2],pch=20)
text(x6G08_MD_588_proj[1],x6G08_MD_588_proj[2],pos=4,label="588",col="green")

x6G08_MD_589_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_589.txt")
x6G08_MD_589_proj<-project.pca(x6G08_MD_589_raw,pcadata)
points(x6G08_MD_589_proj[1],x6G08_MD_589_proj[2],pch=20)
text(x6G08_MD_589_proj[1],x6G08_MD_589_proj[2],pos=4,label="589",col="blue")

x6G08_MD_590_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_590.txt")
x6G08_MD_590_proj<-project.pca(x6G08_MD_590_raw,pcadata)
points(x6G08_MD_590_proj[1],x6G08_MD_590_proj[2],pch=20)
text(x6G08_MD_590_proj[1],x6G08_MD_590_proj[2],pos=4,label="590",col="black")

x6G08_MD_591_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_591.txt")
x6G08_MD_591_proj<-project.pca(x6G08_MD_591_raw,pcadata)
points(x6G08_MD_591_proj[1],x6G08_MD_591_proj[2],pch=20)
text(x6G08_MD_591_proj[1],x6G08_MD_591_proj[2],pos=4,label="591",col="blue")

x6G08_MD_592_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_592.txt")
x6G08_MD_592_proj<-project.pca(x6G08_MD_592_raw,pcadata)
points(x6G08_MD_592_proj[1],x6G08_MD_592_proj[2],pch=20)
text(x6G08_MD_592_proj[1],x6G08_MD_592_proj[2],pos=4,label="592",col="green")

x6G08_MD_593_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_593.txt")
x6G08_MD_593_proj<-project.pca(x6G08_MD_593_raw,pcadata)
points(x6G08_MD_593_proj[1],x6G08_MD_593_proj[2],pch=20)
text(x6G08_MD_593_proj[1],x6G08_MD_593_proj[2],pos=4,label="593",col="blue")

x6G08_MD_594_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_594.txt")
x6G08_MD_594_proj<-project.pca(x6G08_MD_594_raw,pcadata)
points(x6G08_MD_594_proj[1],x6G08_MD_594_proj[2],pch=20)
text(x6G08_MD_594_proj[1],x6G08_MD_594_proj[2],pos=4,label="594",col="blue")

x6G08_MD_595_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_595.txt")
x6G08_MD_595_proj<-project.pca(x6G08_MD_595_raw,pcadata)
points(x6G08_MD_595_proj[1],x6G08_MD_595_proj[2],pch=20)
text(x6G08_MD_595_proj[1],x6G08_MD_595_proj[2],pos=4,label="595",col="black")

x6G08_MD_596_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_596.txt")
x6G08_MD_596_proj<-project.pca(x6G08_MD_596_raw,pcadata)
points(x6G08_MD_596_proj[1],x6G08_MD_596_proj[2],pch=20)
text(x6G08_MD_596_proj[1],x6G08_MD_596_proj[2],pos=4,label="596",col="black")

x6G08_MD_597_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_597.txt")
x6G08_MD_597_proj<-project.pca(x6G08_MD_597_raw,pcadata)
points(x6G08_MD_597_proj[1],x6G08_MD_597_proj[2],pch=20)
text(x6G08_MD_597_proj[1],x6G08_MD_597_proj[2],pos=4,label="597",col="blue")

x6G08_MD_598_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_598.txt")
x6G08_MD_598_proj<-project.pca(x6G08_MD_598_raw,pcadata)
points(x6G08_MD_598_proj[1],x6G08_MD_598_proj[2],pch=20)
text(x6G08_MD_598_proj[1],x6G08_MD_598_proj[2],pos=4,label="598",col="blue")

x6G08_MD_599_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_599.txt")
x6G08_MD_599_proj<-project.pca(x6G08_MD_599_raw,pcadata)
points(x6G08_MD_599_proj[1],x6G08_MD_599_proj[2],pch=20)
text(x6G08_MD_599_proj[1],x6G08_MD_599_proj[2],pos=4,label="599",col="green")

x6G08_MD_600_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_600.txt")
x6G08_MD_600_proj<-project.pca(x6G08_MD_600_raw,pcadata)
points(x6G08_MD_600_proj[1],x6G08_MD_600_proj[2],pch=20)
text(x6G08_MD_600_proj[1],x6G08_MD_600_proj[2],pos=4,label="600",col="green")

x6G08_MD_601_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_601.txt")
x6G08_MD_601_proj<-project.pca(x6G08_MD_601_raw,pcadata)
points(x6G08_MD_601_proj[1],x6G08_MD_601_proj[2],pch=20)
text(x6G08_MD_601_proj[1],x6G08_MD_601_proj[2],pos=4,label="601",col="green")

x6G08_MD_602_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_602.txt")
x6G08_MD_602_proj<-project.pca(x6G08_MD_602_raw,pcadata)
points(x6G08_MD_602_proj[1],x6G08_MD_602_proj[2],pch=20)
text(x6G08_MD_602_proj[1],x6G08_MD_602_proj[2],pos=4,label="602",col="green")

x6G08_MD_603_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_603.txt")
x6G08_MD_603_proj<-project.pca(x6G08_MD_603_raw,pcadata)
points(x6G08_MD_603_proj[1],x6G08_MD_603_proj[2],pch=20)
text(x6G08_MD_603_proj[1],x6G08_MD_603_proj[2],pos=4,label="603",col="blue")

x6G08_MD_604_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_604.txt")
x6G08_MD_604_proj<-project.pca(x6G08_MD_604_raw,pcadata)
points(x6G08_MD_604_proj[1],x6G08_MD_604_proj[2],pch=20)
text(x6G08_MD_604_proj[1],x6G08_MD_604_proj[2],pos=4,label="604",col="black")

x6G08_MD_605_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_605.txt")
x6G08_MD_605_proj<-project.pca(x6G08_MD_605_raw,pcadata)
points(x6G08_MD_605_proj[1],x6G08_MD_605_proj[2],pch=20)
text(x6G08_MD_605_proj[1],x6G08_MD_605_proj[2],pos=4,label="605",col="blue")

x6G08_MD_606_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_606.txt")
x6G08_MD_606_proj<-project.pca(x6G08_MD_606_raw,pcadata)
points(x6G08_MD_606_proj[1],x6G08_MD_606_proj[2],pch=20)
text(x6G08_MD_606_proj[1],x6G08_MD_606_proj[2],pos=4,label="606",col="red")

x6G08_MD_607_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_607.txt")
x6G08_MD_607_proj<-project.pca(x6G08_MD_607_raw,pcadata)
points(x6G08_MD_607_proj[1],x6G08_MD_607_proj[2],pch=20)
text(x6G08_MD_607_proj[1],x6G08_MD_607_proj[2],pos=4,label="607",col="red")

x6G08_MD_608_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_608.txt")
x6G08_MD_608_proj<-project.pca(x6G08_MD_608_raw,pcadata)
points(x6G08_MD_608_proj[1],x6G08_MD_608_proj[2],pch=20)
text(x6G08_MD_608_proj[1],x6G08_MD_608_proj[2],pos=4,label="608",col="red")

x6G08_MD_609_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_609.txt")
x6G08_MD_609_proj<-project.pca(x6G08_MD_609_raw,pcadata)
points(x6G08_MD_609_proj[1],x6G08_MD_609_proj[2],pch=20)
text(x6G08_MD_609_proj[1],x6G08_MD_609_proj[2],pos=4,label="609",col="blue")

x6G08_MD_610_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_610.txt")
x6G08_MD_610_proj<-project.pca(x6G08_MD_610_raw,pcadata)
points(x6G08_MD_610_proj[1],x6G08_MD_610_proj[2],pch=20)
text(x6G08_MD_610_proj[1],x6G08_MD_610_proj[2],pos=4,label="610",col="black")

x6G08_MD_611_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_611.txt")
x6G08_MD_611_proj<-project.pca(x6G08_MD_611_raw,pcadata)
points(x6G08_MD_611_proj[1],x6G08_MD_611_proj[2],pch=20)
text(x6G08_MD_611_proj[1],x6G08_MD_611_proj[2],pos=4,label="611",col="blue")

x6G08_MD_612_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_612.txt")
x6G08_MD_612_proj<-project.pca(x6G08_MD_612_raw,pcadata)
points(x6G08_MD_612_proj[1],x6G08_MD_612_proj[2],pch=20)
text(x6G08_MD_612_proj[1],x6G08_MD_612_proj[2],pos=4,label="612",col="green")

x6G08_MD_613_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_613.txt")
x6G08_MD_613_proj<-project.pca(x6G08_MD_613_raw,pcadata)
points(x6G08_MD_613_proj[1],x6G08_MD_613_proj[2],pch=20)
text(x6G08_MD_613_proj[1],x6G08_MD_613_proj[2],pos=4,label="613",col="red")

x6G08_MD_614_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_614.txt")
x6G08_MD_614_proj<-project.pca(x6G08_MD_614_raw,pcadata)
points(x6G08_MD_614_proj[1],x6G08_MD_614_proj[2],pch=20)
text(x6G08_MD_614_proj[1],x6G08_MD_614_proj[2],pos=4,label="614",col="blue")

x6G08_MD_615_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_615.txt")
x6G08_MD_615_proj<-project.pca(x6G08_MD_615_raw,pcadata)
points(x6G08_MD_615_proj[1],x6G08_MD_615_proj[2],pch=20)
text(x6G08_MD_615_proj[1],x6G08_MD_615_proj[2],pos=4,label="615",col="blue")

x6G08_MD_616_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_616.txt")
x6G08_MD_616_proj<-project.pca(x6G08_MD_616_raw,pcadata)
points(x6G08_MD_616_proj[1],x6G08_MD_616_proj[2],pch=20)
text(x6G08_MD_616_proj[1],x6G08_MD_616_proj[2],pos=4,label="616",col="blue")

x6G08_MD_617_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_617.txt")
x6G08_MD_617_proj<-project.pca(x6G08_MD_617_raw,pcadata)
points(x6G08_MD_617_proj[1],x6G08_MD_617_proj[2],pch=20)
text(x6G08_MD_617_proj[1],x6G08_MD_617_proj[2],pos=4,label="617",col="green")

x6G08_MD_618_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_618.txt")
x6G08_MD_618_proj<-project.pca(x6G08_MD_618_raw,pcadata)
points(x6G08_MD_618_proj[1],x6G08_MD_618_proj[2],pch=20)
text(x6G08_MD_618_proj[1],x6G08_MD_618_proj[2],pos=4,label="618",col="green")

x6G08_MD_619_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_619.txt")
x6G08_MD_619_proj<-project.pca(x6G08_MD_619_raw,pcadata)
points(x6G08_MD_619_proj[1],x6G08_MD_619_proj[2],pch=20)
text(x6G08_MD_619_proj[1],x6G08_MD_619_proj[2],pos=4,label="619",col="blue")

x6G08_MD_620_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_620.txt")
x6G08_MD_620_proj<-project.pca(x6G08_MD_620_raw,pcadata)
points(x6G08_MD_620_proj[1],x6G08_MD_620_proj[2],pch=20)
text(x6G08_MD_620_proj[1],x6G08_MD_620_proj[2],pos=4,label="620",col="green")

x6G08_MD_621_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_621.txt")
x6G08_MD_621_proj<-project.pca(x6G08_MD_621_raw,pcadata)
points(x6G08_MD_621_proj[1],x6G08_MD_621_proj[2],pch=20)
text(x6G08_MD_621_proj[1],x6G08_MD_621_proj[2],pos=4,label="621",col="blue")

x6G08_MD_622_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_622.txt")
x6G08_MD_622_proj<-project.pca(x6G08_MD_622_raw,pcadata)
points(x6G08_MD_622_proj[1],x6G08_MD_622_proj[2],pch=20)
text(x6G08_MD_622_proj[1],x6G08_MD_622_proj[2],pos=4,label="622",col="blue")

x6G08_MD_623_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_623.txt")
x6G08_MD_623_proj<-project.pca(x6G08_MD_623_raw,pcadata)
points(x6G08_MD_623_proj[1],x6G08_MD_623_proj[2],pch=20)
text(x6G08_MD_623_proj[1],x6G08_MD_623_proj[2],pos=4,label="623",col="blue")

x6G08_MD_624_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_624.txt")
x6G08_MD_624_proj<-project.pca(x6G08_MD_624_raw,pcadata)
points(x6G08_MD_624_proj[1],x6G08_MD_624_proj[2],pch=20)
text(x6G08_MD_624_proj[1],x6G08_MD_624_proj[2],pos=4,label="624",col="blue")

x6G08_MD_625_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_625.txt")
x6G08_MD_625_proj<-project.pca(x6G08_MD_625_raw,pcadata)
points(x6G08_MD_625_proj[1],x6G08_MD_625_proj[2],pch=20)
text(x6G08_MD_625_proj[1],x6G08_MD_625_proj[2],pos=4,label="625",col="blue")

x6G08_MD_626_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_626.txt")
x6G08_MD_626_proj<-project.pca(x6G08_MD_626_raw,pcadata)
points(x6G08_MD_626_proj[1],x6G08_MD_626_proj[2],pch=20)
text(x6G08_MD_626_proj[1],x6G08_MD_626_proj[2],pos=4,label="626",col="blue")

x6G08_MD_627_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_627.txt")
x6G08_MD_627_proj<-project.pca(x6G08_MD_627_raw,pcadata)
points(x6G08_MD_627_proj[1],x6G08_MD_627_proj[2],pch=20)
text(x6G08_MD_627_proj[1],x6G08_MD_627_proj[2],pos=4,label="627",col="blue")

x6G08_MD_628_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_628.txt")
x6G08_MD_628_proj<-project.pca(x6G08_MD_628_raw,pcadata)
points(x6G08_MD_628_proj[1],x6G08_MD_628_proj[2],pch=20)
text(x6G08_MD_628_proj[1],x6G08_MD_628_proj[2],pos=4,label="628",col="red")

x6G08_MD_629_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_629.txt")
x6G08_MD_629_proj<-project.pca(x6G08_MD_629_raw,pcadata)
points(x6G08_MD_629_proj[1],x6G08_MD_629_proj[2],pch=20)
text(x6G08_MD_629_proj[1],x6G08_MD_629_proj[2],pos=4,label="629",col="red")

x6G08_MD_630_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_630.txt")
x6G08_MD_630_proj<-project.pca(x6G08_MD_630_raw,pcadata)
points(x6G08_MD_630_proj[1],x6G08_MD_630_proj[2],pch=20)
text(x6G08_MD_630_proj[1],x6G08_MD_630_proj[2],pos=4,label="630",col="blue")

x6G08_MD_631_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_631.txt")
x6G08_MD_631_proj<-project.pca(x6G08_MD_631_raw,pcadata)
points(x6G08_MD_631_proj[1],x6G08_MD_631_proj[2],pch=20)
text(x6G08_MD_631_proj[1],x6G08_MD_631_proj[2],pos=4,label="631",col="green")

x6G08_MD_632_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_632.txt")
x6G08_MD_632_proj<-project.pca(x6G08_MD_632_raw,pcadata)
points(x6G08_MD_632_proj[1],x6G08_MD_632_proj[2],pch=20)
text(x6G08_MD_632_proj[1],x6G08_MD_632_proj[2],pos=4,label="632",col="green")

x6G08_MD_633_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_633.txt")
x6G08_MD_633_proj<-project.pca(x6G08_MD_633_raw,pcadata)
points(x6G08_MD_633_proj[1],x6G08_MD_633_proj[2],pch=20)
text(x6G08_MD_633_proj[1],x6G08_MD_633_proj[2],pos=4,label="633",col="red")

x6G08_MD_634_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_634.txt")
x6G08_MD_634_proj<-project.pca(x6G08_MD_634_raw,pcadata)
points(x6G08_MD_634_proj[1],x6G08_MD_634_proj[2],pch=20)
text(x6G08_MD_634_proj[1],x6G08_MD_634_proj[2],pos=4,label="634",col="green")

x6G08_MD_635_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_635.txt")
x6G08_MD_635_proj<-project.pca(x6G08_MD_635_raw,pcadata)
points(x6G08_MD_635_proj[1],x6G08_MD_635_proj[2],pch=20)
text(x6G08_MD_635_proj[1],x6G08_MD_635_proj[2],pos=4,label="635",col="blue")

x6G08_MD_636_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_636.txt")
x6G08_MD_636_proj<-project.pca(x6G08_MD_636_raw,pcadata)
points(x6G08_MD_636_proj[1],x6G08_MD_636_proj[2],pch=20)
text(x6G08_MD_636_proj[1],x6G08_MD_636_proj[2],pos=4,label="636",col="black")

x6G08_MD_637_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_637.txt")
x6G08_MD_637_proj<-project.pca(x6G08_MD_637_raw,pcadata)
points(x6G08_MD_637_proj[1],x6G08_MD_637_proj[2],pch=20)
text(x6G08_MD_637_proj[1],x6G08_MD_637_proj[2],pos=4,label="637",col="blue")

x6G08_MD_638_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_638.txt")
x6G08_MD_638_proj<-project.pca(x6G08_MD_638_raw,pcadata)
points(x6G08_MD_638_proj[1],x6G08_MD_638_proj[2],pch=20)
text(x6G08_MD_638_proj[1],x6G08_MD_638_proj[2],pos=4,label="638",col="blue")

x6G08_MD_639_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_639.txt")
x6G08_MD_639_proj<-project.pca(x6G08_MD_639_raw,pcadata)
points(x6G08_MD_639_proj[1],x6G08_MD_639_proj[2],pch=20)
text(x6G08_MD_639_proj[1],x6G08_MD_639_proj[2],pos=4,label="639",col="blue")

x6G08_MD_640_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_640.txt")
x6G08_MD_640_proj<-project.pca(x6G08_MD_640_raw,pcadata)
points(x6G08_MD_640_proj[1],x6G08_MD_640_proj[2],pch=20)
text(x6G08_MD_640_proj[1],x6G08_MD_640_proj[2],pos=4,label="640",col="green")

x6G08_MD_641_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_641.txt")
x6G08_MD_641_proj<-project.pca(x6G08_MD_641_raw,pcadata)
points(x6G08_MD_641_proj[1],x6G08_MD_641_proj[2],pch=20)
text(x6G08_MD_641_proj[1],x6G08_MD_641_proj[2],pos=4,label="641",col="green")

x6G08_MD_642_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_642.txt")
x6G08_MD_642_proj<-project.pca(x6G08_MD_642_raw,pcadata)
points(x6G08_MD_642_proj[1],x6G08_MD_642_proj[2],pch=20)
text(x6G08_MD_642_proj[1],x6G08_MD_642_proj[2],pos=4,label="642",col="blue")

x6G08_MD_643_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_643.txt")
x6G08_MD_643_proj<-project.pca(x6G08_MD_643_raw,pcadata)
points(x6G08_MD_643_proj[1],x6G08_MD_643_proj[2],pch=20)
text(x6G08_MD_643_proj[1],x6G08_MD_643_proj[2],pos=4,label="643",col="blue")

x6G08_MD_644_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_644.txt")
x6G08_MD_644_proj<-project.pca(x6G08_MD_644_raw,pcadata)
points(x6G08_MD_644_proj[1],x6G08_MD_644_proj[2],pch=20)
text(x6G08_MD_644_proj[1],x6G08_MD_644_proj[2],pos=4,label="644",col="green")

x6G08_MD_645_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_645.txt")
x6G08_MD_645_proj<-project.pca(x6G08_MD_645_raw,pcadata)
points(x6G08_MD_645_proj[1],x6G08_MD_645_proj[2],pch=20)
text(x6G08_MD_645_proj[1],x6G08_MD_645_proj[2],pos=4,label="645",col="blue")

x6G08_MD_646_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_646.txt")
x6G08_MD_646_proj<-project.pca(x6G08_MD_646_raw,pcadata)
points(x6G08_MD_646_proj[1],x6G08_MD_646_proj[2],pch=20)
text(x6G08_MD_646_proj[1],x6G08_MD_646_proj[2],pos=4,label="646",col="blue")

x6G08_MD_647_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_647.txt")
x6G08_MD_647_proj<-project.pca(x6G08_MD_647_raw,pcadata)
points(x6G08_MD_647_proj[1],x6G08_MD_647_proj[2],pch=20)
text(x6G08_MD_647_proj[1],x6G08_MD_647_proj[2],pos=4,label="647",col="red")

x6G08_MD_648_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_648.txt")
x6G08_MD_648_proj<-project.pca(x6G08_MD_648_raw,pcadata)
points(x6G08_MD_648_proj[1],x6G08_MD_648_proj[2],pch=20)
text(x6G08_MD_648_proj[1],x6G08_MD_648_proj[2],pos=4,label="648",col="red")

x6G08_MD_649_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_649.txt")
x6G08_MD_649_proj<-project.pca(x6G08_MD_649_raw,pcadata)
points(x6G08_MD_649_proj[1],x6G08_MD_649_proj[2],pch=20)
text(x6G08_MD_649_proj[1],x6G08_MD_649_proj[2],pos=4,label="649",col="red")

x6G08_MD_650_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_650.txt")
x6G08_MD_650_proj<-project.pca(x6G08_MD_650_raw,pcadata)
points(x6G08_MD_650_proj[1],x6G08_MD_650_proj[2],pch=20)
text(x6G08_MD_650_proj[1],x6G08_MD_650_proj[2],pos=4,label="650",col="red")

x6G08_MD_651_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_651.txt")
x6G08_MD_651_proj<-project.pca(x6G08_MD_651_raw,pcadata)
points(x6G08_MD_651_proj[1],x6G08_MD_651_proj[2],pch=20)
text(x6G08_MD_651_proj[1],x6G08_MD_651_proj[2],pos=4,label="651",col="red")

x6G08_MD_652_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_652.txt")
x6G08_MD_652_proj<-project.pca(x6G08_MD_652_raw,pcadata)
points(x6G08_MD_652_proj[1],x6G08_MD_652_proj[2],pch=20)
text(x6G08_MD_652_proj[1],x6G08_MD_652_proj[2],pos=4,label="652",col="red")

x6G08_MD_653_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_653.txt")
x6G08_MD_653_proj<-project.pca(x6G08_MD_653_raw,pcadata)
points(x6G08_MD_653_proj[1],x6G08_MD_653_proj[2],pch=20)
text(x6G08_MD_653_proj[1],x6G08_MD_653_proj[2],pos=4,label="653",col="red")

x6G08_MD_654_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_654.txt")
x6G08_MD_654_proj<-project.pca(x6G08_MD_654_raw,pcadata)
points(x6G08_MD_654_proj[1],x6G08_MD_654_proj[2],pch=20)
text(x6G08_MD_654_proj[1],x6G08_MD_654_proj[2],pos=4,label="654",col="red")

x6G08_MD_655_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_655.txt")
x6G08_MD_655_proj<-project.pca(x6G08_MD_655_raw,pcadata)
points(x6G08_MD_655_proj[1],x6G08_MD_655_proj[2],pch=20)
text(x6G08_MD_655_proj[1],x6G08_MD_655_proj[2],pos=4,label="655",col="green")

x6G08_MD_656_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_656.txt")
x6G08_MD_656_proj<-project.pca(x6G08_MD_656_raw,pcadata)
points(x6G08_MD_656_proj[1],x6G08_MD_656_proj[2],pch=20)
text(x6G08_MD_656_proj[1],x6G08_MD_656_proj[2],pos=4,label="656",col="green")

x6G08_MD_657_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_657.txt")
x6G08_MD_657_proj<-project.pca(x6G08_MD_657_raw,pcadata)
points(x6G08_MD_657_proj[1],x6G08_MD_657_proj[2],pch=20)
text(x6G08_MD_657_proj[1],x6G08_MD_657_proj[2],pos=4,label="657",col="blue")

x6G08_MD_658_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_658.txt")
x6G08_MD_658_proj<-project.pca(x6G08_MD_658_raw,pcadata)
points(x6G08_MD_658_proj[1],x6G08_MD_658_proj[2],pch=20)
text(x6G08_MD_658_proj[1],x6G08_MD_658_proj[2],pos=4,label="658",col="green")

x6G08_MD_659_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_659.txt")
x6G08_MD_659_proj<-project.pca(x6G08_MD_659_raw,pcadata)
points(x6G08_MD_659_proj[1],x6G08_MD_659_proj[2],pch=20)
text(x6G08_MD_659_proj[1],x6G08_MD_659_proj[2],pos=4,label="659",col="green")

x6G08_MD_660_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_660.txt")
x6G08_MD_660_proj<-project.pca(x6G08_MD_660_raw,pcadata)
points(x6G08_MD_660_proj[1],x6G08_MD_660_proj[2],pch=20)
text(x6G08_MD_660_proj[1],x6G08_MD_660_proj[2],pos=4,label="660",col="red")

x6G08_MD_661_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_661.txt")
x6G08_MD_661_proj<-project.pca(x6G08_MD_661_raw,pcadata)
points(x6G08_MD_661_proj[1],x6G08_MD_661_proj[2],pch=20)
text(x6G08_MD_661_proj[1],x6G08_MD_661_proj[2],pos=4,label="661",col="blue")

x6G08_MD_662_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_662.txt")
x6G08_MD_662_proj<-project.pca(x6G08_MD_662_raw,pcadata)
points(x6G08_MD_662_proj[1],x6G08_MD_662_proj[2],pch=20)
text(x6G08_MD_662_proj[1],x6G08_MD_662_proj[2],pos=4,label="662",col="blue")

x6G08_MD_663_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_663.txt")
x6G08_MD_663_proj<-project.pca(x6G08_MD_663_raw,pcadata)
points(x6G08_MD_663_proj[1],x6G08_MD_663_proj[2],pch=20)
text(x6G08_MD_663_proj[1],x6G08_MD_663_proj[2],pos=4,label="663",col="green")

x6G08_MD_664_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_664.txt")
x6G08_MD_664_proj<-project.pca(x6G08_MD_664_raw,pcadata)
points(x6G08_MD_664_proj[1],x6G08_MD_664_proj[2],pch=20)
text(x6G08_MD_664_proj[1],x6G08_MD_664_proj[2],pos=4,label="664",col="red")

x6G08_MD_665_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_665.txt")
x6G08_MD_665_proj<-project.pca(x6G08_MD_665_raw,pcadata)
points(x6G08_MD_665_proj[1],x6G08_MD_665_proj[2],pch=20)
text(x6G08_MD_665_proj[1],x6G08_MD_665_proj[2],pos=4,label="665",col="red")

x6G08_MD_666_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_666.txt")
x6G08_MD_666_proj<-project.pca(x6G08_MD_666_raw,pcadata)
points(x6G08_MD_666_proj[1],x6G08_MD_666_proj[2],pch=20)
text(x6G08_MD_666_proj[1],x6G08_MD_666_proj[2],pos=4,label="666",col="red")

x6G08_MD_667_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_667.txt")
x6G08_MD_667_proj<-project.pca(x6G08_MD_667_raw,pcadata)
points(x6G08_MD_667_proj[1],x6G08_MD_667_proj[2],pch=20)
text(x6G08_MD_667_proj[1],x6G08_MD_667_proj[2],pos=4,label="667",col="red")

x6G08_MD_668_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_668.txt")
x6G08_MD_668_proj<-project.pca(x6G08_MD_668_raw,pcadata)
points(x6G08_MD_668_proj[1],x6G08_MD_668_proj[2],pch=20)
text(x6G08_MD_668_proj[1],x6G08_MD_668_proj[2],pos=4,label="668",col="red")

x6G08_MD_669_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_669.txt")
x6G08_MD_669_proj<-project.pca(x6G08_MD_669_raw,pcadata)
points(x6G08_MD_669_proj[1],x6G08_MD_669_proj[2],pch=20)
text(x6G08_MD_669_proj[1],x6G08_MD_669_proj[2],pos=4,label="669",col="red")

x6G08_MD_670_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_670.txt")
x6G08_MD_670_proj<-project.pca(x6G08_MD_670_raw,pcadata)
points(x6G08_MD_670_proj[1],x6G08_MD_670_proj[2],pch=20)
text(x6G08_MD_670_proj[1],x6G08_MD_670_proj[2],pos=4,label="670",col="red")

x6G08_MD_671_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_671.txt")
x6G08_MD_671_proj<-project.pca(x6G08_MD_671_raw,pcadata)
points(x6G08_MD_671_proj[1],x6G08_MD_671_proj[2],pch=20)
text(x6G08_MD_671_proj[1],x6G08_MD_671_proj[2],pos=4,label="671",col="red")

x6G08_MD_672_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_672.txt")
x6G08_MD_672_proj<-project.pca(x6G08_MD_672_raw,pcadata)
points(x6G08_MD_672_proj[1],x6G08_MD_672_proj[2],pch=20)
text(x6G08_MD_672_proj[1],x6G08_MD_672_proj[2],pos=4,label="672",col="red")

x6G08_MD_673_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_673.txt")
x6G08_MD_673_proj<-project.pca(x6G08_MD_673_raw,pcadata)
points(x6G08_MD_673_proj[1],x6G08_MD_673_proj[2],pch=20)
text(x6G08_MD_673_proj[1],x6G08_MD_673_proj[2],pos=4,label="673",col="red")

x6G08_MD_674_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_674.txt")
x6G08_MD_674_proj<-project.pca(x6G08_MD_674_raw,pcadata)
points(x6G08_MD_674_proj[1],x6G08_MD_674_proj[2],pch=20)
text(x6G08_MD_674_proj[1],x6G08_MD_674_proj[2],pos=4,label="674",col="red")

x6G08_MD_675_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_675.txt")
x6G08_MD_675_proj<-project.pca(x6G08_MD_675_raw,pcadata)
points(x6G08_MD_675_proj[1],x6G08_MD_675_proj[2],pch=20)
text(x6G08_MD_675_proj[1],x6G08_MD_675_proj[2],pos=4,label="675",col="red")

x6G08_MD_676_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_676.txt")
x6G08_MD_676_proj<-project.pca(x6G08_MD_676_raw,pcadata)
points(x6G08_MD_676_proj[1],x6G08_MD_676_proj[2],pch=20)
text(x6G08_MD_676_proj[1],x6G08_MD_676_proj[2],pos=4,label="676",col="red")

x6G08_MD_677_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_677.txt")
x6G08_MD_677_proj<-project.pca(x6G08_MD_677_raw,pcadata)
points(x6G08_MD_677_proj[1],x6G08_MD_677_proj[2],pch=20)
text(x6G08_MD_677_proj[1],x6G08_MD_677_proj[2],pos=4,label="677",col="red")

x6G08_MD_678_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_678.txt")
x6G08_MD_678_proj<-project.pca(x6G08_MD_678_raw,pcadata)
points(x6G08_MD_678_proj[1],x6G08_MD_678_proj[2],pch=20)
text(x6G08_MD_678_proj[1],x6G08_MD_678_proj[2],pos=4,label="678",col="red")

x6G08_MD_679_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_679.txt")
x6G08_MD_679_proj<-project.pca(x6G08_MD_679_raw,pcadata)
points(x6G08_MD_679_proj[1],x6G08_MD_679_proj[2],pch=20)
text(x6G08_MD_679_proj[1],x6G08_MD_679_proj[2],pos=4,label="679",col="red")

x6G08_MD_680_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_680.txt")
x6G08_MD_680_proj<-project.pca(x6G08_MD_680_raw,pcadata)
points(x6G08_MD_680_proj[1],x6G08_MD_680_proj[2],pch=20)
text(x6G08_MD_680_proj[1],x6G08_MD_680_proj[2],pos=4,label="680",col="red")

x6G08_MD_681_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_681.txt")
x6G08_MD_681_proj<-project.pca(x6G08_MD_681_raw,pcadata)
points(x6G08_MD_681_proj[1],x6G08_MD_681_proj[2],pch=20)
text(x6G08_MD_681_proj[1],x6G08_MD_681_proj[2],pos=4,label="681",col="red")

x6G08_MD_682_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_682.txt")
x6G08_MD_682_proj<-project.pca(x6G08_MD_682_raw,pcadata)
points(x6G08_MD_682_proj[1],x6G08_MD_682_proj[2],pch=20)
text(x6G08_MD_682_proj[1],x6G08_MD_682_proj[2],pos=4,label="682",col="red")

x6G08_MD_683_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_683.txt")
x6G08_MD_683_proj<-project.pca(x6G08_MD_683_raw,pcadata)
points(x6G08_MD_683_proj[1],x6G08_MD_683_proj[2],pch=20)
text(x6G08_MD_683_proj[1],x6G08_MD_683_proj[2],pos=4,label="683",col="red")

x6G08_MD_684_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_684.txt")
x6G08_MD_684_proj<-project.pca(x6G08_MD_684_raw,pcadata)
points(x6G08_MD_684_proj[1],x6G08_MD_684_proj[2],pch=20)
text(x6G08_MD_684_proj[1],x6G08_MD_684_proj[2],pos=4,label="684",col="red")

x6G08_MD_685_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_685.txt")
x6G08_MD_685_proj<-project.pca(x6G08_MD_685_raw,pcadata)
points(x6G08_MD_685_proj[1],x6G08_MD_685_proj[2],pch=20)
text(x6G08_MD_685_proj[1],x6G08_MD_685_proj[2],pos=4,label="685",col="red")

x6G08_MD_686_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_686.txt")
x6G08_MD_686_proj<-project.pca(x6G08_MD_686_raw,pcadata)
points(x6G08_MD_686_proj[1],x6G08_MD_686_proj[2],pch=20)
text(x6G08_MD_686_proj[1],x6G08_MD_686_proj[2],pos=4,label="686",col="red")

x6G08_MD_687_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_687.txt")
x6G08_MD_687_proj<-project.pca(x6G08_MD_687_raw,pcadata)
points(x6G08_MD_687_proj[1],x6G08_MD_687_proj[2],pch=20)
text(x6G08_MD_687_proj[1],x6G08_MD_687_proj[2],pos=4,label="687",col="red")

x6G08_MD_688_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_688.txt")
x6G08_MD_688_proj<-project.pca(x6G08_MD_688_raw,pcadata)
points(x6G08_MD_688_proj[1],x6G08_MD_688_proj[2],pch=20)
text(x6G08_MD_688_proj[1],x6G08_MD_688_proj[2],pos=4,label="688",col="red")

x6G08_MD_689_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_689.txt")
x6G08_MD_689_proj<-project.pca(x6G08_MD_689_raw,pcadata)
points(x6G08_MD_689_proj[1],x6G08_MD_689_proj[2],pch=20)
text(x6G08_MD_689_proj[1],x6G08_MD_689_proj[2],pos=4,label="689",col="red")

x6G08_MD_690_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_690.txt")
x6G08_MD_690_proj<-project.pca(x6G08_MD_690_raw,pcadata)
points(x6G08_MD_690_proj[1],x6G08_MD_690_proj[2],pch=20)
text(x6G08_MD_690_proj[1],x6G08_MD_690_proj[2],pos=4,label="690",col="red")

x6G08_MD_691_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_691.txt")
x6G08_MD_691_proj<-project.pca(x6G08_MD_691_raw,pcadata)
points(x6G08_MD_691_proj[1],x6G08_MD_691_proj[2],pch=20)
text(x6G08_MD_691_proj[1],x6G08_MD_691_proj[2],pos=4,label="691",col="red")

x6G08_MD_692_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_692.txt")
x6G08_MD_692_proj<-project.pca(x6G08_MD_692_raw,pcadata)
points(x6G08_MD_692_proj[1],x6G08_MD_692_proj[2],pch=20)
text(x6G08_MD_692_proj[1],x6G08_MD_692_proj[2],pos=4,label="692",col="green")

x6G08_MD_693_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_693.txt")
x6G08_MD_693_proj<-project.pca(x6G08_MD_693_raw,pcadata)
points(x6G08_MD_693_proj[1],x6G08_MD_693_proj[2],pch=20)
text(x6G08_MD_693_proj[1],x6G08_MD_693_proj[2],pos=4,label="693",col="red")

x6G08_MD_694_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_694.txt")
x6G08_MD_694_proj<-project.pca(x6G08_MD_694_raw,pcadata)
points(x6G08_MD_694_proj[1],x6G08_MD_694_proj[2],pch=20)
text(x6G08_MD_694_proj[1],x6G08_MD_694_proj[2],pos=4,label="694",col="red")

x6G08_MD_695_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_695.txt")
x6G08_MD_695_proj<-project.pca(x6G08_MD_695_raw,pcadata)
points(x6G08_MD_695_proj[1],x6G08_MD_695_proj[2],pch=20)
text(x6G08_MD_695_proj[1],x6G08_MD_695_proj[2],pos=4,label="695",col="red")

x6G08_MD_696_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_696.txt")
x6G08_MD_696_proj<-project.pca(x6G08_MD_696_raw,pcadata)
points(x6G08_MD_696_proj[1],x6G08_MD_696_proj[2],pch=20)
text(x6G08_MD_696_proj[1],x6G08_MD_696_proj[2],pos=4,label="696",col="red")

x6G08_MD_697_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_697.txt")
x6G08_MD_697_proj<-project.pca(x6G08_MD_697_raw,pcadata)
points(x6G08_MD_697_proj[1],x6G08_MD_697_proj[2],pch=20)
text(x6G08_MD_697_proj[1],x6G08_MD_697_proj[2],pos=4,label="697",col="red")

x6G08_MD_698_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_698.txt")
x6G08_MD_698_proj<-project.pca(x6G08_MD_698_raw,pcadata)
points(x6G08_MD_698_proj[1],x6G08_MD_698_proj[2],pch=20)
text(x6G08_MD_698_proj[1],x6G08_MD_698_proj[2],pos=4,label="698",col="red")

x6G08_MD_699_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_699.txt")
x6G08_MD_699_proj<-project.pca(x6G08_MD_699_raw,pcadata)
points(x6G08_MD_699_proj[1],x6G08_MD_699_proj[2],pch=20)
text(x6G08_MD_699_proj[1],x6G08_MD_699_proj[2],pos=4,label="699",col="red")

x6G08_MD_700_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_700.txt")
x6G08_MD_700_proj<-project.pca(x6G08_MD_700_raw,pcadata)
points(x6G08_MD_700_proj[1],x6G08_MD_700_proj[2],pch=20)
text(x6G08_MD_700_proj[1],x6G08_MD_700_proj[2],pos=4,label="700",col="red")

x6G08_MD_701_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_701.txt")
x6G08_MD_701_proj<-project.pca(x6G08_MD_701_raw,pcadata)
points(x6G08_MD_701_proj[1],x6G08_MD_701_proj[2],pch=20)
text(x6G08_MD_701_proj[1],x6G08_MD_701_proj[2],pos=4,label="701",col="red")

x6G08_MD_702_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_702.txt")
x6G08_MD_702_proj<-project.pca(x6G08_MD_702_raw,pcadata)
points(x6G08_MD_702_proj[1],x6G08_MD_702_proj[2],pch=20)
text(x6G08_MD_702_proj[1],x6G08_MD_702_proj[2],pos=4,label="702",col="green")

x6G08_MD_703_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_703.txt")
x6G08_MD_703_proj<-project.pca(x6G08_MD_703_raw,pcadata)
points(x6G08_MD_703_proj[1],x6G08_MD_703_proj[2],pch=20)
text(x6G08_MD_703_proj[1],x6G08_MD_703_proj[2],pos=4,label="703",col="red")

x6G08_MD_704_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_704.txt")
x6G08_MD_704_proj<-project.pca(x6G08_MD_704_raw,pcadata)
points(x6G08_MD_704_proj[1],x6G08_MD_704_proj[2],pch=20)
text(x6G08_MD_704_proj[1],x6G08_MD_704_proj[2],pos=4,label="704",col="red")

x6G08_MD_705_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_705.txt")
x6G08_MD_705_proj<-project.pca(x6G08_MD_705_raw,pcadata)
points(x6G08_MD_705_proj[1],x6G08_MD_705_proj[2],pch=20)
text(x6G08_MD_705_proj[1],x6G08_MD_705_proj[2],pos=4,label="705",col="red")

x6G08_MD_706_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_706.txt")
x6G08_MD_706_proj<-project.pca(x6G08_MD_706_raw,pcadata)
points(x6G08_MD_706_proj[1],x6G08_MD_706_proj[2],pch=20)
text(x6G08_MD_706_proj[1],x6G08_MD_706_proj[2],pos=4,label="706",col="red")

x6G08_MD_707_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_707.txt")
x6G08_MD_707_proj<-project.pca(x6G08_MD_707_raw,pcadata)
points(x6G08_MD_707_proj[1],x6G08_MD_707_proj[2],pch=20)
text(x6G08_MD_707_proj[1],x6G08_MD_707_proj[2],pos=4,label="707",col="red")

x6G08_MD_708_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_708.txt")
x6G08_MD_708_proj<-project.pca(x6G08_MD_708_raw,pcadata)
points(x6G08_MD_708_proj[1],x6G08_MD_708_proj[2],pch=20)
text(x6G08_MD_708_proj[1],x6G08_MD_708_proj[2],pos=4,label="708",col="red")

x6G08_MD_709_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_709.txt")
x6G08_MD_709_proj<-project.pca(x6G08_MD_709_raw,pcadata)
points(x6G08_MD_709_proj[1],x6G08_MD_709_proj[2],pch=20)
text(x6G08_MD_709_proj[1],x6G08_MD_709_proj[2],pos=4,label="709",col="red")

x6G08_MD_710_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_710.txt")
x6G08_MD_710_proj<-project.pca(x6G08_MD_710_raw,pcadata)
points(x6G08_MD_710_proj[1],x6G08_MD_710_proj[2],pch=20)
text(x6G08_MD_710_proj[1],x6G08_MD_710_proj[2],pos=4,label="710",col="red")

x6G08_MD_711_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_711.txt")
x6G08_MD_711_proj<-project.pca(x6G08_MD_711_raw,pcadata)
points(x6G08_MD_711_proj[1],x6G08_MD_711_proj[2],pch=20)
text(x6G08_MD_711_proj[1],x6G08_MD_711_proj[2],pos=4,label="711",col="red")

x6G08_MD_712_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_712.txt")
x6G08_MD_712_proj<-project.pca(x6G08_MD_712_raw,pcadata)
points(x6G08_MD_712_proj[1],x6G08_MD_712_proj[2],pch=20)
text(x6G08_MD_712_proj[1],x6G08_MD_712_proj[2],pos=4,label="712",col="red")

x6G08_MD_713_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_713.txt")
x6G08_MD_713_proj<-project.pca(x6G08_MD_713_raw,pcadata)
points(x6G08_MD_713_proj[1],x6G08_MD_713_proj[2],pch=20)
text(x6G08_MD_713_proj[1],x6G08_MD_713_proj[2],pos=4,label="713",col="red")

x6G08_MD_714_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_714.txt")
x6G08_MD_714_proj<-project.pca(x6G08_MD_714_raw,pcadata)
points(x6G08_MD_714_proj[1],x6G08_MD_714_proj[2],pch=20)
text(x6G08_MD_714_proj[1],x6G08_MD_714_proj[2],pos=4,label="714",col="red")

x6G08_MD_715_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_715.txt")
x6G08_MD_715_proj<-project.pca(x6G08_MD_715_raw,pcadata)
points(x6G08_MD_715_proj[1],x6G08_MD_715_proj[2],pch=20)
text(x6G08_MD_715_proj[1],x6G08_MD_715_proj[2],pos=4,label="715",col="red")

x6G08_MD_716_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_716.txt")
x6G08_MD_716_proj<-project.pca(x6G08_MD_716_raw,pcadata)
points(x6G08_MD_716_proj[1],x6G08_MD_716_proj[2],pch=20)
text(x6G08_MD_716_proj[1],x6G08_MD_716_proj[2],pos=4,label="716",col="red")

x6G08_MD_717_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_717.txt")
x6G08_MD_717_proj<-project.pca(x6G08_MD_717_raw,pcadata)
points(x6G08_MD_717_proj[1],x6G08_MD_717_proj[2],pch=20)
text(x6G08_MD_717_proj[1],x6G08_MD_717_proj[2],pos=4,label="717",col="red")

x6G08_MD_718_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_718.txt")
x6G08_MD_718_proj<-project.pca(x6G08_MD_718_raw,pcadata)
points(x6G08_MD_718_proj[1],x6G08_MD_718_proj[2],pch=20)
text(x6G08_MD_718_proj[1],x6G08_MD_718_proj[2],pos=4,label="718",col="red")

x6G08_MD_719_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_719.txt")
x6G08_MD_719_proj<-project.pca(x6G08_MD_719_raw,pcadata)
points(x6G08_MD_719_proj[1],x6G08_MD_719_proj[2],pch=20)
text(x6G08_MD_719_proj[1],x6G08_MD_719_proj[2],pos=4,label="719",col="red")

x6G08_MD_720_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_720.txt")
x6G08_MD_720_proj<-project.pca(x6G08_MD_720_raw,pcadata)
points(x6G08_MD_720_proj[1],x6G08_MD_720_proj[2],pch=20)
text(x6G08_MD_720_proj[1],x6G08_MD_720_proj[2],pos=4,label="720",col="red")

x6G08_MD_721_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_721.txt")
x6G08_MD_721_proj<-project.pca(x6G08_MD_721_raw,pcadata)
points(x6G08_MD_721_proj[1],x6G08_MD_721_proj[2],pch=20)
text(x6G08_MD_721_proj[1],x6G08_MD_721_proj[2],pos=4,label="721",col="red")

x6G08_MD_722_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_722.txt")
x6G08_MD_722_proj<-project.pca(x6G08_MD_722_raw,pcadata)
points(x6G08_MD_722_proj[1],x6G08_MD_722_proj[2],pch=20)
text(x6G08_MD_722_proj[1],x6G08_MD_722_proj[2],pos=4,label="722",col="red")

x6G08_MD_723_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_723.txt")
x6G08_MD_723_proj<-project.pca(x6G08_MD_723_raw,pcadata)
points(x6G08_MD_723_proj[1],x6G08_MD_723_proj[2],pch=20)
text(x6G08_MD_723_proj[1],x6G08_MD_723_proj[2],pos=4,label="723",col="red")

x6G08_MD_724_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_724.txt")
x6G08_MD_724_proj<-project.pca(x6G08_MD_724_raw,pcadata)
points(x6G08_MD_724_proj[1],x6G08_MD_724_proj[2],pch=20)
text(x6G08_MD_724_proj[1],x6G08_MD_724_proj[2],pos=4,label="724",col="red")

x6G08_MD_725_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_725.txt")
x6G08_MD_725_proj<-project.pca(x6G08_MD_725_raw,pcadata)
points(x6G08_MD_725_proj[1],x6G08_MD_725_proj[2],pch=20)
text(x6G08_MD_725_proj[1],x6G08_MD_725_proj[2],pos=4,label="725",col="red")

x6G08_MD_726_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_726.txt")
x6G08_MD_726_proj<-project.pca(x6G08_MD_726_raw,pcadata)
points(x6G08_MD_726_proj[1],x6G08_MD_726_proj[2],pch=20)
text(x6G08_MD_726_proj[1],x6G08_MD_726_proj[2],pos=4,label="726",col="red")

x6G08_MD_727_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_727.txt")
x6G08_MD_727_proj<-project.pca(x6G08_MD_727_raw,pcadata)
points(x6G08_MD_727_proj[1],x6G08_MD_727_proj[2],pch=20)
text(x6G08_MD_727_proj[1],x6G08_MD_727_proj[2],pos=4,label="727",col="red")

x6G08_MD_728_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_728.txt")
x6G08_MD_728_proj<-project.pca(x6G08_MD_728_raw,pcadata)
points(x6G08_MD_728_proj[1],x6G08_MD_728_proj[2],pch=20)
text(x6G08_MD_728_proj[1],x6G08_MD_728_proj[2],pos=4,label="728",col="green")

x6G08_MD_729_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_729.txt")
x6G08_MD_729_proj<-project.pca(x6G08_MD_729_raw,pcadata)
points(x6G08_MD_729_proj[1],x6G08_MD_729_proj[2],pch=20)
text(x6G08_MD_729_proj[1],x6G08_MD_729_proj[2],pos=4,label="729",col="red")

x6G08_MD_730_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_730.txt")
x6G08_MD_730_proj<-project.pca(x6G08_MD_730_raw,pcadata)
points(x6G08_MD_730_proj[1],x6G08_MD_730_proj[2],pch=20)
text(x6G08_MD_730_proj[1],x6G08_MD_730_proj[2],pos=4,label="730",col="red")

x6G08_MD_731_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_731.txt")
x6G08_MD_731_proj<-project.pca(x6G08_MD_731_raw,pcadata)
points(x6G08_MD_731_proj[1],x6G08_MD_731_proj[2],pch=20)
text(x6G08_MD_731_proj[1],x6G08_MD_731_proj[2],pos=4,label="731",col="red")

x6G08_MD_732_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_732.txt")
x6G08_MD_732_proj<-project.pca(x6G08_MD_732_raw,pcadata)
points(x6G08_MD_732_proj[1],x6G08_MD_732_proj[2],pch=20)
text(x6G08_MD_732_proj[1],x6G08_MD_732_proj[2],pos=4,label="732",col="red")

x6G08_MD_733_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_733.txt")
x6G08_MD_733_proj<-project.pca(x6G08_MD_733_raw,pcadata)
points(x6G08_MD_733_proj[1],x6G08_MD_733_proj[2],pch=20)
text(x6G08_MD_733_proj[1],x6G08_MD_733_proj[2],pos=4,label="733",col="red")

x6G08_MD_734_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_734.txt")
x6G08_MD_734_proj<-project.pca(x6G08_MD_734_raw,pcadata)
points(x6G08_MD_734_proj[1],x6G08_MD_734_proj[2],pch=20)
text(x6G08_MD_734_proj[1],x6G08_MD_734_proj[2],pos=4,label="734",col="red")

x6G08_MD_735_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_735.txt")
x6G08_MD_735_proj<-project.pca(x6G08_MD_735_raw,pcadata)
points(x6G08_MD_735_proj[1],x6G08_MD_735_proj[2],pch=20)
text(x6G08_MD_735_proj[1],x6G08_MD_735_proj[2],pos=4,label="735",col="red")

x6G08_MD_736_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_736.txt")
x6G08_MD_736_proj<-project.pca(x6G08_MD_736_raw,pcadata)
points(x6G08_MD_736_proj[1],x6G08_MD_736_proj[2],pch=20)
text(x6G08_MD_736_proj[1],x6G08_MD_736_proj[2],pos=4,label="736",col="red")

x6G08_MD_737_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_737.txt")
x6G08_MD_737_proj<-project.pca(x6G08_MD_737_raw,pcadata)
points(x6G08_MD_737_proj[1],x6G08_MD_737_proj[2],pch=20)
text(x6G08_MD_737_proj[1],x6G08_MD_737_proj[2],pos=4,label="737",col="red")

x6G08_MD_738_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_738.txt")
x6G08_MD_738_proj<-project.pca(x6G08_MD_738_raw,pcadata)
points(x6G08_MD_738_proj[1],x6G08_MD_738_proj[2],pch=20)
text(x6G08_MD_738_proj[1],x6G08_MD_738_proj[2],pos=4,label="738",col="red")

x6G08_MD_739_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_739.txt")
x6G08_MD_739_proj<-project.pca(x6G08_MD_739_raw,pcadata)
points(x6G08_MD_739_proj[1],x6G08_MD_739_proj[2],pch=20)
text(x6G08_MD_739_proj[1],x6G08_MD_739_proj[2],pos=4,label="739",col="red")

x6G08_MD_740_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_740.txt")
x6G08_MD_740_proj<-project.pca(x6G08_MD_740_raw,pcadata)
points(x6G08_MD_740_proj[1],x6G08_MD_740_proj[2],pch=20)
text(x6G08_MD_740_proj[1],x6G08_MD_740_proj[2],pos=4,label="740",col="red")

x6G08_MD_741_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_741.txt")
x6G08_MD_741_proj<-project.pca(x6G08_MD_741_raw,pcadata)
points(x6G08_MD_741_proj[1],x6G08_MD_741_proj[2],pch=20)
text(x6G08_MD_741_proj[1],x6G08_MD_741_proj[2],pos=4,label="741",col="red")

x6G08_MD_742_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_742.txt")
x6G08_MD_742_proj<-project.pca(x6G08_MD_742_raw,pcadata)
points(x6G08_MD_742_proj[1],x6G08_MD_742_proj[2],pch=20)
text(x6G08_MD_742_proj[1],x6G08_MD_742_proj[2],pos=4,label="742",col="red")

x6G08_MD_743_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_743.txt")
x6G08_MD_743_proj<-project.pca(x6G08_MD_743_raw,pcadata)
points(x6G08_MD_743_proj[1],x6G08_MD_743_proj[2],pch=20)
text(x6G08_MD_743_proj[1],x6G08_MD_743_proj[2],pos=4,label="743",col="red")

x6G08_MD_744_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_744.txt")
x6G08_MD_744_proj<-project.pca(x6G08_MD_744_raw,pcadata)
points(x6G08_MD_744_proj[1],x6G08_MD_744_proj[2],pch=20)
text(x6G08_MD_744_proj[1],x6G08_MD_744_proj[2],pos=4,label="744",col="red")

x6G08_MD_745_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_745.txt")
x6G08_MD_745_proj<-project.pca(x6G08_MD_745_raw,pcadata)
points(x6G08_MD_745_proj[1],x6G08_MD_745_proj[2],pch=20)
text(x6G08_MD_745_proj[1],x6G08_MD_745_proj[2],pos=4,label="745",col="red")

x6G08_MD_746_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_746.txt")
x6G08_MD_746_proj<-project.pca(x6G08_MD_746_raw,pcadata)
points(x6G08_MD_746_proj[1],x6G08_MD_746_proj[2],pch=20)
text(x6G08_MD_746_proj[1],x6G08_MD_746_proj[2],pos=4,label="746",col="red")

x6G08_MD_747_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_747.txt")
x6G08_MD_747_proj<-project.pca(x6G08_MD_747_raw,pcadata)
points(x6G08_MD_747_proj[1],x6G08_MD_747_proj[2],pch=20)
text(x6G08_MD_747_proj[1],x6G08_MD_747_proj[2],pos=4,label="747",col="red")

x6G08_MD_748_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_748.txt")
x6G08_MD_748_proj<-project.pca(x6G08_MD_748_raw,pcadata)
points(x6G08_MD_748_proj[1],x6G08_MD_748_proj[2],pch=20)
text(x6G08_MD_748_proj[1],x6G08_MD_748_proj[2],pos=4,label="748",col="blue")

x6G08_MD_749_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_749.txt")
x6G08_MD_749_proj<-project.pca(x6G08_MD_749_raw,pcadata)
points(x6G08_MD_749_proj[1],x6G08_MD_749_proj[2],pch=20)
text(x6G08_MD_749_proj[1],x6G08_MD_749_proj[2],pos=4,label="749",col="red")

x6G08_MD_750_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_750.txt")
x6G08_MD_750_proj<-project.pca(x6G08_MD_750_raw,pcadata)
points(x6G08_MD_750_proj[1],x6G08_MD_750_proj[2],pch=20)
text(x6G08_MD_750_proj[1],x6G08_MD_750_proj[2],pos=4,label="750",col="red")

x6G08_MD_751_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_751.txt")
x6G08_MD_751_proj<-project.pca(x6G08_MD_751_raw,pcadata)
points(x6G08_MD_751_proj[1],x6G08_MD_751_proj[2],pch=20)
text(x6G08_MD_751_proj[1],x6G08_MD_751_proj[2],pos=4,label="751",col="red")

x6G08_MD_752_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_752.txt")
x6G08_MD_752_proj<-project.pca(x6G08_MD_752_raw,pcadata)
points(x6G08_MD_752_proj[1],x6G08_MD_752_proj[2],pch=20)
text(x6G08_MD_752_proj[1],x6G08_MD_752_proj[2],pos=4,label="752",col="red")

x6G08_MD_753_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_753.txt")
x6G08_MD_753_proj<-project.pca(x6G08_MD_753_raw,pcadata)
points(x6G08_MD_753_proj[1],x6G08_MD_753_proj[2],pch=20)
text(x6G08_MD_753_proj[1],x6G08_MD_753_proj[2],pos=4,label="753",col="red")

x6G08_MD_754_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_754.txt")
x6G08_MD_754_proj<-project.pca(x6G08_MD_754_raw,pcadata)
points(x6G08_MD_754_proj[1],x6G08_MD_754_proj[2],pch=20)
text(x6G08_MD_754_proj[1],x6G08_MD_754_proj[2],pos=4,label="754",col="green")

x6G08_MD_755_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_755.txt")
x6G08_MD_755_proj<-project.pca(x6G08_MD_755_raw,pcadata)
points(x6G08_MD_755_proj[1],x6G08_MD_755_proj[2],pch=20)
text(x6G08_MD_755_proj[1],x6G08_MD_755_proj[2],pos=4,label="755",col="red")

x6G08_MD_756_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_756.txt")
x6G08_MD_756_proj<-project.pca(x6G08_MD_756_raw,pcadata)
points(x6G08_MD_756_proj[1],x6G08_MD_756_proj[2],pch=20)
text(x6G08_MD_756_proj[1],x6G08_MD_756_proj[2],pos=4,label="756",col="red")

x6G08_MD_757_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_757.txt")
x6G08_MD_757_proj<-project.pca(x6G08_MD_757_raw,pcadata)
points(x6G08_MD_757_proj[1],x6G08_MD_757_proj[2],pch=20)
text(x6G08_MD_757_proj[1],x6G08_MD_757_proj[2],pos=4,label="757",col="red")

x6G08_MD_758_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_758.txt")
x6G08_MD_758_proj<-project.pca(x6G08_MD_758_raw,pcadata)
points(x6G08_MD_758_proj[1],x6G08_MD_758_proj[2],pch=20)
text(x6G08_MD_758_proj[1],x6G08_MD_758_proj[2],pos=4,label="758",col="red")

x6G08_MD_759_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_759.txt")
x6G08_MD_759_proj<-project.pca(x6G08_MD_759_raw,pcadata)
points(x6G08_MD_759_proj[1],x6G08_MD_759_proj[2],pch=20)
text(x6G08_MD_759_proj[1],x6G08_MD_759_proj[2],pos=4,label="759",col="red")

x6G08_MD_760_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_760.txt")
x6G08_MD_760_proj<-project.pca(x6G08_MD_760_raw,pcadata)
points(x6G08_MD_760_proj[1],x6G08_MD_760_proj[2],pch=20)
text(x6G08_MD_760_proj[1],x6G08_MD_760_proj[2],pos=4,label="760",col="red")

x6G08_MD_761_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_761.txt")
x6G08_MD_761_proj<-project.pca(x6G08_MD_761_raw,pcadata)
points(x6G08_MD_761_proj[1],x6G08_MD_761_proj[2],pch=20)
text(x6G08_MD_761_proj[1],x6G08_MD_761_proj[2],pos=4,label="761",col="red")

x6G08_MD_762_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_762.txt")
x6G08_MD_762_proj<-project.pca(x6G08_MD_762_raw,pcadata)
points(x6G08_MD_762_proj[1],x6G08_MD_762_proj[2],pch=20)
text(x6G08_MD_762_proj[1],x6G08_MD_762_proj[2],pos=4,label="762",col="red")

x6G08_MD_763_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_763.txt")
x6G08_MD_763_proj<-project.pca(x6G08_MD_763_raw,pcadata)
points(x6G08_MD_763_proj[1],x6G08_MD_763_proj[2],pch=20)
text(x6G08_MD_763_proj[1],x6G08_MD_763_proj[2],pos=4,label="763",col="red")

x6G08_MD_764_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_764.txt")
x6G08_MD_764_proj<-project.pca(x6G08_MD_764_raw,pcadata)
points(x6G08_MD_764_proj[1],x6G08_MD_764_proj[2],pch=20)
text(x6G08_MD_764_proj[1],x6G08_MD_764_proj[2],pos=4,label="764",col="green")

x6G08_MD_765_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_765.txt")
x6G08_MD_765_proj<-project.pca(x6G08_MD_765_raw,pcadata)
points(x6G08_MD_765_proj[1],x6G08_MD_765_proj[2],pch=20)
text(x6G08_MD_765_proj[1],x6G08_MD_765_proj[2],pos=4,label="765",col="red")

x6G08_MD_766_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_766.txt")
x6G08_MD_766_proj<-project.pca(x6G08_MD_766_raw,pcadata)
points(x6G08_MD_766_proj[1],x6G08_MD_766_proj[2],pch=20)
text(x6G08_MD_766_proj[1],x6G08_MD_766_proj[2],pos=4,label="766",col="red")

x6G08_MD_767_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_767.txt")
x6G08_MD_767_proj<-project.pca(x6G08_MD_767_raw,pcadata)
points(x6G08_MD_767_proj[1],x6G08_MD_767_proj[2],pch=20)
text(x6G08_MD_767_proj[1],x6G08_MD_767_proj[2],pos=4,label="767",col="red")

x6G08_MD_768_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_768.txt")
x6G08_MD_768_proj<-project.pca(x6G08_MD_768_raw,pcadata)
points(x6G08_MD_768_proj[1],x6G08_MD_768_proj[2],pch=20)
text(x6G08_MD_768_proj[1],x6G08_MD_768_proj[2],pos=4,label="768",col="red")

x6G08_MD_769_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_769.txt")
x6G08_MD_769_proj<-project.pca(x6G08_MD_769_raw,pcadata)
points(x6G08_MD_769_proj[1],x6G08_MD_769_proj[2],pch=20)
text(x6G08_MD_769_proj[1],x6G08_MD_769_proj[2],pos=4,label="769",col="red")

x6G08_MD_770_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_770.txt")
x6G08_MD_770_proj<-project.pca(x6G08_MD_770_raw,pcadata)
points(x6G08_MD_770_proj[1],x6G08_MD_770_proj[2],pch=20)
text(x6G08_MD_770_proj[1],x6G08_MD_770_proj[2],pos=4,label="770",col="red")

x6G08_MD_771_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_771.txt")
x6G08_MD_771_proj<-project.pca(x6G08_MD_771_raw,pcadata)
points(x6G08_MD_771_proj[1],x6G08_MD_771_proj[2],pch=20)
text(x6G08_MD_771_proj[1],x6G08_MD_771_proj[2],pos=4,label="771",col="red")

x6G08_MD_772_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_772.txt")
x6G08_MD_772_proj<-project.pca(x6G08_MD_772_raw,pcadata)
points(x6G08_MD_772_proj[1],x6G08_MD_772_proj[2],pch=20)
text(x6G08_MD_772_proj[1],x6G08_MD_772_proj[2],pos=4,label="772",col="red")

x6G08_MD_773_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_773.txt")
x6G08_MD_773_proj<-project.pca(x6G08_MD_773_raw,pcadata)
points(x6G08_MD_773_proj[1],x6G08_MD_773_proj[2],pch=20)
text(x6G08_MD_773_proj[1],x6G08_MD_773_proj[2],pos=4,label="773",col="red")

x6G08_MD_774_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_774.txt")
x6G08_MD_774_proj<-project.pca(x6G08_MD_774_raw,pcadata)
points(x6G08_MD_774_proj[1],x6G08_MD_774_proj[2],pch=20)
text(x6G08_MD_774_proj[1],x6G08_MD_774_proj[2],pos=4,label="774",col="red")

x6G08_MD_775_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_775.txt")
x6G08_MD_775_proj<-project.pca(x6G08_MD_775_raw,pcadata)
points(x6G08_MD_775_proj[1],x6G08_MD_775_proj[2],pch=20)
text(x6G08_MD_775_proj[1],x6G08_MD_775_proj[2],pos=4,label="775",col="red")

x6G08_MD_776_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_776.txt")
x6G08_MD_776_proj<-project.pca(x6G08_MD_776_raw,pcadata)
points(x6G08_MD_776_proj[1],x6G08_MD_776_proj[2],pch=20)
text(x6G08_MD_776_proj[1],x6G08_MD_776_proj[2],pos=4,label="776",col="red")

x6G08_MD_777_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_777.txt")
x6G08_MD_777_proj<-project.pca(x6G08_MD_777_raw,pcadata)
points(x6G08_MD_777_proj[1],x6G08_MD_777_proj[2],pch=20)
text(x6G08_MD_777_proj[1],x6G08_MD_777_proj[2],pos=4,label="777",col="green")

x6G08_MD_778_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_778.txt")
x6G08_MD_778_proj<-project.pca(x6G08_MD_778_raw,pcadata)
points(x6G08_MD_778_proj[1],x6G08_MD_778_proj[2],pch=20)
text(x6G08_MD_778_proj[1],x6G08_MD_778_proj[2],pos=4,label="778",col="green")

x6G08_MD_779_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_779.txt")
x6G08_MD_779_proj<-project.pca(x6G08_MD_779_raw,pcadata)
points(x6G08_MD_779_proj[1],x6G08_MD_779_proj[2],pch=20)
text(x6G08_MD_779_proj[1],x6G08_MD_779_proj[2],pos=4,label="779",col="green")

x6G08_MD_780_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_780.txt")
x6G08_MD_780_proj<-project.pca(x6G08_MD_780_raw,pcadata)
points(x6G08_MD_780_proj[1],x6G08_MD_780_proj[2],pch=20)
text(x6G08_MD_780_proj[1],x6G08_MD_780_proj[2],pos=4,label="780",col="red")

x6G08_MD_781_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_781.txt")
x6G08_MD_781_proj<-project.pca(x6G08_MD_781_raw,pcadata)
points(x6G08_MD_781_proj[1],x6G08_MD_781_proj[2],pch=20)
text(x6G08_MD_781_proj[1],x6G08_MD_781_proj[2],pos=4,label="781",col="red")

x6G08_MD_782_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_782.txt")
x6G08_MD_782_proj<-project.pca(x6G08_MD_782_raw,pcadata)
points(x6G08_MD_782_proj[1],x6G08_MD_782_proj[2],pch=20)
text(x6G08_MD_782_proj[1],x6G08_MD_782_proj[2],pos=4,label="782",col="red")

x6G08_MD_783_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_783.txt")
x6G08_MD_783_proj<-project.pca(x6G08_MD_783_raw,pcadata)
points(x6G08_MD_783_proj[1],x6G08_MD_783_proj[2],pch=20)
text(x6G08_MD_783_proj[1],x6G08_MD_783_proj[2],pos=4,label="783",col="green")

x6G08_MD_784_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_784.txt")
x6G08_MD_784_proj<-project.pca(x6G08_MD_784_raw,pcadata)
points(x6G08_MD_784_proj[1],x6G08_MD_784_proj[2],pch=20)
text(x6G08_MD_784_proj[1],x6G08_MD_784_proj[2],pos=4,label="784",col="red")

x6G08_MD_785_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_785.txt")
x6G08_MD_785_proj<-project.pca(x6G08_MD_785_raw,pcadata)
points(x6G08_MD_785_proj[1],x6G08_MD_785_proj[2],pch=20)
text(x6G08_MD_785_proj[1],x6G08_MD_785_proj[2],pos=4,label="785",col="red")

x6G08_MD_786_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_786.txt")
x6G08_MD_786_proj<-project.pca(x6G08_MD_786_raw,pcadata)
points(x6G08_MD_786_proj[1],x6G08_MD_786_proj[2],pch=20)
text(x6G08_MD_786_proj[1],x6G08_MD_786_proj[2],pos=4,label="786",col="red")

x6G08_MD_787_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_787.txt")
x6G08_MD_787_proj<-project.pca(x6G08_MD_787_raw,pcadata)
points(x6G08_MD_787_proj[1],x6G08_MD_787_proj[2],pch=20)
text(x6G08_MD_787_proj[1],x6G08_MD_787_proj[2],pos=4,label="787",col="red")

x6G08_MD_788_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_788.txt")
x6G08_MD_788_proj<-project.pca(x6G08_MD_788_raw,pcadata)
points(x6G08_MD_788_proj[1],x6G08_MD_788_proj[2],pch=20)
text(x6G08_MD_788_proj[1],x6G08_MD_788_proj[2],pos=4,label="788",col="red")

x6G08_MD_789_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_789.txt")
x6G08_MD_789_proj<-project.pca(x6G08_MD_789_raw,pcadata)
points(x6G08_MD_789_proj[1],x6G08_MD_789_proj[2],pch=20)
text(x6G08_MD_789_proj[1],x6G08_MD_789_proj[2],pos=4,label="789",col="red")

x6G08_MD_790_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_790.txt")
x6G08_MD_790_proj<-project.pca(x6G08_MD_790_raw,pcadata)
points(x6G08_MD_790_proj[1],x6G08_MD_790_proj[2],pch=20)
text(x6G08_MD_790_proj[1],x6G08_MD_790_proj[2],pos=4,label="790",col="red")

x6G08_MD_791_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_791.txt")
x6G08_MD_791_proj<-project.pca(x6G08_MD_791_raw,pcadata)
points(x6G08_MD_791_proj[1],x6G08_MD_791_proj[2],pch=20)
text(x6G08_MD_791_proj[1],x6G08_MD_791_proj[2],pos=4,label="791",col="red")

x6G08_MD_792_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_792.txt")
x6G08_MD_792_proj<-project.pca(x6G08_MD_792_raw,pcadata)
points(x6G08_MD_792_proj[1],x6G08_MD_792_proj[2],pch=20)
text(x6G08_MD_792_proj[1],x6G08_MD_792_proj[2],pos=4,label="792",col="red")

x6G08_MD_793_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_793.txt")
x6G08_MD_793_proj<-project.pca(x6G08_MD_793_raw,pcadata)
points(x6G08_MD_793_proj[1],x6G08_MD_793_proj[2],pch=20)
text(x6G08_MD_793_proj[1],x6G08_MD_793_proj[2],pos=4,label="793",col="red")

x6G08_MD_794_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_794.txt")
x6G08_MD_794_proj<-project.pca(x6G08_MD_794_raw,pcadata)
points(x6G08_MD_794_proj[1],x6G08_MD_794_proj[2],pch=20)
text(x6G08_MD_794_proj[1],x6G08_MD_794_proj[2],pos=4,label="794",col="red")

x6G08_MD_795_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_795.txt")
x6G08_MD_795_proj<-project.pca(x6G08_MD_795_raw,pcadata)
points(x6G08_MD_795_proj[1],x6G08_MD_795_proj[2],pch=20)
text(x6G08_MD_795_proj[1],x6G08_MD_795_proj[2],pos=4,label="795",col="red")

x6G08_MD_796_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_796.txt")
x6G08_MD_796_proj<-project.pca(x6G08_MD_796_raw,pcadata)
points(x6G08_MD_796_proj[1],x6G08_MD_796_proj[2],pch=20)
text(x6G08_MD_796_proj[1],x6G08_MD_796_proj[2],pos=4,label="796",col="red")

x6G08_MD_797_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_797.txt")
x6G08_MD_797_proj<-project.pca(x6G08_MD_797_raw,pcadata)
points(x6G08_MD_797_proj[1],x6G08_MD_797_proj[2],pch=20)
text(x6G08_MD_797_proj[1],x6G08_MD_797_proj[2],pos=4,label="797",col="red")

x6G08_MD_798_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_798.txt")
x6G08_MD_798_proj<-project.pca(x6G08_MD_798_raw,pcadata)
points(x6G08_MD_798_proj[1],x6G08_MD_798_proj[2],pch=20)
text(x6G08_MD_798_proj[1],x6G08_MD_798_proj[2],pos=4,label="798",col="red")

x6G08_MD_799_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_799.txt")
x6G08_MD_799_proj<-project.pca(x6G08_MD_799_raw,pcadata)
points(x6G08_MD_799_proj[1],x6G08_MD_799_proj[2],pch=20)
text(x6G08_MD_799_proj[1],x6G08_MD_799_proj[2],pos=4,label="799",col="red")

x6G08_MD_800_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_800.txt")
x6G08_MD_800_proj<-project.pca(x6G08_MD_800_raw,pcadata)
points(x6G08_MD_800_proj[1],x6G08_MD_800_proj[2],pch=20)
text(x6G08_MD_800_proj[1],x6G08_MD_800_proj[2],pos=4,label="800",col="red")

x6G08_MD_801_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_801.txt")
x6G08_MD_801_proj<-project.pca(x6G08_MD_801_raw,pcadata)
points(x6G08_MD_801_proj[1],x6G08_MD_801_proj[2],pch=20)
text(x6G08_MD_801_proj[1],x6G08_MD_801_proj[2],pos=4,label="801",col="red")

x6G08_MD_802_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_802.txt")
x6G08_MD_802_proj<-project.pca(x6G08_MD_802_raw,pcadata)
points(x6G08_MD_802_proj[1],x6G08_MD_802_proj[2],pch=20)
text(x6G08_MD_802_proj[1],x6G08_MD_802_proj[2],pos=4,label="802",col="red")

x6G08_MD_803_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_803.txt")
x6G08_MD_803_proj<-project.pca(x6G08_MD_803_raw,pcadata)
points(x6G08_MD_803_proj[1],x6G08_MD_803_proj[2],pch=20)
text(x6G08_MD_803_proj[1],x6G08_MD_803_proj[2],pos=4,label="803",col="red")

x6G08_MD_804_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_804.txt")
x6G08_MD_804_proj<-project.pca(x6G08_MD_804_raw,pcadata)
points(x6G08_MD_804_proj[1],x6G08_MD_804_proj[2],pch=20)
text(x6G08_MD_804_proj[1],x6G08_MD_804_proj[2],pos=4,label="804",col="red")

x6G08_MD_805_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_805.txt")
x6G08_MD_805_proj<-project.pca(x6G08_MD_805_raw,pcadata)
points(x6G08_MD_805_proj[1],x6G08_MD_805_proj[2],pch=20)
text(x6G08_MD_805_proj[1],x6G08_MD_805_proj[2],pos=4,label="805",col="red")

x6G08_MD_806_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_806.txt")
x6G08_MD_806_proj<-project.pca(x6G08_MD_806_raw,pcadata)
points(x6G08_MD_806_proj[1],x6G08_MD_806_proj[2],pch=20)
text(x6G08_MD_806_proj[1],x6G08_MD_806_proj[2],pos=4,label="806",col="red")

x6G08_MD_807_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_807.txt")
x6G08_MD_807_proj<-project.pca(x6G08_MD_807_raw,pcadata)
points(x6G08_MD_807_proj[1],x6G08_MD_807_proj[2],pch=20)
text(x6G08_MD_807_proj[1],x6G08_MD_807_proj[2],pos=4,label="807",col="red")

x6G08_MD_808_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_808.txt")
x6G08_MD_808_proj<-project.pca(x6G08_MD_808_raw,pcadata)
points(x6G08_MD_808_proj[1],x6G08_MD_808_proj[2],pch=20)
text(x6G08_MD_808_proj[1],x6G08_MD_808_proj[2],pos=4,label="808",col="red")

x6G08_MD_809_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_809.txt")
x6G08_MD_809_proj<-project.pca(x6G08_MD_809_raw,pcadata)
points(x6G08_MD_809_proj[1],x6G08_MD_809_proj[2],pch=20)
text(x6G08_MD_809_proj[1],x6G08_MD_809_proj[2],pos=4,label="809",col="red")

x6G08_MD_810_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_810.txt")
x6G08_MD_810_proj<-project.pca(x6G08_MD_810_raw,pcadata)
points(x6G08_MD_810_proj[1],x6G08_MD_810_proj[2],pch=20)
text(x6G08_MD_810_proj[1],x6G08_MD_810_proj[2],pos=4,label="810",col="red")

x6G08_MD_811_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_811.txt")
x6G08_MD_811_proj<-project.pca(x6G08_MD_811_raw,pcadata)
points(x6G08_MD_811_proj[1],x6G08_MD_811_proj[2],pch=20)
text(x6G08_MD_811_proj[1],x6G08_MD_811_proj[2],pos=4,label="811",col="red")

x6G08_MD_812_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_812.txt")
x6G08_MD_812_proj<-project.pca(x6G08_MD_812_raw,pcadata)
points(x6G08_MD_812_proj[1],x6G08_MD_812_proj[2],pch=20)
text(x6G08_MD_812_proj[1],x6G08_MD_812_proj[2],pos=4,label="812",col="red")

x6G08_MD_813_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_813.txt")
x6G08_MD_813_proj<-project.pca(x6G08_MD_813_raw,pcadata)
points(x6G08_MD_813_proj[1],x6G08_MD_813_proj[2],pch=20)
text(x6G08_MD_813_proj[1],x6G08_MD_813_proj[2],pos=4,label="813",col="red")

x6G08_MD_814_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_814.txt")
x6G08_MD_814_proj<-project.pca(x6G08_MD_814_raw,pcadata)
points(x6G08_MD_814_proj[1],x6G08_MD_814_proj[2],pch=20)
text(x6G08_MD_814_proj[1],x6G08_MD_814_proj[2],pos=4,label="814",col="red")

x6G08_MD_815_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_815.txt")
x6G08_MD_815_proj<-project.pca(x6G08_MD_815_raw,pcadata)
points(x6G08_MD_815_proj[1],x6G08_MD_815_proj[2],pch=20)
text(x6G08_MD_815_proj[1],x6G08_MD_815_proj[2],pos=4,label="815",col="red")

x6G08_MD_816_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_816.txt")
x6G08_MD_816_proj<-project.pca(x6G08_MD_816_raw,pcadata)
points(x6G08_MD_816_proj[1],x6G08_MD_816_proj[2],pch=20)
text(x6G08_MD_816_proj[1],x6G08_MD_816_proj[2],pos=4,label="816",col="blue")

x6G08_MD_817_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_817.txt")
x6G08_MD_817_proj<-project.pca(x6G08_MD_817_raw,pcadata)
points(x6G08_MD_817_proj[1],x6G08_MD_817_proj[2],pch=20)
text(x6G08_MD_817_proj[1],x6G08_MD_817_proj[2],pos=4,label="817",col="red")

x6G08_MD_818_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_818.txt")
x6G08_MD_818_proj<-project.pca(x6G08_MD_818_raw,pcadata)
points(x6G08_MD_818_proj[1],x6G08_MD_818_proj[2],pch=20)
text(x6G08_MD_818_proj[1],x6G08_MD_818_proj[2],pos=4,label="818",col="blue")

x6G08_MD_819_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_819.txt")
x6G08_MD_819_proj<-project.pca(x6G08_MD_819_raw,pcadata)
points(x6G08_MD_819_proj[1],x6G08_MD_819_proj[2],pch=20)
text(x6G08_MD_819_proj[1],x6G08_MD_819_proj[2],pos=4,label="819",col="blue")

x6G08_MD_820_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_820.txt")
x6G08_MD_820_proj<-project.pca(x6G08_MD_820_raw,pcadata)
points(x6G08_MD_820_proj[1],x6G08_MD_820_proj[2],pch=20)
text(x6G08_MD_820_proj[1],x6G08_MD_820_proj[2],pos=4,label="820",col="blue")

x6G08_MD_821_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_821.txt")
x6G08_MD_821_proj<-project.pca(x6G08_MD_821_raw,pcadata)
points(x6G08_MD_821_proj[1],x6G08_MD_821_proj[2],pch=20)
text(x6G08_MD_821_proj[1],x6G08_MD_821_proj[2],pos=4,label="821",col="red")

x6G08_MD_822_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_822.txt")
x6G08_MD_822_proj<-project.pca(x6G08_MD_822_raw,pcadata)
points(x6G08_MD_822_proj[1],x6G08_MD_822_proj[2],pch=20)
text(x6G08_MD_822_proj[1],x6G08_MD_822_proj[2],pos=4,label="822",col="red")

x6G08_MD_823_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_823.txt")
x6G08_MD_823_proj<-project.pca(x6G08_MD_823_raw,pcadata)
points(x6G08_MD_823_proj[1],x6G08_MD_823_proj[2],pch=20)
text(x6G08_MD_823_proj[1],x6G08_MD_823_proj[2],pos=4,label="823",col="red")

x6G08_MD_824_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_824.txt")
x6G08_MD_824_proj<-project.pca(x6G08_MD_824_raw,pcadata)
points(x6G08_MD_824_proj[1],x6G08_MD_824_proj[2],pch=20)
text(x6G08_MD_824_proj[1],x6G08_MD_824_proj[2],pos=4,label="824",col="red")

x6G08_MD_825_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_825.txt")
x6G08_MD_825_proj<-project.pca(x6G08_MD_825_raw,pcadata)
points(x6G08_MD_825_proj[1],x6G08_MD_825_proj[2],pch=20)
text(x6G08_MD_825_proj[1],x6G08_MD_825_proj[2],pos=4,label="825",col="red")

x6G08_MD_826_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_826.txt")
x6G08_MD_826_proj<-project.pca(x6G08_MD_826_raw,pcadata)
points(x6G08_MD_826_proj[1],x6G08_MD_826_proj[2],pch=20)
text(x6G08_MD_826_proj[1],x6G08_MD_826_proj[2],pos=4,label="826",col="red")

x6G08_MD_827_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_827.txt")
x6G08_MD_827_proj<-project.pca(x6G08_MD_827_raw,pcadata)
points(x6G08_MD_827_proj[1],x6G08_MD_827_proj[2],pch=20)
text(x6G08_MD_827_proj[1],x6G08_MD_827_proj[2],pos=4,label="827",col="red")

x6G08_MD_828_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_828.txt")
x6G08_MD_828_proj<-project.pca(x6G08_MD_828_raw,pcadata)
points(x6G08_MD_828_proj[1],x6G08_MD_828_proj[2],pch=20)
text(x6G08_MD_828_proj[1],x6G08_MD_828_proj[2],pos=4,label="828",col="red")

x6G08_MD_829_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_829.txt")
x6G08_MD_829_proj<-project.pca(x6G08_MD_829_raw,pcadata)
points(x6G08_MD_829_proj[1],x6G08_MD_829_proj[2],pch=20)
text(x6G08_MD_829_proj[1],x6G08_MD_829_proj[2],pos=4,label="829",col="red")

x6G08_MD_830_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_830.txt")
x6G08_MD_830_proj<-project.pca(x6G08_MD_830_raw,pcadata)
points(x6G08_MD_830_proj[1],x6G08_MD_830_proj[2],pch=20)
text(x6G08_MD_830_proj[1],x6G08_MD_830_proj[2],pos=4,label="830",col="red")

x6G08_MD_831_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_831.txt")
x6G08_MD_831_proj<-project.pca(x6G08_MD_831_raw,pcadata)
points(x6G08_MD_831_proj[1],x6G08_MD_831_proj[2],pch=20)
text(x6G08_MD_831_proj[1],x6G08_MD_831_proj[2],pos=4,label="831",col="red")

x6G08_MD_832_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_832.txt")
x6G08_MD_832_proj<-project.pca(x6G08_MD_832_raw,pcadata)
points(x6G08_MD_832_proj[1],x6G08_MD_832_proj[2],pch=20)
text(x6G08_MD_832_proj[1],x6G08_MD_832_proj[2],pos=4,label="832",col="red")

x6G08_MD_833_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_833.txt")
x6G08_MD_833_proj<-project.pca(x6G08_MD_833_raw,pcadata)
points(x6G08_MD_833_proj[1],x6G08_MD_833_proj[2],pch=20)
text(x6G08_MD_833_proj[1],x6G08_MD_833_proj[2],pos=4,label="833",col="red")

x6G08_MD_834_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_834.txt")
x6G08_MD_834_proj<-project.pca(x6G08_MD_834_raw,pcadata)
points(x6G08_MD_834_proj[1],x6G08_MD_834_proj[2],pch=20)
text(x6G08_MD_834_proj[1],x6G08_MD_834_proj[2],pos=4,label="834",col="red")

x6G08_MD_835_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_835.txt")
x6G08_MD_835_proj<-project.pca(x6G08_MD_835_raw,pcadata)
points(x6G08_MD_835_proj[1],x6G08_MD_835_proj[2],pch=20)
text(x6G08_MD_835_proj[1],x6G08_MD_835_proj[2],pos=4,label="835",col="red")

x6G08_MD_836_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_836.txt")
x6G08_MD_836_proj<-project.pca(x6G08_MD_836_raw,pcadata)
points(x6G08_MD_836_proj[1],x6G08_MD_836_proj[2],pch=20)
text(x6G08_MD_836_proj[1],x6G08_MD_836_proj[2],pos=4,label="836",col="red")

x6G08_MD_837_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_837.txt")
x6G08_MD_837_proj<-project.pca(x6G08_MD_837_raw,pcadata)
points(x6G08_MD_837_proj[1],x6G08_MD_837_proj[2],pch=20)
text(x6G08_MD_837_proj[1],x6G08_MD_837_proj[2],pos=4,label="837",col="red")

x6G08_MD_838_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_838.txt")
x6G08_MD_838_proj<-project.pca(x6G08_MD_838_raw,pcadata)
points(x6G08_MD_838_proj[1],x6G08_MD_838_proj[2],pch=20)
text(x6G08_MD_838_proj[1],x6G08_MD_838_proj[2],pos=4,label="838",col="red")

x6G08_MD_839_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_839.txt")
x6G08_MD_839_proj<-project.pca(x6G08_MD_839_raw,pcadata)
points(x6G08_MD_839_proj[1],x6G08_MD_839_proj[2],pch=20)
text(x6G08_MD_839_proj[1],x6G08_MD_839_proj[2],pos=4,label="839",col="red")

x6G08_MD_840_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_840.txt")
x6G08_MD_840_proj<-project.pca(x6G08_MD_840_raw,pcadata)
points(x6G08_MD_840_proj[1],x6G08_MD_840_proj[2],pch=20)
text(x6G08_MD_840_proj[1],x6G08_MD_840_proj[2],pos=4,label="840",col="red")

x6G08_MD_841_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_841.txt")
x6G08_MD_841_proj<-project.pca(x6G08_MD_841_raw,pcadata)
points(x6G08_MD_841_proj[1],x6G08_MD_841_proj[2],pch=20)
text(x6G08_MD_841_proj[1],x6G08_MD_841_proj[2],pos=4,label="841",col="red")

x6G08_MD_842_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_842.txt")
x6G08_MD_842_proj<-project.pca(x6G08_MD_842_raw,pcadata)
points(x6G08_MD_842_proj[1],x6G08_MD_842_proj[2],pch=20)
text(x6G08_MD_842_proj[1],x6G08_MD_842_proj[2],pos=4,label="842",col="green")

x6G08_MD_843_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_843.txt")
x6G08_MD_843_proj<-project.pca(x6G08_MD_843_raw,pcadata)
points(x6G08_MD_843_proj[1],x6G08_MD_843_proj[2],pch=20)
text(x6G08_MD_843_proj[1],x6G08_MD_843_proj[2],pos=4,label="843",col="red")

x6G08_MD_844_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_844.txt")
x6G08_MD_844_proj<-project.pca(x6G08_MD_844_raw,pcadata)
points(x6G08_MD_844_proj[1],x6G08_MD_844_proj[2],pch=20)
text(x6G08_MD_844_proj[1],x6G08_MD_844_proj[2],pos=4,label="844",col="red")

x6G08_MD_845_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_845.txt")
x6G08_MD_845_proj<-project.pca(x6G08_MD_845_raw,pcadata)
points(x6G08_MD_845_proj[1],x6G08_MD_845_proj[2],pch=20)
text(x6G08_MD_845_proj[1],x6G08_MD_845_proj[2],pos=4,label="845",col="red")

x6G08_MD_846_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_846.txt")
x6G08_MD_846_proj<-project.pca(x6G08_MD_846_raw,pcadata)
points(x6G08_MD_846_proj[1],x6G08_MD_846_proj[2],pch=20)
text(x6G08_MD_846_proj[1],x6G08_MD_846_proj[2],pos=4,label="846",col="red")

x6G08_MD_847_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_847.txt")
x6G08_MD_847_proj<-project.pca(x6G08_MD_847_raw,pcadata)
points(x6G08_MD_847_proj[1],x6G08_MD_847_proj[2],pch=20)
text(x6G08_MD_847_proj[1],x6G08_MD_847_proj[2],pos=4,label="847",col="red")

x6G08_MD_848_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_848.txt")
x6G08_MD_848_proj<-project.pca(x6G08_MD_848_raw,pcadata)
points(x6G08_MD_848_proj[1],x6G08_MD_848_proj[2],pch=20)
text(x6G08_MD_848_proj[1],x6G08_MD_848_proj[2],pos=4,label="848",col="red")

x6G08_MD_849_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_849.txt")
x6G08_MD_849_proj<-project.pca(x6G08_MD_849_raw,pcadata)
points(x6G08_MD_849_proj[1],x6G08_MD_849_proj[2],pch=20)
text(x6G08_MD_849_proj[1],x6G08_MD_849_proj[2],pos=4,label="849",col="red")

x6G08_MD_850_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_850.txt")
x6G08_MD_850_proj<-project.pca(x6G08_MD_850_raw,pcadata)
points(x6G08_MD_850_proj[1],x6G08_MD_850_proj[2],pch=20)
text(x6G08_MD_850_proj[1],x6G08_MD_850_proj[2],pos=4,label="850",col="red")

x6G08_MD_851_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_851.txt")
x6G08_MD_851_proj<-project.pca(x6G08_MD_851_raw,pcadata)
points(x6G08_MD_851_proj[1],x6G08_MD_851_proj[2],pch=20)
text(x6G08_MD_851_proj[1],x6G08_MD_851_proj[2],pos=4,label="851",col="red")

x6G08_MD_852_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_852.txt")
x6G08_MD_852_proj<-project.pca(x6G08_MD_852_raw,pcadata)
points(x6G08_MD_852_proj[1],x6G08_MD_852_proj[2],pch=20)
text(x6G08_MD_852_proj[1],x6G08_MD_852_proj[2],pos=4,label="852",col="red")

x6G08_MD_853_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_853.txt")
x6G08_MD_853_proj<-project.pca(x6G08_MD_853_raw,pcadata)
points(x6G08_MD_853_proj[1],x6G08_MD_853_proj[2],pch=20)
text(x6G08_MD_853_proj[1],x6G08_MD_853_proj[2],pos=4,label="853",col="red")

x6G08_MD_854_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_854.txt")
x6G08_MD_854_proj<-project.pca(x6G08_MD_854_raw,pcadata)
points(x6G08_MD_854_proj[1],x6G08_MD_854_proj[2],pch=20)
text(x6G08_MD_854_proj[1],x6G08_MD_854_proj[2],pos=4,label="854",col="red")

x6G08_MD_855_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_855.txt")
x6G08_MD_855_proj<-project.pca(x6G08_MD_855_raw,pcadata)
points(x6G08_MD_855_proj[1],x6G08_MD_855_proj[2],pch=20)
text(x6G08_MD_855_proj[1],x6G08_MD_855_proj[2],pos=4,label="855",col="green")

x6G08_MD_856_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_856.txt")
x6G08_MD_856_proj<-project.pca(x6G08_MD_856_raw,pcadata)
points(x6G08_MD_856_proj[1],x6G08_MD_856_proj[2],pch=20)
text(x6G08_MD_856_proj[1],x6G08_MD_856_proj[2],pos=4,label="856",col="red")

x6G08_MD_857_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_857.txt")
x6G08_MD_857_proj<-project.pca(x6G08_MD_857_raw,pcadata)
points(x6G08_MD_857_proj[1],x6G08_MD_857_proj[2],pch=20)
text(x6G08_MD_857_proj[1],x6G08_MD_857_proj[2],pos=4,label="857",col="red")

x6G08_MD_858_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_858.txt")
x6G08_MD_858_proj<-project.pca(x6G08_MD_858_raw,pcadata)
points(x6G08_MD_858_proj[1],x6G08_MD_858_proj[2],pch=20)
text(x6G08_MD_858_proj[1],x6G08_MD_858_proj[2],pos=4,label="858",col="red")

x6G08_MD_859_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_859.txt")
x6G08_MD_859_proj<-project.pca(x6G08_MD_859_raw,pcadata)
points(x6G08_MD_859_proj[1],x6G08_MD_859_proj[2],pch=20)
text(x6G08_MD_859_proj[1],x6G08_MD_859_proj[2],pos=4,label="859",col="red")

x6G08_MD_860_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_860.txt")
x6G08_MD_860_proj<-project.pca(x6G08_MD_860_raw,pcadata)
points(x6G08_MD_860_proj[1],x6G08_MD_860_proj[2],pch=20)
text(x6G08_MD_860_proj[1],x6G08_MD_860_proj[2],pos=4,label="860",col="red")

x6G08_MD_861_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_861.txt")
x6G08_MD_861_proj<-project.pca(x6G08_MD_861_raw,pcadata)
points(x6G08_MD_861_proj[1],x6G08_MD_861_proj[2],pch=20)
text(x6G08_MD_861_proj[1],x6G08_MD_861_proj[2],pos=4,label="861",col="red")

x6G08_MD_862_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_862.txt")
x6G08_MD_862_proj<-project.pca(x6G08_MD_862_raw,pcadata)
points(x6G08_MD_862_proj[1],x6G08_MD_862_proj[2],pch=20)
text(x6G08_MD_862_proj[1],x6G08_MD_862_proj[2],pos=4,label="862",col="red")

x6G08_MD_863_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_863.txt")
x6G08_MD_863_proj<-project.pca(x6G08_MD_863_raw,pcadata)
points(x6G08_MD_863_proj[1],x6G08_MD_863_proj[2],pch=20)
text(x6G08_MD_863_proj[1],x6G08_MD_863_proj[2],pos=4,label="863",col="red")

x6G08_MD_864_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_864.txt")
x6G08_MD_864_proj<-project.pca(x6G08_MD_864_raw,pcadata)
points(x6G08_MD_864_proj[1],x6G08_MD_864_proj[2],pch=20)
text(x6G08_MD_864_proj[1],x6G08_MD_864_proj[2],pos=4,label="864",col="red")

x6G08_MD_865_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_865.txt")
x6G08_MD_865_proj<-project.pca(x6G08_MD_865_raw,pcadata)
points(x6G08_MD_865_proj[1],x6G08_MD_865_proj[2],pch=20)
text(x6G08_MD_865_proj[1],x6G08_MD_865_proj[2],pos=4,label="865",col="red")

x6G08_MD_866_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_866.txt")
x6G08_MD_866_proj<-project.pca(x6G08_MD_866_raw,pcadata)
points(x6G08_MD_866_proj[1],x6G08_MD_866_proj[2],pch=20)
text(x6G08_MD_866_proj[1],x6G08_MD_866_proj[2],pos=4,label="866",col="red")

x6G08_MD_867_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_867.txt")
x6G08_MD_867_proj<-project.pca(x6G08_MD_867_raw,pcadata)
points(x6G08_MD_867_proj[1],x6G08_MD_867_proj[2],pch=20)
text(x6G08_MD_867_proj[1],x6G08_MD_867_proj[2],pos=4,label="867",col="red")

x6G08_MD_868_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_868.txt")
x6G08_MD_868_proj<-project.pca(x6G08_MD_868_raw,pcadata)
points(x6G08_MD_868_proj[1],x6G08_MD_868_proj[2],pch=20)
text(x6G08_MD_868_proj[1],x6G08_MD_868_proj[2],pos=4,label="868",col="red")

x6G08_MD_869_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_869.txt")
x6G08_MD_869_proj<-project.pca(x6G08_MD_869_raw,pcadata)
points(x6G08_MD_869_proj[1],x6G08_MD_869_proj[2],pch=20)
text(x6G08_MD_869_proj[1],x6G08_MD_869_proj[2],pos=4,label="869",col="red")

x6G08_MD_870_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_870.txt")
x6G08_MD_870_proj<-project.pca(x6G08_MD_870_raw,pcadata)
points(x6G08_MD_870_proj[1],x6G08_MD_870_proj[2],pch=20)
text(x6G08_MD_870_proj[1],x6G08_MD_870_proj[2],pos=4,label="870",col="red")

x6G08_MD_871_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_871.txt")
x6G08_MD_871_proj<-project.pca(x6G08_MD_871_raw,pcadata)
points(x6G08_MD_871_proj[1],x6G08_MD_871_proj[2],pch=20)
text(x6G08_MD_871_proj[1],x6G08_MD_871_proj[2],pos=4,label="871",col="red")

x6G08_MD_872_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_872.txt")
x6G08_MD_872_proj<-project.pca(x6G08_MD_872_raw,pcadata)
points(x6G08_MD_872_proj[1],x6G08_MD_872_proj[2],pch=20)
text(x6G08_MD_872_proj[1],x6G08_MD_872_proj[2],pos=4,label="872",col="red")

x6G08_MD_873_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_873.txt")
x6G08_MD_873_proj<-project.pca(x6G08_MD_873_raw,pcadata)
points(x6G08_MD_873_proj[1],x6G08_MD_873_proj[2],pch=20)
text(x6G08_MD_873_proj[1],x6G08_MD_873_proj[2],pos=4,label="873",col="red")

x6G08_MD_874_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_874.txt")
x6G08_MD_874_proj<-project.pca(x6G08_MD_874_raw,pcadata)
points(x6G08_MD_874_proj[1],x6G08_MD_874_proj[2],pch=20)
text(x6G08_MD_874_proj[1],x6G08_MD_874_proj[2],pos=4,label="874",col="red")

x6G08_MD_875_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_875.txt")
x6G08_MD_875_proj<-project.pca(x6G08_MD_875_raw,pcadata)
points(x6G08_MD_875_proj[1],x6G08_MD_875_proj[2],pch=20)
text(x6G08_MD_875_proj[1],x6G08_MD_875_proj[2],pos=4,label="875",col="red")

x6G08_MD_876_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_876.txt")
x6G08_MD_876_proj<-project.pca(x6G08_MD_876_raw,pcadata)
points(x6G08_MD_876_proj[1],x6G08_MD_876_proj[2],pch=20)
text(x6G08_MD_876_proj[1],x6G08_MD_876_proj[2],pos=4,label="876",col="red")

x6G08_MD_877_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_877.txt")
x6G08_MD_877_proj<-project.pca(x6G08_MD_877_raw,pcadata)
points(x6G08_MD_877_proj[1],x6G08_MD_877_proj[2],pch=20)
text(x6G08_MD_877_proj[1],x6G08_MD_877_proj[2],pos=4,label="877",col="red")

x6G08_MD_878_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_878.txt")
x6G08_MD_878_proj<-project.pca(x6G08_MD_878_raw,pcadata)
points(x6G08_MD_878_proj[1],x6G08_MD_878_proj[2],pch=20)
text(x6G08_MD_878_proj[1],x6G08_MD_878_proj[2],pos=4,label="878",col="red")

x6G08_MD_879_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_879.txt")
x6G08_MD_879_proj<-project.pca(x6G08_MD_879_raw,pcadata)
points(x6G08_MD_879_proj[1],x6G08_MD_879_proj[2],pch=20)
text(x6G08_MD_879_proj[1],x6G08_MD_879_proj[2],pos=4,label="879",col="red")

x6G08_MD_880_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_880.txt")
x6G08_MD_880_proj<-project.pca(x6G08_MD_880_raw,pcadata)
points(x6G08_MD_880_proj[1],x6G08_MD_880_proj[2],pch=20)
text(x6G08_MD_880_proj[1],x6G08_MD_880_proj[2],pos=4,label="880",col="red")

x6G08_MD_881_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_881.txt")
x6G08_MD_881_proj<-project.pca(x6G08_MD_881_raw,pcadata)
points(x6G08_MD_881_proj[1],x6G08_MD_881_proj[2],pch=20)
text(x6G08_MD_881_proj[1],x6G08_MD_881_proj[2],pos=4,label="881",col="red")

x6G08_MD_882_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_882.txt")
x6G08_MD_882_proj<-project.pca(x6G08_MD_882_raw,pcadata)
points(x6G08_MD_882_proj[1],x6G08_MD_882_proj[2],pch=20)
text(x6G08_MD_882_proj[1],x6G08_MD_882_proj[2],pos=4,label="882",col="red")

x6G08_MD_883_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_883.txt")
x6G08_MD_883_proj<-project.pca(x6G08_MD_883_raw,pcadata)
points(x6G08_MD_883_proj[1],x6G08_MD_883_proj[2],pch=20)
text(x6G08_MD_883_proj[1],x6G08_MD_883_proj[2],pos=4,label="883",col="red")

x6G08_MD_884_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_884.txt")
x6G08_MD_884_proj<-project.pca(x6G08_MD_884_raw,pcadata)
points(x6G08_MD_884_proj[1],x6G08_MD_884_proj[2],pch=20)
text(x6G08_MD_884_proj[1],x6G08_MD_884_proj[2],pos=4,label="884",col="red")

x6G08_MD_885_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_885.txt")
x6G08_MD_885_proj<-project.pca(x6G08_MD_885_raw,pcadata)
points(x6G08_MD_885_proj[1],x6G08_MD_885_proj[2],pch=20)
text(x6G08_MD_885_proj[1],x6G08_MD_885_proj[2],pos=4,label="885",col="red")

x6G08_MD_886_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_886.txt")
x6G08_MD_886_proj<-project.pca(x6G08_MD_886_raw,pcadata)
points(x6G08_MD_886_proj[1],x6G08_MD_886_proj[2],pch=20)
text(x6G08_MD_886_proj[1],x6G08_MD_886_proj[2],pos=4,label="886",col="red")

x6G08_MD_887_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_887.txt")
x6G08_MD_887_proj<-project.pca(x6G08_MD_887_raw,pcadata)
points(x6G08_MD_887_proj[1],x6G08_MD_887_proj[2],pch=20)
text(x6G08_MD_887_proj[1],x6G08_MD_887_proj[2],pos=4,label="887",col="red")

x6G08_MD_888_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_888.txt")
x6G08_MD_888_proj<-project.pca(x6G08_MD_888_raw,pcadata)
points(x6G08_MD_888_proj[1],x6G08_MD_888_proj[2],pch=20)
text(x6G08_MD_888_proj[1],x6G08_MD_888_proj[2],pos=4,label="888",col="red")

x6G08_MD_889_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_889.txt")
x6G08_MD_889_proj<-project.pca(x6G08_MD_889_raw,pcadata)
points(x6G08_MD_889_proj[1],x6G08_MD_889_proj[2],pch=20)
text(x6G08_MD_889_proj[1],x6G08_MD_889_proj[2],pos=4,label="889",col="red")

x6G08_MD_890_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_890.txt")
x6G08_MD_890_proj<-project.pca(x6G08_MD_890_raw,pcadata)
points(x6G08_MD_890_proj[1],x6G08_MD_890_proj[2],pch=20)
text(x6G08_MD_890_proj[1],x6G08_MD_890_proj[2],pos=4,label="890",col="red")

x6G08_MD_891_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_891.txt")
x6G08_MD_891_proj<-project.pca(x6G08_MD_891_raw,pcadata)
points(x6G08_MD_891_proj[1],x6G08_MD_891_proj[2],pch=20)
text(x6G08_MD_891_proj[1],x6G08_MD_891_proj[2],pos=4,label="891",col="red")

x6G08_MD_892_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_892.txt")
x6G08_MD_892_proj<-project.pca(x6G08_MD_892_raw,pcadata)
points(x6G08_MD_892_proj[1],x6G08_MD_892_proj[2],pch=20)
text(x6G08_MD_892_proj[1],x6G08_MD_892_proj[2],pos=4,label="892",col="red")

x6G08_MD_893_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_893.txt")
x6G08_MD_893_proj<-project.pca(x6G08_MD_893_raw,pcadata)
points(x6G08_MD_893_proj[1],x6G08_MD_893_proj[2],pch=20)
text(x6G08_MD_893_proj[1],x6G08_MD_893_proj[2],pos=4,label="893",col="red")

x6G08_MD_894_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_894.txt")
x6G08_MD_894_proj<-project.pca(x6G08_MD_894_raw,pcadata)
points(x6G08_MD_894_proj[1],x6G08_MD_894_proj[2],pch=20)
text(x6G08_MD_894_proj[1],x6G08_MD_894_proj[2],pos=4,label="894",col="red")

x6G08_MD_895_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_895.txt")
x6G08_MD_895_proj<-project.pca(x6G08_MD_895_raw,pcadata)
points(x6G08_MD_895_proj[1],x6G08_MD_895_proj[2],pch=20)
text(x6G08_MD_895_proj[1],x6G08_MD_895_proj[2],pos=4,label="895",col="green")

x6G08_MD_896_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_896.txt")
x6G08_MD_896_proj<-project.pca(x6G08_MD_896_raw,pcadata)
points(x6G08_MD_896_proj[1],x6G08_MD_896_proj[2],pch=20)
text(x6G08_MD_896_proj[1],x6G08_MD_896_proj[2],pos=4,label="896",col="red")

x6G08_MD_897_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_897.txt")
x6G08_MD_897_proj<-project.pca(x6G08_MD_897_raw,pcadata)
points(x6G08_MD_897_proj[1],x6G08_MD_897_proj[2],pch=20)
text(x6G08_MD_897_proj[1],x6G08_MD_897_proj[2],pos=4,label="897",col="green")

x6G08_MD_898_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_898.txt")
x6G08_MD_898_proj<-project.pca(x6G08_MD_898_raw,pcadata)
points(x6G08_MD_898_proj[1],x6G08_MD_898_proj[2],pch=20)
text(x6G08_MD_898_proj[1],x6G08_MD_898_proj[2],pos=4,label="898",col="green")

x6G08_MD_899_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_899.txt")
x6G08_MD_899_proj<-project.pca(x6G08_MD_899_raw,pcadata)
points(x6G08_MD_899_proj[1],x6G08_MD_899_proj[2],pch=20)
text(x6G08_MD_899_proj[1],x6G08_MD_899_proj[2],pos=4,label="899",col="red")

x6G08_MD_900_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_900.txt")
x6G08_MD_900_proj<-project.pca(x6G08_MD_900_raw,pcadata)
points(x6G08_MD_900_proj[1],x6G08_MD_900_proj[2],pch=20)
text(x6G08_MD_900_proj[1],x6G08_MD_900_proj[2],pos=4,label="900",col="red")

x6G08_MD_901_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_901.txt")
x6G08_MD_901_proj<-project.pca(x6G08_MD_901_raw,pcadata)
points(x6G08_MD_901_proj[1],x6G08_MD_901_proj[2],pch=20)
text(x6G08_MD_901_proj[1],x6G08_MD_901_proj[2],pos=4,label="901",col="red")

x6G08_MD_902_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_902.txt")
x6G08_MD_902_proj<-project.pca(x6G08_MD_902_raw,pcadata)
points(x6G08_MD_902_proj[1],x6G08_MD_902_proj[2],pch=20)
text(x6G08_MD_902_proj[1],x6G08_MD_902_proj[2],pos=4,label="902",col="red")

x6G08_MD_903_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_903.txt")
x6G08_MD_903_proj<-project.pca(x6G08_MD_903_raw,pcadata)
points(x6G08_MD_903_proj[1],x6G08_MD_903_proj[2],pch=20)
text(x6G08_MD_903_proj[1],x6G08_MD_903_proj[2],pos=4,label="903",col="red")

x6G08_MD_904_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_904.txt")
x6G08_MD_904_proj<-project.pca(x6G08_MD_904_raw,pcadata)
points(x6G08_MD_904_proj[1],x6G08_MD_904_proj[2],pch=20)
text(x6G08_MD_904_proj[1],x6G08_MD_904_proj[2],pos=4,label="904",col="red")

x6G08_MD_905_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_905.txt")
x6G08_MD_905_proj<-project.pca(x6G08_MD_905_raw,pcadata)
points(x6G08_MD_905_proj[1],x6G08_MD_905_proj[2],pch=20)
text(x6G08_MD_905_proj[1],x6G08_MD_905_proj[2],pos=4,label="905",col="red")

x6G08_MD_906_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_906.txt")
x6G08_MD_906_proj<-project.pca(x6G08_MD_906_raw,pcadata)
points(x6G08_MD_906_proj[1],x6G08_MD_906_proj[2],pch=20)
text(x6G08_MD_906_proj[1],x6G08_MD_906_proj[2],pos=4,label="906",col="red")

x6G08_MD_907_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_907.txt")
x6G08_MD_907_proj<-project.pca(x6G08_MD_907_raw,pcadata)
points(x6G08_MD_907_proj[1],x6G08_MD_907_proj[2],pch=20)
text(x6G08_MD_907_proj[1],x6G08_MD_907_proj[2],pos=4,label="907",col="red")

x6G08_MD_908_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_908.txt")
x6G08_MD_908_proj<-project.pca(x6G08_MD_908_raw,pcadata)
points(x6G08_MD_908_proj[1],x6G08_MD_908_proj[2],pch=20)
text(x6G08_MD_908_proj[1],x6G08_MD_908_proj[2],pos=4,label="908",col="red")

x6G08_MD_909_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_909.txt")
x6G08_MD_909_proj<-project.pca(x6G08_MD_909_raw,pcadata)
points(x6G08_MD_909_proj[1],x6G08_MD_909_proj[2],pch=20)
text(x6G08_MD_909_proj[1],x6G08_MD_909_proj[2],pos=4,label="909",col="red")

x6G08_MD_910_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_910.txt")
x6G08_MD_910_proj<-project.pca(x6G08_MD_910_raw,pcadata)
points(x6G08_MD_910_proj[1],x6G08_MD_910_proj[2],pch=20)
text(x6G08_MD_910_proj[1],x6G08_MD_910_proj[2],pos=4,label="910",col="red")

x6G08_MD_911_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_911.txt")
x6G08_MD_911_proj<-project.pca(x6G08_MD_911_raw,pcadata)
points(x6G08_MD_911_proj[1],x6G08_MD_911_proj[2],pch=20)
text(x6G08_MD_911_proj[1],x6G08_MD_911_proj[2],pos=4,label="911",col="red")

x6G08_MD_912_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_912.txt")
x6G08_MD_912_proj<-project.pca(x6G08_MD_912_raw,pcadata)
points(x6G08_MD_912_proj[1],x6G08_MD_912_proj[2],pch=20)
text(x6G08_MD_912_proj[1],x6G08_MD_912_proj[2],pos=4,label="912",col="red")

x6G08_MD_913_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_913.txt")
x6G08_MD_913_proj<-project.pca(x6G08_MD_913_raw,pcadata)
points(x6G08_MD_913_proj[1],x6G08_MD_913_proj[2],pch=20)
text(x6G08_MD_913_proj[1],x6G08_MD_913_proj[2],pos=4,label="913",col="green")

x6G08_MD_914_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_914.txt")
x6G08_MD_914_proj<-project.pca(x6G08_MD_914_raw,pcadata)
points(x6G08_MD_914_proj[1],x6G08_MD_914_proj[2],pch=20)
text(x6G08_MD_914_proj[1],x6G08_MD_914_proj[2],pos=4,label="914",col="red")

x6G08_MD_915_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_915.txt")
x6G08_MD_915_proj<-project.pca(x6G08_MD_915_raw,pcadata)
points(x6G08_MD_915_proj[1],x6G08_MD_915_proj[2],pch=20)
text(x6G08_MD_915_proj[1],x6G08_MD_915_proj[2],pos=4,label="915",col="red")

x6G08_MD_916_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_916.txt")
x6G08_MD_916_proj<-project.pca(x6G08_MD_916_raw,pcadata)
points(x6G08_MD_916_proj[1],x6G08_MD_916_proj[2],pch=20)
text(x6G08_MD_916_proj[1],x6G08_MD_916_proj[2],pos=4,label="916",col="red")

x6G08_MD_917_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_917.txt")
x6G08_MD_917_proj<-project.pca(x6G08_MD_917_raw,pcadata)
points(x6G08_MD_917_proj[1],x6G08_MD_917_proj[2],pch=20)
text(x6G08_MD_917_proj[1],x6G08_MD_917_proj[2],pos=4,label="917",col="red")

x6G08_MD_918_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_918.txt")
x6G08_MD_918_proj<-project.pca(x6G08_MD_918_raw,pcadata)
points(x6G08_MD_918_proj[1],x6G08_MD_918_proj[2],pch=20)
text(x6G08_MD_918_proj[1],x6G08_MD_918_proj[2],pos=4,label="918",col="red")

x6G08_MD_919_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_919.txt")
x6G08_MD_919_proj<-project.pca(x6G08_MD_919_raw,pcadata)
points(x6G08_MD_919_proj[1],x6G08_MD_919_proj[2],pch=20)
text(x6G08_MD_919_proj[1],x6G08_MD_919_proj[2],pos=4,label="919",col="red")

x6G08_MD_920_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_920.txt")
x6G08_MD_920_proj<-project.pca(x6G08_MD_920_raw,pcadata)
points(x6G08_MD_920_proj[1],x6G08_MD_920_proj[2],pch=20)
text(x6G08_MD_920_proj[1],x6G08_MD_920_proj[2],pos=4,label="920",col="red")

x6G08_MD_921_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_921.txt")
x6G08_MD_921_proj<-project.pca(x6G08_MD_921_raw,pcadata)
points(x6G08_MD_921_proj[1],x6G08_MD_921_proj[2],pch=20)
text(x6G08_MD_921_proj[1],x6G08_MD_921_proj[2],pos=4,label="921",col="red")

x6G08_MD_922_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_922.txt")
x6G08_MD_922_proj<-project.pca(x6G08_MD_922_raw,pcadata)
points(x6G08_MD_922_proj[1],x6G08_MD_922_proj[2],pch=20)
text(x6G08_MD_922_proj[1],x6G08_MD_922_proj[2],pos=4,label="922",col="red")

x6G08_MD_923_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_923.txt")
x6G08_MD_923_proj<-project.pca(x6G08_MD_923_raw,pcadata)
points(x6G08_MD_923_proj[1],x6G08_MD_923_proj[2],pch=20)
text(x6G08_MD_923_proj[1],x6G08_MD_923_proj[2],pos=4,label="923",col="red")

x6G08_MD_924_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_924.txt")
x6G08_MD_924_proj<-project.pca(x6G08_MD_924_raw,pcadata)
points(x6G08_MD_924_proj[1],x6G08_MD_924_proj[2],pch=20)
text(x6G08_MD_924_proj[1],x6G08_MD_924_proj[2],pos=4,label="924",col="red")

x6G08_MD_925_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_925.txt")
x6G08_MD_925_proj<-project.pca(x6G08_MD_925_raw,pcadata)
points(x6G08_MD_925_proj[1],x6G08_MD_925_proj[2],pch=20)
text(x6G08_MD_925_proj[1],x6G08_MD_925_proj[2],pos=4,label="925",col="red")

x6G08_MD_926_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_926.txt")
x6G08_MD_926_proj<-project.pca(x6G08_MD_926_raw,pcadata)
points(x6G08_MD_926_proj[1],x6G08_MD_926_proj[2],pch=20)
text(x6G08_MD_926_proj[1],x6G08_MD_926_proj[2],pos=4,label="926",col="red")

x6G08_MD_927_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_927.txt")
x6G08_MD_927_proj<-project.pca(x6G08_MD_927_raw,pcadata)
points(x6G08_MD_927_proj[1],x6G08_MD_927_proj[2],pch=20)
text(x6G08_MD_927_proj[1],x6G08_MD_927_proj[2],pos=4,label="927",col="red")

x6G08_MD_928_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_928.txt")
x6G08_MD_928_proj<-project.pca(x6G08_MD_928_raw,pcadata)
points(x6G08_MD_928_proj[1],x6G08_MD_928_proj[2],pch=20)
text(x6G08_MD_928_proj[1],x6G08_MD_928_proj[2],pos=4,label="928",col="red")

x6G08_MD_929_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_929.txt")
x6G08_MD_929_proj<-project.pca(x6G08_MD_929_raw,pcadata)
points(x6G08_MD_929_proj[1],x6G08_MD_929_proj[2],pch=20)
text(x6G08_MD_929_proj[1],x6G08_MD_929_proj[2],pos=4,label="929",col="red")

x6G08_MD_930_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_930.txt")
x6G08_MD_930_proj<-project.pca(x6G08_MD_930_raw,pcadata)
points(x6G08_MD_930_proj[1],x6G08_MD_930_proj[2],pch=20)
text(x6G08_MD_930_proj[1],x6G08_MD_930_proj[2],pos=4,label="930",col="red")

x6G08_MD_931_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_931.txt")
x6G08_MD_931_proj<-project.pca(x6G08_MD_931_raw,pcadata)
points(x6G08_MD_931_proj[1],x6G08_MD_931_proj[2],pch=20)
text(x6G08_MD_931_proj[1],x6G08_MD_931_proj[2],pos=4,label="931",col="red")

x6G08_MD_932_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_932.txt")
x6G08_MD_932_proj<-project.pca(x6G08_MD_932_raw,pcadata)
points(x6G08_MD_932_proj[1],x6G08_MD_932_proj[2],pch=20)
text(x6G08_MD_932_proj[1],x6G08_MD_932_proj[2],pos=4,label="932",col="red")

x6G08_MD_933_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_933.txt")
x6G08_MD_933_proj<-project.pca(x6G08_MD_933_raw,pcadata)
points(x6G08_MD_933_proj[1],x6G08_MD_933_proj[2],pch=20)
text(x6G08_MD_933_proj[1],x6G08_MD_933_proj[2],pos=4,label="933",col="red")

x6G08_MD_934_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_934.txt")
x6G08_MD_934_proj<-project.pca(x6G08_MD_934_raw,pcadata)
points(x6G08_MD_934_proj[1],x6G08_MD_934_proj[2],pch=20)
text(x6G08_MD_934_proj[1],x6G08_MD_934_proj[2],pos=4,label="934",col="red")

x6G08_MD_935_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_935.txt")
x6G08_MD_935_proj<-project.pca(x6G08_MD_935_raw,pcadata)
points(x6G08_MD_935_proj[1],x6G08_MD_935_proj[2],pch=20)
text(x6G08_MD_935_proj[1],x6G08_MD_935_proj[2],pos=4,label="935",col="red")

x6G08_MD_936_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_936.txt")
x6G08_MD_936_proj<-project.pca(x6G08_MD_936_raw,pcadata)
points(x6G08_MD_936_proj[1],x6G08_MD_936_proj[2],pch=20)
text(x6G08_MD_936_proj[1],x6G08_MD_936_proj[2],pos=4,label="936",col="red")

x6G08_MD_937_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_937.txt")
x6G08_MD_937_proj<-project.pca(x6G08_MD_937_raw,pcadata)
points(x6G08_MD_937_proj[1],x6G08_MD_937_proj[2],pch=20)
text(x6G08_MD_937_proj[1],x6G08_MD_937_proj[2],pos=4,label="937",col="green")

x6G08_MD_938_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_938.txt")
x6G08_MD_938_proj<-project.pca(x6G08_MD_938_raw,pcadata)
points(x6G08_MD_938_proj[1],x6G08_MD_938_proj[2],pch=20)
text(x6G08_MD_938_proj[1],x6G08_MD_938_proj[2],pos=4,label="938",col="green")

x6G08_MD_939_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_939.txt")
x6G08_MD_939_proj<-project.pca(x6G08_MD_939_raw,pcadata)
points(x6G08_MD_939_proj[1],x6G08_MD_939_proj[2],pch=20)
text(x6G08_MD_939_proj[1],x6G08_MD_939_proj[2],pos=4,label="939",col="red")

x6G08_MD_940_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_940.txt")
x6G08_MD_940_proj<-project.pca(x6G08_MD_940_raw,pcadata)
points(x6G08_MD_940_proj[1],x6G08_MD_940_proj[2],pch=20)
text(x6G08_MD_940_proj[1],x6G08_MD_940_proj[2],pos=4,label="940",col="red")

x6G08_MD_941_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_941.txt")
x6G08_MD_941_proj<-project.pca(x6G08_MD_941_raw,pcadata)
points(x6G08_MD_941_proj[1],x6G08_MD_941_proj[2],pch=20)
text(x6G08_MD_941_proj[1],x6G08_MD_941_proj[2],pos=4,label="941",col="red")

x6G08_MD_942_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_942.txt")
x6G08_MD_942_proj<-project.pca(x6G08_MD_942_raw,pcadata)
points(x6G08_MD_942_proj[1],x6G08_MD_942_proj[2],pch=20)
text(x6G08_MD_942_proj[1],x6G08_MD_942_proj[2],pos=4,label="942",col="red")

x6G08_MD_943_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_943.txt")
x6G08_MD_943_proj<-project.pca(x6G08_MD_943_raw,pcadata)
points(x6G08_MD_943_proj[1],x6G08_MD_943_proj[2],pch=20)
text(x6G08_MD_943_proj[1],x6G08_MD_943_proj[2],pos=4,label="943",col="red")

x6G08_MD_944_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_944.txt")
x6G08_MD_944_proj<-project.pca(x6G08_MD_944_raw,pcadata)
points(x6G08_MD_944_proj[1],x6G08_MD_944_proj[2],pch=20)
text(x6G08_MD_944_proj[1],x6G08_MD_944_proj[2],pos=4,label="944",col="red")

x6G08_MD_945_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_945.txt")
x6G08_MD_945_proj<-project.pca(x6G08_MD_945_raw,pcadata)
points(x6G08_MD_945_proj[1],x6G08_MD_945_proj[2],pch=20)
text(x6G08_MD_945_proj[1],x6G08_MD_945_proj[2],pos=4,label="945",col="red")

x6G08_MD_946_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_946.txt")
x6G08_MD_946_proj<-project.pca(x6G08_MD_946_raw,pcadata)
points(x6G08_MD_946_proj[1],x6G08_MD_946_proj[2],pch=20)
text(x6G08_MD_946_proj[1],x6G08_MD_946_proj[2],pos=4,label="946",col="red")

x6G08_MD_947_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_947.txt")
x6G08_MD_947_proj<-project.pca(x6G08_MD_947_raw,pcadata)
points(x6G08_MD_947_proj[1],x6G08_MD_947_proj[2],pch=20)
text(x6G08_MD_947_proj[1],x6G08_MD_947_proj[2],pos=4,label="947",col="red")

x6G08_MD_948_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_948.txt")
x6G08_MD_948_proj<-project.pca(x6G08_MD_948_raw,pcadata)
points(x6G08_MD_948_proj[1],x6G08_MD_948_proj[2],pch=20)
text(x6G08_MD_948_proj[1],x6G08_MD_948_proj[2],pos=4,label="948",col="red")

x6G08_MD_949_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_949.txt")
x6G08_MD_949_proj<-project.pca(x6G08_MD_949_raw,pcadata)
points(x6G08_MD_949_proj[1],x6G08_MD_949_proj[2],pch=20)
text(x6G08_MD_949_proj[1],x6G08_MD_949_proj[2],pos=4,label="949",col="red")

x6G08_MD_950_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_950.txt")
x6G08_MD_950_proj<-project.pca(x6G08_MD_950_raw,pcadata)
points(x6G08_MD_950_proj[1],x6G08_MD_950_proj[2],pch=20)
text(x6G08_MD_950_proj[1],x6G08_MD_950_proj[2],pos=4,label="950",col="red")

x6G08_MD_951_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_951.txt")
x6G08_MD_951_proj<-project.pca(x6G08_MD_951_raw,pcadata)
points(x6G08_MD_951_proj[1],x6G08_MD_951_proj[2],pch=20)
text(x6G08_MD_951_proj[1],x6G08_MD_951_proj[2],pos=4,label="951",col="red")

x6G08_MD_952_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_952.txt")
x6G08_MD_952_proj<-project.pca(x6G08_MD_952_raw,pcadata)
points(x6G08_MD_952_proj[1],x6G08_MD_952_proj[2],pch=20)
text(x6G08_MD_952_proj[1],x6G08_MD_952_proj[2],pos=4,label="952",col="red")

x6G08_MD_953_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_953.txt")
x6G08_MD_953_proj<-project.pca(x6G08_MD_953_raw,pcadata)
points(x6G08_MD_953_proj[1],x6G08_MD_953_proj[2],pch=20)
text(x6G08_MD_953_proj[1],x6G08_MD_953_proj[2],pos=4,label="953",col="red")

x6G08_MD_954_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_954.txt")
x6G08_MD_954_proj<-project.pca(x6G08_MD_954_raw,pcadata)
points(x6G08_MD_954_proj[1],x6G08_MD_954_proj[2],pch=20)
text(x6G08_MD_954_proj[1],x6G08_MD_954_proj[2],pos=4,label="954",col="red")

x6G08_MD_955_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_955.txt")
x6G08_MD_955_proj<-project.pca(x6G08_MD_955_raw,pcadata)
points(x6G08_MD_955_proj[1],x6G08_MD_955_proj[2],pch=20)
text(x6G08_MD_955_proj[1],x6G08_MD_955_proj[2],pos=4,label="955",col="red")

x6G08_MD_956_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_956.txt")
x6G08_MD_956_proj<-project.pca(x6G08_MD_956_raw,pcadata)
points(x6G08_MD_956_proj[1],x6G08_MD_956_proj[2],pch=20)
text(x6G08_MD_956_proj[1],x6G08_MD_956_proj[2],pos=4,label="956",col="red")

x6G08_MD_957_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_957.txt")
x6G08_MD_957_proj<-project.pca(x6G08_MD_957_raw,pcadata)
points(x6G08_MD_957_proj[1],x6G08_MD_957_proj[2],pch=20)
text(x6G08_MD_957_proj[1],x6G08_MD_957_proj[2],pos=4,label="957",col="red")

x6G08_MD_958_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_958.txt")
x6G08_MD_958_proj<-project.pca(x6G08_MD_958_raw,pcadata)
points(x6G08_MD_958_proj[1],x6G08_MD_958_proj[2],pch=20)
text(x6G08_MD_958_proj[1],x6G08_MD_958_proj[2],pos=4,label="958",col="red")

x6G08_MD_959_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_959.txt")
x6G08_MD_959_proj<-project.pca(x6G08_MD_959_raw,pcadata)
points(x6G08_MD_959_proj[1],x6G08_MD_959_proj[2],pch=20)
text(x6G08_MD_959_proj[1],x6G08_MD_959_proj[2],pos=4,label="959",col="red")

x6G08_MD_960_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_960.txt")
x6G08_MD_960_proj<-project.pca(x6G08_MD_960_raw,pcadata)
points(x6G08_MD_960_proj[1],x6G08_MD_960_proj[2],pch=20)
text(x6G08_MD_960_proj[1],x6G08_MD_960_proj[2],pos=4,label="960",col="red")

x6G08_MD_961_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_961.txt")
x6G08_MD_961_proj<-project.pca(x6G08_MD_961_raw,pcadata)
points(x6G08_MD_961_proj[1],x6G08_MD_961_proj[2],pch=20)
text(x6G08_MD_961_proj[1],x6G08_MD_961_proj[2],pos=4,label="961",col="red")

x6G08_MD_962_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_962.txt")
x6G08_MD_962_proj<-project.pca(x6G08_MD_962_raw,pcadata)
points(x6G08_MD_962_proj[1],x6G08_MD_962_proj[2],pch=20)
text(x6G08_MD_962_proj[1],x6G08_MD_962_proj[2],pos=4,label="962",col="red")

x6G08_MD_963_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_963.txt")
x6G08_MD_963_proj<-project.pca(x6G08_MD_963_raw,pcadata)
points(x6G08_MD_963_proj[1],x6G08_MD_963_proj[2],pch=20)
text(x6G08_MD_963_proj[1],x6G08_MD_963_proj[2],pos=4,label="963",col="red")

x6G08_MD_964_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_964.txt")
x6G08_MD_964_proj<-project.pca(x6G08_MD_964_raw,pcadata)
points(x6G08_MD_964_proj[1],x6G08_MD_964_proj[2],pch=20)
text(x6G08_MD_964_proj[1],x6G08_MD_964_proj[2],pos=4,label="964",col="red")

x6G08_MD_965_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_965.txt")
x6G08_MD_965_proj<-project.pca(x6G08_MD_965_raw,pcadata)
points(x6G08_MD_965_proj[1],x6G08_MD_965_proj[2],pch=20)
text(x6G08_MD_965_proj[1],x6G08_MD_965_proj[2],pos=4,label="965",col="red")

x6G08_MD_966_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_966.txt")
x6G08_MD_966_proj<-project.pca(x6G08_MD_966_raw,pcadata)
points(x6G08_MD_966_proj[1],x6G08_MD_966_proj[2],pch=20)
text(x6G08_MD_966_proj[1],x6G08_MD_966_proj[2],pos=4,label="966",col="red")

x6G08_MD_967_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_967.txt")
x6G08_MD_967_proj<-project.pca(x6G08_MD_967_raw,pcadata)
points(x6G08_MD_967_proj[1],x6G08_MD_967_proj[2],pch=20)
text(x6G08_MD_967_proj[1],x6G08_MD_967_proj[2],pos=4,label="967",col="red")

x6G08_MD_968_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_968.txt")
x6G08_MD_968_proj<-project.pca(x6G08_MD_968_raw,pcadata)
points(x6G08_MD_968_proj[1],x6G08_MD_968_proj[2],pch=20)
text(x6G08_MD_968_proj[1],x6G08_MD_968_proj[2],pos=4,label="968",col="red")

x6G08_MD_969_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_969.txt")
x6G08_MD_969_proj<-project.pca(x6G08_MD_969_raw,pcadata)
points(x6G08_MD_969_proj[1],x6G08_MD_969_proj[2],pch=20)
text(x6G08_MD_969_proj[1],x6G08_MD_969_proj[2],pos=4,label="969",col="red")

x6G08_MD_970_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_970.txt")
x6G08_MD_970_proj<-project.pca(x6G08_MD_970_raw,pcadata)
points(x6G08_MD_970_proj[1],x6G08_MD_970_proj[2],pch=20)
text(x6G08_MD_970_proj[1],x6G08_MD_970_proj[2],pos=4,label="970",col="red")

x6G08_MD_971_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_971.txt")
x6G08_MD_971_proj<-project.pca(x6G08_MD_971_raw,pcadata)
points(x6G08_MD_971_proj[1],x6G08_MD_971_proj[2],pch=20)
text(x6G08_MD_971_proj[1],x6G08_MD_971_proj[2],pos=4,label="971",col="red")

x6G08_MD_972_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_972.txt")
x6G08_MD_972_proj<-project.pca(x6G08_MD_972_raw,pcadata)
points(x6G08_MD_972_proj[1],x6G08_MD_972_proj[2],pch=20)
text(x6G08_MD_972_proj[1],x6G08_MD_972_proj[2],pos=4,label="972",col="red")

x6G08_MD_973_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_973.txt")
x6G08_MD_973_proj<-project.pca(x6G08_MD_973_raw,pcadata)
points(x6G08_MD_973_proj[1],x6G08_MD_973_proj[2],pch=20)
text(x6G08_MD_973_proj[1],x6G08_MD_973_proj[2],pos=4,label="973",col="red")

x6G08_MD_974_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_974.txt")
x6G08_MD_974_proj<-project.pca(x6G08_MD_974_raw,pcadata)
points(x6G08_MD_974_proj[1],x6G08_MD_974_proj[2],pch=20)
text(x6G08_MD_974_proj[1],x6G08_MD_974_proj[2],pos=4,label="974",col="red")

x6G08_MD_975_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_975.txt")
x6G08_MD_975_proj<-project.pca(x6G08_MD_975_raw,pcadata)
points(x6G08_MD_975_proj[1],x6G08_MD_975_proj[2],pch=20)
text(x6G08_MD_975_proj[1],x6G08_MD_975_proj[2],pos=4,label="975",col="red")

x6G08_MD_976_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_976.txt")
x6G08_MD_976_proj<-project.pca(x6G08_MD_976_raw,pcadata)
points(x6G08_MD_976_proj[1],x6G08_MD_976_proj[2],pch=20)
text(x6G08_MD_976_proj[1],x6G08_MD_976_proj[2],pos=4,label="976",col="red")

x6G08_MD_977_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_977.txt")
x6G08_MD_977_proj<-project.pca(x6G08_MD_977_raw,pcadata)
points(x6G08_MD_977_proj[1],x6G08_MD_977_proj[2],pch=20)
text(x6G08_MD_977_proj[1],x6G08_MD_977_proj[2],pos=4,label="977",col="red")

x6G08_MD_978_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_978.txt")
x6G08_MD_978_proj<-project.pca(x6G08_MD_978_raw,pcadata)
points(x6G08_MD_978_proj[1],x6G08_MD_978_proj[2],pch=20)
text(x6G08_MD_978_proj[1],x6G08_MD_978_proj[2],pos=4,label="978",col="red")

x6G08_MD_979_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_979.txt")
x6G08_MD_979_proj<-project.pca(x6G08_MD_979_raw,pcadata)
points(x6G08_MD_979_proj[1],x6G08_MD_979_proj[2],pch=20)
text(x6G08_MD_979_proj[1],x6G08_MD_979_proj[2],pos=4,label="979",col="red")

x6G08_MD_980_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_980.txt")
x6G08_MD_980_proj<-project.pca(x6G08_MD_980_raw,pcadata)
points(x6G08_MD_980_proj[1],x6G08_MD_980_proj[2],pch=20)
text(x6G08_MD_980_proj[1],x6G08_MD_980_proj[2],pos=4,label="980",col="green")

x6G08_MD_981_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_981.txt")
x6G08_MD_981_proj<-project.pca(x6G08_MD_981_raw,pcadata)
points(x6G08_MD_981_proj[1],x6G08_MD_981_proj[2],pch=20)
text(x6G08_MD_981_proj[1],x6G08_MD_981_proj[2],pos=4,label="981",col="red")

x6G08_MD_982_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_982.txt")
x6G08_MD_982_proj<-project.pca(x6G08_MD_982_raw,pcadata)
points(x6G08_MD_982_proj[1],x6G08_MD_982_proj[2],pch=20)
text(x6G08_MD_982_proj[1],x6G08_MD_982_proj[2],pos=4,label="982",col="red")

x6G08_MD_983_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_983.txt")
x6G08_MD_983_proj<-project.pca(x6G08_MD_983_raw,pcadata)
points(x6G08_MD_983_proj[1],x6G08_MD_983_proj[2],pch=20)
text(x6G08_MD_983_proj[1],x6G08_MD_983_proj[2],pos=4,label="983",col="green")

x6G08_MD_984_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_984.txt")
x6G08_MD_984_proj<-project.pca(x6G08_MD_984_raw,pcadata)
points(x6G08_MD_984_proj[1],x6G08_MD_984_proj[2],pch=20)
text(x6G08_MD_984_proj[1],x6G08_MD_984_proj[2],pos=4,label="984",col="green")

x6G08_MD_985_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_985.txt")
x6G08_MD_985_proj<-project.pca(x6G08_MD_985_raw,pcadata)
points(x6G08_MD_985_proj[1],x6G08_MD_985_proj[2],pch=20)
text(x6G08_MD_985_proj[1],x6G08_MD_985_proj[2],pos=4,label="985",col="red")

x6G08_MD_986_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_986.txt")
x6G08_MD_986_proj<-project.pca(x6G08_MD_986_raw,pcadata)
points(x6G08_MD_986_proj[1],x6G08_MD_986_proj[2],pch=20)
text(x6G08_MD_986_proj[1],x6G08_MD_986_proj[2],pos=4,label="986",col="red")

x6G08_MD_987_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_987.txt")
x6G08_MD_987_proj<-project.pca(x6G08_MD_987_raw,pcadata)
points(x6G08_MD_987_proj[1],x6G08_MD_987_proj[2],pch=20)
text(x6G08_MD_987_proj[1],x6G08_MD_987_proj[2],pos=4,label="987",col="green")

x6G08_MD_988_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_988.txt")
x6G08_MD_988_proj<-project.pca(x6G08_MD_988_raw,pcadata)
points(x6G08_MD_988_proj[1],x6G08_MD_988_proj[2],pch=20)
text(x6G08_MD_988_proj[1],x6G08_MD_988_proj[2],pos=4,label="988",col="green")

x6G08_MD_989_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_989.txt")
x6G08_MD_989_proj<-project.pca(x6G08_MD_989_raw,pcadata)
points(x6G08_MD_989_proj[1],x6G08_MD_989_proj[2],pch=20)
text(x6G08_MD_989_proj[1],x6G08_MD_989_proj[2],pos=4,label="989",col="green")

x6G08_MD_990_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_990.txt")
x6G08_MD_990_proj<-project.pca(x6G08_MD_990_raw,pcadata)
points(x6G08_MD_990_proj[1],x6G08_MD_990_proj[2],pch=20)
text(x6G08_MD_990_proj[1],x6G08_MD_990_proj[2],pos=4,label="990",col="red")

x6G08_MD_991_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_991.txt")
x6G08_MD_991_proj<-project.pca(x6G08_MD_991_raw,pcadata)
points(x6G08_MD_991_proj[1],x6G08_MD_991_proj[2],pch=20)
text(x6G08_MD_991_proj[1],x6G08_MD_991_proj[2],pos=4,label="991",col="red")

x6G08_MD_992_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_992.txt")
x6G08_MD_992_proj<-project.pca(x6G08_MD_992_raw,pcadata)
points(x6G08_MD_992_proj[1],x6G08_MD_992_proj[2],pch=20)
text(x6G08_MD_992_proj[1],x6G08_MD_992_proj[2],pos=4,label="992",col="green")

x6G08_MD_993_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_993.txt")
x6G08_MD_993_proj<-project.pca(x6G08_MD_993_raw,pcadata)
points(x6G08_MD_993_proj[1],x6G08_MD_993_proj[2],pch=20)
text(x6G08_MD_993_proj[1],x6G08_MD_993_proj[2],pos=4,label="993",col="green")

x6G08_MD_994_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_994.txt")
x6G08_MD_994_proj<-project.pca(x6G08_MD_994_raw,pcadata)
points(x6G08_MD_994_proj[1],x6G08_MD_994_proj[2],pch=20)
text(x6G08_MD_994_proj[1],x6G08_MD_994_proj[2],pos=4,label="994",col="red")

x6G08_MD_995_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_995.txt")
x6G08_MD_995_proj<-project.pca(x6G08_MD_995_raw,pcadata)
points(x6G08_MD_995_proj[1],x6G08_MD_995_proj[2],pch=20)
text(x6G08_MD_995_proj[1],x6G08_MD_995_proj[2],pos=4,label="995",col="green")

x6G08_MD_996_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_996.txt")
x6G08_MD_996_proj<-project.pca(x6G08_MD_996_raw,pcadata)
points(x6G08_MD_996_proj[1],x6G08_MD_996_proj[2],pch=20)
text(x6G08_MD_996_proj[1],x6G08_MD_996_proj[2],pos=4,label="996",col="green")

x6G08_MD_997_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_997.txt")
x6G08_MD_997_proj<-project.pca(x6G08_MD_997_raw,pcadata)
points(x6G08_MD_997_proj[1],x6G08_MD_997_proj[2],pch=20)
text(x6G08_MD_997_proj[1],x6G08_MD_997_proj[2],pos=4,label="997",col="green")

x6G08_MD_998_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_998.txt")
x6G08_MD_998_proj<-project.pca(x6G08_MD_998_raw,pcadata)
points(x6G08_MD_998_proj[1],x6G08_MD_998_proj[2],pch=20)
text(x6G08_MD_998_proj[1],x6G08_MD_998_proj[2],pos=4,label="998",col="blue")

x6G08_MD_999_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_999.txt")
x6G08_MD_999_proj<-project.pca(x6G08_MD_999_raw,pcadata)
points(x6G08_MD_999_proj[1],x6G08_MD_999_proj[2],pch=20)
text(x6G08_MD_999_proj[1],x6G08_MD_999_proj[2],pos=4,label="999",col="green")

x6G08_MD_1000_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_1000.txt")
x6G08_MD_1000_proj<-project.pca(x6G08_MD_1000_raw,pcadata)
points(x6G08_MD_1000_proj[1],x6G08_MD_1000_proj[2],pch=20)
text(x6G08_MD_1000_proj[1],x6G08_MD_1000_proj[2],pos=4,label="1000",col="red")

plot(PC1,PC3,xlim=c(-150,200),ylim=c(-30,50))

xCA_6G08_fab_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/CA_6G08_fab.txt")
xCA_6G08_fab_proj<-project.pca(xCA_6G08_fab_raw,pcadata)
points(xCA_6G08_fab_proj[1],xCA_6G08_fab_proj[3],pch=20)
text(xCA_6G08_fab_proj[1],xCA_6G08_fab_proj[3],pos=4,label="6G08",col="purple")

x6G08_MD_1_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_1.txt")
x6G08_MD_1_proj<-project.pca(x6G08_MD_1_raw,pcadata)
points(x6G08_MD_1_proj[1],x6G08_MD_1_proj[3],pch=20)
text(x6G08_MD_1_proj[1],x6G08_MD_1_proj[3],pos=4,label="1",col="blue")

x6G08_MD_2_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_2.txt")
x6G08_MD_2_proj<-project.pca(x6G08_MD_2_raw,pcadata)
points(x6G08_MD_2_proj[1],x6G08_MD_2_proj[3],pch=20)
text(x6G08_MD_2_proj[1],x6G08_MD_2_proj[3],pos=4,label="2",col="blue")

x6G08_MD_3_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_3.txt")
x6G08_MD_3_proj<-project.pca(x6G08_MD_3_raw,pcadata)
points(x6G08_MD_3_proj[1],x6G08_MD_3_proj[3],pch=20)
text(x6G08_MD_3_proj[1],x6G08_MD_3_proj[3],pos=4,label="3",col="blue")

x6G08_MD_4_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_4.txt")
x6G08_MD_4_proj<-project.pca(x6G08_MD_4_raw,pcadata)
points(x6G08_MD_4_proj[1],x6G08_MD_4_proj[3],pch=20)
text(x6G08_MD_4_proj[1],x6G08_MD_4_proj[3],pos=4,label="4",col="blue")

x6G08_MD_5_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_5.txt")
x6G08_MD_5_proj<-project.pca(x6G08_MD_5_raw,pcadata)
points(x6G08_MD_5_proj[1],x6G08_MD_5_proj[3],pch=20)
text(x6G08_MD_5_proj[1],x6G08_MD_5_proj[3],pos=4,label="5",col="blue")

x6G08_MD_6_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_6.txt")
x6G08_MD_6_proj<-project.pca(x6G08_MD_6_raw,pcadata)
points(x6G08_MD_6_proj[1],x6G08_MD_6_proj[3],pch=20)
text(x6G08_MD_6_proj[1],x6G08_MD_6_proj[3],pos=4,label="6",col="black")

x6G08_MD_7_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_7.txt")
x6G08_MD_7_proj<-project.pca(x6G08_MD_7_raw,pcadata)
points(x6G08_MD_7_proj[1],x6G08_MD_7_proj[3],pch=20)
text(x6G08_MD_7_proj[1],x6G08_MD_7_proj[3],pos=4,label="7",col="blue")

x6G08_MD_8_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_8.txt")
x6G08_MD_8_proj<-project.pca(x6G08_MD_8_raw,pcadata)
points(x6G08_MD_8_proj[1],x6G08_MD_8_proj[3],pch=20)
text(x6G08_MD_8_proj[1],x6G08_MD_8_proj[3],pos=4,label="8",col="blue")

x6G08_MD_9_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_9.txt")
x6G08_MD_9_proj<-project.pca(x6G08_MD_9_raw,pcadata)
points(x6G08_MD_9_proj[1],x6G08_MD_9_proj[3],pch=20)
text(x6G08_MD_9_proj[1],x6G08_MD_9_proj[3],pos=4,label="9",col="green")

x6G08_MD_10_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_10.txt")
x6G08_MD_10_proj<-project.pca(x6G08_MD_10_raw,pcadata)
points(x6G08_MD_10_proj[1],x6G08_MD_10_proj[3],pch=20)
text(x6G08_MD_10_proj[1],x6G08_MD_10_proj[3],pos=4,label="10",col="red")

x6G08_MD_11_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_11.txt")
x6G08_MD_11_proj<-project.pca(x6G08_MD_11_raw,pcadata)
points(x6G08_MD_11_proj[1],x6G08_MD_11_proj[3],pch=20)
text(x6G08_MD_11_proj[1],x6G08_MD_11_proj[3],pos=4,label="11",col="blue")

x6G08_MD_12_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_12.txt")
x6G08_MD_12_proj<-project.pca(x6G08_MD_12_raw,pcadata)
points(x6G08_MD_12_proj[1],x6G08_MD_12_proj[3],pch=20)
text(x6G08_MD_12_proj[1],x6G08_MD_12_proj[3],pos=4,label="12",col="black")

x6G08_MD_13_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_13.txt")
x6G08_MD_13_proj<-project.pca(x6G08_MD_13_raw,pcadata)
points(x6G08_MD_13_proj[1],x6G08_MD_13_proj[3],pch=20)
text(x6G08_MD_13_proj[1],x6G08_MD_13_proj[3],pos=4,label="13",col="black")

x6G08_MD_14_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_14.txt")
x6G08_MD_14_proj<-project.pca(x6G08_MD_14_raw,pcadata)
points(x6G08_MD_14_proj[1],x6G08_MD_14_proj[3],pch=20)
text(x6G08_MD_14_proj[1],x6G08_MD_14_proj[3],pos=4,label="14",col="blue")

x6G08_MD_15_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_15.txt")
x6G08_MD_15_proj<-project.pca(x6G08_MD_15_raw,pcadata)
points(x6G08_MD_15_proj[1],x6G08_MD_15_proj[3],pch=20)
text(x6G08_MD_15_proj[1],x6G08_MD_15_proj[3],pos=4,label="15",col="black")

x6G08_MD_16_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_16.txt")
x6G08_MD_16_proj<-project.pca(x6G08_MD_16_raw,pcadata)
points(x6G08_MD_16_proj[1],x6G08_MD_16_proj[3],pch=20)
text(x6G08_MD_16_proj[1],x6G08_MD_16_proj[3],pos=4,label="16",col="black")

x6G08_MD_17_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_17.txt")
x6G08_MD_17_proj<-project.pca(x6G08_MD_17_raw,pcadata)
points(x6G08_MD_17_proj[1],x6G08_MD_17_proj[3],pch=20)
text(x6G08_MD_17_proj[1],x6G08_MD_17_proj[3],pos=4,label="17",col="blue")

x6G08_MD_18_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_18.txt")
x6G08_MD_18_proj<-project.pca(x6G08_MD_18_raw,pcadata)
points(x6G08_MD_18_proj[1],x6G08_MD_18_proj[3],pch=20)
text(x6G08_MD_18_proj[1],x6G08_MD_18_proj[3],pos=4,label="18",col="black")

x6G08_MD_19_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_19.txt")
x6G08_MD_19_proj<-project.pca(x6G08_MD_19_raw,pcadata)
points(x6G08_MD_19_proj[1],x6G08_MD_19_proj[3],pch=20)
text(x6G08_MD_19_proj[1],x6G08_MD_19_proj[3],pos=4,label="19",col="black")

x6G08_MD_20_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_20.txt")
x6G08_MD_20_proj<-project.pca(x6G08_MD_20_raw,pcadata)
points(x6G08_MD_20_proj[1],x6G08_MD_20_proj[3],pch=20)
text(x6G08_MD_20_proj[1],x6G08_MD_20_proj[3],pos=4,label="20",col="black")

x6G08_MD_21_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_21.txt")
x6G08_MD_21_proj<-project.pca(x6G08_MD_21_raw,pcadata)
points(x6G08_MD_21_proj[1],x6G08_MD_21_proj[3],pch=20)
text(x6G08_MD_21_proj[1],x6G08_MD_21_proj[3],pos=4,label="21",col="blue")

x6G08_MD_22_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_22.txt")
x6G08_MD_22_proj<-project.pca(x6G08_MD_22_raw,pcadata)
points(x6G08_MD_22_proj[1],x6G08_MD_22_proj[3],pch=20)
text(x6G08_MD_22_proj[1],x6G08_MD_22_proj[3],pos=4,label="22",col="blue")

x6G08_MD_23_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_23.txt")
x6G08_MD_23_proj<-project.pca(x6G08_MD_23_raw,pcadata)
points(x6G08_MD_23_proj[1],x6G08_MD_23_proj[3],pch=20)
text(x6G08_MD_23_proj[1],x6G08_MD_23_proj[3],pos=4,label="23",col="blue")

x6G08_MD_24_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_24.txt")
x6G08_MD_24_proj<-project.pca(x6G08_MD_24_raw,pcadata)
points(x6G08_MD_24_proj[1],x6G08_MD_24_proj[3],pch=20)
text(x6G08_MD_24_proj[1],x6G08_MD_24_proj[3],pos=4,label="24",col="blue")

x6G08_MD_25_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_25.txt")
x6G08_MD_25_proj<-project.pca(x6G08_MD_25_raw,pcadata)
points(x6G08_MD_25_proj[1],x6G08_MD_25_proj[3],pch=20)
text(x6G08_MD_25_proj[1],x6G08_MD_25_proj[3],pos=4,label="25",col="blue")

x6G08_MD_26_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_26.txt")
x6G08_MD_26_proj<-project.pca(x6G08_MD_26_raw,pcadata)
points(x6G08_MD_26_proj[1],x6G08_MD_26_proj[3],pch=20)
text(x6G08_MD_26_proj[1],x6G08_MD_26_proj[3],pos=4,label="26",col="blue")

x6G08_MD_27_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_27.txt")
x6G08_MD_27_proj<-project.pca(x6G08_MD_27_raw,pcadata)
points(x6G08_MD_27_proj[1],x6G08_MD_27_proj[3],pch=20)
text(x6G08_MD_27_proj[1],x6G08_MD_27_proj[3],pos=4,label="27",col="blue")

x6G08_MD_28_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_28.txt")
x6G08_MD_28_proj<-project.pca(x6G08_MD_28_raw,pcadata)
points(x6G08_MD_28_proj[1],x6G08_MD_28_proj[3],pch=20)
text(x6G08_MD_28_proj[1],x6G08_MD_28_proj[3],pos=4,label="28",col="blue")

x6G08_MD_29_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_29.txt")
x6G08_MD_29_proj<-project.pca(x6G08_MD_29_raw,pcadata)
points(x6G08_MD_29_proj[1],x6G08_MD_29_proj[3],pch=20)
text(x6G08_MD_29_proj[1],x6G08_MD_29_proj[3],pos=4,label="29",col="blue")

x6G08_MD_30_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_30.txt")
x6G08_MD_30_proj<-project.pca(x6G08_MD_30_raw,pcadata)
points(x6G08_MD_30_proj[1],x6G08_MD_30_proj[3],pch=20)
text(x6G08_MD_30_proj[1],x6G08_MD_30_proj[3],pos=4,label="30",col="black")

x6G08_MD_31_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_31.txt")
x6G08_MD_31_proj<-project.pca(x6G08_MD_31_raw,pcadata)
points(x6G08_MD_31_proj[1],x6G08_MD_31_proj[3],pch=20)
text(x6G08_MD_31_proj[1],x6G08_MD_31_proj[3],pos=4,label="31",col="black")

x6G08_MD_32_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_32.txt")
x6G08_MD_32_proj<-project.pca(x6G08_MD_32_raw,pcadata)
points(x6G08_MD_32_proj[1],x6G08_MD_32_proj[3],pch=20)
text(x6G08_MD_32_proj[1],x6G08_MD_32_proj[3],pos=4,label="32",col="black")

x6G08_MD_33_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_33.txt")
x6G08_MD_33_proj<-project.pca(x6G08_MD_33_raw,pcadata)
points(x6G08_MD_33_proj[1],x6G08_MD_33_proj[3],pch=20)
text(x6G08_MD_33_proj[1],x6G08_MD_33_proj[3],pos=4,label="33",col="black")

x6G08_MD_34_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_34.txt")
x6G08_MD_34_proj<-project.pca(x6G08_MD_34_raw,pcadata)
points(x6G08_MD_34_proj[1],x6G08_MD_34_proj[3],pch=20)
text(x6G08_MD_34_proj[1],x6G08_MD_34_proj[3],pos=4,label="34",col="blue")

x6G08_MD_35_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_35.txt")
x6G08_MD_35_proj<-project.pca(x6G08_MD_35_raw,pcadata)
points(x6G08_MD_35_proj[1],x6G08_MD_35_proj[3],pch=20)
text(x6G08_MD_35_proj[1],x6G08_MD_35_proj[3],pos=4,label="35",col="black")

x6G08_MD_36_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_36.txt")
x6G08_MD_36_proj<-project.pca(x6G08_MD_36_raw,pcadata)
points(x6G08_MD_36_proj[1],x6G08_MD_36_proj[3],pch=20)
text(x6G08_MD_36_proj[1],x6G08_MD_36_proj[3],pos=4,label="36",col="black")

x6G08_MD_37_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_37.txt")
x6G08_MD_37_proj<-project.pca(x6G08_MD_37_raw,pcadata)
points(x6G08_MD_37_proj[1],x6G08_MD_37_proj[3],pch=20)
text(x6G08_MD_37_proj[1],x6G08_MD_37_proj[3],pos=4,label="37",col="blue")

x6G08_MD_38_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_38.txt")
x6G08_MD_38_proj<-project.pca(x6G08_MD_38_raw,pcadata)
points(x6G08_MD_38_proj[1],x6G08_MD_38_proj[3],pch=20)
text(x6G08_MD_38_proj[1],x6G08_MD_38_proj[3],pos=4,label="38",col="blue")

x6G08_MD_39_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_39.txt")
x6G08_MD_39_proj<-project.pca(x6G08_MD_39_raw,pcadata)
points(x6G08_MD_39_proj[1],x6G08_MD_39_proj[3],pch=20)
text(x6G08_MD_39_proj[1],x6G08_MD_39_proj[3],pos=4,label="39",col="blue")

x6G08_MD_40_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_40.txt")
x6G08_MD_40_proj<-project.pca(x6G08_MD_40_raw,pcadata)
points(x6G08_MD_40_proj[1],x6G08_MD_40_proj[3],pch=20)
text(x6G08_MD_40_proj[1],x6G08_MD_40_proj[3],pos=4,label="40",col="black")

x6G08_MD_41_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_41.txt")
x6G08_MD_41_proj<-project.pca(x6G08_MD_41_raw,pcadata)
points(x6G08_MD_41_proj[1],x6G08_MD_41_proj[3],pch=20)
text(x6G08_MD_41_proj[1],x6G08_MD_41_proj[3],pos=4,label="41",col="black")

x6G08_MD_42_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_42.txt")
x6G08_MD_42_proj<-project.pca(x6G08_MD_42_raw,pcadata)
points(x6G08_MD_42_proj[1],x6G08_MD_42_proj[3],pch=20)
text(x6G08_MD_42_proj[1],x6G08_MD_42_proj[3],pos=4,label="42",col="black")

x6G08_MD_43_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_43.txt")
x6G08_MD_43_proj<-project.pca(x6G08_MD_43_raw,pcadata)
points(x6G08_MD_43_proj[1],x6G08_MD_43_proj[3],pch=20)
text(x6G08_MD_43_proj[1],x6G08_MD_43_proj[3],pos=4,label="43",col="black")

x6G08_MD_44_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_44.txt")
x6G08_MD_44_proj<-project.pca(x6G08_MD_44_raw,pcadata)
points(x6G08_MD_44_proj[1],x6G08_MD_44_proj[3],pch=20)
text(x6G08_MD_44_proj[1],x6G08_MD_44_proj[3],pos=4,label="44",col="black")

x6G08_MD_45_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_45.txt")
x6G08_MD_45_proj<-project.pca(x6G08_MD_45_raw,pcadata)
points(x6G08_MD_45_proj[1],x6G08_MD_45_proj[3],pch=20)
text(x6G08_MD_45_proj[1],x6G08_MD_45_proj[3],pos=4,label="45",col="black")

x6G08_MD_46_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_46.txt")
x6G08_MD_46_proj<-project.pca(x6G08_MD_46_raw,pcadata)
points(x6G08_MD_46_proj[1],x6G08_MD_46_proj[3],pch=20)
text(x6G08_MD_46_proj[1],x6G08_MD_46_proj[3],pos=4,label="46",col="black")

x6G08_MD_47_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_47.txt")
x6G08_MD_47_proj<-project.pca(x6G08_MD_47_raw,pcadata)
points(x6G08_MD_47_proj[1],x6G08_MD_47_proj[3],pch=20)
text(x6G08_MD_47_proj[1],x6G08_MD_47_proj[3],pos=4,label="47",col="black")

x6G08_MD_48_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_48.txt")
x6G08_MD_48_proj<-project.pca(x6G08_MD_48_raw,pcadata)
points(x6G08_MD_48_proj[1],x6G08_MD_48_proj[3],pch=20)
text(x6G08_MD_48_proj[1],x6G08_MD_48_proj[3],pos=4,label="48",col="black")

x6G08_MD_49_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_49.txt")
x6G08_MD_49_proj<-project.pca(x6G08_MD_49_raw,pcadata)
points(x6G08_MD_49_proj[1],x6G08_MD_49_proj[3],pch=20)
text(x6G08_MD_49_proj[1],x6G08_MD_49_proj[3],pos=4,label="49",col="black")

x6G08_MD_50_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_50.txt")
x6G08_MD_50_proj<-project.pca(x6G08_MD_50_raw,pcadata)
points(x6G08_MD_50_proj[1],x6G08_MD_50_proj[3],pch=20)
text(x6G08_MD_50_proj[1],x6G08_MD_50_proj[3],pos=4,label="50",col="black")

x6G08_MD_51_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_51.txt")
x6G08_MD_51_proj<-project.pca(x6G08_MD_51_raw,pcadata)
points(x6G08_MD_51_proj[1],x6G08_MD_51_proj[3],pch=20)
text(x6G08_MD_51_proj[1],x6G08_MD_51_proj[3],pos=4,label="51",col="black")

x6G08_MD_52_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_52.txt")
x6G08_MD_52_proj<-project.pca(x6G08_MD_52_raw,pcadata)
points(x6G08_MD_52_proj[1],x6G08_MD_52_proj[3],pch=20)
text(x6G08_MD_52_proj[1],x6G08_MD_52_proj[3],pos=4,label="52",col="black")

x6G08_MD_53_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_53.txt")
x6G08_MD_53_proj<-project.pca(x6G08_MD_53_raw,pcadata)
points(x6G08_MD_53_proj[1],x6G08_MD_53_proj[3],pch=20)
text(x6G08_MD_53_proj[1],x6G08_MD_53_proj[3],pos=4,label="53",col="black")

x6G08_MD_54_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_54.txt")
x6G08_MD_54_proj<-project.pca(x6G08_MD_54_raw,pcadata)
points(x6G08_MD_54_proj[1],x6G08_MD_54_proj[3],pch=20)
text(x6G08_MD_54_proj[1],x6G08_MD_54_proj[3],pos=4,label="54",col="black")

x6G08_MD_55_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_55.txt")
x6G08_MD_55_proj<-project.pca(x6G08_MD_55_raw,pcadata)
points(x6G08_MD_55_proj[1],x6G08_MD_55_proj[3],pch=20)
text(x6G08_MD_55_proj[1],x6G08_MD_55_proj[3],pos=4,label="55",col="black")

x6G08_MD_56_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_56.txt")
x6G08_MD_56_proj<-project.pca(x6G08_MD_56_raw,pcadata)
points(x6G08_MD_56_proj[1],x6G08_MD_56_proj[3],pch=20)
text(x6G08_MD_56_proj[1],x6G08_MD_56_proj[3],pos=4,label="56",col="black")

x6G08_MD_57_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_57.txt")
x6G08_MD_57_proj<-project.pca(x6G08_MD_57_raw,pcadata)
points(x6G08_MD_57_proj[1],x6G08_MD_57_proj[3],pch=20)
text(x6G08_MD_57_proj[1],x6G08_MD_57_proj[3],pos=4,label="57",col="black")

x6G08_MD_58_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_58.txt")
x6G08_MD_58_proj<-project.pca(x6G08_MD_58_raw,pcadata)
points(x6G08_MD_58_proj[1],x6G08_MD_58_proj[3],pch=20)
text(x6G08_MD_58_proj[1],x6G08_MD_58_proj[3],pos=4,label="58",col="blue")

x6G08_MD_59_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_59.txt")
x6G08_MD_59_proj<-project.pca(x6G08_MD_59_raw,pcadata)
points(x6G08_MD_59_proj[1],x6G08_MD_59_proj[3],pch=20)
text(x6G08_MD_59_proj[1],x6G08_MD_59_proj[3],pos=4,label="59",col="black")

x6G08_MD_60_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_60.txt")
x6G08_MD_60_proj<-project.pca(x6G08_MD_60_raw,pcadata)
points(x6G08_MD_60_proj[1],x6G08_MD_60_proj[3],pch=20)
text(x6G08_MD_60_proj[1],x6G08_MD_60_proj[3],pos=4,label="60",col="black")

x6G08_MD_61_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_61.txt")
x6G08_MD_61_proj<-project.pca(x6G08_MD_61_raw,pcadata)
points(x6G08_MD_61_proj[1],x6G08_MD_61_proj[3],pch=20)
text(x6G08_MD_61_proj[1],x6G08_MD_61_proj[3],pos=4,label="61",col="black")

x6G08_MD_62_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_62.txt")
x6G08_MD_62_proj<-project.pca(x6G08_MD_62_raw,pcadata)
points(x6G08_MD_62_proj[1],x6G08_MD_62_proj[3],pch=20)
text(x6G08_MD_62_proj[1],x6G08_MD_62_proj[3],pos=4,label="62",col="black")

x6G08_MD_63_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_63.txt")
x6G08_MD_63_proj<-project.pca(x6G08_MD_63_raw,pcadata)
points(x6G08_MD_63_proj[1],x6G08_MD_63_proj[3],pch=20)
text(x6G08_MD_63_proj[1],x6G08_MD_63_proj[3],pos=4,label="63",col="black")

x6G08_MD_64_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_64.txt")
x6G08_MD_64_proj<-project.pca(x6G08_MD_64_raw,pcadata)
points(x6G08_MD_64_proj[1],x6G08_MD_64_proj[3],pch=20)
text(x6G08_MD_64_proj[1],x6G08_MD_64_proj[3],pos=4,label="64",col="black")

x6G08_MD_65_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_65.txt")
x6G08_MD_65_proj<-project.pca(x6G08_MD_65_raw,pcadata)
points(x6G08_MD_65_proj[1],x6G08_MD_65_proj[3],pch=20)
text(x6G08_MD_65_proj[1],x6G08_MD_65_proj[3],pos=4,label="65",col="black")

x6G08_MD_66_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_66.txt")
x6G08_MD_66_proj<-project.pca(x6G08_MD_66_raw,pcadata)
points(x6G08_MD_66_proj[1],x6G08_MD_66_proj[3],pch=20)
text(x6G08_MD_66_proj[1],x6G08_MD_66_proj[3],pos=4,label="66",col="black")

x6G08_MD_67_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_67.txt")
x6G08_MD_67_proj<-project.pca(x6G08_MD_67_raw,pcadata)
points(x6G08_MD_67_proj[1],x6G08_MD_67_proj[3],pch=20)
text(x6G08_MD_67_proj[1],x6G08_MD_67_proj[3],pos=4,label="67",col="green")

x6G08_MD_68_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_68.txt")
x6G08_MD_68_proj<-project.pca(x6G08_MD_68_raw,pcadata)
points(x6G08_MD_68_proj[1],x6G08_MD_68_proj[3],pch=20)
text(x6G08_MD_68_proj[1],x6G08_MD_68_proj[3],pos=4,label="68",col="blue")

x6G08_MD_69_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_69.txt")
x6G08_MD_69_proj<-project.pca(x6G08_MD_69_raw,pcadata)
points(x6G08_MD_69_proj[1],x6G08_MD_69_proj[3],pch=20)
text(x6G08_MD_69_proj[1],x6G08_MD_69_proj[3],pos=4,label="69",col="blue")

x6G08_MD_70_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_70.txt")
x6G08_MD_70_proj<-project.pca(x6G08_MD_70_raw,pcadata)
points(x6G08_MD_70_proj[1],x6G08_MD_70_proj[3],pch=20)
text(x6G08_MD_70_proj[1],x6G08_MD_70_proj[3],pos=4,label="70",col="black")

x6G08_MD_71_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_71.txt")
x6G08_MD_71_proj<-project.pca(x6G08_MD_71_raw,pcadata)
points(x6G08_MD_71_proj[1],x6G08_MD_71_proj[3],pch=20)
text(x6G08_MD_71_proj[1],x6G08_MD_71_proj[3],pos=4,label="71",col="black")

x6G08_MD_72_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_72.txt")
x6G08_MD_72_proj<-project.pca(x6G08_MD_72_raw,pcadata)
points(x6G08_MD_72_proj[1],x6G08_MD_72_proj[3],pch=20)
text(x6G08_MD_72_proj[1],x6G08_MD_72_proj[3],pos=4,label="72",col="black")

x6G08_MD_73_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_73.txt")
x6G08_MD_73_proj<-project.pca(x6G08_MD_73_raw,pcadata)
points(x6G08_MD_73_proj[1],x6G08_MD_73_proj[3],pch=20)
text(x6G08_MD_73_proj[1],x6G08_MD_73_proj[3],pos=4,label="73",col="blue")

x6G08_MD_74_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_74.txt")
x6G08_MD_74_proj<-project.pca(x6G08_MD_74_raw,pcadata)
points(x6G08_MD_74_proj[1],x6G08_MD_74_proj[3],pch=20)
text(x6G08_MD_74_proj[1],x6G08_MD_74_proj[3],pos=4,label="74",col="blue")

x6G08_MD_75_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_75.txt")
x6G08_MD_75_proj<-project.pca(x6G08_MD_75_raw,pcadata)
points(x6G08_MD_75_proj[1],x6G08_MD_75_proj[3],pch=20)
text(x6G08_MD_75_proj[1],x6G08_MD_75_proj[3],pos=4,label="75",col="black")

x6G08_MD_76_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_76.txt")
x6G08_MD_76_proj<-project.pca(x6G08_MD_76_raw,pcadata)
points(x6G08_MD_76_proj[1],x6G08_MD_76_proj[3],pch=20)
text(x6G08_MD_76_proj[1],x6G08_MD_76_proj[3],pos=4,label="76",col="black")

x6G08_MD_77_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_77.txt")
x6G08_MD_77_proj<-project.pca(x6G08_MD_77_raw,pcadata)
points(x6G08_MD_77_proj[1],x6G08_MD_77_proj[3],pch=20)
text(x6G08_MD_77_proj[1],x6G08_MD_77_proj[3],pos=4,label="77",col="black")

x6G08_MD_78_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_78.txt")
x6G08_MD_78_proj<-project.pca(x6G08_MD_78_raw,pcadata)
points(x6G08_MD_78_proj[1],x6G08_MD_78_proj[3],pch=20)
text(x6G08_MD_78_proj[1],x6G08_MD_78_proj[3],pos=4,label="78",col="blue")

x6G08_MD_79_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_79.txt")
x6G08_MD_79_proj<-project.pca(x6G08_MD_79_raw,pcadata)
points(x6G08_MD_79_proj[1],x6G08_MD_79_proj[3],pch=20)
text(x6G08_MD_79_proj[1],x6G08_MD_79_proj[3],pos=4,label="79",col="green")

x6G08_MD_80_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_80.txt")
x6G08_MD_80_proj<-project.pca(x6G08_MD_80_raw,pcadata)
points(x6G08_MD_80_proj[1],x6G08_MD_80_proj[3],pch=20)
text(x6G08_MD_80_proj[1],x6G08_MD_80_proj[3],pos=4,label="80",col="blue")

x6G08_MD_81_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_81.txt")
x6G08_MD_81_proj<-project.pca(x6G08_MD_81_raw,pcadata)
points(x6G08_MD_81_proj[1],x6G08_MD_81_proj[3],pch=20)
text(x6G08_MD_81_proj[1],x6G08_MD_81_proj[3],pos=4,label="81",col="black")

x6G08_MD_82_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_82.txt")
x6G08_MD_82_proj<-project.pca(x6G08_MD_82_raw,pcadata)
points(x6G08_MD_82_proj[1],x6G08_MD_82_proj[3],pch=20)
text(x6G08_MD_82_proj[1],x6G08_MD_82_proj[3],pos=4,label="82",col="black")

x6G08_MD_83_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_83.txt")
x6G08_MD_83_proj<-project.pca(x6G08_MD_83_raw,pcadata)
points(x6G08_MD_83_proj[1],x6G08_MD_83_proj[3],pch=20)
text(x6G08_MD_83_proj[1],x6G08_MD_83_proj[3],pos=4,label="83",col="green")

x6G08_MD_84_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_84.txt")
x6G08_MD_84_proj<-project.pca(x6G08_MD_84_raw,pcadata)
points(x6G08_MD_84_proj[1],x6G08_MD_84_proj[3],pch=20)
text(x6G08_MD_84_proj[1],x6G08_MD_84_proj[3],pos=4,label="84",col="blue")

x6G08_MD_85_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_85.txt")
x6G08_MD_85_proj<-project.pca(x6G08_MD_85_raw,pcadata)
points(x6G08_MD_85_proj[1],x6G08_MD_85_proj[3],pch=20)
text(x6G08_MD_85_proj[1],x6G08_MD_85_proj[3],pos=4,label="85",col="black")

x6G08_MD_86_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_86.txt")
x6G08_MD_86_proj<-project.pca(x6G08_MD_86_raw,pcadata)
points(x6G08_MD_86_proj[1],x6G08_MD_86_proj[3],pch=20)
text(x6G08_MD_86_proj[1],x6G08_MD_86_proj[3],pos=4,label="86",col="blue")

x6G08_MD_87_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_87.txt")
x6G08_MD_87_proj<-project.pca(x6G08_MD_87_raw,pcadata)
points(x6G08_MD_87_proj[1],x6G08_MD_87_proj[3],pch=20)
text(x6G08_MD_87_proj[1],x6G08_MD_87_proj[3],pos=4,label="87",col="blue")

x6G08_MD_88_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_88.txt")
x6G08_MD_88_proj<-project.pca(x6G08_MD_88_raw,pcadata)
points(x6G08_MD_88_proj[1],x6G08_MD_88_proj[3],pch=20)
text(x6G08_MD_88_proj[1],x6G08_MD_88_proj[3],pos=4,label="88",col="green")

x6G08_MD_89_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_89.txt")
x6G08_MD_89_proj<-project.pca(x6G08_MD_89_raw,pcadata)
points(x6G08_MD_89_proj[1],x6G08_MD_89_proj[3],pch=20)
text(x6G08_MD_89_proj[1],x6G08_MD_89_proj[3],pos=4,label="89",col="green")

x6G08_MD_90_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_90.txt")
x6G08_MD_90_proj<-project.pca(x6G08_MD_90_raw,pcadata)
points(x6G08_MD_90_proj[1],x6G08_MD_90_proj[3],pch=20)
text(x6G08_MD_90_proj[1],x6G08_MD_90_proj[3],pos=4,label="90",col="green")

x6G08_MD_91_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_91.txt")
x6G08_MD_91_proj<-project.pca(x6G08_MD_91_raw,pcadata)
points(x6G08_MD_91_proj[1],x6G08_MD_91_proj[3],pch=20)
text(x6G08_MD_91_proj[1],x6G08_MD_91_proj[3],pos=4,label="91",col="green")

x6G08_MD_92_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_92.txt")
x6G08_MD_92_proj<-project.pca(x6G08_MD_92_raw,pcadata)
points(x6G08_MD_92_proj[1],x6G08_MD_92_proj[3],pch=20)
text(x6G08_MD_92_proj[1],x6G08_MD_92_proj[3],pos=4,label="92",col="green")

x6G08_MD_93_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_93.txt")
x6G08_MD_93_proj<-project.pca(x6G08_MD_93_raw,pcadata)
points(x6G08_MD_93_proj[1],x6G08_MD_93_proj[3],pch=20)
text(x6G08_MD_93_proj[1],x6G08_MD_93_proj[3],pos=4,label="93",col="blue")

x6G08_MD_94_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_94.txt")
x6G08_MD_94_proj<-project.pca(x6G08_MD_94_raw,pcadata)
points(x6G08_MD_94_proj[1],x6G08_MD_94_proj[3],pch=20)
text(x6G08_MD_94_proj[1],x6G08_MD_94_proj[3],pos=4,label="94",col="black")

x6G08_MD_95_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_95.txt")
x6G08_MD_95_proj<-project.pca(x6G08_MD_95_raw,pcadata)
points(x6G08_MD_95_proj[1],x6G08_MD_95_proj[3],pch=20)
text(x6G08_MD_95_proj[1],x6G08_MD_95_proj[3],pos=4,label="95",col="green")

x6G08_MD_96_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_96.txt")
x6G08_MD_96_proj<-project.pca(x6G08_MD_96_raw,pcadata)
points(x6G08_MD_96_proj[1],x6G08_MD_96_proj[3],pch=20)
text(x6G08_MD_96_proj[1],x6G08_MD_96_proj[3],pos=4,label="96",col="green")

x6G08_MD_97_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_97.txt")
x6G08_MD_97_proj<-project.pca(x6G08_MD_97_raw,pcadata)
points(x6G08_MD_97_proj[1],x6G08_MD_97_proj[3],pch=20)
text(x6G08_MD_97_proj[1],x6G08_MD_97_proj[3],pos=4,label="97",col="green")

x6G08_MD_98_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_98.txt")
x6G08_MD_98_proj<-project.pca(x6G08_MD_98_raw,pcadata)
points(x6G08_MD_98_proj[1],x6G08_MD_98_proj[3],pch=20)
text(x6G08_MD_98_proj[1],x6G08_MD_98_proj[3],pos=4,label="98",col="black")

x6G08_MD_99_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_99.txt")
x6G08_MD_99_proj<-project.pca(x6G08_MD_99_raw,pcadata)
points(x6G08_MD_99_proj[1],x6G08_MD_99_proj[3],pch=20)
text(x6G08_MD_99_proj[1],x6G08_MD_99_proj[3],pos=4,label="99",col="blue")

x6G08_MD_100_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_100.txt")
x6G08_MD_100_proj<-project.pca(x6G08_MD_100_raw,pcadata)
points(x6G08_MD_100_proj[1],x6G08_MD_100_proj[3],pch=20)
text(x6G08_MD_100_proj[1],x6G08_MD_100_proj[3],pos=4,label="100",col="black")

x6G08_MD_101_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_101.txt")
x6G08_MD_101_proj<-project.pca(x6G08_MD_101_raw,pcadata)
points(x6G08_MD_101_proj[1],x6G08_MD_101_proj[3],pch=20)
text(x6G08_MD_101_proj[1],x6G08_MD_101_proj[3],pos=4,label="101",col="black")

x6G08_MD_102_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_102.txt")
x6G08_MD_102_proj<-project.pca(x6G08_MD_102_raw,pcadata)
points(x6G08_MD_102_proj[1],x6G08_MD_102_proj[3],pch=20)
text(x6G08_MD_102_proj[1],x6G08_MD_102_proj[3],pos=4,label="102",col="green")

x6G08_MD_103_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_103.txt")
x6G08_MD_103_proj<-project.pca(x6G08_MD_103_raw,pcadata)
points(x6G08_MD_103_proj[1],x6G08_MD_103_proj[3],pch=20)
text(x6G08_MD_103_proj[1],x6G08_MD_103_proj[3],pos=4,label="103",col="blue")

x6G08_MD_104_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_104.txt")
x6G08_MD_104_proj<-project.pca(x6G08_MD_104_raw,pcadata)
points(x6G08_MD_104_proj[1],x6G08_MD_104_proj[3],pch=20)
text(x6G08_MD_104_proj[1],x6G08_MD_104_proj[3],pos=4,label="104",col="blue")

x6G08_MD_105_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_105.txt")
x6G08_MD_105_proj<-project.pca(x6G08_MD_105_raw,pcadata)
points(x6G08_MD_105_proj[1],x6G08_MD_105_proj[3],pch=20)
text(x6G08_MD_105_proj[1],x6G08_MD_105_proj[3],pos=4,label="105",col="blue")

x6G08_MD_106_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_106.txt")
x6G08_MD_106_proj<-project.pca(x6G08_MD_106_raw,pcadata)
points(x6G08_MD_106_proj[1],x6G08_MD_106_proj[3],pch=20)
text(x6G08_MD_106_proj[1],x6G08_MD_106_proj[3],pos=4,label="106",col="black")

x6G08_MD_107_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_107.txt")
x6G08_MD_107_proj<-project.pca(x6G08_MD_107_raw,pcadata)
points(x6G08_MD_107_proj[1],x6G08_MD_107_proj[3],pch=20)
text(x6G08_MD_107_proj[1],x6G08_MD_107_proj[3],pos=4,label="107",col="green")

x6G08_MD_108_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_108.txt")
x6G08_MD_108_proj<-project.pca(x6G08_MD_108_raw,pcadata)
points(x6G08_MD_108_proj[1],x6G08_MD_108_proj[3],pch=20)
text(x6G08_MD_108_proj[1],x6G08_MD_108_proj[3],pos=4,label="108",col="black")

x6G08_MD_109_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_109.txt")
x6G08_MD_109_proj<-project.pca(x6G08_MD_109_raw,pcadata)
points(x6G08_MD_109_proj[1],x6G08_MD_109_proj[3],pch=20)
text(x6G08_MD_109_proj[1],x6G08_MD_109_proj[3],pos=4,label="109",col="black")

x6G08_MD_110_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_110.txt")
x6G08_MD_110_proj<-project.pca(x6G08_MD_110_raw,pcadata)
points(x6G08_MD_110_proj[1],x6G08_MD_110_proj[3],pch=20)
text(x6G08_MD_110_proj[1],x6G08_MD_110_proj[3],pos=4,label="110",col="black")

x6G08_MD_111_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_111.txt")
x6G08_MD_111_proj<-project.pca(x6G08_MD_111_raw,pcadata)
points(x6G08_MD_111_proj[1],x6G08_MD_111_proj[3],pch=20)
text(x6G08_MD_111_proj[1],x6G08_MD_111_proj[3],pos=4,label="111",col="black")

x6G08_MD_112_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_112.txt")
x6G08_MD_112_proj<-project.pca(x6G08_MD_112_raw,pcadata)
points(x6G08_MD_112_proj[1],x6G08_MD_112_proj[3],pch=20)
text(x6G08_MD_112_proj[1],x6G08_MD_112_proj[3],pos=4,label="112",col="blue")

x6G08_MD_113_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_113.txt")
x6G08_MD_113_proj<-project.pca(x6G08_MD_113_raw,pcadata)
points(x6G08_MD_113_proj[1],x6G08_MD_113_proj[3],pch=20)
text(x6G08_MD_113_proj[1],x6G08_MD_113_proj[3],pos=4,label="113",col="blue")

x6G08_MD_114_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_114.txt")
x6G08_MD_114_proj<-project.pca(x6G08_MD_114_raw,pcadata)
points(x6G08_MD_114_proj[1],x6G08_MD_114_proj[3],pch=20)
text(x6G08_MD_114_proj[1],x6G08_MD_114_proj[3],pos=4,label="114",col="blue")

x6G08_MD_115_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_115.txt")
x6G08_MD_115_proj<-project.pca(x6G08_MD_115_raw,pcadata)
points(x6G08_MD_115_proj[1],x6G08_MD_115_proj[3],pch=20)
text(x6G08_MD_115_proj[1],x6G08_MD_115_proj[3],pos=4,label="115",col="black")

x6G08_MD_116_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_116.txt")
x6G08_MD_116_proj<-project.pca(x6G08_MD_116_raw,pcadata)
points(x6G08_MD_116_proj[1],x6G08_MD_116_proj[3],pch=20)
text(x6G08_MD_116_proj[1],x6G08_MD_116_proj[3],pos=4,label="116",col="black")

x6G08_MD_117_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_117.txt")
x6G08_MD_117_proj<-project.pca(x6G08_MD_117_raw,pcadata)
points(x6G08_MD_117_proj[1],x6G08_MD_117_proj[3],pch=20)
text(x6G08_MD_117_proj[1],x6G08_MD_117_proj[3],pos=4,label="117",col="blue")

x6G08_MD_118_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_118.txt")
x6G08_MD_118_proj<-project.pca(x6G08_MD_118_raw,pcadata)
points(x6G08_MD_118_proj[1],x6G08_MD_118_proj[3],pch=20)
text(x6G08_MD_118_proj[1],x6G08_MD_118_proj[3],pos=4,label="118",col="green")

x6G08_MD_119_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_119.txt")
x6G08_MD_119_proj<-project.pca(x6G08_MD_119_raw,pcadata)
points(x6G08_MD_119_proj[1],x6G08_MD_119_proj[3],pch=20)
text(x6G08_MD_119_proj[1],x6G08_MD_119_proj[3],pos=4,label="119",col="blue")

x6G08_MD_120_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_120.txt")
x6G08_MD_120_proj<-project.pca(x6G08_MD_120_raw,pcadata)
points(x6G08_MD_120_proj[1],x6G08_MD_120_proj[3],pch=20)
text(x6G08_MD_120_proj[1],x6G08_MD_120_proj[3],pos=4,label="120",col="black")

x6G08_MD_121_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_121.txt")
x6G08_MD_121_proj<-project.pca(x6G08_MD_121_raw,pcadata)
points(x6G08_MD_121_proj[1],x6G08_MD_121_proj[3],pch=20)
text(x6G08_MD_121_proj[1],x6G08_MD_121_proj[3],pos=4,label="121",col="black")

x6G08_MD_122_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_122.txt")
x6G08_MD_122_proj<-project.pca(x6G08_MD_122_raw,pcadata)
points(x6G08_MD_122_proj[1],x6G08_MD_122_proj[3],pch=20)
text(x6G08_MD_122_proj[1],x6G08_MD_122_proj[3],pos=4,label="122",col="black")

x6G08_MD_123_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_123.txt")
x6G08_MD_123_proj<-project.pca(x6G08_MD_123_raw,pcadata)
points(x6G08_MD_123_proj[1],x6G08_MD_123_proj[3],pch=20)
text(x6G08_MD_123_proj[1],x6G08_MD_123_proj[3],pos=4,label="123",col="black")

x6G08_MD_124_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_124.txt")
x6G08_MD_124_proj<-project.pca(x6G08_MD_124_raw,pcadata)
points(x6G08_MD_124_proj[1],x6G08_MD_124_proj[3],pch=20)
text(x6G08_MD_124_proj[1],x6G08_MD_124_proj[3],pos=4,label="124",col="black")

x6G08_MD_125_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_125.txt")
x6G08_MD_125_proj<-project.pca(x6G08_MD_125_raw,pcadata)
points(x6G08_MD_125_proj[1],x6G08_MD_125_proj[3],pch=20)
text(x6G08_MD_125_proj[1],x6G08_MD_125_proj[3],pos=4,label="125",col="blue")

x6G08_MD_126_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_126.txt")
x6G08_MD_126_proj<-project.pca(x6G08_MD_126_raw,pcadata)
points(x6G08_MD_126_proj[1],x6G08_MD_126_proj[3],pch=20)
text(x6G08_MD_126_proj[1],x6G08_MD_126_proj[3],pos=4,label="126",col="black")

x6G08_MD_127_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_127.txt")
x6G08_MD_127_proj<-project.pca(x6G08_MD_127_raw,pcadata)
points(x6G08_MD_127_proj[1],x6G08_MD_127_proj[3],pch=20)
text(x6G08_MD_127_proj[1],x6G08_MD_127_proj[3],pos=4,label="127",col="black")

x6G08_MD_128_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_128.txt")
x6G08_MD_128_proj<-project.pca(x6G08_MD_128_raw,pcadata)
points(x6G08_MD_128_proj[1],x6G08_MD_128_proj[3],pch=20)
text(x6G08_MD_128_proj[1],x6G08_MD_128_proj[3],pos=4,label="128",col="black")

x6G08_MD_129_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_129.txt")
x6G08_MD_129_proj<-project.pca(x6G08_MD_129_raw,pcadata)
points(x6G08_MD_129_proj[1],x6G08_MD_129_proj[3],pch=20)
text(x6G08_MD_129_proj[1],x6G08_MD_129_proj[3],pos=4,label="129",col="black")

x6G08_MD_130_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_130.txt")
x6G08_MD_130_proj<-project.pca(x6G08_MD_130_raw,pcadata)
points(x6G08_MD_130_proj[1],x6G08_MD_130_proj[3],pch=20)
text(x6G08_MD_130_proj[1],x6G08_MD_130_proj[3],pos=4,label="130",col="black")

x6G08_MD_131_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_131.txt")
x6G08_MD_131_proj<-project.pca(x6G08_MD_131_raw,pcadata)
points(x6G08_MD_131_proj[1],x6G08_MD_131_proj[3],pch=20)
text(x6G08_MD_131_proj[1],x6G08_MD_131_proj[3],pos=4,label="131",col="black")

x6G08_MD_132_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_132.txt")
x6G08_MD_132_proj<-project.pca(x6G08_MD_132_raw,pcadata)
points(x6G08_MD_132_proj[1],x6G08_MD_132_proj[3],pch=20)
text(x6G08_MD_132_proj[1],x6G08_MD_132_proj[3],pos=4,label="132",col="black")

x6G08_MD_133_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_133.txt")
x6G08_MD_133_proj<-project.pca(x6G08_MD_133_raw,pcadata)
points(x6G08_MD_133_proj[1],x6G08_MD_133_proj[3],pch=20)
text(x6G08_MD_133_proj[1],x6G08_MD_133_proj[3],pos=4,label="133",col="black")

x6G08_MD_134_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_134.txt")
x6G08_MD_134_proj<-project.pca(x6G08_MD_134_raw,pcadata)
points(x6G08_MD_134_proj[1],x6G08_MD_134_proj[3],pch=20)
text(x6G08_MD_134_proj[1],x6G08_MD_134_proj[3],pos=4,label="134",col="black")

x6G08_MD_135_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_135.txt")
x6G08_MD_135_proj<-project.pca(x6G08_MD_135_raw,pcadata)
points(x6G08_MD_135_proj[1],x6G08_MD_135_proj[3],pch=20)
text(x6G08_MD_135_proj[1],x6G08_MD_135_proj[3],pos=4,label="135",col="black")

x6G08_MD_136_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_136.txt")
x6G08_MD_136_proj<-project.pca(x6G08_MD_136_raw,pcadata)
points(x6G08_MD_136_proj[1],x6G08_MD_136_proj[3],pch=20)
text(x6G08_MD_136_proj[1],x6G08_MD_136_proj[3],pos=4,label="136",col="black")

x6G08_MD_137_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_137.txt")
x6G08_MD_137_proj<-project.pca(x6G08_MD_137_raw,pcadata)
points(x6G08_MD_137_proj[1],x6G08_MD_137_proj[3],pch=20)
text(x6G08_MD_137_proj[1],x6G08_MD_137_proj[3],pos=4,label="137",col="black")

x6G08_MD_138_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_138.txt")
x6G08_MD_138_proj<-project.pca(x6G08_MD_138_raw,pcadata)
points(x6G08_MD_138_proj[1],x6G08_MD_138_proj[3],pch=20)
text(x6G08_MD_138_proj[1],x6G08_MD_138_proj[3],pos=4,label="138",col="black")

x6G08_MD_139_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_139.txt")
x6G08_MD_139_proj<-project.pca(x6G08_MD_139_raw,pcadata)
points(x6G08_MD_139_proj[1],x6G08_MD_139_proj[3],pch=20)
text(x6G08_MD_139_proj[1],x6G08_MD_139_proj[3],pos=4,label="139",col="black")

x6G08_MD_140_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_140.txt")
x6G08_MD_140_proj<-project.pca(x6G08_MD_140_raw,pcadata)
points(x6G08_MD_140_proj[1],x6G08_MD_140_proj[3],pch=20)
text(x6G08_MD_140_proj[1],x6G08_MD_140_proj[3],pos=4,label="140",col="black")

x6G08_MD_141_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_141.txt")
x6G08_MD_141_proj<-project.pca(x6G08_MD_141_raw,pcadata)
points(x6G08_MD_141_proj[1],x6G08_MD_141_proj[3],pch=20)
text(x6G08_MD_141_proj[1],x6G08_MD_141_proj[3],pos=4,label="141",col="green")

x6G08_MD_142_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_142.txt")
x6G08_MD_142_proj<-project.pca(x6G08_MD_142_raw,pcadata)
points(x6G08_MD_142_proj[1],x6G08_MD_142_proj[3],pch=20)
text(x6G08_MD_142_proj[1],x6G08_MD_142_proj[3],pos=4,label="142",col="black")

x6G08_MD_143_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_143.txt")
x6G08_MD_143_proj<-project.pca(x6G08_MD_143_raw,pcadata)
points(x6G08_MD_143_proj[1],x6G08_MD_143_proj[3],pch=20)
text(x6G08_MD_143_proj[1],x6G08_MD_143_proj[3],pos=4,label="143",col="black")

x6G08_MD_144_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_144.txt")
x6G08_MD_144_proj<-project.pca(x6G08_MD_144_raw,pcadata)
points(x6G08_MD_144_proj[1],x6G08_MD_144_proj[3],pch=20)
text(x6G08_MD_144_proj[1],x6G08_MD_144_proj[3],pos=4,label="144",col="blue")

x6G08_MD_145_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_145.txt")
x6G08_MD_145_proj<-project.pca(x6G08_MD_145_raw,pcadata)
points(x6G08_MD_145_proj[1],x6G08_MD_145_proj[3],pch=20)
text(x6G08_MD_145_proj[1],x6G08_MD_145_proj[3],pos=4,label="145",col="black")

x6G08_MD_146_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_146.txt")
x6G08_MD_146_proj<-project.pca(x6G08_MD_146_raw,pcadata)
points(x6G08_MD_146_proj[1],x6G08_MD_146_proj[3],pch=20)
text(x6G08_MD_146_proj[1],x6G08_MD_146_proj[3],pos=4,label="146",col="black")

x6G08_MD_147_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_147.txt")
x6G08_MD_147_proj<-project.pca(x6G08_MD_147_raw,pcadata)
points(x6G08_MD_147_proj[1],x6G08_MD_147_proj[3],pch=20)
text(x6G08_MD_147_proj[1],x6G08_MD_147_proj[3],pos=4,label="147",col="green")

x6G08_MD_148_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_148.txt")
x6G08_MD_148_proj<-project.pca(x6G08_MD_148_raw,pcadata)
points(x6G08_MD_148_proj[1],x6G08_MD_148_proj[3],pch=20)
text(x6G08_MD_148_proj[1],x6G08_MD_148_proj[3],pos=4,label="148",col="green")

x6G08_MD_149_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_149.txt")
x6G08_MD_149_proj<-project.pca(x6G08_MD_149_raw,pcadata)
points(x6G08_MD_149_proj[1],x6G08_MD_149_proj[3],pch=20)
text(x6G08_MD_149_proj[1],x6G08_MD_149_proj[3],pos=4,label="149",col="blue")

x6G08_MD_150_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_150.txt")
x6G08_MD_150_proj<-project.pca(x6G08_MD_150_raw,pcadata)
points(x6G08_MD_150_proj[1],x6G08_MD_150_proj[3],pch=20)
text(x6G08_MD_150_proj[1],x6G08_MD_150_proj[3],pos=4,label="150",col="red")

x6G08_MD_151_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_151.txt")
x6G08_MD_151_proj<-project.pca(x6G08_MD_151_raw,pcadata)
points(x6G08_MD_151_proj[1],x6G08_MD_151_proj[3],pch=20)
text(x6G08_MD_151_proj[1],x6G08_MD_151_proj[3],pos=4,label="151",col="black")

x6G08_MD_152_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_152.txt")
x6G08_MD_152_proj<-project.pca(x6G08_MD_152_raw,pcadata)
points(x6G08_MD_152_proj[1],x6G08_MD_152_proj[3],pch=20)
text(x6G08_MD_152_proj[1],x6G08_MD_152_proj[3],pos=4,label="152",col="black")

x6G08_MD_153_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_153.txt")
x6G08_MD_153_proj<-project.pca(x6G08_MD_153_raw,pcadata)
points(x6G08_MD_153_proj[1],x6G08_MD_153_proj[3],pch=20)
text(x6G08_MD_153_proj[1],x6G08_MD_153_proj[3],pos=4,label="153",col="black")

x6G08_MD_154_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_154.txt")
x6G08_MD_154_proj<-project.pca(x6G08_MD_154_raw,pcadata)
points(x6G08_MD_154_proj[1],x6G08_MD_154_proj[3],pch=20)
text(x6G08_MD_154_proj[1],x6G08_MD_154_proj[3],pos=4,label="154",col="black")

x6G08_MD_155_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_155.txt")
x6G08_MD_155_proj<-project.pca(x6G08_MD_155_raw,pcadata)
points(x6G08_MD_155_proj[1],x6G08_MD_155_proj[3],pch=20)
text(x6G08_MD_155_proj[1],x6G08_MD_155_proj[3],pos=4,label="155",col="black")

x6G08_MD_156_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_156.txt")
x6G08_MD_156_proj<-project.pca(x6G08_MD_156_raw,pcadata)
points(x6G08_MD_156_proj[1],x6G08_MD_156_proj[3],pch=20)
text(x6G08_MD_156_proj[1],x6G08_MD_156_proj[3],pos=4,label="156",col="black")

x6G08_MD_157_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_157.txt")
x6G08_MD_157_proj<-project.pca(x6G08_MD_157_raw,pcadata)
points(x6G08_MD_157_proj[1],x6G08_MD_157_proj[3],pch=20)
text(x6G08_MD_157_proj[1],x6G08_MD_157_proj[3],pos=4,label="157",col="black")

x6G08_MD_158_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_158.txt")
x6G08_MD_158_proj<-project.pca(x6G08_MD_158_raw,pcadata)
points(x6G08_MD_158_proj[1],x6G08_MD_158_proj[3],pch=20)
text(x6G08_MD_158_proj[1],x6G08_MD_158_proj[3],pos=4,label="158",col="black")

x6G08_MD_159_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_159.txt")
x6G08_MD_159_proj<-project.pca(x6G08_MD_159_raw,pcadata)
points(x6G08_MD_159_proj[1],x6G08_MD_159_proj[3],pch=20)
text(x6G08_MD_159_proj[1],x6G08_MD_159_proj[3],pos=4,label="159",col="black")

x6G08_MD_160_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_160.txt")
x6G08_MD_160_proj<-project.pca(x6G08_MD_160_raw,pcadata)
points(x6G08_MD_160_proj[1],x6G08_MD_160_proj[3],pch=20)
text(x6G08_MD_160_proj[1],x6G08_MD_160_proj[3],pos=4,label="160",col="black")

x6G08_MD_161_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_161.txt")
x6G08_MD_161_proj<-project.pca(x6G08_MD_161_raw,pcadata)
points(x6G08_MD_161_proj[1],x6G08_MD_161_proj[3],pch=20)
text(x6G08_MD_161_proj[1],x6G08_MD_161_proj[3],pos=4,label="161",col="black")

x6G08_MD_162_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_162.txt")
x6G08_MD_162_proj<-project.pca(x6G08_MD_162_raw,pcadata)
points(x6G08_MD_162_proj[1],x6G08_MD_162_proj[3],pch=20)
text(x6G08_MD_162_proj[1],x6G08_MD_162_proj[3],pos=4,label="162",col="black")

x6G08_MD_163_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_163.txt")
x6G08_MD_163_proj<-project.pca(x6G08_MD_163_raw,pcadata)
points(x6G08_MD_163_proj[1],x6G08_MD_163_proj[3],pch=20)
text(x6G08_MD_163_proj[1],x6G08_MD_163_proj[3],pos=4,label="163",col="black")

x6G08_MD_164_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_164.txt")
x6G08_MD_164_proj<-project.pca(x6G08_MD_164_raw,pcadata)
points(x6G08_MD_164_proj[1],x6G08_MD_164_proj[3],pch=20)
text(x6G08_MD_164_proj[1],x6G08_MD_164_proj[3],pos=4,label="164",col="black")

x6G08_MD_165_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_165.txt")
x6G08_MD_165_proj<-project.pca(x6G08_MD_165_raw,pcadata)
points(x6G08_MD_165_proj[1],x6G08_MD_165_proj[3],pch=20)
text(x6G08_MD_165_proj[1],x6G08_MD_165_proj[3],pos=4,label="165",col="green")

x6G08_MD_166_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_166.txt")
x6G08_MD_166_proj<-project.pca(x6G08_MD_166_raw,pcadata)
points(x6G08_MD_166_proj[1],x6G08_MD_166_proj[3],pch=20)
text(x6G08_MD_166_proj[1],x6G08_MD_166_proj[3],pos=4,label="166",col="black")

x6G08_MD_167_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_167.txt")
x6G08_MD_167_proj<-project.pca(x6G08_MD_167_raw,pcadata)
points(x6G08_MD_167_proj[1],x6G08_MD_167_proj[3],pch=20)
text(x6G08_MD_167_proj[1],x6G08_MD_167_proj[3],pos=4,label="167",col="black")

x6G08_MD_168_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_168.txt")
x6G08_MD_168_proj<-project.pca(x6G08_MD_168_raw,pcadata)
points(x6G08_MD_168_proj[1],x6G08_MD_168_proj[3],pch=20)
text(x6G08_MD_168_proj[1],x6G08_MD_168_proj[3],pos=4,label="168",col="black")

x6G08_MD_169_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_169.txt")
x6G08_MD_169_proj<-project.pca(x6G08_MD_169_raw,pcadata)
points(x6G08_MD_169_proj[1],x6G08_MD_169_proj[3],pch=20)
text(x6G08_MD_169_proj[1],x6G08_MD_169_proj[3],pos=4,label="169",col="black")

x6G08_MD_170_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_170.txt")
x6G08_MD_170_proj<-project.pca(x6G08_MD_170_raw,pcadata)
points(x6G08_MD_170_proj[1],x6G08_MD_170_proj[3],pch=20)
text(x6G08_MD_170_proj[1],x6G08_MD_170_proj[3],pos=4,label="170",col="red")

x6G08_MD_171_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_171.txt")
x6G08_MD_171_proj<-project.pca(x6G08_MD_171_raw,pcadata)
points(x6G08_MD_171_proj[1],x6G08_MD_171_proj[3],pch=20)
text(x6G08_MD_171_proj[1],x6G08_MD_171_proj[3],pos=4,label="171",col="red")

x6G08_MD_172_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_172.txt")
x6G08_MD_172_proj<-project.pca(x6G08_MD_172_raw,pcadata)
points(x6G08_MD_172_proj[1],x6G08_MD_172_proj[3],pch=20)
text(x6G08_MD_172_proj[1],x6G08_MD_172_proj[3],pos=4,label="172",col="blue")

x6G08_MD_173_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_173.txt")
x6G08_MD_173_proj<-project.pca(x6G08_MD_173_raw,pcadata)
points(x6G08_MD_173_proj[1],x6G08_MD_173_proj[3],pch=20)
text(x6G08_MD_173_proj[1],x6G08_MD_173_proj[3],pos=4,label="173",col="black")

x6G08_MD_174_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_174.txt")
x6G08_MD_174_proj<-project.pca(x6G08_MD_174_raw,pcadata)
points(x6G08_MD_174_proj[1],x6G08_MD_174_proj[3],pch=20)
text(x6G08_MD_174_proj[1],x6G08_MD_174_proj[3],pos=4,label="174",col="blue")

x6G08_MD_175_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_175.txt")
x6G08_MD_175_proj<-project.pca(x6G08_MD_175_raw,pcadata)
points(x6G08_MD_175_proj[1],x6G08_MD_175_proj[3],pch=20)
text(x6G08_MD_175_proj[1],x6G08_MD_175_proj[3],pos=4,label="175",col="blue")

x6G08_MD_176_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_176.txt")
x6G08_MD_176_proj<-project.pca(x6G08_MD_176_raw,pcadata)
points(x6G08_MD_176_proj[1],x6G08_MD_176_proj[3],pch=20)
text(x6G08_MD_176_proj[1],x6G08_MD_176_proj[3],pos=4,label="176",col="blue")

x6G08_MD_177_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_177.txt")
x6G08_MD_177_proj<-project.pca(x6G08_MD_177_raw,pcadata)
points(x6G08_MD_177_proj[1],x6G08_MD_177_proj[3],pch=20)
text(x6G08_MD_177_proj[1],x6G08_MD_177_proj[3],pos=4,label="177",col="blue")

x6G08_MD_178_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_178.txt")
x6G08_MD_178_proj<-project.pca(x6G08_MD_178_raw,pcadata)
points(x6G08_MD_178_proj[1],x6G08_MD_178_proj[3],pch=20)
text(x6G08_MD_178_proj[1],x6G08_MD_178_proj[3],pos=4,label="178",col="blue")

x6G08_MD_179_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_179.txt")
x6G08_MD_179_proj<-project.pca(x6G08_MD_179_raw,pcadata)
points(x6G08_MD_179_proj[1],x6G08_MD_179_proj[3],pch=20)
text(x6G08_MD_179_proj[1],x6G08_MD_179_proj[3],pos=4,label="179",col="blue")

x6G08_MD_180_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_180.txt")
x6G08_MD_180_proj<-project.pca(x6G08_MD_180_raw,pcadata)
points(x6G08_MD_180_proj[1],x6G08_MD_180_proj[3],pch=20)
text(x6G08_MD_180_proj[1],x6G08_MD_180_proj[3],pos=4,label="180",col="black")

x6G08_MD_181_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_181.txt")
x6G08_MD_181_proj<-project.pca(x6G08_MD_181_raw,pcadata)
points(x6G08_MD_181_proj[1],x6G08_MD_181_proj[3],pch=20)
text(x6G08_MD_181_proj[1],x6G08_MD_181_proj[3],pos=4,label="181",col="black")

x6G08_MD_182_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_182.txt")
x6G08_MD_182_proj<-project.pca(x6G08_MD_182_raw,pcadata)
points(x6G08_MD_182_proj[1],x6G08_MD_182_proj[3],pch=20)
text(x6G08_MD_182_proj[1],x6G08_MD_182_proj[3],pos=4,label="182",col="black")

x6G08_MD_183_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_183.txt")
x6G08_MD_183_proj<-project.pca(x6G08_MD_183_raw,pcadata)
points(x6G08_MD_183_proj[1],x6G08_MD_183_proj[3],pch=20)
text(x6G08_MD_183_proj[1],x6G08_MD_183_proj[3],pos=4,label="183",col="black")

x6G08_MD_184_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_184.txt")
x6G08_MD_184_proj<-project.pca(x6G08_MD_184_raw,pcadata)
points(x6G08_MD_184_proj[1],x6G08_MD_184_proj[3],pch=20)
text(x6G08_MD_184_proj[1],x6G08_MD_184_proj[3],pos=4,label="184",col="green")

x6G08_MD_185_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_185.txt")
x6G08_MD_185_proj<-project.pca(x6G08_MD_185_raw,pcadata)
points(x6G08_MD_185_proj[1],x6G08_MD_185_proj[3],pch=20)
text(x6G08_MD_185_proj[1],x6G08_MD_185_proj[3],pos=4,label="185",col="blue")

x6G08_MD_186_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_186.txt")
x6G08_MD_186_proj<-project.pca(x6G08_MD_186_raw,pcadata)
points(x6G08_MD_186_proj[1],x6G08_MD_186_proj[3],pch=20)
text(x6G08_MD_186_proj[1],x6G08_MD_186_proj[3],pos=4,label="186",col="red")

x6G08_MD_187_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_187.txt")
x6G08_MD_187_proj<-project.pca(x6G08_MD_187_raw,pcadata)
points(x6G08_MD_187_proj[1],x6G08_MD_187_proj[3],pch=20)
text(x6G08_MD_187_proj[1],x6G08_MD_187_proj[3],pos=4,label="187",col="blue")

x6G08_MD_188_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_188.txt")
x6G08_MD_188_proj<-project.pca(x6G08_MD_188_raw,pcadata)
points(x6G08_MD_188_proj[1],x6G08_MD_188_proj[3],pch=20)
text(x6G08_MD_188_proj[1],x6G08_MD_188_proj[3],pos=4,label="188",col="blue")

x6G08_MD_189_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_189.txt")
x6G08_MD_189_proj<-project.pca(x6G08_MD_189_raw,pcadata)
points(x6G08_MD_189_proj[1],x6G08_MD_189_proj[3],pch=20)
text(x6G08_MD_189_proj[1],x6G08_MD_189_proj[3],pos=4,label="189",col="blue")

x6G08_MD_190_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_190.txt")
x6G08_MD_190_proj<-project.pca(x6G08_MD_190_raw,pcadata)
points(x6G08_MD_190_proj[1],x6G08_MD_190_proj[3],pch=20)
text(x6G08_MD_190_proj[1],x6G08_MD_190_proj[3],pos=4,label="190",col="blue")

x6G08_MD_191_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_191.txt")
x6G08_MD_191_proj<-project.pca(x6G08_MD_191_raw,pcadata)
points(x6G08_MD_191_proj[1],x6G08_MD_191_proj[3],pch=20)
text(x6G08_MD_191_proj[1],x6G08_MD_191_proj[3],pos=4,label="191",col="black")

x6G08_MD_192_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_192.txt")
x6G08_MD_192_proj<-project.pca(x6G08_MD_192_raw,pcadata)
points(x6G08_MD_192_proj[1],x6G08_MD_192_proj[3],pch=20)
text(x6G08_MD_192_proj[1],x6G08_MD_192_proj[3],pos=4,label="192",col="black")

x6G08_MD_193_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_193.txt")
x6G08_MD_193_proj<-project.pca(x6G08_MD_193_raw,pcadata)
points(x6G08_MD_193_proj[1],x6G08_MD_193_proj[3],pch=20)
text(x6G08_MD_193_proj[1],x6G08_MD_193_proj[3],pos=4,label="193",col="blue")

x6G08_MD_194_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_194.txt")
x6G08_MD_194_proj<-project.pca(x6G08_MD_194_raw,pcadata)
points(x6G08_MD_194_proj[1],x6G08_MD_194_proj[3],pch=20)
text(x6G08_MD_194_proj[1],x6G08_MD_194_proj[3],pos=4,label="194",col="blue")

x6G08_MD_195_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_195.txt")
x6G08_MD_195_proj<-project.pca(x6G08_MD_195_raw,pcadata)
points(x6G08_MD_195_proj[1],x6G08_MD_195_proj[3],pch=20)
text(x6G08_MD_195_proj[1],x6G08_MD_195_proj[3],pos=4,label="195",col="blue")

x6G08_MD_196_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_196.txt")
x6G08_MD_196_proj<-project.pca(x6G08_MD_196_raw,pcadata)
points(x6G08_MD_196_proj[1],x6G08_MD_196_proj[3],pch=20)
text(x6G08_MD_196_proj[1],x6G08_MD_196_proj[3],pos=4,label="196",col="black")

x6G08_MD_197_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_197.txt")
x6G08_MD_197_proj<-project.pca(x6G08_MD_197_raw,pcadata)
points(x6G08_MD_197_proj[1],x6G08_MD_197_proj[3],pch=20)
text(x6G08_MD_197_proj[1],x6G08_MD_197_proj[3],pos=4,label="197",col="black")

x6G08_MD_198_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_198.txt")
x6G08_MD_198_proj<-project.pca(x6G08_MD_198_raw,pcadata)
points(x6G08_MD_198_proj[1],x6G08_MD_198_proj[3],pch=20)
text(x6G08_MD_198_proj[1],x6G08_MD_198_proj[3],pos=4,label="198",col="black")

x6G08_MD_199_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_199.txt")
x6G08_MD_199_proj<-project.pca(x6G08_MD_199_raw,pcadata)
points(x6G08_MD_199_proj[1],x6G08_MD_199_proj[3],pch=20)
text(x6G08_MD_199_proj[1],x6G08_MD_199_proj[3],pos=4,label="199",col="black")

x6G08_MD_200_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_200.txt")
x6G08_MD_200_proj<-project.pca(x6G08_MD_200_raw,pcadata)
points(x6G08_MD_200_proj[1],x6G08_MD_200_proj[3],pch=20)
text(x6G08_MD_200_proj[1],x6G08_MD_200_proj[3],pos=4,label="200",col="black")

x6G08_MD_201_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_201.txt")
x6G08_MD_201_proj<-project.pca(x6G08_MD_201_raw,pcadata)
points(x6G08_MD_201_proj[1],x6G08_MD_201_proj[3],pch=20)
text(x6G08_MD_201_proj[1],x6G08_MD_201_proj[3],pos=4,label="201",col="blue")

x6G08_MD_202_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_202.txt")
x6G08_MD_202_proj<-project.pca(x6G08_MD_202_raw,pcadata)
points(x6G08_MD_202_proj[1],x6G08_MD_202_proj[3],pch=20)
text(x6G08_MD_202_proj[1],x6G08_MD_202_proj[3],pos=4,label="202",col="black")

x6G08_MD_203_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_203.txt")
x6G08_MD_203_proj<-project.pca(x6G08_MD_203_raw,pcadata)
points(x6G08_MD_203_proj[1],x6G08_MD_203_proj[3],pch=20)
text(x6G08_MD_203_proj[1],x6G08_MD_203_proj[3],pos=4,label="203",col="blue")

x6G08_MD_204_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_204.txt")
x6G08_MD_204_proj<-project.pca(x6G08_MD_204_raw,pcadata)
points(x6G08_MD_204_proj[1],x6G08_MD_204_proj[3],pch=20)
text(x6G08_MD_204_proj[1],x6G08_MD_204_proj[3],pos=4,label="204",col="blue")

x6G08_MD_205_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_205.txt")
x6G08_MD_205_proj<-project.pca(x6G08_MD_205_raw,pcadata)
points(x6G08_MD_205_proj[1],x6G08_MD_205_proj[3],pch=20)
text(x6G08_MD_205_proj[1],x6G08_MD_205_proj[3],pos=4,label="205",col="blue")

x6G08_MD_206_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_206.txt")
x6G08_MD_206_proj<-project.pca(x6G08_MD_206_raw,pcadata)
points(x6G08_MD_206_proj[1],x6G08_MD_206_proj[3],pch=20)
text(x6G08_MD_206_proj[1],x6G08_MD_206_proj[3],pos=4,label="206",col="red")

x6G08_MD_207_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_207.txt")
x6G08_MD_207_proj<-project.pca(x6G08_MD_207_raw,pcadata)
points(x6G08_MD_207_proj[1],x6G08_MD_207_proj[3],pch=20)
text(x6G08_MD_207_proj[1],x6G08_MD_207_proj[3],pos=4,label="207",col="blue")

x6G08_MD_208_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_208.txt")
x6G08_MD_208_proj<-project.pca(x6G08_MD_208_raw,pcadata)
points(x6G08_MD_208_proj[1],x6G08_MD_208_proj[3],pch=20)
text(x6G08_MD_208_proj[1],x6G08_MD_208_proj[3],pos=4,label="208",col="blue")

x6G08_MD_209_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_209.txt")
x6G08_MD_209_proj<-project.pca(x6G08_MD_209_raw,pcadata)
points(x6G08_MD_209_proj[1],x6G08_MD_209_proj[3],pch=20)
text(x6G08_MD_209_proj[1],x6G08_MD_209_proj[3],pos=4,label="209",col="black")

x6G08_MD_210_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_210.txt")
x6G08_MD_210_proj<-project.pca(x6G08_MD_210_raw,pcadata)
points(x6G08_MD_210_proj[1],x6G08_MD_210_proj[3],pch=20)
text(x6G08_MD_210_proj[1],x6G08_MD_210_proj[3],pos=4,label="210",col="blue")

x6G08_MD_211_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_211.txt")
x6G08_MD_211_proj<-project.pca(x6G08_MD_211_raw,pcadata)
points(x6G08_MD_211_proj[1],x6G08_MD_211_proj[3],pch=20)
text(x6G08_MD_211_proj[1],x6G08_MD_211_proj[3],pos=4,label="211",col="blue")

x6G08_MD_212_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_212.txt")
x6G08_MD_212_proj<-project.pca(x6G08_MD_212_raw,pcadata)
points(x6G08_MD_212_proj[1],x6G08_MD_212_proj[3],pch=20)
text(x6G08_MD_212_proj[1],x6G08_MD_212_proj[3],pos=4,label="212",col="blue")

x6G08_MD_213_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_213.txt")
x6G08_MD_213_proj<-project.pca(x6G08_MD_213_raw,pcadata)
points(x6G08_MD_213_proj[1],x6G08_MD_213_proj[3],pch=20)
text(x6G08_MD_213_proj[1],x6G08_MD_213_proj[3],pos=4,label="213",col="blue")

x6G08_MD_214_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_214.txt")
x6G08_MD_214_proj<-project.pca(x6G08_MD_214_raw,pcadata)
points(x6G08_MD_214_proj[1],x6G08_MD_214_proj[3],pch=20)
text(x6G08_MD_214_proj[1],x6G08_MD_214_proj[3],pos=4,label="214",col="blue")

x6G08_MD_215_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_215.txt")
x6G08_MD_215_proj<-project.pca(x6G08_MD_215_raw,pcadata)
points(x6G08_MD_215_proj[1],x6G08_MD_215_proj[3],pch=20)
text(x6G08_MD_215_proj[1],x6G08_MD_215_proj[3],pos=4,label="215",col="blue")

x6G08_MD_216_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_216.txt")
x6G08_MD_216_proj<-project.pca(x6G08_MD_216_raw,pcadata)
points(x6G08_MD_216_proj[1],x6G08_MD_216_proj[3],pch=20)
text(x6G08_MD_216_proj[1],x6G08_MD_216_proj[3],pos=4,label="216",col="blue")

x6G08_MD_217_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_217.txt")
x6G08_MD_217_proj<-project.pca(x6G08_MD_217_raw,pcadata)
points(x6G08_MD_217_proj[1],x6G08_MD_217_proj[3],pch=20)
text(x6G08_MD_217_proj[1],x6G08_MD_217_proj[3],pos=4,label="217",col="blue")

x6G08_MD_218_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_218.txt")
x6G08_MD_218_proj<-project.pca(x6G08_MD_218_raw,pcadata)
points(x6G08_MD_218_proj[1],x6G08_MD_218_proj[3],pch=20)
text(x6G08_MD_218_proj[1],x6G08_MD_218_proj[3],pos=4,label="218",col="green")

x6G08_MD_219_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_219.txt")
x6G08_MD_219_proj<-project.pca(x6G08_MD_219_raw,pcadata)
points(x6G08_MD_219_proj[1],x6G08_MD_219_proj[3],pch=20)
text(x6G08_MD_219_proj[1],x6G08_MD_219_proj[3],pos=4,label="219",col="red")

x6G08_MD_220_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_220.txt")
x6G08_MD_220_proj<-project.pca(x6G08_MD_220_raw,pcadata)
points(x6G08_MD_220_proj[1],x6G08_MD_220_proj[3],pch=20)
text(x6G08_MD_220_proj[1],x6G08_MD_220_proj[3],pos=4,label="220",col="red")

x6G08_MD_221_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_221.txt")
x6G08_MD_221_proj<-project.pca(x6G08_MD_221_raw,pcadata)
points(x6G08_MD_221_proj[1],x6G08_MD_221_proj[3],pch=20)
text(x6G08_MD_221_proj[1],x6G08_MD_221_proj[3],pos=4,label="221",col="blue")

x6G08_MD_222_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_222.txt")
x6G08_MD_222_proj<-project.pca(x6G08_MD_222_raw,pcadata)
points(x6G08_MD_222_proj[1],x6G08_MD_222_proj[3],pch=20)
text(x6G08_MD_222_proj[1],x6G08_MD_222_proj[3],pos=4,label="222",col="green")

x6G08_MD_223_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_223.txt")
x6G08_MD_223_proj<-project.pca(x6G08_MD_223_raw,pcadata)
points(x6G08_MD_223_proj[1],x6G08_MD_223_proj[3],pch=20)
text(x6G08_MD_223_proj[1],x6G08_MD_223_proj[3],pos=4,label="223",col="red")

x6G08_MD_224_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_224.txt")
x6G08_MD_224_proj<-project.pca(x6G08_MD_224_raw,pcadata)
points(x6G08_MD_224_proj[1],x6G08_MD_224_proj[3],pch=20)
text(x6G08_MD_224_proj[1],x6G08_MD_224_proj[3],pos=4,label="224",col="red")

x6G08_MD_225_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_225.txt")
x6G08_MD_225_proj<-project.pca(x6G08_MD_225_raw,pcadata)
points(x6G08_MD_225_proj[1],x6G08_MD_225_proj[3],pch=20)
text(x6G08_MD_225_proj[1],x6G08_MD_225_proj[3],pos=4,label="225",col="red")

x6G08_MD_226_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_226.txt")
x6G08_MD_226_proj<-project.pca(x6G08_MD_226_raw,pcadata)
points(x6G08_MD_226_proj[1],x6G08_MD_226_proj[3],pch=20)
text(x6G08_MD_226_proj[1],x6G08_MD_226_proj[3],pos=4,label="226",col="blue")

x6G08_MD_227_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_227.txt")
x6G08_MD_227_proj<-project.pca(x6G08_MD_227_raw,pcadata)
points(x6G08_MD_227_proj[1],x6G08_MD_227_proj[3],pch=20)
text(x6G08_MD_227_proj[1],x6G08_MD_227_proj[3],pos=4,label="227",col="red")

x6G08_MD_228_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_228.txt")
x6G08_MD_228_proj<-project.pca(x6G08_MD_228_raw,pcadata)
points(x6G08_MD_228_proj[1],x6G08_MD_228_proj[3],pch=20)
text(x6G08_MD_228_proj[1],x6G08_MD_228_proj[3],pos=4,label="228",col="green")

x6G08_MD_229_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_229.txt")
x6G08_MD_229_proj<-project.pca(x6G08_MD_229_raw,pcadata)
points(x6G08_MD_229_proj[1],x6G08_MD_229_proj[3],pch=20)
text(x6G08_MD_229_proj[1],x6G08_MD_229_proj[3],pos=4,label="229",col="red")

x6G08_MD_230_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_230.txt")
x6G08_MD_230_proj<-project.pca(x6G08_MD_230_raw,pcadata)
points(x6G08_MD_230_proj[1],x6G08_MD_230_proj[3],pch=20)
text(x6G08_MD_230_proj[1],x6G08_MD_230_proj[3],pos=4,label="230",col="green")

x6G08_MD_231_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_231.txt")
x6G08_MD_231_proj<-project.pca(x6G08_MD_231_raw,pcadata)
points(x6G08_MD_231_proj[1],x6G08_MD_231_proj[3],pch=20)
text(x6G08_MD_231_proj[1],x6G08_MD_231_proj[3],pos=4,label="231",col="green")

x6G08_MD_232_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_232.txt")
x6G08_MD_232_proj<-project.pca(x6G08_MD_232_raw,pcadata)
points(x6G08_MD_232_proj[1],x6G08_MD_232_proj[3],pch=20)
text(x6G08_MD_232_proj[1],x6G08_MD_232_proj[3],pos=4,label="232",col="red")

x6G08_MD_233_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_233.txt")
x6G08_MD_233_proj<-project.pca(x6G08_MD_233_raw,pcadata)
points(x6G08_MD_233_proj[1],x6G08_MD_233_proj[3],pch=20)
text(x6G08_MD_233_proj[1],x6G08_MD_233_proj[3],pos=4,label="233",col="blue")

x6G08_MD_234_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_234.txt")
x6G08_MD_234_proj<-project.pca(x6G08_MD_234_raw,pcadata)
points(x6G08_MD_234_proj[1],x6G08_MD_234_proj[3],pch=20)
text(x6G08_MD_234_proj[1],x6G08_MD_234_proj[3],pos=4,label="234",col="green")

x6G08_MD_235_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_235.txt")
x6G08_MD_235_proj<-project.pca(x6G08_MD_235_raw,pcadata)
points(x6G08_MD_235_proj[1],x6G08_MD_235_proj[3],pch=20)
text(x6G08_MD_235_proj[1],x6G08_MD_235_proj[3],pos=4,label="235",col="green")

x6G08_MD_236_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_236.txt")
x6G08_MD_236_proj<-project.pca(x6G08_MD_236_raw,pcadata)
points(x6G08_MD_236_proj[1],x6G08_MD_236_proj[3],pch=20)
text(x6G08_MD_236_proj[1],x6G08_MD_236_proj[3],pos=4,label="236",col="blue")

x6G08_MD_237_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_237.txt")
x6G08_MD_237_proj<-project.pca(x6G08_MD_237_raw,pcadata)
points(x6G08_MD_237_proj[1],x6G08_MD_237_proj[3],pch=20)
text(x6G08_MD_237_proj[1],x6G08_MD_237_proj[3],pos=4,label="237",col="green")

x6G08_MD_238_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_238.txt")
x6G08_MD_238_proj<-project.pca(x6G08_MD_238_raw,pcadata)
points(x6G08_MD_238_proj[1],x6G08_MD_238_proj[3],pch=20)
text(x6G08_MD_238_proj[1],x6G08_MD_238_proj[3],pos=4,label="238",col="green")

x6G08_MD_239_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_239.txt")
x6G08_MD_239_proj<-project.pca(x6G08_MD_239_raw,pcadata)
points(x6G08_MD_239_proj[1],x6G08_MD_239_proj[3],pch=20)
text(x6G08_MD_239_proj[1],x6G08_MD_239_proj[3],pos=4,label="239",col="red")

x6G08_MD_240_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_240.txt")
x6G08_MD_240_proj<-project.pca(x6G08_MD_240_raw,pcadata)
points(x6G08_MD_240_proj[1],x6G08_MD_240_proj[3],pch=20)
text(x6G08_MD_240_proj[1],x6G08_MD_240_proj[3],pos=4,label="240",col="green")

x6G08_MD_241_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_241.txt")
x6G08_MD_241_proj<-project.pca(x6G08_MD_241_raw,pcadata)
points(x6G08_MD_241_proj[1],x6G08_MD_241_proj[3],pch=20)
text(x6G08_MD_241_proj[1],x6G08_MD_241_proj[3],pos=4,label="241",col="green")

x6G08_MD_242_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_242.txt")
x6G08_MD_242_proj<-project.pca(x6G08_MD_242_raw,pcadata)
points(x6G08_MD_242_proj[1],x6G08_MD_242_proj[3],pch=20)
text(x6G08_MD_242_proj[1],x6G08_MD_242_proj[3],pos=4,label="242",col="green")

x6G08_MD_243_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_243.txt")
x6G08_MD_243_proj<-project.pca(x6G08_MD_243_raw,pcadata)
points(x6G08_MD_243_proj[1],x6G08_MD_243_proj[3],pch=20)
text(x6G08_MD_243_proj[1],x6G08_MD_243_proj[3],pos=4,label="243",col="green")

x6G08_MD_244_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_244.txt")
x6G08_MD_244_proj<-project.pca(x6G08_MD_244_raw,pcadata)
points(x6G08_MD_244_proj[1],x6G08_MD_244_proj[3],pch=20)
text(x6G08_MD_244_proj[1],x6G08_MD_244_proj[3],pos=4,label="244",col="blue")

x6G08_MD_245_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_245.txt")
x6G08_MD_245_proj<-project.pca(x6G08_MD_245_raw,pcadata)
points(x6G08_MD_245_proj[1],x6G08_MD_245_proj[3],pch=20)
text(x6G08_MD_245_proj[1],x6G08_MD_245_proj[3],pos=4,label="245",col="green")

x6G08_MD_246_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_246.txt")
x6G08_MD_246_proj<-project.pca(x6G08_MD_246_raw,pcadata)
points(x6G08_MD_246_proj[1],x6G08_MD_246_proj[3],pch=20)
text(x6G08_MD_246_proj[1],x6G08_MD_246_proj[3],pos=4,label="246",col="green")

x6G08_MD_247_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_247.txt")
x6G08_MD_247_proj<-project.pca(x6G08_MD_247_raw,pcadata)
points(x6G08_MD_247_proj[1],x6G08_MD_247_proj[3],pch=20)
text(x6G08_MD_247_proj[1],x6G08_MD_247_proj[3],pos=4,label="247",col="blue")

x6G08_MD_248_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_248.txt")
x6G08_MD_248_proj<-project.pca(x6G08_MD_248_raw,pcadata)
points(x6G08_MD_248_proj[1],x6G08_MD_248_proj[3],pch=20)
text(x6G08_MD_248_proj[1],x6G08_MD_248_proj[3],pos=4,label="248",col="black")

x6G08_MD_249_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_249.txt")
x6G08_MD_249_proj<-project.pca(x6G08_MD_249_raw,pcadata)
points(x6G08_MD_249_proj[1],x6G08_MD_249_proj[3],pch=20)
text(x6G08_MD_249_proj[1],x6G08_MD_249_proj[3],pos=4,label="249",col="blue")

x6G08_MD_250_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_250.txt")
x6G08_MD_250_proj<-project.pca(x6G08_MD_250_raw,pcadata)
points(x6G08_MD_250_proj[1],x6G08_MD_250_proj[3],pch=20)
text(x6G08_MD_250_proj[1],x6G08_MD_250_proj[3],pos=4,label="250",col="black")

x6G08_MD_251_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_251.txt")
x6G08_MD_251_proj<-project.pca(x6G08_MD_251_raw,pcadata)
points(x6G08_MD_251_proj[1],x6G08_MD_251_proj[3],pch=20)
text(x6G08_MD_251_proj[1],x6G08_MD_251_proj[3],pos=4,label="251",col="black")

x6G08_MD_252_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_252.txt")
x6G08_MD_252_proj<-project.pca(x6G08_MD_252_raw,pcadata)
points(x6G08_MD_252_proj[1],x6G08_MD_252_proj[3],pch=20)
text(x6G08_MD_252_proj[1],x6G08_MD_252_proj[3],pos=4,label="252",col="blue")

x6G08_MD_253_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_253.txt")
x6G08_MD_253_proj<-project.pca(x6G08_MD_253_raw,pcadata)
points(x6G08_MD_253_proj[1],x6G08_MD_253_proj[3],pch=20)
text(x6G08_MD_253_proj[1],x6G08_MD_253_proj[3],pos=4,label="253",col="blue")

x6G08_MD_254_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_254.txt")
x6G08_MD_254_proj<-project.pca(x6G08_MD_254_raw,pcadata)
points(x6G08_MD_254_proj[1],x6G08_MD_254_proj[3],pch=20)
text(x6G08_MD_254_proj[1],x6G08_MD_254_proj[3],pos=4,label="254",col="black")

x6G08_MD_255_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_255.txt")
x6G08_MD_255_proj<-project.pca(x6G08_MD_255_raw,pcadata)
points(x6G08_MD_255_proj[1],x6G08_MD_255_proj[3],pch=20)
text(x6G08_MD_255_proj[1],x6G08_MD_255_proj[3],pos=4,label="255",col="blue")

x6G08_MD_256_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_256.txt")
x6G08_MD_256_proj<-project.pca(x6G08_MD_256_raw,pcadata)
points(x6G08_MD_256_proj[1],x6G08_MD_256_proj[3],pch=20)
text(x6G08_MD_256_proj[1],x6G08_MD_256_proj[3],pos=4,label="256",col="blue")

x6G08_MD_257_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_257.txt")
x6G08_MD_257_proj<-project.pca(x6G08_MD_257_raw,pcadata)
points(x6G08_MD_257_proj[1],x6G08_MD_257_proj[3],pch=20)
text(x6G08_MD_257_proj[1],x6G08_MD_257_proj[3],pos=4,label="257",col="green")

x6G08_MD_258_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_258.txt")
x6G08_MD_258_proj<-project.pca(x6G08_MD_258_raw,pcadata)
points(x6G08_MD_258_proj[1],x6G08_MD_258_proj[3],pch=20)
text(x6G08_MD_258_proj[1],x6G08_MD_258_proj[3],pos=4,label="258",col="green")

x6G08_MD_259_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_259.txt")
x6G08_MD_259_proj<-project.pca(x6G08_MD_259_raw,pcadata)
points(x6G08_MD_259_proj[1],x6G08_MD_259_proj[3],pch=20)
text(x6G08_MD_259_proj[1],x6G08_MD_259_proj[3],pos=4,label="259",col="green")

x6G08_MD_260_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_260.txt")
x6G08_MD_260_proj<-project.pca(x6G08_MD_260_raw,pcadata)
points(x6G08_MD_260_proj[1],x6G08_MD_260_proj[3],pch=20)
text(x6G08_MD_260_proj[1],x6G08_MD_260_proj[3],pos=4,label="260",col="blue")

x6G08_MD_261_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_261.txt")
x6G08_MD_261_proj<-project.pca(x6G08_MD_261_raw,pcadata)
points(x6G08_MD_261_proj[1],x6G08_MD_261_proj[3],pch=20)
text(x6G08_MD_261_proj[1],x6G08_MD_261_proj[3],pos=4,label="261",col="blue")

x6G08_MD_262_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_262.txt")
x6G08_MD_262_proj<-project.pca(x6G08_MD_262_raw,pcadata)
points(x6G08_MD_262_proj[1],x6G08_MD_262_proj[3],pch=20)
text(x6G08_MD_262_proj[1],x6G08_MD_262_proj[3],pos=4,label="262",col="blue")

x6G08_MD_263_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_263.txt")
x6G08_MD_263_proj<-project.pca(x6G08_MD_263_raw,pcadata)
points(x6G08_MD_263_proj[1],x6G08_MD_263_proj[3],pch=20)
text(x6G08_MD_263_proj[1],x6G08_MD_263_proj[3],pos=4,label="263",col="blue")

x6G08_MD_264_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_264.txt")
x6G08_MD_264_proj<-project.pca(x6G08_MD_264_raw,pcadata)
points(x6G08_MD_264_proj[1],x6G08_MD_264_proj[3],pch=20)
text(x6G08_MD_264_proj[1],x6G08_MD_264_proj[3],pos=4,label="264",col="blue")

x6G08_MD_265_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_265.txt")
x6G08_MD_265_proj<-project.pca(x6G08_MD_265_raw,pcadata)
points(x6G08_MD_265_proj[1],x6G08_MD_265_proj[3],pch=20)
text(x6G08_MD_265_proj[1],x6G08_MD_265_proj[3],pos=4,label="265",col="blue")

x6G08_MD_266_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_266.txt")
x6G08_MD_266_proj<-project.pca(x6G08_MD_266_raw,pcadata)
points(x6G08_MD_266_proj[1],x6G08_MD_266_proj[3],pch=20)
text(x6G08_MD_266_proj[1],x6G08_MD_266_proj[3],pos=4,label="266",col="blue")

x6G08_MD_267_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_267.txt")
x6G08_MD_267_proj<-project.pca(x6G08_MD_267_raw,pcadata)
points(x6G08_MD_267_proj[1],x6G08_MD_267_proj[3],pch=20)
text(x6G08_MD_267_proj[1],x6G08_MD_267_proj[3],pos=4,label="267",col="blue")

x6G08_MD_268_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_268.txt")
x6G08_MD_268_proj<-project.pca(x6G08_MD_268_raw,pcadata)
points(x6G08_MD_268_proj[1],x6G08_MD_268_proj[3],pch=20)
text(x6G08_MD_268_proj[1],x6G08_MD_268_proj[3],pos=4,label="268",col="blue")

x6G08_MD_269_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_269.txt")
x6G08_MD_269_proj<-project.pca(x6G08_MD_269_raw,pcadata)
points(x6G08_MD_269_proj[1],x6G08_MD_269_proj[3],pch=20)
text(x6G08_MD_269_proj[1],x6G08_MD_269_proj[3],pos=4,label="269",col="blue")

x6G08_MD_270_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_270.txt")
x6G08_MD_270_proj<-project.pca(x6G08_MD_270_raw,pcadata)
points(x6G08_MD_270_proj[1],x6G08_MD_270_proj[3],pch=20)
text(x6G08_MD_270_proj[1],x6G08_MD_270_proj[3],pos=4,label="270",col="black")

x6G08_MD_271_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_271.txt")
x6G08_MD_271_proj<-project.pca(x6G08_MD_271_raw,pcadata)
points(x6G08_MD_271_proj[1],x6G08_MD_271_proj[3],pch=20)
text(x6G08_MD_271_proj[1],x6G08_MD_271_proj[3],pos=4,label="271",col="green")

x6G08_MD_272_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_272.txt")
x6G08_MD_272_proj<-project.pca(x6G08_MD_272_raw,pcadata)
points(x6G08_MD_272_proj[1],x6G08_MD_272_proj[3],pch=20)
text(x6G08_MD_272_proj[1],x6G08_MD_272_proj[3],pos=4,label="272",col="blue")

x6G08_MD_273_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_273.txt")
x6G08_MD_273_proj<-project.pca(x6G08_MD_273_raw,pcadata)
points(x6G08_MD_273_proj[1],x6G08_MD_273_proj[3],pch=20)
text(x6G08_MD_273_proj[1],x6G08_MD_273_proj[3],pos=4,label="273",col="blue")

x6G08_MD_274_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_274.txt")
x6G08_MD_274_proj<-project.pca(x6G08_MD_274_raw,pcadata)
points(x6G08_MD_274_proj[1],x6G08_MD_274_proj[3],pch=20)
text(x6G08_MD_274_proj[1],x6G08_MD_274_proj[3],pos=4,label="274",col="blue")

x6G08_MD_275_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_275.txt")
x6G08_MD_275_proj<-project.pca(x6G08_MD_275_raw,pcadata)
points(x6G08_MD_275_proj[1],x6G08_MD_275_proj[3],pch=20)
text(x6G08_MD_275_proj[1],x6G08_MD_275_proj[3],pos=4,label="275",col="green")

x6G08_MD_276_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_276.txt")
x6G08_MD_276_proj<-project.pca(x6G08_MD_276_raw,pcadata)
points(x6G08_MD_276_proj[1],x6G08_MD_276_proj[3],pch=20)
text(x6G08_MD_276_proj[1],x6G08_MD_276_proj[3],pos=4,label="276",col="blue")

x6G08_MD_277_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_277.txt")
x6G08_MD_277_proj<-project.pca(x6G08_MD_277_raw,pcadata)
points(x6G08_MD_277_proj[1],x6G08_MD_277_proj[3],pch=20)
text(x6G08_MD_277_proj[1],x6G08_MD_277_proj[3],pos=4,label="277",col="green")

x6G08_MD_278_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_278.txt")
x6G08_MD_278_proj<-project.pca(x6G08_MD_278_raw,pcadata)
points(x6G08_MD_278_proj[1],x6G08_MD_278_proj[3],pch=20)
text(x6G08_MD_278_proj[1],x6G08_MD_278_proj[3],pos=4,label="278",col="red")

x6G08_MD_279_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_279.txt")
x6G08_MD_279_proj<-project.pca(x6G08_MD_279_raw,pcadata)
points(x6G08_MD_279_proj[1],x6G08_MD_279_proj[3],pch=20)
text(x6G08_MD_279_proj[1],x6G08_MD_279_proj[3],pos=4,label="279",col="green")

x6G08_MD_280_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_280.txt")
x6G08_MD_280_proj<-project.pca(x6G08_MD_280_raw,pcadata)
points(x6G08_MD_280_proj[1],x6G08_MD_280_proj[3],pch=20)
text(x6G08_MD_280_proj[1],x6G08_MD_280_proj[3],pos=4,label="280",col="blue")

x6G08_MD_281_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_281.txt")
x6G08_MD_281_proj<-project.pca(x6G08_MD_281_raw,pcadata)
points(x6G08_MD_281_proj[1],x6G08_MD_281_proj[3],pch=20)
text(x6G08_MD_281_proj[1],x6G08_MD_281_proj[3],pos=4,label="281",col="green")

x6G08_MD_282_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_282.txt")
x6G08_MD_282_proj<-project.pca(x6G08_MD_282_raw,pcadata)
points(x6G08_MD_282_proj[1],x6G08_MD_282_proj[3],pch=20)
text(x6G08_MD_282_proj[1],x6G08_MD_282_proj[3],pos=4,label="282",col="green")

x6G08_MD_283_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_283.txt")
x6G08_MD_283_proj<-project.pca(x6G08_MD_283_raw,pcadata)
points(x6G08_MD_283_proj[1],x6G08_MD_283_proj[3],pch=20)
text(x6G08_MD_283_proj[1],x6G08_MD_283_proj[3],pos=4,label="283",col="green")

x6G08_MD_284_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_284.txt")
x6G08_MD_284_proj<-project.pca(x6G08_MD_284_raw,pcadata)
points(x6G08_MD_284_proj[1],x6G08_MD_284_proj[3],pch=20)
text(x6G08_MD_284_proj[1],x6G08_MD_284_proj[3],pos=4,label="284",col="green")

x6G08_MD_285_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_285.txt")
x6G08_MD_285_proj<-project.pca(x6G08_MD_285_raw,pcadata)
points(x6G08_MD_285_proj[1],x6G08_MD_285_proj[3],pch=20)
text(x6G08_MD_285_proj[1],x6G08_MD_285_proj[3],pos=4,label="285",col="red")

x6G08_MD_286_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_286.txt")
x6G08_MD_286_proj<-project.pca(x6G08_MD_286_raw,pcadata)
points(x6G08_MD_286_proj[1],x6G08_MD_286_proj[3],pch=20)
text(x6G08_MD_286_proj[1],x6G08_MD_286_proj[3],pos=4,label="286",col="red")

x6G08_MD_287_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_287.txt")
x6G08_MD_287_proj<-project.pca(x6G08_MD_287_raw,pcadata)
points(x6G08_MD_287_proj[1],x6G08_MD_287_proj[3],pch=20)
text(x6G08_MD_287_proj[1],x6G08_MD_287_proj[3],pos=4,label="287",col="red")

x6G08_MD_288_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_288.txt")
x6G08_MD_288_proj<-project.pca(x6G08_MD_288_raw,pcadata)
points(x6G08_MD_288_proj[1],x6G08_MD_288_proj[3],pch=20)
text(x6G08_MD_288_proj[1],x6G08_MD_288_proj[3],pos=4,label="288",col="green")

x6G08_MD_289_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_289.txt")
x6G08_MD_289_proj<-project.pca(x6G08_MD_289_raw,pcadata)
points(x6G08_MD_289_proj[1],x6G08_MD_289_proj[3],pch=20)
text(x6G08_MD_289_proj[1],x6G08_MD_289_proj[3],pos=4,label="289",col="blue")

x6G08_MD_290_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_290.txt")
x6G08_MD_290_proj<-project.pca(x6G08_MD_290_raw,pcadata)
points(x6G08_MD_290_proj[1],x6G08_MD_290_proj[3],pch=20)
text(x6G08_MD_290_proj[1],x6G08_MD_290_proj[3],pos=4,label="290",col="black")

x6G08_MD_291_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_291.txt")
x6G08_MD_291_proj<-project.pca(x6G08_MD_291_raw,pcadata)
points(x6G08_MD_291_proj[1],x6G08_MD_291_proj[3],pch=20)
text(x6G08_MD_291_proj[1],x6G08_MD_291_proj[3],pos=4,label="291",col="blue")

x6G08_MD_292_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_292.txt")
x6G08_MD_292_proj<-project.pca(x6G08_MD_292_raw,pcadata)
points(x6G08_MD_292_proj[1],x6G08_MD_292_proj[3],pch=20)
text(x6G08_MD_292_proj[1],x6G08_MD_292_proj[3],pos=4,label="292",col="black")

x6G08_MD_293_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_293.txt")
x6G08_MD_293_proj<-project.pca(x6G08_MD_293_raw,pcadata)
points(x6G08_MD_293_proj[1],x6G08_MD_293_proj[3],pch=20)
text(x6G08_MD_293_proj[1],x6G08_MD_293_proj[3],pos=4,label="293",col="blue")

x6G08_MD_294_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_294.txt")
x6G08_MD_294_proj<-project.pca(x6G08_MD_294_raw,pcadata)
points(x6G08_MD_294_proj[1],x6G08_MD_294_proj[3],pch=20)
text(x6G08_MD_294_proj[1],x6G08_MD_294_proj[3],pos=4,label="294",col="red")

x6G08_MD_295_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_295.txt")
x6G08_MD_295_proj<-project.pca(x6G08_MD_295_raw,pcadata)
points(x6G08_MD_295_proj[1],x6G08_MD_295_proj[3],pch=20)
text(x6G08_MD_295_proj[1],x6G08_MD_295_proj[3],pos=4,label="295",col="green")

x6G08_MD_296_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_296.txt")
x6G08_MD_296_proj<-project.pca(x6G08_MD_296_raw,pcadata)
points(x6G08_MD_296_proj[1],x6G08_MD_296_proj[3],pch=20)
text(x6G08_MD_296_proj[1],x6G08_MD_296_proj[3],pos=4,label="296",col="blue")

x6G08_MD_297_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_297.txt")
x6G08_MD_297_proj<-project.pca(x6G08_MD_297_raw,pcadata)
points(x6G08_MD_297_proj[1],x6G08_MD_297_proj[3],pch=20)
text(x6G08_MD_297_proj[1],x6G08_MD_297_proj[3],pos=4,label="297",col="red")

x6G08_MD_298_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_298.txt")
x6G08_MD_298_proj<-project.pca(x6G08_MD_298_raw,pcadata)
points(x6G08_MD_298_proj[1],x6G08_MD_298_proj[3],pch=20)
text(x6G08_MD_298_proj[1],x6G08_MD_298_proj[3],pos=4,label="298",col="red")

x6G08_MD_299_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_299.txt")
x6G08_MD_299_proj<-project.pca(x6G08_MD_299_raw,pcadata)
points(x6G08_MD_299_proj[1],x6G08_MD_299_proj[3],pch=20)
text(x6G08_MD_299_proj[1],x6G08_MD_299_proj[3],pos=4,label="299",col="green")

x6G08_MD_300_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_300.txt")
x6G08_MD_300_proj<-project.pca(x6G08_MD_300_raw,pcadata)
points(x6G08_MD_300_proj[1],x6G08_MD_300_proj[3],pch=20)
text(x6G08_MD_300_proj[1],x6G08_MD_300_proj[3],pos=4,label="300",col="green")

x6G08_MD_301_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_301.txt")
x6G08_MD_301_proj<-project.pca(x6G08_MD_301_raw,pcadata)
points(x6G08_MD_301_proj[1],x6G08_MD_301_proj[3],pch=20)
text(x6G08_MD_301_proj[1],x6G08_MD_301_proj[3],pos=4,label="301",col="blue")

x6G08_MD_302_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_302.txt")
x6G08_MD_302_proj<-project.pca(x6G08_MD_302_raw,pcadata)
points(x6G08_MD_302_proj[1],x6G08_MD_302_proj[3],pch=20)
text(x6G08_MD_302_proj[1],x6G08_MD_302_proj[3],pos=4,label="302",col="green")

x6G08_MD_303_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_303.txt")
x6G08_MD_303_proj<-project.pca(x6G08_MD_303_raw,pcadata)
points(x6G08_MD_303_proj[1],x6G08_MD_303_proj[3],pch=20)
text(x6G08_MD_303_proj[1],x6G08_MD_303_proj[3],pos=4,label="303",col="green")

x6G08_MD_304_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_304.txt")
x6G08_MD_304_proj<-project.pca(x6G08_MD_304_raw,pcadata)
points(x6G08_MD_304_proj[1],x6G08_MD_304_proj[3],pch=20)
text(x6G08_MD_304_proj[1],x6G08_MD_304_proj[3],pos=4,label="304",col="green")

x6G08_MD_305_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_305.txt")
x6G08_MD_305_proj<-project.pca(x6G08_MD_305_raw,pcadata)
points(x6G08_MD_305_proj[1],x6G08_MD_305_proj[3],pch=20)
text(x6G08_MD_305_proj[1],x6G08_MD_305_proj[3],pos=4,label="305",col="red")

x6G08_MD_306_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_306.txt")
x6G08_MD_306_proj<-project.pca(x6G08_MD_306_raw,pcadata)
points(x6G08_MD_306_proj[1],x6G08_MD_306_proj[3],pch=20)
text(x6G08_MD_306_proj[1],x6G08_MD_306_proj[3],pos=4,label="306",col="green")

x6G08_MD_307_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_307.txt")
x6G08_MD_307_proj<-project.pca(x6G08_MD_307_raw,pcadata)
points(x6G08_MD_307_proj[1],x6G08_MD_307_proj[3],pch=20)
text(x6G08_MD_307_proj[1],x6G08_MD_307_proj[3],pos=4,label="307",col="green")

x6G08_MD_308_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_308.txt")
x6G08_MD_308_proj<-project.pca(x6G08_MD_308_raw,pcadata)
points(x6G08_MD_308_proj[1],x6G08_MD_308_proj[3],pch=20)
text(x6G08_MD_308_proj[1],x6G08_MD_308_proj[3],pos=4,label="308",col="red")

x6G08_MD_309_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_309.txt")
x6G08_MD_309_proj<-project.pca(x6G08_MD_309_raw,pcadata)
points(x6G08_MD_309_proj[1],x6G08_MD_309_proj[3],pch=20)
text(x6G08_MD_309_proj[1],x6G08_MD_309_proj[3],pos=4,label="309",col="green")

x6G08_MD_310_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_310.txt")
x6G08_MD_310_proj<-project.pca(x6G08_MD_310_raw,pcadata)
points(x6G08_MD_310_proj[1],x6G08_MD_310_proj[3],pch=20)
text(x6G08_MD_310_proj[1],x6G08_MD_310_proj[3],pos=4,label="310",col="green")

x6G08_MD_311_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_311.txt")
x6G08_MD_311_proj<-project.pca(x6G08_MD_311_raw,pcadata)
points(x6G08_MD_311_proj[1],x6G08_MD_311_proj[3],pch=20)
text(x6G08_MD_311_proj[1],x6G08_MD_311_proj[3],pos=4,label="311",col="blue")

x6G08_MD_312_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_312.txt")
x6G08_MD_312_proj<-project.pca(x6G08_MD_312_raw,pcadata)
points(x6G08_MD_312_proj[1],x6G08_MD_312_proj[3],pch=20)
text(x6G08_MD_312_proj[1],x6G08_MD_312_proj[3],pos=4,label="312",col="blue")

x6G08_MD_313_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_313.txt")
x6G08_MD_313_proj<-project.pca(x6G08_MD_313_raw,pcadata)
points(x6G08_MD_313_proj[1],x6G08_MD_313_proj[3],pch=20)
text(x6G08_MD_313_proj[1],x6G08_MD_313_proj[3],pos=4,label="313",col="blue")

x6G08_MD_314_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_314.txt")
x6G08_MD_314_proj<-project.pca(x6G08_MD_314_raw,pcadata)
points(x6G08_MD_314_proj[1],x6G08_MD_314_proj[3],pch=20)
text(x6G08_MD_314_proj[1],x6G08_MD_314_proj[3],pos=4,label="314",col="green")

x6G08_MD_315_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_315.txt")
x6G08_MD_315_proj<-project.pca(x6G08_MD_315_raw,pcadata)
points(x6G08_MD_315_proj[1],x6G08_MD_315_proj[3],pch=20)
text(x6G08_MD_315_proj[1],x6G08_MD_315_proj[3],pos=4,label="315",col="blue")

x6G08_MD_316_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_316.txt")
x6G08_MD_316_proj<-project.pca(x6G08_MD_316_raw,pcadata)
points(x6G08_MD_316_proj[1],x6G08_MD_316_proj[3],pch=20)
text(x6G08_MD_316_proj[1],x6G08_MD_316_proj[3],pos=4,label="316",col="blue")

x6G08_MD_317_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_317.txt")
x6G08_MD_317_proj<-project.pca(x6G08_MD_317_raw,pcadata)
points(x6G08_MD_317_proj[1],x6G08_MD_317_proj[3],pch=20)
text(x6G08_MD_317_proj[1],x6G08_MD_317_proj[3],pos=4,label="317",col="green")

x6G08_MD_318_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_318.txt")
x6G08_MD_318_proj<-project.pca(x6G08_MD_318_raw,pcadata)
points(x6G08_MD_318_proj[1],x6G08_MD_318_proj[3],pch=20)
text(x6G08_MD_318_proj[1],x6G08_MD_318_proj[3],pos=4,label="318",col="green")

x6G08_MD_319_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_319.txt")
x6G08_MD_319_proj<-project.pca(x6G08_MD_319_raw,pcadata)
points(x6G08_MD_319_proj[1],x6G08_MD_319_proj[3],pch=20)
text(x6G08_MD_319_proj[1],x6G08_MD_319_proj[3],pos=4,label="319",col="green")

x6G08_MD_320_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_320.txt")
x6G08_MD_320_proj<-project.pca(x6G08_MD_320_raw,pcadata)
points(x6G08_MD_320_proj[1],x6G08_MD_320_proj[3],pch=20)
text(x6G08_MD_320_proj[1],x6G08_MD_320_proj[3],pos=4,label="320",col="green")

x6G08_MD_321_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_321.txt")
x6G08_MD_321_proj<-project.pca(x6G08_MD_321_raw,pcadata)
points(x6G08_MD_321_proj[1],x6G08_MD_321_proj[3],pch=20)
text(x6G08_MD_321_proj[1],x6G08_MD_321_proj[3],pos=4,label="321",col="blue")

x6G08_MD_322_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_322.txt")
x6G08_MD_322_proj<-project.pca(x6G08_MD_322_raw,pcadata)
points(x6G08_MD_322_proj[1],x6G08_MD_322_proj[3],pch=20)
text(x6G08_MD_322_proj[1],x6G08_MD_322_proj[3],pos=4,label="322",col="green")

x6G08_MD_323_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_323.txt")
x6G08_MD_323_proj<-project.pca(x6G08_MD_323_raw,pcadata)
points(x6G08_MD_323_proj[1],x6G08_MD_323_proj[3],pch=20)
text(x6G08_MD_323_proj[1],x6G08_MD_323_proj[3],pos=4,label="323",col="red")

x6G08_MD_324_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_324.txt")
x6G08_MD_324_proj<-project.pca(x6G08_MD_324_raw,pcadata)
points(x6G08_MD_324_proj[1],x6G08_MD_324_proj[3],pch=20)
text(x6G08_MD_324_proj[1],x6G08_MD_324_proj[3],pos=4,label="324",col="green")

x6G08_MD_325_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_325.txt")
x6G08_MD_325_proj<-project.pca(x6G08_MD_325_raw,pcadata)
points(x6G08_MD_325_proj[1],x6G08_MD_325_proj[3],pch=20)
text(x6G08_MD_325_proj[1],x6G08_MD_325_proj[3],pos=4,label="325",col="red")

x6G08_MD_326_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_326.txt")
x6G08_MD_326_proj<-project.pca(x6G08_MD_326_raw,pcadata)
points(x6G08_MD_326_proj[1],x6G08_MD_326_proj[3],pch=20)
text(x6G08_MD_326_proj[1],x6G08_MD_326_proj[3],pos=4,label="326",col="red")

x6G08_MD_327_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_327.txt")
x6G08_MD_327_proj<-project.pca(x6G08_MD_327_raw,pcadata)
points(x6G08_MD_327_proj[1],x6G08_MD_327_proj[3],pch=20)
text(x6G08_MD_327_proj[1],x6G08_MD_327_proj[3],pos=4,label="327",col="blue")

x6G08_MD_328_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_328.txt")
x6G08_MD_328_proj<-project.pca(x6G08_MD_328_raw,pcadata)
points(x6G08_MD_328_proj[1],x6G08_MD_328_proj[3],pch=20)
text(x6G08_MD_328_proj[1],x6G08_MD_328_proj[3],pos=4,label="328",col="red")

x6G08_MD_329_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_329.txt")
x6G08_MD_329_proj<-project.pca(x6G08_MD_329_raw,pcadata)
points(x6G08_MD_329_proj[1],x6G08_MD_329_proj[3],pch=20)
text(x6G08_MD_329_proj[1],x6G08_MD_329_proj[3],pos=4,label="329",col="red")

x6G08_MD_330_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_330.txt")
x6G08_MD_330_proj<-project.pca(x6G08_MD_330_raw,pcadata)
points(x6G08_MD_330_proj[1],x6G08_MD_330_proj[3],pch=20)
text(x6G08_MD_330_proj[1],x6G08_MD_330_proj[3],pos=4,label="330",col="green")

x6G08_MD_331_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_331.txt")
x6G08_MD_331_proj<-project.pca(x6G08_MD_331_raw,pcadata)
points(x6G08_MD_331_proj[1],x6G08_MD_331_proj[3],pch=20)
text(x6G08_MD_331_proj[1],x6G08_MD_331_proj[3],pos=4,label="331",col="green")

x6G08_MD_332_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_332.txt")
x6G08_MD_332_proj<-project.pca(x6G08_MD_332_raw,pcadata)
points(x6G08_MD_332_proj[1],x6G08_MD_332_proj[3],pch=20)
text(x6G08_MD_332_proj[1],x6G08_MD_332_proj[3],pos=4,label="332",col="red")

x6G08_MD_333_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_333.txt")
x6G08_MD_333_proj<-project.pca(x6G08_MD_333_raw,pcadata)
points(x6G08_MD_333_proj[1],x6G08_MD_333_proj[3],pch=20)
text(x6G08_MD_333_proj[1],x6G08_MD_333_proj[3],pos=4,label="333",col="green")

x6G08_MD_334_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_334.txt")
x6G08_MD_334_proj<-project.pca(x6G08_MD_334_raw,pcadata)
points(x6G08_MD_334_proj[1],x6G08_MD_334_proj[3],pch=20)
text(x6G08_MD_334_proj[1],x6G08_MD_334_proj[3],pos=4,label="334",col="green")

x6G08_MD_335_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_335.txt")
x6G08_MD_335_proj<-project.pca(x6G08_MD_335_raw,pcadata)
points(x6G08_MD_335_proj[1],x6G08_MD_335_proj[3],pch=20)
text(x6G08_MD_335_proj[1],x6G08_MD_335_proj[3],pos=4,label="335",col="green")

x6G08_MD_336_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_336.txt")
x6G08_MD_336_proj<-project.pca(x6G08_MD_336_raw,pcadata)
points(x6G08_MD_336_proj[1],x6G08_MD_336_proj[3],pch=20)
text(x6G08_MD_336_proj[1],x6G08_MD_336_proj[3],pos=4,label="336",col="green")

x6G08_MD_337_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_337.txt")
x6G08_MD_337_proj<-project.pca(x6G08_MD_337_raw,pcadata)
points(x6G08_MD_337_proj[1],x6G08_MD_337_proj[3],pch=20)
text(x6G08_MD_337_proj[1],x6G08_MD_337_proj[3],pos=4,label="337",col="green")

x6G08_MD_338_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_338.txt")
x6G08_MD_338_proj<-project.pca(x6G08_MD_338_raw,pcadata)
points(x6G08_MD_338_proj[1],x6G08_MD_338_proj[3],pch=20)
text(x6G08_MD_338_proj[1],x6G08_MD_338_proj[3],pos=4,label="338",col="blue")

x6G08_MD_339_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_339.txt")
x6G08_MD_339_proj<-project.pca(x6G08_MD_339_raw,pcadata)
points(x6G08_MD_339_proj[1],x6G08_MD_339_proj[3],pch=20)
text(x6G08_MD_339_proj[1],x6G08_MD_339_proj[3],pos=4,label="339",col="black")

x6G08_MD_340_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_340.txt")
x6G08_MD_340_proj<-project.pca(x6G08_MD_340_raw,pcadata)
points(x6G08_MD_340_proj[1],x6G08_MD_340_proj[3],pch=20)
text(x6G08_MD_340_proj[1],x6G08_MD_340_proj[3],pos=4,label="340",col="blue")

x6G08_MD_341_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_341.txt")
x6G08_MD_341_proj<-project.pca(x6G08_MD_341_raw,pcadata)
points(x6G08_MD_341_proj[1],x6G08_MD_341_proj[3],pch=20)
text(x6G08_MD_341_proj[1],x6G08_MD_341_proj[3],pos=4,label="341",col="blue")

x6G08_MD_342_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_342.txt")
x6G08_MD_342_proj<-project.pca(x6G08_MD_342_raw,pcadata)
points(x6G08_MD_342_proj[1],x6G08_MD_342_proj[3],pch=20)
text(x6G08_MD_342_proj[1],x6G08_MD_342_proj[3],pos=4,label="342",col="green")

x6G08_MD_343_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_343.txt")
x6G08_MD_343_proj<-project.pca(x6G08_MD_343_raw,pcadata)
points(x6G08_MD_343_proj[1],x6G08_MD_343_proj[3],pch=20)
text(x6G08_MD_343_proj[1],x6G08_MD_343_proj[3],pos=4,label="343",col="black")

x6G08_MD_344_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_344.txt")
x6G08_MD_344_proj<-project.pca(x6G08_MD_344_raw,pcadata)
points(x6G08_MD_344_proj[1],x6G08_MD_344_proj[3],pch=20)
text(x6G08_MD_344_proj[1],x6G08_MD_344_proj[3],pos=4,label="344",col="blue")

x6G08_MD_345_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_345.txt")
x6G08_MD_345_proj<-project.pca(x6G08_MD_345_raw,pcadata)
points(x6G08_MD_345_proj[1],x6G08_MD_345_proj[3],pch=20)
text(x6G08_MD_345_proj[1],x6G08_MD_345_proj[3],pos=4,label="345",col="blue")

x6G08_MD_346_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_346.txt")
x6G08_MD_346_proj<-project.pca(x6G08_MD_346_raw,pcadata)
points(x6G08_MD_346_proj[1],x6G08_MD_346_proj[3],pch=20)
text(x6G08_MD_346_proj[1],x6G08_MD_346_proj[3],pos=4,label="346",col="black")

x6G08_MD_347_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_347.txt")
x6G08_MD_347_proj<-project.pca(x6G08_MD_347_raw,pcadata)
points(x6G08_MD_347_proj[1],x6G08_MD_347_proj[3],pch=20)
text(x6G08_MD_347_proj[1],x6G08_MD_347_proj[3],pos=4,label="347",col="blue")

x6G08_MD_348_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_348.txt")
x6G08_MD_348_proj<-project.pca(x6G08_MD_348_raw,pcadata)
points(x6G08_MD_348_proj[1],x6G08_MD_348_proj[3],pch=20)
text(x6G08_MD_348_proj[1],x6G08_MD_348_proj[3],pos=4,label="348",col="black")

x6G08_MD_349_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_349.txt")
x6G08_MD_349_proj<-project.pca(x6G08_MD_349_raw,pcadata)
points(x6G08_MD_349_proj[1],x6G08_MD_349_proj[3],pch=20)
text(x6G08_MD_349_proj[1],x6G08_MD_349_proj[3],pos=4,label="349",col="black")

x6G08_MD_350_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_350.txt")
x6G08_MD_350_proj<-project.pca(x6G08_MD_350_raw,pcadata)
points(x6G08_MD_350_proj[1],x6G08_MD_350_proj[3],pch=20)
text(x6G08_MD_350_proj[1],x6G08_MD_350_proj[3],pos=4,label="350",col="black")

x6G08_MD_351_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_351.txt")
x6G08_MD_351_proj<-project.pca(x6G08_MD_351_raw,pcadata)
points(x6G08_MD_351_proj[1],x6G08_MD_351_proj[3],pch=20)
text(x6G08_MD_351_proj[1],x6G08_MD_351_proj[3],pos=4,label="351",col="blue")

x6G08_MD_352_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_352.txt")
x6G08_MD_352_proj<-project.pca(x6G08_MD_352_raw,pcadata)
points(x6G08_MD_352_proj[1],x6G08_MD_352_proj[3],pch=20)
text(x6G08_MD_352_proj[1],x6G08_MD_352_proj[3],pos=4,label="352",col="black")

x6G08_MD_353_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_353.txt")
x6G08_MD_353_proj<-project.pca(x6G08_MD_353_raw,pcadata)
points(x6G08_MD_353_proj[1],x6G08_MD_353_proj[3],pch=20)
text(x6G08_MD_353_proj[1],x6G08_MD_353_proj[3],pos=4,label="353",col="blue")

x6G08_MD_354_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_354.txt")
x6G08_MD_354_proj<-project.pca(x6G08_MD_354_raw,pcadata)
points(x6G08_MD_354_proj[1],x6G08_MD_354_proj[3],pch=20)
text(x6G08_MD_354_proj[1],x6G08_MD_354_proj[3],pos=4,label="354",col="black")

x6G08_MD_355_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_355.txt")
x6G08_MD_355_proj<-project.pca(x6G08_MD_355_raw,pcadata)
points(x6G08_MD_355_proj[1],x6G08_MD_355_proj[3],pch=20)
text(x6G08_MD_355_proj[1],x6G08_MD_355_proj[3],pos=4,label="355",col="black")

x6G08_MD_356_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_356.txt")
x6G08_MD_356_proj<-project.pca(x6G08_MD_356_raw,pcadata)
points(x6G08_MD_356_proj[1],x6G08_MD_356_proj[3],pch=20)
text(x6G08_MD_356_proj[1],x6G08_MD_356_proj[3],pos=4,label="356",col="green")

x6G08_MD_357_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_357.txt")
x6G08_MD_357_proj<-project.pca(x6G08_MD_357_raw,pcadata)
points(x6G08_MD_357_proj[1],x6G08_MD_357_proj[3],pch=20)
text(x6G08_MD_357_proj[1],x6G08_MD_357_proj[3],pos=4,label="357",col="blue")

x6G08_MD_358_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_358.txt")
x6G08_MD_358_proj<-project.pca(x6G08_MD_358_raw,pcadata)
points(x6G08_MD_358_proj[1],x6G08_MD_358_proj[3],pch=20)
text(x6G08_MD_358_proj[1],x6G08_MD_358_proj[3],pos=4,label="358",col="black")

x6G08_MD_359_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_359.txt")
x6G08_MD_359_proj<-project.pca(x6G08_MD_359_raw,pcadata)
points(x6G08_MD_359_proj[1],x6G08_MD_359_proj[3],pch=20)
text(x6G08_MD_359_proj[1],x6G08_MD_359_proj[3],pos=4,label="359",col="black")

x6G08_MD_360_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_360.txt")
x6G08_MD_360_proj<-project.pca(x6G08_MD_360_raw,pcadata)
points(x6G08_MD_360_proj[1],x6G08_MD_360_proj[3],pch=20)
text(x6G08_MD_360_proj[1],x6G08_MD_360_proj[3],pos=4,label="360",col="blue")

x6G08_MD_361_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_361.txt")
x6G08_MD_361_proj<-project.pca(x6G08_MD_361_raw,pcadata)
points(x6G08_MD_361_proj[1],x6G08_MD_361_proj[3],pch=20)
text(x6G08_MD_361_proj[1],x6G08_MD_361_proj[3],pos=4,label="361",col="red")

x6G08_MD_362_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_362.txt")
x6G08_MD_362_proj<-project.pca(x6G08_MD_362_raw,pcadata)
points(x6G08_MD_362_proj[1],x6G08_MD_362_proj[3],pch=20)
text(x6G08_MD_362_proj[1],x6G08_MD_362_proj[3],pos=4,label="362",col="red")

x6G08_MD_363_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_363.txt")
x6G08_MD_363_proj<-project.pca(x6G08_MD_363_raw,pcadata)
points(x6G08_MD_363_proj[1],x6G08_MD_363_proj[3],pch=20)
text(x6G08_MD_363_proj[1],x6G08_MD_363_proj[3],pos=4,label="363",col="red")

x6G08_MD_364_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_364.txt")
x6G08_MD_364_proj<-project.pca(x6G08_MD_364_raw,pcadata)
points(x6G08_MD_364_proj[1],x6G08_MD_364_proj[3],pch=20)
text(x6G08_MD_364_proj[1],x6G08_MD_364_proj[3],pos=4,label="364",col="red")

x6G08_MD_365_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_365.txt")
x6G08_MD_365_proj<-project.pca(x6G08_MD_365_raw,pcadata)
points(x6G08_MD_365_proj[1],x6G08_MD_365_proj[3],pch=20)
text(x6G08_MD_365_proj[1],x6G08_MD_365_proj[3],pos=4,label="365",col="red")

x6G08_MD_366_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_366.txt")
x6G08_MD_366_proj<-project.pca(x6G08_MD_366_raw,pcadata)
points(x6G08_MD_366_proj[1],x6G08_MD_366_proj[3],pch=20)
text(x6G08_MD_366_proj[1],x6G08_MD_366_proj[3],pos=4,label="366",col="red")

x6G08_MD_367_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_367.txt")
x6G08_MD_367_proj<-project.pca(x6G08_MD_367_raw,pcadata)
points(x6G08_MD_367_proj[1],x6G08_MD_367_proj[3],pch=20)
text(x6G08_MD_367_proj[1],x6G08_MD_367_proj[3],pos=4,label="367",col="green")

x6G08_MD_368_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_368.txt")
x6G08_MD_368_proj<-project.pca(x6G08_MD_368_raw,pcadata)
points(x6G08_MD_368_proj[1],x6G08_MD_368_proj[3],pch=20)
text(x6G08_MD_368_proj[1],x6G08_MD_368_proj[3],pos=4,label="368",col="green")

x6G08_MD_369_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_369.txt")
x6G08_MD_369_proj<-project.pca(x6G08_MD_369_raw,pcadata)
points(x6G08_MD_369_proj[1],x6G08_MD_369_proj[3],pch=20)
text(x6G08_MD_369_proj[1],x6G08_MD_369_proj[3],pos=4,label="369",col="red")

x6G08_MD_370_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_370.txt")
x6G08_MD_370_proj<-project.pca(x6G08_MD_370_raw,pcadata)
points(x6G08_MD_370_proj[1],x6G08_MD_370_proj[3],pch=20)
text(x6G08_MD_370_proj[1],x6G08_MD_370_proj[3],pos=4,label="370",col="green")

x6G08_MD_371_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_371.txt")
x6G08_MD_371_proj<-project.pca(x6G08_MD_371_raw,pcadata)
points(x6G08_MD_371_proj[1],x6G08_MD_371_proj[3],pch=20)
text(x6G08_MD_371_proj[1],x6G08_MD_371_proj[3],pos=4,label="371",col="blue")

x6G08_MD_372_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_372.txt")
x6G08_MD_372_proj<-project.pca(x6G08_MD_372_raw,pcadata)
points(x6G08_MD_372_proj[1],x6G08_MD_372_proj[3],pch=20)
text(x6G08_MD_372_proj[1],x6G08_MD_372_proj[3],pos=4,label="372",col="blue")

x6G08_MD_373_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_373.txt")
x6G08_MD_373_proj<-project.pca(x6G08_MD_373_raw,pcadata)
points(x6G08_MD_373_proj[1],x6G08_MD_373_proj[3],pch=20)
text(x6G08_MD_373_proj[1],x6G08_MD_373_proj[3],pos=4,label="373",col="blue")

x6G08_MD_374_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_374.txt")
x6G08_MD_374_proj<-project.pca(x6G08_MD_374_raw,pcadata)
points(x6G08_MD_374_proj[1],x6G08_MD_374_proj[3],pch=20)
text(x6G08_MD_374_proj[1],x6G08_MD_374_proj[3],pos=4,label="374",col="blue")

x6G08_MD_375_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_375.txt")
x6G08_MD_375_proj<-project.pca(x6G08_MD_375_raw,pcadata)
points(x6G08_MD_375_proj[1],x6G08_MD_375_proj[3],pch=20)
text(x6G08_MD_375_proj[1],x6G08_MD_375_proj[3],pos=4,label="375",col="red")

x6G08_MD_376_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_376.txt")
x6G08_MD_376_proj<-project.pca(x6G08_MD_376_raw,pcadata)
points(x6G08_MD_376_proj[1],x6G08_MD_376_proj[3],pch=20)
text(x6G08_MD_376_proj[1],x6G08_MD_376_proj[3],pos=4,label="376",col="blue")

x6G08_MD_377_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_377.txt")
x6G08_MD_377_proj<-project.pca(x6G08_MD_377_raw,pcadata)
points(x6G08_MD_377_proj[1],x6G08_MD_377_proj[3],pch=20)
text(x6G08_MD_377_proj[1],x6G08_MD_377_proj[3],pos=4,label="377",col="blue")

x6G08_MD_378_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_378.txt")
x6G08_MD_378_proj<-project.pca(x6G08_MD_378_raw,pcadata)
points(x6G08_MD_378_proj[1],x6G08_MD_378_proj[3],pch=20)
text(x6G08_MD_378_proj[1],x6G08_MD_378_proj[3],pos=4,label="378",col="black")

x6G08_MD_379_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_379.txt")
x6G08_MD_379_proj<-project.pca(x6G08_MD_379_raw,pcadata)
points(x6G08_MD_379_proj[1],x6G08_MD_379_proj[3],pch=20)
text(x6G08_MD_379_proj[1],x6G08_MD_379_proj[3],pos=4,label="379",col="blue")

x6G08_MD_380_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_380.txt")
x6G08_MD_380_proj<-project.pca(x6G08_MD_380_raw,pcadata)
points(x6G08_MD_380_proj[1],x6G08_MD_380_proj[3],pch=20)
text(x6G08_MD_380_proj[1],x6G08_MD_380_proj[3],pos=4,label="380",col="blue")

x6G08_MD_381_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_381.txt")
x6G08_MD_381_proj<-project.pca(x6G08_MD_381_raw,pcadata)
points(x6G08_MD_381_proj[1],x6G08_MD_381_proj[3],pch=20)
text(x6G08_MD_381_proj[1],x6G08_MD_381_proj[3],pos=4,label="381",col="black")

x6G08_MD_382_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_382.txt")
x6G08_MD_382_proj<-project.pca(x6G08_MD_382_raw,pcadata)
points(x6G08_MD_382_proj[1],x6G08_MD_382_proj[3],pch=20)
text(x6G08_MD_382_proj[1],x6G08_MD_382_proj[3],pos=4,label="382",col="blue")

x6G08_MD_383_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_383.txt")
x6G08_MD_383_proj<-project.pca(x6G08_MD_383_raw,pcadata)
points(x6G08_MD_383_proj[1],x6G08_MD_383_proj[3],pch=20)
text(x6G08_MD_383_proj[1],x6G08_MD_383_proj[3],pos=4,label="383",col="blue")

x6G08_MD_384_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_384.txt")
x6G08_MD_384_proj<-project.pca(x6G08_MD_384_raw,pcadata)
points(x6G08_MD_384_proj[1],x6G08_MD_384_proj[3],pch=20)
text(x6G08_MD_384_proj[1],x6G08_MD_384_proj[3],pos=4,label="384",col="black")

x6G08_MD_385_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_385.txt")
x6G08_MD_385_proj<-project.pca(x6G08_MD_385_raw,pcadata)
points(x6G08_MD_385_proj[1],x6G08_MD_385_proj[3],pch=20)
text(x6G08_MD_385_proj[1],x6G08_MD_385_proj[3],pos=4,label="385",col="blue")

x6G08_MD_386_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_386.txt")
x6G08_MD_386_proj<-project.pca(x6G08_MD_386_raw,pcadata)
points(x6G08_MD_386_proj[1],x6G08_MD_386_proj[3],pch=20)
text(x6G08_MD_386_proj[1],x6G08_MD_386_proj[3],pos=4,label="386",col="blue")

x6G08_MD_387_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_387.txt")
x6G08_MD_387_proj<-project.pca(x6G08_MD_387_raw,pcadata)
points(x6G08_MD_387_proj[1],x6G08_MD_387_proj[3],pch=20)
text(x6G08_MD_387_proj[1],x6G08_MD_387_proj[3],pos=4,label="387",col="blue")

x6G08_MD_388_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_388.txt")
x6G08_MD_388_proj<-project.pca(x6G08_MD_388_raw,pcadata)
points(x6G08_MD_388_proj[1],x6G08_MD_388_proj[3],pch=20)
text(x6G08_MD_388_proj[1],x6G08_MD_388_proj[3],pos=4,label="388",col="red")

x6G08_MD_389_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_389.txt")
x6G08_MD_389_proj<-project.pca(x6G08_MD_389_raw,pcadata)
points(x6G08_MD_389_proj[1],x6G08_MD_389_proj[3],pch=20)
text(x6G08_MD_389_proj[1],x6G08_MD_389_proj[3],pos=4,label="389",col="green")

x6G08_MD_390_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_390.txt")
x6G08_MD_390_proj<-project.pca(x6G08_MD_390_raw,pcadata)
points(x6G08_MD_390_proj[1],x6G08_MD_390_proj[3],pch=20)
text(x6G08_MD_390_proj[1],x6G08_MD_390_proj[3],pos=4,label="390",col="blue")

x6G08_MD_391_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_391.txt")
x6G08_MD_391_proj<-project.pca(x6G08_MD_391_raw,pcadata)
points(x6G08_MD_391_proj[1],x6G08_MD_391_proj[3],pch=20)
text(x6G08_MD_391_proj[1],x6G08_MD_391_proj[3],pos=4,label="391",col="green")

x6G08_MD_392_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_392.txt")
x6G08_MD_392_proj<-project.pca(x6G08_MD_392_raw,pcadata)
points(x6G08_MD_392_proj[1],x6G08_MD_392_proj[3],pch=20)
text(x6G08_MD_392_proj[1],x6G08_MD_392_proj[3],pos=4,label="392",col="green")

x6G08_MD_393_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_393.txt")
x6G08_MD_393_proj<-project.pca(x6G08_MD_393_raw,pcadata)
points(x6G08_MD_393_proj[1],x6G08_MD_393_proj[3],pch=20)
text(x6G08_MD_393_proj[1],x6G08_MD_393_proj[3],pos=4,label="393",col="green")

x6G08_MD_394_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_394.txt")
x6G08_MD_394_proj<-project.pca(x6G08_MD_394_raw,pcadata)
points(x6G08_MD_394_proj[1],x6G08_MD_394_proj[3],pch=20)
text(x6G08_MD_394_proj[1],x6G08_MD_394_proj[3],pos=4,label="394",col="green")

x6G08_MD_395_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_395.txt")
x6G08_MD_395_proj<-project.pca(x6G08_MD_395_raw,pcadata)
points(x6G08_MD_395_proj[1],x6G08_MD_395_proj[3],pch=20)
text(x6G08_MD_395_proj[1],x6G08_MD_395_proj[3],pos=4,label="395",col="green")

x6G08_MD_396_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_396.txt")
x6G08_MD_396_proj<-project.pca(x6G08_MD_396_raw,pcadata)
points(x6G08_MD_396_proj[1],x6G08_MD_396_proj[3],pch=20)
text(x6G08_MD_396_proj[1],x6G08_MD_396_proj[3],pos=4,label="396",col="green")

x6G08_MD_397_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_397.txt")
x6G08_MD_397_proj<-project.pca(x6G08_MD_397_raw,pcadata)
points(x6G08_MD_397_proj[1],x6G08_MD_397_proj[3],pch=20)
text(x6G08_MD_397_proj[1],x6G08_MD_397_proj[3],pos=4,label="397",col="blue")

x6G08_MD_398_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_398.txt")
x6G08_MD_398_proj<-project.pca(x6G08_MD_398_raw,pcadata)
points(x6G08_MD_398_proj[1],x6G08_MD_398_proj[3],pch=20)
text(x6G08_MD_398_proj[1],x6G08_MD_398_proj[3],pos=4,label="398",col="blue")

x6G08_MD_399_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_399.txt")
x6G08_MD_399_proj<-project.pca(x6G08_MD_399_raw,pcadata)
points(x6G08_MD_399_proj[1],x6G08_MD_399_proj[3],pch=20)
text(x6G08_MD_399_proj[1],x6G08_MD_399_proj[3],pos=4,label="399",col="green")

x6G08_MD_400_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_400.txt")
x6G08_MD_400_proj<-project.pca(x6G08_MD_400_raw,pcadata)
points(x6G08_MD_400_proj[1],x6G08_MD_400_proj[3],pch=20)
text(x6G08_MD_400_proj[1],x6G08_MD_400_proj[3],pos=4,label="400",col="blue")

x6G08_MD_401_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_401.txt")
x6G08_MD_401_proj<-project.pca(x6G08_MD_401_raw,pcadata)
points(x6G08_MD_401_proj[1],x6G08_MD_401_proj[3],pch=20)
text(x6G08_MD_401_proj[1],x6G08_MD_401_proj[3],pos=4,label="401",col="green")

x6G08_MD_402_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_402.txt")
x6G08_MD_402_proj<-project.pca(x6G08_MD_402_raw,pcadata)
points(x6G08_MD_402_proj[1],x6G08_MD_402_proj[3],pch=20)
text(x6G08_MD_402_proj[1],x6G08_MD_402_proj[3],pos=4,label="402",col="green")

x6G08_MD_403_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_403.txt")
x6G08_MD_403_proj<-project.pca(x6G08_MD_403_raw,pcadata)
points(x6G08_MD_403_proj[1],x6G08_MD_403_proj[3],pch=20)
text(x6G08_MD_403_proj[1],x6G08_MD_403_proj[3],pos=4,label="403",col="red")

x6G08_MD_404_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_404.txt")
x6G08_MD_404_proj<-project.pca(x6G08_MD_404_raw,pcadata)
points(x6G08_MD_404_proj[1],x6G08_MD_404_proj[3],pch=20)
text(x6G08_MD_404_proj[1],x6G08_MD_404_proj[3],pos=4,label="404",col="blue")

x6G08_MD_405_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_405.txt")
x6G08_MD_405_proj<-project.pca(x6G08_MD_405_raw,pcadata)
points(x6G08_MD_405_proj[1],x6G08_MD_405_proj[3],pch=20)
text(x6G08_MD_405_proj[1],x6G08_MD_405_proj[3],pos=4,label="405",col="red")

x6G08_MD_406_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_406.txt")
x6G08_MD_406_proj<-project.pca(x6G08_MD_406_raw,pcadata)
points(x6G08_MD_406_proj[1],x6G08_MD_406_proj[3],pch=20)
text(x6G08_MD_406_proj[1],x6G08_MD_406_proj[3],pos=4,label="406",col="red")

x6G08_MD_407_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_407.txt")
x6G08_MD_407_proj<-project.pca(x6G08_MD_407_raw,pcadata)
points(x6G08_MD_407_proj[1],x6G08_MD_407_proj[3],pch=20)
text(x6G08_MD_407_proj[1],x6G08_MD_407_proj[3],pos=4,label="407",col="red")

x6G08_MD_408_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_408.txt")
x6G08_MD_408_proj<-project.pca(x6G08_MD_408_raw,pcadata)
points(x6G08_MD_408_proj[1],x6G08_MD_408_proj[3],pch=20)
text(x6G08_MD_408_proj[1],x6G08_MD_408_proj[3],pos=4,label="408",col="red")

x6G08_MD_409_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_409.txt")
x6G08_MD_409_proj<-project.pca(x6G08_MD_409_raw,pcadata)
points(x6G08_MD_409_proj[1],x6G08_MD_409_proj[3],pch=20)
text(x6G08_MD_409_proj[1],x6G08_MD_409_proj[3],pos=4,label="409",col="red")

x6G08_MD_410_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_410.txt")
x6G08_MD_410_proj<-project.pca(x6G08_MD_410_raw,pcadata)
points(x6G08_MD_410_proj[1],x6G08_MD_410_proj[3],pch=20)
text(x6G08_MD_410_proj[1],x6G08_MD_410_proj[3],pos=4,label="410",col="green")

x6G08_MD_411_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_411.txt")
x6G08_MD_411_proj<-project.pca(x6G08_MD_411_raw,pcadata)
points(x6G08_MD_411_proj[1],x6G08_MD_411_proj[3],pch=20)
text(x6G08_MD_411_proj[1],x6G08_MD_411_proj[3],pos=4,label="411",col="red")

x6G08_MD_412_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_412.txt")
x6G08_MD_412_proj<-project.pca(x6G08_MD_412_raw,pcadata)
points(x6G08_MD_412_proj[1],x6G08_MD_412_proj[3],pch=20)
text(x6G08_MD_412_proj[1],x6G08_MD_412_proj[3],pos=4,label="412",col="red")

x6G08_MD_413_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_413.txt")
x6G08_MD_413_proj<-project.pca(x6G08_MD_413_raw,pcadata)
points(x6G08_MD_413_proj[1],x6G08_MD_413_proj[3],pch=20)
text(x6G08_MD_413_proj[1],x6G08_MD_413_proj[3],pos=4,label="413",col="green")

x6G08_MD_414_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_414.txt")
x6G08_MD_414_proj<-project.pca(x6G08_MD_414_raw,pcadata)
points(x6G08_MD_414_proj[1],x6G08_MD_414_proj[3],pch=20)
text(x6G08_MD_414_proj[1],x6G08_MD_414_proj[3],pos=4,label="414",col="green")

x6G08_MD_415_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_415.txt")
x6G08_MD_415_proj<-project.pca(x6G08_MD_415_raw,pcadata)
points(x6G08_MD_415_proj[1],x6G08_MD_415_proj[3],pch=20)
text(x6G08_MD_415_proj[1],x6G08_MD_415_proj[3],pos=4,label="415",col="green")

x6G08_MD_416_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_416.txt")
x6G08_MD_416_proj<-project.pca(x6G08_MD_416_raw,pcadata)
points(x6G08_MD_416_proj[1],x6G08_MD_416_proj[3],pch=20)
text(x6G08_MD_416_proj[1],x6G08_MD_416_proj[3],pos=4,label="416",col="red")

x6G08_MD_417_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_417.txt")
x6G08_MD_417_proj<-project.pca(x6G08_MD_417_raw,pcadata)
points(x6G08_MD_417_proj[1],x6G08_MD_417_proj[3],pch=20)
text(x6G08_MD_417_proj[1],x6G08_MD_417_proj[3],pos=4,label="417",col="red")

x6G08_MD_418_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_418.txt")
x6G08_MD_418_proj<-project.pca(x6G08_MD_418_raw,pcadata)
points(x6G08_MD_418_proj[1],x6G08_MD_418_proj[3],pch=20)
text(x6G08_MD_418_proj[1],x6G08_MD_418_proj[3],pos=4,label="418",col="red")

x6G08_MD_419_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_419.txt")
x6G08_MD_419_proj<-project.pca(x6G08_MD_419_raw,pcadata)
points(x6G08_MD_419_proj[1],x6G08_MD_419_proj[3],pch=20)
text(x6G08_MD_419_proj[1],x6G08_MD_419_proj[3],pos=4,label="419",col="green")

x6G08_MD_420_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_420.txt")
x6G08_MD_420_proj<-project.pca(x6G08_MD_420_raw,pcadata)
points(x6G08_MD_420_proj[1],x6G08_MD_420_proj[3],pch=20)
text(x6G08_MD_420_proj[1],x6G08_MD_420_proj[3],pos=4,label="420",col="red")

x6G08_MD_421_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_421.txt")
x6G08_MD_421_proj<-project.pca(x6G08_MD_421_raw,pcadata)
points(x6G08_MD_421_proj[1],x6G08_MD_421_proj[3],pch=20)
text(x6G08_MD_421_proj[1],x6G08_MD_421_proj[3],pos=4,label="421",col="red")

x6G08_MD_422_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_422.txt")
x6G08_MD_422_proj<-project.pca(x6G08_MD_422_raw,pcadata)
points(x6G08_MD_422_proj[1],x6G08_MD_422_proj[3],pch=20)
text(x6G08_MD_422_proj[1],x6G08_MD_422_proj[3],pos=4,label="422",col="red")

x6G08_MD_423_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_423.txt")
x6G08_MD_423_proj<-project.pca(x6G08_MD_423_raw,pcadata)
points(x6G08_MD_423_proj[1],x6G08_MD_423_proj[3],pch=20)
text(x6G08_MD_423_proj[1],x6G08_MD_423_proj[3],pos=4,label="423",col="blue")

x6G08_MD_424_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_424.txt")
x6G08_MD_424_proj<-project.pca(x6G08_MD_424_raw,pcadata)
points(x6G08_MD_424_proj[1],x6G08_MD_424_proj[3],pch=20)
text(x6G08_MD_424_proj[1],x6G08_MD_424_proj[3],pos=4,label="424",col="red")

x6G08_MD_425_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_425.txt")
x6G08_MD_425_proj<-project.pca(x6G08_MD_425_raw,pcadata)
points(x6G08_MD_425_proj[1],x6G08_MD_425_proj[3],pch=20)
text(x6G08_MD_425_proj[1],x6G08_MD_425_proj[3],pos=4,label="425",col="red")

x6G08_MD_426_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_426.txt")
x6G08_MD_426_proj<-project.pca(x6G08_MD_426_raw,pcadata)
points(x6G08_MD_426_proj[1],x6G08_MD_426_proj[3],pch=20)
text(x6G08_MD_426_proj[1],x6G08_MD_426_proj[3],pos=4,label="426",col="red")

x6G08_MD_427_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_427.txt")
x6G08_MD_427_proj<-project.pca(x6G08_MD_427_raw,pcadata)
points(x6G08_MD_427_proj[1],x6G08_MD_427_proj[3],pch=20)
text(x6G08_MD_427_proj[1],x6G08_MD_427_proj[3],pos=4,label="427",col="green")

x6G08_MD_428_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_428.txt")
x6G08_MD_428_proj<-project.pca(x6G08_MD_428_raw,pcadata)
points(x6G08_MD_428_proj[1],x6G08_MD_428_proj[3],pch=20)
text(x6G08_MD_428_proj[1],x6G08_MD_428_proj[3],pos=4,label="428",col="red")

x6G08_MD_429_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_429.txt")
x6G08_MD_429_proj<-project.pca(x6G08_MD_429_raw,pcadata)
points(x6G08_MD_429_proj[1],x6G08_MD_429_proj[3],pch=20)
text(x6G08_MD_429_proj[1],x6G08_MD_429_proj[3],pos=4,label="429",col="green")

x6G08_MD_430_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_430.txt")
x6G08_MD_430_proj<-project.pca(x6G08_MD_430_raw,pcadata)
points(x6G08_MD_430_proj[1],x6G08_MD_430_proj[3],pch=20)
text(x6G08_MD_430_proj[1],x6G08_MD_430_proj[3],pos=4,label="430",col="red")

x6G08_MD_431_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_431.txt")
x6G08_MD_431_proj<-project.pca(x6G08_MD_431_raw,pcadata)
points(x6G08_MD_431_proj[1],x6G08_MD_431_proj[3],pch=20)
text(x6G08_MD_431_proj[1],x6G08_MD_431_proj[3],pos=4,label="431",col="green")

x6G08_MD_432_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_432.txt")
x6G08_MD_432_proj<-project.pca(x6G08_MD_432_raw,pcadata)
points(x6G08_MD_432_proj[1],x6G08_MD_432_proj[3],pch=20)
text(x6G08_MD_432_proj[1],x6G08_MD_432_proj[3],pos=4,label="432",col="green")

x6G08_MD_433_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_433.txt")
x6G08_MD_433_proj<-project.pca(x6G08_MD_433_raw,pcadata)
points(x6G08_MD_433_proj[1],x6G08_MD_433_proj[3],pch=20)
text(x6G08_MD_433_proj[1],x6G08_MD_433_proj[3],pos=4,label="433",col="green")

x6G08_MD_434_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_434.txt")
x6G08_MD_434_proj<-project.pca(x6G08_MD_434_raw,pcadata)
points(x6G08_MD_434_proj[1],x6G08_MD_434_proj[3],pch=20)
text(x6G08_MD_434_proj[1],x6G08_MD_434_proj[3],pos=4,label="434",col="red")

x6G08_MD_435_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_435.txt")
x6G08_MD_435_proj<-project.pca(x6G08_MD_435_raw,pcadata)
points(x6G08_MD_435_proj[1],x6G08_MD_435_proj[3],pch=20)
text(x6G08_MD_435_proj[1],x6G08_MD_435_proj[3],pos=4,label="435",col="green")

x6G08_MD_436_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_436.txt")
x6G08_MD_436_proj<-project.pca(x6G08_MD_436_raw,pcadata)
points(x6G08_MD_436_proj[1],x6G08_MD_436_proj[3],pch=20)
text(x6G08_MD_436_proj[1],x6G08_MD_436_proj[3],pos=4,label="436",col="red")

x6G08_MD_437_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_437.txt")
x6G08_MD_437_proj<-project.pca(x6G08_MD_437_raw,pcadata)
points(x6G08_MD_437_proj[1],x6G08_MD_437_proj[3],pch=20)
text(x6G08_MD_437_proj[1],x6G08_MD_437_proj[3],pos=4,label="437",col="green")

x6G08_MD_438_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_438.txt")
x6G08_MD_438_proj<-project.pca(x6G08_MD_438_raw,pcadata)
points(x6G08_MD_438_proj[1],x6G08_MD_438_proj[3],pch=20)
text(x6G08_MD_438_proj[1],x6G08_MD_438_proj[3],pos=4,label="438",col="red")

x6G08_MD_439_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_439.txt")
x6G08_MD_439_proj<-project.pca(x6G08_MD_439_raw,pcadata)
points(x6G08_MD_439_proj[1],x6G08_MD_439_proj[3],pch=20)
text(x6G08_MD_439_proj[1],x6G08_MD_439_proj[3],pos=4,label="439",col="red")

x6G08_MD_440_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_440.txt")
x6G08_MD_440_proj<-project.pca(x6G08_MD_440_raw,pcadata)
points(x6G08_MD_440_proj[1],x6G08_MD_440_proj[3],pch=20)
text(x6G08_MD_440_proj[1],x6G08_MD_440_proj[3],pos=4,label="440",col="red")

x6G08_MD_441_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_441.txt")
x6G08_MD_441_proj<-project.pca(x6G08_MD_441_raw,pcadata)
points(x6G08_MD_441_proj[1],x6G08_MD_441_proj[3],pch=20)
text(x6G08_MD_441_proj[1],x6G08_MD_441_proj[3],pos=4,label="441",col="red")

x6G08_MD_442_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_442.txt")
x6G08_MD_442_proj<-project.pca(x6G08_MD_442_raw,pcadata)
points(x6G08_MD_442_proj[1],x6G08_MD_442_proj[3],pch=20)
text(x6G08_MD_442_proj[1],x6G08_MD_442_proj[3],pos=4,label="442",col="red")

x6G08_MD_443_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_443.txt")
x6G08_MD_443_proj<-project.pca(x6G08_MD_443_raw,pcadata)
points(x6G08_MD_443_proj[1],x6G08_MD_443_proj[3],pch=20)
text(x6G08_MD_443_proj[1],x6G08_MD_443_proj[3],pos=4,label="443",col="red")

x6G08_MD_444_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_444.txt")
x6G08_MD_444_proj<-project.pca(x6G08_MD_444_raw,pcadata)
points(x6G08_MD_444_proj[1],x6G08_MD_444_proj[3],pch=20)
text(x6G08_MD_444_proj[1],x6G08_MD_444_proj[3],pos=4,label="444",col="red")

x6G08_MD_445_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_445.txt")
x6G08_MD_445_proj<-project.pca(x6G08_MD_445_raw,pcadata)
points(x6G08_MD_445_proj[1],x6G08_MD_445_proj[3],pch=20)
text(x6G08_MD_445_proj[1],x6G08_MD_445_proj[3],pos=4,label="445",col="red")

x6G08_MD_446_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_446.txt")
x6G08_MD_446_proj<-project.pca(x6G08_MD_446_raw,pcadata)
points(x6G08_MD_446_proj[1],x6G08_MD_446_proj[3],pch=20)
text(x6G08_MD_446_proj[1],x6G08_MD_446_proj[3],pos=4,label="446",col="red")

x6G08_MD_447_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_447.txt")
x6G08_MD_447_proj<-project.pca(x6G08_MD_447_raw,pcadata)
points(x6G08_MD_447_proj[1],x6G08_MD_447_proj[3],pch=20)
text(x6G08_MD_447_proj[1],x6G08_MD_447_proj[3],pos=4,label="447",col="red")

x6G08_MD_448_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_448.txt")
x6G08_MD_448_proj<-project.pca(x6G08_MD_448_raw,pcadata)
points(x6G08_MD_448_proj[1],x6G08_MD_448_proj[3],pch=20)
text(x6G08_MD_448_proj[1],x6G08_MD_448_proj[3],pos=4,label="448",col="red")

x6G08_MD_449_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_449.txt")
x6G08_MD_449_proj<-project.pca(x6G08_MD_449_raw,pcadata)
points(x6G08_MD_449_proj[1],x6G08_MD_449_proj[3],pch=20)
text(x6G08_MD_449_proj[1],x6G08_MD_449_proj[3],pos=4,label="449",col="red")

x6G08_MD_450_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_450.txt")
x6G08_MD_450_proj<-project.pca(x6G08_MD_450_raw,pcadata)
points(x6G08_MD_450_proj[1],x6G08_MD_450_proj[3],pch=20)
text(x6G08_MD_450_proj[1],x6G08_MD_450_proj[3],pos=4,label="450",col="red")

x6G08_MD_451_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_451.txt")
x6G08_MD_451_proj<-project.pca(x6G08_MD_451_raw,pcadata)
points(x6G08_MD_451_proj[1],x6G08_MD_451_proj[3],pch=20)
text(x6G08_MD_451_proj[1],x6G08_MD_451_proj[3],pos=4,label="451",col="red")

x6G08_MD_452_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_452.txt")
x6G08_MD_452_proj<-project.pca(x6G08_MD_452_raw,pcadata)
points(x6G08_MD_452_proj[1],x6G08_MD_452_proj[3],pch=20)
text(x6G08_MD_452_proj[1],x6G08_MD_452_proj[3],pos=4,label="452",col="red")

x6G08_MD_453_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_453.txt")
x6G08_MD_453_proj<-project.pca(x6G08_MD_453_raw,pcadata)
points(x6G08_MD_453_proj[1],x6G08_MD_453_proj[3],pch=20)
text(x6G08_MD_453_proj[1],x6G08_MD_453_proj[3],pos=4,label="453",col="red")

x6G08_MD_454_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_454.txt")
x6G08_MD_454_proj<-project.pca(x6G08_MD_454_raw,pcadata)
points(x6G08_MD_454_proj[1],x6G08_MD_454_proj[3],pch=20)
text(x6G08_MD_454_proj[1],x6G08_MD_454_proj[3],pos=4,label="454",col="green")

x6G08_MD_455_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_455.txt")
x6G08_MD_455_proj<-project.pca(x6G08_MD_455_raw,pcadata)
points(x6G08_MD_455_proj[1],x6G08_MD_455_proj[3],pch=20)
text(x6G08_MD_455_proj[1],x6G08_MD_455_proj[3],pos=4,label="455",col="green")

x6G08_MD_456_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_456.txt")
x6G08_MD_456_proj<-project.pca(x6G08_MD_456_raw,pcadata)
points(x6G08_MD_456_proj[1],x6G08_MD_456_proj[3],pch=20)
text(x6G08_MD_456_proj[1],x6G08_MD_456_proj[3],pos=4,label="456",col="red")

x6G08_MD_457_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_457.txt")
x6G08_MD_457_proj<-project.pca(x6G08_MD_457_raw,pcadata)
points(x6G08_MD_457_proj[1],x6G08_MD_457_proj[3],pch=20)
text(x6G08_MD_457_proj[1],x6G08_MD_457_proj[3],pos=4,label="457",col="red")

x6G08_MD_458_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_458.txt")
x6G08_MD_458_proj<-project.pca(x6G08_MD_458_raw,pcadata)
points(x6G08_MD_458_proj[1],x6G08_MD_458_proj[3],pch=20)
text(x6G08_MD_458_proj[1],x6G08_MD_458_proj[3],pos=4,label="458",col="red")

x6G08_MD_459_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_459.txt")
x6G08_MD_459_proj<-project.pca(x6G08_MD_459_raw,pcadata)
points(x6G08_MD_459_proj[1],x6G08_MD_459_proj[3],pch=20)
text(x6G08_MD_459_proj[1],x6G08_MD_459_proj[3],pos=4,label="459",col="red")

x6G08_MD_460_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_460.txt")
x6G08_MD_460_proj<-project.pca(x6G08_MD_460_raw,pcadata)
points(x6G08_MD_460_proj[1],x6G08_MD_460_proj[3],pch=20)
text(x6G08_MD_460_proj[1],x6G08_MD_460_proj[3],pos=4,label="460",col="blue")

x6G08_MD_461_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_461.txt")
x6G08_MD_461_proj<-project.pca(x6G08_MD_461_raw,pcadata)
points(x6G08_MD_461_proj[1],x6G08_MD_461_proj[3],pch=20)
text(x6G08_MD_461_proj[1],x6G08_MD_461_proj[3],pos=4,label="461",col="blue")

x6G08_MD_462_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_462.txt")
x6G08_MD_462_proj<-project.pca(x6G08_MD_462_raw,pcadata)
points(x6G08_MD_462_proj[1],x6G08_MD_462_proj[3],pch=20)
text(x6G08_MD_462_proj[1],x6G08_MD_462_proj[3],pos=4,label="462",col="green")

x6G08_MD_463_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_463.txt")
x6G08_MD_463_proj<-project.pca(x6G08_MD_463_raw,pcadata)
points(x6G08_MD_463_proj[1],x6G08_MD_463_proj[3],pch=20)
text(x6G08_MD_463_proj[1],x6G08_MD_463_proj[3],pos=4,label="463",col="blue")

x6G08_MD_464_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_464.txt")
x6G08_MD_464_proj<-project.pca(x6G08_MD_464_raw,pcadata)
points(x6G08_MD_464_proj[1],x6G08_MD_464_proj[3],pch=20)
text(x6G08_MD_464_proj[1],x6G08_MD_464_proj[3],pos=4,label="464",col="green")

x6G08_MD_465_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_465.txt")
x6G08_MD_465_proj<-project.pca(x6G08_MD_465_raw,pcadata)
points(x6G08_MD_465_proj[1],x6G08_MD_465_proj[3],pch=20)
text(x6G08_MD_465_proj[1],x6G08_MD_465_proj[3],pos=4,label="465",col="green")

x6G08_MD_466_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_466.txt")
x6G08_MD_466_proj<-project.pca(x6G08_MD_466_raw,pcadata)
points(x6G08_MD_466_proj[1],x6G08_MD_466_proj[3],pch=20)
text(x6G08_MD_466_proj[1],x6G08_MD_466_proj[3],pos=4,label="466",col="blue")

x6G08_MD_467_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_467.txt")
x6G08_MD_467_proj<-project.pca(x6G08_MD_467_raw,pcadata)
points(x6G08_MD_467_proj[1],x6G08_MD_467_proj[3],pch=20)
text(x6G08_MD_467_proj[1],x6G08_MD_467_proj[3],pos=4,label="467",col="blue")

x6G08_MD_468_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_468.txt")
x6G08_MD_468_proj<-project.pca(x6G08_MD_468_raw,pcadata)
points(x6G08_MD_468_proj[1],x6G08_MD_468_proj[3],pch=20)
text(x6G08_MD_468_proj[1],x6G08_MD_468_proj[3],pos=4,label="468",col="blue")

x6G08_MD_469_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_469.txt")
x6G08_MD_469_proj<-project.pca(x6G08_MD_469_raw,pcadata)
points(x6G08_MD_469_proj[1],x6G08_MD_469_proj[3],pch=20)
text(x6G08_MD_469_proj[1],x6G08_MD_469_proj[3],pos=4,label="469",col="black")

x6G08_MD_470_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_470.txt")
x6G08_MD_470_proj<-project.pca(x6G08_MD_470_raw,pcadata)
points(x6G08_MD_470_proj[1],x6G08_MD_470_proj[3],pch=20)
text(x6G08_MD_470_proj[1],x6G08_MD_470_proj[3],pos=4,label="470",col="blue")

x6G08_MD_471_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_471.txt")
x6G08_MD_471_proj<-project.pca(x6G08_MD_471_raw,pcadata)
points(x6G08_MD_471_proj[1],x6G08_MD_471_proj[3],pch=20)
text(x6G08_MD_471_proj[1],x6G08_MD_471_proj[3],pos=4,label="471",col="green")

x6G08_MD_472_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_472.txt")
x6G08_MD_472_proj<-project.pca(x6G08_MD_472_raw,pcadata)
points(x6G08_MD_472_proj[1],x6G08_MD_472_proj[3],pch=20)
text(x6G08_MD_472_proj[1],x6G08_MD_472_proj[3],pos=4,label="472",col="blue")

x6G08_MD_473_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_473.txt")
x6G08_MD_473_proj<-project.pca(x6G08_MD_473_raw,pcadata)
points(x6G08_MD_473_proj[1],x6G08_MD_473_proj[3],pch=20)
text(x6G08_MD_473_proj[1],x6G08_MD_473_proj[3],pos=4,label="473",col="red")

x6G08_MD_474_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_474.txt")
x6G08_MD_474_proj<-project.pca(x6G08_MD_474_raw,pcadata)
points(x6G08_MD_474_proj[1],x6G08_MD_474_proj[3],pch=20)
text(x6G08_MD_474_proj[1],x6G08_MD_474_proj[3],pos=4,label="474",col="green")

x6G08_MD_475_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_475.txt")
x6G08_MD_475_proj<-project.pca(x6G08_MD_475_raw,pcadata)
points(x6G08_MD_475_proj[1],x6G08_MD_475_proj[3],pch=20)
text(x6G08_MD_475_proj[1],x6G08_MD_475_proj[3],pos=4,label="475",col="green")

x6G08_MD_476_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_476.txt")
x6G08_MD_476_proj<-project.pca(x6G08_MD_476_raw,pcadata)
points(x6G08_MD_476_proj[1],x6G08_MD_476_proj[3],pch=20)
text(x6G08_MD_476_proj[1],x6G08_MD_476_proj[3],pos=4,label="476",col="red")

x6G08_MD_477_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_477.txt")
x6G08_MD_477_proj<-project.pca(x6G08_MD_477_raw,pcadata)
points(x6G08_MD_477_proj[1],x6G08_MD_477_proj[3],pch=20)
text(x6G08_MD_477_proj[1],x6G08_MD_477_proj[3],pos=4,label="477",col="green")

x6G08_MD_478_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_478.txt")
x6G08_MD_478_proj<-project.pca(x6G08_MD_478_raw,pcadata)
points(x6G08_MD_478_proj[1],x6G08_MD_478_proj[3],pch=20)
text(x6G08_MD_478_proj[1],x6G08_MD_478_proj[3],pos=4,label="478",col="blue")

x6G08_MD_479_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_479.txt")
x6G08_MD_479_proj<-project.pca(x6G08_MD_479_raw,pcadata)
points(x6G08_MD_479_proj[1],x6G08_MD_479_proj[3],pch=20)
text(x6G08_MD_479_proj[1],x6G08_MD_479_proj[3],pos=4,label="479",col="black")

x6G08_MD_480_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_480.txt")
x6G08_MD_480_proj<-project.pca(x6G08_MD_480_raw,pcadata)
points(x6G08_MD_480_proj[1],x6G08_MD_480_proj[3],pch=20)
text(x6G08_MD_480_proj[1],x6G08_MD_480_proj[3],pos=4,label="480",col="blue")

x6G08_MD_481_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_481.txt")
x6G08_MD_481_proj<-project.pca(x6G08_MD_481_raw,pcadata)
points(x6G08_MD_481_proj[1],x6G08_MD_481_proj[3],pch=20)
text(x6G08_MD_481_proj[1],x6G08_MD_481_proj[3],pos=4,label="481",col="blue")

x6G08_MD_482_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_482.txt")
x6G08_MD_482_proj<-project.pca(x6G08_MD_482_raw,pcadata)
points(x6G08_MD_482_proj[1],x6G08_MD_482_proj[3],pch=20)
text(x6G08_MD_482_proj[1],x6G08_MD_482_proj[3],pos=4,label="482",col="black")

x6G08_MD_483_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_483.txt")
x6G08_MD_483_proj<-project.pca(x6G08_MD_483_raw,pcadata)
points(x6G08_MD_483_proj[1],x6G08_MD_483_proj[3],pch=20)
text(x6G08_MD_483_proj[1],x6G08_MD_483_proj[3],pos=4,label="483",col="blue")

x6G08_MD_484_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_484.txt")
x6G08_MD_484_proj<-project.pca(x6G08_MD_484_raw,pcadata)
points(x6G08_MD_484_proj[1],x6G08_MD_484_proj[3],pch=20)
text(x6G08_MD_484_proj[1],x6G08_MD_484_proj[3],pos=4,label="484",col="black")

x6G08_MD_485_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_485.txt")
x6G08_MD_485_proj<-project.pca(x6G08_MD_485_raw,pcadata)
points(x6G08_MD_485_proj[1],x6G08_MD_485_proj[3],pch=20)
text(x6G08_MD_485_proj[1],x6G08_MD_485_proj[3],pos=4,label="485",col="blue")

x6G08_MD_486_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_486.txt")
x6G08_MD_486_proj<-project.pca(x6G08_MD_486_raw,pcadata)
points(x6G08_MD_486_proj[1],x6G08_MD_486_proj[3],pch=20)
text(x6G08_MD_486_proj[1],x6G08_MD_486_proj[3],pos=4,label="486",col="black")

x6G08_MD_487_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_487.txt")
x6G08_MD_487_proj<-project.pca(x6G08_MD_487_raw,pcadata)
points(x6G08_MD_487_proj[1],x6G08_MD_487_proj[3],pch=20)
text(x6G08_MD_487_proj[1],x6G08_MD_487_proj[3],pos=4,label="487",col="blue")

x6G08_MD_488_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_488.txt")
x6G08_MD_488_proj<-project.pca(x6G08_MD_488_raw,pcadata)
points(x6G08_MD_488_proj[1],x6G08_MD_488_proj[3],pch=20)
text(x6G08_MD_488_proj[1],x6G08_MD_488_proj[3],pos=4,label="488",col="blue")

x6G08_MD_489_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_489.txt")
x6G08_MD_489_proj<-project.pca(x6G08_MD_489_raw,pcadata)
points(x6G08_MD_489_proj[1],x6G08_MD_489_proj[3],pch=20)
text(x6G08_MD_489_proj[1],x6G08_MD_489_proj[3],pos=4,label="489",col="blue")

x6G08_MD_490_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_490.txt")
x6G08_MD_490_proj<-project.pca(x6G08_MD_490_raw,pcadata)
points(x6G08_MD_490_proj[1],x6G08_MD_490_proj[3],pch=20)
text(x6G08_MD_490_proj[1],x6G08_MD_490_proj[3],pos=4,label="490",col="blue")

x6G08_MD_491_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_491.txt")
x6G08_MD_491_proj<-project.pca(x6G08_MD_491_raw,pcadata)
points(x6G08_MD_491_proj[1],x6G08_MD_491_proj[3],pch=20)
text(x6G08_MD_491_proj[1],x6G08_MD_491_proj[3],pos=4,label="491",col="blue")

x6G08_MD_492_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_492.txt")
x6G08_MD_492_proj<-project.pca(x6G08_MD_492_raw,pcadata)
points(x6G08_MD_492_proj[1],x6G08_MD_492_proj[3],pch=20)
text(x6G08_MD_492_proj[1],x6G08_MD_492_proj[3],pos=4,label="492",col="blue")

x6G08_MD_493_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_493.txt")
x6G08_MD_493_proj<-project.pca(x6G08_MD_493_raw,pcadata)
points(x6G08_MD_493_proj[1],x6G08_MD_493_proj[3],pch=20)
text(x6G08_MD_493_proj[1],x6G08_MD_493_proj[3],pos=4,label="493",col="green")

x6G08_MD_494_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_494.txt")
x6G08_MD_494_proj<-project.pca(x6G08_MD_494_raw,pcadata)
points(x6G08_MD_494_proj[1],x6G08_MD_494_proj[3],pch=20)
text(x6G08_MD_494_proj[1],x6G08_MD_494_proj[3],pos=4,label="494",col="green")

x6G08_MD_495_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_495.txt")
x6G08_MD_495_proj<-project.pca(x6G08_MD_495_raw,pcadata)
points(x6G08_MD_495_proj[1],x6G08_MD_495_proj[3],pch=20)
text(x6G08_MD_495_proj[1],x6G08_MD_495_proj[3],pos=4,label="495",col="red")

x6G08_MD_496_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_496.txt")
x6G08_MD_496_proj<-project.pca(x6G08_MD_496_raw,pcadata)
points(x6G08_MD_496_proj[1],x6G08_MD_496_proj[3],pch=20)
text(x6G08_MD_496_proj[1],x6G08_MD_496_proj[3],pos=4,label="496",col="green")

x6G08_MD_497_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_497.txt")
x6G08_MD_497_proj<-project.pca(x6G08_MD_497_raw,pcadata)
points(x6G08_MD_497_proj[1],x6G08_MD_497_proj[3],pch=20)
text(x6G08_MD_497_proj[1],x6G08_MD_497_proj[3],pos=4,label="497",col="green")

x6G08_MD_498_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_498.txt")
x6G08_MD_498_proj<-project.pca(x6G08_MD_498_raw,pcadata)
points(x6G08_MD_498_proj[1],x6G08_MD_498_proj[3],pch=20)
text(x6G08_MD_498_proj[1],x6G08_MD_498_proj[3],pos=4,label="498",col="black")

x6G08_MD_499_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_499.txt")
x6G08_MD_499_proj<-project.pca(x6G08_MD_499_raw,pcadata)
points(x6G08_MD_499_proj[1],x6G08_MD_499_proj[3],pch=20)
text(x6G08_MD_499_proj[1],x6G08_MD_499_proj[3],pos=4,label="499",col="green")

x6G08_MD_500_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_500.txt")
x6G08_MD_500_proj<-project.pca(x6G08_MD_500_raw,pcadata)
points(x6G08_MD_500_proj[1],x6G08_MD_500_proj[3],pch=20)
text(x6G08_MD_500_proj[1],x6G08_MD_500_proj[3],pos=4,label="500",col="red")

x6G08_MD_501_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_501.txt")
x6G08_MD_501_proj<-project.pca(x6G08_MD_501_raw,pcadata)
points(x6G08_MD_501_proj[1],x6G08_MD_501_proj[3],pch=20)
text(x6G08_MD_501_proj[1],x6G08_MD_501_proj[3],pos=4,label="501",col="red")

x6G08_MD_502_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_502.txt")
x6G08_MD_502_proj<-project.pca(x6G08_MD_502_raw,pcadata)
points(x6G08_MD_502_proj[1],x6G08_MD_502_proj[3],pch=20)
text(x6G08_MD_502_proj[1],x6G08_MD_502_proj[3],pos=4,label="502",col="blue")

x6G08_MD_503_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_503.txt")
x6G08_MD_503_proj<-project.pca(x6G08_MD_503_raw,pcadata)
points(x6G08_MD_503_proj[1],x6G08_MD_503_proj[3],pch=20)
text(x6G08_MD_503_proj[1],x6G08_MD_503_proj[3],pos=4,label="503",col="black")

x6G08_MD_504_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_504.txt")
x6G08_MD_504_proj<-project.pca(x6G08_MD_504_raw,pcadata)
points(x6G08_MD_504_proj[1],x6G08_MD_504_proj[3],pch=20)
text(x6G08_MD_504_proj[1],x6G08_MD_504_proj[3],pos=4,label="504",col="blue")

x6G08_MD_505_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_505.txt")
x6G08_MD_505_proj<-project.pca(x6G08_MD_505_raw,pcadata)
points(x6G08_MD_505_proj[1],x6G08_MD_505_proj[3],pch=20)
text(x6G08_MD_505_proj[1],x6G08_MD_505_proj[3],pos=4,label="505",col="green")

x6G08_MD_506_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_506.txt")
x6G08_MD_506_proj<-project.pca(x6G08_MD_506_raw,pcadata)
points(x6G08_MD_506_proj[1],x6G08_MD_506_proj[3],pch=20)
text(x6G08_MD_506_proj[1],x6G08_MD_506_proj[3],pos=4,label="506",col="blue")

x6G08_MD_507_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_507.txt")
x6G08_MD_507_proj<-project.pca(x6G08_MD_507_raw,pcadata)
points(x6G08_MD_507_proj[1],x6G08_MD_507_proj[3],pch=20)
text(x6G08_MD_507_proj[1],x6G08_MD_507_proj[3],pos=4,label="507",col="blue")

x6G08_MD_508_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_508.txt")
x6G08_MD_508_proj<-project.pca(x6G08_MD_508_raw,pcadata)
points(x6G08_MD_508_proj[1],x6G08_MD_508_proj[3],pch=20)
text(x6G08_MD_508_proj[1],x6G08_MD_508_proj[3],pos=4,label="508",col="blue")

x6G08_MD_509_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_509.txt")
x6G08_MD_509_proj<-project.pca(x6G08_MD_509_raw,pcadata)
points(x6G08_MD_509_proj[1],x6G08_MD_509_proj[3],pch=20)
text(x6G08_MD_509_proj[1],x6G08_MD_509_proj[3],pos=4,label="509",col="blue")

x6G08_MD_510_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_510.txt")
x6G08_MD_510_proj<-project.pca(x6G08_MD_510_raw,pcadata)
points(x6G08_MD_510_proj[1],x6G08_MD_510_proj[3],pch=20)
text(x6G08_MD_510_proj[1],x6G08_MD_510_proj[3],pos=4,label="510",col="blue")

x6G08_MD_511_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_511.txt")
x6G08_MD_511_proj<-project.pca(x6G08_MD_511_raw,pcadata)
points(x6G08_MD_511_proj[1],x6G08_MD_511_proj[3],pch=20)
text(x6G08_MD_511_proj[1],x6G08_MD_511_proj[3],pos=4,label="511",col="blue")

x6G08_MD_512_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_512.txt")
x6G08_MD_512_proj<-project.pca(x6G08_MD_512_raw,pcadata)
points(x6G08_MD_512_proj[1],x6G08_MD_512_proj[3],pch=20)
text(x6G08_MD_512_proj[1],x6G08_MD_512_proj[3],pos=4,label="512",col="blue")

x6G08_MD_513_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_513.txt")
x6G08_MD_513_proj<-project.pca(x6G08_MD_513_raw,pcadata)
points(x6G08_MD_513_proj[1],x6G08_MD_513_proj[3],pch=20)
text(x6G08_MD_513_proj[1],x6G08_MD_513_proj[3],pos=4,label="513",col="black")

x6G08_MD_514_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_514.txt")
x6G08_MD_514_proj<-project.pca(x6G08_MD_514_raw,pcadata)
points(x6G08_MD_514_proj[1],x6G08_MD_514_proj[3],pch=20)
text(x6G08_MD_514_proj[1],x6G08_MD_514_proj[3],pos=4,label="514",col="black")

x6G08_MD_515_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_515.txt")
x6G08_MD_515_proj<-project.pca(x6G08_MD_515_raw,pcadata)
points(x6G08_MD_515_proj[1],x6G08_MD_515_proj[3],pch=20)
text(x6G08_MD_515_proj[1],x6G08_MD_515_proj[3],pos=4,label="515",col="blue")

x6G08_MD_516_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_516.txt")
x6G08_MD_516_proj<-project.pca(x6G08_MD_516_raw,pcadata)
points(x6G08_MD_516_proj[1],x6G08_MD_516_proj[3],pch=20)
text(x6G08_MD_516_proj[1],x6G08_MD_516_proj[3],pos=4,label="516",col="black")

x6G08_MD_517_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_517.txt")
x6G08_MD_517_proj<-project.pca(x6G08_MD_517_raw,pcadata)
points(x6G08_MD_517_proj[1],x6G08_MD_517_proj[3],pch=20)
text(x6G08_MD_517_proj[1],x6G08_MD_517_proj[3],pos=4,label="517",col="black")

x6G08_MD_518_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_518.txt")
x6G08_MD_518_proj<-project.pca(x6G08_MD_518_raw,pcadata)
points(x6G08_MD_518_proj[1],x6G08_MD_518_proj[3],pch=20)
text(x6G08_MD_518_proj[1],x6G08_MD_518_proj[3],pos=4,label="518",col="black")

x6G08_MD_519_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_519.txt")
x6G08_MD_519_proj<-project.pca(x6G08_MD_519_raw,pcadata)
points(x6G08_MD_519_proj[1],x6G08_MD_519_proj[3],pch=20)
text(x6G08_MD_519_proj[1],x6G08_MD_519_proj[3],pos=4,label="519",col="black")

x6G08_MD_520_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_520.txt")
x6G08_MD_520_proj<-project.pca(x6G08_MD_520_raw,pcadata)
points(x6G08_MD_520_proj[1],x6G08_MD_520_proj[3],pch=20)
text(x6G08_MD_520_proj[1],x6G08_MD_520_proj[3],pos=4,label="520",col="blue")

x6G08_MD_521_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_521.txt")
x6G08_MD_521_proj<-project.pca(x6G08_MD_521_raw,pcadata)
points(x6G08_MD_521_proj[1],x6G08_MD_521_proj[3],pch=20)
text(x6G08_MD_521_proj[1],x6G08_MD_521_proj[3],pos=4,label="521",col="blue")

x6G08_MD_522_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_522.txt")
x6G08_MD_522_proj<-project.pca(x6G08_MD_522_raw,pcadata)
points(x6G08_MD_522_proj[1],x6G08_MD_522_proj[3],pch=20)
text(x6G08_MD_522_proj[1],x6G08_MD_522_proj[3],pos=4,label="522",col="black")

x6G08_MD_523_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_523.txt")
x6G08_MD_523_proj<-project.pca(x6G08_MD_523_raw,pcadata)
points(x6G08_MD_523_proj[1],x6G08_MD_523_proj[3],pch=20)
text(x6G08_MD_523_proj[1],x6G08_MD_523_proj[3],pos=4,label="523",col="black")

x6G08_MD_524_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_524.txt")
x6G08_MD_524_proj<-project.pca(x6G08_MD_524_raw,pcadata)
points(x6G08_MD_524_proj[1],x6G08_MD_524_proj[3],pch=20)
text(x6G08_MD_524_proj[1],x6G08_MD_524_proj[3],pos=4,label="524",col="black")

x6G08_MD_525_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_525.txt")
x6G08_MD_525_proj<-project.pca(x6G08_MD_525_raw,pcadata)
points(x6G08_MD_525_proj[1],x6G08_MD_525_proj[3],pch=20)
text(x6G08_MD_525_proj[1],x6G08_MD_525_proj[3],pos=4,label="525",col="blue")

x6G08_MD_526_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_526.txt")
x6G08_MD_526_proj<-project.pca(x6G08_MD_526_raw,pcadata)
points(x6G08_MD_526_proj[1],x6G08_MD_526_proj[3],pch=20)
text(x6G08_MD_526_proj[1],x6G08_MD_526_proj[3],pos=4,label="526",col="green")

x6G08_MD_527_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_527.txt")
x6G08_MD_527_proj<-project.pca(x6G08_MD_527_raw,pcadata)
points(x6G08_MD_527_proj[1],x6G08_MD_527_proj[3],pch=20)
text(x6G08_MD_527_proj[1],x6G08_MD_527_proj[3],pos=4,label="527",col="blue")

x6G08_MD_528_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_528.txt")
x6G08_MD_528_proj<-project.pca(x6G08_MD_528_raw,pcadata)
points(x6G08_MD_528_proj[1],x6G08_MD_528_proj[3],pch=20)
text(x6G08_MD_528_proj[1],x6G08_MD_528_proj[3],pos=4,label="528",col="green")

x6G08_MD_529_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_529.txt")
x6G08_MD_529_proj<-project.pca(x6G08_MD_529_raw,pcadata)
points(x6G08_MD_529_proj[1],x6G08_MD_529_proj[3],pch=20)
text(x6G08_MD_529_proj[1],x6G08_MD_529_proj[3],pos=4,label="529",col="blue")

x6G08_MD_530_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_530.txt")
x6G08_MD_530_proj<-project.pca(x6G08_MD_530_raw,pcadata)
points(x6G08_MD_530_proj[1],x6G08_MD_530_proj[3],pch=20)
text(x6G08_MD_530_proj[1],x6G08_MD_530_proj[3],pos=4,label="530",col="blue")

x6G08_MD_531_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_531.txt")
x6G08_MD_531_proj<-project.pca(x6G08_MD_531_raw,pcadata)
points(x6G08_MD_531_proj[1],x6G08_MD_531_proj[3],pch=20)
text(x6G08_MD_531_proj[1],x6G08_MD_531_proj[3],pos=4,label="531",col="blue")

x6G08_MD_532_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_532.txt")
x6G08_MD_532_proj<-project.pca(x6G08_MD_532_raw,pcadata)
points(x6G08_MD_532_proj[1],x6G08_MD_532_proj[3],pch=20)
text(x6G08_MD_532_proj[1],x6G08_MD_532_proj[3],pos=4,label="532",col="red")

x6G08_MD_533_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_533.txt")
x6G08_MD_533_proj<-project.pca(x6G08_MD_533_raw,pcadata)
points(x6G08_MD_533_proj[1],x6G08_MD_533_proj[3],pch=20)
text(x6G08_MD_533_proj[1],x6G08_MD_533_proj[3],pos=4,label="533",col="blue")

x6G08_MD_534_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_534.txt")
x6G08_MD_534_proj<-project.pca(x6G08_MD_534_raw,pcadata)
points(x6G08_MD_534_proj[1],x6G08_MD_534_proj[3],pch=20)
text(x6G08_MD_534_proj[1],x6G08_MD_534_proj[3],pos=4,label="534",col="blue")

x6G08_MD_535_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_535.txt")
x6G08_MD_535_proj<-project.pca(x6G08_MD_535_raw,pcadata)
points(x6G08_MD_535_proj[1],x6G08_MD_535_proj[3],pch=20)
text(x6G08_MD_535_proj[1],x6G08_MD_535_proj[3],pos=4,label="535",col="red")

x6G08_MD_536_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_536.txt")
x6G08_MD_536_proj<-project.pca(x6G08_MD_536_raw,pcadata)
points(x6G08_MD_536_proj[1],x6G08_MD_536_proj[3],pch=20)
text(x6G08_MD_536_proj[1],x6G08_MD_536_proj[3],pos=4,label="536",col="black")

x6G08_MD_537_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_537.txt")
x6G08_MD_537_proj<-project.pca(x6G08_MD_537_raw,pcadata)
points(x6G08_MD_537_proj[1],x6G08_MD_537_proj[3],pch=20)
text(x6G08_MD_537_proj[1],x6G08_MD_537_proj[3],pos=4,label="537",col="blue")

x6G08_MD_538_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_538.txt")
x6G08_MD_538_proj<-project.pca(x6G08_MD_538_raw,pcadata)
points(x6G08_MD_538_proj[1],x6G08_MD_538_proj[3],pch=20)
text(x6G08_MD_538_proj[1],x6G08_MD_538_proj[3],pos=4,label="538",col="red")

x6G08_MD_539_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_539.txt")
x6G08_MD_539_proj<-project.pca(x6G08_MD_539_raw,pcadata)
points(x6G08_MD_539_proj[1],x6G08_MD_539_proj[3],pch=20)
text(x6G08_MD_539_proj[1],x6G08_MD_539_proj[3],pos=4,label="539",col="green")

x6G08_MD_540_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_540.txt")
x6G08_MD_540_proj<-project.pca(x6G08_MD_540_raw,pcadata)
points(x6G08_MD_540_proj[1],x6G08_MD_540_proj[3],pch=20)
text(x6G08_MD_540_proj[1],x6G08_MD_540_proj[3],pos=4,label="540",col="red")

x6G08_MD_541_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_541.txt")
x6G08_MD_541_proj<-project.pca(x6G08_MD_541_raw,pcadata)
points(x6G08_MD_541_proj[1],x6G08_MD_541_proj[3],pch=20)
text(x6G08_MD_541_proj[1],x6G08_MD_541_proj[3],pos=4,label="541",col="red")

x6G08_MD_542_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_542.txt")
x6G08_MD_542_proj<-project.pca(x6G08_MD_542_raw,pcadata)
points(x6G08_MD_542_proj[1],x6G08_MD_542_proj[3],pch=20)
text(x6G08_MD_542_proj[1],x6G08_MD_542_proj[3],pos=4,label="542",col="green")

x6G08_MD_543_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_543.txt")
x6G08_MD_543_proj<-project.pca(x6G08_MD_543_raw,pcadata)
points(x6G08_MD_543_proj[1],x6G08_MD_543_proj[3],pch=20)
text(x6G08_MD_543_proj[1],x6G08_MD_543_proj[3],pos=4,label="543",col="green")

x6G08_MD_544_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_544.txt")
x6G08_MD_544_proj<-project.pca(x6G08_MD_544_raw,pcadata)
points(x6G08_MD_544_proj[1],x6G08_MD_544_proj[3],pch=20)
text(x6G08_MD_544_proj[1],x6G08_MD_544_proj[3],pos=4,label="544",col="green")

x6G08_MD_545_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_545.txt")
x6G08_MD_545_proj<-project.pca(x6G08_MD_545_raw,pcadata)
points(x6G08_MD_545_proj[1],x6G08_MD_545_proj[3],pch=20)
text(x6G08_MD_545_proj[1],x6G08_MD_545_proj[3],pos=4,label="545",col="blue")

x6G08_MD_546_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_546.txt")
x6G08_MD_546_proj<-project.pca(x6G08_MD_546_raw,pcadata)
points(x6G08_MD_546_proj[1],x6G08_MD_546_proj[3],pch=20)
text(x6G08_MD_546_proj[1],x6G08_MD_546_proj[3],pos=4,label="546",col="red")

x6G08_MD_547_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_547.txt")
x6G08_MD_547_proj<-project.pca(x6G08_MD_547_raw,pcadata)
points(x6G08_MD_547_proj[1],x6G08_MD_547_proj[3],pch=20)
text(x6G08_MD_547_proj[1],x6G08_MD_547_proj[3],pos=4,label="547",col="green")

x6G08_MD_548_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_548.txt")
x6G08_MD_548_proj<-project.pca(x6G08_MD_548_raw,pcadata)
points(x6G08_MD_548_proj[1],x6G08_MD_548_proj[3],pch=20)
text(x6G08_MD_548_proj[1],x6G08_MD_548_proj[3],pos=4,label="548",col="blue")

x6G08_MD_549_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_549.txt")
x6G08_MD_549_proj<-project.pca(x6G08_MD_549_raw,pcadata)
points(x6G08_MD_549_proj[1],x6G08_MD_549_proj[3],pch=20)
text(x6G08_MD_549_proj[1],x6G08_MD_549_proj[3],pos=4,label="549",col="red")

x6G08_MD_550_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_550.txt")
x6G08_MD_550_proj<-project.pca(x6G08_MD_550_raw,pcadata)
points(x6G08_MD_550_proj[1],x6G08_MD_550_proj[3],pch=20)
text(x6G08_MD_550_proj[1],x6G08_MD_550_proj[3],pos=4,label="550",col="green")

x6G08_MD_551_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_551.txt")
x6G08_MD_551_proj<-project.pca(x6G08_MD_551_raw,pcadata)
points(x6G08_MD_551_proj[1],x6G08_MD_551_proj[3],pch=20)
text(x6G08_MD_551_proj[1],x6G08_MD_551_proj[3],pos=4,label="551",col="red")

x6G08_MD_552_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_552.txt")
x6G08_MD_552_proj<-project.pca(x6G08_MD_552_raw,pcadata)
points(x6G08_MD_552_proj[1],x6G08_MD_552_proj[3],pch=20)
text(x6G08_MD_552_proj[1],x6G08_MD_552_proj[3],pos=4,label="552",col="red")

x6G08_MD_553_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_553.txt")
x6G08_MD_553_proj<-project.pca(x6G08_MD_553_raw,pcadata)
points(x6G08_MD_553_proj[1],x6G08_MD_553_proj[3],pch=20)
text(x6G08_MD_553_proj[1],x6G08_MD_553_proj[3],pos=4,label="553",col="green")

x6G08_MD_554_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_554.txt")
x6G08_MD_554_proj<-project.pca(x6G08_MD_554_raw,pcadata)
points(x6G08_MD_554_proj[1],x6G08_MD_554_proj[3],pch=20)
text(x6G08_MD_554_proj[1],x6G08_MD_554_proj[3],pos=4,label="554",col="green")

x6G08_MD_555_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_555.txt")
x6G08_MD_555_proj<-project.pca(x6G08_MD_555_raw,pcadata)
points(x6G08_MD_555_proj[1],x6G08_MD_555_proj[3],pch=20)
text(x6G08_MD_555_proj[1],x6G08_MD_555_proj[3],pos=4,label="555",col="blue")

x6G08_MD_556_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_556.txt")
x6G08_MD_556_proj<-project.pca(x6G08_MD_556_raw,pcadata)
points(x6G08_MD_556_proj[1],x6G08_MD_556_proj[3],pch=20)
text(x6G08_MD_556_proj[1],x6G08_MD_556_proj[3],pos=4,label="556",col="green")

x6G08_MD_557_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_557.txt")
x6G08_MD_557_proj<-project.pca(x6G08_MD_557_raw,pcadata)
points(x6G08_MD_557_proj[1],x6G08_MD_557_proj[3],pch=20)
text(x6G08_MD_557_proj[1],x6G08_MD_557_proj[3],pos=4,label="557",col="blue")

x6G08_MD_558_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_558.txt")
x6G08_MD_558_proj<-project.pca(x6G08_MD_558_raw,pcadata)
points(x6G08_MD_558_proj[1],x6G08_MD_558_proj[3],pch=20)
text(x6G08_MD_558_proj[1],x6G08_MD_558_proj[3],pos=4,label="558",col="blue")

x6G08_MD_559_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_559.txt")
x6G08_MD_559_proj<-project.pca(x6G08_MD_559_raw,pcadata)
points(x6G08_MD_559_proj[1],x6G08_MD_559_proj[3],pch=20)
text(x6G08_MD_559_proj[1],x6G08_MD_559_proj[3],pos=4,label="559",col="blue")

x6G08_MD_560_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_560.txt")
x6G08_MD_560_proj<-project.pca(x6G08_MD_560_raw,pcadata)
points(x6G08_MD_560_proj[1],x6G08_MD_560_proj[3],pch=20)
text(x6G08_MD_560_proj[1],x6G08_MD_560_proj[3],pos=4,label="560",col="blue")

x6G08_MD_561_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_561.txt")
x6G08_MD_561_proj<-project.pca(x6G08_MD_561_raw,pcadata)
points(x6G08_MD_561_proj[1],x6G08_MD_561_proj[3],pch=20)
text(x6G08_MD_561_proj[1],x6G08_MD_561_proj[3],pos=4,label="561",col="red")

x6G08_MD_562_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_562.txt")
x6G08_MD_562_proj<-project.pca(x6G08_MD_562_raw,pcadata)
points(x6G08_MD_562_proj[1],x6G08_MD_562_proj[3],pch=20)
text(x6G08_MD_562_proj[1],x6G08_MD_562_proj[3],pos=4,label="562",col="green")

x6G08_MD_563_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_563.txt")
x6G08_MD_563_proj<-project.pca(x6G08_MD_563_raw,pcadata)
points(x6G08_MD_563_proj[1],x6G08_MD_563_proj[3],pch=20)
text(x6G08_MD_563_proj[1],x6G08_MD_563_proj[3],pos=4,label="563",col="red")

x6G08_MD_564_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_564.txt")
x6G08_MD_564_proj<-project.pca(x6G08_MD_564_raw,pcadata)
points(x6G08_MD_564_proj[1],x6G08_MD_564_proj[3],pch=20)
text(x6G08_MD_564_proj[1],x6G08_MD_564_proj[3],pos=4,label="564",col="blue")

x6G08_MD_565_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_565.txt")
x6G08_MD_565_proj<-project.pca(x6G08_MD_565_raw,pcadata)
points(x6G08_MD_565_proj[1],x6G08_MD_565_proj[3],pch=20)
text(x6G08_MD_565_proj[1],x6G08_MD_565_proj[3],pos=4,label="565",col="blue")

x6G08_MD_566_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_566.txt")
x6G08_MD_566_proj<-project.pca(x6G08_MD_566_raw,pcadata)
points(x6G08_MD_566_proj[1],x6G08_MD_566_proj[3],pch=20)
text(x6G08_MD_566_proj[1],x6G08_MD_566_proj[3],pos=4,label="566",col="red")

x6G08_MD_567_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_567.txt")
x6G08_MD_567_proj<-project.pca(x6G08_MD_567_raw,pcadata)
points(x6G08_MD_567_proj[1],x6G08_MD_567_proj[3],pch=20)
text(x6G08_MD_567_proj[1],x6G08_MD_567_proj[3],pos=4,label="567",col="red")

x6G08_MD_568_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_568.txt")
x6G08_MD_568_proj<-project.pca(x6G08_MD_568_raw,pcadata)
points(x6G08_MD_568_proj[1],x6G08_MD_568_proj[3],pch=20)
text(x6G08_MD_568_proj[1],x6G08_MD_568_proj[3],pos=4,label="568",col="blue")

x6G08_MD_569_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_569.txt")
x6G08_MD_569_proj<-project.pca(x6G08_MD_569_raw,pcadata)
points(x6G08_MD_569_proj[1],x6G08_MD_569_proj[3],pch=20)
text(x6G08_MD_569_proj[1],x6G08_MD_569_proj[3],pos=4,label="569",col="green")

x6G08_MD_570_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_570.txt")
x6G08_MD_570_proj<-project.pca(x6G08_MD_570_raw,pcadata)
points(x6G08_MD_570_proj[1],x6G08_MD_570_proj[3],pch=20)
text(x6G08_MD_570_proj[1],x6G08_MD_570_proj[3],pos=4,label="570",col="green")

x6G08_MD_571_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_571.txt")
x6G08_MD_571_proj<-project.pca(x6G08_MD_571_raw,pcadata)
points(x6G08_MD_571_proj[1],x6G08_MD_571_proj[3],pch=20)
text(x6G08_MD_571_proj[1],x6G08_MD_571_proj[3],pos=4,label="571",col="red")

x6G08_MD_572_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_572.txt")
x6G08_MD_572_proj<-project.pca(x6G08_MD_572_raw,pcadata)
points(x6G08_MD_572_proj[1],x6G08_MD_572_proj[3],pch=20)
text(x6G08_MD_572_proj[1],x6G08_MD_572_proj[3],pos=4,label="572",col="blue")

x6G08_MD_573_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_573.txt")
x6G08_MD_573_proj<-project.pca(x6G08_MD_573_raw,pcadata)
points(x6G08_MD_573_proj[1],x6G08_MD_573_proj[3],pch=20)
text(x6G08_MD_573_proj[1],x6G08_MD_573_proj[3],pos=4,label="573",col="black")

x6G08_MD_574_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_574.txt")
x6G08_MD_574_proj<-project.pca(x6G08_MD_574_raw,pcadata)
points(x6G08_MD_574_proj[1],x6G08_MD_574_proj[3],pch=20)
text(x6G08_MD_574_proj[1],x6G08_MD_574_proj[3],pos=4,label="574",col="black")

x6G08_MD_575_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_575.txt")
x6G08_MD_575_proj<-project.pca(x6G08_MD_575_raw,pcadata)
points(x6G08_MD_575_proj[1],x6G08_MD_575_proj[3],pch=20)
text(x6G08_MD_575_proj[1],x6G08_MD_575_proj[3],pos=4,label="575",col="red")

x6G08_MD_576_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_576.txt")
x6G08_MD_576_proj<-project.pca(x6G08_MD_576_raw,pcadata)
points(x6G08_MD_576_proj[1],x6G08_MD_576_proj[3],pch=20)
text(x6G08_MD_576_proj[1],x6G08_MD_576_proj[3],pos=4,label="576",col="blue")

x6G08_MD_577_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_577.txt")
x6G08_MD_577_proj<-project.pca(x6G08_MD_577_raw,pcadata)
points(x6G08_MD_577_proj[1],x6G08_MD_577_proj[3],pch=20)
text(x6G08_MD_577_proj[1],x6G08_MD_577_proj[3],pos=4,label="577",col="black")

x6G08_MD_578_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_578.txt")
x6G08_MD_578_proj<-project.pca(x6G08_MD_578_raw,pcadata)
points(x6G08_MD_578_proj[1],x6G08_MD_578_proj[3],pch=20)
text(x6G08_MD_578_proj[1],x6G08_MD_578_proj[3],pos=4,label="578",col="blue")

x6G08_MD_579_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_579.txt")
x6G08_MD_579_proj<-project.pca(x6G08_MD_579_raw,pcadata)
points(x6G08_MD_579_proj[1],x6G08_MD_579_proj[3],pch=20)
text(x6G08_MD_579_proj[1],x6G08_MD_579_proj[3],pos=4,label="579",col="blue")

x6G08_MD_580_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_580.txt")
x6G08_MD_580_proj<-project.pca(x6G08_MD_580_raw,pcadata)
points(x6G08_MD_580_proj[1],x6G08_MD_580_proj[3],pch=20)
text(x6G08_MD_580_proj[1],x6G08_MD_580_proj[3],pos=4,label="580",col="green")

x6G08_MD_581_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_581.txt")
x6G08_MD_581_proj<-project.pca(x6G08_MD_581_raw,pcadata)
points(x6G08_MD_581_proj[1],x6G08_MD_581_proj[3],pch=20)
text(x6G08_MD_581_proj[1],x6G08_MD_581_proj[3],pos=4,label="581",col="green")

x6G08_MD_582_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_582.txt")
x6G08_MD_582_proj<-project.pca(x6G08_MD_582_raw,pcadata)
points(x6G08_MD_582_proj[1],x6G08_MD_582_proj[3],pch=20)
text(x6G08_MD_582_proj[1],x6G08_MD_582_proj[3],pos=4,label="582",col="green")

x6G08_MD_583_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_583.txt")
x6G08_MD_583_proj<-project.pca(x6G08_MD_583_raw,pcadata)
points(x6G08_MD_583_proj[1],x6G08_MD_583_proj[3],pch=20)
text(x6G08_MD_583_proj[1],x6G08_MD_583_proj[3],pos=4,label="583",col="blue")

x6G08_MD_584_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_584.txt")
x6G08_MD_584_proj<-project.pca(x6G08_MD_584_raw,pcadata)
points(x6G08_MD_584_proj[1],x6G08_MD_584_proj[3],pch=20)
text(x6G08_MD_584_proj[1],x6G08_MD_584_proj[3],pos=4,label="584",col="red")

x6G08_MD_585_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_585.txt")
x6G08_MD_585_proj<-project.pca(x6G08_MD_585_raw,pcadata)
points(x6G08_MD_585_proj[1],x6G08_MD_585_proj[3],pch=20)
text(x6G08_MD_585_proj[1],x6G08_MD_585_proj[3],pos=4,label="585",col="blue")

x6G08_MD_586_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_586.txt")
x6G08_MD_586_proj<-project.pca(x6G08_MD_586_raw,pcadata)
points(x6G08_MD_586_proj[1],x6G08_MD_586_proj[3],pch=20)
text(x6G08_MD_586_proj[1],x6G08_MD_586_proj[3],pos=4,label="586",col="blue")

x6G08_MD_587_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_587.txt")
x6G08_MD_587_proj<-project.pca(x6G08_MD_587_raw,pcadata)
points(x6G08_MD_587_proj[1],x6G08_MD_587_proj[3],pch=20)
text(x6G08_MD_587_proj[1],x6G08_MD_587_proj[3],pos=4,label="587",col="black")

x6G08_MD_588_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_588.txt")
x6G08_MD_588_proj<-project.pca(x6G08_MD_588_raw,pcadata)
points(x6G08_MD_588_proj[1],x6G08_MD_588_proj[3],pch=20)
text(x6G08_MD_588_proj[1],x6G08_MD_588_proj[3],pos=4,label="588",col="green")

x6G08_MD_589_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_589.txt")
x6G08_MD_589_proj<-project.pca(x6G08_MD_589_raw,pcadata)
points(x6G08_MD_589_proj[1],x6G08_MD_589_proj[3],pch=20)
text(x6G08_MD_589_proj[1],x6G08_MD_589_proj[3],pos=4,label="589",col="blue")

x6G08_MD_590_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_590.txt")
x6G08_MD_590_proj<-project.pca(x6G08_MD_590_raw,pcadata)
points(x6G08_MD_590_proj[1],x6G08_MD_590_proj[3],pch=20)
text(x6G08_MD_590_proj[1],x6G08_MD_590_proj[3],pos=4,label="590",col="black")

x6G08_MD_591_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_591.txt")
x6G08_MD_591_proj<-project.pca(x6G08_MD_591_raw,pcadata)
points(x6G08_MD_591_proj[1],x6G08_MD_591_proj[3],pch=20)
text(x6G08_MD_591_proj[1],x6G08_MD_591_proj[3],pos=4,label="591",col="blue")

x6G08_MD_592_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_592.txt")
x6G08_MD_592_proj<-project.pca(x6G08_MD_592_raw,pcadata)
points(x6G08_MD_592_proj[1],x6G08_MD_592_proj[3],pch=20)
text(x6G08_MD_592_proj[1],x6G08_MD_592_proj[3],pos=4,label="592",col="green")

x6G08_MD_593_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_593.txt")
x6G08_MD_593_proj<-project.pca(x6G08_MD_593_raw,pcadata)
points(x6G08_MD_593_proj[1],x6G08_MD_593_proj[3],pch=20)
text(x6G08_MD_593_proj[1],x6G08_MD_593_proj[3],pos=4,label="593",col="blue")

x6G08_MD_594_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_594.txt")
x6G08_MD_594_proj<-project.pca(x6G08_MD_594_raw,pcadata)
points(x6G08_MD_594_proj[1],x6G08_MD_594_proj[3],pch=20)
text(x6G08_MD_594_proj[1],x6G08_MD_594_proj[3],pos=4,label="594",col="blue")

x6G08_MD_595_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_595.txt")
x6G08_MD_595_proj<-project.pca(x6G08_MD_595_raw,pcadata)
points(x6G08_MD_595_proj[1],x6G08_MD_595_proj[3],pch=20)
text(x6G08_MD_595_proj[1],x6G08_MD_595_proj[3],pos=4,label="595",col="black")

x6G08_MD_596_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_596.txt")
x6G08_MD_596_proj<-project.pca(x6G08_MD_596_raw,pcadata)
points(x6G08_MD_596_proj[1],x6G08_MD_596_proj[3],pch=20)
text(x6G08_MD_596_proj[1],x6G08_MD_596_proj[3],pos=4,label="596",col="black")

x6G08_MD_597_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_597.txt")
x6G08_MD_597_proj<-project.pca(x6G08_MD_597_raw,pcadata)
points(x6G08_MD_597_proj[1],x6G08_MD_597_proj[3],pch=20)
text(x6G08_MD_597_proj[1],x6G08_MD_597_proj[3],pos=4,label="597",col="blue")

x6G08_MD_598_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_598.txt")
x6G08_MD_598_proj<-project.pca(x6G08_MD_598_raw,pcadata)
points(x6G08_MD_598_proj[1],x6G08_MD_598_proj[3],pch=20)
text(x6G08_MD_598_proj[1],x6G08_MD_598_proj[3],pos=4,label="598",col="blue")

x6G08_MD_599_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_599.txt")
x6G08_MD_599_proj<-project.pca(x6G08_MD_599_raw,pcadata)
points(x6G08_MD_599_proj[1],x6G08_MD_599_proj[3],pch=20)
text(x6G08_MD_599_proj[1],x6G08_MD_599_proj[3],pos=4,label="599",col="green")

x6G08_MD_600_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_600.txt")
x6G08_MD_600_proj<-project.pca(x6G08_MD_600_raw,pcadata)
points(x6G08_MD_600_proj[1],x6G08_MD_600_proj[3],pch=20)
text(x6G08_MD_600_proj[1],x6G08_MD_600_proj[3],pos=4,label="600",col="green")

x6G08_MD_601_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_601.txt")
x6G08_MD_601_proj<-project.pca(x6G08_MD_601_raw,pcadata)
points(x6G08_MD_601_proj[1],x6G08_MD_601_proj[3],pch=20)
text(x6G08_MD_601_proj[1],x6G08_MD_601_proj[3],pos=4,label="601",col="green")

x6G08_MD_602_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_602.txt")
x6G08_MD_602_proj<-project.pca(x6G08_MD_602_raw,pcadata)
points(x6G08_MD_602_proj[1],x6G08_MD_602_proj[3],pch=20)
text(x6G08_MD_602_proj[1],x6G08_MD_602_proj[3],pos=4,label="602",col="green")

x6G08_MD_603_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_603.txt")
x6G08_MD_603_proj<-project.pca(x6G08_MD_603_raw,pcadata)
points(x6G08_MD_603_proj[1],x6G08_MD_603_proj[3],pch=20)
text(x6G08_MD_603_proj[1],x6G08_MD_603_proj[3],pos=4,label="603",col="blue")

x6G08_MD_604_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_604.txt")
x6G08_MD_604_proj<-project.pca(x6G08_MD_604_raw,pcadata)
points(x6G08_MD_604_proj[1],x6G08_MD_604_proj[3],pch=20)
text(x6G08_MD_604_proj[1],x6G08_MD_604_proj[3],pos=4,label="604",col="black")

x6G08_MD_605_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_605.txt")
x6G08_MD_605_proj<-project.pca(x6G08_MD_605_raw,pcadata)
points(x6G08_MD_605_proj[1],x6G08_MD_605_proj[3],pch=20)
text(x6G08_MD_605_proj[1],x6G08_MD_605_proj[3],pos=4,label="605",col="blue")

x6G08_MD_606_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_606.txt")
x6G08_MD_606_proj<-project.pca(x6G08_MD_606_raw,pcadata)
points(x6G08_MD_606_proj[1],x6G08_MD_606_proj[3],pch=20)
text(x6G08_MD_606_proj[1],x6G08_MD_606_proj[3],pos=4,label="606",col="red")

x6G08_MD_607_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_607.txt")
x6G08_MD_607_proj<-project.pca(x6G08_MD_607_raw,pcadata)
points(x6G08_MD_607_proj[1],x6G08_MD_607_proj[3],pch=20)
text(x6G08_MD_607_proj[1],x6G08_MD_607_proj[3],pos=4,label="607",col="red")

x6G08_MD_608_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_608.txt")
x6G08_MD_608_proj<-project.pca(x6G08_MD_608_raw,pcadata)
points(x6G08_MD_608_proj[1],x6G08_MD_608_proj[3],pch=20)
text(x6G08_MD_608_proj[1],x6G08_MD_608_proj[3],pos=4,label="608",col="red")

x6G08_MD_609_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_609.txt")
x6G08_MD_609_proj<-project.pca(x6G08_MD_609_raw,pcadata)
points(x6G08_MD_609_proj[1],x6G08_MD_609_proj[3],pch=20)
text(x6G08_MD_609_proj[1],x6G08_MD_609_proj[3],pos=4,label="609",col="blue")

x6G08_MD_610_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_610.txt")
x6G08_MD_610_proj<-project.pca(x6G08_MD_610_raw,pcadata)
points(x6G08_MD_610_proj[1],x6G08_MD_610_proj[3],pch=20)
text(x6G08_MD_610_proj[1],x6G08_MD_610_proj[3],pos=4,label="610",col="black")

x6G08_MD_611_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_611.txt")
x6G08_MD_611_proj<-project.pca(x6G08_MD_611_raw,pcadata)
points(x6G08_MD_611_proj[1],x6G08_MD_611_proj[3],pch=20)
text(x6G08_MD_611_proj[1],x6G08_MD_611_proj[3],pos=4,label="611",col="blue")

x6G08_MD_612_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_612.txt")
x6G08_MD_612_proj<-project.pca(x6G08_MD_612_raw,pcadata)
points(x6G08_MD_612_proj[1],x6G08_MD_612_proj[3],pch=20)
text(x6G08_MD_612_proj[1],x6G08_MD_612_proj[3],pos=4,label="612",col="green")

x6G08_MD_613_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_613.txt")
x6G08_MD_613_proj<-project.pca(x6G08_MD_613_raw,pcadata)
points(x6G08_MD_613_proj[1],x6G08_MD_613_proj[3],pch=20)
text(x6G08_MD_613_proj[1],x6G08_MD_613_proj[3],pos=4,label="613",col="red")

x6G08_MD_614_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_614.txt")
x6G08_MD_614_proj<-project.pca(x6G08_MD_614_raw,pcadata)
points(x6G08_MD_614_proj[1],x6G08_MD_614_proj[3],pch=20)
text(x6G08_MD_614_proj[1],x6G08_MD_614_proj[3],pos=4,label="614",col="blue")

x6G08_MD_615_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_615.txt")
x6G08_MD_615_proj<-project.pca(x6G08_MD_615_raw,pcadata)
points(x6G08_MD_615_proj[1],x6G08_MD_615_proj[3],pch=20)
text(x6G08_MD_615_proj[1],x6G08_MD_615_proj[3],pos=4,label="615",col="blue")

x6G08_MD_616_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_616.txt")
x6G08_MD_616_proj<-project.pca(x6G08_MD_616_raw,pcadata)
points(x6G08_MD_616_proj[1],x6G08_MD_616_proj[3],pch=20)
text(x6G08_MD_616_proj[1],x6G08_MD_616_proj[3],pos=4,label="616",col="blue")

x6G08_MD_617_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_617.txt")
x6G08_MD_617_proj<-project.pca(x6G08_MD_617_raw,pcadata)
points(x6G08_MD_617_proj[1],x6G08_MD_617_proj[3],pch=20)
text(x6G08_MD_617_proj[1],x6G08_MD_617_proj[3],pos=4,label="617",col="green")

x6G08_MD_618_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_618.txt")
x6G08_MD_618_proj<-project.pca(x6G08_MD_618_raw,pcadata)
points(x6G08_MD_618_proj[1],x6G08_MD_618_proj[3],pch=20)
text(x6G08_MD_618_proj[1],x6G08_MD_618_proj[3],pos=4,label="618",col="green")

x6G08_MD_619_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_619.txt")
x6G08_MD_619_proj<-project.pca(x6G08_MD_619_raw,pcadata)
points(x6G08_MD_619_proj[1],x6G08_MD_619_proj[3],pch=20)
text(x6G08_MD_619_proj[1],x6G08_MD_619_proj[3],pos=4,label="619",col="blue")

x6G08_MD_620_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_620.txt")
x6G08_MD_620_proj<-project.pca(x6G08_MD_620_raw,pcadata)
points(x6G08_MD_620_proj[1],x6G08_MD_620_proj[3],pch=20)
text(x6G08_MD_620_proj[1],x6G08_MD_620_proj[3],pos=4,label="620",col="green")

x6G08_MD_621_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_621.txt")
x6G08_MD_621_proj<-project.pca(x6G08_MD_621_raw,pcadata)
points(x6G08_MD_621_proj[1],x6G08_MD_621_proj[3],pch=20)
text(x6G08_MD_621_proj[1],x6G08_MD_621_proj[3],pos=4,label="621",col="blue")

x6G08_MD_622_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_622.txt")
x6G08_MD_622_proj<-project.pca(x6G08_MD_622_raw,pcadata)
points(x6G08_MD_622_proj[1],x6G08_MD_622_proj[3],pch=20)
text(x6G08_MD_622_proj[1],x6G08_MD_622_proj[3],pos=4,label="622",col="blue")

x6G08_MD_623_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_623.txt")
x6G08_MD_623_proj<-project.pca(x6G08_MD_623_raw,pcadata)
points(x6G08_MD_623_proj[1],x6G08_MD_623_proj[3],pch=20)
text(x6G08_MD_623_proj[1],x6G08_MD_623_proj[3],pos=4,label="623",col="blue")

x6G08_MD_624_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_624.txt")
x6G08_MD_624_proj<-project.pca(x6G08_MD_624_raw,pcadata)
points(x6G08_MD_624_proj[1],x6G08_MD_624_proj[3],pch=20)
text(x6G08_MD_624_proj[1],x6G08_MD_624_proj[3],pos=4,label="624",col="blue")

x6G08_MD_625_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_625.txt")
x6G08_MD_625_proj<-project.pca(x6G08_MD_625_raw,pcadata)
points(x6G08_MD_625_proj[1],x6G08_MD_625_proj[3],pch=20)
text(x6G08_MD_625_proj[1],x6G08_MD_625_proj[3],pos=4,label="625",col="blue")

x6G08_MD_626_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_626.txt")
x6G08_MD_626_proj<-project.pca(x6G08_MD_626_raw,pcadata)
points(x6G08_MD_626_proj[1],x6G08_MD_626_proj[3],pch=20)
text(x6G08_MD_626_proj[1],x6G08_MD_626_proj[3],pos=4,label="626",col="blue")

x6G08_MD_627_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_627.txt")
x6G08_MD_627_proj<-project.pca(x6G08_MD_627_raw,pcadata)
points(x6G08_MD_627_proj[1],x6G08_MD_627_proj[3],pch=20)
text(x6G08_MD_627_proj[1],x6G08_MD_627_proj[3],pos=4,label="627",col="blue")

x6G08_MD_628_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_628.txt")
x6G08_MD_628_proj<-project.pca(x6G08_MD_628_raw,pcadata)
points(x6G08_MD_628_proj[1],x6G08_MD_628_proj[3],pch=20)
text(x6G08_MD_628_proj[1],x6G08_MD_628_proj[3],pos=4,label="628",col="red")

x6G08_MD_629_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_629.txt")
x6G08_MD_629_proj<-project.pca(x6G08_MD_629_raw,pcadata)
points(x6G08_MD_629_proj[1],x6G08_MD_629_proj[3],pch=20)
text(x6G08_MD_629_proj[1],x6G08_MD_629_proj[3],pos=4,label="629",col="red")

x6G08_MD_630_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_630.txt")
x6G08_MD_630_proj<-project.pca(x6G08_MD_630_raw,pcadata)
points(x6G08_MD_630_proj[1],x6G08_MD_630_proj[3],pch=20)
text(x6G08_MD_630_proj[1],x6G08_MD_630_proj[3],pos=4,label="630",col="blue")

x6G08_MD_631_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_631.txt")
x6G08_MD_631_proj<-project.pca(x6G08_MD_631_raw,pcadata)
points(x6G08_MD_631_proj[1],x6G08_MD_631_proj[3],pch=20)
text(x6G08_MD_631_proj[1],x6G08_MD_631_proj[3],pos=4,label="631",col="green")

x6G08_MD_632_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_632.txt")
x6G08_MD_632_proj<-project.pca(x6G08_MD_632_raw,pcadata)
points(x6G08_MD_632_proj[1],x6G08_MD_632_proj[3],pch=20)
text(x6G08_MD_632_proj[1],x6G08_MD_632_proj[3],pos=4,label="632",col="green")

x6G08_MD_633_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_633.txt")
x6G08_MD_633_proj<-project.pca(x6G08_MD_633_raw,pcadata)
points(x6G08_MD_633_proj[1],x6G08_MD_633_proj[3],pch=20)
text(x6G08_MD_633_proj[1],x6G08_MD_633_proj[3],pos=4,label="633",col="red")

x6G08_MD_634_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_634.txt")
x6G08_MD_634_proj<-project.pca(x6G08_MD_634_raw,pcadata)
points(x6G08_MD_634_proj[1],x6G08_MD_634_proj[3],pch=20)
text(x6G08_MD_634_proj[1],x6G08_MD_634_proj[3],pos=4,label="634",col="green")

x6G08_MD_635_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_635.txt")
x6G08_MD_635_proj<-project.pca(x6G08_MD_635_raw,pcadata)
points(x6G08_MD_635_proj[1],x6G08_MD_635_proj[3],pch=20)
text(x6G08_MD_635_proj[1],x6G08_MD_635_proj[3],pos=4,label="635",col="blue")

x6G08_MD_636_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_636.txt")
x6G08_MD_636_proj<-project.pca(x6G08_MD_636_raw,pcadata)
points(x6G08_MD_636_proj[1],x6G08_MD_636_proj[3],pch=20)
text(x6G08_MD_636_proj[1],x6G08_MD_636_proj[3],pos=4,label="636",col="black")

x6G08_MD_637_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_637.txt")
x6G08_MD_637_proj<-project.pca(x6G08_MD_637_raw,pcadata)
points(x6G08_MD_637_proj[1],x6G08_MD_637_proj[3],pch=20)
text(x6G08_MD_637_proj[1],x6G08_MD_637_proj[3],pos=4,label="637",col="blue")

x6G08_MD_638_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_638.txt")
x6G08_MD_638_proj<-project.pca(x6G08_MD_638_raw,pcadata)
points(x6G08_MD_638_proj[1],x6G08_MD_638_proj[3],pch=20)
text(x6G08_MD_638_proj[1],x6G08_MD_638_proj[3],pos=4,label="638",col="blue")

x6G08_MD_639_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_639.txt")
x6G08_MD_639_proj<-project.pca(x6G08_MD_639_raw,pcadata)
points(x6G08_MD_639_proj[1],x6G08_MD_639_proj[3],pch=20)
text(x6G08_MD_639_proj[1],x6G08_MD_639_proj[3],pos=4,label="639",col="blue")

x6G08_MD_640_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_640.txt")
x6G08_MD_640_proj<-project.pca(x6G08_MD_640_raw,pcadata)
points(x6G08_MD_640_proj[1],x6G08_MD_640_proj[3],pch=20)
text(x6G08_MD_640_proj[1],x6G08_MD_640_proj[3],pos=4,label="640",col="green")

x6G08_MD_641_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_641.txt")
x6G08_MD_641_proj<-project.pca(x6G08_MD_641_raw,pcadata)
points(x6G08_MD_641_proj[1],x6G08_MD_641_proj[3],pch=20)
text(x6G08_MD_641_proj[1],x6G08_MD_641_proj[3],pos=4,label="641",col="green")

x6G08_MD_642_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_642.txt")
x6G08_MD_642_proj<-project.pca(x6G08_MD_642_raw,pcadata)
points(x6G08_MD_642_proj[1],x6G08_MD_642_proj[3],pch=20)
text(x6G08_MD_642_proj[1],x6G08_MD_642_proj[3],pos=4,label="642",col="blue")

x6G08_MD_643_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_643.txt")
x6G08_MD_643_proj<-project.pca(x6G08_MD_643_raw,pcadata)
points(x6G08_MD_643_proj[1],x6G08_MD_643_proj[3],pch=20)
text(x6G08_MD_643_proj[1],x6G08_MD_643_proj[3],pos=4,label="643",col="blue")

x6G08_MD_644_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_644.txt")
x6G08_MD_644_proj<-project.pca(x6G08_MD_644_raw,pcadata)
points(x6G08_MD_644_proj[1],x6G08_MD_644_proj[3],pch=20)
text(x6G08_MD_644_proj[1],x6G08_MD_644_proj[3],pos=4,label="644",col="green")

x6G08_MD_645_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_645.txt")
x6G08_MD_645_proj<-project.pca(x6G08_MD_645_raw,pcadata)
points(x6G08_MD_645_proj[1],x6G08_MD_645_proj[3],pch=20)
text(x6G08_MD_645_proj[1],x6G08_MD_645_proj[3],pos=4,label="645",col="blue")

x6G08_MD_646_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_646.txt")
x6G08_MD_646_proj<-project.pca(x6G08_MD_646_raw,pcadata)
points(x6G08_MD_646_proj[1],x6G08_MD_646_proj[3],pch=20)
text(x6G08_MD_646_proj[1],x6G08_MD_646_proj[3],pos=4,label="646",col="blue")

x6G08_MD_647_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_647.txt")
x6G08_MD_647_proj<-project.pca(x6G08_MD_647_raw,pcadata)
points(x6G08_MD_647_proj[1],x6G08_MD_647_proj[3],pch=20)
text(x6G08_MD_647_proj[1],x6G08_MD_647_proj[3],pos=4,label="647",col="red")

x6G08_MD_648_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_648.txt")
x6G08_MD_648_proj<-project.pca(x6G08_MD_648_raw,pcadata)
points(x6G08_MD_648_proj[1],x6G08_MD_648_proj[3],pch=20)
text(x6G08_MD_648_proj[1],x6G08_MD_648_proj[3],pos=4,label="648",col="red")

x6G08_MD_649_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_649.txt")
x6G08_MD_649_proj<-project.pca(x6G08_MD_649_raw,pcadata)
points(x6G08_MD_649_proj[1],x6G08_MD_649_proj[3],pch=20)
text(x6G08_MD_649_proj[1],x6G08_MD_649_proj[3],pos=4,label="649",col="red")

x6G08_MD_650_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_650.txt")
x6G08_MD_650_proj<-project.pca(x6G08_MD_650_raw,pcadata)
points(x6G08_MD_650_proj[1],x6G08_MD_650_proj[3],pch=20)
text(x6G08_MD_650_proj[1],x6G08_MD_650_proj[3],pos=4,label="650",col="red")

x6G08_MD_651_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_651.txt")
x6G08_MD_651_proj<-project.pca(x6G08_MD_651_raw,pcadata)
points(x6G08_MD_651_proj[1],x6G08_MD_651_proj[3],pch=20)
text(x6G08_MD_651_proj[1],x6G08_MD_651_proj[3],pos=4,label="651",col="red")

x6G08_MD_652_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_652.txt")
x6G08_MD_652_proj<-project.pca(x6G08_MD_652_raw,pcadata)
points(x6G08_MD_652_proj[1],x6G08_MD_652_proj[3],pch=20)
text(x6G08_MD_652_proj[1],x6G08_MD_652_proj[3],pos=4,label="652",col="red")

x6G08_MD_653_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_653.txt")
x6G08_MD_653_proj<-project.pca(x6G08_MD_653_raw,pcadata)
points(x6G08_MD_653_proj[1],x6G08_MD_653_proj[3],pch=20)
text(x6G08_MD_653_proj[1],x6G08_MD_653_proj[3],pos=4,label="653",col="red")

x6G08_MD_654_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_654.txt")
x6G08_MD_654_proj<-project.pca(x6G08_MD_654_raw,pcadata)
points(x6G08_MD_654_proj[1],x6G08_MD_654_proj[3],pch=20)
text(x6G08_MD_654_proj[1],x6G08_MD_654_proj[3],pos=4,label="654",col="red")

x6G08_MD_655_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_655.txt")
x6G08_MD_655_proj<-project.pca(x6G08_MD_655_raw,pcadata)
points(x6G08_MD_655_proj[1],x6G08_MD_655_proj[3],pch=20)
text(x6G08_MD_655_proj[1],x6G08_MD_655_proj[3],pos=4,label="655",col="green")

x6G08_MD_656_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_656.txt")
x6G08_MD_656_proj<-project.pca(x6G08_MD_656_raw,pcadata)
points(x6G08_MD_656_proj[1],x6G08_MD_656_proj[3],pch=20)
text(x6G08_MD_656_proj[1],x6G08_MD_656_proj[3],pos=4,label="656",col="green")

x6G08_MD_657_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_657.txt")
x6G08_MD_657_proj<-project.pca(x6G08_MD_657_raw,pcadata)
points(x6G08_MD_657_proj[1],x6G08_MD_657_proj[3],pch=20)
text(x6G08_MD_657_proj[1],x6G08_MD_657_proj[3],pos=4,label="657",col="blue")

x6G08_MD_658_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_658.txt")
x6G08_MD_658_proj<-project.pca(x6G08_MD_658_raw,pcadata)
points(x6G08_MD_658_proj[1],x6G08_MD_658_proj[3],pch=20)
text(x6G08_MD_658_proj[1],x6G08_MD_658_proj[3],pos=4,label="658",col="green")

x6G08_MD_659_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_659.txt")
x6G08_MD_659_proj<-project.pca(x6G08_MD_659_raw,pcadata)
points(x6G08_MD_659_proj[1],x6G08_MD_659_proj[3],pch=20)
text(x6G08_MD_659_proj[1],x6G08_MD_659_proj[3],pos=4,label="659",col="green")

x6G08_MD_660_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_660.txt")
x6G08_MD_660_proj<-project.pca(x6G08_MD_660_raw,pcadata)
points(x6G08_MD_660_proj[1],x6G08_MD_660_proj[3],pch=20)
text(x6G08_MD_660_proj[1],x6G08_MD_660_proj[3],pos=4,label="660",col="red")

x6G08_MD_661_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_661.txt")
x6G08_MD_661_proj<-project.pca(x6G08_MD_661_raw,pcadata)
points(x6G08_MD_661_proj[1],x6G08_MD_661_proj[3],pch=20)
text(x6G08_MD_661_proj[1],x6G08_MD_661_proj[3],pos=4,label="661",col="blue")

x6G08_MD_662_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_662.txt")
x6G08_MD_662_proj<-project.pca(x6G08_MD_662_raw,pcadata)
points(x6G08_MD_662_proj[1],x6G08_MD_662_proj[3],pch=20)
text(x6G08_MD_662_proj[1],x6G08_MD_662_proj[3],pos=4,label="662",col="blue")

x6G08_MD_663_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_663.txt")
x6G08_MD_663_proj<-project.pca(x6G08_MD_663_raw,pcadata)
points(x6G08_MD_663_proj[1],x6G08_MD_663_proj[3],pch=20)
text(x6G08_MD_663_proj[1],x6G08_MD_663_proj[3],pos=4,label="663",col="green")

x6G08_MD_664_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_664.txt")
x6G08_MD_664_proj<-project.pca(x6G08_MD_664_raw,pcadata)
points(x6G08_MD_664_proj[1],x6G08_MD_664_proj[3],pch=20)
text(x6G08_MD_664_proj[1],x6G08_MD_664_proj[3],pos=4,label="664",col="red")

x6G08_MD_665_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_665.txt")
x6G08_MD_665_proj<-project.pca(x6G08_MD_665_raw,pcadata)
points(x6G08_MD_665_proj[1],x6G08_MD_665_proj[3],pch=20)
text(x6G08_MD_665_proj[1],x6G08_MD_665_proj[3],pos=4,label="665",col="red")

x6G08_MD_666_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_666.txt")
x6G08_MD_666_proj<-project.pca(x6G08_MD_666_raw,pcadata)
points(x6G08_MD_666_proj[1],x6G08_MD_666_proj[3],pch=20)
text(x6G08_MD_666_proj[1],x6G08_MD_666_proj[3],pos=4,label="666",col="red")

x6G08_MD_667_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_667.txt")
x6G08_MD_667_proj<-project.pca(x6G08_MD_667_raw,pcadata)
points(x6G08_MD_667_proj[1],x6G08_MD_667_proj[3],pch=20)
text(x6G08_MD_667_proj[1],x6G08_MD_667_proj[3],pos=4,label="667",col="red")

x6G08_MD_668_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_668.txt")
x6G08_MD_668_proj<-project.pca(x6G08_MD_668_raw,pcadata)
points(x6G08_MD_668_proj[1],x6G08_MD_668_proj[3],pch=20)
text(x6G08_MD_668_proj[1],x6G08_MD_668_proj[3],pos=4,label="668",col="red")

x6G08_MD_669_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_669.txt")
x6G08_MD_669_proj<-project.pca(x6G08_MD_669_raw,pcadata)
points(x6G08_MD_669_proj[1],x6G08_MD_669_proj[3],pch=20)
text(x6G08_MD_669_proj[1],x6G08_MD_669_proj[3],pos=4,label="669",col="red")

x6G08_MD_670_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_670.txt")
x6G08_MD_670_proj<-project.pca(x6G08_MD_670_raw,pcadata)
points(x6G08_MD_670_proj[1],x6G08_MD_670_proj[3],pch=20)
text(x6G08_MD_670_proj[1],x6G08_MD_670_proj[3],pos=4,label="670",col="red")

x6G08_MD_671_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_671.txt")
x6G08_MD_671_proj<-project.pca(x6G08_MD_671_raw,pcadata)
points(x6G08_MD_671_proj[1],x6G08_MD_671_proj[3],pch=20)
text(x6G08_MD_671_proj[1],x6G08_MD_671_proj[3],pos=4,label="671",col="red")

x6G08_MD_672_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_672.txt")
x6G08_MD_672_proj<-project.pca(x6G08_MD_672_raw,pcadata)
points(x6G08_MD_672_proj[1],x6G08_MD_672_proj[3],pch=20)
text(x6G08_MD_672_proj[1],x6G08_MD_672_proj[3],pos=4,label="672",col="red")

x6G08_MD_673_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_673.txt")
x6G08_MD_673_proj<-project.pca(x6G08_MD_673_raw,pcadata)
points(x6G08_MD_673_proj[1],x6G08_MD_673_proj[3],pch=20)
text(x6G08_MD_673_proj[1],x6G08_MD_673_proj[3],pos=4,label="673",col="red")

x6G08_MD_674_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_674.txt")
x6G08_MD_674_proj<-project.pca(x6G08_MD_674_raw,pcadata)
points(x6G08_MD_674_proj[1],x6G08_MD_674_proj[3],pch=20)
text(x6G08_MD_674_proj[1],x6G08_MD_674_proj[3],pos=4,label="674",col="red")

x6G08_MD_675_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_675.txt")
x6G08_MD_675_proj<-project.pca(x6G08_MD_675_raw,pcadata)
points(x6G08_MD_675_proj[1],x6G08_MD_675_proj[3],pch=20)
text(x6G08_MD_675_proj[1],x6G08_MD_675_proj[3],pos=4,label="675",col="red")

x6G08_MD_676_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_676.txt")
x6G08_MD_676_proj<-project.pca(x6G08_MD_676_raw,pcadata)
points(x6G08_MD_676_proj[1],x6G08_MD_676_proj[3],pch=20)
text(x6G08_MD_676_proj[1],x6G08_MD_676_proj[3],pos=4,label="676",col="red")

x6G08_MD_677_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_677.txt")
x6G08_MD_677_proj<-project.pca(x6G08_MD_677_raw,pcadata)
points(x6G08_MD_677_proj[1],x6G08_MD_677_proj[3],pch=20)
text(x6G08_MD_677_proj[1],x6G08_MD_677_proj[3],pos=4,label="677",col="red")

x6G08_MD_678_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_678.txt")
x6G08_MD_678_proj<-project.pca(x6G08_MD_678_raw,pcadata)
points(x6G08_MD_678_proj[1],x6G08_MD_678_proj[3],pch=20)
text(x6G08_MD_678_proj[1],x6G08_MD_678_proj[3],pos=4,label="678",col="red")

x6G08_MD_679_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_679.txt")
x6G08_MD_679_proj<-project.pca(x6G08_MD_679_raw,pcadata)
points(x6G08_MD_679_proj[1],x6G08_MD_679_proj[3],pch=20)
text(x6G08_MD_679_proj[1],x6G08_MD_679_proj[3],pos=4,label="679",col="red")

x6G08_MD_680_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_680.txt")
x6G08_MD_680_proj<-project.pca(x6G08_MD_680_raw,pcadata)
points(x6G08_MD_680_proj[1],x6G08_MD_680_proj[3],pch=20)
text(x6G08_MD_680_proj[1],x6G08_MD_680_proj[3],pos=4,label="680",col="red")

x6G08_MD_681_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_681.txt")
x6G08_MD_681_proj<-project.pca(x6G08_MD_681_raw,pcadata)
points(x6G08_MD_681_proj[1],x6G08_MD_681_proj[3],pch=20)
text(x6G08_MD_681_proj[1],x6G08_MD_681_proj[3],pos=4,label="681",col="red")

x6G08_MD_682_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_682.txt")
x6G08_MD_682_proj<-project.pca(x6G08_MD_682_raw,pcadata)
points(x6G08_MD_682_proj[1],x6G08_MD_682_proj[3],pch=20)
text(x6G08_MD_682_proj[1],x6G08_MD_682_proj[3],pos=4,label="682",col="red")

x6G08_MD_683_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_683.txt")
x6G08_MD_683_proj<-project.pca(x6G08_MD_683_raw,pcadata)
points(x6G08_MD_683_proj[1],x6G08_MD_683_proj[3],pch=20)
text(x6G08_MD_683_proj[1],x6G08_MD_683_proj[3],pos=4,label="683",col="red")

x6G08_MD_684_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_684.txt")
x6G08_MD_684_proj<-project.pca(x6G08_MD_684_raw,pcadata)
points(x6G08_MD_684_proj[1],x6G08_MD_684_proj[3],pch=20)
text(x6G08_MD_684_proj[1],x6G08_MD_684_proj[3],pos=4,label="684",col="red")

x6G08_MD_685_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_685.txt")
x6G08_MD_685_proj<-project.pca(x6G08_MD_685_raw,pcadata)
points(x6G08_MD_685_proj[1],x6G08_MD_685_proj[3],pch=20)
text(x6G08_MD_685_proj[1],x6G08_MD_685_proj[3],pos=4,label="685",col="red")

x6G08_MD_686_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_686.txt")
x6G08_MD_686_proj<-project.pca(x6G08_MD_686_raw,pcadata)
points(x6G08_MD_686_proj[1],x6G08_MD_686_proj[3],pch=20)
text(x6G08_MD_686_proj[1],x6G08_MD_686_proj[3],pos=4,label="686",col="red")

x6G08_MD_687_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_687.txt")
x6G08_MD_687_proj<-project.pca(x6G08_MD_687_raw,pcadata)
points(x6G08_MD_687_proj[1],x6G08_MD_687_proj[3],pch=20)
text(x6G08_MD_687_proj[1],x6G08_MD_687_proj[3],pos=4,label="687",col="red")

x6G08_MD_688_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_688.txt")
x6G08_MD_688_proj<-project.pca(x6G08_MD_688_raw,pcadata)
points(x6G08_MD_688_proj[1],x6G08_MD_688_proj[3],pch=20)
text(x6G08_MD_688_proj[1],x6G08_MD_688_proj[3],pos=4,label="688",col="red")

x6G08_MD_689_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_689.txt")
x6G08_MD_689_proj<-project.pca(x6G08_MD_689_raw,pcadata)
points(x6G08_MD_689_proj[1],x6G08_MD_689_proj[3],pch=20)
text(x6G08_MD_689_proj[1],x6G08_MD_689_proj[3],pos=4,label="689",col="red")

x6G08_MD_690_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_690.txt")
x6G08_MD_690_proj<-project.pca(x6G08_MD_690_raw,pcadata)
points(x6G08_MD_690_proj[1],x6G08_MD_690_proj[3],pch=20)
text(x6G08_MD_690_proj[1],x6G08_MD_690_proj[3],pos=4,label="690",col="red")

x6G08_MD_691_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_691.txt")
x6G08_MD_691_proj<-project.pca(x6G08_MD_691_raw,pcadata)
points(x6G08_MD_691_proj[1],x6G08_MD_691_proj[3],pch=20)
text(x6G08_MD_691_proj[1],x6G08_MD_691_proj[3],pos=4,label="691",col="red")

x6G08_MD_692_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_692.txt")
x6G08_MD_692_proj<-project.pca(x6G08_MD_692_raw,pcadata)
points(x6G08_MD_692_proj[1],x6G08_MD_692_proj[3],pch=20)
text(x6G08_MD_692_proj[1],x6G08_MD_692_proj[3],pos=4,label="692",col="green")

x6G08_MD_693_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_693.txt")
x6G08_MD_693_proj<-project.pca(x6G08_MD_693_raw,pcadata)
points(x6G08_MD_693_proj[1],x6G08_MD_693_proj[3],pch=20)
text(x6G08_MD_693_proj[1],x6G08_MD_693_proj[3],pos=4,label="693",col="red")

x6G08_MD_694_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_694.txt")
x6G08_MD_694_proj<-project.pca(x6G08_MD_694_raw,pcadata)
points(x6G08_MD_694_proj[1],x6G08_MD_694_proj[3],pch=20)
text(x6G08_MD_694_proj[1],x6G08_MD_694_proj[3],pos=4,label="694",col="red")

x6G08_MD_695_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_695.txt")
x6G08_MD_695_proj<-project.pca(x6G08_MD_695_raw,pcadata)
points(x6G08_MD_695_proj[1],x6G08_MD_695_proj[3],pch=20)
text(x6G08_MD_695_proj[1],x6G08_MD_695_proj[3],pos=4,label="695",col="red")

x6G08_MD_696_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_696.txt")
x6G08_MD_696_proj<-project.pca(x6G08_MD_696_raw,pcadata)
points(x6G08_MD_696_proj[1],x6G08_MD_696_proj[3],pch=20)
text(x6G08_MD_696_proj[1],x6G08_MD_696_proj[3],pos=4,label="696",col="red")

x6G08_MD_697_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_697.txt")
x6G08_MD_697_proj<-project.pca(x6G08_MD_697_raw,pcadata)
points(x6G08_MD_697_proj[1],x6G08_MD_697_proj[3],pch=20)
text(x6G08_MD_697_proj[1],x6G08_MD_697_proj[3],pos=4,label="697",col="red")

x6G08_MD_698_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_698.txt")
x6G08_MD_698_proj<-project.pca(x6G08_MD_698_raw,pcadata)
points(x6G08_MD_698_proj[1],x6G08_MD_698_proj[3],pch=20)
text(x6G08_MD_698_proj[1],x6G08_MD_698_proj[3],pos=4,label="698",col="red")

x6G08_MD_699_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_699.txt")
x6G08_MD_699_proj<-project.pca(x6G08_MD_699_raw,pcadata)
points(x6G08_MD_699_proj[1],x6G08_MD_699_proj[3],pch=20)
text(x6G08_MD_699_proj[1],x6G08_MD_699_proj[3],pos=4,label="699",col="red")

x6G08_MD_700_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_700.txt")
x6G08_MD_700_proj<-project.pca(x6G08_MD_700_raw,pcadata)
points(x6G08_MD_700_proj[1],x6G08_MD_700_proj[3],pch=20)
text(x6G08_MD_700_proj[1],x6G08_MD_700_proj[3],pos=4,label="700",col="red")

x6G08_MD_701_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_701.txt")
x6G08_MD_701_proj<-project.pca(x6G08_MD_701_raw,pcadata)
points(x6G08_MD_701_proj[1],x6G08_MD_701_proj[3],pch=20)
text(x6G08_MD_701_proj[1],x6G08_MD_701_proj[3],pos=4,label="701",col="red")

x6G08_MD_702_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_702.txt")
x6G08_MD_702_proj<-project.pca(x6G08_MD_702_raw,pcadata)
points(x6G08_MD_702_proj[1],x6G08_MD_702_proj[3],pch=20)
text(x6G08_MD_702_proj[1],x6G08_MD_702_proj[3],pos=4,label="702",col="green")

x6G08_MD_703_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_703.txt")
x6G08_MD_703_proj<-project.pca(x6G08_MD_703_raw,pcadata)
points(x6G08_MD_703_proj[1],x6G08_MD_703_proj[3],pch=20)
text(x6G08_MD_703_proj[1],x6G08_MD_703_proj[3],pos=4,label="703",col="red")

x6G08_MD_704_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_704.txt")
x6G08_MD_704_proj<-project.pca(x6G08_MD_704_raw,pcadata)
points(x6G08_MD_704_proj[1],x6G08_MD_704_proj[3],pch=20)
text(x6G08_MD_704_proj[1],x6G08_MD_704_proj[3],pos=4,label="704",col="red")

x6G08_MD_705_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_705.txt")
x6G08_MD_705_proj<-project.pca(x6G08_MD_705_raw,pcadata)
points(x6G08_MD_705_proj[1],x6G08_MD_705_proj[3],pch=20)
text(x6G08_MD_705_proj[1],x6G08_MD_705_proj[3],pos=4,label="705",col="red")

x6G08_MD_706_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_706.txt")
x6G08_MD_706_proj<-project.pca(x6G08_MD_706_raw,pcadata)
points(x6G08_MD_706_proj[1],x6G08_MD_706_proj[3],pch=20)
text(x6G08_MD_706_proj[1],x6G08_MD_706_proj[3],pos=4,label="706",col="red")

x6G08_MD_707_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_707.txt")
x6G08_MD_707_proj<-project.pca(x6G08_MD_707_raw,pcadata)
points(x6G08_MD_707_proj[1],x6G08_MD_707_proj[3],pch=20)
text(x6G08_MD_707_proj[1],x6G08_MD_707_proj[3],pos=4,label="707",col="red")

x6G08_MD_708_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_708.txt")
x6G08_MD_708_proj<-project.pca(x6G08_MD_708_raw,pcadata)
points(x6G08_MD_708_proj[1],x6G08_MD_708_proj[3],pch=20)
text(x6G08_MD_708_proj[1],x6G08_MD_708_proj[3],pos=4,label="708",col="red")

x6G08_MD_709_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_709.txt")
x6G08_MD_709_proj<-project.pca(x6G08_MD_709_raw,pcadata)
points(x6G08_MD_709_proj[1],x6G08_MD_709_proj[3],pch=20)
text(x6G08_MD_709_proj[1],x6G08_MD_709_proj[3],pos=4,label="709",col="red")

x6G08_MD_710_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_710.txt")
x6G08_MD_710_proj<-project.pca(x6G08_MD_710_raw,pcadata)
points(x6G08_MD_710_proj[1],x6G08_MD_710_proj[3],pch=20)
text(x6G08_MD_710_proj[1],x6G08_MD_710_proj[3],pos=4,label="710",col="red")

x6G08_MD_711_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_711.txt")
x6G08_MD_711_proj<-project.pca(x6G08_MD_711_raw,pcadata)
points(x6G08_MD_711_proj[1],x6G08_MD_711_proj[3],pch=20)
text(x6G08_MD_711_proj[1],x6G08_MD_711_proj[3],pos=4,label="711",col="red")

x6G08_MD_712_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_712.txt")
x6G08_MD_712_proj<-project.pca(x6G08_MD_712_raw,pcadata)
points(x6G08_MD_712_proj[1],x6G08_MD_712_proj[3],pch=20)
text(x6G08_MD_712_proj[1],x6G08_MD_712_proj[3],pos=4,label="712",col="red")

x6G08_MD_713_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_713.txt")
x6G08_MD_713_proj<-project.pca(x6G08_MD_713_raw,pcadata)
points(x6G08_MD_713_proj[1],x6G08_MD_713_proj[3],pch=20)
text(x6G08_MD_713_proj[1],x6G08_MD_713_proj[3],pos=4,label="713",col="red")

x6G08_MD_714_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_714.txt")
x6G08_MD_714_proj<-project.pca(x6G08_MD_714_raw,pcadata)
points(x6G08_MD_714_proj[1],x6G08_MD_714_proj[3],pch=20)
text(x6G08_MD_714_proj[1],x6G08_MD_714_proj[3],pos=4,label="714",col="red")

x6G08_MD_715_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_715.txt")
x6G08_MD_715_proj<-project.pca(x6G08_MD_715_raw,pcadata)
points(x6G08_MD_715_proj[1],x6G08_MD_715_proj[3],pch=20)
text(x6G08_MD_715_proj[1],x6G08_MD_715_proj[3],pos=4,label="715",col="red")

x6G08_MD_716_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_716.txt")
x6G08_MD_716_proj<-project.pca(x6G08_MD_716_raw,pcadata)
points(x6G08_MD_716_proj[1],x6G08_MD_716_proj[3],pch=20)
text(x6G08_MD_716_proj[1],x6G08_MD_716_proj[3],pos=4,label="716",col="red")

x6G08_MD_717_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_717.txt")
x6G08_MD_717_proj<-project.pca(x6G08_MD_717_raw,pcadata)
points(x6G08_MD_717_proj[1],x6G08_MD_717_proj[3],pch=20)
text(x6G08_MD_717_proj[1],x6G08_MD_717_proj[3],pos=4,label="717",col="red")

x6G08_MD_718_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_718.txt")
x6G08_MD_718_proj<-project.pca(x6G08_MD_718_raw,pcadata)
points(x6G08_MD_718_proj[1],x6G08_MD_718_proj[3],pch=20)
text(x6G08_MD_718_proj[1],x6G08_MD_718_proj[3],pos=4,label="718",col="red")

x6G08_MD_719_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_719.txt")
x6G08_MD_719_proj<-project.pca(x6G08_MD_719_raw,pcadata)
points(x6G08_MD_719_proj[1],x6G08_MD_719_proj[3],pch=20)
text(x6G08_MD_719_proj[1],x6G08_MD_719_proj[3],pos=4,label="719",col="red")

x6G08_MD_720_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_720.txt")
x6G08_MD_720_proj<-project.pca(x6G08_MD_720_raw,pcadata)
points(x6G08_MD_720_proj[1],x6G08_MD_720_proj[3],pch=20)
text(x6G08_MD_720_proj[1],x6G08_MD_720_proj[3],pos=4,label="720",col="red")

x6G08_MD_721_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_721.txt")
x6G08_MD_721_proj<-project.pca(x6G08_MD_721_raw,pcadata)
points(x6G08_MD_721_proj[1],x6G08_MD_721_proj[3],pch=20)
text(x6G08_MD_721_proj[1],x6G08_MD_721_proj[3],pos=4,label="721",col="red")

x6G08_MD_722_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_722.txt")
x6G08_MD_722_proj<-project.pca(x6G08_MD_722_raw,pcadata)
points(x6G08_MD_722_proj[1],x6G08_MD_722_proj[3],pch=20)
text(x6G08_MD_722_proj[1],x6G08_MD_722_proj[3],pos=4,label="722",col="red")

x6G08_MD_723_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_723.txt")
x6G08_MD_723_proj<-project.pca(x6G08_MD_723_raw,pcadata)
points(x6G08_MD_723_proj[1],x6G08_MD_723_proj[3],pch=20)
text(x6G08_MD_723_proj[1],x6G08_MD_723_proj[3],pos=4,label="723",col="red")

x6G08_MD_724_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_724.txt")
x6G08_MD_724_proj<-project.pca(x6G08_MD_724_raw,pcadata)
points(x6G08_MD_724_proj[1],x6G08_MD_724_proj[3],pch=20)
text(x6G08_MD_724_proj[1],x6G08_MD_724_proj[3],pos=4,label="724",col="red")

x6G08_MD_725_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_725.txt")
x6G08_MD_725_proj<-project.pca(x6G08_MD_725_raw,pcadata)
points(x6G08_MD_725_proj[1],x6G08_MD_725_proj[3],pch=20)
text(x6G08_MD_725_proj[1],x6G08_MD_725_proj[3],pos=4,label="725",col="red")

x6G08_MD_726_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_726.txt")
x6G08_MD_726_proj<-project.pca(x6G08_MD_726_raw,pcadata)
points(x6G08_MD_726_proj[1],x6G08_MD_726_proj[3],pch=20)
text(x6G08_MD_726_proj[1],x6G08_MD_726_proj[3],pos=4,label="726",col="red")

x6G08_MD_727_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_727.txt")
x6G08_MD_727_proj<-project.pca(x6G08_MD_727_raw,pcadata)
points(x6G08_MD_727_proj[1],x6G08_MD_727_proj[3],pch=20)
text(x6G08_MD_727_proj[1],x6G08_MD_727_proj[3],pos=4,label="727",col="red")

x6G08_MD_728_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_728.txt")
x6G08_MD_728_proj<-project.pca(x6G08_MD_728_raw,pcadata)
points(x6G08_MD_728_proj[1],x6G08_MD_728_proj[3],pch=20)
text(x6G08_MD_728_proj[1],x6G08_MD_728_proj[3],pos=4,label="728",col="green")

x6G08_MD_729_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_729.txt")
x6G08_MD_729_proj<-project.pca(x6G08_MD_729_raw,pcadata)
points(x6G08_MD_729_proj[1],x6G08_MD_729_proj[3],pch=20)
text(x6G08_MD_729_proj[1],x6G08_MD_729_proj[3],pos=4,label="729",col="red")

x6G08_MD_730_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_730.txt")
x6G08_MD_730_proj<-project.pca(x6G08_MD_730_raw,pcadata)
points(x6G08_MD_730_proj[1],x6G08_MD_730_proj[3],pch=20)
text(x6G08_MD_730_proj[1],x6G08_MD_730_proj[3],pos=4,label="730",col="red")

x6G08_MD_731_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_731.txt")
x6G08_MD_731_proj<-project.pca(x6G08_MD_731_raw,pcadata)
points(x6G08_MD_731_proj[1],x6G08_MD_731_proj[3],pch=20)
text(x6G08_MD_731_proj[1],x6G08_MD_731_proj[3],pos=4,label="731",col="red")

x6G08_MD_732_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_732.txt")
x6G08_MD_732_proj<-project.pca(x6G08_MD_732_raw,pcadata)
points(x6G08_MD_732_proj[1],x6G08_MD_732_proj[3],pch=20)
text(x6G08_MD_732_proj[1],x6G08_MD_732_proj[3],pos=4,label="732",col="red")

x6G08_MD_733_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_733.txt")
x6G08_MD_733_proj<-project.pca(x6G08_MD_733_raw,pcadata)
points(x6G08_MD_733_proj[1],x6G08_MD_733_proj[3],pch=20)
text(x6G08_MD_733_proj[1],x6G08_MD_733_proj[3],pos=4,label="733",col="red")

x6G08_MD_734_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_734.txt")
x6G08_MD_734_proj<-project.pca(x6G08_MD_734_raw,pcadata)
points(x6G08_MD_734_proj[1],x6G08_MD_734_proj[3],pch=20)
text(x6G08_MD_734_proj[1],x6G08_MD_734_proj[3],pos=4,label="734",col="red")

x6G08_MD_735_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_735.txt")
x6G08_MD_735_proj<-project.pca(x6G08_MD_735_raw,pcadata)
points(x6G08_MD_735_proj[1],x6G08_MD_735_proj[3],pch=20)
text(x6G08_MD_735_proj[1],x6G08_MD_735_proj[3],pos=4,label="735",col="red")

x6G08_MD_736_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_736.txt")
x6G08_MD_736_proj<-project.pca(x6G08_MD_736_raw,pcadata)
points(x6G08_MD_736_proj[1],x6G08_MD_736_proj[3],pch=20)
text(x6G08_MD_736_proj[1],x6G08_MD_736_proj[3],pos=4,label="736",col="red")

x6G08_MD_737_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_737.txt")
x6G08_MD_737_proj<-project.pca(x6G08_MD_737_raw,pcadata)
points(x6G08_MD_737_proj[1],x6G08_MD_737_proj[3],pch=20)
text(x6G08_MD_737_proj[1],x6G08_MD_737_proj[3],pos=4,label="737",col="red")

x6G08_MD_738_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_738.txt")
x6G08_MD_738_proj<-project.pca(x6G08_MD_738_raw,pcadata)
points(x6G08_MD_738_proj[1],x6G08_MD_738_proj[3],pch=20)
text(x6G08_MD_738_proj[1],x6G08_MD_738_proj[3],pos=4,label="738",col="red")

x6G08_MD_739_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_739.txt")
x6G08_MD_739_proj<-project.pca(x6G08_MD_739_raw,pcadata)
points(x6G08_MD_739_proj[1],x6G08_MD_739_proj[3],pch=20)
text(x6G08_MD_739_proj[1],x6G08_MD_739_proj[3],pos=4,label="739",col="red")

x6G08_MD_740_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_740.txt")
x6G08_MD_740_proj<-project.pca(x6G08_MD_740_raw,pcadata)
points(x6G08_MD_740_proj[1],x6G08_MD_740_proj[3],pch=20)
text(x6G08_MD_740_proj[1],x6G08_MD_740_proj[3],pos=4,label="740",col="red")

x6G08_MD_741_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_741.txt")
x6G08_MD_741_proj<-project.pca(x6G08_MD_741_raw,pcadata)
points(x6G08_MD_741_proj[1],x6G08_MD_741_proj[3],pch=20)
text(x6G08_MD_741_proj[1],x6G08_MD_741_proj[3],pos=4,label="741",col="red")

x6G08_MD_742_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_742.txt")
x6G08_MD_742_proj<-project.pca(x6G08_MD_742_raw,pcadata)
points(x6G08_MD_742_proj[1],x6G08_MD_742_proj[3],pch=20)
text(x6G08_MD_742_proj[1],x6G08_MD_742_proj[3],pos=4,label="742",col="red")

x6G08_MD_743_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_743.txt")
x6G08_MD_743_proj<-project.pca(x6G08_MD_743_raw,pcadata)
points(x6G08_MD_743_proj[1],x6G08_MD_743_proj[3],pch=20)
text(x6G08_MD_743_proj[1],x6G08_MD_743_proj[3],pos=4,label="743",col="red")

x6G08_MD_744_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_744.txt")
x6G08_MD_744_proj<-project.pca(x6G08_MD_744_raw,pcadata)
points(x6G08_MD_744_proj[1],x6G08_MD_744_proj[3],pch=20)
text(x6G08_MD_744_proj[1],x6G08_MD_744_proj[3],pos=4,label="744",col="red")

x6G08_MD_745_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_745.txt")
x6G08_MD_745_proj<-project.pca(x6G08_MD_745_raw,pcadata)
points(x6G08_MD_745_proj[1],x6G08_MD_745_proj[3],pch=20)
text(x6G08_MD_745_proj[1],x6G08_MD_745_proj[3],pos=4,label="745",col="red")

x6G08_MD_746_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_746.txt")
x6G08_MD_746_proj<-project.pca(x6G08_MD_746_raw,pcadata)
points(x6G08_MD_746_proj[1],x6G08_MD_746_proj[3],pch=20)
text(x6G08_MD_746_proj[1],x6G08_MD_746_proj[3],pos=4,label="746",col="red")

x6G08_MD_747_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_747.txt")
x6G08_MD_747_proj<-project.pca(x6G08_MD_747_raw,pcadata)
points(x6G08_MD_747_proj[1],x6G08_MD_747_proj[3],pch=20)
text(x6G08_MD_747_proj[1],x6G08_MD_747_proj[3],pos=4,label="747",col="red")

x6G08_MD_748_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_748.txt")
x6G08_MD_748_proj<-project.pca(x6G08_MD_748_raw,pcadata)
points(x6G08_MD_748_proj[1],x6G08_MD_748_proj[3],pch=20)
text(x6G08_MD_748_proj[1],x6G08_MD_748_proj[3],pos=4,label="748",col="blue")

x6G08_MD_749_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_749.txt")
x6G08_MD_749_proj<-project.pca(x6G08_MD_749_raw,pcadata)
points(x6G08_MD_749_proj[1],x6G08_MD_749_proj[3],pch=20)
text(x6G08_MD_749_proj[1],x6G08_MD_749_proj[3],pos=4,label="749",col="red")

x6G08_MD_750_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_750.txt")
x6G08_MD_750_proj<-project.pca(x6G08_MD_750_raw,pcadata)
points(x6G08_MD_750_proj[1],x6G08_MD_750_proj[3],pch=20)
text(x6G08_MD_750_proj[1],x6G08_MD_750_proj[3],pos=4,label="750",col="red")

x6G08_MD_751_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_751.txt")
x6G08_MD_751_proj<-project.pca(x6G08_MD_751_raw,pcadata)
points(x6G08_MD_751_proj[1],x6G08_MD_751_proj[3],pch=20)
text(x6G08_MD_751_proj[1],x6G08_MD_751_proj[3],pos=4,label="751",col="red")

x6G08_MD_752_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_752.txt")
x6G08_MD_752_proj<-project.pca(x6G08_MD_752_raw,pcadata)
points(x6G08_MD_752_proj[1],x6G08_MD_752_proj[3],pch=20)
text(x6G08_MD_752_proj[1],x6G08_MD_752_proj[3],pos=4,label="752",col="red")

x6G08_MD_753_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_753.txt")
x6G08_MD_753_proj<-project.pca(x6G08_MD_753_raw,pcadata)
points(x6G08_MD_753_proj[1],x6G08_MD_753_proj[3],pch=20)
text(x6G08_MD_753_proj[1],x6G08_MD_753_proj[3],pos=4,label="753",col="red")

x6G08_MD_754_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_754.txt")
x6G08_MD_754_proj<-project.pca(x6G08_MD_754_raw,pcadata)
points(x6G08_MD_754_proj[1],x6G08_MD_754_proj[3],pch=20)
text(x6G08_MD_754_proj[1],x6G08_MD_754_proj[3],pos=4,label="754",col="green")

x6G08_MD_755_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_755.txt")
x6G08_MD_755_proj<-project.pca(x6G08_MD_755_raw,pcadata)
points(x6G08_MD_755_proj[1],x6G08_MD_755_proj[3],pch=20)
text(x6G08_MD_755_proj[1],x6G08_MD_755_proj[3],pos=4,label="755",col="red")

x6G08_MD_756_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_756.txt")
x6G08_MD_756_proj<-project.pca(x6G08_MD_756_raw,pcadata)
points(x6G08_MD_756_proj[1],x6G08_MD_756_proj[3],pch=20)
text(x6G08_MD_756_proj[1],x6G08_MD_756_proj[3],pos=4,label="756",col="red")

x6G08_MD_757_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_757.txt")
x6G08_MD_757_proj<-project.pca(x6G08_MD_757_raw,pcadata)
points(x6G08_MD_757_proj[1],x6G08_MD_757_proj[3],pch=20)
text(x6G08_MD_757_proj[1],x6G08_MD_757_proj[3],pos=4,label="757",col="red")

x6G08_MD_758_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_758.txt")
x6G08_MD_758_proj<-project.pca(x6G08_MD_758_raw,pcadata)
points(x6G08_MD_758_proj[1],x6G08_MD_758_proj[3],pch=20)
text(x6G08_MD_758_proj[1],x6G08_MD_758_proj[3],pos=4,label="758",col="red")

x6G08_MD_759_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_759.txt")
x6G08_MD_759_proj<-project.pca(x6G08_MD_759_raw,pcadata)
points(x6G08_MD_759_proj[1],x6G08_MD_759_proj[3],pch=20)
text(x6G08_MD_759_proj[1],x6G08_MD_759_proj[3],pos=4,label="759",col="red")

x6G08_MD_760_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_760.txt")
x6G08_MD_760_proj<-project.pca(x6G08_MD_760_raw,pcadata)
points(x6G08_MD_760_proj[1],x6G08_MD_760_proj[3],pch=20)
text(x6G08_MD_760_proj[1],x6G08_MD_760_proj[3],pos=4,label="760",col="red")

x6G08_MD_761_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_761.txt")
x6G08_MD_761_proj<-project.pca(x6G08_MD_761_raw,pcadata)
points(x6G08_MD_761_proj[1],x6G08_MD_761_proj[3],pch=20)
text(x6G08_MD_761_proj[1],x6G08_MD_761_proj[3],pos=4,label="761",col="red")

x6G08_MD_762_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_762.txt")
x6G08_MD_762_proj<-project.pca(x6G08_MD_762_raw,pcadata)
points(x6G08_MD_762_proj[1],x6G08_MD_762_proj[3],pch=20)
text(x6G08_MD_762_proj[1],x6G08_MD_762_proj[3],pos=4,label="762",col="red")

x6G08_MD_763_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_763.txt")
x6G08_MD_763_proj<-project.pca(x6G08_MD_763_raw,pcadata)
points(x6G08_MD_763_proj[1],x6G08_MD_763_proj[3],pch=20)
text(x6G08_MD_763_proj[1],x6G08_MD_763_proj[3],pos=4,label="763",col="red")

x6G08_MD_764_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_764.txt")
x6G08_MD_764_proj<-project.pca(x6G08_MD_764_raw,pcadata)
points(x6G08_MD_764_proj[1],x6G08_MD_764_proj[3],pch=20)
text(x6G08_MD_764_proj[1],x6G08_MD_764_proj[3],pos=4,label="764",col="green")

x6G08_MD_765_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_765.txt")
x6G08_MD_765_proj<-project.pca(x6G08_MD_765_raw,pcadata)
points(x6G08_MD_765_proj[1],x6G08_MD_765_proj[3],pch=20)
text(x6G08_MD_765_proj[1],x6G08_MD_765_proj[3],pos=4,label="765",col="red")

x6G08_MD_766_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_766.txt")
x6G08_MD_766_proj<-project.pca(x6G08_MD_766_raw,pcadata)
points(x6G08_MD_766_proj[1],x6G08_MD_766_proj[3],pch=20)
text(x6G08_MD_766_proj[1],x6G08_MD_766_proj[3],pos=4,label="766",col="red")

x6G08_MD_767_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_767.txt")
x6G08_MD_767_proj<-project.pca(x6G08_MD_767_raw,pcadata)
points(x6G08_MD_767_proj[1],x6G08_MD_767_proj[3],pch=20)
text(x6G08_MD_767_proj[1],x6G08_MD_767_proj[3],pos=4,label="767",col="red")

x6G08_MD_768_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_768.txt")
x6G08_MD_768_proj<-project.pca(x6G08_MD_768_raw,pcadata)
points(x6G08_MD_768_proj[1],x6G08_MD_768_proj[3],pch=20)
text(x6G08_MD_768_proj[1],x6G08_MD_768_proj[3],pos=4,label="768",col="red")

x6G08_MD_769_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_769.txt")
x6G08_MD_769_proj<-project.pca(x6G08_MD_769_raw,pcadata)
points(x6G08_MD_769_proj[1],x6G08_MD_769_proj[3],pch=20)
text(x6G08_MD_769_proj[1],x6G08_MD_769_proj[3],pos=4,label="769",col="red")

x6G08_MD_770_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_770.txt")
x6G08_MD_770_proj<-project.pca(x6G08_MD_770_raw,pcadata)
points(x6G08_MD_770_proj[1],x6G08_MD_770_proj[3],pch=20)
text(x6G08_MD_770_proj[1],x6G08_MD_770_proj[3],pos=4,label="770",col="red")

x6G08_MD_771_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_771.txt")
x6G08_MD_771_proj<-project.pca(x6G08_MD_771_raw,pcadata)
points(x6G08_MD_771_proj[1],x6G08_MD_771_proj[3],pch=20)
text(x6G08_MD_771_proj[1],x6G08_MD_771_proj[3],pos=4,label="771",col="red")

x6G08_MD_772_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_772.txt")
x6G08_MD_772_proj<-project.pca(x6G08_MD_772_raw,pcadata)
points(x6G08_MD_772_proj[1],x6G08_MD_772_proj[3],pch=20)
text(x6G08_MD_772_proj[1],x6G08_MD_772_proj[3],pos=4,label="772",col="red")

x6G08_MD_773_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_773.txt")
x6G08_MD_773_proj<-project.pca(x6G08_MD_773_raw,pcadata)
points(x6G08_MD_773_proj[1],x6G08_MD_773_proj[3],pch=20)
text(x6G08_MD_773_proj[1],x6G08_MD_773_proj[3],pos=4,label="773",col="red")

x6G08_MD_774_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_774.txt")
x6G08_MD_774_proj<-project.pca(x6G08_MD_774_raw,pcadata)
points(x6G08_MD_774_proj[1],x6G08_MD_774_proj[3],pch=20)
text(x6G08_MD_774_proj[1],x6G08_MD_774_proj[3],pos=4,label="774",col="red")

x6G08_MD_775_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_775.txt")
x6G08_MD_775_proj<-project.pca(x6G08_MD_775_raw,pcadata)
points(x6G08_MD_775_proj[1],x6G08_MD_775_proj[3],pch=20)
text(x6G08_MD_775_proj[1],x6G08_MD_775_proj[3],pos=4,label="775",col="red")

x6G08_MD_776_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_776.txt")
x6G08_MD_776_proj<-project.pca(x6G08_MD_776_raw,pcadata)
points(x6G08_MD_776_proj[1],x6G08_MD_776_proj[3],pch=20)
text(x6G08_MD_776_proj[1],x6G08_MD_776_proj[3],pos=4,label="776",col="red")

x6G08_MD_777_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_777.txt")
x6G08_MD_777_proj<-project.pca(x6G08_MD_777_raw,pcadata)
points(x6G08_MD_777_proj[1],x6G08_MD_777_proj[3],pch=20)
text(x6G08_MD_777_proj[1],x6G08_MD_777_proj[3],pos=4,label="777",col="green")

x6G08_MD_778_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_778.txt")
x6G08_MD_778_proj<-project.pca(x6G08_MD_778_raw,pcadata)
points(x6G08_MD_778_proj[1],x6G08_MD_778_proj[3],pch=20)
text(x6G08_MD_778_proj[1],x6G08_MD_778_proj[3],pos=4,label="778",col="green")

x6G08_MD_779_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_779.txt")
x6G08_MD_779_proj<-project.pca(x6G08_MD_779_raw,pcadata)
points(x6G08_MD_779_proj[1],x6G08_MD_779_proj[3],pch=20)
text(x6G08_MD_779_proj[1],x6G08_MD_779_proj[3],pos=4,label="779",col="green")

x6G08_MD_780_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_780.txt")
x6G08_MD_780_proj<-project.pca(x6G08_MD_780_raw,pcadata)
points(x6G08_MD_780_proj[1],x6G08_MD_780_proj[3],pch=20)
text(x6G08_MD_780_proj[1],x6G08_MD_780_proj[3],pos=4,label="780",col="red")

x6G08_MD_781_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_781.txt")
x6G08_MD_781_proj<-project.pca(x6G08_MD_781_raw,pcadata)
points(x6G08_MD_781_proj[1],x6G08_MD_781_proj[3],pch=20)
text(x6G08_MD_781_proj[1],x6G08_MD_781_proj[3],pos=4,label="781",col="red")

x6G08_MD_782_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_782.txt")
x6G08_MD_782_proj<-project.pca(x6G08_MD_782_raw,pcadata)
points(x6G08_MD_782_proj[1],x6G08_MD_782_proj[3],pch=20)
text(x6G08_MD_782_proj[1],x6G08_MD_782_proj[3],pos=4,label="782",col="red")

x6G08_MD_783_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_783.txt")
x6G08_MD_783_proj<-project.pca(x6G08_MD_783_raw,pcadata)
points(x6G08_MD_783_proj[1],x6G08_MD_783_proj[3],pch=20)
text(x6G08_MD_783_proj[1],x6G08_MD_783_proj[3],pos=4,label="783",col="green")

x6G08_MD_784_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_784.txt")
x6G08_MD_784_proj<-project.pca(x6G08_MD_784_raw,pcadata)
points(x6G08_MD_784_proj[1],x6G08_MD_784_proj[3],pch=20)
text(x6G08_MD_784_proj[1],x6G08_MD_784_proj[3],pos=4,label="784",col="red")

x6G08_MD_785_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_785.txt")
x6G08_MD_785_proj<-project.pca(x6G08_MD_785_raw,pcadata)
points(x6G08_MD_785_proj[1],x6G08_MD_785_proj[3],pch=20)
text(x6G08_MD_785_proj[1],x6G08_MD_785_proj[3],pos=4,label="785",col="red")

x6G08_MD_786_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_786.txt")
x6G08_MD_786_proj<-project.pca(x6G08_MD_786_raw,pcadata)
points(x6G08_MD_786_proj[1],x6G08_MD_786_proj[3],pch=20)
text(x6G08_MD_786_proj[1],x6G08_MD_786_proj[3],pos=4,label="786",col="red")

x6G08_MD_787_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_787.txt")
x6G08_MD_787_proj<-project.pca(x6G08_MD_787_raw,pcadata)
points(x6G08_MD_787_proj[1],x6G08_MD_787_proj[3],pch=20)
text(x6G08_MD_787_proj[1],x6G08_MD_787_proj[3],pos=4,label="787",col="red")

x6G08_MD_788_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_788.txt")
x6G08_MD_788_proj<-project.pca(x6G08_MD_788_raw,pcadata)
points(x6G08_MD_788_proj[1],x6G08_MD_788_proj[3],pch=20)
text(x6G08_MD_788_proj[1],x6G08_MD_788_proj[3],pos=4,label="788",col="red")

x6G08_MD_789_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_789.txt")
x6G08_MD_789_proj<-project.pca(x6G08_MD_789_raw,pcadata)
points(x6G08_MD_789_proj[1],x6G08_MD_789_proj[3],pch=20)
text(x6G08_MD_789_proj[1],x6G08_MD_789_proj[3],pos=4,label="789",col="red")

x6G08_MD_790_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_790.txt")
x6G08_MD_790_proj<-project.pca(x6G08_MD_790_raw,pcadata)
points(x6G08_MD_790_proj[1],x6G08_MD_790_proj[3],pch=20)
text(x6G08_MD_790_proj[1],x6G08_MD_790_proj[3],pos=4,label="790",col="red")

x6G08_MD_791_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_791.txt")
x6G08_MD_791_proj<-project.pca(x6G08_MD_791_raw,pcadata)
points(x6G08_MD_791_proj[1],x6G08_MD_791_proj[3],pch=20)
text(x6G08_MD_791_proj[1],x6G08_MD_791_proj[3],pos=4,label="791",col="red")

x6G08_MD_792_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_792.txt")
x6G08_MD_792_proj<-project.pca(x6G08_MD_792_raw,pcadata)
points(x6G08_MD_792_proj[1],x6G08_MD_792_proj[3],pch=20)
text(x6G08_MD_792_proj[1],x6G08_MD_792_proj[3],pos=4,label="792",col="red")

x6G08_MD_793_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_793.txt")
x6G08_MD_793_proj<-project.pca(x6G08_MD_793_raw,pcadata)
points(x6G08_MD_793_proj[1],x6G08_MD_793_proj[3],pch=20)
text(x6G08_MD_793_proj[1],x6G08_MD_793_proj[3],pos=4,label="793",col="red")

x6G08_MD_794_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_794.txt")
x6G08_MD_794_proj<-project.pca(x6G08_MD_794_raw,pcadata)
points(x6G08_MD_794_proj[1],x6G08_MD_794_proj[3],pch=20)
text(x6G08_MD_794_proj[1],x6G08_MD_794_proj[3],pos=4,label="794",col="red")

x6G08_MD_795_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_795.txt")
x6G08_MD_795_proj<-project.pca(x6G08_MD_795_raw,pcadata)
points(x6G08_MD_795_proj[1],x6G08_MD_795_proj[3],pch=20)
text(x6G08_MD_795_proj[1],x6G08_MD_795_proj[3],pos=4,label="795",col="red")

x6G08_MD_796_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_796.txt")
x6G08_MD_796_proj<-project.pca(x6G08_MD_796_raw,pcadata)
points(x6G08_MD_796_proj[1],x6G08_MD_796_proj[3],pch=20)
text(x6G08_MD_796_proj[1],x6G08_MD_796_proj[3],pos=4,label="796",col="red")

x6G08_MD_797_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_797.txt")
x6G08_MD_797_proj<-project.pca(x6G08_MD_797_raw,pcadata)
points(x6G08_MD_797_proj[1],x6G08_MD_797_proj[3],pch=20)
text(x6G08_MD_797_proj[1],x6G08_MD_797_proj[3],pos=4,label="797",col="red")

x6G08_MD_798_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_798.txt")
x6G08_MD_798_proj<-project.pca(x6G08_MD_798_raw,pcadata)
points(x6G08_MD_798_proj[1],x6G08_MD_798_proj[3],pch=20)
text(x6G08_MD_798_proj[1],x6G08_MD_798_proj[3],pos=4,label="798",col="red")

x6G08_MD_799_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_799.txt")
x6G08_MD_799_proj<-project.pca(x6G08_MD_799_raw,pcadata)
points(x6G08_MD_799_proj[1],x6G08_MD_799_proj[3],pch=20)
text(x6G08_MD_799_proj[1],x6G08_MD_799_proj[3],pos=4,label="799",col="red")

x6G08_MD_800_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_800.txt")
x6G08_MD_800_proj<-project.pca(x6G08_MD_800_raw,pcadata)
points(x6G08_MD_800_proj[1],x6G08_MD_800_proj[3],pch=20)
text(x6G08_MD_800_proj[1],x6G08_MD_800_proj[3],pos=4,label="800",col="red")

x6G08_MD_801_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_801.txt")
x6G08_MD_801_proj<-project.pca(x6G08_MD_801_raw,pcadata)
points(x6G08_MD_801_proj[1],x6G08_MD_801_proj[3],pch=20)
text(x6G08_MD_801_proj[1],x6G08_MD_801_proj[3],pos=4,label="801",col="red")

x6G08_MD_802_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_802.txt")
x6G08_MD_802_proj<-project.pca(x6G08_MD_802_raw,pcadata)
points(x6G08_MD_802_proj[1],x6G08_MD_802_proj[3],pch=20)
text(x6G08_MD_802_proj[1],x6G08_MD_802_proj[3],pos=4,label="802",col="red")

x6G08_MD_803_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_803.txt")
x6G08_MD_803_proj<-project.pca(x6G08_MD_803_raw,pcadata)
points(x6G08_MD_803_proj[1],x6G08_MD_803_proj[3],pch=20)
text(x6G08_MD_803_proj[1],x6G08_MD_803_proj[3],pos=4,label="803",col="red")

x6G08_MD_804_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_804.txt")
x6G08_MD_804_proj<-project.pca(x6G08_MD_804_raw,pcadata)
points(x6G08_MD_804_proj[1],x6G08_MD_804_proj[3],pch=20)
text(x6G08_MD_804_proj[1],x6G08_MD_804_proj[3],pos=4,label="804",col="red")

x6G08_MD_805_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_805.txt")
x6G08_MD_805_proj<-project.pca(x6G08_MD_805_raw,pcadata)
points(x6G08_MD_805_proj[1],x6G08_MD_805_proj[3],pch=20)
text(x6G08_MD_805_proj[1],x6G08_MD_805_proj[3],pos=4,label="805",col="red")

x6G08_MD_806_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_806.txt")
x6G08_MD_806_proj<-project.pca(x6G08_MD_806_raw,pcadata)
points(x6G08_MD_806_proj[1],x6G08_MD_806_proj[3],pch=20)
text(x6G08_MD_806_proj[1],x6G08_MD_806_proj[3],pos=4,label="806",col="red")

x6G08_MD_807_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_807.txt")
x6G08_MD_807_proj<-project.pca(x6G08_MD_807_raw,pcadata)
points(x6G08_MD_807_proj[1],x6G08_MD_807_proj[3],pch=20)
text(x6G08_MD_807_proj[1],x6G08_MD_807_proj[3],pos=4,label="807",col="red")

x6G08_MD_808_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_808.txt")
x6G08_MD_808_proj<-project.pca(x6G08_MD_808_raw,pcadata)
points(x6G08_MD_808_proj[1],x6G08_MD_808_proj[3],pch=20)
text(x6G08_MD_808_proj[1],x6G08_MD_808_proj[3],pos=4,label="808",col="red")

x6G08_MD_809_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_809.txt")
x6G08_MD_809_proj<-project.pca(x6G08_MD_809_raw,pcadata)
points(x6G08_MD_809_proj[1],x6G08_MD_809_proj[3],pch=20)
text(x6G08_MD_809_proj[1],x6G08_MD_809_proj[3],pos=4,label="809",col="red")

x6G08_MD_810_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_810.txt")
x6G08_MD_810_proj<-project.pca(x6G08_MD_810_raw,pcadata)
points(x6G08_MD_810_proj[1],x6G08_MD_810_proj[3],pch=20)
text(x6G08_MD_810_proj[1],x6G08_MD_810_proj[3],pos=4,label="810",col="red")

x6G08_MD_811_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_811.txt")
x6G08_MD_811_proj<-project.pca(x6G08_MD_811_raw,pcadata)
points(x6G08_MD_811_proj[1],x6G08_MD_811_proj[3],pch=20)
text(x6G08_MD_811_proj[1],x6G08_MD_811_proj[3],pos=4,label="811",col="red")

x6G08_MD_812_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_812.txt")
x6G08_MD_812_proj<-project.pca(x6G08_MD_812_raw,pcadata)
points(x6G08_MD_812_proj[1],x6G08_MD_812_proj[3],pch=20)
text(x6G08_MD_812_proj[1],x6G08_MD_812_proj[3],pos=4,label="812",col="red")

x6G08_MD_813_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_813.txt")
x6G08_MD_813_proj<-project.pca(x6G08_MD_813_raw,pcadata)
points(x6G08_MD_813_proj[1],x6G08_MD_813_proj[3],pch=20)
text(x6G08_MD_813_proj[1],x6G08_MD_813_proj[3],pos=4,label="813",col="red")

x6G08_MD_814_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_814.txt")
x6G08_MD_814_proj<-project.pca(x6G08_MD_814_raw,pcadata)
points(x6G08_MD_814_proj[1],x6G08_MD_814_proj[3],pch=20)
text(x6G08_MD_814_proj[1],x6G08_MD_814_proj[3],pos=4,label="814",col="red")

x6G08_MD_815_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_815.txt")
x6G08_MD_815_proj<-project.pca(x6G08_MD_815_raw,pcadata)
points(x6G08_MD_815_proj[1],x6G08_MD_815_proj[3],pch=20)
text(x6G08_MD_815_proj[1],x6G08_MD_815_proj[3],pos=4,label="815",col="red")

x6G08_MD_816_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_816.txt")
x6G08_MD_816_proj<-project.pca(x6G08_MD_816_raw,pcadata)
points(x6G08_MD_816_proj[1],x6G08_MD_816_proj[3],pch=20)
text(x6G08_MD_816_proj[1],x6G08_MD_816_proj[3],pos=4,label="816",col="blue")

x6G08_MD_817_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_817.txt")
x6G08_MD_817_proj<-project.pca(x6G08_MD_817_raw,pcadata)
points(x6G08_MD_817_proj[1],x6G08_MD_817_proj[3],pch=20)
text(x6G08_MD_817_proj[1],x6G08_MD_817_proj[3],pos=4,label="817",col="red")

x6G08_MD_818_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_818.txt")
x6G08_MD_818_proj<-project.pca(x6G08_MD_818_raw,pcadata)
points(x6G08_MD_818_proj[1],x6G08_MD_818_proj[3],pch=20)
text(x6G08_MD_818_proj[1],x6G08_MD_818_proj[3],pos=4,label="818",col="blue")

x6G08_MD_819_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_819.txt")
x6G08_MD_819_proj<-project.pca(x6G08_MD_819_raw,pcadata)
points(x6G08_MD_819_proj[1],x6G08_MD_819_proj[3],pch=20)
text(x6G08_MD_819_proj[1],x6G08_MD_819_proj[3],pos=4,label="819",col="blue")

x6G08_MD_820_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_820.txt")
x6G08_MD_820_proj<-project.pca(x6G08_MD_820_raw,pcadata)
points(x6G08_MD_820_proj[1],x6G08_MD_820_proj[3],pch=20)
text(x6G08_MD_820_proj[1],x6G08_MD_820_proj[3],pos=4,label="820",col="blue")

x6G08_MD_821_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_821.txt")
x6G08_MD_821_proj<-project.pca(x6G08_MD_821_raw,pcadata)
points(x6G08_MD_821_proj[1],x6G08_MD_821_proj[3],pch=20)
text(x6G08_MD_821_proj[1],x6G08_MD_821_proj[3],pos=4,label="821",col="red")

x6G08_MD_822_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_822.txt")
x6G08_MD_822_proj<-project.pca(x6G08_MD_822_raw,pcadata)
points(x6G08_MD_822_proj[1],x6G08_MD_822_proj[3],pch=20)
text(x6G08_MD_822_proj[1],x6G08_MD_822_proj[3],pos=4,label="822",col="red")

x6G08_MD_823_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_823.txt")
x6G08_MD_823_proj<-project.pca(x6G08_MD_823_raw,pcadata)
points(x6G08_MD_823_proj[1],x6G08_MD_823_proj[3],pch=20)
text(x6G08_MD_823_proj[1],x6G08_MD_823_proj[3],pos=4,label="823",col="red")

x6G08_MD_824_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_824.txt")
x6G08_MD_824_proj<-project.pca(x6G08_MD_824_raw,pcadata)
points(x6G08_MD_824_proj[1],x6G08_MD_824_proj[3],pch=20)
text(x6G08_MD_824_proj[1],x6G08_MD_824_proj[3],pos=4,label="824",col="red")

x6G08_MD_825_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_825.txt")
x6G08_MD_825_proj<-project.pca(x6G08_MD_825_raw,pcadata)
points(x6G08_MD_825_proj[1],x6G08_MD_825_proj[3],pch=20)
text(x6G08_MD_825_proj[1],x6G08_MD_825_proj[3],pos=4,label="825",col="red")

x6G08_MD_826_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_826.txt")
x6G08_MD_826_proj<-project.pca(x6G08_MD_826_raw,pcadata)
points(x6G08_MD_826_proj[1],x6G08_MD_826_proj[3],pch=20)
text(x6G08_MD_826_proj[1],x6G08_MD_826_proj[3],pos=4,label="826",col="red")

x6G08_MD_827_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_827.txt")
x6G08_MD_827_proj<-project.pca(x6G08_MD_827_raw,pcadata)
points(x6G08_MD_827_proj[1],x6G08_MD_827_proj[3],pch=20)
text(x6G08_MD_827_proj[1],x6G08_MD_827_proj[3],pos=4,label="827",col="red")

x6G08_MD_828_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_828.txt")
x6G08_MD_828_proj<-project.pca(x6G08_MD_828_raw,pcadata)
points(x6G08_MD_828_proj[1],x6G08_MD_828_proj[3],pch=20)
text(x6G08_MD_828_proj[1],x6G08_MD_828_proj[3],pos=4,label="828",col="red")

x6G08_MD_829_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_829.txt")
x6G08_MD_829_proj<-project.pca(x6G08_MD_829_raw,pcadata)
points(x6G08_MD_829_proj[1],x6G08_MD_829_proj[3],pch=20)
text(x6G08_MD_829_proj[1],x6G08_MD_829_proj[3],pos=4,label="829",col="red")

x6G08_MD_830_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_830.txt")
x6G08_MD_830_proj<-project.pca(x6G08_MD_830_raw,pcadata)
points(x6G08_MD_830_proj[1],x6G08_MD_830_proj[3],pch=20)
text(x6G08_MD_830_proj[1],x6G08_MD_830_proj[3],pos=4,label="830",col="red")

x6G08_MD_831_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_831.txt")
x6G08_MD_831_proj<-project.pca(x6G08_MD_831_raw,pcadata)
points(x6G08_MD_831_proj[1],x6G08_MD_831_proj[3],pch=20)
text(x6G08_MD_831_proj[1],x6G08_MD_831_proj[3],pos=4,label="831",col="red")

x6G08_MD_832_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_832.txt")
x6G08_MD_832_proj<-project.pca(x6G08_MD_832_raw,pcadata)
points(x6G08_MD_832_proj[1],x6G08_MD_832_proj[3],pch=20)
text(x6G08_MD_832_proj[1],x6G08_MD_832_proj[3],pos=4,label="832",col="red")

x6G08_MD_833_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_833.txt")
x6G08_MD_833_proj<-project.pca(x6G08_MD_833_raw,pcadata)
points(x6G08_MD_833_proj[1],x6G08_MD_833_proj[3],pch=20)
text(x6G08_MD_833_proj[1],x6G08_MD_833_proj[3],pos=4,label="833",col="red")

x6G08_MD_834_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_834.txt")
x6G08_MD_834_proj<-project.pca(x6G08_MD_834_raw,pcadata)
points(x6G08_MD_834_proj[1],x6G08_MD_834_proj[3],pch=20)
text(x6G08_MD_834_proj[1],x6G08_MD_834_proj[3],pos=4,label="834",col="red")

x6G08_MD_835_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_835.txt")
x6G08_MD_835_proj<-project.pca(x6G08_MD_835_raw,pcadata)
points(x6G08_MD_835_proj[1],x6G08_MD_835_proj[3],pch=20)
text(x6G08_MD_835_proj[1],x6G08_MD_835_proj[3],pos=4,label="835",col="red")

x6G08_MD_836_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_836.txt")
x6G08_MD_836_proj<-project.pca(x6G08_MD_836_raw,pcadata)
points(x6G08_MD_836_proj[1],x6G08_MD_836_proj[3],pch=20)
text(x6G08_MD_836_proj[1],x6G08_MD_836_proj[3],pos=4,label="836",col="red")

x6G08_MD_837_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_837.txt")
x6G08_MD_837_proj<-project.pca(x6G08_MD_837_raw,pcadata)
points(x6G08_MD_837_proj[1],x6G08_MD_837_proj[3],pch=20)
text(x6G08_MD_837_proj[1],x6G08_MD_837_proj[3],pos=4,label="837",col="red")

x6G08_MD_838_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_838.txt")
x6G08_MD_838_proj<-project.pca(x6G08_MD_838_raw,pcadata)
points(x6G08_MD_838_proj[1],x6G08_MD_838_proj[3],pch=20)
text(x6G08_MD_838_proj[1],x6G08_MD_838_proj[3],pos=4,label="838",col="red")

x6G08_MD_839_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_839.txt")
x6G08_MD_839_proj<-project.pca(x6G08_MD_839_raw,pcadata)
points(x6G08_MD_839_proj[1],x6G08_MD_839_proj[3],pch=20)
text(x6G08_MD_839_proj[1],x6G08_MD_839_proj[3],pos=4,label="839",col="red")

x6G08_MD_840_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_840.txt")
x6G08_MD_840_proj<-project.pca(x6G08_MD_840_raw,pcadata)
points(x6G08_MD_840_proj[1],x6G08_MD_840_proj[3],pch=20)
text(x6G08_MD_840_proj[1],x6G08_MD_840_proj[3],pos=4,label="840",col="red")

x6G08_MD_841_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_841.txt")
x6G08_MD_841_proj<-project.pca(x6G08_MD_841_raw,pcadata)
points(x6G08_MD_841_proj[1],x6G08_MD_841_proj[3],pch=20)
text(x6G08_MD_841_proj[1],x6G08_MD_841_proj[3],pos=4,label="841",col="red")

x6G08_MD_842_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_842.txt")
x6G08_MD_842_proj<-project.pca(x6G08_MD_842_raw,pcadata)
points(x6G08_MD_842_proj[1],x6G08_MD_842_proj[3],pch=20)
text(x6G08_MD_842_proj[1],x6G08_MD_842_proj[3],pos=4,label="842",col="green")

x6G08_MD_843_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_843.txt")
x6G08_MD_843_proj<-project.pca(x6G08_MD_843_raw,pcadata)
points(x6G08_MD_843_proj[1],x6G08_MD_843_proj[3],pch=20)
text(x6G08_MD_843_proj[1],x6G08_MD_843_proj[3],pos=4,label="843",col="red")

x6G08_MD_844_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_844.txt")
x6G08_MD_844_proj<-project.pca(x6G08_MD_844_raw,pcadata)
points(x6G08_MD_844_proj[1],x6G08_MD_844_proj[3],pch=20)
text(x6G08_MD_844_proj[1],x6G08_MD_844_proj[3],pos=4,label="844",col="red")

x6G08_MD_845_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_845.txt")
x6G08_MD_845_proj<-project.pca(x6G08_MD_845_raw,pcadata)
points(x6G08_MD_845_proj[1],x6G08_MD_845_proj[3],pch=20)
text(x6G08_MD_845_proj[1],x6G08_MD_845_proj[3],pos=4,label="845",col="red")

x6G08_MD_846_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_846.txt")
x6G08_MD_846_proj<-project.pca(x6G08_MD_846_raw,pcadata)
points(x6G08_MD_846_proj[1],x6G08_MD_846_proj[3],pch=20)
text(x6G08_MD_846_proj[1],x6G08_MD_846_proj[3],pos=4,label="846",col="red")

x6G08_MD_847_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_847.txt")
x6G08_MD_847_proj<-project.pca(x6G08_MD_847_raw,pcadata)
points(x6G08_MD_847_proj[1],x6G08_MD_847_proj[3],pch=20)
text(x6G08_MD_847_proj[1],x6G08_MD_847_proj[3],pos=4,label="847",col="red")

x6G08_MD_848_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_848.txt")
x6G08_MD_848_proj<-project.pca(x6G08_MD_848_raw,pcadata)
points(x6G08_MD_848_proj[1],x6G08_MD_848_proj[3],pch=20)
text(x6G08_MD_848_proj[1],x6G08_MD_848_proj[3],pos=4,label="848",col="red")

x6G08_MD_849_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_849.txt")
x6G08_MD_849_proj<-project.pca(x6G08_MD_849_raw,pcadata)
points(x6G08_MD_849_proj[1],x6G08_MD_849_proj[3],pch=20)
text(x6G08_MD_849_proj[1],x6G08_MD_849_proj[3],pos=4,label="849",col="red")

x6G08_MD_850_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_850.txt")
x6G08_MD_850_proj<-project.pca(x6G08_MD_850_raw,pcadata)
points(x6G08_MD_850_proj[1],x6G08_MD_850_proj[3],pch=20)
text(x6G08_MD_850_proj[1],x6G08_MD_850_proj[3],pos=4,label="850",col="red")

x6G08_MD_851_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_851.txt")
x6G08_MD_851_proj<-project.pca(x6G08_MD_851_raw,pcadata)
points(x6G08_MD_851_proj[1],x6G08_MD_851_proj[3],pch=20)
text(x6G08_MD_851_proj[1],x6G08_MD_851_proj[3],pos=4,label="851",col="red")

x6G08_MD_852_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_852.txt")
x6G08_MD_852_proj<-project.pca(x6G08_MD_852_raw,pcadata)
points(x6G08_MD_852_proj[1],x6G08_MD_852_proj[3],pch=20)
text(x6G08_MD_852_proj[1],x6G08_MD_852_proj[3],pos=4,label="852",col="red")

x6G08_MD_853_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_853.txt")
x6G08_MD_853_proj<-project.pca(x6G08_MD_853_raw,pcadata)
points(x6G08_MD_853_proj[1],x6G08_MD_853_proj[3],pch=20)
text(x6G08_MD_853_proj[1],x6G08_MD_853_proj[3],pos=4,label="853",col="red")

x6G08_MD_854_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_854.txt")
x6G08_MD_854_proj<-project.pca(x6G08_MD_854_raw,pcadata)
points(x6G08_MD_854_proj[1],x6G08_MD_854_proj[3],pch=20)
text(x6G08_MD_854_proj[1],x6G08_MD_854_proj[3],pos=4,label="854",col="red")

x6G08_MD_855_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_855.txt")
x6G08_MD_855_proj<-project.pca(x6G08_MD_855_raw,pcadata)
points(x6G08_MD_855_proj[1],x6G08_MD_855_proj[3],pch=20)
text(x6G08_MD_855_proj[1],x6G08_MD_855_proj[3],pos=4,label="855",col="green")

x6G08_MD_856_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_856.txt")
x6G08_MD_856_proj<-project.pca(x6G08_MD_856_raw,pcadata)
points(x6G08_MD_856_proj[1],x6G08_MD_856_proj[3],pch=20)
text(x6G08_MD_856_proj[1],x6G08_MD_856_proj[3],pos=4,label="856",col="red")

x6G08_MD_857_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_857.txt")
x6G08_MD_857_proj<-project.pca(x6G08_MD_857_raw,pcadata)
points(x6G08_MD_857_proj[1],x6G08_MD_857_proj[3],pch=20)
text(x6G08_MD_857_proj[1],x6G08_MD_857_proj[3],pos=4,label="857",col="red")

x6G08_MD_858_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_858.txt")
x6G08_MD_858_proj<-project.pca(x6G08_MD_858_raw,pcadata)
points(x6G08_MD_858_proj[1],x6G08_MD_858_proj[3],pch=20)
text(x6G08_MD_858_proj[1],x6G08_MD_858_proj[3],pos=4,label="858",col="red")

x6G08_MD_859_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_859.txt")
x6G08_MD_859_proj<-project.pca(x6G08_MD_859_raw,pcadata)
points(x6G08_MD_859_proj[1],x6G08_MD_859_proj[3],pch=20)
text(x6G08_MD_859_proj[1],x6G08_MD_859_proj[3],pos=4,label="859",col="red")

x6G08_MD_860_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_860.txt")
x6G08_MD_860_proj<-project.pca(x6G08_MD_860_raw,pcadata)
points(x6G08_MD_860_proj[1],x6G08_MD_860_proj[3],pch=20)
text(x6G08_MD_860_proj[1],x6G08_MD_860_proj[3],pos=4,label="860",col="red")

x6G08_MD_861_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_861.txt")
x6G08_MD_861_proj<-project.pca(x6G08_MD_861_raw,pcadata)
points(x6G08_MD_861_proj[1],x6G08_MD_861_proj[3],pch=20)
text(x6G08_MD_861_proj[1],x6G08_MD_861_proj[3],pos=4,label="861",col="red")

x6G08_MD_862_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_862.txt")
x6G08_MD_862_proj<-project.pca(x6G08_MD_862_raw,pcadata)
points(x6G08_MD_862_proj[1],x6G08_MD_862_proj[3],pch=20)
text(x6G08_MD_862_proj[1],x6G08_MD_862_proj[3],pos=4,label="862",col="red")

x6G08_MD_863_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_863.txt")
x6G08_MD_863_proj<-project.pca(x6G08_MD_863_raw,pcadata)
points(x6G08_MD_863_proj[1],x6G08_MD_863_proj[3],pch=20)
text(x6G08_MD_863_proj[1],x6G08_MD_863_proj[3],pos=4,label="863",col="red")

x6G08_MD_864_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_864.txt")
x6G08_MD_864_proj<-project.pca(x6G08_MD_864_raw,pcadata)
points(x6G08_MD_864_proj[1],x6G08_MD_864_proj[3],pch=20)
text(x6G08_MD_864_proj[1],x6G08_MD_864_proj[3],pos=4,label="864",col="red")

x6G08_MD_865_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_865.txt")
x6G08_MD_865_proj<-project.pca(x6G08_MD_865_raw,pcadata)
points(x6G08_MD_865_proj[1],x6G08_MD_865_proj[3],pch=20)
text(x6G08_MD_865_proj[1],x6G08_MD_865_proj[3],pos=4,label="865",col="red")

x6G08_MD_866_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_866.txt")
x6G08_MD_866_proj<-project.pca(x6G08_MD_866_raw,pcadata)
points(x6G08_MD_866_proj[1],x6G08_MD_866_proj[3],pch=20)
text(x6G08_MD_866_proj[1],x6G08_MD_866_proj[3],pos=4,label="866",col="red")

x6G08_MD_867_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_867.txt")
x6G08_MD_867_proj<-project.pca(x6G08_MD_867_raw,pcadata)
points(x6G08_MD_867_proj[1],x6G08_MD_867_proj[3],pch=20)
text(x6G08_MD_867_proj[1],x6G08_MD_867_proj[3],pos=4,label="867",col="red")

x6G08_MD_868_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_868.txt")
x6G08_MD_868_proj<-project.pca(x6G08_MD_868_raw,pcadata)
points(x6G08_MD_868_proj[1],x6G08_MD_868_proj[3],pch=20)
text(x6G08_MD_868_proj[1],x6G08_MD_868_proj[3],pos=4,label="868",col="red")

x6G08_MD_869_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_869.txt")
x6G08_MD_869_proj<-project.pca(x6G08_MD_869_raw,pcadata)
points(x6G08_MD_869_proj[1],x6G08_MD_869_proj[3],pch=20)
text(x6G08_MD_869_proj[1],x6G08_MD_869_proj[3],pos=4,label="869",col="red")

x6G08_MD_870_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_870.txt")
x6G08_MD_870_proj<-project.pca(x6G08_MD_870_raw,pcadata)
points(x6G08_MD_870_proj[1],x6G08_MD_870_proj[3],pch=20)
text(x6G08_MD_870_proj[1],x6G08_MD_870_proj[3],pos=4,label="870",col="red")

x6G08_MD_871_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_871.txt")
x6G08_MD_871_proj<-project.pca(x6G08_MD_871_raw,pcadata)
points(x6G08_MD_871_proj[1],x6G08_MD_871_proj[3],pch=20)
text(x6G08_MD_871_proj[1],x6G08_MD_871_proj[3],pos=4,label="871",col="red")

x6G08_MD_872_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_872.txt")
x6G08_MD_872_proj<-project.pca(x6G08_MD_872_raw,pcadata)
points(x6G08_MD_872_proj[1],x6G08_MD_872_proj[3],pch=20)
text(x6G08_MD_872_proj[1],x6G08_MD_872_proj[3],pos=4,label="872",col="red")

x6G08_MD_873_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_873.txt")
x6G08_MD_873_proj<-project.pca(x6G08_MD_873_raw,pcadata)
points(x6G08_MD_873_proj[1],x6G08_MD_873_proj[3],pch=20)
text(x6G08_MD_873_proj[1],x6G08_MD_873_proj[3],pos=4,label="873",col="red")

x6G08_MD_874_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_874.txt")
x6G08_MD_874_proj<-project.pca(x6G08_MD_874_raw,pcadata)
points(x6G08_MD_874_proj[1],x6G08_MD_874_proj[3],pch=20)
text(x6G08_MD_874_proj[1],x6G08_MD_874_proj[3],pos=4,label="874",col="red")

x6G08_MD_875_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_875.txt")
x6G08_MD_875_proj<-project.pca(x6G08_MD_875_raw,pcadata)
points(x6G08_MD_875_proj[1],x6G08_MD_875_proj[3],pch=20)
text(x6G08_MD_875_proj[1],x6G08_MD_875_proj[3],pos=4,label="875",col="red")

x6G08_MD_876_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_876.txt")
x6G08_MD_876_proj<-project.pca(x6G08_MD_876_raw,pcadata)
points(x6G08_MD_876_proj[1],x6G08_MD_876_proj[3],pch=20)
text(x6G08_MD_876_proj[1],x6G08_MD_876_proj[3],pos=4,label="876",col="red")

x6G08_MD_877_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_877.txt")
x6G08_MD_877_proj<-project.pca(x6G08_MD_877_raw,pcadata)
points(x6G08_MD_877_proj[1],x6G08_MD_877_proj[3],pch=20)
text(x6G08_MD_877_proj[1],x6G08_MD_877_proj[3],pos=4,label="877",col="red")

x6G08_MD_878_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_878.txt")
x6G08_MD_878_proj<-project.pca(x6G08_MD_878_raw,pcadata)
points(x6G08_MD_878_proj[1],x6G08_MD_878_proj[3],pch=20)
text(x6G08_MD_878_proj[1],x6G08_MD_878_proj[3],pos=4,label="878",col="red")

x6G08_MD_879_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_879.txt")
x6G08_MD_879_proj<-project.pca(x6G08_MD_879_raw,pcadata)
points(x6G08_MD_879_proj[1],x6G08_MD_879_proj[3],pch=20)
text(x6G08_MD_879_proj[1],x6G08_MD_879_proj[3],pos=4,label="879",col="red")

x6G08_MD_880_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_880.txt")
x6G08_MD_880_proj<-project.pca(x6G08_MD_880_raw,pcadata)
points(x6G08_MD_880_proj[1],x6G08_MD_880_proj[3],pch=20)
text(x6G08_MD_880_proj[1],x6G08_MD_880_proj[3],pos=4,label="880",col="red")

x6G08_MD_881_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_881.txt")
x6G08_MD_881_proj<-project.pca(x6G08_MD_881_raw,pcadata)
points(x6G08_MD_881_proj[1],x6G08_MD_881_proj[3],pch=20)
text(x6G08_MD_881_proj[1],x6G08_MD_881_proj[3],pos=4,label="881",col="red")

x6G08_MD_882_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_882.txt")
x6G08_MD_882_proj<-project.pca(x6G08_MD_882_raw,pcadata)
points(x6G08_MD_882_proj[1],x6G08_MD_882_proj[3],pch=20)
text(x6G08_MD_882_proj[1],x6G08_MD_882_proj[3],pos=4,label="882",col="red")

x6G08_MD_883_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_883.txt")
x6G08_MD_883_proj<-project.pca(x6G08_MD_883_raw,pcadata)
points(x6G08_MD_883_proj[1],x6G08_MD_883_proj[3],pch=20)
text(x6G08_MD_883_proj[1],x6G08_MD_883_proj[3],pos=4,label="883",col="red")

x6G08_MD_884_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_884.txt")
x6G08_MD_884_proj<-project.pca(x6G08_MD_884_raw,pcadata)
points(x6G08_MD_884_proj[1],x6G08_MD_884_proj[3],pch=20)
text(x6G08_MD_884_proj[1],x6G08_MD_884_proj[3],pos=4,label="884",col="red")

x6G08_MD_885_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_885.txt")
x6G08_MD_885_proj<-project.pca(x6G08_MD_885_raw,pcadata)
points(x6G08_MD_885_proj[1],x6G08_MD_885_proj[3],pch=20)
text(x6G08_MD_885_proj[1],x6G08_MD_885_proj[3],pos=4,label="885",col="red")

x6G08_MD_886_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_886.txt")
x6G08_MD_886_proj<-project.pca(x6G08_MD_886_raw,pcadata)
points(x6G08_MD_886_proj[1],x6G08_MD_886_proj[3],pch=20)
text(x6G08_MD_886_proj[1],x6G08_MD_886_proj[3],pos=4,label="886",col="red")

x6G08_MD_887_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_887.txt")
x6G08_MD_887_proj<-project.pca(x6G08_MD_887_raw,pcadata)
points(x6G08_MD_887_proj[1],x6G08_MD_887_proj[3],pch=20)
text(x6G08_MD_887_proj[1],x6G08_MD_887_proj[3],pos=4,label="887",col="red")

x6G08_MD_888_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_888.txt")
x6G08_MD_888_proj<-project.pca(x6G08_MD_888_raw,pcadata)
points(x6G08_MD_888_proj[1],x6G08_MD_888_proj[3],pch=20)
text(x6G08_MD_888_proj[1],x6G08_MD_888_proj[3],pos=4,label="888",col="red")

x6G08_MD_889_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_889.txt")
x6G08_MD_889_proj<-project.pca(x6G08_MD_889_raw,pcadata)
points(x6G08_MD_889_proj[1],x6G08_MD_889_proj[3],pch=20)
text(x6G08_MD_889_proj[1],x6G08_MD_889_proj[3],pos=4,label="889",col="red")

x6G08_MD_890_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_890.txt")
x6G08_MD_890_proj<-project.pca(x6G08_MD_890_raw,pcadata)
points(x6G08_MD_890_proj[1],x6G08_MD_890_proj[3],pch=20)
text(x6G08_MD_890_proj[1],x6G08_MD_890_proj[3],pos=4,label="890",col="red")

x6G08_MD_891_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_891.txt")
x6G08_MD_891_proj<-project.pca(x6G08_MD_891_raw,pcadata)
points(x6G08_MD_891_proj[1],x6G08_MD_891_proj[3],pch=20)
text(x6G08_MD_891_proj[1],x6G08_MD_891_proj[3],pos=4,label="891",col="red")

x6G08_MD_892_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_892.txt")
x6G08_MD_892_proj<-project.pca(x6G08_MD_892_raw,pcadata)
points(x6G08_MD_892_proj[1],x6G08_MD_892_proj[3],pch=20)
text(x6G08_MD_892_proj[1],x6G08_MD_892_proj[3],pos=4,label="892",col="red")

x6G08_MD_893_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_893.txt")
x6G08_MD_893_proj<-project.pca(x6G08_MD_893_raw,pcadata)
points(x6G08_MD_893_proj[1],x6G08_MD_893_proj[3],pch=20)
text(x6G08_MD_893_proj[1],x6G08_MD_893_proj[3],pos=4,label="893",col="red")

x6G08_MD_894_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_894.txt")
x6G08_MD_894_proj<-project.pca(x6G08_MD_894_raw,pcadata)
points(x6G08_MD_894_proj[1],x6G08_MD_894_proj[3],pch=20)
text(x6G08_MD_894_proj[1],x6G08_MD_894_proj[3],pos=4,label="894",col="red")

x6G08_MD_895_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_895.txt")
x6G08_MD_895_proj<-project.pca(x6G08_MD_895_raw,pcadata)
points(x6G08_MD_895_proj[1],x6G08_MD_895_proj[3],pch=20)
text(x6G08_MD_895_proj[1],x6G08_MD_895_proj[3],pos=4,label="895",col="green")

x6G08_MD_896_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_896.txt")
x6G08_MD_896_proj<-project.pca(x6G08_MD_896_raw,pcadata)
points(x6G08_MD_896_proj[1],x6G08_MD_896_proj[3],pch=20)
text(x6G08_MD_896_proj[1],x6G08_MD_896_proj[3],pos=4,label="896",col="red")

x6G08_MD_897_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_897.txt")
x6G08_MD_897_proj<-project.pca(x6G08_MD_897_raw,pcadata)
points(x6G08_MD_897_proj[1],x6G08_MD_897_proj[3],pch=20)
text(x6G08_MD_897_proj[1],x6G08_MD_897_proj[3],pos=4,label="897",col="green")

x6G08_MD_898_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_898.txt")
x6G08_MD_898_proj<-project.pca(x6G08_MD_898_raw,pcadata)
points(x6G08_MD_898_proj[1],x6G08_MD_898_proj[3],pch=20)
text(x6G08_MD_898_proj[1],x6G08_MD_898_proj[3],pos=4,label="898",col="green")

x6G08_MD_899_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_899.txt")
x6G08_MD_899_proj<-project.pca(x6G08_MD_899_raw,pcadata)
points(x6G08_MD_899_proj[1],x6G08_MD_899_proj[3],pch=20)
text(x6G08_MD_899_proj[1],x6G08_MD_899_proj[3],pos=4,label="899",col="red")

x6G08_MD_900_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_900.txt")
x6G08_MD_900_proj<-project.pca(x6G08_MD_900_raw,pcadata)
points(x6G08_MD_900_proj[1],x6G08_MD_900_proj[3],pch=20)
text(x6G08_MD_900_proj[1],x6G08_MD_900_proj[3],pos=4,label="900",col="red")

x6G08_MD_901_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_901.txt")
x6G08_MD_901_proj<-project.pca(x6G08_MD_901_raw,pcadata)
points(x6G08_MD_901_proj[1],x6G08_MD_901_proj[3],pch=20)
text(x6G08_MD_901_proj[1],x6G08_MD_901_proj[3],pos=4,label="901",col="red")

x6G08_MD_902_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_902.txt")
x6G08_MD_902_proj<-project.pca(x6G08_MD_902_raw,pcadata)
points(x6G08_MD_902_proj[1],x6G08_MD_902_proj[3],pch=20)
text(x6G08_MD_902_proj[1],x6G08_MD_902_proj[3],pos=4,label="902",col="red")

x6G08_MD_903_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_903.txt")
x6G08_MD_903_proj<-project.pca(x6G08_MD_903_raw,pcadata)
points(x6G08_MD_903_proj[1],x6G08_MD_903_proj[3],pch=20)
text(x6G08_MD_903_proj[1],x6G08_MD_903_proj[3],pos=4,label="903",col="red")

x6G08_MD_904_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_904.txt")
x6G08_MD_904_proj<-project.pca(x6G08_MD_904_raw,pcadata)
points(x6G08_MD_904_proj[1],x6G08_MD_904_proj[3],pch=20)
text(x6G08_MD_904_proj[1],x6G08_MD_904_proj[3],pos=4,label="904",col="red")

x6G08_MD_905_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_905.txt")
x6G08_MD_905_proj<-project.pca(x6G08_MD_905_raw,pcadata)
points(x6G08_MD_905_proj[1],x6G08_MD_905_proj[3],pch=20)
text(x6G08_MD_905_proj[1],x6G08_MD_905_proj[3],pos=4,label="905",col="red")

x6G08_MD_906_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_906.txt")
x6G08_MD_906_proj<-project.pca(x6G08_MD_906_raw,pcadata)
points(x6G08_MD_906_proj[1],x6G08_MD_906_proj[3],pch=20)
text(x6G08_MD_906_proj[1],x6G08_MD_906_proj[3],pos=4,label="906",col="red")

x6G08_MD_907_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_907.txt")
x6G08_MD_907_proj<-project.pca(x6G08_MD_907_raw,pcadata)
points(x6G08_MD_907_proj[1],x6G08_MD_907_proj[3],pch=20)
text(x6G08_MD_907_proj[1],x6G08_MD_907_proj[3],pos=4,label="907",col="red")

x6G08_MD_908_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_908.txt")
x6G08_MD_908_proj<-project.pca(x6G08_MD_908_raw,pcadata)
points(x6G08_MD_908_proj[1],x6G08_MD_908_proj[3],pch=20)
text(x6G08_MD_908_proj[1],x6G08_MD_908_proj[3],pos=4,label="908",col="red")

x6G08_MD_909_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_909.txt")
x6G08_MD_909_proj<-project.pca(x6G08_MD_909_raw,pcadata)
points(x6G08_MD_909_proj[1],x6G08_MD_909_proj[3],pch=20)
text(x6G08_MD_909_proj[1],x6G08_MD_909_proj[3],pos=4,label="909",col="red")

x6G08_MD_910_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_910.txt")
x6G08_MD_910_proj<-project.pca(x6G08_MD_910_raw,pcadata)
points(x6G08_MD_910_proj[1],x6G08_MD_910_proj[3],pch=20)
text(x6G08_MD_910_proj[1],x6G08_MD_910_proj[3],pos=4,label="910",col="red")

x6G08_MD_911_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_911.txt")
x6G08_MD_911_proj<-project.pca(x6G08_MD_911_raw,pcadata)
points(x6G08_MD_911_proj[1],x6G08_MD_911_proj[3],pch=20)
text(x6G08_MD_911_proj[1],x6G08_MD_911_proj[3],pos=4,label="911",col="red")

x6G08_MD_912_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_912.txt")
x6G08_MD_912_proj<-project.pca(x6G08_MD_912_raw,pcadata)
points(x6G08_MD_912_proj[1],x6G08_MD_912_proj[3],pch=20)
text(x6G08_MD_912_proj[1],x6G08_MD_912_proj[3],pos=4,label="912",col="red")

x6G08_MD_913_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_913.txt")
x6G08_MD_913_proj<-project.pca(x6G08_MD_913_raw,pcadata)
points(x6G08_MD_913_proj[1],x6G08_MD_913_proj[3],pch=20)
text(x6G08_MD_913_proj[1],x6G08_MD_913_proj[3],pos=4,label="913",col="green")

x6G08_MD_914_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_914.txt")
x6G08_MD_914_proj<-project.pca(x6G08_MD_914_raw,pcadata)
points(x6G08_MD_914_proj[1],x6G08_MD_914_proj[3],pch=20)
text(x6G08_MD_914_proj[1],x6G08_MD_914_proj[3],pos=4,label="914",col="red")

x6G08_MD_915_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_915.txt")
x6G08_MD_915_proj<-project.pca(x6G08_MD_915_raw,pcadata)
points(x6G08_MD_915_proj[1],x6G08_MD_915_proj[3],pch=20)
text(x6G08_MD_915_proj[1],x6G08_MD_915_proj[3],pos=4,label="915",col="red")

x6G08_MD_916_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_916.txt")
x6G08_MD_916_proj<-project.pca(x6G08_MD_916_raw,pcadata)
points(x6G08_MD_916_proj[1],x6G08_MD_916_proj[3],pch=20)
text(x6G08_MD_916_proj[1],x6G08_MD_916_proj[3],pos=4,label="916",col="red")

x6G08_MD_917_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_917.txt")
x6G08_MD_917_proj<-project.pca(x6G08_MD_917_raw,pcadata)
points(x6G08_MD_917_proj[1],x6G08_MD_917_proj[3],pch=20)
text(x6G08_MD_917_proj[1],x6G08_MD_917_proj[3],pos=4,label="917",col="red")

x6G08_MD_918_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_918.txt")
x6G08_MD_918_proj<-project.pca(x6G08_MD_918_raw,pcadata)
points(x6G08_MD_918_proj[1],x6G08_MD_918_proj[3],pch=20)
text(x6G08_MD_918_proj[1],x6G08_MD_918_proj[3],pos=4,label="918",col="red")

x6G08_MD_919_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_919.txt")
x6G08_MD_919_proj<-project.pca(x6G08_MD_919_raw,pcadata)
points(x6G08_MD_919_proj[1],x6G08_MD_919_proj[3],pch=20)
text(x6G08_MD_919_proj[1],x6G08_MD_919_proj[3],pos=4,label="919",col="red")

x6G08_MD_920_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_920.txt")
x6G08_MD_920_proj<-project.pca(x6G08_MD_920_raw,pcadata)
points(x6G08_MD_920_proj[1],x6G08_MD_920_proj[3],pch=20)
text(x6G08_MD_920_proj[1],x6G08_MD_920_proj[3],pos=4,label="920",col="red")

x6G08_MD_921_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_921.txt")
x6G08_MD_921_proj<-project.pca(x6G08_MD_921_raw,pcadata)
points(x6G08_MD_921_proj[1],x6G08_MD_921_proj[3],pch=20)
text(x6G08_MD_921_proj[1],x6G08_MD_921_proj[3],pos=4,label="921",col="red")

x6G08_MD_922_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_922.txt")
x6G08_MD_922_proj<-project.pca(x6G08_MD_922_raw,pcadata)
points(x6G08_MD_922_proj[1],x6G08_MD_922_proj[3],pch=20)
text(x6G08_MD_922_proj[1],x6G08_MD_922_proj[3],pos=4,label="922",col="red")

x6G08_MD_923_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_923.txt")
x6G08_MD_923_proj<-project.pca(x6G08_MD_923_raw,pcadata)
points(x6G08_MD_923_proj[1],x6G08_MD_923_proj[3],pch=20)
text(x6G08_MD_923_proj[1],x6G08_MD_923_proj[3],pos=4,label="923",col="red")

x6G08_MD_924_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_924.txt")
x6G08_MD_924_proj<-project.pca(x6G08_MD_924_raw,pcadata)
points(x6G08_MD_924_proj[1],x6G08_MD_924_proj[3],pch=20)
text(x6G08_MD_924_proj[1],x6G08_MD_924_proj[3],pos=4,label="924",col="red")

x6G08_MD_925_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_925.txt")
x6G08_MD_925_proj<-project.pca(x6G08_MD_925_raw,pcadata)
points(x6G08_MD_925_proj[1],x6G08_MD_925_proj[3],pch=20)
text(x6G08_MD_925_proj[1],x6G08_MD_925_proj[3],pos=4,label="925",col="red")

x6G08_MD_926_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_926.txt")
x6G08_MD_926_proj<-project.pca(x6G08_MD_926_raw,pcadata)
points(x6G08_MD_926_proj[1],x6G08_MD_926_proj[3],pch=20)
text(x6G08_MD_926_proj[1],x6G08_MD_926_proj[3],pos=4,label="926",col="red")

x6G08_MD_927_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_927.txt")
x6G08_MD_927_proj<-project.pca(x6G08_MD_927_raw,pcadata)
points(x6G08_MD_927_proj[1],x6G08_MD_927_proj[3],pch=20)
text(x6G08_MD_927_proj[1],x6G08_MD_927_proj[3],pos=4,label="927",col="red")

x6G08_MD_928_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_928.txt")
x6G08_MD_928_proj<-project.pca(x6G08_MD_928_raw,pcadata)
points(x6G08_MD_928_proj[1],x6G08_MD_928_proj[3],pch=20)
text(x6G08_MD_928_proj[1],x6G08_MD_928_proj[3],pos=4,label="928",col="red")

x6G08_MD_929_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_929.txt")
x6G08_MD_929_proj<-project.pca(x6G08_MD_929_raw,pcadata)
points(x6G08_MD_929_proj[1],x6G08_MD_929_proj[3],pch=20)
text(x6G08_MD_929_proj[1],x6G08_MD_929_proj[3],pos=4,label="929",col="red")

x6G08_MD_930_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_930.txt")
x6G08_MD_930_proj<-project.pca(x6G08_MD_930_raw,pcadata)
points(x6G08_MD_930_proj[1],x6G08_MD_930_proj[3],pch=20)
text(x6G08_MD_930_proj[1],x6G08_MD_930_proj[3],pos=4,label="930",col="red")

x6G08_MD_931_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_931.txt")
x6G08_MD_931_proj<-project.pca(x6G08_MD_931_raw,pcadata)
points(x6G08_MD_931_proj[1],x6G08_MD_931_proj[3],pch=20)
text(x6G08_MD_931_proj[1],x6G08_MD_931_proj[3],pos=4,label="931",col="red")

x6G08_MD_932_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_932.txt")
x6G08_MD_932_proj<-project.pca(x6G08_MD_932_raw,pcadata)
points(x6G08_MD_932_proj[1],x6G08_MD_932_proj[3],pch=20)
text(x6G08_MD_932_proj[1],x6G08_MD_932_proj[3],pos=4,label="932",col="red")

x6G08_MD_933_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_933.txt")
x6G08_MD_933_proj<-project.pca(x6G08_MD_933_raw,pcadata)
points(x6G08_MD_933_proj[1],x6G08_MD_933_proj[3],pch=20)
text(x6G08_MD_933_proj[1],x6G08_MD_933_proj[3],pos=4,label="933",col="red")

x6G08_MD_934_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_934.txt")
x6G08_MD_934_proj<-project.pca(x6G08_MD_934_raw,pcadata)
points(x6G08_MD_934_proj[1],x6G08_MD_934_proj[3],pch=20)
text(x6G08_MD_934_proj[1],x6G08_MD_934_proj[3],pos=4,label="934",col="red")

x6G08_MD_935_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_935.txt")
x6G08_MD_935_proj<-project.pca(x6G08_MD_935_raw,pcadata)
points(x6G08_MD_935_proj[1],x6G08_MD_935_proj[3],pch=20)
text(x6G08_MD_935_proj[1],x6G08_MD_935_proj[3],pos=4,label="935",col="red")

x6G08_MD_936_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_936.txt")
x6G08_MD_936_proj<-project.pca(x6G08_MD_936_raw,pcadata)
points(x6G08_MD_936_proj[1],x6G08_MD_936_proj[3],pch=20)
text(x6G08_MD_936_proj[1],x6G08_MD_936_proj[3],pos=4,label="936",col="red")

x6G08_MD_937_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_937.txt")
x6G08_MD_937_proj<-project.pca(x6G08_MD_937_raw,pcadata)
points(x6G08_MD_937_proj[1],x6G08_MD_937_proj[3],pch=20)
text(x6G08_MD_937_proj[1],x6G08_MD_937_proj[3],pos=4,label="937",col="green")

x6G08_MD_938_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_938.txt")
x6G08_MD_938_proj<-project.pca(x6G08_MD_938_raw,pcadata)
points(x6G08_MD_938_proj[1],x6G08_MD_938_proj[3],pch=20)
text(x6G08_MD_938_proj[1],x6G08_MD_938_proj[3],pos=4,label="938",col="green")

x6G08_MD_939_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_939.txt")
x6G08_MD_939_proj<-project.pca(x6G08_MD_939_raw,pcadata)
points(x6G08_MD_939_proj[1],x6G08_MD_939_proj[3],pch=20)
text(x6G08_MD_939_proj[1],x6G08_MD_939_proj[3],pos=4,label="939",col="red")

x6G08_MD_940_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_940.txt")
x6G08_MD_940_proj<-project.pca(x6G08_MD_940_raw,pcadata)
points(x6G08_MD_940_proj[1],x6G08_MD_940_proj[3],pch=20)
text(x6G08_MD_940_proj[1],x6G08_MD_940_proj[3],pos=4,label="940",col="red")

x6G08_MD_941_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_941.txt")
x6G08_MD_941_proj<-project.pca(x6G08_MD_941_raw,pcadata)
points(x6G08_MD_941_proj[1],x6G08_MD_941_proj[3],pch=20)
text(x6G08_MD_941_proj[1],x6G08_MD_941_proj[3],pos=4,label="941",col="red")

x6G08_MD_942_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_942.txt")
x6G08_MD_942_proj<-project.pca(x6G08_MD_942_raw,pcadata)
points(x6G08_MD_942_proj[1],x6G08_MD_942_proj[3],pch=20)
text(x6G08_MD_942_proj[1],x6G08_MD_942_proj[3],pos=4,label="942",col="red")

x6G08_MD_943_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_943.txt")
x6G08_MD_943_proj<-project.pca(x6G08_MD_943_raw,pcadata)
points(x6G08_MD_943_proj[1],x6G08_MD_943_proj[3],pch=20)
text(x6G08_MD_943_proj[1],x6G08_MD_943_proj[3],pos=4,label="943",col="red")

x6G08_MD_944_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_944.txt")
x6G08_MD_944_proj<-project.pca(x6G08_MD_944_raw,pcadata)
points(x6G08_MD_944_proj[1],x6G08_MD_944_proj[3],pch=20)
text(x6G08_MD_944_proj[1],x6G08_MD_944_proj[3],pos=4,label="944",col="red")

x6G08_MD_945_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_945.txt")
x6G08_MD_945_proj<-project.pca(x6G08_MD_945_raw,pcadata)
points(x6G08_MD_945_proj[1],x6G08_MD_945_proj[3],pch=20)
text(x6G08_MD_945_proj[1],x6G08_MD_945_proj[3],pos=4,label="945",col="red")

x6G08_MD_946_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_946.txt")
x6G08_MD_946_proj<-project.pca(x6G08_MD_946_raw,pcadata)
points(x6G08_MD_946_proj[1],x6G08_MD_946_proj[3],pch=20)
text(x6G08_MD_946_proj[1],x6G08_MD_946_proj[3],pos=4,label="946",col="red")

x6G08_MD_947_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_947.txt")
x6G08_MD_947_proj<-project.pca(x6G08_MD_947_raw,pcadata)
points(x6G08_MD_947_proj[1],x6G08_MD_947_proj[3],pch=20)
text(x6G08_MD_947_proj[1],x6G08_MD_947_proj[3],pos=4,label="947",col="red")

x6G08_MD_948_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_948.txt")
x6G08_MD_948_proj<-project.pca(x6G08_MD_948_raw,pcadata)
points(x6G08_MD_948_proj[1],x6G08_MD_948_proj[3],pch=20)
text(x6G08_MD_948_proj[1],x6G08_MD_948_proj[3],pos=4,label="948",col="red")

x6G08_MD_949_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_949.txt")
x6G08_MD_949_proj<-project.pca(x6G08_MD_949_raw,pcadata)
points(x6G08_MD_949_proj[1],x6G08_MD_949_proj[3],pch=20)
text(x6G08_MD_949_proj[1],x6G08_MD_949_proj[3],pos=4,label="949",col="red")

x6G08_MD_950_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_950.txt")
x6G08_MD_950_proj<-project.pca(x6G08_MD_950_raw,pcadata)
points(x6G08_MD_950_proj[1],x6G08_MD_950_proj[3],pch=20)
text(x6G08_MD_950_proj[1],x6G08_MD_950_proj[3],pos=4,label="950",col="red")

x6G08_MD_951_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_951.txt")
x6G08_MD_951_proj<-project.pca(x6G08_MD_951_raw,pcadata)
points(x6G08_MD_951_proj[1],x6G08_MD_951_proj[3],pch=20)
text(x6G08_MD_951_proj[1],x6G08_MD_951_proj[3],pos=4,label="951",col="red")

x6G08_MD_952_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_952.txt")
x6G08_MD_952_proj<-project.pca(x6G08_MD_952_raw,pcadata)
points(x6G08_MD_952_proj[1],x6G08_MD_952_proj[3],pch=20)
text(x6G08_MD_952_proj[1],x6G08_MD_952_proj[3],pos=4,label="952",col="red")

x6G08_MD_953_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_953.txt")
x6G08_MD_953_proj<-project.pca(x6G08_MD_953_raw,pcadata)
points(x6G08_MD_953_proj[1],x6G08_MD_953_proj[3],pch=20)
text(x6G08_MD_953_proj[1],x6G08_MD_953_proj[3],pos=4,label="953",col="red")

x6G08_MD_954_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_954.txt")
x6G08_MD_954_proj<-project.pca(x6G08_MD_954_raw,pcadata)
points(x6G08_MD_954_proj[1],x6G08_MD_954_proj[3],pch=20)
text(x6G08_MD_954_proj[1],x6G08_MD_954_proj[3],pos=4,label="954",col="red")

x6G08_MD_955_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_955.txt")
x6G08_MD_955_proj<-project.pca(x6G08_MD_955_raw,pcadata)
points(x6G08_MD_955_proj[1],x6G08_MD_955_proj[3],pch=20)
text(x6G08_MD_955_proj[1],x6G08_MD_955_proj[3],pos=4,label="955",col="red")

x6G08_MD_956_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_956.txt")
x6G08_MD_956_proj<-project.pca(x6G08_MD_956_raw,pcadata)
points(x6G08_MD_956_proj[1],x6G08_MD_956_proj[3],pch=20)
text(x6G08_MD_956_proj[1],x6G08_MD_956_proj[3],pos=4,label="956",col="red")

x6G08_MD_957_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_957.txt")
x6G08_MD_957_proj<-project.pca(x6G08_MD_957_raw,pcadata)
points(x6G08_MD_957_proj[1],x6G08_MD_957_proj[3],pch=20)
text(x6G08_MD_957_proj[1],x6G08_MD_957_proj[3],pos=4,label="957",col="red")

x6G08_MD_958_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_958.txt")
x6G08_MD_958_proj<-project.pca(x6G08_MD_958_raw,pcadata)
points(x6G08_MD_958_proj[1],x6G08_MD_958_proj[3],pch=20)
text(x6G08_MD_958_proj[1],x6G08_MD_958_proj[3],pos=4,label="958",col="red")

x6G08_MD_959_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_959.txt")
x6G08_MD_959_proj<-project.pca(x6G08_MD_959_raw,pcadata)
points(x6G08_MD_959_proj[1],x6G08_MD_959_proj[3],pch=20)
text(x6G08_MD_959_proj[1],x6G08_MD_959_proj[3],pos=4,label="959",col="red")

x6G08_MD_960_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_960.txt")
x6G08_MD_960_proj<-project.pca(x6G08_MD_960_raw,pcadata)
points(x6G08_MD_960_proj[1],x6G08_MD_960_proj[3],pch=20)
text(x6G08_MD_960_proj[1],x6G08_MD_960_proj[3],pos=4,label="960",col="red")

x6G08_MD_961_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_961.txt")
x6G08_MD_961_proj<-project.pca(x6G08_MD_961_raw,pcadata)
points(x6G08_MD_961_proj[1],x6G08_MD_961_proj[3],pch=20)
text(x6G08_MD_961_proj[1],x6G08_MD_961_proj[3],pos=4,label="961",col="red")

x6G08_MD_962_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_962.txt")
x6G08_MD_962_proj<-project.pca(x6G08_MD_962_raw,pcadata)
points(x6G08_MD_962_proj[1],x6G08_MD_962_proj[3],pch=20)
text(x6G08_MD_962_proj[1],x6G08_MD_962_proj[3],pos=4,label="962",col="red")

x6G08_MD_963_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_963.txt")
x6G08_MD_963_proj<-project.pca(x6G08_MD_963_raw,pcadata)
points(x6G08_MD_963_proj[1],x6G08_MD_963_proj[3],pch=20)
text(x6G08_MD_963_proj[1],x6G08_MD_963_proj[3],pos=4,label="963",col="red")

x6G08_MD_964_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_964.txt")
x6G08_MD_964_proj<-project.pca(x6G08_MD_964_raw,pcadata)
points(x6G08_MD_964_proj[1],x6G08_MD_964_proj[3],pch=20)
text(x6G08_MD_964_proj[1],x6G08_MD_964_proj[3],pos=4,label="964",col="red")

x6G08_MD_965_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_965.txt")
x6G08_MD_965_proj<-project.pca(x6G08_MD_965_raw,pcadata)
points(x6G08_MD_965_proj[1],x6G08_MD_965_proj[3],pch=20)
text(x6G08_MD_965_proj[1],x6G08_MD_965_proj[3],pos=4,label="965",col="red")

x6G08_MD_966_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_966.txt")
x6G08_MD_966_proj<-project.pca(x6G08_MD_966_raw,pcadata)
points(x6G08_MD_966_proj[1],x6G08_MD_966_proj[3],pch=20)
text(x6G08_MD_966_proj[1],x6G08_MD_966_proj[3],pos=4,label="966",col="red")

x6G08_MD_967_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_967.txt")
x6G08_MD_967_proj<-project.pca(x6G08_MD_967_raw,pcadata)
points(x6G08_MD_967_proj[1],x6G08_MD_967_proj[3],pch=20)
text(x6G08_MD_967_proj[1],x6G08_MD_967_proj[3],pos=4,label="967",col="red")

x6G08_MD_968_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_968.txt")
x6G08_MD_968_proj<-project.pca(x6G08_MD_968_raw,pcadata)
points(x6G08_MD_968_proj[1],x6G08_MD_968_proj[3],pch=20)
text(x6G08_MD_968_proj[1],x6G08_MD_968_proj[3],pos=4,label="968",col="red")

x6G08_MD_969_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_969.txt")
x6G08_MD_969_proj<-project.pca(x6G08_MD_969_raw,pcadata)
points(x6G08_MD_969_proj[1],x6G08_MD_969_proj[3],pch=20)
text(x6G08_MD_969_proj[1],x6G08_MD_969_proj[3],pos=4,label="969",col="red")

x6G08_MD_970_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_970.txt")
x6G08_MD_970_proj<-project.pca(x6G08_MD_970_raw,pcadata)
points(x6G08_MD_970_proj[1],x6G08_MD_970_proj[3],pch=20)
text(x6G08_MD_970_proj[1],x6G08_MD_970_proj[3],pos=4,label="970",col="red")

x6G08_MD_971_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_971.txt")
x6G08_MD_971_proj<-project.pca(x6G08_MD_971_raw,pcadata)
points(x6G08_MD_971_proj[1],x6G08_MD_971_proj[3],pch=20)
text(x6G08_MD_971_proj[1],x6G08_MD_971_proj[3],pos=4,label="971",col="red")

x6G08_MD_972_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_972.txt")
x6G08_MD_972_proj<-project.pca(x6G08_MD_972_raw,pcadata)
points(x6G08_MD_972_proj[1],x6G08_MD_972_proj[3],pch=20)
text(x6G08_MD_972_proj[1],x6G08_MD_972_proj[3],pos=4,label="972",col="red")

x6G08_MD_973_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_973.txt")
x6G08_MD_973_proj<-project.pca(x6G08_MD_973_raw,pcadata)
points(x6G08_MD_973_proj[1],x6G08_MD_973_proj[3],pch=20)
text(x6G08_MD_973_proj[1],x6G08_MD_973_proj[3],pos=4,label="973",col="red")

x6G08_MD_974_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_974.txt")
x6G08_MD_974_proj<-project.pca(x6G08_MD_974_raw,pcadata)
points(x6G08_MD_974_proj[1],x6G08_MD_974_proj[3],pch=20)
text(x6G08_MD_974_proj[1],x6G08_MD_974_proj[3],pos=4,label="974",col="red")

x6G08_MD_975_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_975.txt")
x6G08_MD_975_proj<-project.pca(x6G08_MD_975_raw,pcadata)
points(x6G08_MD_975_proj[1],x6G08_MD_975_proj[3],pch=20)
text(x6G08_MD_975_proj[1],x6G08_MD_975_proj[3],pos=4,label="975",col="red")

x6G08_MD_976_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_976.txt")
x6G08_MD_976_proj<-project.pca(x6G08_MD_976_raw,pcadata)
points(x6G08_MD_976_proj[1],x6G08_MD_976_proj[3],pch=20)
text(x6G08_MD_976_proj[1],x6G08_MD_976_proj[3],pos=4,label="976",col="red")

x6G08_MD_977_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_977.txt")
x6G08_MD_977_proj<-project.pca(x6G08_MD_977_raw,pcadata)
points(x6G08_MD_977_proj[1],x6G08_MD_977_proj[3],pch=20)
text(x6G08_MD_977_proj[1],x6G08_MD_977_proj[3],pos=4,label="977",col="red")

x6G08_MD_978_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_978.txt")
x6G08_MD_978_proj<-project.pca(x6G08_MD_978_raw,pcadata)
points(x6G08_MD_978_proj[1],x6G08_MD_978_proj[3],pch=20)
text(x6G08_MD_978_proj[1],x6G08_MD_978_proj[3],pos=4,label="978",col="red")

x6G08_MD_979_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_979.txt")
x6G08_MD_979_proj<-project.pca(x6G08_MD_979_raw,pcadata)
points(x6G08_MD_979_proj[1],x6G08_MD_979_proj[3],pch=20)
text(x6G08_MD_979_proj[1],x6G08_MD_979_proj[3],pos=4,label="979",col="red")

x6G08_MD_980_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_980.txt")
x6G08_MD_980_proj<-project.pca(x6G08_MD_980_raw,pcadata)
points(x6G08_MD_980_proj[1],x6G08_MD_980_proj[3],pch=20)
text(x6G08_MD_980_proj[1],x6G08_MD_980_proj[3],pos=4,label="980",col="green")

x6G08_MD_981_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_981.txt")
x6G08_MD_981_proj<-project.pca(x6G08_MD_981_raw,pcadata)
points(x6G08_MD_981_proj[1],x6G08_MD_981_proj[3],pch=20)
text(x6G08_MD_981_proj[1],x6G08_MD_981_proj[3],pos=4,label="981",col="red")

x6G08_MD_982_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_982.txt")
x6G08_MD_982_proj<-project.pca(x6G08_MD_982_raw,pcadata)
points(x6G08_MD_982_proj[1],x6G08_MD_982_proj[3],pch=20)
text(x6G08_MD_982_proj[1],x6G08_MD_982_proj[3],pos=4,label="982",col="red")

x6G08_MD_983_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_983.txt")
x6G08_MD_983_proj<-project.pca(x6G08_MD_983_raw,pcadata)
points(x6G08_MD_983_proj[1],x6G08_MD_983_proj[3],pch=20)
text(x6G08_MD_983_proj[1],x6G08_MD_983_proj[3],pos=4,label="983",col="green")

x6G08_MD_984_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_984.txt")
x6G08_MD_984_proj<-project.pca(x6G08_MD_984_raw,pcadata)
points(x6G08_MD_984_proj[1],x6G08_MD_984_proj[3],pch=20)
text(x6G08_MD_984_proj[1],x6G08_MD_984_proj[3],pos=4,label="984",col="green")

x6G08_MD_985_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_985.txt")
x6G08_MD_985_proj<-project.pca(x6G08_MD_985_raw,pcadata)
points(x6G08_MD_985_proj[1],x6G08_MD_985_proj[3],pch=20)
text(x6G08_MD_985_proj[1],x6G08_MD_985_proj[3],pos=4,label="985",col="red")

x6G08_MD_986_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_986.txt")
x6G08_MD_986_proj<-project.pca(x6G08_MD_986_raw,pcadata)
points(x6G08_MD_986_proj[1],x6G08_MD_986_proj[3],pch=20)
text(x6G08_MD_986_proj[1],x6G08_MD_986_proj[3],pos=4,label="986",col="red")

x6G08_MD_987_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_987.txt")
x6G08_MD_987_proj<-project.pca(x6G08_MD_987_raw,pcadata)
points(x6G08_MD_987_proj[1],x6G08_MD_987_proj[3],pch=20)
text(x6G08_MD_987_proj[1],x6G08_MD_987_proj[3],pos=4,label="987",col="green")

x6G08_MD_988_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_988.txt")
x6G08_MD_988_proj<-project.pca(x6G08_MD_988_raw,pcadata)
points(x6G08_MD_988_proj[1],x6G08_MD_988_proj[3],pch=20)
text(x6G08_MD_988_proj[1],x6G08_MD_988_proj[3],pos=4,label="988",col="green")

x6G08_MD_989_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_989.txt")
x6G08_MD_989_proj<-project.pca(x6G08_MD_989_raw,pcadata)
points(x6G08_MD_989_proj[1],x6G08_MD_989_proj[3],pch=20)
text(x6G08_MD_989_proj[1],x6G08_MD_989_proj[3],pos=4,label="989",col="green")

x6G08_MD_990_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_990.txt")
x6G08_MD_990_proj<-project.pca(x6G08_MD_990_raw,pcadata)
points(x6G08_MD_990_proj[1],x6G08_MD_990_proj[3],pch=20)
text(x6G08_MD_990_proj[1],x6G08_MD_990_proj[3],pos=4,label="990",col="red")

x6G08_MD_991_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_991.txt")
x6G08_MD_991_proj<-project.pca(x6G08_MD_991_raw,pcadata)
points(x6G08_MD_991_proj[1],x6G08_MD_991_proj[3],pch=20)
text(x6G08_MD_991_proj[1],x6G08_MD_991_proj[3],pos=4,label="991",col="red")

x6G08_MD_992_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_992.txt")
x6G08_MD_992_proj<-project.pca(x6G08_MD_992_raw,pcadata)
points(x6G08_MD_992_proj[1],x6G08_MD_992_proj[3],pch=20)
text(x6G08_MD_992_proj[1],x6G08_MD_992_proj[3],pos=4,label="992",col="green")

x6G08_MD_993_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_993.txt")
x6G08_MD_993_proj<-project.pca(x6G08_MD_993_raw,pcadata)
points(x6G08_MD_993_proj[1],x6G08_MD_993_proj[3],pch=20)
text(x6G08_MD_993_proj[1],x6G08_MD_993_proj[3],pos=4,label="993",col="green")

x6G08_MD_994_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_994.txt")
x6G08_MD_994_proj<-project.pca(x6G08_MD_994_raw,pcadata)
points(x6G08_MD_994_proj[1],x6G08_MD_994_proj[3],pch=20)
text(x6G08_MD_994_proj[1],x6G08_MD_994_proj[3],pos=4,label="994",col="red")

x6G08_MD_995_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_995.txt")
x6G08_MD_995_proj<-project.pca(x6G08_MD_995_raw,pcadata)
points(x6G08_MD_995_proj[1],x6G08_MD_995_proj[3],pch=20)
text(x6G08_MD_995_proj[1],x6G08_MD_995_proj[3],pos=4,label="995",col="green")

x6G08_MD_996_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_996.txt")
x6G08_MD_996_proj<-project.pca(x6G08_MD_996_raw,pcadata)
points(x6G08_MD_996_proj[1],x6G08_MD_996_proj[3],pch=20)
text(x6G08_MD_996_proj[1],x6G08_MD_996_proj[3],pos=4,label="996",col="green")

x6G08_MD_997_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_997.txt")
x6G08_MD_997_proj<-project.pca(x6G08_MD_997_raw,pcadata)
points(x6G08_MD_997_proj[1],x6G08_MD_997_proj[3],pch=20)
text(x6G08_MD_997_proj[1],x6G08_MD_997_proj[3],pos=4,label="997",col="green")

x6G08_MD_998_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_998.txt")
x6G08_MD_998_proj<-project.pca(x6G08_MD_998_raw,pcadata)
points(x6G08_MD_998_proj[1],x6G08_MD_998_proj[3],pch=20)
text(x6G08_MD_998_proj[1],x6G08_MD_998_proj[3],pos=4,label="998",col="blue")

x6G08_MD_999_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_999.txt")
x6G08_MD_999_proj<-project.pca(x6G08_MD_999_raw,pcadata)
points(x6G08_MD_999_proj[1],x6G08_MD_999_proj[3],pch=20)
text(x6G08_MD_999_proj[1],x6G08_MD_999_proj[3],pos=4,label="999",col="green")

x6G08_MD_1000_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_1000.txt")
x6G08_MD_1000_proj<-project.pca(x6G08_MD_1000_raw,pcadata)
points(x6G08_MD_1000_proj[1],x6G08_MD_1000_proj[3],pch=20)
text(x6G08_MD_1000_proj[1],x6G08_MD_1000_proj[3],pos=4,label="1000",col="red")


plot(PC3,PC2,xlim=c(-40,50),ylim=c(-40,50))

xCA_6G08_fab_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/CA_6G08_fab.txt")
xCA_6G08_fab_proj<-project.pca(xCA_6G08_fab_raw,pcadata)
points(xCA_6G08_fab_proj[3],xCA_6G08_fab_proj[2],pch=20)
text(xCA_6G08_fab_proj[3],xCA_6G08_fab_proj[2],pos=4,label="6G08",col="purple")

x6G08_MD_1_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_1.txt")
x6G08_MD_1_proj<-project.pca(x6G08_MD_1_raw,pcadata)
points(x6G08_MD_1_proj[3],x6G08_MD_1_proj[2],pch=20)
text(x6G08_MD_1_proj[3],x6G08_MD_1_proj[2],pos=4,label="1",col="blue")

x6G08_MD_2_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_2.txt")
x6G08_MD_2_proj<-project.pca(x6G08_MD_2_raw,pcadata)
points(x6G08_MD_2_proj[3],x6G08_MD_2_proj[2],pch=20)
text(x6G08_MD_2_proj[3],x6G08_MD_2_proj[2],pos=4,label="2",col="blue")

x6G08_MD_3_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_3.txt")
x6G08_MD_3_proj<-project.pca(x6G08_MD_3_raw,pcadata)
points(x6G08_MD_3_proj[3],x6G08_MD_3_proj[2],pch=20)
text(x6G08_MD_3_proj[3],x6G08_MD_3_proj[2],pos=4,label="3",col="blue")

x6G08_MD_4_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_4.txt")
x6G08_MD_4_proj<-project.pca(x6G08_MD_4_raw,pcadata)
points(x6G08_MD_4_proj[3],x6G08_MD_4_proj[2],pch=20)
text(x6G08_MD_4_proj[3],x6G08_MD_4_proj[2],pos=4,label="4",col="blue")

x6G08_MD_5_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_5.txt")
x6G08_MD_5_proj<-project.pca(x6G08_MD_5_raw,pcadata)
points(x6G08_MD_5_proj[3],x6G08_MD_5_proj[2],pch=20)
text(x6G08_MD_5_proj[3],x6G08_MD_5_proj[2],pos=4,label="5",col="blue")

x6G08_MD_6_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_6.txt")
x6G08_MD_6_proj<-project.pca(x6G08_MD_6_raw,pcadata)
points(x6G08_MD_6_proj[3],x6G08_MD_6_proj[2],pch=20)
text(x6G08_MD_6_proj[3],x6G08_MD_6_proj[2],pos=4,label="6",col="black")

x6G08_MD_7_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_7.txt")
x6G08_MD_7_proj<-project.pca(x6G08_MD_7_raw,pcadata)
points(x6G08_MD_7_proj[3],x6G08_MD_7_proj[2],pch=20)
text(x6G08_MD_7_proj[3],x6G08_MD_7_proj[2],pos=4,label="7",col="blue")

x6G08_MD_8_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_8.txt")
x6G08_MD_8_proj<-project.pca(x6G08_MD_8_raw,pcadata)
points(x6G08_MD_8_proj[3],x6G08_MD_8_proj[2],pch=20)
text(x6G08_MD_8_proj[3],x6G08_MD_8_proj[2],pos=4,label="8",col="blue")

x6G08_MD_9_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_9.txt")
x6G08_MD_9_proj<-project.pca(x6G08_MD_9_raw,pcadata)
points(x6G08_MD_9_proj[3],x6G08_MD_9_proj[2],pch=20)
text(x6G08_MD_9_proj[3],x6G08_MD_9_proj[2],pos=4,label="9",col="green")

x6G08_MD_10_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_10.txt")
x6G08_MD_10_proj<-project.pca(x6G08_MD_10_raw,pcadata)
points(x6G08_MD_10_proj[3],x6G08_MD_10_proj[2],pch=20)
text(x6G08_MD_10_proj[3],x6G08_MD_10_proj[2],pos=4,label="10",col="red")

x6G08_MD_11_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_11.txt")
x6G08_MD_11_proj<-project.pca(x6G08_MD_11_raw,pcadata)
points(x6G08_MD_11_proj[3],x6G08_MD_11_proj[2],pch=20)
text(x6G08_MD_11_proj[3],x6G08_MD_11_proj[2],pos=4,label="11",col="blue")

x6G08_MD_12_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_12.txt")
x6G08_MD_12_proj<-project.pca(x6G08_MD_12_raw,pcadata)
points(x6G08_MD_12_proj[3],x6G08_MD_12_proj[2],pch=20)
text(x6G08_MD_12_proj[3],x6G08_MD_12_proj[2],pos=4,label="12",col="black")

x6G08_MD_13_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_13.txt")
x6G08_MD_13_proj<-project.pca(x6G08_MD_13_raw,pcadata)
points(x6G08_MD_13_proj[3],x6G08_MD_13_proj[2],pch=20)
text(x6G08_MD_13_proj[3],x6G08_MD_13_proj[2],pos=4,label="13",col="black")

x6G08_MD_14_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_14.txt")
x6G08_MD_14_proj<-project.pca(x6G08_MD_14_raw,pcadata)
points(x6G08_MD_14_proj[3],x6G08_MD_14_proj[2],pch=20)
text(x6G08_MD_14_proj[3],x6G08_MD_14_proj[2],pos=4,label="14",col="blue")

x6G08_MD_15_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_15.txt")
x6G08_MD_15_proj<-project.pca(x6G08_MD_15_raw,pcadata)
points(x6G08_MD_15_proj[3],x6G08_MD_15_proj[2],pch=20)
text(x6G08_MD_15_proj[3],x6G08_MD_15_proj[2],pos=4,label="15",col="black")

x6G08_MD_16_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_16.txt")
x6G08_MD_16_proj<-project.pca(x6G08_MD_16_raw,pcadata)
points(x6G08_MD_16_proj[3],x6G08_MD_16_proj[2],pch=20)
text(x6G08_MD_16_proj[3],x6G08_MD_16_proj[2],pos=4,label="16",col="black")

x6G08_MD_17_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_17.txt")
x6G08_MD_17_proj<-project.pca(x6G08_MD_17_raw,pcadata)
points(x6G08_MD_17_proj[3],x6G08_MD_17_proj[2],pch=20)
text(x6G08_MD_17_proj[3],x6G08_MD_17_proj[2],pos=4,label="17",col="blue")

x6G08_MD_18_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_18.txt")
x6G08_MD_18_proj<-project.pca(x6G08_MD_18_raw,pcadata)
points(x6G08_MD_18_proj[3],x6G08_MD_18_proj[2],pch=20)
text(x6G08_MD_18_proj[3],x6G08_MD_18_proj[2],pos=4,label="18",col="black")

x6G08_MD_19_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_19.txt")
x6G08_MD_19_proj<-project.pca(x6G08_MD_19_raw,pcadata)
points(x6G08_MD_19_proj[3],x6G08_MD_19_proj[2],pch=20)
text(x6G08_MD_19_proj[3],x6G08_MD_19_proj[2],pos=4,label="19",col="black")

x6G08_MD_20_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_20.txt")
x6G08_MD_20_proj<-project.pca(x6G08_MD_20_raw,pcadata)
points(x6G08_MD_20_proj[3],x6G08_MD_20_proj[2],pch=20)
text(x6G08_MD_20_proj[3],x6G08_MD_20_proj[2],pos=4,label="20",col="black")

x6G08_MD_21_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_21.txt")
x6G08_MD_21_proj<-project.pca(x6G08_MD_21_raw,pcadata)
points(x6G08_MD_21_proj[3],x6G08_MD_21_proj[2],pch=20)
text(x6G08_MD_21_proj[3],x6G08_MD_21_proj[2],pos=4,label="21",col="blue")

x6G08_MD_22_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_22.txt")
x6G08_MD_22_proj<-project.pca(x6G08_MD_22_raw,pcadata)
points(x6G08_MD_22_proj[3],x6G08_MD_22_proj[2],pch=20)
text(x6G08_MD_22_proj[3],x6G08_MD_22_proj[2],pos=4,label="22",col="blue")

x6G08_MD_23_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_23.txt")
x6G08_MD_23_proj<-project.pca(x6G08_MD_23_raw,pcadata)
points(x6G08_MD_23_proj[3],x6G08_MD_23_proj[2],pch=20)
text(x6G08_MD_23_proj[3],x6G08_MD_23_proj[2],pos=4,label="23",col="blue")

x6G08_MD_24_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_24.txt")
x6G08_MD_24_proj<-project.pca(x6G08_MD_24_raw,pcadata)
points(x6G08_MD_24_proj[3],x6G08_MD_24_proj[2],pch=20)
text(x6G08_MD_24_proj[3],x6G08_MD_24_proj[2],pos=4,label="24",col="blue")

x6G08_MD_25_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_25.txt")
x6G08_MD_25_proj<-project.pca(x6G08_MD_25_raw,pcadata)
points(x6G08_MD_25_proj[3],x6G08_MD_25_proj[2],pch=20)
text(x6G08_MD_25_proj[3],x6G08_MD_25_proj[2],pos=4,label="25",col="blue")

x6G08_MD_26_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_26.txt")
x6G08_MD_26_proj<-project.pca(x6G08_MD_26_raw,pcadata)
points(x6G08_MD_26_proj[3],x6G08_MD_26_proj[2],pch=20)
text(x6G08_MD_26_proj[3],x6G08_MD_26_proj[2],pos=4,label="26",col="blue")

x6G08_MD_27_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_27.txt")
x6G08_MD_27_proj<-project.pca(x6G08_MD_27_raw,pcadata)
points(x6G08_MD_27_proj[3],x6G08_MD_27_proj[2],pch=20)
text(x6G08_MD_27_proj[3],x6G08_MD_27_proj[2],pos=4,label="27",col="blue")

x6G08_MD_28_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_28.txt")
x6G08_MD_28_proj<-project.pca(x6G08_MD_28_raw,pcadata)
points(x6G08_MD_28_proj[3],x6G08_MD_28_proj[2],pch=20)
text(x6G08_MD_28_proj[3],x6G08_MD_28_proj[2],pos=4,label="28",col="blue")

x6G08_MD_29_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_29.txt")
x6G08_MD_29_proj<-project.pca(x6G08_MD_29_raw,pcadata)
points(x6G08_MD_29_proj[3],x6G08_MD_29_proj[2],pch=20)
text(x6G08_MD_29_proj[3],x6G08_MD_29_proj[2],pos=4,label="29",col="blue")

x6G08_MD_30_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_30.txt")
x6G08_MD_30_proj<-project.pca(x6G08_MD_30_raw,pcadata)
points(x6G08_MD_30_proj[3],x6G08_MD_30_proj[2],pch=20)
text(x6G08_MD_30_proj[3],x6G08_MD_30_proj[2],pos=4,label="30",col="black")

x6G08_MD_31_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_31.txt")
x6G08_MD_31_proj<-project.pca(x6G08_MD_31_raw,pcadata)
points(x6G08_MD_31_proj[3],x6G08_MD_31_proj[2],pch=20)
text(x6G08_MD_31_proj[3],x6G08_MD_31_proj[2],pos=4,label="31",col="black")

x6G08_MD_32_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_32.txt")
x6G08_MD_32_proj<-project.pca(x6G08_MD_32_raw,pcadata)
points(x6G08_MD_32_proj[3],x6G08_MD_32_proj[2],pch=20)
text(x6G08_MD_32_proj[3],x6G08_MD_32_proj[2],pos=4,label="32",col="black")

x6G08_MD_33_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_33.txt")
x6G08_MD_33_proj<-project.pca(x6G08_MD_33_raw,pcadata)
points(x6G08_MD_33_proj[3],x6G08_MD_33_proj[2],pch=20)
text(x6G08_MD_33_proj[3],x6G08_MD_33_proj[2],pos=4,label="33",col="black")

x6G08_MD_34_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_34.txt")
x6G08_MD_34_proj<-project.pca(x6G08_MD_34_raw,pcadata)
points(x6G08_MD_34_proj[3],x6G08_MD_34_proj[2],pch=20)
text(x6G08_MD_34_proj[3],x6G08_MD_34_proj[2],pos=4,label="34",col="blue")

x6G08_MD_35_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_35.txt")
x6G08_MD_35_proj<-project.pca(x6G08_MD_35_raw,pcadata)
points(x6G08_MD_35_proj[3],x6G08_MD_35_proj[2],pch=20)
text(x6G08_MD_35_proj[3],x6G08_MD_35_proj[2],pos=4,label="35",col="black")

x6G08_MD_36_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_36.txt")
x6G08_MD_36_proj<-project.pca(x6G08_MD_36_raw,pcadata)
points(x6G08_MD_36_proj[3],x6G08_MD_36_proj[2],pch=20)
text(x6G08_MD_36_proj[3],x6G08_MD_36_proj[2],pos=4,label="36",col="black")

x6G08_MD_37_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_37.txt")
x6G08_MD_37_proj<-project.pca(x6G08_MD_37_raw,pcadata)
points(x6G08_MD_37_proj[3],x6G08_MD_37_proj[2],pch=20)
text(x6G08_MD_37_proj[3],x6G08_MD_37_proj[2],pos=4,label="37",col="blue")

x6G08_MD_38_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_38.txt")
x6G08_MD_38_proj<-project.pca(x6G08_MD_38_raw,pcadata)
points(x6G08_MD_38_proj[3],x6G08_MD_38_proj[2],pch=20)
text(x6G08_MD_38_proj[3],x6G08_MD_38_proj[2],pos=4,label="38",col="blue")

x6G08_MD_39_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_39.txt")
x6G08_MD_39_proj<-project.pca(x6G08_MD_39_raw,pcadata)
points(x6G08_MD_39_proj[3],x6G08_MD_39_proj[2],pch=20)
text(x6G08_MD_39_proj[3],x6G08_MD_39_proj[2],pos=4,label="39",col="blue")

x6G08_MD_40_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_40.txt")
x6G08_MD_40_proj<-project.pca(x6G08_MD_40_raw,pcadata)
points(x6G08_MD_40_proj[3],x6G08_MD_40_proj[2],pch=20)
text(x6G08_MD_40_proj[3],x6G08_MD_40_proj[2],pos=4,label="40",col="black")

x6G08_MD_41_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_41.txt")
x6G08_MD_41_proj<-project.pca(x6G08_MD_41_raw,pcadata)
points(x6G08_MD_41_proj[3],x6G08_MD_41_proj[2],pch=20)
text(x6G08_MD_41_proj[3],x6G08_MD_41_proj[2],pos=4,label="41",col="black")

x6G08_MD_42_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_42.txt")
x6G08_MD_42_proj<-project.pca(x6G08_MD_42_raw,pcadata)
points(x6G08_MD_42_proj[3],x6G08_MD_42_proj[2],pch=20)
text(x6G08_MD_42_proj[3],x6G08_MD_42_proj[2],pos=4,label="42",col="black")

x6G08_MD_43_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_43.txt")
x6G08_MD_43_proj<-project.pca(x6G08_MD_43_raw,pcadata)
points(x6G08_MD_43_proj[3],x6G08_MD_43_proj[2],pch=20)
text(x6G08_MD_43_proj[3],x6G08_MD_43_proj[2],pos=4,label="43",col="black")

x6G08_MD_44_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_44.txt")
x6G08_MD_44_proj<-project.pca(x6G08_MD_44_raw,pcadata)
points(x6G08_MD_44_proj[3],x6G08_MD_44_proj[2],pch=20)
text(x6G08_MD_44_proj[3],x6G08_MD_44_proj[2],pos=4,label="44",col="black")

x6G08_MD_45_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_45.txt")
x6G08_MD_45_proj<-project.pca(x6G08_MD_45_raw,pcadata)
points(x6G08_MD_45_proj[3],x6G08_MD_45_proj[2],pch=20)
text(x6G08_MD_45_proj[3],x6G08_MD_45_proj[2],pos=4,label="45",col="black")

x6G08_MD_46_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_46.txt")
x6G08_MD_46_proj<-project.pca(x6G08_MD_46_raw,pcadata)
points(x6G08_MD_46_proj[3],x6G08_MD_46_proj[2],pch=20)
text(x6G08_MD_46_proj[3],x6G08_MD_46_proj[2],pos=4,label="46",col="black")

x6G08_MD_47_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_47.txt")
x6G08_MD_47_proj<-project.pca(x6G08_MD_47_raw,pcadata)
points(x6G08_MD_47_proj[3],x6G08_MD_47_proj[2],pch=20)
text(x6G08_MD_47_proj[3],x6G08_MD_47_proj[2],pos=4,label="47",col="black")

x6G08_MD_48_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_48.txt")
x6G08_MD_48_proj<-project.pca(x6G08_MD_48_raw,pcadata)
points(x6G08_MD_48_proj[3],x6G08_MD_48_proj[2],pch=20)
text(x6G08_MD_48_proj[3],x6G08_MD_48_proj[2],pos=4,label="48",col="black")

x6G08_MD_49_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_49.txt")
x6G08_MD_49_proj<-project.pca(x6G08_MD_49_raw,pcadata)
points(x6G08_MD_49_proj[3],x6G08_MD_49_proj[2],pch=20)
text(x6G08_MD_49_proj[3],x6G08_MD_49_proj[2],pos=4,label="49",col="black")

x6G08_MD_50_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_50.txt")
x6G08_MD_50_proj<-project.pca(x6G08_MD_50_raw,pcadata)
points(x6G08_MD_50_proj[3],x6G08_MD_50_proj[2],pch=20)
text(x6G08_MD_50_proj[3],x6G08_MD_50_proj[2],pos=4,label="50",col="black")

x6G08_MD_51_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_51.txt")
x6G08_MD_51_proj<-project.pca(x6G08_MD_51_raw,pcadata)
points(x6G08_MD_51_proj[3],x6G08_MD_51_proj[2],pch=20)
text(x6G08_MD_51_proj[3],x6G08_MD_51_proj[2],pos=4,label="51",col="black")

x6G08_MD_52_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_52.txt")
x6G08_MD_52_proj<-project.pca(x6G08_MD_52_raw,pcadata)
points(x6G08_MD_52_proj[3],x6G08_MD_52_proj[2],pch=20)
text(x6G08_MD_52_proj[3],x6G08_MD_52_proj[2],pos=4,label="52",col="black")

x6G08_MD_53_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_53.txt")
x6G08_MD_53_proj<-project.pca(x6G08_MD_53_raw,pcadata)
points(x6G08_MD_53_proj[3],x6G08_MD_53_proj[2],pch=20)
text(x6G08_MD_53_proj[3],x6G08_MD_53_proj[2],pos=4,label="53",col="black")

x6G08_MD_54_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_54.txt")
x6G08_MD_54_proj<-project.pca(x6G08_MD_54_raw,pcadata)
points(x6G08_MD_54_proj[3],x6G08_MD_54_proj[2],pch=20)
text(x6G08_MD_54_proj[3],x6G08_MD_54_proj[2],pos=4,label="54",col="black")

x6G08_MD_55_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_55.txt")
x6G08_MD_55_proj<-project.pca(x6G08_MD_55_raw,pcadata)
points(x6G08_MD_55_proj[3],x6G08_MD_55_proj[2],pch=20)
text(x6G08_MD_55_proj[3],x6G08_MD_55_proj[2],pos=4,label="55",col="black")

x6G08_MD_56_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_56.txt")
x6G08_MD_56_proj<-project.pca(x6G08_MD_56_raw,pcadata)
points(x6G08_MD_56_proj[3],x6G08_MD_56_proj[2],pch=20)
text(x6G08_MD_56_proj[3],x6G08_MD_56_proj[2],pos=4,label="56",col="black")

x6G08_MD_57_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_57.txt")
x6G08_MD_57_proj<-project.pca(x6G08_MD_57_raw,pcadata)
points(x6G08_MD_57_proj[3],x6G08_MD_57_proj[2],pch=20)
text(x6G08_MD_57_proj[3],x6G08_MD_57_proj[2],pos=4,label="57",col="black")

x6G08_MD_58_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_58.txt")
x6G08_MD_58_proj<-project.pca(x6G08_MD_58_raw,pcadata)
points(x6G08_MD_58_proj[3],x6G08_MD_58_proj[2],pch=20)
text(x6G08_MD_58_proj[3],x6G08_MD_58_proj[2],pos=4,label="58",col="blue")

x6G08_MD_59_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_59.txt")
x6G08_MD_59_proj<-project.pca(x6G08_MD_59_raw,pcadata)
points(x6G08_MD_59_proj[3],x6G08_MD_59_proj[2],pch=20)
text(x6G08_MD_59_proj[3],x6G08_MD_59_proj[2],pos=4,label="59",col="black")

x6G08_MD_60_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_60.txt")
x6G08_MD_60_proj<-project.pca(x6G08_MD_60_raw,pcadata)
points(x6G08_MD_60_proj[3],x6G08_MD_60_proj[2],pch=20)
text(x6G08_MD_60_proj[3],x6G08_MD_60_proj[2],pos=4,label="60",col="black")

x6G08_MD_61_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_61.txt")
x6G08_MD_61_proj<-project.pca(x6G08_MD_61_raw,pcadata)
points(x6G08_MD_61_proj[3],x6G08_MD_61_proj[2],pch=20)
text(x6G08_MD_61_proj[3],x6G08_MD_61_proj[2],pos=4,label="61",col="black")

x6G08_MD_62_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_62.txt")
x6G08_MD_62_proj<-project.pca(x6G08_MD_62_raw,pcadata)
points(x6G08_MD_62_proj[3],x6G08_MD_62_proj[2],pch=20)
text(x6G08_MD_62_proj[3],x6G08_MD_62_proj[2],pos=4,label="62",col="black")

x6G08_MD_63_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_63.txt")
x6G08_MD_63_proj<-project.pca(x6G08_MD_63_raw,pcadata)
points(x6G08_MD_63_proj[3],x6G08_MD_63_proj[2],pch=20)
text(x6G08_MD_63_proj[3],x6G08_MD_63_proj[2],pos=4,label="63",col="black")

x6G08_MD_64_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_64.txt")
x6G08_MD_64_proj<-project.pca(x6G08_MD_64_raw,pcadata)
points(x6G08_MD_64_proj[3],x6G08_MD_64_proj[2],pch=20)
text(x6G08_MD_64_proj[3],x6G08_MD_64_proj[2],pos=4,label="64",col="black")

x6G08_MD_65_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_65.txt")
x6G08_MD_65_proj<-project.pca(x6G08_MD_65_raw,pcadata)
points(x6G08_MD_65_proj[3],x6G08_MD_65_proj[2],pch=20)
text(x6G08_MD_65_proj[3],x6G08_MD_65_proj[2],pos=4,label="65",col="black")

x6G08_MD_66_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_66.txt")
x6G08_MD_66_proj<-project.pca(x6G08_MD_66_raw,pcadata)
points(x6G08_MD_66_proj[3],x6G08_MD_66_proj[2],pch=20)
text(x6G08_MD_66_proj[3],x6G08_MD_66_proj[2],pos=4,label="66",col="black")

x6G08_MD_67_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_67.txt")
x6G08_MD_67_proj<-project.pca(x6G08_MD_67_raw,pcadata)
points(x6G08_MD_67_proj[3],x6G08_MD_67_proj[2],pch=20)
text(x6G08_MD_67_proj[3],x6G08_MD_67_proj[2],pos=4,label="67",col="green")

x6G08_MD_68_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_68.txt")
x6G08_MD_68_proj<-project.pca(x6G08_MD_68_raw,pcadata)
points(x6G08_MD_68_proj[3],x6G08_MD_68_proj[2],pch=20)
text(x6G08_MD_68_proj[3],x6G08_MD_68_proj[2],pos=4,label="68",col="blue")

x6G08_MD_69_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_69.txt")
x6G08_MD_69_proj<-project.pca(x6G08_MD_69_raw,pcadata)
points(x6G08_MD_69_proj[3],x6G08_MD_69_proj[2],pch=20)
text(x6G08_MD_69_proj[3],x6G08_MD_69_proj[2],pos=4,label="69",col="blue")

x6G08_MD_70_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_70.txt")
x6G08_MD_70_proj<-project.pca(x6G08_MD_70_raw,pcadata)
points(x6G08_MD_70_proj[3],x6G08_MD_70_proj[2],pch=20)
text(x6G08_MD_70_proj[3],x6G08_MD_70_proj[2],pos=4,label="70",col="black")

x6G08_MD_71_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_71.txt")
x6G08_MD_71_proj<-project.pca(x6G08_MD_71_raw,pcadata)
points(x6G08_MD_71_proj[3],x6G08_MD_71_proj[2],pch=20)
text(x6G08_MD_71_proj[3],x6G08_MD_71_proj[2],pos=4,label="71",col="black")

x6G08_MD_72_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_72.txt")
x6G08_MD_72_proj<-project.pca(x6G08_MD_72_raw,pcadata)
points(x6G08_MD_72_proj[3],x6G08_MD_72_proj[2],pch=20)
text(x6G08_MD_72_proj[3],x6G08_MD_72_proj[2],pos=4,label="72",col="black")

x6G08_MD_73_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_73.txt")
x6G08_MD_73_proj<-project.pca(x6G08_MD_73_raw,pcadata)
points(x6G08_MD_73_proj[3],x6G08_MD_73_proj[2],pch=20)
text(x6G08_MD_73_proj[3],x6G08_MD_73_proj[2],pos=4,label="73",col="blue")

x6G08_MD_74_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_74.txt")
x6G08_MD_74_proj<-project.pca(x6G08_MD_74_raw,pcadata)
points(x6G08_MD_74_proj[3],x6G08_MD_74_proj[2],pch=20)
text(x6G08_MD_74_proj[3],x6G08_MD_74_proj[2],pos=4,label="74",col="blue")

x6G08_MD_75_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_75.txt")
x6G08_MD_75_proj<-project.pca(x6G08_MD_75_raw,pcadata)
points(x6G08_MD_75_proj[3],x6G08_MD_75_proj[2],pch=20)
text(x6G08_MD_75_proj[3],x6G08_MD_75_proj[2],pos=4,label="75",col="black")

x6G08_MD_76_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_76.txt")
x6G08_MD_76_proj<-project.pca(x6G08_MD_76_raw,pcadata)
points(x6G08_MD_76_proj[3],x6G08_MD_76_proj[2],pch=20)
text(x6G08_MD_76_proj[3],x6G08_MD_76_proj[2],pos=4,label="76",col="black")

x6G08_MD_77_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_77.txt")
x6G08_MD_77_proj<-project.pca(x6G08_MD_77_raw,pcadata)
points(x6G08_MD_77_proj[3],x6G08_MD_77_proj[2],pch=20)
text(x6G08_MD_77_proj[3],x6G08_MD_77_proj[2],pos=4,label="77",col="black")

x6G08_MD_78_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_78.txt")
x6G08_MD_78_proj<-project.pca(x6G08_MD_78_raw,pcadata)
points(x6G08_MD_78_proj[3],x6G08_MD_78_proj[2],pch=20)
text(x6G08_MD_78_proj[3],x6G08_MD_78_proj[2],pos=4,label="78",col="blue")

x6G08_MD_79_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_79.txt")
x6G08_MD_79_proj<-project.pca(x6G08_MD_79_raw,pcadata)
points(x6G08_MD_79_proj[3],x6G08_MD_79_proj[2],pch=20)
text(x6G08_MD_79_proj[3],x6G08_MD_79_proj[2],pos=4,label="79",col="green")

x6G08_MD_80_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_80.txt")
x6G08_MD_80_proj<-project.pca(x6G08_MD_80_raw,pcadata)
points(x6G08_MD_80_proj[3],x6G08_MD_80_proj[2],pch=20)
text(x6G08_MD_80_proj[3],x6G08_MD_80_proj[2],pos=4,label="80",col="blue")

x6G08_MD_81_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_81.txt")
x6G08_MD_81_proj<-project.pca(x6G08_MD_81_raw,pcadata)
points(x6G08_MD_81_proj[3],x6G08_MD_81_proj[2],pch=20)
text(x6G08_MD_81_proj[3],x6G08_MD_81_proj[2],pos=4,label="81",col="black")

x6G08_MD_82_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_82.txt")
x6G08_MD_82_proj<-project.pca(x6G08_MD_82_raw,pcadata)
points(x6G08_MD_82_proj[3],x6G08_MD_82_proj[2],pch=20)
text(x6G08_MD_82_proj[3],x6G08_MD_82_proj[2],pos=4,label="82",col="black")

x6G08_MD_83_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_83.txt")
x6G08_MD_83_proj<-project.pca(x6G08_MD_83_raw,pcadata)
points(x6G08_MD_83_proj[3],x6G08_MD_83_proj[2],pch=20)
text(x6G08_MD_83_proj[3],x6G08_MD_83_proj[2],pos=4,label="83",col="green")

x6G08_MD_84_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_84.txt")
x6G08_MD_84_proj<-project.pca(x6G08_MD_84_raw,pcadata)
points(x6G08_MD_84_proj[3],x6G08_MD_84_proj[2],pch=20)
text(x6G08_MD_84_proj[3],x6G08_MD_84_proj[2],pos=4,label="84",col="blue")

x6G08_MD_85_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_85.txt")
x6G08_MD_85_proj<-project.pca(x6G08_MD_85_raw,pcadata)
points(x6G08_MD_85_proj[3],x6G08_MD_85_proj[2],pch=20)
text(x6G08_MD_85_proj[3],x6G08_MD_85_proj[2],pos=4,label="85",col="black")

x6G08_MD_86_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_86.txt")
x6G08_MD_86_proj<-project.pca(x6G08_MD_86_raw,pcadata)
points(x6G08_MD_86_proj[3],x6G08_MD_86_proj[2],pch=20)
text(x6G08_MD_86_proj[3],x6G08_MD_86_proj[2],pos=4,label="86",col="blue")

x6G08_MD_87_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_87.txt")
x6G08_MD_87_proj<-project.pca(x6G08_MD_87_raw,pcadata)
points(x6G08_MD_87_proj[3],x6G08_MD_87_proj[2],pch=20)
text(x6G08_MD_87_proj[3],x6G08_MD_87_proj[2],pos=4,label="87",col="blue")

x6G08_MD_88_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_88.txt")
x6G08_MD_88_proj<-project.pca(x6G08_MD_88_raw,pcadata)
points(x6G08_MD_88_proj[3],x6G08_MD_88_proj[2],pch=20)
text(x6G08_MD_88_proj[3],x6G08_MD_88_proj[2],pos=4,label="88",col="green")

x6G08_MD_89_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_89.txt")
x6G08_MD_89_proj<-project.pca(x6G08_MD_89_raw,pcadata)
points(x6G08_MD_89_proj[3],x6G08_MD_89_proj[2],pch=20)
text(x6G08_MD_89_proj[3],x6G08_MD_89_proj[2],pos=4,label="89",col="green")

x6G08_MD_90_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_90.txt")
x6G08_MD_90_proj<-project.pca(x6G08_MD_90_raw,pcadata)
points(x6G08_MD_90_proj[3],x6G08_MD_90_proj[2],pch=20)
text(x6G08_MD_90_proj[3],x6G08_MD_90_proj[2],pos=4,label="90",col="green")

x6G08_MD_91_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_91.txt")
x6G08_MD_91_proj<-project.pca(x6G08_MD_91_raw,pcadata)
points(x6G08_MD_91_proj[3],x6G08_MD_91_proj[2],pch=20)
text(x6G08_MD_91_proj[3],x6G08_MD_91_proj[2],pos=4,label="91",col="green")

x6G08_MD_92_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_92.txt")
x6G08_MD_92_proj<-project.pca(x6G08_MD_92_raw,pcadata)
points(x6G08_MD_92_proj[3],x6G08_MD_92_proj[2],pch=20)
text(x6G08_MD_92_proj[3],x6G08_MD_92_proj[2],pos=4,label="92",col="green")

x6G08_MD_93_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_93.txt")
x6G08_MD_93_proj<-project.pca(x6G08_MD_93_raw,pcadata)
points(x6G08_MD_93_proj[3],x6G08_MD_93_proj[2],pch=20)
text(x6G08_MD_93_proj[3],x6G08_MD_93_proj[2],pos=4,label="93",col="blue")

x6G08_MD_94_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_94.txt")
x6G08_MD_94_proj<-project.pca(x6G08_MD_94_raw,pcadata)
points(x6G08_MD_94_proj[3],x6G08_MD_94_proj[2],pch=20)
text(x6G08_MD_94_proj[3],x6G08_MD_94_proj[2],pos=4,label="94",col="black")

x6G08_MD_95_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_95.txt")
x6G08_MD_95_proj<-project.pca(x6G08_MD_95_raw,pcadata)
points(x6G08_MD_95_proj[3],x6G08_MD_95_proj[2],pch=20)
text(x6G08_MD_95_proj[3],x6G08_MD_95_proj[2],pos=4,label="95",col="green")

x6G08_MD_96_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_96.txt")
x6G08_MD_96_proj<-project.pca(x6G08_MD_96_raw,pcadata)
points(x6G08_MD_96_proj[3],x6G08_MD_96_proj[2],pch=20)
text(x6G08_MD_96_proj[3],x6G08_MD_96_proj[2],pos=4,label="96",col="green")

x6G08_MD_97_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_97.txt")
x6G08_MD_97_proj<-project.pca(x6G08_MD_97_raw,pcadata)
points(x6G08_MD_97_proj[3],x6G08_MD_97_proj[2],pch=20)
text(x6G08_MD_97_proj[3],x6G08_MD_97_proj[2],pos=4,label="97",col="green")

x6G08_MD_98_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_98.txt")
x6G08_MD_98_proj<-project.pca(x6G08_MD_98_raw,pcadata)
points(x6G08_MD_98_proj[3],x6G08_MD_98_proj[2],pch=20)
text(x6G08_MD_98_proj[3],x6G08_MD_98_proj[2],pos=4,label="98",col="black")

x6G08_MD_99_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_99.txt")
x6G08_MD_99_proj<-project.pca(x6G08_MD_99_raw,pcadata)
points(x6G08_MD_99_proj[3],x6G08_MD_99_proj[2],pch=20)
text(x6G08_MD_99_proj[3],x6G08_MD_99_proj[2],pos=4,label="99",col="blue")

x6G08_MD_100_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_100.txt")
x6G08_MD_100_proj<-project.pca(x6G08_MD_100_raw,pcadata)
points(x6G08_MD_100_proj[3],x6G08_MD_100_proj[2],pch=20)
text(x6G08_MD_100_proj[3],x6G08_MD_100_proj[2],pos=4,label="100",col="black")

x6G08_MD_101_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_101.txt")
x6G08_MD_101_proj<-project.pca(x6G08_MD_101_raw,pcadata)
points(x6G08_MD_101_proj[3],x6G08_MD_101_proj[2],pch=20)
text(x6G08_MD_101_proj[3],x6G08_MD_101_proj[2],pos=4,label="101",col="black")

x6G08_MD_102_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_102.txt")
x6G08_MD_102_proj<-project.pca(x6G08_MD_102_raw,pcadata)
points(x6G08_MD_102_proj[3],x6G08_MD_102_proj[2],pch=20)
text(x6G08_MD_102_proj[3],x6G08_MD_102_proj[2],pos=4,label="102",col="green")

x6G08_MD_103_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_103.txt")
x6G08_MD_103_proj<-project.pca(x6G08_MD_103_raw,pcadata)
points(x6G08_MD_103_proj[3],x6G08_MD_103_proj[2],pch=20)
text(x6G08_MD_103_proj[3],x6G08_MD_103_proj[2],pos=4,label="103",col="blue")

x6G08_MD_104_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_104.txt")
x6G08_MD_104_proj<-project.pca(x6G08_MD_104_raw,pcadata)
points(x6G08_MD_104_proj[3],x6G08_MD_104_proj[2],pch=20)
text(x6G08_MD_104_proj[3],x6G08_MD_104_proj[2],pos=4,label="104",col="blue")

x6G08_MD_105_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_105.txt")
x6G08_MD_105_proj<-project.pca(x6G08_MD_105_raw,pcadata)
points(x6G08_MD_105_proj[3],x6G08_MD_105_proj[2],pch=20)
text(x6G08_MD_105_proj[3],x6G08_MD_105_proj[2],pos=4,label="105",col="blue")

x6G08_MD_106_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_106.txt")
x6G08_MD_106_proj<-project.pca(x6G08_MD_106_raw,pcadata)
points(x6G08_MD_106_proj[3],x6G08_MD_106_proj[2],pch=20)
text(x6G08_MD_106_proj[3],x6G08_MD_106_proj[2],pos=4,label="106",col="black")

x6G08_MD_107_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_107.txt")
x6G08_MD_107_proj<-project.pca(x6G08_MD_107_raw,pcadata)
points(x6G08_MD_107_proj[3],x6G08_MD_107_proj[2],pch=20)
text(x6G08_MD_107_proj[3],x6G08_MD_107_proj[2],pos=4,label="107",col="green")

x6G08_MD_108_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_108.txt")
x6G08_MD_108_proj<-project.pca(x6G08_MD_108_raw,pcadata)
points(x6G08_MD_108_proj[3],x6G08_MD_108_proj[2],pch=20)
text(x6G08_MD_108_proj[3],x6G08_MD_108_proj[2],pos=4,label="108",col="black")

x6G08_MD_109_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_109.txt")
x6G08_MD_109_proj<-project.pca(x6G08_MD_109_raw,pcadata)
points(x6G08_MD_109_proj[3],x6G08_MD_109_proj[2],pch=20)
text(x6G08_MD_109_proj[3],x6G08_MD_109_proj[2],pos=4,label="109",col="black")

x6G08_MD_110_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_110.txt")
x6G08_MD_110_proj<-project.pca(x6G08_MD_110_raw,pcadata)
points(x6G08_MD_110_proj[3],x6G08_MD_110_proj[2],pch=20)
text(x6G08_MD_110_proj[3],x6G08_MD_110_proj[2],pos=4,label="110",col="black")

x6G08_MD_111_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_111.txt")
x6G08_MD_111_proj<-project.pca(x6G08_MD_111_raw,pcadata)
points(x6G08_MD_111_proj[3],x6G08_MD_111_proj[2],pch=20)
text(x6G08_MD_111_proj[3],x6G08_MD_111_proj[2],pos=4,label="111",col="black")

x6G08_MD_112_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_112.txt")
x6G08_MD_112_proj<-project.pca(x6G08_MD_112_raw,pcadata)
points(x6G08_MD_112_proj[3],x6G08_MD_112_proj[2],pch=20)
text(x6G08_MD_112_proj[3],x6G08_MD_112_proj[2],pos=4,label="112",col="blue")

x6G08_MD_113_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_113.txt")
x6G08_MD_113_proj<-project.pca(x6G08_MD_113_raw,pcadata)
points(x6G08_MD_113_proj[3],x6G08_MD_113_proj[2],pch=20)
text(x6G08_MD_113_proj[3],x6G08_MD_113_proj[2],pos=4,label="113",col="blue")

x6G08_MD_114_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_114.txt")
x6G08_MD_114_proj<-project.pca(x6G08_MD_114_raw,pcadata)
points(x6G08_MD_114_proj[3],x6G08_MD_114_proj[2],pch=20)
text(x6G08_MD_114_proj[3],x6G08_MD_114_proj[2],pos=4,label="114",col="blue")

x6G08_MD_115_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_115.txt")
x6G08_MD_115_proj<-project.pca(x6G08_MD_115_raw,pcadata)
points(x6G08_MD_115_proj[3],x6G08_MD_115_proj[2],pch=20)
text(x6G08_MD_115_proj[3],x6G08_MD_115_proj[2],pos=4,label="115",col="black")

x6G08_MD_116_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_116.txt")
x6G08_MD_116_proj<-project.pca(x6G08_MD_116_raw,pcadata)
points(x6G08_MD_116_proj[3],x6G08_MD_116_proj[2],pch=20)
text(x6G08_MD_116_proj[3],x6G08_MD_116_proj[2],pos=4,label="116",col="black")

x6G08_MD_117_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_117.txt")
x6G08_MD_117_proj<-project.pca(x6G08_MD_117_raw,pcadata)
points(x6G08_MD_117_proj[3],x6G08_MD_117_proj[2],pch=20)
text(x6G08_MD_117_proj[3],x6G08_MD_117_proj[2],pos=4,label="117",col="blue")

x6G08_MD_118_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_118.txt")
x6G08_MD_118_proj<-project.pca(x6G08_MD_118_raw,pcadata)
points(x6G08_MD_118_proj[3],x6G08_MD_118_proj[2],pch=20)
text(x6G08_MD_118_proj[3],x6G08_MD_118_proj[2],pos=4,label="118",col="green")

x6G08_MD_119_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_119.txt")
x6G08_MD_119_proj<-project.pca(x6G08_MD_119_raw,pcadata)
points(x6G08_MD_119_proj[3],x6G08_MD_119_proj[2],pch=20)
text(x6G08_MD_119_proj[3],x6G08_MD_119_proj[2],pos=4,label="119",col="blue")

x6G08_MD_120_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_120.txt")
x6G08_MD_120_proj<-project.pca(x6G08_MD_120_raw,pcadata)
points(x6G08_MD_120_proj[3],x6G08_MD_120_proj[2],pch=20)
text(x6G08_MD_120_proj[3],x6G08_MD_120_proj[2],pos=4,label="120",col="black")

x6G08_MD_121_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_121.txt")
x6G08_MD_121_proj<-project.pca(x6G08_MD_121_raw,pcadata)
points(x6G08_MD_121_proj[3],x6G08_MD_121_proj[2],pch=20)
text(x6G08_MD_121_proj[3],x6G08_MD_121_proj[2],pos=4,label="121",col="black")

x6G08_MD_122_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_122.txt")
x6G08_MD_122_proj<-project.pca(x6G08_MD_122_raw,pcadata)
points(x6G08_MD_122_proj[3],x6G08_MD_122_proj[2],pch=20)
text(x6G08_MD_122_proj[3],x6G08_MD_122_proj[2],pos=4,label="122",col="black")

x6G08_MD_123_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_123.txt")
x6G08_MD_123_proj<-project.pca(x6G08_MD_123_raw,pcadata)
points(x6G08_MD_123_proj[3],x6G08_MD_123_proj[2],pch=20)
text(x6G08_MD_123_proj[3],x6G08_MD_123_proj[2],pos=4,label="123",col="black")

x6G08_MD_124_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_124.txt")
x6G08_MD_124_proj<-project.pca(x6G08_MD_124_raw,pcadata)
points(x6G08_MD_124_proj[3],x6G08_MD_124_proj[2],pch=20)
text(x6G08_MD_124_proj[3],x6G08_MD_124_proj[2],pos=4,label="124",col="black")

x6G08_MD_125_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_125.txt")
x6G08_MD_125_proj<-project.pca(x6G08_MD_125_raw,pcadata)
points(x6G08_MD_125_proj[3],x6G08_MD_125_proj[2],pch=20)
text(x6G08_MD_125_proj[3],x6G08_MD_125_proj[2],pos=4,label="125",col="blue")

x6G08_MD_126_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_126.txt")
x6G08_MD_126_proj<-project.pca(x6G08_MD_126_raw,pcadata)
points(x6G08_MD_126_proj[3],x6G08_MD_126_proj[2],pch=20)
text(x6G08_MD_126_proj[3],x6G08_MD_126_proj[2],pos=4,label="126",col="black")

x6G08_MD_127_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_127.txt")
x6G08_MD_127_proj<-project.pca(x6G08_MD_127_raw,pcadata)
points(x6G08_MD_127_proj[3],x6G08_MD_127_proj[2],pch=20)
text(x6G08_MD_127_proj[3],x6G08_MD_127_proj[2],pos=4,label="127",col="black")

x6G08_MD_128_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_128.txt")
x6G08_MD_128_proj<-project.pca(x6G08_MD_128_raw,pcadata)
points(x6G08_MD_128_proj[3],x6G08_MD_128_proj[2],pch=20)
text(x6G08_MD_128_proj[3],x6G08_MD_128_proj[2],pos=4,label="128",col="black")

x6G08_MD_129_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_129.txt")
x6G08_MD_129_proj<-project.pca(x6G08_MD_129_raw,pcadata)
points(x6G08_MD_129_proj[3],x6G08_MD_129_proj[2],pch=20)
text(x6G08_MD_129_proj[3],x6G08_MD_129_proj[2],pos=4,label="129",col="black")

x6G08_MD_130_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_130.txt")
x6G08_MD_130_proj<-project.pca(x6G08_MD_130_raw,pcadata)
points(x6G08_MD_130_proj[3],x6G08_MD_130_proj[2],pch=20)
text(x6G08_MD_130_proj[3],x6G08_MD_130_proj[2],pos=4,label="130",col="black")

x6G08_MD_131_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_131.txt")
x6G08_MD_131_proj<-project.pca(x6G08_MD_131_raw,pcadata)
points(x6G08_MD_131_proj[3],x6G08_MD_131_proj[2],pch=20)
text(x6G08_MD_131_proj[3],x6G08_MD_131_proj[2],pos=4,label="131",col="black")

x6G08_MD_132_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_132.txt")
x6G08_MD_132_proj<-project.pca(x6G08_MD_132_raw,pcadata)
points(x6G08_MD_132_proj[3],x6G08_MD_132_proj[2],pch=20)
text(x6G08_MD_132_proj[3],x6G08_MD_132_proj[2],pos=4,label="132",col="black")

x6G08_MD_133_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_133.txt")
x6G08_MD_133_proj<-project.pca(x6G08_MD_133_raw,pcadata)
points(x6G08_MD_133_proj[3],x6G08_MD_133_proj[2],pch=20)
text(x6G08_MD_133_proj[3],x6G08_MD_133_proj[2],pos=4,label="133",col="black")

x6G08_MD_134_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_134.txt")
x6G08_MD_134_proj<-project.pca(x6G08_MD_134_raw,pcadata)
points(x6G08_MD_134_proj[3],x6G08_MD_134_proj[2],pch=20)
text(x6G08_MD_134_proj[3],x6G08_MD_134_proj[2],pos=4,label="134",col="black")

x6G08_MD_135_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_135.txt")
x6G08_MD_135_proj<-project.pca(x6G08_MD_135_raw,pcadata)
points(x6G08_MD_135_proj[3],x6G08_MD_135_proj[2],pch=20)
text(x6G08_MD_135_proj[3],x6G08_MD_135_proj[2],pos=4,label="135",col="black")

x6G08_MD_136_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_136.txt")
x6G08_MD_136_proj<-project.pca(x6G08_MD_136_raw,pcadata)
points(x6G08_MD_136_proj[3],x6G08_MD_136_proj[2],pch=20)
text(x6G08_MD_136_proj[3],x6G08_MD_136_proj[2],pos=4,label="136",col="black")

x6G08_MD_137_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_137.txt")
x6G08_MD_137_proj<-project.pca(x6G08_MD_137_raw,pcadata)
points(x6G08_MD_137_proj[3],x6G08_MD_137_proj[2],pch=20)
text(x6G08_MD_137_proj[3],x6G08_MD_137_proj[2],pos=4,label="137",col="black")

x6G08_MD_138_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_138.txt")
x6G08_MD_138_proj<-project.pca(x6G08_MD_138_raw,pcadata)
points(x6G08_MD_138_proj[3],x6G08_MD_138_proj[2],pch=20)
text(x6G08_MD_138_proj[3],x6G08_MD_138_proj[2],pos=4,label="138",col="black")

x6G08_MD_139_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_139.txt")
x6G08_MD_139_proj<-project.pca(x6G08_MD_139_raw,pcadata)
points(x6G08_MD_139_proj[3],x6G08_MD_139_proj[2],pch=20)
text(x6G08_MD_139_proj[3],x6G08_MD_139_proj[2],pos=4,label="139",col="black")

x6G08_MD_140_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_140.txt")
x6G08_MD_140_proj<-project.pca(x6G08_MD_140_raw,pcadata)
points(x6G08_MD_140_proj[3],x6G08_MD_140_proj[2],pch=20)
text(x6G08_MD_140_proj[3],x6G08_MD_140_proj[2],pos=4,label="140",col="black")

x6G08_MD_141_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_141.txt")
x6G08_MD_141_proj<-project.pca(x6G08_MD_141_raw,pcadata)
points(x6G08_MD_141_proj[3],x6G08_MD_141_proj[2],pch=20)
text(x6G08_MD_141_proj[3],x6G08_MD_141_proj[2],pos=4,label="141",col="green")

x6G08_MD_142_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_142.txt")
x6G08_MD_142_proj<-project.pca(x6G08_MD_142_raw,pcadata)
points(x6G08_MD_142_proj[3],x6G08_MD_142_proj[2],pch=20)
text(x6G08_MD_142_proj[3],x6G08_MD_142_proj[2],pos=4,label="142",col="black")

x6G08_MD_143_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_143.txt")
x6G08_MD_143_proj<-project.pca(x6G08_MD_143_raw,pcadata)
points(x6G08_MD_143_proj[3],x6G08_MD_143_proj[2],pch=20)
text(x6G08_MD_143_proj[3],x6G08_MD_143_proj[2],pos=4,label="143",col="black")

x6G08_MD_144_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_144.txt")
x6G08_MD_144_proj<-project.pca(x6G08_MD_144_raw,pcadata)
points(x6G08_MD_144_proj[3],x6G08_MD_144_proj[2],pch=20)
text(x6G08_MD_144_proj[3],x6G08_MD_144_proj[2],pos=4,label="144",col="blue")

x6G08_MD_145_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_145.txt")
x6G08_MD_145_proj<-project.pca(x6G08_MD_145_raw,pcadata)
points(x6G08_MD_145_proj[3],x6G08_MD_145_proj[2],pch=20)
text(x6G08_MD_145_proj[3],x6G08_MD_145_proj[2],pos=4,label="145",col="black")

x6G08_MD_146_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_146.txt")
x6G08_MD_146_proj<-project.pca(x6G08_MD_146_raw,pcadata)
points(x6G08_MD_146_proj[3],x6G08_MD_146_proj[2],pch=20)
text(x6G08_MD_146_proj[3],x6G08_MD_146_proj[2],pos=4,label="146",col="black")

x6G08_MD_147_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_147.txt")
x6G08_MD_147_proj<-project.pca(x6G08_MD_147_raw,pcadata)
points(x6G08_MD_147_proj[3],x6G08_MD_147_proj[2],pch=20)
text(x6G08_MD_147_proj[3],x6G08_MD_147_proj[2],pos=4,label="147",col="green")

x6G08_MD_148_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_148.txt")
x6G08_MD_148_proj<-project.pca(x6G08_MD_148_raw,pcadata)
points(x6G08_MD_148_proj[3],x6G08_MD_148_proj[2],pch=20)
text(x6G08_MD_148_proj[3],x6G08_MD_148_proj[2],pos=4,label="148",col="green")

x6G08_MD_149_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_149.txt")
x6G08_MD_149_proj<-project.pca(x6G08_MD_149_raw,pcadata)
points(x6G08_MD_149_proj[3],x6G08_MD_149_proj[2],pch=20)
text(x6G08_MD_149_proj[3],x6G08_MD_149_proj[2],pos=4,label="149",col="blue")

x6G08_MD_150_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_150.txt")
x6G08_MD_150_proj<-project.pca(x6G08_MD_150_raw,pcadata)
points(x6G08_MD_150_proj[3],x6G08_MD_150_proj[2],pch=20)
text(x6G08_MD_150_proj[3],x6G08_MD_150_proj[2],pos=4,label="150",col="red")

x6G08_MD_151_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_151.txt")
x6G08_MD_151_proj<-project.pca(x6G08_MD_151_raw,pcadata)
points(x6G08_MD_151_proj[3],x6G08_MD_151_proj[2],pch=20)
text(x6G08_MD_151_proj[3],x6G08_MD_151_proj[2],pos=4,label="151",col="black")

x6G08_MD_152_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_152.txt")
x6G08_MD_152_proj<-project.pca(x6G08_MD_152_raw,pcadata)
points(x6G08_MD_152_proj[3],x6G08_MD_152_proj[2],pch=20)
text(x6G08_MD_152_proj[3],x6G08_MD_152_proj[2],pos=4,label="152",col="black")

x6G08_MD_153_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_153.txt")
x6G08_MD_153_proj<-project.pca(x6G08_MD_153_raw,pcadata)
points(x6G08_MD_153_proj[3],x6G08_MD_153_proj[2],pch=20)
text(x6G08_MD_153_proj[3],x6G08_MD_153_proj[2],pos=4,label="153",col="black")

x6G08_MD_154_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_154.txt")
x6G08_MD_154_proj<-project.pca(x6G08_MD_154_raw,pcadata)
points(x6G08_MD_154_proj[3],x6G08_MD_154_proj[2],pch=20)
text(x6G08_MD_154_proj[3],x6G08_MD_154_proj[2],pos=4,label="154",col="black")

x6G08_MD_155_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_155.txt")
x6G08_MD_155_proj<-project.pca(x6G08_MD_155_raw,pcadata)
points(x6G08_MD_155_proj[3],x6G08_MD_155_proj[2],pch=20)
text(x6G08_MD_155_proj[3],x6G08_MD_155_proj[2],pos=4,label="155",col="black")

x6G08_MD_156_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_156.txt")
x6G08_MD_156_proj<-project.pca(x6G08_MD_156_raw,pcadata)
points(x6G08_MD_156_proj[3],x6G08_MD_156_proj[2],pch=20)
text(x6G08_MD_156_proj[3],x6G08_MD_156_proj[2],pos=4,label="156",col="black")

x6G08_MD_157_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_157.txt")
x6G08_MD_157_proj<-project.pca(x6G08_MD_157_raw,pcadata)
points(x6G08_MD_157_proj[3],x6G08_MD_157_proj[2],pch=20)
text(x6G08_MD_157_proj[3],x6G08_MD_157_proj[2],pos=4,label="157",col="black")

x6G08_MD_158_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_158.txt")
x6G08_MD_158_proj<-project.pca(x6G08_MD_158_raw,pcadata)
points(x6G08_MD_158_proj[3],x6G08_MD_158_proj[2],pch=20)
text(x6G08_MD_158_proj[3],x6G08_MD_158_proj[2],pos=4,label="158",col="black")

x6G08_MD_159_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_159.txt")
x6G08_MD_159_proj<-project.pca(x6G08_MD_159_raw,pcadata)
points(x6G08_MD_159_proj[3],x6G08_MD_159_proj[2],pch=20)
text(x6G08_MD_159_proj[3],x6G08_MD_159_proj[2],pos=4,label="159",col="black")

x6G08_MD_160_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_160.txt")
x6G08_MD_160_proj<-project.pca(x6G08_MD_160_raw,pcadata)
points(x6G08_MD_160_proj[3],x6G08_MD_160_proj[2],pch=20)
text(x6G08_MD_160_proj[3],x6G08_MD_160_proj[2],pos=4,label="160",col="black")

x6G08_MD_161_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_161.txt")
x6G08_MD_161_proj<-project.pca(x6G08_MD_161_raw,pcadata)
points(x6G08_MD_161_proj[3],x6G08_MD_161_proj[2],pch=20)
text(x6G08_MD_161_proj[3],x6G08_MD_161_proj[2],pos=4,label="161",col="black")

x6G08_MD_162_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_162.txt")
x6G08_MD_162_proj<-project.pca(x6G08_MD_162_raw,pcadata)
points(x6G08_MD_162_proj[3],x6G08_MD_162_proj[2],pch=20)
text(x6G08_MD_162_proj[3],x6G08_MD_162_proj[2],pos=4,label="162",col="black")

x6G08_MD_163_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_163.txt")
x6G08_MD_163_proj<-project.pca(x6G08_MD_163_raw,pcadata)
points(x6G08_MD_163_proj[3],x6G08_MD_163_proj[2],pch=20)
text(x6G08_MD_163_proj[3],x6G08_MD_163_proj[2],pos=4,label="163",col="black")

x6G08_MD_164_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_164.txt")
x6G08_MD_164_proj<-project.pca(x6G08_MD_164_raw,pcadata)
points(x6G08_MD_164_proj[3],x6G08_MD_164_proj[2],pch=20)
text(x6G08_MD_164_proj[3],x6G08_MD_164_proj[2],pos=4,label="164",col="black")

x6G08_MD_165_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_165.txt")
x6G08_MD_165_proj<-project.pca(x6G08_MD_165_raw,pcadata)
points(x6G08_MD_165_proj[3],x6G08_MD_165_proj[2],pch=20)
text(x6G08_MD_165_proj[3],x6G08_MD_165_proj[2],pos=4,label="165",col="green")

x6G08_MD_166_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_166.txt")
x6G08_MD_166_proj<-project.pca(x6G08_MD_166_raw,pcadata)
points(x6G08_MD_166_proj[3],x6G08_MD_166_proj[2],pch=20)
text(x6G08_MD_166_proj[3],x6G08_MD_166_proj[2],pos=4,label="166",col="black")

x6G08_MD_167_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_167.txt")
x6G08_MD_167_proj<-project.pca(x6G08_MD_167_raw,pcadata)
points(x6G08_MD_167_proj[3],x6G08_MD_167_proj[2],pch=20)
text(x6G08_MD_167_proj[3],x6G08_MD_167_proj[2],pos=4,label="167",col="black")

x6G08_MD_168_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_168.txt")
x6G08_MD_168_proj<-project.pca(x6G08_MD_168_raw,pcadata)
points(x6G08_MD_168_proj[3],x6G08_MD_168_proj[2],pch=20)
text(x6G08_MD_168_proj[3],x6G08_MD_168_proj[2],pos=4,label="168",col="black")

x6G08_MD_169_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_169.txt")
x6G08_MD_169_proj<-project.pca(x6G08_MD_169_raw,pcadata)
points(x6G08_MD_169_proj[3],x6G08_MD_169_proj[2],pch=20)
text(x6G08_MD_169_proj[3],x6G08_MD_169_proj[2],pos=4,label="169",col="black")

x6G08_MD_170_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_170.txt")
x6G08_MD_170_proj<-project.pca(x6G08_MD_170_raw,pcadata)
points(x6G08_MD_170_proj[3],x6G08_MD_170_proj[2],pch=20)
text(x6G08_MD_170_proj[3],x6G08_MD_170_proj[2],pos=4,label="170",col="red")

x6G08_MD_171_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_171.txt")
x6G08_MD_171_proj<-project.pca(x6G08_MD_171_raw,pcadata)
points(x6G08_MD_171_proj[3],x6G08_MD_171_proj[2],pch=20)
text(x6G08_MD_171_proj[3],x6G08_MD_171_proj[2],pos=4,label="171",col="red")

x6G08_MD_172_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_172.txt")
x6G08_MD_172_proj<-project.pca(x6G08_MD_172_raw,pcadata)
points(x6G08_MD_172_proj[3],x6G08_MD_172_proj[2],pch=20)
text(x6G08_MD_172_proj[3],x6G08_MD_172_proj[2],pos=4,label="172",col="blue")

x6G08_MD_173_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_173.txt")
x6G08_MD_173_proj<-project.pca(x6G08_MD_173_raw,pcadata)
points(x6G08_MD_173_proj[3],x6G08_MD_173_proj[2],pch=20)
text(x6G08_MD_173_proj[3],x6G08_MD_173_proj[2],pos=4,label="173",col="black")

x6G08_MD_174_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_174.txt")
x6G08_MD_174_proj<-project.pca(x6G08_MD_174_raw,pcadata)
points(x6G08_MD_174_proj[3],x6G08_MD_174_proj[2],pch=20)
text(x6G08_MD_174_proj[3],x6G08_MD_174_proj[2],pos=4,label="174",col="blue")

x6G08_MD_175_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_175.txt")
x6G08_MD_175_proj<-project.pca(x6G08_MD_175_raw,pcadata)
points(x6G08_MD_175_proj[3],x6G08_MD_175_proj[2],pch=20)
text(x6G08_MD_175_proj[3],x6G08_MD_175_proj[2],pos=4,label="175",col="blue")

x6G08_MD_176_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_176.txt")
x6G08_MD_176_proj<-project.pca(x6G08_MD_176_raw,pcadata)
points(x6G08_MD_176_proj[3],x6G08_MD_176_proj[2],pch=20)
text(x6G08_MD_176_proj[3],x6G08_MD_176_proj[2],pos=4,label="176",col="blue")

x6G08_MD_177_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_177.txt")
x6G08_MD_177_proj<-project.pca(x6G08_MD_177_raw,pcadata)
points(x6G08_MD_177_proj[3],x6G08_MD_177_proj[2],pch=20)
text(x6G08_MD_177_proj[3],x6G08_MD_177_proj[2],pos=4,label="177",col="blue")

x6G08_MD_178_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_178.txt")
x6G08_MD_178_proj<-project.pca(x6G08_MD_178_raw,pcadata)
points(x6G08_MD_178_proj[3],x6G08_MD_178_proj[2],pch=20)
text(x6G08_MD_178_proj[3],x6G08_MD_178_proj[2],pos=4,label="178",col="blue")

x6G08_MD_179_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_179.txt")
x6G08_MD_179_proj<-project.pca(x6G08_MD_179_raw,pcadata)
points(x6G08_MD_179_proj[3],x6G08_MD_179_proj[2],pch=20)
text(x6G08_MD_179_proj[3],x6G08_MD_179_proj[2],pos=4,label="179",col="blue")

x6G08_MD_180_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_180.txt")
x6G08_MD_180_proj<-project.pca(x6G08_MD_180_raw,pcadata)
points(x6G08_MD_180_proj[3],x6G08_MD_180_proj[2],pch=20)
text(x6G08_MD_180_proj[3],x6G08_MD_180_proj[2],pos=4,label="180",col="black")

x6G08_MD_181_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_181.txt")
x6G08_MD_181_proj<-project.pca(x6G08_MD_181_raw,pcadata)
points(x6G08_MD_181_proj[3],x6G08_MD_181_proj[2],pch=20)
text(x6G08_MD_181_proj[3],x6G08_MD_181_proj[2],pos=4,label="181",col="black")

x6G08_MD_182_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_182.txt")
x6G08_MD_182_proj<-project.pca(x6G08_MD_182_raw,pcadata)
points(x6G08_MD_182_proj[3],x6G08_MD_182_proj[2],pch=20)
text(x6G08_MD_182_proj[3],x6G08_MD_182_proj[2],pos=4,label="182",col="black")

x6G08_MD_183_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_183.txt")
x6G08_MD_183_proj<-project.pca(x6G08_MD_183_raw,pcadata)
points(x6G08_MD_183_proj[3],x6G08_MD_183_proj[2],pch=20)
text(x6G08_MD_183_proj[3],x6G08_MD_183_proj[2],pos=4,label="183",col="black")

x6G08_MD_184_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_184.txt")
x6G08_MD_184_proj<-project.pca(x6G08_MD_184_raw,pcadata)
points(x6G08_MD_184_proj[3],x6G08_MD_184_proj[2],pch=20)
text(x6G08_MD_184_proj[3],x6G08_MD_184_proj[2],pos=4,label="184",col="green")

x6G08_MD_185_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_185.txt")
x6G08_MD_185_proj<-project.pca(x6G08_MD_185_raw,pcadata)
points(x6G08_MD_185_proj[3],x6G08_MD_185_proj[2],pch=20)
text(x6G08_MD_185_proj[3],x6G08_MD_185_proj[2],pos=4,label="185",col="blue")

x6G08_MD_186_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_186.txt")
x6G08_MD_186_proj<-project.pca(x6G08_MD_186_raw,pcadata)
points(x6G08_MD_186_proj[3],x6G08_MD_186_proj[2],pch=20)
text(x6G08_MD_186_proj[3],x6G08_MD_186_proj[2],pos=4,label="186",col="red")

x6G08_MD_187_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_187.txt")
x6G08_MD_187_proj<-project.pca(x6G08_MD_187_raw,pcadata)
points(x6G08_MD_187_proj[3],x6G08_MD_187_proj[2],pch=20)
text(x6G08_MD_187_proj[3],x6G08_MD_187_proj[2],pos=4,label="187",col="blue")

x6G08_MD_188_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_188.txt")
x6G08_MD_188_proj<-project.pca(x6G08_MD_188_raw,pcadata)
points(x6G08_MD_188_proj[3],x6G08_MD_188_proj[2],pch=20)
text(x6G08_MD_188_proj[3],x6G08_MD_188_proj[2],pos=4,label="188",col="blue")

x6G08_MD_189_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_189.txt")
x6G08_MD_189_proj<-project.pca(x6G08_MD_189_raw,pcadata)
points(x6G08_MD_189_proj[3],x6G08_MD_189_proj[2],pch=20)
text(x6G08_MD_189_proj[3],x6G08_MD_189_proj[2],pos=4,label="189",col="blue")

x6G08_MD_190_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_190.txt")
x6G08_MD_190_proj<-project.pca(x6G08_MD_190_raw,pcadata)
points(x6G08_MD_190_proj[3],x6G08_MD_190_proj[2],pch=20)
text(x6G08_MD_190_proj[3],x6G08_MD_190_proj[2],pos=4,label="190",col="blue")

x6G08_MD_191_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_191.txt")
x6G08_MD_191_proj<-project.pca(x6G08_MD_191_raw,pcadata)
points(x6G08_MD_191_proj[3],x6G08_MD_191_proj[2],pch=20)
text(x6G08_MD_191_proj[3],x6G08_MD_191_proj[2],pos=4,label="191",col="black")

x6G08_MD_192_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_192.txt")
x6G08_MD_192_proj<-project.pca(x6G08_MD_192_raw,pcadata)
points(x6G08_MD_192_proj[3],x6G08_MD_192_proj[2],pch=20)
text(x6G08_MD_192_proj[3],x6G08_MD_192_proj[2],pos=4,label="192",col="black")

x6G08_MD_193_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_193.txt")
x6G08_MD_193_proj<-project.pca(x6G08_MD_193_raw,pcadata)
points(x6G08_MD_193_proj[3],x6G08_MD_193_proj[2],pch=20)
text(x6G08_MD_193_proj[3],x6G08_MD_193_proj[2],pos=4,label="193",col="blue")

x6G08_MD_194_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_194.txt")
x6G08_MD_194_proj<-project.pca(x6G08_MD_194_raw,pcadata)
points(x6G08_MD_194_proj[3],x6G08_MD_194_proj[2],pch=20)
text(x6G08_MD_194_proj[3],x6G08_MD_194_proj[2],pos=4,label="194",col="blue")

x6G08_MD_195_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_195.txt")
x6G08_MD_195_proj<-project.pca(x6G08_MD_195_raw,pcadata)
points(x6G08_MD_195_proj[3],x6G08_MD_195_proj[2],pch=20)
text(x6G08_MD_195_proj[3],x6G08_MD_195_proj[2],pos=4,label="195",col="blue")

x6G08_MD_196_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_196.txt")
x6G08_MD_196_proj<-project.pca(x6G08_MD_196_raw,pcadata)
points(x6G08_MD_196_proj[3],x6G08_MD_196_proj[2],pch=20)
text(x6G08_MD_196_proj[3],x6G08_MD_196_proj[2],pos=4,label="196",col="black")

x6G08_MD_197_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_197.txt")
x6G08_MD_197_proj<-project.pca(x6G08_MD_197_raw,pcadata)
points(x6G08_MD_197_proj[3],x6G08_MD_197_proj[2],pch=20)
text(x6G08_MD_197_proj[3],x6G08_MD_197_proj[2],pos=4,label="197",col="black")

x6G08_MD_198_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_198.txt")
x6G08_MD_198_proj<-project.pca(x6G08_MD_198_raw,pcadata)
points(x6G08_MD_198_proj[3],x6G08_MD_198_proj[2],pch=20)
text(x6G08_MD_198_proj[3],x6G08_MD_198_proj[2],pos=4,label="198",col="black")

x6G08_MD_199_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_199.txt")
x6G08_MD_199_proj<-project.pca(x6G08_MD_199_raw,pcadata)
points(x6G08_MD_199_proj[3],x6G08_MD_199_proj[2],pch=20)
text(x6G08_MD_199_proj[3],x6G08_MD_199_proj[2],pos=4,label="199",col="black")

x6G08_MD_200_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_200.txt")
x6G08_MD_200_proj<-project.pca(x6G08_MD_200_raw,pcadata)
points(x6G08_MD_200_proj[3],x6G08_MD_200_proj[2],pch=20)
text(x6G08_MD_200_proj[3],x6G08_MD_200_proj[2],pos=4,label="200",col="black")

x6G08_MD_201_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_201.txt")
x6G08_MD_201_proj<-project.pca(x6G08_MD_201_raw,pcadata)
points(x6G08_MD_201_proj[3],x6G08_MD_201_proj[2],pch=20)
text(x6G08_MD_201_proj[3],x6G08_MD_201_proj[2],pos=4,label="201",col="blue")

x6G08_MD_202_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_202.txt")
x6G08_MD_202_proj<-project.pca(x6G08_MD_202_raw,pcadata)
points(x6G08_MD_202_proj[3],x6G08_MD_202_proj[2],pch=20)
text(x6G08_MD_202_proj[3],x6G08_MD_202_proj[2],pos=4,label="202",col="black")

x6G08_MD_203_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_203.txt")
x6G08_MD_203_proj<-project.pca(x6G08_MD_203_raw,pcadata)
points(x6G08_MD_203_proj[3],x6G08_MD_203_proj[2],pch=20)
text(x6G08_MD_203_proj[3],x6G08_MD_203_proj[2],pos=4,label="203",col="blue")

x6G08_MD_204_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_204.txt")
x6G08_MD_204_proj<-project.pca(x6G08_MD_204_raw,pcadata)
points(x6G08_MD_204_proj[3],x6G08_MD_204_proj[2],pch=20)
text(x6G08_MD_204_proj[3],x6G08_MD_204_proj[2],pos=4,label="204",col="blue")

x6G08_MD_205_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_205.txt")
x6G08_MD_205_proj<-project.pca(x6G08_MD_205_raw,pcadata)
points(x6G08_MD_205_proj[3],x6G08_MD_205_proj[2],pch=20)
text(x6G08_MD_205_proj[3],x6G08_MD_205_proj[2],pos=4,label="205",col="blue")

x6G08_MD_206_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_206.txt")
x6G08_MD_206_proj<-project.pca(x6G08_MD_206_raw,pcadata)
points(x6G08_MD_206_proj[3],x6G08_MD_206_proj[2],pch=20)
text(x6G08_MD_206_proj[3],x6G08_MD_206_proj[2],pos=4,label="206",col="red")

x6G08_MD_207_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_207.txt")
x6G08_MD_207_proj<-project.pca(x6G08_MD_207_raw,pcadata)
points(x6G08_MD_207_proj[3],x6G08_MD_207_proj[2],pch=20)
text(x6G08_MD_207_proj[3],x6G08_MD_207_proj[2],pos=4,label="207",col="blue")

x6G08_MD_208_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_208.txt")
x6G08_MD_208_proj<-project.pca(x6G08_MD_208_raw,pcadata)
points(x6G08_MD_208_proj[3],x6G08_MD_208_proj[2],pch=20)
text(x6G08_MD_208_proj[3],x6G08_MD_208_proj[2],pos=4,label="208",col="blue")

x6G08_MD_209_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_209.txt")
x6G08_MD_209_proj<-project.pca(x6G08_MD_209_raw,pcadata)
points(x6G08_MD_209_proj[3],x6G08_MD_209_proj[2],pch=20)
text(x6G08_MD_209_proj[3],x6G08_MD_209_proj[2],pos=4,label="209",col="black")

x6G08_MD_210_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_210.txt")
x6G08_MD_210_proj<-project.pca(x6G08_MD_210_raw,pcadata)
points(x6G08_MD_210_proj[3],x6G08_MD_210_proj[2],pch=20)
text(x6G08_MD_210_proj[3],x6G08_MD_210_proj[2],pos=4,label="210",col="blue")

x6G08_MD_211_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_211.txt")
x6G08_MD_211_proj<-project.pca(x6G08_MD_211_raw,pcadata)
points(x6G08_MD_211_proj[3],x6G08_MD_211_proj[2],pch=20)
text(x6G08_MD_211_proj[3],x6G08_MD_211_proj[2],pos=4,label="211",col="blue")

x6G08_MD_212_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_212.txt")
x6G08_MD_212_proj<-project.pca(x6G08_MD_212_raw,pcadata)
points(x6G08_MD_212_proj[3],x6G08_MD_212_proj[2],pch=20)
text(x6G08_MD_212_proj[3],x6G08_MD_212_proj[2],pos=4,label="212",col="blue")

x6G08_MD_213_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_213.txt")
x6G08_MD_213_proj<-project.pca(x6G08_MD_213_raw,pcadata)
points(x6G08_MD_213_proj[3],x6G08_MD_213_proj[2],pch=20)
text(x6G08_MD_213_proj[3],x6G08_MD_213_proj[2],pos=4,label="213",col="blue")

x6G08_MD_214_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_214.txt")
x6G08_MD_214_proj<-project.pca(x6G08_MD_214_raw,pcadata)
points(x6G08_MD_214_proj[3],x6G08_MD_214_proj[2],pch=20)
text(x6G08_MD_214_proj[3],x6G08_MD_214_proj[2],pos=4,label="214",col="blue")

x6G08_MD_215_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_215.txt")
x6G08_MD_215_proj<-project.pca(x6G08_MD_215_raw,pcadata)
points(x6G08_MD_215_proj[3],x6G08_MD_215_proj[2],pch=20)
text(x6G08_MD_215_proj[3],x6G08_MD_215_proj[2],pos=4,label="215",col="blue")

x6G08_MD_216_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_216.txt")
x6G08_MD_216_proj<-project.pca(x6G08_MD_216_raw,pcadata)
points(x6G08_MD_216_proj[3],x6G08_MD_216_proj[2],pch=20)
text(x6G08_MD_216_proj[3],x6G08_MD_216_proj[2],pos=4,label="216",col="blue")

x6G08_MD_217_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_217.txt")
x6G08_MD_217_proj<-project.pca(x6G08_MD_217_raw,pcadata)
points(x6G08_MD_217_proj[3],x6G08_MD_217_proj[2],pch=20)
text(x6G08_MD_217_proj[3],x6G08_MD_217_proj[2],pos=4,label="217",col="blue")

x6G08_MD_218_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_218.txt")
x6G08_MD_218_proj<-project.pca(x6G08_MD_218_raw,pcadata)
points(x6G08_MD_218_proj[3],x6G08_MD_218_proj[2],pch=20)
text(x6G08_MD_218_proj[3],x6G08_MD_218_proj[2],pos=4,label="218",col="green")

x6G08_MD_219_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_219.txt")
x6G08_MD_219_proj<-project.pca(x6G08_MD_219_raw,pcadata)
points(x6G08_MD_219_proj[3],x6G08_MD_219_proj[2],pch=20)
text(x6G08_MD_219_proj[3],x6G08_MD_219_proj[2],pos=4,label="219",col="red")

x6G08_MD_220_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_220.txt")
x6G08_MD_220_proj<-project.pca(x6G08_MD_220_raw,pcadata)
points(x6G08_MD_220_proj[3],x6G08_MD_220_proj[2],pch=20)
text(x6G08_MD_220_proj[3],x6G08_MD_220_proj[2],pos=4,label="220",col="red")

x6G08_MD_221_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_221.txt")
x6G08_MD_221_proj<-project.pca(x6G08_MD_221_raw,pcadata)
points(x6G08_MD_221_proj[3],x6G08_MD_221_proj[2],pch=20)
text(x6G08_MD_221_proj[3],x6G08_MD_221_proj[2],pos=4,label="221",col="blue")

x6G08_MD_222_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_222.txt")
x6G08_MD_222_proj<-project.pca(x6G08_MD_222_raw,pcadata)
points(x6G08_MD_222_proj[3],x6G08_MD_222_proj[2],pch=20)
text(x6G08_MD_222_proj[3],x6G08_MD_222_proj[2],pos=4,label="222",col="green")

x6G08_MD_223_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_223.txt")
x6G08_MD_223_proj<-project.pca(x6G08_MD_223_raw,pcadata)
points(x6G08_MD_223_proj[3],x6G08_MD_223_proj[2],pch=20)
text(x6G08_MD_223_proj[3],x6G08_MD_223_proj[2],pos=4,label="223",col="red")

x6G08_MD_224_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_224.txt")
x6G08_MD_224_proj<-project.pca(x6G08_MD_224_raw,pcadata)
points(x6G08_MD_224_proj[3],x6G08_MD_224_proj[2],pch=20)
text(x6G08_MD_224_proj[3],x6G08_MD_224_proj[2],pos=4,label="224",col="red")

x6G08_MD_225_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_225.txt")
x6G08_MD_225_proj<-project.pca(x6G08_MD_225_raw,pcadata)
points(x6G08_MD_225_proj[3],x6G08_MD_225_proj[2],pch=20)
text(x6G08_MD_225_proj[3],x6G08_MD_225_proj[2],pos=4,label="225",col="red")

x6G08_MD_226_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_226.txt")
x6G08_MD_226_proj<-project.pca(x6G08_MD_226_raw,pcadata)
points(x6G08_MD_226_proj[3],x6G08_MD_226_proj[2],pch=20)
text(x6G08_MD_226_proj[3],x6G08_MD_226_proj[2],pos=4,label="226",col="blue")

x6G08_MD_227_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_227.txt")
x6G08_MD_227_proj<-project.pca(x6G08_MD_227_raw,pcadata)
points(x6G08_MD_227_proj[3],x6G08_MD_227_proj[2],pch=20)
text(x6G08_MD_227_proj[3],x6G08_MD_227_proj[2],pos=4,label="227",col="red")

x6G08_MD_228_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_228.txt")
x6G08_MD_228_proj<-project.pca(x6G08_MD_228_raw,pcadata)
points(x6G08_MD_228_proj[3],x6G08_MD_228_proj[2],pch=20)
text(x6G08_MD_228_proj[3],x6G08_MD_228_proj[2],pos=4,label="228",col="green")

x6G08_MD_229_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_229.txt")
x6G08_MD_229_proj<-project.pca(x6G08_MD_229_raw,pcadata)
points(x6G08_MD_229_proj[3],x6G08_MD_229_proj[2],pch=20)
text(x6G08_MD_229_proj[3],x6G08_MD_229_proj[2],pos=4,label="229",col="red")

x6G08_MD_230_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_230.txt")
x6G08_MD_230_proj<-project.pca(x6G08_MD_230_raw,pcadata)
points(x6G08_MD_230_proj[3],x6G08_MD_230_proj[2],pch=20)
text(x6G08_MD_230_proj[3],x6G08_MD_230_proj[2],pos=4,label="230",col="green")

x6G08_MD_231_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_231.txt")
x6G08_MD_231_proj<-project.pca(x6G08_MD_231_raw,pcadata)
points(x6G08_MD_231_proj[3],x6G08_MD_231_proj[2],pch=20)
text(x6G08_MD_231_proj[3],x6G08_MD_231_proj[2],pos=4,label="231",col="green")

x6G08_MD_232_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_232.txt")
x6G08_MD_232_proj<-project.pca(x6G08_MD_232_raw,pcadata)
points(x6G08_MD_232_proj[3],x6G08_MD_232_proj[2],pch=20)
text(x6G08_MD_232_proj[3],x6G08_MD_232_proj[2],pos=4,label="232",col="red")

x6G08_MD_233_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_233.txt")
x6G08_MD_233_proj<-project.pca(x6G08_MD_233_raw,pcadata)
points(x6G08_MD_233_proj[3],x6G08_MD_233_proj[2],pch=20)
text(x6G08_MD_233_proj[3],x6G08_MD_233_proj[2],pos=4,label="233",col="blue")

x6G08_MD_234_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_234.txt")
x6G08_MD_234_proj<-project.pca(x6G08_MD_234_raw,pcadata)
points(x6G08_MD_234_proj[3],x6G08_MD_234_proj[2],pch=20)
text(x6G08_MD_234_proj[3],x6G08_MD_234_proj[2],pos=4,label="234",col="green")

x6G08_MD_235_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_235.txt")
x6G08_MD_235_proj<-project.pca(x6G08_MD_235_raw,pcadata)
points(x6G08_MD_235_proj[3],x6G08_MD_235_proj[2],pch=20)
text(x6G08_MD_235_proj[3],x6G08_MD_235_proj[2],pos=4,label="235",col="green")

x6G08_MD_236_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_236.txt")
x6G08_MD_236_proj<-project.pca(x6G08_MD_236_raw,pcadata)
points(x6G08_MD_236_proj[3],x6G08_MD_236_proj[2],pch=20)
text(x6G08_MD_236_proj[3],x6G08_MD_236_proj[2],pos=4,label="236",col="blue")

x6G08_MD_237_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_237.txt")
x6G08_MD_237_proj<-project.pca(x6G08_MD_237_raw,pcadata)
points(x6G08_MD_237_proj[3],x6G08_MD_237_proj[2],pch=20)
text(x6G08_MD_237_proj[3],x6G08_MD_237_proj[2],pos=4,label="237",col="green")

x6G08_MD_238_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_238.txt")
x6G08_MD_238_proj<-project.pca(x6G08_MD_238_raw,pcadata)
points(x6G08_MD_238_proj[3],x6G08_MD_238_proj[2],pch=20)
text(x6G08_MD_238_proj[3],x6G08_MD_238_proj[2],pos=4,label="238",col="green")

x6G08_MD_239_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_239.txt")
x6G08_MD_239_proj<-project.pca(x6G08_MD_239_raw,pcadata)
points(x6G08_MD_239_proj[3],x6G08_MD_239_proj[2],pch=20)
text(x6G08_MD_239_proj[3],x6G08_MD_239_proj[2],pos=4,label="239",col="red")

x6G08_MD_240_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_240.txt")
x6G08_MD_240_proj<-project.pca(x6G08_MD_240_raw,pcadata)
points(x6G08_MD_240_proj[3],x6G08_MD_240_proj[2],pch=20)
text(x6G08_MD_240_proj[3],x6G08_MD_240_proj[2],pos=4,label="240",col="green")

x6G08_MD_241_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_241.txt")
x6G08_MD_241_proj<-project.pca(x6G08_MD_241_raw,pcadata)
points(x6G08_MD_241_proj[3],x6G08_MD_241_proj[2],pch=20)
text(x6G08_MD_241_proj[3],x6G08_MD_241_proj[2],pos=4,label="241",col="green")

x6G08_MD_242_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_242.txt")
x6G08_MD_242_proj<-project.pca(x6G08_MD_242_raw,pcadata)
points(x6G08_MD_242_proj[3],x6G08_MD_242_proj[2],pch=20)
text(x6G08_MD_242_proj[3],x6G08_MD_242_proj[2],pos=4,label="242",col="green")

x6G08_MD_243_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_243.txt")
x6G08_MD_243_proj<-project.pca(x6G08_MD_243_raw,pcadata)
points(x6G08_MD_243_proj[3],x6G08_MD_243_proj[2],pch=20)
text(x6G08_MD_243_proj[3],x6G08_MD_243_proj[2],pos=4,label="243",col="green")

x6G08_MD_244_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_244.txt")
x6G08_MD_244_proj<-project.pca(x6G08_MD_244_raw,pcadata)
points(x6G08_MD_244_proj[3],x6G08_MD_244_proj[2],pch=20)
text(x6G08_MD_244_proj[3],x6G08_MD_244_proj[2],pos=4,label="244",col="blue")

x6G08_MD_245_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_245.txt")
x6G08_MD_245_proj<-project.pca(x6G08_MD_245_raw,pcadata)
points(x6G08_MD_245_proj[3],x6G08_MD_245_proj[2],pch=20)
text(x6G08_MD_245_proj[3],x6G08_MD_245_proj[2],pos=4,label="245",col="green")

x6G08_MD_246_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_246.txt")
x6G08_MD_246_proj<-project.pca(x6G08_MD_246_raw,pcadata)
points(x6G08_MD_246_proj[3],x6G08_MD_246_proj[2],pch=20)
text(x6G08_MD_246_proj[3],x6G08_MD_246_proj[2],pos=4,label="246",col="green")

x6G08_MD_247_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_247.txt")
x6G08_MD_247_proj<-project.pca(x6G08_MD_247_raw,pcadata)
points(x6G08_MD_247_proj[3],x6G08_MD_247_proj[2],pch=20)
text(x6G08_MD_247_proj[3],x6G08_MD_247_proj[2],pos=4,label="247",col="blue")

x6G08_MD_248_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_248.txt")
x6G08_MD_248_proj<-project.pca(x6G08_MD_248_raw,pcadata)
points(x6G08_MD_248_proj[3],x6G08_MD_248_proj[2],pch=20)
text(x6G08_MD_248_proj[3],x6G08_MD_248_proj[2],pos=4,label="248",col="black")

x6G08_MD_249_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_249.txt")
x6G08_MD_249_proj<-project.pca(x6G08_MD_249_raw,pcadata)
points(x6G08_MD_249_proj[3],x6G08_MD_249_proj[2],pch=20)
text(x6G08_MD_249_proj[3],x6G08_MD_249_proj[2],pos=4,label="249",col="blue")

x6G08_MD_250_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_250.txt")
x6G08_MD_250_proj<-project.pca(x6G08_MD_250_raw,pcadata)
points(x6G08_MD_250_proj[3],x6G08_MD_250_proj[2],pch=20)
text(x6G08_MD_250_proj[3],x6G08_MD_250_proj[2],pos=4,label="250",col="black")

x6G08_MD_251_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_251.txt")
x6G08_MD_251_proj<-project.pca(x6G08_MD_251_raw,pcadata)
points(x6G08_MD_251_proj[3],x6G08_MD_251_proj[2],pch=20)
text(x6G08_MD_251_proj[3],x6G08_MD_251_proj[2],pos=4,label="251",col="black")

x6G08_MD_252_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_252.txt")
x6G08_MD_252_proj<-project.pca(x6G08_MD_252_raw,pcadata)
points(x6G08_MD_252_proj[3],x6G08_MD_252_proj[2],pch=20)
text(x6G08_MD_252_proj[3],x6G08_MD_252_proj[2],pos=4,label="252",col="blue")

x6G08_MD_253_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_253.txt")
x6G08_MD_253_proj<-project.pca(x6G08_MD_253_raw,pcadata)
points(x6G08_MD_253_proj[3],x6G08_MD_253_proj[2],pch=20)
text(x6G08_MD_253_proj[3],x6G08_MD_253_proj[2],pos=4,label="253",col="blue")

x6G08_MD_254_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_254.txt")
x6G08_MD_254_proj<-project.pca(x6G08_MD_254_raw,pcadata)
points(x6G08_MD_254_proj[3],x6G08_MD_254_proj[2],pch=20)
text(x6G08_MD_254_proj[3],x6G08_MD_254_proj[2],pos=4,label="254",col="black")

x6G08_MD_255_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_255.txt")
x6G08_MD_255_proj<-project.pca(x6G08_MD_255_raw,pcadata)
points(x6G08_MD_255_proj[3],x6G08_MD_255_proj[2],pch=20)
text(x6G08_MD_255_proj[3],x6G08_MD_255_proj[2],pos=4,label="255",col="blue")

x6G08_MD_256_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_256.txt")
x6G08_MD_256_proj<-project.pca(x6G08_MD_256_raw,pcadata)
points(x6G08_MD_256_proj[3],x6G08_MD_256_proj[2],pch=20)
text(x6G08_MD_256_proj[3],x6G08_MD_256_proj[2],pos=4,label="256",col="blue")

x6G08_MD_257_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_257.txt")
x6G08_MD_257_proj<-project.pca(x6G08_MD_257_raw,pcadata)
points(x6G08_MD_257_proj[3],x6G08_MD_257_proj[2],pch=20)
text(x6G08_MD_257_proj[3],x6G08_MD_257_proj[2],pos=4,label="257",col="green")

x6G08_MD_258_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_258.txt")
x6G08_MD_258_proj<-project.pca(x6G08_MD_258_raw,pcadata)
points(x6G08_MD_258_proj[3],x6G08_MD_258_proj[2],pch=20)
text(x6G08_MD_258_proj[3],x6G08_MD_258_proj[2],pos=4,label="258",col="green")

x6G08_MD_259_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_259.txt")
x6G08_MD_259_proj<-project.pca(x6G08_MD_259_raw,pcadata)
points(x6G08_MD_259_proj[3],x6G08_MD_259_proj[2],pch=20)
text(x6G08_MD_259_proj[3],x6G08_MD_259_proj[2],pos=4,label="259",col="green")

x6G08_MD_260_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_260.txt")
x6G08_MD_260_proj<-project.pca(x6G08_MD_260_raw,pcadata)
points(x6G08_MD_260_proj[3],x6G08_MD_260_proj[2],pch=20)
text(x6G08_MD_260_proj[3],x6G08_MD_260_proj[2],pos=4,label="260",col="blue")

x6G08_MD_261_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_261.txt")
x6G08_MD_261_proj<-project.pca(x6G08_MD_261_raw,pcadata)
points(x6G08_MD_261_proj[3],x6G08_MD_261_proj[2],pch=20)
text(x6G08_MD_261_proj[3],x6G08_MD_261_proj[2],pos=4,label="261",col="blue")

x6G08_MD_262_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_262.txt")
x6G08_MD_262_proj<-project.pca(x6G08_MD_262_raw,pcadata)
points(x6G08_MD_262_proj[3],x6G08_MD_262_proj[2],pch=20)
text(x6G08_MD_262_proj[3],x6G08_MD_262_proj[2],pos=4,label="262",col="blue")

x6G08_MD_263_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_263.txt")
x6G08_MD_263_proj<-project.pca(x6G08_MD_263_raw,pcadata)
points(x6G08_MD_263_proj[3],x6G08_MD_263_proj[2],pch=20)
text(x6G08_MD_263_proj[3],x6G08_MD_263_proj[2],pos=4,label="263",col="blue")

x6G08_MD_264_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_264.txt")
x6G08_MD_264_proj<-project.pca(x6G08_MD_264_raw,pcadata)
points(x6G08_MD_264_proj[3],x6G08_MD_264_proj[2],pch=20)
text(x6G08_MD_264_proj[3],x6G08_MD_264_proj[2],pos=4,label="264",col="blue")

x6G08_MD_265_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_265.txt")
x6G08_MD_265_proj<-project.pca(x6G08_MD_265_raw,pcadata)
points(x6G08_MD_265_proj[3],x6G08_MD_265_proj[2],pch=20)
text(x6G08_MD_265_proj[3],x6G08_MD_265_proj[2],pos=4,label="265",col="blue")

x6G08_MD_266_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_266.txt")
x6G08_MD_266_proj<-project.pca(x6G08_MD_266_raw,pcadata)
points(x6G08_MD_266_proj[3],x6G08_MD_266_proj[2],pch=20)
text(x6G08_MD_266_proj[3],x6G08_MD_266_proj[2],pos=4,label="266",col="blue")

x6G08_MD_267_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_267.txt")
x6G08_MD_267_proj<-project.pca(x6G08_MD_267_raw,pcadata)
points(x6G08_MD_267_proj[3],x6G08_MD_267_proj[2],pch=20)
text(x6G08_MD_267_proj[3],x6G08_MD_267_proj[2],pos=4,label="267",col="blue")

x6G08_MD_268_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_268.txt")
x6G08_MD_268_proj<-project.pca(x6G08_MD_268_raw,pcadata)
points(x6G08_MD_268_proj[3],x6G08_MD_268_proj[2],pch=20)
text(x6G08_MD_268_proj[3],x6G08_MD_268_proj[2],pos=4,label="268",col="blue")

x6G08_MD_269_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_269.txt")
x6G08_MD_269_proj<-project.pca(x6G08_MD_269_raw,pcadata)
points(x6G08_MD_269_proj[3],x6G08_MD_269_proj[2],pch=20)
text(x6G08_MD_269_proj[3],x6G08_MD_269_proj[2],pos=4,label="269",col="blue")

x6G08_MD_270_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_270.txt")
x6G08_MD_270_proj<-project.pca(x6G08_MD_270_raw,pcadata)
points(x6G08_MD_270_proj[3],x6G08_MD_270_proj[2],pch=20)
text(x6G08_MD_270_proj[3],x6G08_MD_270_proj[2],pos=4,label="270",col="black")

x6G08_MD_271_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_271.txt")
x6G08_MD_271_proj<-project.pca(x6G08_MD_271_raw,pcadata)
points(x6G08_MD_271_proj[3],x6G08_MD_271_proj[2],pch=20)
text(x6G08_MD_271_proj[3],x6G08_MD_271_proj[2],pos=4,label="271",col="green")

x6G08_MD_272_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_272.txt")
x6G08_MD_272_proj<-project.pca(x6G08_MD_272_raw,pcadata)
points(x6G08_MD_272_proj[3],x6G08_MD_272_proj[2],pch=20)
text(x6G08_MD_272_proj[3],x6G08_MD_272_proj[2],pos=4,label="272",col="blue")

x6G08_MD_273_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_273.txt")
x6G08_MD_273_proj<-project.pca(x6G08_MD_273_raw,pcadata)
points(x6G08_MD_273_proj[3],x6G08_MD_273_proj[2],pch=20)
text(x6G08_MD_273_proj[3],x6G08_MD_273_proj[2],pos=4,label="273",col="blue")

x6G08_MD_274_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_274.txt")
x6G08_MD_274_proj<-project.pca(x6G08_MD_274_raw,pcadata)
points(x6G08_MD_274_proj[3],x6G08_MD_274_proj[2],pch=20)
text(x6G08_MD_274_proj[3],x6G08_MD_274_proj[2],pos=4,label="274",col="blue")

x6G08_MD_275_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_275.txt")
x6G08_MD_275_proj<-project.pca(x6G08_MD_275_raw,pcadata)
points(x6G08_MD_275_proj[3],x6G08_MD_275_proj[2],pch=20)
text(x6G08_MD_275_proj[3],x6G08_MD_275_proj[2],pos=4,label="275",col="green")

x6G08_MD_276_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_276.txt")
x6G08_MD_276_proj<-project.pca(x6G08_MD_276_raw,pcadata)
points(x6G08_MD_276_proj[3],x6G08_MD_276_proj[2],pch=20)
text(x6G08_MD_276_proj[3],x6G08_MD_276_proj[2],pos=4,label="276",col="blue")

x6G08_MD_277_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_277.txt")
x6G08_MD_277_proj<-project.pca(x6G08_MD_277_raw,pcadata)
points(x6G08_MD_277_proj[3],x6G08_MD_277_proj[2],pch=20)
text(x6G08_MD_277_proj[3],x6G08_MD_277_proj[2],pos=4,label="277",col="green")

x6G08_MD_278_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_278.txt")
x6G08_MD_278_proj<-project.pca(x6G08_MD_278_raw,pcadata)
points(x6G08_MD_278_proj[3],x6G08_MD_278_proj[2],pch=20)
text(x6G08_MD_278_proj[3],x6G08_MD_278_proj[2],pos=4,label="278",col="red")

x6G08_MD_279_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_279.txt")
x6G08_MD_279_proj<-project.pca(x6G08_MD_279_raw,pcadata)
points(x6G08_MD_279_proj[3],x6G08_MD_279_proj[2],pch=20)
text(x6G08_MD_279_proj[3],x6G08_MD_279_proj[2],pos=4,label="279",col="green")

x6G08_MD_280_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_280.txt")
x6G08_MD_280_proj<-project.pca(x6G08_MD_280_raw,pcadata)
points(x6G08_MD_280_proj[3],x6G08_MD_280_proj[2],pch=20)
text(x6G08_MD_280_proj[3],x6G08_MD_280_proj[2],pos=4,label="280",col="blue")

x6G08_MD_281_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_281.txt")
x6G08_MD_281_proj<-project.pca(x6G08_MD_281_raw,pcadata)
points(x6G08_MD_281_proj[3],x6G08_MD_281_proj[2],pch=20)
text(x6G08_MD_281_proj[3],x6G08_MD_281_proj[2],pos=4,label="281",col="green")

x6G08_MD_282_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_282.txt")
x6G08_MD_282_proj<-project.pca(x6G08_MD_282_raw,pcadata)
points(x6G08_MD_282_proj[3],x6G08_MD_282_proj[2],pch=20)
text(x6G08_MD_282_proj[3],x6G08_MD_282_proj[2],pos=4,label="282",col="green")

x6G08_MD_283_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_283.txt")
x6G08_MD_283_proj<-project.pca(x6G08_MD_283_raw,pcadata)
points(x6G08_MD_283_proj[3],x6G08_MD_283_proj[2],pch=20)
text(x6G08_MD_283_proj[3],x6G08_MD_283_proj[2],pos=4,label="283",col="green")

x6G08_MD_284_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_284.txt")
x6G08_MD_284_proj<-project.pca(x6G08_MD_284_raw,pcadata)
points(x6G08_MD_284_proj[3],x6G08_MD_284_proj[2],pch=20)
text(x6G08_MD_284_proj[3],x6G08_MD_284_proj[2],pos=4,label="284",col="green")

x6G08_MD_285_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_285.txt")
x6G08_MD_285_proj<-project.pca(x6G08_MD_285_raw,pcadata)
points(x6G08_MD_285_proj[3],x6G08_MD_285_proj[2],pch=20)
text(x6G08_MD_285_proj[3],x6G08_MD_285_proj[2],pos=4,label="285",col="red")

x6G08_MD_286_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_286.txt")
x6G08_MD_286_proj<-project.pca(x6G08_MD_286_raw,pcadata)
points(x6G08_MD_286_proj[3],x6G08_MD_286_proj[2],pch=20)
text(x6G08_MD_286_proj[3],x6G08_MD_286_proj[2],pos=4,label="286",col="red")

x6G08_MD_287_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_287.txt")
x6G08_MD_287_proj<-project.pca(x6G08_MD_287_raw,pcadata)
points(x6G08_MD_287_proj[3],x6G08_MD_287_proj[2],pch=20)
text(x6G08_MD_287_proj[3],x6G08_MD_287_proj[2],pos=4,label="287",col="red")

x6G08_MD_288_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_288.txt")
x6G08_MD_288_proj<-project.pca(x6G08_MD_288_raw,pcadata)
points(x6G08_MD_288_proj[3],x6G08_MD_288_proj[2],pch=20)
text(x6G08_MD_288_proj[3],x6G08_MD_288_proj[2],pos=4,label="288",col="green")

x6G08_MD_289_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_289.txt")
x6G08_MD_289_proj<-project.pca(x6G08_MD_289_raw,pcadata)
points(x6G08_MD_289_proj[3],x6G08_MD_289_proj[2],pch=20)
text(x6G08_MD_289_proj[3],x6G08_MD_289_proj[2],pos=4,label="289",col="blue")

x6G08_MD_290_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_290.txt")
x6G08_MD_290_proj<-project.pca(x6G08_MD_290_raw,pcadata)
points(x6G08_MD_290_proj[3],x6G08_MD_290_proj[2],pch=20)
text(x6G08_MD_290_proj[3],x6G08_MD_290_proj[2],pos=4,label="290",col="black")

x6G08_MD_291_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_291.txt")
x6G08_MD_291_proj<-project.pca(x6G08_MD_291_raw,pcadata)
points(x6G08_MD_291_proj[3],x6G08_MD_291_proj[2],pch=20)
text(x6G08_MD_291_proj[3],x6G08_MD_291_proj[2],pos=4,label="291",col="blue")

x6G08_MD_292_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_292.txt")
x6G08_MD_292_proj<-project.pca(x6G08_MD_292_raw,pcadata)
points(x6G08_MD_292_proj[3],x6G08_MD_292_proj[2],pch=20)
text(x6G08_MD_292_proj[3],x6G08_MD_292_proj[2],pos=4,label="292",col="black")

x6G08_MD_293_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_293.txt")
x6G08_MD_293_proj<-project.pca(x6G08_MD_293_raw,pcadata)
points(x6G08_MD_293_proj[3],x6G08_MD_293_proj[2],pch=20)
text(x6G08_MD_293_proj[3],x6G08_MD_293_proj[2],pos=4,label="293",col="blue")

x6G08_MD_294_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_294.txt")
x6G08_MD_294_proj<-project.pca(x6G08_MD_294_raw,pcadata)
points(x6G08_MD_294_proj[3],x6G08_MD_294_proj[2],pch=20)
text(x6G08_MD_294_proj[3],x6G08_MD_294_proj[2],pos=4,label="294",col="red")

x6G08_MD_295_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_295.txt")
x6G08_MD_295_proj<-project.pca(x6G08_MD_295_raw,pcadata)
points(x6G08_MD_295_proj[3],x6G08_MD_295_proj[2],pch=20)
text(x6G08_MD_295_proj[3],x6G08_MD_295_proj[2],pos=4,label="295",col="green")

x6G08_MD_296_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_296.txt")
x6G08_MD_296_proj<-project.pca(x6G08_MD_296_raw,pcadata)
points(x6G08_MD_296_proj[3],x6G08_MD_296_proj[2],pch=20)
text(x6G08_MD_296_proj[3],x6G08_MD_296_proj[2],pos=4,label="296",col="blue")

x6G08_MD_297_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_297.txt")
x6G08_MD_297_proj<-project.pca(x6G08_MD_297_raw,pcadata)
points(x6G08_MD_297_proj[3],x6G08_MD_297_proj[2],pch=20)
text(x6G08_MD_297_proj[3],x6G08_MD_297_proj[2],pos=4,label="297",col="red")

x6G08_MD_298_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_298.txt")
x6G08_MD_298_proj<-project.pca(x6G08_MD_298_raw,pcadata)
points(x6G08_MD_298_proj[3],x6G08_MD_298_proj[2],pch=20)
text(x6G08_MD_298_proj[3],x6G08_MD_298_proj[2],pos=4,label="298",col="red")

x6G08_MD_299_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_299.txt")
x6G08_MD_299_proj<-project.pca(x6G08_MD_299_raw,pcadata)
points(x6G08_MD_299_proj[3],x6G08_MD_299_proj[2],pch=20)
text(x6G08_MD_299_proj[3],x6G08_MD_299_proj[2],pos=4,label="299",col="green")

x6G08_MD_300_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_300.txt")
x6G08_MD_300_proj<-project.pca(x6G08_MD_300_raw,pcadata)
points(x6G08_MD_300_proj[3],x6G08_MD_300_proj[2],pch=20)
text(x6G08_MD_300_proj[3],x6G08_MD_300_proj[2],pos=4,label="300",col="green")

x6G08_MD_301_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_301.txt")
x6G08_MD_301_proj<-project.pca(x6G08_MD_301_raw,pcadata)
points(x6G08_MD_301_proj[3],x6G08_MD_301_proj[2],pch=20)
text(x6G08_MD_301_proj[3],x6G08_MD_301_proj[2],pos=4,label="301",col="blue")

x6G08_MD_302_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_302.txt")
x6G08_MD_302_proj<-project.pca(x6G08_MD_302_raw,pcadata)
points(x6G08_MD_302_proj[3],x6G08_MD_302_proj[2],pch=20)
text(x6G08_MD_302_proj[3],x6G08_MD_302_proj[2],pos=4,label="302",col="green")

x6G08_MD_303_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_303.txt")
x6G08_MD_303_proj<-project.pca(x6G08_MD_303_raw,pcadata)
points(x6G08_MD_303_proj[3],x6G08_MD_303_proj[2],pch=20)
text(x6G08_MD_303_proj[3],x6G08_MD_303_proj[2],pos=4,label="303",col="green")

x6G08_MD_304_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_304.txt")
x6G08_MD_304_proj<-project.pca(x6G08_MD_304_raw,pcadata)
points(x6G08_MD_304_proj[3],x6G08_MD_304_proj[2],pch=20)
text(x6G08_MD_304_proj[3],x6G08_MD_304_proj[2],pos=4,label="304",col="green")

x6G08_MD_305_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_305.txt")
x6G08_MD_305_proj<-project.pca(x6G08_MD_305_raw,pcadata)
points(x6G08_MD_305_proj[3],x6G08_MD_305_proj[2],pch=20)
text(x6G08_MD_305_proj[3],x6G08_MD_305_proj[2],pos=4,label="305",col="red")

x6G08_MD_306_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_306.txt")
x6G08_MD_306_proj<-project.pca(x6G08_MD_306_raw,pcadata)
points(x6G08_MD_306_proj[3],x6G08_MD_306_proj[2],pch=20)
text(x6G08_MD_306_proj[3],x6G08_MD_306_proj[2],pos=4,label="306",col="green")

x6G08_MD_307_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_307.txt")
x6G08_MD_307_proj<-project.pca(x6G08_MD_307_raw,pcadata)
points(x6G08_MD_307_proj[3],x6G08_MD_307_proj[2],pch=20)
text(x6G08_MD_307_proj[3],x6G08_MD_307_proj[2],pos=4,label="307",col="green")

x6G08_MD_308_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_308.txt")
x6G08_MD_308_proj<-project.pca(x6G08_MD_308_raw,pcadata)
points(x6G08_MD_308_proj[3],x6G08_MD_308_proj[2],pch=20)
text(x6G08_MD_308_proj[3],x6G08_MD_308_proj[2],pos=4,label="308",col="red")

x6G08_MD_309_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_309.txt")
x6G08_MD_309_proj<-project.pca(x6G08_MD_309_raw,pcadata)
points(x6G08_MD_309_proj[3],x6G08_MD_309_proj[2],pch=20)
text(x6G08_MD_309_proj[3],x6G08_MD_309_proj[2],pos=4,label="309",col="green")

x6G08_MD_310_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_310.txt")
x6G08_MD_310_proj<-project.pca(x6G08_MD_310_raw,pcadata)
points(x6G08_MD_310_proj[3],x6G08_MD_310_proj[2],pch=20)
text(x6G08_MD_310_proj[3],x6G08_MD_310_proj[2],pos=4,label="310",col="green")

x6G08_MD_311_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_311.txt")
x6G08_MD_311_proj<-project.pca(x6G08_MD_311_raw,pcadata)
points(x6G08_MD_311_proj[3],x6G08_MD_311_proj[2],pch=20)
text(x6G08_MD_311_proj[3],x6G08_MD_311_proj[2],pos=4,label="311",col="blue")

x6G08_MD_312_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_312.txt")
x6G08_MD_312_proj<-project.pca(x6G08_MD_312_raw,pcadata)
points(x6G08_MD_312_proj[3],x6G08_MD_312_proj[2],pch=20)
text(x6G08_MD_312_proj[3],x6G08_MD_312_proj[2],pos=4,label="312",col="blue")

x6G08_MD_313_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_313.txt")
x6G08_MD_313_proj<-project.pca(x6G08_MD_313_raw,pcadata)
points(x6G08_MD_313_proj[3],x6G08_MD_313_proj[2],pch=20)
text(x6G08_MD_313_proj[3],x6G08_MD_313_proj[2],pos=4,label="313",col="blue")

x6G08_MD_314_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_314.txt")
x6G08_MD_314_proj<-project.pca(x6G08_MD_314_raw,pcadata)
points(x6G08_MD_314_proj[3],x6G08_MD_314_proj[2],pch=20)
text(x6G08_MD_314_proj[3],x6G08_MD_314_proj[2],pos=4,label="314",col="green")

x6G08_MD_315_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_315.txt")
x6G08_MD_315_proj<-project.pca(x6G08_MD_315_raw,pcadata)
points(x6G08_MD_315_proj[3],x6G08_MD_315_proj[2],pch=20)
text(x6G08_MD_315_proj[3],x6G08_MD_315_proj[2],pos=4,label="315",col="blue")

x6G08_MD_316_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_316.txt")
x6G08_MD_316_proj<-project.pca(x6G08_MD_316_raw,pcadata)
points(x6G08_MD_316_proj[3],x6G08_MD_316_proj[2],pch=20)
text(x6G08_MD_316_proj[3],x6G08_MD_316_proj[2],pos=4,label="316",col="blue")

x6G08_MD_317_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_317.txt")
x6G08_MD_317_proj<-project.pca(x6G08_MD_317_raw,pcadata)
points(x6G08_MD_317_proj[3],x6G08_MD_317_proj[2],pch=20)
text(x6G08_MD_317_proj[3],x6G08_MD_317_proj[2],pos=4,label="317",col="green")

x6G08_MD_318_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_318.txt")
x6G08_MD_318_proj<-project.pca(x6G08_MD_318_raw,pcadata)
points(x6G08_MD_318_proj[3],x6G08_MD_318_proj[2],pch=20)
text(x6G08_MD_318_proj[3],x6G08_MD_318_proj[2],pos=4,label="318",col="green")

x6G08_MD_319_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_319.txt")
x6G08_MD_319_proj<-project.pca(x6G08_MD_319_raw,pcadata)
points(x6G08_MD_319_proj[3],x6G08_MD_319_proj[2],pch=20)
text(x6G08_MD_319_proj[3],x6G08_MD_319_proj[2],pos=4,label="319",col="green")

x6G08_MD_320_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_320.txt")
x6G08_MD_320_proj<-project.pca(x6G08_MD_320_raw,pcadata)
points(x6G08_MD_320_proj[3],x6G08_MD_320_proj[2],pch=20)
text(x6G08_MD_320_proj[3],x6G08_MD_320_proj[2],pos=4,label="320",col="green")

x6G08_MD_321_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_321.txt")
x6G08_MD_321_proj<-project.pca(x6G08_MD_321_raw,pcadata)
points(x6G08_MD_321_proj[3],x6G08_MD_321_proj[2],pch=20)
text(x6G08_MD_321_proj[3],x6G08_MD_321_proj[2],pos=4,label="321",col="blue")

x6G08_MD_322_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_322.txt")
x6G08_MD_322_proj<-project.pca(x6G08_MD_322_raw,pcadata)
points(x6G08_MD_322_proj[3],x6G08_MD_322_proj[2],pch=20)
text(x6G08_MD_322_proj[3],x6G08_MD_322_proj[2],pos=4,label="322",col="green")

x6G08_MD_323_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_323.txt")
x6G08_MD_323_proj<-project.pca(x6G08_MD_323_raw,pcadata)
points(x6G08_MD_323_proj[3],x6G08_MD_323_proj[2],pch=20)
text(x6G08_MD_323_proj[3],x6G08_MD_323_proj[2],pos=4,label="323",col="red")

x6G08_MD_324_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_324.txt")
x6G08_MD_324_proj<-project.pca(x6G08_MD_324_raw,pcadata)
points(x6G08_MD_324_proj[3],x6G08_MD_324_proj[2],pch=20)
text(x6G08_MD_324_proj[3],x6G08_MD_324_proj[2],pos=4,label="324",col="green")

x6G08_MD_325_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_325.txt")
x6G08_MD_325_proj<-project.pca(x6G08_MD_325_raw,pcadata)
points(x6G08_MD_325_proj[3],x6G08_MD_325_proj[2],pch=20)
text(x6G08_MD_325_proj[3],x6G08_MD_325_proj[2],pos=4,label="325",col="red")

x6G08_MD_326_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_326.txt")
x6G08_MD_326_proj<-project.pca(x6G08_MD_326_raw,pcadata)
points(x6G08_MD_326_proj[3],x6G08_MD_326_proj[2],pch=20)
text(x6G08_MD_326_proj[3],x6G08_MD_326_proj[2],pos=4,label="326",col="red")

x6G08_MD_327_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_327.txt")
x6G08_MD_327_proj<-project.pca(x6G08_MD_327_raw,pcadata)
points(x6G08_MD_327_proj[3],x6G08_MD_327_proj[2],pch=20)
text(x6G08_MD_327_proj[3],x6G08_MD_327_proj[2],pos=4,label="327",col="blue")

x6G08_MD_328_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_328.txt")
x6G08_MD_328_proj<-project.pca(x6G08_MD_328_raw,pcadata)
points(x6G08_MD_328_proj[3],x6G08_MD_328_proj[2],pch=20)
text(x6G08_MD_328_proj[3],x6G08_MD_328_proj[2],pos=4,label="328",col="red")

x6G08_MD_329_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_329.txt")
x6G08_MD_329_proj<-project.pca(x6G08_MD_329_raw,pcadata)
points(x6G08_MD_329_proj[3],x6G08_MD_329_proj[2],pch=20)
text(x6G08_MD_329_proj[3],x6G08_MD_329_proj[2],pos=4,label="329",col="red")

x6G08_MD_330_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_330.txt")
x6G08_MD_330_proj<-project.pca(x6G08_MD_330_raw,pcadata)
points(x6G08_MD_330_proj[3],x6G08_MD_330_proj[2],pch=20)
text(x6G08_MD_330_proj[3],x6G08_MD_330_proj[2],pos=4,label="330",col="green")

x6G08_MD_331_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_331.txt")
x6G08_MD_331_proj<-project.pca(x6G08_MD_331_raw,pcadata)
points(x6G08_MD_331_proj[3],x6G08_MD_331_proj[2],pch=20)
text(x6G08_MD_331_proj[3],x6G08_MD_331_proj[2],pos=4,label="331",col="green")

x6G08_MD_332_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_332.txt")
x6G08_MD_332_proj<-project.pca(x6G08_MD_332_raw,pcadata)
points(x6G08_MD_332_proj[3],x6G08_MD_332_proj[2],pch=20)
text(x6G08_MD_332_proj[3],x6G08_MD_332_proj[2],pos=4,label="332",col="red")

x6G08_MD_333_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_333.txt")
x6G08_MD_333_proj<-project.pca(x6G08_MD_333_raw,pcadata)
points(x6G08_MD_333_proj[3],x6G08_MD_333_proj[2],pch=20)
text(x6G08_MD_333_proj[3],x6G08_MD_333_proj[2],pos=4,label="333",col="green")

x6G08_MD_334_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_334.txt")
x6G08_MD_334_proj<-project.pca(x6G08_MD_334_raw,pcadata)
points(x6G08_MD_334_proj[3],x6G08_MD_334_proj[2],pch=20)
text(x6G08_MD_334_proj[3],x6G08_MD_334_proj[2],pos=4,label="334",col="green")

x6G08_MD_335_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_335.txt")
x6G08_MD_335_proj<-project.pca(x6G08_MD_335_raw,pcadata)
points(x6G08_MD_335_proj[3],x6G08_MD_335_proj[2],pch=20)
text(x6G08_MD_335_proj[3],x6G08_MD_335_proj[2],pos=4,label="335",col="green")

x6G08_MD_336_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_336.txt")
x6G08_MD_336_proj<-project.pca(x6G08_MD_336_raw,pcadata)
points(x6G08_MD_336_proj[3],x6G08_MD_336_proj[2],pch=20)
text(x6G08_MD_336_proj[3],x6G08_MD_336_proj[2],pos=4,label="336",col="green")

x6G08_MD_337_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_337.txt")
x6G08_MD_337_proj<-project.pca(x6G08_MD_337_raw,pcadata)
points(x6G08_MD_337_proj[3],x6G08_MD_337_proj[2],pch=20)
text(x6G08_MD_337_proj[3],x6G08_MD_337_proj[2],pos=4,label="337",col="green")

x6G08_MD_338_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_338.txt")
x6G08_MD_338_proj<-project.pca(x6G08_MD_338_raw,pcadata)
points(x6G08_MD_338_proj[3],x6G08_MD_338_proj[2],pch=20)
text(x6G08_MD_338_proj[3],x6G08_MD_338_proj[2],pos=4,label="338",col="blue")

x6G08_MD_339_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_339.txt")
x6G08_MD_339_proj<-project.pca(x6G08_MD_339_raw,pcadata)
points(x6G08_MD_339_proj[3],x6G08_MD_339_proj[2],pch=20)
text(x6G08_MD_339_proj[3],x6G08_MD_339_proj[2],pos=4,label="339",col="black")

x6G08_MD_340_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_340.txt")
x6G08_MD_340_proj<-project.pca(x6G08_MD_340_raw,pcadata)
points(x6G08_MD_340_proj[3],x6G08_MD_340_proj[2],pch=20)
text(x6G08_MD_340_proj[3],x6G08_MD_340_proj[2],pos=4,label="340",col="blue")

x6G08_MD_341_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_341.txt")
x6G08_MD_341_proj<-project.pca(x6G08_MD_341_raw,pcadata)
points(x6G08_MD_341_proj[3],x6G08_MD_341_proj[2],pch=20)
text(x6G08_MD_341_proj[3],x6G08_MD_341_proj[2],pos=4,label="341",col="blue")

x6G08_MD_342_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_342.txt")
x6G08_MD_342_proj<-project.pca(x6G08_MD_342_raw,pcadata)
points(x6G08_MD_342_proj[3],x6G08_MD_342_proj[2],pch=20)
text(x6G08_MD_342_proj[3],x6G08_MD_342_proj[2],pos=4,label="342",col="green")

x6G08_MD_343_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_343.txt")
x6G08_MD_343_proj<-project.pca(x6G08_MD_343_raw,pcadata)
points(x6G08_MD_343_proj[3],x6G08_MD_343_proj[2],pch=20)
text(x6G08_MD_343_proj[3],x6G08_MD_343_proj[2],pos=4,label="343",col="black")

x6G08_MD_344_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_344.txt")
x6G08_MD_344_proj<-project.pca(x6G08_MD_344_raw,pcadata)
points(x6G08_MD_344_proj[3],x6G08_MD_344_proj[2],pch=20)
text(x6G08_MD_344_proj[3],x6G08_MD_344_proj[2],pos=4,label="344",col="blue")

x6G08_MD_345_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_345.txt")
x6G08_MD_345_proj<-project.pca(x6G08_MD_345_raw,pcadata)
points(x6G08_MD_345_proj[3],x6G08_MD_345_proj[2],pch=20)
text(x6G08_MD_345_proj[3],x6G08_MD_345_proj[2],pos=4,label="345",col="blue")

x6G08_MD_346_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_346.txt")
x6G08_MD_346_proj<-project.pca(x6G08_MD_346_raw,pcadata)
points(x6G08_MD_346_proj[3],x6G08_MD_346_proj[2],pch=20)
text(x6G08_MD_346_proj[3],x6G08_MD_346_proj[2],pos=4,label="346",col="black")

x6G08_MD_347_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_347.txt")
x6G08_MD_347_proj<-project.pca(x6G08_MD_347_raw,pcadata)
points(x6G08_MD_347_proj[3],x6G08_MD_347_proj[2],pch=20)
text(x6G08_MD_347_proj[3],x6G08_MD_347_proj[2],pos=4,label="347",col="blue")

x6G08_MD_348_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_348.txt")
x6G08_MD_348_proj<-project.pca(x6G08_MD_348_raw,pcadata)
points(x6G08_MD_348_proj[3],x6G08_MD_348_proj[2],pch=20)
text(x6G08_MD_348_proj[3],x6G08_MD_348_proj[2],pos=4,label="348",col="black")

x6G08_MD_349_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_349.txt")
x6G08_MD_349_proj<-project.pca(x6G08_MD_349_raw,pcadata)
points(x6G08_MD_349_proj[3],x6G08_MD_349_proj[2],pch=20)
text(x6G08_MD_349_proj[3],x6G08_MD_349_proj[2],pos=4,label="349",col="black")

x6G08_MD_350_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_350.txt")
x6G08_MD_350_proj<-project.pca(x6G08_MD_350_raw,pcadata)
points(x6G08_MD_350_proj[3],x6G08_MD_350_proj[2],pch=20)
text(x6G08_MD_350_proj[3],x6G08_MD_350_proj[2],pos=4,label="350",col="black")

x6G08_MD_351_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_351.txt")
x6G08_MD_351_proj<-project.pca(x6G08_MD_351_raw,pcadata)
points(x6G08_MD_351_proj[3],x6G08_MD_351_proj[2],pch=20)
text(x6G08_MD_351_proj[3],x6G08_MD_351_proj[2],pos=4,label="351",col="blue")

x6G08_MD_352_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_352.txt")
x6G08_MD_352_proj<-project.pca(x6G08_MD_352_raw,pcadata)
points(x6G08_MD_352_proj[3],x6G08_MD_352_proj[2],pch=20)
text(x6G08_MD_352_proj[3],x6G08_MD_352_proj[2],pos=4,label="352",col="black")

x6G08_MD_353_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_353.txt")
x6G08_MD_353_proj<-project.pca(x6G08_MD_353_raw,pcadata)
points(x6G08_MD_353_proj[3],x6G08_MD_353_proj[2],pch=20)
text(x6G08_MD_353_proj[3],x6G08_MD_353_proj[2],pos=4,label="353",col="blue")

x6G08_MD_354_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_354.txt")
x6G08_MD_354_proj<-project.pca(x6G08_MD_354_raw,pcadata)
points(x6G08_MD_354_proj[3],x6G08_MD_354_proj[2],pch=20)
text(x6G08_MD_354_proj[3],x6G08_MD_354_proj[2],pos=4,label="354",col="black")

x6G08_MD_355_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_355.txt")
x6G08_MD_355_proj<-project.pca(x6G08_MD_355_raw,pcadata)
points(x6G08_MD_355_proj[3],x6G08_MD_355_proj[2],pch=20)
text(x6G08_MD_355_proj[3],x6G08_MD_355_proj[2],pos=4,label="355",col="black")

x6G08_MD_356_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_356.txt")
x6G08_MD_356_proj<-project.pca(x6G08_MD_356_raw,pcadata)
points(x6G08_MD_356_proj[3],x6G08_MD_356_proj[2],pch=20)
text(x6G08_MD_356_proj[3],x6G08_MD_356_proj[2],pos=4,label="356",col="green")

x6G08_MD_357_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_357.txt")
x6G08_MD_357_proj<-project.pca(x6G08_MD_357_raw,pcadata)
points(x6G08_MD_357_proj[3],x6G08_MD_357_proj[2],pch=20)
text(x6G08_MD_357_proj[3],x6G08_MD_357_proj[2],pos=4,label="357",col="blue")

x6G08_MD_358_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_358.txt")
x6G08_MD_358_proj<-project.pca(x6G08_MD_358_raw,pcadata)
points(x6G08_MD_358_proj[3],x6G08_MD_358_proj[2],pch=20)
text(x6G08_MD_358_proj[3],x6G08_MD_358_proj[2],pos=4,label="358",col="black")

x6G08_MD_359_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_359.txt")
x6G08_MD_359_proj<-project.pca(x6G08_MD_359_raw,pcadata)
points(x6G08_MD_359_proj[3],x6G08_MD_359_proj[2],pch=20)
text(x6G08_MD_359_proj[3],x6G08_MD_359_proj[2],pos=4,label="359",col="black")

x6G08_MD_360_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_360.txt")
x6G08_MD_360_proj<-project.pca(x6G08_MD_360_raw,pcadata)
points(x6G08_MD_360_proj[3],x6G08_MD_360_proj[2],pch=20)
text(x6G08_MD_360_proj[3],x6G08_MD_360_proj[2],pos=4,label="360",col="blue")

x6G08_MD_361_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_361.txt")
x6G08_MD_361_proj<-project.pca(x6G08_MD_361_raw,pcadata)
points(x6G08_MD_361_proj[3],x6G08_MD_361_proj[2],pch=20)
text(x6G08_MD_361_proj[3],x6G08_MD_361_proj[2],pos=4,label="361",col="red")

x6G08_MD_362_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_362.txt")
x6G08_MD_362_proj<-project.pca(x6G08_MD_362_raw,pcadata)
points(x6G08_MD_362_proj[3],x6G08_MD_362_proj[2],pch=20)
text(x6G08_MD_362_proj[3],x6G08_MD_362_proj[2],pos=4,label="362",col="red")

x6G08_MD_363_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_363.txt")
x6G08_MD_363_proj<-project.pca(x6G08_MD_363_raw,pcadata)
points(x6G08_MD_363_proj[3],x6G08_MD_363_proj[2],pch=20)
text(x6G08_MD_363_proj[3],x6G08_MD_363_proj[2],pos=4,label="363",col="red")

x6G08_MD_364_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_364.txt")
x6G08_MD_364_proj<-project.pca(x6G08_MD_364_raw,pcadata)
points(x6G08_MD_364_proj[3],x6G08_MD_364_proj[2],pch=20)
text(x6G08_MD_364_proj[3],x6G08_MD_364_proj[2],pos=4,label="364",col="red")

x6G08_MD_365_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_365.txt")
x6G08_MD_365_proj<-project.pca(x6G08_MD_365_raw,pcadata)
points(x6G08_MD_365_proj[3],x6G08_MD_365_proj[2],pch=20)
text(x6G08_MD_365_proj[3],x6G08_MD_365_proj[2],pos=4,label="365",col="red")

x6G08_MD_366_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_366.txt")
x6G08_MD_366_proj<-project.pca(x6G08_MD_366_raw,pcadata)
points(x6G08_MD_366_proj[3],x6G08_MD_366_proj[2],pch=20)
text(x6G08_MD_366_proj[3],x6G08_MD_366_proj[2],pos=4,label="366",col="red")

x6G08_MD_367_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_367.txt")
x6G08_MD_367_proj<-project.pca(x6G08_MD_367_raw,pcadata)
points(x6G08_MD_367_proj[3],x6G08_MD_367_proj[2],pch=20)
text(x6G08_MD_367_proj[3],x6G08_MD_367_proj[2],pos=4,label="367",col="green")

x6G08_MD_368_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_368.txt")
x6G08_MD_368_proj<-project.pca(x6G08_MD_368_raw,pcadata)
points(x6G08_MD_368_proj[3],x6G08_MD_368_proj[2],pch=20)
text(x6G08_MD_368_proj[3],x6G08_MD_368_proj[2],pos=4,label="368",col="green")

x6G08_MD_369_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_369.txt")
x6G08_MD_369_proj<-project.pca(x6G08_MD_369_raw,pcadata)
points(x6G08_MD_369_proj[3],x6G08_MD_369_proj[2],pch=20)
text(x6G08_MD_369_proj[3],x6G08_MD_369_proj[2],pos=4,label="369",col="red")

x6G08_MD_370_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_370.txt")
x6G08_MD_370_proj<-project.pca(x6G08_MD_370_raw,pcadata)
points(x6G08_MD_370_proj[3],x6G08_MD_370_proj[2],pch=20)
text(x6G08_MD_370_proj[3],x6G08_MD_370_proj[2],pos=4,label="370",col="green")

x6G08_MD_371_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_371.txt")
x6G08_MD_371_proj<-project.pca(x6G08_MD_371_raw,pcadata)
points(x6G08_MD_371_proj[3],x6G08_MD_371_proj[2],pch=20)
text(x6G08_MD_371_proj[3],x6G08_MD_371_proj[2],pos=4,label="371",col="blue")

x6G08_MD_372_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_372.txt")
x6G08_MD_372_proj<-project.pca(x6G08_MD_372_raw,pcadata)
points(x6G08_MD_372_proj[3],x6G08_MD_372_proj[2],pch=20)
text(x6G08_MD_372_proj[3],x6G08_MD_372_proj[2],pos=4,label="372",col="blue")

x6G08_MD_373_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_373.txt")
x6G08_MD_373_proj<-project.pca(x6G08_MD_373_raw,pcadata)
points(x6G08_MD_373_proj[3],x6G08_MD_373_proj[2],pch=20)
text(x6G08_MD_373_proj[3],x6G08_MD_373_proj[2],pos=4,label="373",col="blue")

x6G08_MD_374_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_374.txt")
x6G08_MD_374_proj<-project.pca(x6G08_MD_374_raw,pcadata)
points(x6G08_MD_374_proj[3],x6G08_MD_374_proj[2],pch=20)
text(x6G08_MD_374_proj[3],x6G08_MD_374_proj[2],pos=4,label="374",col="blue")

x6G08_MD_375_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_375.txt")
x6G08_MD_375_proj<-project.pca(x6G08_MD_375_raw,pcadata)
points(x6G08_MD_375_proj[3],x6G08_MD_375_proj[2],pch=20)
text(x6G08_MD_375_proj[3],x6G08_MD_375_proj[2],pos=4,label="375",col="red")

x6G08_MD_376_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_376.txt")
x6G08_MD_376_proj<-project.pca(x6G08_MD_376_raw,pcadata)
points(x6G08_MD_376_proj[3],x6G08_MD_376_proj[2],pch=20)
text(x6G08_MD_376_proj[3],x6G08_MD_376_proj[2],pos=4,label="376",col="blue")

x6G08_MD_377_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_377.txt")
x6G08_MD_377_proj<-project.pca(x6G08_MD_377_raw,pcadata)
points(x6G08_MD_377_proj[3],x6G08_MD_377_proj[2],pch=20)
text(x6G08_MD_377_proj[3],x6G08_MD_377_proj[2],pos=4,label="377",col="blue")

x6G08_MD_378_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_378.txt")
x6G08_MD_378_proj<-project.pca(x6G08_MD_378_raw,pcadata)
points(x6G08_MD_378_proj[3],x6G08_MD_378_proj[2],pch=20)
text(x6G08_MD_378_proj[3],x6G08_MD_378_proj[2],pos=4,label="378",col="black")

x6G08_MD_379_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_379.txt")
x6G08_MD_379_proj<-project.pca(x6G08_MD_379_raw,pcadata)
points(x6G08_MD_379_proj[3],x6G08_MD_379_proj[2],pch=20)
text(x6G08_MD_379_proj[3],x6G08_MD_379_proj[2],pos=4,label="379",col="blue")

x6G08_MD_380_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_380.txt")
x6G08_MD_380_proj<-project.pca(x6G08_MD_380_raw,pcadata)
points(x6G08_MD_380_proj[3],x6G08_MD_380_proj[2],pch=20)
text(x6G08_MD_380_proj[3],x6G08_MD_380_proj[2],pos=4,label="380",col="blue")

x6G08_MD_381_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_381.txt")
x6G08_MD_381_proj<-project.pca(x6G08_MD_381_raw,pcadata)
points(x6G08_MD_381_proj[3],x6G08_MD_381_proj[2],pch=20)
text(x6G08_MD_381_proj[3],x6G08_MD_381_proj[2],pos=4,label="381",col="black")

x6G08_MD_382_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_382.txt")
x6G08_MD_382_proj<-project.pca(x6G08_MD_382_raw,pcadata)
points(x6G08_MD_382_proj[3],x6G08_MD_382_proj[2],pch=20)
text(x6G08_MD_382_proj[3],x6G08_MD_382_proj[2],pos=4,label="382",col="blue")

x6G08_MD_383_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_383.txt")
x6G08_MD_383_proj<-project.pca(x6G08_MD_383_raw,pcadata)
points(x6G08_MD_383_proj[3],x6G08_MD_383_proj[2],pch=20)
text(x6G08_MD_383_proj[3],x6G08_MD_383_proj[2],pos=4,label="383",col="blue")

x6G08_MD_384_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_384.txt")
x6G08_MD_384_proj<-project.pca(x6G08_MD_384_raw,pcadata)
points(x6G08_MD_384_proj[3],x6G08_MD_384_proj[2],pch=20)
text(x6G08_MD_384_proj[3],x6G08_MD_384_proj[2],pos=4,label="384",col="black")

x6G08_MD_385_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_385.txt")
x6G08_MD_385_proj<-project.pca(x6G08_MD_385_raw,pcadata)
points(x6G08_MD_385_proj[3],x6G08_MD_385_proj[2],pch=20)
text(x6G08_MD_385_proj[3],x6G08_MD_385_proj[2],pos=4,label="385",col="blue")

x6G08_MD_386_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_386.txt")
x6G08_MD_386_proj<-project.pca(x6G08_MD_386_raw,pcadata)
points(x6G08_MD_386_proj[3],x6G08_MD_386_proj[2],pch=20)
text(x6G08_MD_386_proj[3],x6G08_MD_386_proj[2],pos=4,label="386",col="blue")

x6G08_MD_387_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_387.txt")
x6G08_MD_387_proj<-project.pca(x6G08_MD_387_raw,pcadata)
points(x6G08_MD_387_proj[3],x6G08_MD_387_proj[2],pch=20)
text(x6G08_MD_387_proj[3],x6G08_MD_387_proj[2],pos=4,label="387",col="blue")

x6G08_MD_388_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_388.txt")
x6G08_MD_388_proj<-project.pca(x6G08_MD_388_raw,pcadata)
points(x6G08_MD_388_proj[3],x6G08_MD_388_proj[2],pch=20)
text(x6G08_MD_388_proj[3],x6G08_MD_388_proj[2],pos=4,label="388",col="red")

x6G08_MD_389_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_389.txt")
x6G08_MD_389_proj<-project.pca(x6G08_MD_389_raw,pcadata)
points(x6G08_MD_389_proj[3],x6G08_MD_389_proj[2],pch=20)
text(x6G08_MD_389_proj[3],x6G08_MD_389_proj[2],pos=4,label="389",col="green")

x6G08_MD_390_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_390.txt")
x6G08_MD_390_proj<-project.pca(x6G08_MD_390_raw,pcadata)
points(x6G08_MD_390_proj[3],x6G08_MD_390_proj[2],pch=20)
text(x6G08_MD_390_proj[3],x6G08_MD_390_proj[2],pos=4,label="390",col="blue")

x6G08_MD_391_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_391.txt")
x6G08_MD_391_proj<-project.pca(x6G08_MD_391_raw,pcadata)
points(x6G08_MD_391_proj[3],x6G08_MD_391_proj[2],pch=20)
text(x6G08_MD_391_proj[3],x6G08_MD_391_proj[2],pos=4,label="391",col="green")

x6G08_MD_392_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_392.txt")
x6G08_MD_392_proj<-project.pca(x6G08_MD_392_raw,pcadata)
points(x6G08_MD_392_proj[3],x6G08_MD_392_proj[2],pch=20)
text(x6G08_MD_392_proj[3],x6G08_MD_392_proj[2],pos=4,label="392",col="green")

x6G08_MD_393_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_393.txt")
x6G08_MD_393_proj<-project.pca(x6G08_MD_393_raw,pcadata)
points(x6G08_MD_393_proj[3],x6G08_MD_393_proj[2],pch=20)
text(x6G08_MD_393_proj[3],x6G08_MD_393_proj[2],pos=4,label="393",col="green")

x6G08_MD_394_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_394.txt")
x6G08_MD_394_proj<-project.pca(x6G08_MD_394_raw,pcadata)
points(x6G08_MD_394_proj[3],x6G08_MD_394_proj[2],pch=20)
text(x6G08_MD_394_proj[3],x6G08_MD_394_proj[2],pos=4,label="394",col="green")

x6G08_MD_395_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_395.txt")
x6G08_MD_395_proj<-project.pca(x6G08_MD_395_raw,pcadata)
points(x6G08_MD_395_proj[3],x6G08_MD_395_proj[2],pch=20)
text(x6G08_MD_395_proj[3],x6G08_MD_395_proj[2],pos=4,label="395",col="green")

x6G08_MD_396_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_396.txt")
x6G08_MD_396_proj<-project.pca(x6G08_MD_396_raw,pcadata)
points(x6G08_MD_396_proj[3],x6G08_MD_396_proj[2],pch=20)
text(x6G08_MD_396_proj[3],x6G08_MD_396_proj[2],pos=4,label="396",col="green")

x6G08_MD_397_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_397.txt")
x6G08_MD_397_proj<-project.pca(x6G08_MD_397_raw,pcadata)
points(x6G08_MD_397_proj[3],x6G08_MD_397_proj[2],pch=20)
text(x6G08_MD_397_proj[3],x6G08_MD_397_proj[2],pos=4,label="397",col="blue")

x6G08_MD_398_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_398.txt")
x6G08_MD_398_proj<-project.pca(x6G08_MD_398_raw,pcadata)
points(x6G08_MD_398_proj[3],x6G08_MD_398_proj[2],pch=20)
text(x6G08_MD_398_proj[3],x6G08_MD_398_proj[2],pos=4,label="398",col="blue")

x6G08_MD_399_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_399.txt")
x6G08_MD_399_proj<-project.pca(x6G08_MD_399_raw,pcadata)
points(x6G08_MD_399_proj[3],x6G08_MD_399_proj[2],pch=20)
text(x6G08_MD_399_proj[3],x6G08_MD_399_proj[2],pos=4,label="399",col="green")

x6G08_MD_400_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_400.txt")
x6G08_MD_400_proj<-project.pca(x6G08_MD_400_raw,pcadata)
points(x6G08_MD_400_proj[3],x6G08_MD_400_proj[2],pch=20)
text(x6G08_MD_400_proj[3],x6G08_MD_400_proj[2],pos=4,label="400",col="blue")

x6G08_MD_401_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_401.txt")
x6G08_MD_401_proj<-project.pca(x6G08_MD_401_raw,pcadata)
points(x6G08_MD_401_proj[3],x6G08_MD_401_proj[2],pch=20)
text(x6G08_MD_401_proj[3],x6G08_MD_401_proj[2],pos=4,label="401",col="green")

x6G08_MD_402_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_402.txt")
x6G08_MD_402_proj<-project.pca(x6G08_MD_402_raw,pcadata)
points(x6G08_MD_402_proj[3],x6G08_MD_402_proj[2],pch=20)
text(x6G08_MD_402_proj[3],x6G08_MD_402_proj[2],pos=4,label="402",col="green")

x6G08_MD_403_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_403.txt")
x6G08_MD_403_proj<-project.pca(x6G08_MD_403_raw,pcadata)
points(x6G08_MD_403_proj[3],x6G08_MD_403_proj[2],pch=20)
text(x6G08_MD_403_proj[3],x6G08_MD_403_proj[2],pos=4,label="403",col="red")

x6G08_MD_404_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_404.txt")
x6G08_MD_404_proj<-project.pca(x6G08_MD_404_raw,pcadata)
points(x6G08_MD_404_proj[3],x6G08_MD_404_proj[2],pch=20)
text(x6G08_MD_404_proj[3],x6G08_MD_404_proj[2],pos=4,label="404",col="blue")

x6G08_MD_405_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_405.txt")
x6G08_MD_405_proj<-project.pca(x6G08_MD_405_raw,pcadata)
points(x6G08_MD_405_proj[3],x6G08_MD_405_proj[2],pch=20)
text(x6G08_MD_405_proj[3],x6G08_MD_405_proj[2],pos=4,label="405",col="red")

x6G08_MD_406_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_406.txt")
x6G08_MD_406_proj<-project.pca(x6G08_MD_406_raw,pcadata)
points(x6G08_MD_406_proj[3],x6G08_MD_406_proj[2],pch=20)
text(x6G08_MD_406_proj[3],x6G08_MD_406_proj[2],pos=4,label="406",col="red")

x6G08_MD_407_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_407.txt")
x6G08_MD_407_proj<-project.pca(x6G08_MD_407_raw,pcadata)
points(x6G08_MD_407_proj[3],x6G08_MD_407_proj[2],pch=20)
text(x6G08_MD_407_proj[3],x6G08_MD_407_proj[2],pos=4,label="407",col="red")

x6G08_MD_408_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_408.txt")
x6G08_MD_408_proj<-project.pca(x6G08_MD_408_raw,pcadata)
points(x6G08_MD_408_proj[3],x6G08_MD_408_proj[2],pch=20)
text(x6G08_MD_408_proj[3],x6G08_MD_408_proj[2],pos=4,label="408",col="red")

x6G08_MD_409_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_409.txt")
x6G08_MD_409_proj<-project.pca(x6G08_MD_409_raw,pcadata)
points(x6G08_MD_409_proj[3],x6G08_MD_409_proj[2],pch=20)
text(x6G08_MD_409_proj[3],x6G08_MD_409_proj[2],pos=4,label="409",col="red")

x6G08_MD_410_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_410.txt")
x6G08_MD_410_proj<-project.pca(x6G08_MD_410_raw,pcadata)
points(x6G08_MD_410_proj[3],x6G08_MD_410_proj[2],pch=20)
text(x6G08_MD_410_proj[3],x6G08_MD_410_proj[2],pos=4,label="410",col="green")

x6G08_MD_411_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_411.txt")
x6G08_MD_411_proj<-project.pca(x6G08_MD_411_raw,pcadata)
points(x6G08_MD_411_proj[3],x6G08_MD_411_proj[2],pch=20)
text(x6G08_MD_411_proj[3],x6G08_MD_411_proj[2],pos=4,label="411",col="red")

x6G08_MD_412_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_412.txt")
x6G08_MD_412_proj<-project.pca(x6G08_MD_412_raw,pcadata)
points(x6G08_MD_412_proj[3],x6G08_MD_412_proj[2],pch=20)
text(x6G08_MD_412_proj[3],x6G08_MD_412_proj[2],pos=4,label="412",col="red")

x6G08_MD_413_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_413.txt")
x6G08_MD_413_proj<-project.pca(x6G08_MD_413_raw,pcadata)
points(x6G08_MD_413_proj[3],x6G08_MD_413_proj[2],pch=20)
text(x6G08_MD_413_proj[3],x6G08_MD_413_proj[2],pos=4,label="413",col="green")

x6G08_MD_414_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_414.txt")
x6G08_MD_414_proj<-project.pca(x6G08_MD_414_raw,pcadata)
points(x6G08_MD_414_proj[3],x6G08_MD_414_proj[2],pch=20)
text(x6G08_MD_414_proj[3],x6G08_MD_414_proj[2],pos=4,label="414",col="green")

x6G08_MD_415_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_415.txt")
x6G08_MD_415_proj<-project.pca(x6G08_MD_415_raw,pcadata)
points(x6G08_MD_415_proj[3],x6G08_MD_415_proj[2],pch=20)
text(x6G08_MD_415_proj[3],x6G08_MD_415_proj[2],pos=4,label="415",col="green")

x6G08_MD_416_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_416.txt")
x6G08_MD_416_proj<-project.pca(x6G08_MD_416_raw,pcadata)
points(x6G08_MD_416_proj[3],x6G08_MD_416_proj[2],pch=20)
text(x6G08_MD_416_proj[3],x6G08_MD_416_proj[2],pos=4,label="416",col="red")

x6G08_MD_417_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_417.txt")
x6G08_MD_417_proj<-project.pca(x6G08_MD_417_raw,pcadata)
points(x6G08_MD_417_proj[3],x6G08_MD_417_proj[2],pch=20)
text(x6G08_MD_417_proj[3],x6G08_MD_417_proj[2],pos=4,label="417",col="red")

x6G08_MD_418_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_418.txt")
x6G08_MD_418_proj<-project.pca(x6G08_MD_418_raw,pcadata)
points(x6G08_MD_418_proj[3],x6G08_MD_418_proj[2],pch=20)
text(x6G08_MD_418_proj[3],x6G08_MD_418_proj[2],pos=4,label="418",col="red")

x6G08_MD_419_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_419.txt")
x6G08_MD_419_proj<-project.pca(x6G08_MD_419_raw,pcadata)
points(x6G08_MD_419_proj[3],x6G08_MD_419_proj[2],pch=20)
text(x6G08_MD_419_proj[3],x6G08_MD_419_proj[2],pos=4,label="419",col="green")

x6G08_MD_420_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_420.txt")
x6G08_MD_420_proj<-project.pca(x6G08_MD_420_raw,pcadata)
points(x6G08_MD_420_proj[3],x6G08_MD_420_proj[2],pch=20)
text(x6G08_MD_420_proj[3],x6G08_MD_420_proj[2],pos=4,label="420",col="red")

x6G08_MD_421_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_421.txt")
x6G08_MD_421_proj<-project.pca(x6G08_MD_421_raw,pcadata)
points(x6G08_MD_421_proj[3],x6G08_MD_421_proj[2],pch=20)
text(x6G08_MD_421_proj[3],x6G08_MD_421_proj[2],pos=4,label="421",col="red")

x6G08_MD_422_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_422.txt")
x6G08_MD_422_proj<-project.pca(x6G08_MD_422_raw,pcadata)
points(x6G08_MD_422_proj[3],x6G08_MD_422_proj[2],pch=20)
text(x6G08_MD_422_proj[3],x6G08_MD_422_proj[2],pos=4,label="422",col="red")

x6G08_MD_423_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_423.txt")
x6G08_MD_423_proj<-project.pca(x6G08_MD_423_raw,pcadata)
points(x6G08_MD_423_proj[3],x6G08_MD_423_proj[2],pch=20)
text(x6G08_MD_423_proj[3],x6G08_MD_423_proj[2],pos=4,label="423",col="blue")

x6G08_MD_424_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_424.txt")
x6G08_MD_424_proj<-project.pca(x6G08_MD_424_raw,pcadata)
points(x6G08_MD_424_proj[3],x6G08_MD_424_proj[2],pch=20)
text(x6G08_MD_424_proj[3],x6G08_MD_424_proj[2],pos=4,label="424",col="red")

x6G08_MD_425_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_425.txt")
x6G08_MD_425_proj<-project.pca(x6G08_MD_425_raw,pcadata)
points(x6G08_MD_425_proj[3],x6G08_MD_425_proj[2],pch=20)
text(x6G08_MD_425_proj[3],x6G08_MD_425_proj[2],pos=4,label="425",col="red")

x6G08_MD_426_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_426.txt")
x6G08_MD_426_proj<-project.pca(x6G08_MD_426_raw,pcadata)
points(x6G08_MD_426_proj[3],x6G08_MD_426_proj[2],pch=20)
text(x6G08_MD_426_proj[3],x6G08_MD_426_proj[2],pos=4,label="426",col="red")

x6G08_MD_427_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_427.txt")
x6G08_MD_427_proj<-project.pca(x6G08_MD_427_raw,pcadata)
points(x6G08_MD_427_proj[3],x6G08_MD_427_proj[2],pch=20)
text(x6G08_MD_427_proj[3],x6G08_MD_427_proj[2],pos=4,label="427",col="green")

x6G08_MD_428_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_428.txt")
x6G08_MD_428_proj<-project.pca(x6G08_MD_428_raw,pcadata)
points(x6G08_MD_428_proj[3],x6G08_MD_428_proj[2],pch=20)
text(x6G08_MD_428_proj[3],x6G08_MD_428_proj[2],pos=4,label="428",col="red")

x6G08_MD_429_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_429.txt")
x6G08_MD_429_proj<-project.pca(x6G08_MD_429_raw,pcadata)
points(x6G08_MD_429_proj[3],x6G08_MD_429_proj[2],pch=20)
text(x6G08_MD_429_proj[3],x6G08_MD_429_proj[2],pos=4,label="429",col="green")

x6G08_MD_430_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_430.txt")
x6G08_MD_430_proj<-project.pca(x6G08_MD_430_raw,pcadata)
points(x6G08_MD_430_proj[3],x6G08_MD_430_proj[2],pch=20)
text(x6G08_MD_430_proj[3],x6G08_MD_430_proj[2],pos=4,label="430",col="red")

x6G08_MD_431_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_431.txt")
x6G08_MD_431_proj<-project.pca(x6G08_MD_431_raw,pcadata)
points(x6G08_MD_431_proj[3],x6G08_MD_431_proj[2],pch=20)
text(x6G08_MD_431_proj[3],x6G08_MD_431_proj[2],pos=4,label="431",col="green")

x6G08_MD_432_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_432.txt")
x6G08_MD_432_proj<-project.pca(x6G08_MD_432_raw,pcadata)
points(x6G08_MD_432_proj[3],x6G08_MD_432_proj[2],pch=20)
text(x6G08_MD_432_proj[3],x6G08_MD_432_proj[2],pos=4,label="432",col="green")

x6G08_MD_433_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_433.txt")
x6G08_MD_433_proj<-project.pca(x6G08_MD_433_raw,pcadata)
points(x6G08_MD_433_proj[3],x6G08_MD_433_proj[2],pch=20)
text(x6G08_MD_433_proj[3],x6G08_MD_433_proj[2],pos=4,label="433",col="green")

x6G08_MD_434_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_434.txt")
x6G08_MD_434_proj<-project.pca(x6G08_MD_434_raw,pcadata)
points(x6G08_MD_434_proj[3],x6G08_MD_434_proj[2],pch=20)
text(x6G08_MD_434_proj[3],x6G08_MD_434_proj[2],pos=4,label="434",col="red")

x6G08_MD_435_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_435.txt")
x6G08_MD_435_proj<-project.pca(x6G08_MD_435_raw,pcadata)
points(x6G08_MD_435_proj[3],x6G08_MD_435_proj[2],pch=20)
text(x6G08_MD_435_proj[3],x6G08_MD_435_proj[2],pos=4,label="435",col="green")

x6G08_MD_436_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_436.txt")
x6G08_MD_436_proj<-project.pca(x6G08_MD_436_raw,pcadata)
points(x6G08_MD_436_proj[3],x6G08_MD_436_proj[2],pch=20)
text(x6G08_MD_436_proj[3],x6G08_MD_436_proj[2],pos=4,label="436",col="red")

x6G08_MD_437_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_437.txt")
x6G08_MD_437_proj<-project.pca(x6G08_MD_437_raw,pcadata)
points(x6G08_MD_437_proj[3],x6G08_MD_437_proj[2],pch=20)
text(x6G08_MD_437_proj[3],x6G08_MD_437_proj[2],pos=4,label="437",col="green")

x6G08_MD_438_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_438.txt")
x6G08_MD_438_proj<-project.pca(x6G08_MD_438_raw,pcadata)
points(x6G08_MD_438_proj[3],x6G08_MD_438_proj[2],pch=20)
text(x6G08_MD_438_proj[3],x6G08_MD_438_proj[2],pos=4,label="438",col="red")

x6G08_MD_439_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_439.txt")
x6G08_MD_439_proj<-project.pca(x6G08_MD_439_raw,pcadata)
points(x6G08_MD_439_proj[3],x6G08_MD_439_proj[2],pch=20)
text(x6G08_MD_439_proj[3],x6G08_MD_439_proj[2],pos=4,label="439",col="red")

x6G08_MD_440_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_440.txt")
x6G08_MD_440_proj<-project.pca(x6G08_MD_440_raw,pcadata)
points(x6G08_MD_440_proj[3],x6G08_MD_440_proj[2],pch=20)
text(x6G08_MD_440_proj[3],x6G08_MD_440_proj[2],pos=4,label="440",col="red")

x6G08_MD_441_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_441.txt")
x6G08_MD_441_proj<-project.pca(x6G08_MD_441_raw,pcadata)
points(x6G08_MD_441_proj[3],x6G08_MD_441_proj[2],pch=20)
text(x6G08_MD_441_proj[3],x6G08_MD_441_proj[2],pos=4,label="441",col="red")

x6G08_MD_442_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_442.txt")
x6G08_MD_442_proj<-project.pca(x6G08_MD_442_raw,pcadata)
points(x6G08_MD_442_proj[3],x6G08_MD_442_proj[2],pch=20)
text(x6G08_MD_442_proj[3],x6G08_MD_442_proj[2],pos=4,label="442",col="red")

x6G08_MD_443_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_443.txt")
x6G08_MD_443_proj<-project.pca(x6G08_MD_443_raw,pcadata)
points(x6G08_MD_443_proj[3],x6G08_MD_443_proj[2],pch=20)
text(x6G08_MD_443_proj[3],x6G08_MD_443_proj[2],pos=4,label="443",col="red")

x6G08_MD_444_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_444.txt")
x6G08_MD_444_proj<-project.pca(x6G08_MD_444_raw,pcadata)
points(x6G08_MD_444_proj[3],x6G08_MD_444_proj[2],pch=20)
text(x6G08_MD_444_proj[3],x6G08_MD_444_proj[2],pos=4,label="444",col="red")

x6G08_MD_445_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_445.txt")
x6G08_MD_445_proj<-project.pca(x6G08_MD_445_raw,pcadata)
points(x6G08_MD_445_proj[3],x6G08_MD_445_proj[2],pch=20)
text(x6G08_MD_445_proj[3],x6G08_MD_445_proj[2],pos=4,label="445",col="red")

x6G08_MD_446_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_446.txt")
x6G08_MD_446_proj<-project.pca(x6G08_MD_446_raw,pcadata)
points(x6G08_MD_446_proj[3],x6G08_MD_446_proj[2],pch=20)
text(x6G08_MD_446_proj[3],x6G08_MD_446_proj[2],pos=4,label="446",col="red")

x6G08_MD_447_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_447.txt")
x6G08_MD_447_proj<-project.pca(x6G08_MD_447_raw,pcadata)
points(x6G08_MD_447_proj[3],x6G08_MD_447_proj[2],pch=20)
text(x6G08_MD_447_proj[3],x6G08_MD_447_proj[2],pos=4,label="447",col="red")

x6G08_MD_448_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_448.txt")
x6G08_MD_448_proj<-project.pca(x6G08_MD_448_raw,pcadata)
points(x6G08_MD_448_proj[3],x6G08_MD_448_proj[2],pch=20)
text(x6G08_MD_448_proj[3],x6G08_MD_448_proj[2],pos=4,label="448",col="red")

x6G08_MD_449_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_449.txt")
x6G08_MD_449_proj<-project.pca(x6G08_MD_449_raw,pcadata)
points(x6G08_MD_449_proj[3],x6G08_MD_449_proj[2],pch=20)
text(x6G08_MD_449_proj[3],x6G08_MD_449_proj[2],pos=4,label="449",col="red")

x6G08_MD_450_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_450.txt")
x6G08_MD_450_proj<-project.pca(x6G08_MD_450_raw,pcadata)
points(x6G08_MD_450_proj[3],x6G08_MD_450_proj[2],pch=20)
text(x6G08_MD_450_proj[3],x6G08_MD_450_proj[2],pos=4,label="450",col="red")

x6G08_MD_451_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_451.txt")
x6G08_MD_451_proj<-project.pca(x6G08_MD_451_raw,pcadata)
points(x6G08_MD_451_proj[3],x6G08_MD_451_proj[2],pch=20)
text(x6G08_MD_451_proj[3],x6G08_MD_451_proj[2],pos=4,label="451",col="red")

x6G08_MD_452_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_452.txt")
x6G08_MD_452_proj<-project.pca(x6G08_MD_452_raw,pcadata)
points(x6G08_MD_452_proj[3],x6G08_MD_452_proj[2],pch=20)
text(x6G08_MD_452_proj[3],x6G08_MD_452_proj[2],pos=4,label="452",col="red")

x6G08_MD_453_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_453.txt")
x6G08_MD_453_proj<-project.pca(x6G08_MD_453_raw,pcadata)
points(x6G08_MD_453_proj[3],x6G08_MD_453_proj[2],pch=20)
text(x6G08_MD_453_proj[3],x6G08_MD_453_proj[2],pos=4,label="453",col="red")

x6G08_MD_454_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_454.txt")
x6G08_MD_454_proj<-project.pca(x6G08_MD_454_raw,pcadata)
points(x6G08_MD_454_proj[3],x6G08_MD_454_proj[2],pch=20)
text(x6G08_MD_454_proj[3],x6G08_MD_454_proj[2],pos=4,label="454",col="green")

x6G08_MD_455_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_455.txt")
x6G08_MD_455_proj<-project.pca(x6G08_MD_455_raw,pcadata)
points(x6G08_MD_455_proj[3],x6G08_MD_455_proj[2],pch=20)
text(x6G08_MD_455_proj[3],x6G08_MD_455_proj[2],pos=4,label="455",col="green")

x6G08_MD_456_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_456.txt")
x6G08_MD_456_proj<-project.pca(x6G08_MD_456_raw,pcadata)
points(x6G08_MD_456_proj[3],x6G08_MD_456_proj[2],pch=20)
text(x6G08_MD_456_proj[3],x6G08_MD_456_proj[2],pos=4,label="456",col="red")

x6G08_MD_457_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_457.txt")
x6G08_MD_457_proj<-project.pca(x6G08_MD_457_raw,pcadata)
points(x6G08_MD_457_proj[3],x6G08_MD_457_proj[2],pch=20)
text(x6G08_MD_457_proj[3],x6G08_MD_457_proj[2],pos=4,label="457",col="red")

x6G08_MD_458_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_458.txt")
x6G08_MD_458_proj<-project.pca(x6G08_MD_458_raw,pcadata)
points(x6G08_MD_458_proj[3],x6G08_MD_458_proj[2],pch=20)
text(x6G08_MD_458_proj[3],x6G08_MD_458_proj[2],pos=4,label="458",col="red")

x6G08_MD_459_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_459.txt")
x6G08_MD_459_proj<-project.pca(x6G08_MD_459_raw,pcadata)
points(x6G08_MD_459_proj[3],x6G08_MD_459_proj[2],pch=20)
text(x6G08_MD_459_proj[3],x6G08_MD_459_proj[2],pos=4,label="459",col="red")

x6G08_MD_460_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_460.txt")
x6G08_MD_460_proj<-project.pca(x6G08_MD_460_raw,pcadata)
points(x6G08_MD_460_proj[3],x6G08_MD_460_proj[2],pch=20)
text(x6G08_MD_460_proj[3],x6G08_MD_460_proj[2],pos=4,label="460",col="blue")

x6G08_MD_461_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_461.txt")
x6G08_MD_461_proj<-project.pca(x6G08_MD_461_raw,pcadata)
points(x6G08_MD_461_proj[3],x6G08_MD_461_proj[2],pch=20)
text(x6G08_MD_461_proj[3],x6G08_MD_461_proj[2],pos=4,label="461",col="blue")

x6G08_MD_462_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_462.txt")
x6G08_MD_462_proj<-project.pca(x6G08_MD_462_raw,pcadata)
points(x6G08_MD_462_proj[3],x6G08_MD_462_proj[2],pch=20)
text(x6G08_MD_462_proj[3],x6G08_MD_462_proj[2],pos=4,label="462",col="green")

x6G08_MD_463_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_463.txt")
x6G08_MD_463_proj<-project.pca(x6G08_MD_463_raw,pcadata)
points(x6G08_MD_463_proj[3],x6G08_MD_463_proj[2],pch=20)
text(x6G08_MD_463_proj[3],x6G08_MD_463_proj[2],pos=4,label="463",col="blue")

x6G08_MD_464_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_464.txt")
x6G08_MD_464_proj<-project.pca(x6G08_MD_464_raw,pcadata)
points(x6G08_MD_464_proj[3],x6G08_MD_464_proj[2],pch=20)
text(x6G08_MD_464_proj[3],x6G08_MD_464_proj[2],pos=4,label="464",col="green")

x6G08_MD_465_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_465.txt")
x6G08_MD_465_proj<-project.pca(x6G08_MD_465_raw,pcadata)
points(x6G08_MD_465_proj[3],x6G08_MD_465_proj[2],pch=20)
text(x6G08_MD_465_proj[3],x6G08_MD_465_proj[2],pos=4,label="465",col="green")

x6G08_MD_466_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_466.txt")
x6G08_MD_466_proj<-project.pca(x6G08_MD_466_raw,pcadata)
points(x6G08_MD_466_proj[3],x6G08_MD_466_proj[2],pch=20)
text(x6G08_MD_466_proj[3],x6G08_MD_466_proj[2],pos=4,label="466",col="blue")

x6G08_MD_467_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_467.txt")
x6G08_MD_467_proj<-project.pca(x6G08_MD_467_raw,pcadata)
points(x6G08_MD_467_proj[3],x6G08_MD_467_proj[2],pch=20)
text(x6G08_MD_467_proj[3],x6G08_MD_467_proj[2],pos=4,label="467",col="blue")

x6G08_MD_468_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_468.txt")
x6G08_MD_468_proj<-project.pca(x6G08_MD_468_raw,pcadata)
points(x6G08_MD_468_proj[3],x6G08_MD_468_proj[2],pch=20)
text(x6G08_MD_468_proj[3],x6G08_MD_468_proj[2],pos=4,label="468",col="blue")

x6G08_MD_469_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_469.txt")
x6G08_MD_469_proj<-project.pca(x6G08_MD_469_raw,pcadata)
points(x6G08_MD_469_proj[3],x6G08_MD_469_proj[2],pch=20)
text(x6G08_MD_469_proj[3],x6G08_MD_469_proj[2],pos=4,label="469",col="black")

x6G08_MD_470_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_470.txt")
x6G08_MD_470_proj<-project.pca(x6G08_MD_470_raw,pcadata)
points(x6G08_MD_470_proj[3],x6G08_MD_470_proj[2],pch=20)
text(x6G08_MD_470_proj[3],x6G08_MD_470_proj[2],pos=4,label="470",col="blue")

x6G08_MD_471_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_471.txt")
x6G08_MD_471_proj<-project.pca(x6G08_MD_471_raw,pcadata)
points(x6G08_MD_471_proj[3],x6G08_MD_471_proj[2],pch=20)
text(x6G08_MD_471_proj[3],x6G08_MD_471_proj[2],pos=4,label="471",col="green")

x6G08_MD_472_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_472.txt")
x6G08_MD_472_proj<-project.pca(x6G08_MD_472_raw,pcadata)
points(x6G08_MD_472_proj[3],x6G08_MD_472_proj[2],pch=20)
text(x6G08_MD_472_proj[3],x6G08_MD_472_proj[2],pos=4,label="472",col="blue")

x6G08_MD_473_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_473.txt")
x6G08_MD_473_proj<-project.pca(x6G08_MD_473_raw,pcadata)
points(x6G08_MD_473_proj[3],x6G08_MD_473_proj[2],pch=20)
text(x6G08_MD_473_proj[3],x6G08_MD_473_proj[2],pos=4,label="473",col="red")

x6G08_MD_474_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_474.txt")
x6G08_MD_474_proj<-project.pca(x6G08_MD_474_raw,pcadata)
points(x6G08_MD_474_proj[3],x6G08_MD_474_proj[2],pch=20)
text(x6G08_MD_474_proj[3],x6G08_MD_474_proj[2],pos=4,label="474",col="green")

x6G08_MD_475_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_475.txt")
x6G08_MD_475_proj<-project.pca(x6G08_MD_475_raw,pcadata)
points(x6G08_MD_475_proj[3],x6G08_MD_475_proj[2],pch=20)
text(x6G08_MD_475_proj[3],x6G08_MD_475_proj[2],pos=4,label="475",col="green")

x6G08_MD_476_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_476.txt")
x6G08_MD_476_proj<-project.pca(x6G08_MD_476_raw,pcadata)
points(x6G08_MD_476_proj[3],x6G08_MD_476_proj[2],pch=20)
text(x6G08_MD_476_proj[3],x6G08_MD_476_proj[2],pos=4,label="476",col="red")

x6G08_MD_477_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_477.txt")
x6G08_MD_477_proj<-project.pca(x6G08_MD_477_raw,pcadata)
points(x6G08_MD_477_proj[3],x6G08_MD_477_proj[2],pch=20)
text(x6G08_MD_477_proj[3],x6G08_MD_477_proj[2],pos=4,label="477",col="green")

x6G08_MD_478_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_478.txt")
x6G08_MD_478_proj<-project.pca(x6G08_MD_478_raw,pcadata)
points(x6G08_MD_478_proj[3],x6G08_MD_478_proj[2],pch=20)
text(x6G08_MD_478_proj[3],x6G08_MD_478_proj[2],pos=4,label="478",col="blue")

x6G08_MD_479_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_479.txt")
x6G08_MD_479_proj<-project.pca(x6G08_MD_479_raw,pcadata)
points(x6G08_MD_479_proj[3],x6G08_MD_479_proj[2],pch=20)
text(x6G08_MD_479_proj[3],x6G08_MD_479_proj[2],pos=4,label="479",col="black")

x6G08_MD_480_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_480.txt")
x6G08_MD_480_proj<-project.pca(x6G08_MD_480_raw,pcadata)
points(x6G08_MD_480_proj[3],x6G08_MD_480_proj[2],pch=20)
text(x6G08_MD_480_proj[3],x6G08_MD_480_proj[2],pos=4,label="480",col="blue")

x6G08_MD_481_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_481.txt")
x6G08_MD_481_proj<-project.pca(x6G08_MD_481_raw,pcadata)
points(x6G08_MD_481_proj[3],x6G08_MD_481_proj[2],pch=20)
text(x6G08_MD_481_proj[3],x6G08_MD_481_proj[2],pos=4,label="481",col="blue")

x6G08_MD_482_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_482.txt")
x6G08_MD_482_proj<-project.pca(x6G08_MD_482_raw,pcadata)
points(x6G08_MD_482_proj[3],x6G08_MD_482_proj[2],pch=20)
text(x6G08_MD_482_proj[3],x6G08_MD_482_proj[2],pos=4,label="482",col="black")

x6G08_MD_483_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_483.txt")
x6G08_MD_483_proj<-project.pca(x6G08_MD_483_raw,pcadata)
points(x6G08_MD_483_proj[3],x6G08_MD_483_proj[2],pch=20)
text(x6G08_MD_483_proj[3],x6G08_MD_483_proj[2],pos=4,label="483",col="blue")

x6G08_MD_484_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_484.txt")
x6G08_MD_484_proj<-project.pca(x6G08_MD_484_raw,pcadata)
points(x6G08_MD_484_proj[3],x6G08_MD_484_proj[2],pch=20)
text(x6G08_MD_484_proj[3],x6G08_MD_484_proj[2],pos=4,label="484",col="black")

x6G08_MD_485_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_485.txt")
x6G08_MD_485_proj<-project.pca(x6G08_MD_485_raw,pcadata)
points(x6G08_MD_485_proj[3],x6G08_MD_485_proj[2],pch=20)
text(x6G08_MD_485_proj[3],x6G08_MD_485_proj[2],pos=4,label="485",col="blue")

x6G08_MD_486_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_486.txt")
x6G08_MD_486_proj<-project.pca(x6G08_MD_486_raw,pcadata)
points(x6G08_MD_486_proj[3],x6G08_MD_486_proj[2],pch=20)
text(x6G08_MD_486_proj[3],x6G08_MD_486_proj[2],pos=4,label="486",col="black")

x6G08_MD_487_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_487.txt")
x6G08_MD_487_proj<-project.pca(x6G08_MD_487_raw,pcadata)
points(x6G08_MD_487_proj[3],x6G08_MD_487_proj[2],pch=20)
text(x6G08_MD_487_proj[3],x6G08_MD_487_proj[2],pos=4,label="487",col="blue")

x6G08_MD_488_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_488.txt")
x6G08_MD_488_proj<-project.pca(x6G08_MD_488_raw,pcadata)
points(x6G08_MD_488_proj[3],x6G08_MD_488_proj[2],pch=20)
text(x6G08_MD_488_proj[3],x6G08_MD_488_proj[2],pos=4,label="488",col="blue")

x6G08_MD_489_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_489.txt")
x6G08_MD_489_proj<-project.pca(x6G08_MD_489_raw,pcadata)
points(x6G08_MD_489_proj[3],x6G08_MD_489_proj[2],pch=20)
text(x6G08_MD_489_proj[3],x6G08_MD_489_proj[2],pos=4,label="489",col="blue")

x6G08_MD_490_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_490.txt")
x6G08_MD_490_proj<-project.pca(x6G08_MD_490_raw,pcadata)
points(x6G08_MD_490_proj[3],x6G08_MD_490_proj[2],pch=20)
text(x6G08_MD_490_proj[3],x6G08_MD_490_proj[2],pos=4,label="490",col="blue")

x6G08_MD_491_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_491.txt")
x6G08_MD_491_proj<-project.pca(x6G08_MD_491_raw,pcadata)
points(x6G08_MD_491_proj[3],x6G08_MD_491_proj[2],pch=20)
text(x6G08_MD_491_proj[3],x6G08_MD_491_proj[2],pos=4,label="491",col="blue")

x6G08_MD_492_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_492.txt")
x6G08_MD_492_proj<-project.pca(x6G08_MD_492_raw,pcadata)
points(x6G08_MD_492_proj[3],x6G08_MD_492_proj[2],pch=20)
text(x6G08_MD_492_proj[3],x6G08_MD_492_proj[2],pos=4,label="492",col="blue")

x6G08_MD_493_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_493.txt")
x6G08_MD_493_proj<-project.pca(x6G08_MD_493_raw,pcadata)
points(x6G08_MD_493_proj[3],x6G08_MD_493_proj[2],pch=20)
text(x6G08_MD_493_proj[3],x6G08_MD_493_proj[2],pos=4,label="493",col="green")

x6G08_MD_494_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_494.txt")
x6G08_MD_494_proj<-project.pca(x6G08_MD_494_raw,pcadata)
points(x6G08_MD_494_proj[3],x6G08_MD_494_proj[2],pch=20)
text(x6G08_MD_494_proj[3],x6G08_MD_494_proj[2],pos=4,label="494",col="green")

x6G08_MD_495_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_495.txt")
x6G08_MD_495_proj<-project.pca(x6G08_MD_495_raw,pcadata)
points(x6G08_MD_495_proj[3],x6G08_MD_495_proj[2],pch=20)
text(x6G08_MD_495_proj[3],x6G08_MD_495_proj[2],pos=4,label="495",col="red")

x6G08_MD_496_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_496.txt")
x6G08_MD_496_proj<-project.pca(x6G08_MD_496_raw,pcadata)
points(x6G08_MD_496_proj[3],x6G08_MD_496_proj[2],pch=20)
text(x6G08_MD_496_proj[3],x6G08_MD_496_proj[2],pos=4,label="496",col="green")

x6G08_MD_497_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_497.txt")
x6G08_MD_497_proj<-project.pca(x6G08_MD_497_raw,pcadata)
points(x6G08_MD_497_proj[3],x6G08_MD_497_proj[2],pch=20)
text(x6G08_MD_497_proj[3],x6G08_MD_497_proj[2],pos=4,label="497",col="green")

x6G08_MD_498_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_498.txt")
x6G08_MD_498_proj<-project.pca(x6G08_MD_498_raw,pcadata)
points(x6G08_MD_498_proj[3],x6G08_MD_498_proj[2],pch=20)
text(x6G08_MD_498_proj[3],x6G08_MD_498_proj[2],pos=4,label="498",col="black")

x6G08_MD_499_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_499.txt")
x6G08_MD_499_proj<-project.pca(x6G08_MD_499_raw,pcadata)
points(x6G08_MD_499_proj[3],x6G08_MD_499_proj[2],pch=20)
text(x6G08_MD_499_proj[3],x6G08_MD_499_proj[2],pos=4,label="499",col="green")

x6G08_MD_500_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_500.txt")
x6G08_MD_500_proj<-project.pca(x6G08_MD_500_raw,pcadata)
points(x6G08_MD_500_proj[3],x6G08_MD_500_proj[2],pch=20)
text(x6G08_MD_500_proj[3],x6G08_MD_500_proj[2],pos=4,label="500",col="red")

x6G08_MD_501_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_501.txt")
x6G08_MD_501_proj<-project.pca(x6G08_MD_501_raw,pcadata)
points(x6G08_MD_501_proj[3],x6G08_MD_501_proj[2],pch=20)
text(x6G08_MD_501_proj[3],x6G08_MD_501_proj[2],pos=4,label="501",col="red")

x6G08_MD_502_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_502.txt")
x6G08_MD_502_proj<-project.pca(x6G08_MD_502_raw,pcadata)
points(x6G08_MD_502_proj[3],x6G08_MD_502_proj[2],pch=20)
text(x6G08_MD_502_proj[3],x6G08_MD_502_proj[2],pos=4,label="502",col="blue")

x6G08_MD_503_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_503.txt")
x6G08_MD_503_proj<-project.pca(x6G08_MD_503_raw,pcadata)
points(x6G08_MD_503_proj[3],x6G08_MD_503_proj[2],pch=20)
text(x6G08_MD_503_proj[3],x6G08_MD_503_proj[2],pos=4,label="503",col="black")

x6G08_MD_504_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_504.txt")
x6G08_MD_504_proj<-project.pca(x6G08_MD_504_raw,pcadata)
points(x6G08_MD_504_proj[3],x6G08_MD_504_proj[2],pch=20)
text(x6G08_MD_504_proj[3],x6G08_MD_504_proj[2],pos=4,label="504",col="blue")

x6G08_MD_505_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_505.txt")
x6G08_MD_505_proj<-project.pca(x6G08_MD_505_raw,pcadata)
points(x6G08_MD_505_proj[3],x6G08_MD_505_proj[2],pch=20)
text(x6G08_MD_505_proj[3],x6G08_MD_505_proj[2],pos=4,label="505",col="green")

x6G08_MD_506_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_506.txt")
x6G08_MD_506_proj<-project.pca(x6G08_MD_506_raw,pcadata)
points(x6G08_MD_506_proj[3],x6G08_MD_506_proj[2],pch=20)
text(x6G08_MD_506_proj[3],x6G08_MD_506_proj[2],pos=4,label="506",col="blue")

x6G08_MD_507_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_507.txt")
x6G08_MD_507_proj<-project.pca(x6G08_MD_507_raw,pcadata)
points(x6G08_MD_507_proj[3],x6G08_MD_507_proj[2],pch=20)
text(x6G08_MD_507_proj[3],x6G08_MD_507_proj[2],pos=4,label="507",col="blue")

x6G08_MD_508_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_508.txt")
x6G08_MD_508_proj<-project.pca(x6G08_MD_508_raw,pcadata)
points(x6G08_MD_508_proj[3],x6G08_MD_508_proj[2],pch=20)
text(x6G08_MD_508_proj[3],x6G08_MD_508_proj[2],pos=4,label="508",col="blue")

x6G08_MD_509_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_509.txt")
x6G08_MD_509_proj<-project.pca(x6G08_MD_509_raw,pcadata)
points(x6G08_MD_509_proj[3],x6G08_MD_509_proj[2],pch=20)
text(x6G08_MD_509_proj[3],x6G08_MD_509_proj[2],pos=4,label="509",col="blue")

x6G08_MD_510_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_510.txt")
x6G08_MD_510_proj<-project.pca(x6G08_MD_510_raw,pcadata)
points(x6G08_MD_510_proj[3],x6G08_MD_510_proj[2],pch=20)
text(x6G08_MD_510_proj[3],x6G08_MD_510_proj[2],pos=4,label="510",col="blue")

x6G08_MD_511_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_511.txt")
x6G08_MD_511_proj<-project.pca(x6G08_MD_511_raw,pcadata)
points(x6G08_MD_511_proj[3],x6G08_MD_511_proj[2],pch=20)
text(x6G08_MD_511_proj[3],x6G08_MD_511_proj[2],pos=4,label="511",col="blue")

x6G08_MD_512_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_512.txt")
x6G08_MD_512_proj<-project.pca(x6G08_MD_512_raw,pcadata)
points(x6G08_MD_512_proj[3],x6G08_MD_512_proj[2],pch=20)
text(x6G08_MD_512_proj[3],x6G08_MD_512_proj[2],pos=4,label="512",col="blue")

x6G08_MD_513_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_513.txt")
x6G08_MD_513_proj<-project.pca(x6G08_MD_513_raw,pcadata)
points(x6G08_MD_513_proj[3],x6G08_MD_513_proj[2],pch=20)
text(x6G08_MD_513_proj[3],x6G08_MD_513_proj[2],pos=4,label="513",col="black")

x6G08_MD_514_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_514.txt")
x6G08_MD_514_proj<-project.pca(x6G08_MD_514_raw,pcadata)
points(x6G08_MD_514_proj[3],x6G08_MD_514_proj[2],pch=20)
text(x6G08_MD_514_proj[3],x6G08_MD_514_proj[2],pos=4,label="514",col="black")

x6G08_MD_515_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_515.txt")
x6G08_MD_515_proj<-project.pca(x6G08_MD_515_raw,pcadata)
points(x6G08_MD_515_proj[3],x6G08_MD_515_proj[2],pch=20)
text(x6G08_MD_515_proj[3],x6G08_MD_515_proj[2],pos=4,label="515",col="blue")

x6G08_MD_516_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_516.txt")
x6G08_MD_516_proj<-project.pca(x6G08_MD_516_raw,pcadata)
points(x6G08_MD_516_proj[3],x6G08_MD_516_proj[2],pch=20)
text(x6G08_MD_516_proj[3],x6G08_MD_516_proj[2],pos=4,label="516",col="black")

x6G08_MD_517_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_517.txt")
x6G08_MD_517_proj<-project.pca(x6G08_MD_517_raw,pcadata)
points(x6G08_MD_517_proj[3],x6G08_MD_517_proj[2],pch=20)
text(x6G08_MD_517_proj[3],x6G08_MD_517_proj[2],pos=4,label="517",col="black")

x6G08_MD_518_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_518.txt")
x6G08_MD_518_proj<-project.pca(x6G08_MD_518_raw,pcadata)
points(x6G08_MD_518_proj[3],x6G08_MD_518_proj[2],pch=20)
text(x6G08_MD_518_proj[3],x6G08_MD_518_proj[2],pos=4,label="518",col="black")

x6G08_MD_519_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_519.txt")
x6G08_MD_519_proj<-project.pca(x6G08_MD_519_raw,pcadata)
points(x6G08_MD_519_proj[3],x6G08_MD_519_proj[2],pch=20)
text(x6G08_MD_519_proj[3],x6G08_MD_519_proj[2],pos=4,label="519",col="black")

x6G08_MD_520_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_520.txt")
x6G08_MD_520_proj<-project.pca(x6G08_MD_520_raw,pcadata)
points(x6G08_MD_520_proj[3],x6G08_MD_520_proj[2],pch=20)
text(x6G08_MD_520_proj[3],x6G08_MD_520_proj[2],pos=4,label="520",col="blue")

x6G08_MD_521_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_521.txt")
x6G08_MD_521_proj<-project.pca(x6G08_MD_521_raw,pcadata)
points(x6G08_MD_521_proj[3],x6G08_MD_521_proj[2],pch=20)
text(x6G08_MD_521_proj[3],x6G08_MD_521_proj[2],pos=4,label="521",col="blue")

x6G08_MD_522_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_522.txt")
x6G08_MD_522_proj<-project.pca(x6G08_MD_522_raw,pcadata)
points(x6G08_MD_522_proj[3],x6G08_MD_522_proj[2],pch=20)
text(x6G08_MD_522_proj[3],x6G08_MD_522_proj[2],pos=4,label="522",col="black")

x6G08_MD_523_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_523.txt")
x6G08_MD_523_proj<-project.pca(x6G08_MD_523_raw,pcadata)
points(x6G08_MD_523_proj[3],x6G08_MD_523_proj[2],pch=20)
text(x6G08_MD_523_proj[3],x6G08_MD_523_proj[2],pos=4,label="523",col="black")

x6G08_MD_524_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_524.txt")
x6G08_MD_524_proj<-project.pca(x6G08_MD_524_raw,pcadata)
points(x6G08_MD_524_proj[3],x6G08_MD_524_proj[2],pch=20)
text(x6G08_MD_524_proj[3],x6G08_MD_524_proj[2],pos=4,label="524",col="black")

x6G08_MD_525_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_525.txt")
x6G08_MD_525_proj<-project.pca(x6G08_MD_525_raw,pcadata)
points(x6G08_MD_525_proj[3],x6G08_MD_525_proj[2],pch=20)
text(x6G08_MD_525_proj[3],x6G08_MD_525_proj[2],pos=4,label="525",col="blue")

x6G08_MD_526_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_526.txt")
x6G08_MD_526_proj<-project.pca(x6G08_MD_526_raw,pcadata)
points(x6G08_MD_526_proj[3],x6G08_MD_526_proj[2],pch=20)
text(x6G08_MD_526_proj[3],x6G08_MD_526_proj[2],pos=4,label="526",col="green")

x6G08_MD_527_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_527.txt")
x6G08_MD_527_proj<-project.pca(x6G08_MD_527_raw,pcadata)
points(x6G08_MD_527_proj[3],x6G08_MD_527_proj[2],pch=20)
text(x6G08_MD_527_proj[3],x6G08_MD_527_proj[2],pos=4,label="527",col="blue")

x6G08_MD_528_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_528.txt")
x6G08_MD_528_proj<-project.pca(x6G08_MD_528_raw,pcadata)
points(x6G08_MD_528_proj[3],x6G08_MD_528_proj[2],pch=20)
text(x6G08_MD_528_proj[3],x6G08_MD_528_proj[2],pos=4,label="528",col="green")

x6G08_MD_529_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_529.txt")
x6G08_MD_529_proj<-project.pca(x6G08_MD_529_raw,pcadata)
points(x6G08_MD_529_proj[3],x6G08_MD_529_proj[2],pch=20)
text(x6G08_MD_529_proj[3],x6G08_MD_529_proj[2],pos=4,label="529",col="blue")

x6G08_MD_530_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_530.txt")
x6G08_MD_530_proj<-project.pca(x6G08_MD_530_raw,pcadata)
points(x6G08_MD_530_proj[3],x6G08_MD_530_proj[2],pch=20)
text(x6G08_MD_530_proj[3],x6G08_MD_530_proj[2],pos=4,label="530",col="blue")

x6G08_MD_531_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_531.txt")
x6G08_MD_531_proj<-project.pca(x6G08_MD_531_raw,pcadata)
points(x6G08_MD_531_proj[3],x6G08_MD_531_proj[2],pch=20)
text(x6G08_MD_531_proj[3],x6G08_MD_531_proj[2],pos=4,label="531",col="blue")

x6G08_MD_532_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_532.txt")
x6G08_MD_532_proj<-project.pca(x6G08_MD_532_raw,pcadata)
points(x6G08_MD_532_proj[3],x6G08_MD_532_proj[2],pch=20)
text(x6G08_MD_532_proj[3],x6G08_MD_532_proj[2],pos=4,label="532",col="red")

x6G08_MD_533_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_533.txt")
x6G08_MD_533_proj<-project.pca(x6G08_MD_533_raw,pcadata)
points(x6G08_MD_533_proj[3],x6G08_MD_533_proj[2],pch=20)
text(x6G08_MD_533_proj[3],x6G08_MD_533_proj[2],pos=4,label="533",col="blue")

x6G08_MD_534_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_534.txt")
x6G08_MD_534_proj<-project.pca(x6G08_MD_534_raw,pcadata)
points(x6G08_MD_534_proj[3],x6G08_MD_534_proj[2],pch=20)
text(x6G08_MD_534_proj[3],x6G08_MD_534_proj[2],pos=4,label="534",col="blue")

x6G08_MD_535_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_535.txt")
x6G08_MD_535_proj<-project.pca(x6G08_MD_535_raw,pcadata)
points(x6G08_MD_535_proj[3],x6G08_MD_535_proj[2],pch=20)
text(x6G08_MD_535_proj[3],x6G08_MD_535_proj[2],pos=4,label="535",col="red")

x6G08_MD_536_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_536.txt")
x6G08_MD_536_proj<-project.pca(x6G08_MD_536_raw,pcadata)
points(x6G08_MD_536_proj[3],x6G08_MD_536_proj[2],pch=20)
text(x6G08_MD_536_proj[3],x6G08_MD_536_proj[2],pos=4,label="536",col="black")

x6G08_MD_537_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_537.txt")
x6G08_MD_537_proj<-project.pca(x6G08_MD_537_raw,pcadata)
points(x6G08_MD_537_proj[3],x6G08_MD_537_proj[2],pch=20)
text(x6G08_MD_537_proj[3],x6G08_MD_537_proj[2],pos=4,label="537",col="blue")

x6G08_MD_538_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_538.txt")
x6G08_MD_538_proj<-project.pca(x6G08_MD_538_raw,pcadata)
points(x6G08_MD_538_proj[3],x6G08_MD_538_proj[2],pch=20)
text(x6G08_MD_538_proj[3],x6G08_MD_538_proj[2],pos=4,label="538",col="red")

x6G08_MD_539_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_539.txt")
x6G08_MD_539_proj<-project.pca(x6G08_MD_539_raw,pcadata)
points(x6G08_MD_539_proj[3],x6G08_MD_539_proj[2],pch=20)
text(x6G08_MD_539_proj[3],x6G08_MD_539_proj[2],pos=4,label="539",col="green")

x6G08_MD_540_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_540.txt")
x6G08_MD_540_proj<-project.pca(x6G08_MD_540_raw,pcadata)
points(x6G08_MD_540_proj[3],x6G08_MD_540_proj[2],pch=20)
text(x6G08_MD_540_proj[3],x6G08_MD_540_proj[2],pos=4,label="540",col="red")

x6G08_MD_541_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_541.txt")
x6G08_MD_541_proj<-project.pca(x6G08_MD_541_raw,pcadata)
points(x6G08_MD_541_proj[3],x6G08_MD_541_proj[2],pch=20)
text(x6G08_MD_541_proj[3],x6G08_MD_541_proj[2],pos=4,label="541",col="red")

x6G08_MD_542_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_542.txt")
x6G08_MD_542_proj<-project.pca(x6G08_MD_542_raw,pcadata)
points(x6G08_MD_542_proj[3],x6G08_MD_542_proj[2],pch=20)
text(x6G08_MD_542_proj[3],x6G08_MD_542_proj[2],pos=4,label="542",col="green")

x6G08_MD_543_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_543.txt")
x6G08_MD_543_proj<-project.pca(x6G08_MD_543_raw,pcadata)
points(x6G08_MD_543_proj[3],x6G08_MD_543_proj[2],pch=20)
text(x6G08_MD_543_proj[3],x6G08_MD_543_proj[2],pos=4,label="543",col="green")

x6G08_MD_544_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_544.txt")
x6G08_MD_544_proj<-project.pca(x6G08_MD_544_raw,pcadata)
points(x6G08_MD_544_proj[3],x6G08_MD_544_proj[2],pch=20)
text(x6G08_MD_544_proj[3],x6G08_MD_544_proj[2],pos=4,label="544",col="green")

x6G08_MD_545_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_545.txt")
x6G08_MD_545_proj<-project.pca(x6G08_MD_545_raw,pcadata)
points(x6G08_MD_545_proj[3],x6G08_MD_545_proj[2],pch=20)
text(x6G08_MD_545_proj[3],x6G08_MD_545_proj[2],pos=4,label="545",col="blue")

x6G08_MD_546_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_546.txt")
x6G08_MD_546_proj<-project.pca(x6G08_MD_546_raw,pcadata)
points(x6G08_MD_546_proj[3],x6G08_MD_546_proj[2],pch=20)
text(x6G08_MD_546_proj[3],x6G08_MD_546_proj[2],pos=4,label="546",col="red")

x6G08_MD_547_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_547.txt")
x6G08_MD_547_proj<-project.pca(x6G08_MD_547_raw,pcadata)
points(x6G08_MD_547_proj[3],x6G08_MD_547_proj[2],pch=20)
text(x6G08_MD_547_proj[3],x6G08_MD_547_proj[2],pos=4,label="547",col="green")

x6G08_MD_548_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_548.txt")
x6G08_MD_548_proj<-project.pca(x6G08_MD_548_raw,pcadata)
points(x6G08_MD_548_proj[3],x6G08_MD_548_proj[2],pch=20)
text(x6G08_MD_548_proj[3],x6G08_MD_548_proj[2],pos=4,label="548",col="blue")

x6G08_MD_549_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_549.txt")
x6G08_MD_549_proj<-project.pca(x6G08_MD_549_raw,pcadata)
points(x6G08_MD_549_proj[3],x6G08_MD_549_proj[2],pch=20)
text(x6G08_MD_549_proj[3],x6G08_MD_549_proj[2],pos=4,label="549",col="red")

x6G08_MD_550_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_550.txt")
x6G08_MD_550_proj<-project.pca(x6G08_MD_550_raw,pcadata)
points(x6G08_MD_550_proj[3],x6G08_MD_550_proj[2],pch=20)
text(x6G08_MD_550_proj[3],x6G08_MD_550_proj[2],pos=4,label="550",col="green")

x6G08_MD_551_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_551.txt")
x6G08_MD_551_proj<-project.pca(x6G08_MD_551_raw,pcadata)
points(x6G08_MD_551_proj[3],x6G08_MD_551_proj[2],pch=20)
text(x6G08_MD_551_proj[3],x6G08_MD_551_proj[2],pos=4,label="551",col="red")

x6G08_MD_552_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_552.txt")
x6G08_MD_552_proj<-project.pca(x6G08_MD_552_raw,pcadata)
points(x6G08_MD_552_proj[3],x6G08_MD_552_proj[2],pch=20)
text(x6G08_MD_552_proj[3],x6G08_MD_552_proj[2],pos=4,label="552",col="red")

x6G08_MD_553_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_553.txt")
x6G08_MD_553_proj<-project.pca(x6G08_MD_553_raw,pcadata)
points(x6G08_MD_553_proj[3],x6G08_MD_553_proj[2],pch=20)
text(x6G08_MD_553_proj[3],x6G08_MD_553_proj[2],pos=4,label="553",col="green")

x6G08_MD_554_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_554.txt")
x6G08_MD_554_proj<-project.pca(x6G08_MD_554_raw,pcadata)
points(x6G08_MD_554_proj[3],x6G08_MD_554_proj[2],pch=20)
text(x6G08_MD_554_proj[3],x6G08_MD_554_proj[2],pos=4,label="554",col="green")

x6G08_MD_555_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_555.txt")
x6G08_MD_555_proj<-project.pca(x6G08_MD_555_raw,pcadata)
points(x6G08_MD_555_proj[3],x6G08_MD_555_proj[2],pch=20)
text(x6G08_MD_555_proj[3],x6G08_MD_555_proj[2],pos=4,label="555",col="blue")

x6G08_MD_556_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_556.txt")
x6G08_MD_556_proj<-project.pca(x6G08_MD_556_raw,pcadata)
points(x6G08_MD_556_proj[3],x6G08_MD_556_proj[2],pch=20)
text(x6G08_MD_556_proj[3],x6G08_MD_556_proj[2],pos=4,label="556",col="green")

x6G08_MD_557_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_557.txt")
x6G08_MD_557_proj<-project.pca(x6G08_MD_557_raw,pcadata)
points(x6G08_MD_557_proj[3],x6G08_MD_557_proj[2],pch=20)
text(x6G08_MD_557_proj[3],x6G08_MD_557_proj[2],pos=4,label="557",col="blue")

x6G08_MD_558_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_558.txt")
x6G08_MD_558_proj<-project.pca(x6G08_MD_558_raw,pcadata)
points(x6G08_MD_558_proj[3],x6G08_MD_558_proj[2],pch=20)
text(x6G08_MD_558_proj[3],x6G08_MD_558_proj[2],pos=4,label="558",col="blue")

x6G08_MD_559_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_559.txt")
x6G08_MD_559_proj<-project.pca(x6G08_MD_559_raw,pcadata)
points(x6G08_MD_559_proj[3],x6G08_MD_559_proj[2],pch=20)
text(x6G08_MD_559_proj[3],x6G08_MD_559_proj[2],pos=4,label="559",col="blue")

x6G08_MD_560_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_560.txt")
x6G08_MD_560_proj<-project.pca(x6G08_MD_560_raw,pcadata)
points(x6G08_MD_560_proj[3],x6G08_MD_560_proj[2],pch=20)
text(x6G08_MD_560_proj[3],x6G08_MD_560_proj[2],pos=4,label="560",col="blue")

x6G08_MD_561_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_561.txt")
x6G08_MD_561_proj<-project.pca(x6G08_MD_561_raw,pcadata)
points(x6G08_MD_561_proj[3],x6G08_MD_561_proj[2],pch=20)
text(x6G08_MD_561_proj[3],x6G08_MD_561_proj[2],pos=4,label="561",col="red")

x6G08_MD_562_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_562.txt")
x6G08_MD_562_proj<-project.pca(x6G08_MD_562_raw,pcadata)
points(x6G08_MD_562_proj[3],x6G08_MD_562_proj[2],pch=20)
text(x6G08_MD_562_proj[3],x6G08_MD_562_proj[2],pos=4,label="562",col="green")

x6G08_MD_563_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_563.txt")
x6G08_MD_563_proj<-project.pca(x6G08_MD_563_raw,pcadata)
points(x6G08_MD_563_proj[3],x6G08_MD_563_proj[2],pch=20)
text(x6G08_MD_563_proj[3],x6G08_MD_563_proj[2],pos=4,label="563",col="red")

x6G08_MD_564_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_564.txt")
x6G08_MD_564_proj<-project.pca(x6G08_MD_564_raw,pcadata)
points(x6G08_MD_564_proj[3],x6G08_MD_564_proj[2],pch=20)
text(x6G08_MD_564_proj[3],x6G08_MD_564_proj[2],pos=4,label="564",col="blue")

x6G08_MD_565_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_565.txt")
x6G08_MD_565_proj<-project.pca(x6G08_MD_565_raw,pcadata)
points(x6G08_MD_565_proj[3],x6G08_MD_565_proj[2],pch=20)
text(x6G08_MD_565_proj[3],x6G08_MD_565_proj[2],pos=4,label="565",col="blue")

x6G08_MD_566_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_566.txt")
x6G08_MD_566_proj<-project.pca(x6G08_MD_566_raw,pcadata)
points(x6G08_MD_566_proj[3],x6G08_MD_566_proj[2],pch=20)
text(x6G08_MD_566_proj[3],x6G08_MD_566_proj[2],pos=4,label="566",col="red")

x6G08_MD_567_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_567.txt")
x6G08_MD_567_proj<-project.pca(x6G08_MD_567_raw,pcadata)
points(x6G08_MD_567_proj[3],x6G08_MD_567_proj[2],pch=20)
text(x6G08_MD_567_proj[3],x6G08_MD_567_proj[2],pos=4,label="567",col="red")

x6G08_MD_568_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_568.txt")
x6G08_MD_568_proj<-project.pca(x6G08_MD_568_raw,pcadata)
points(x6G08_MD_568_proj[3],x6G08_MD_568_proj[2],pch=20)
text(x6G08_MD_568_proj[3],x6G08_MD_568_proj[2],pos=4,label="568",col="blue")

x6G08_MD_569_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_569.txt")
x6G08_MD_569_proj<-project.pca(x6G08_MD_569_raw,pcadata)
points(x6G08_MD_569_proj[3],x6G08_MD_569_proj[2],pch=20)
text(x6G08_MD_569_proj[3],x6G08_MD_569_proj[2],pos=4,label="569",col="green")

x6G08_MD_570_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_570.txt")
x6G08_MD_570_proj<-project.pca(x6G08_MD_570_raw,pcadata)
points(x6G08_MD_570_proj[3],x6G08_MD_570_proj[2],pch=20)
text(x6G08_MD_570_proj[3],x6G08_MD_570_proj[2],pos=4,label="570",col="green")

x6G08_MD_571_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_571.txt")
x6G08_MD_571_proj<-project.pca(x6G08_MD_571_raw,pcadata)
points(x6G08_MD_571_proj[3],x6G08_MD_571_proj[2],pch=20)
text(x6G08_MD_571_proj[3],x6G08_MD_571_proj[2],pos=4,label="571",col="red")

x6G08_MD_572_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_572.txt")
x6G08_MD_572_proj<-project.pca(x6G08_MD_572_raw,pcadata)
points(x6G08_MD_572_proj[3],x6G08_MD_572_proj[2],pch=20)
text(x6G08_MD_572_proj[3],x6G08_MD_572_proj[2],pos=4,label="572",col="blue")

x6G08_MD_573_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_573.txt")
x6G08_MD_573_proj<-project.pca(x6G08_MD_573_raw,pcadata)
points(x6G08_MD_573_proj[3],x6G08_MD_573_proj[2],pch=20)
text(x6G08_MD_573_proj[3],x6G08_MD_573_proj[2],pos=4,label="573",col="black")

x6G08_MD_574_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_574.txt")
x6G08_MD_574_proj<-project.pca(x6G08_MD_574_raw,pcadata)
points(x6G08_MD_574_proj[3],x6G08_MD_574_proj[2],pch=20)
text(x6G08_MD_574_proj[3],x6G08_MD_574_proj[2],pos=4,label="574",col="black")

x6G08_MD_575_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_575.txt")
x6G08_MD_575_proj<-project.pca(x6G08_MD_575_raw,pcadata)
points(x6G08_MD_575_proj[3],x6G08_MD_575_proj[2],pch=20)
text(x6G08_MD_575_proj[3],x6G08_MD_575_proj[2],pos=4,label="575",col="red")

x6G08_MD_576_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_576.txt")
x6G08_MD_576_proj<-project.pca(x6G08_MD_576_raw,pcadata)
points(x6G08_MD_576_proj[3],x6G08_MD_576_proj[2],pch=20)
text(x6G08_MD_576_proj[3],x6G08_MD_576_proj[2],pos=4,label="576",col="blue")

x6G08_MD_577_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_577.txt")
x6G08_MD_577_proj<-project.pca(x6G08_MD_577_raw,pcadata)
points(x6G08_MD_577_proj[3],x6G08_MD_577_proj[2],pch=20)
text(x6G08_MD_577_proj[3],x6G08_MD_577_proj[2],pos=4,label="577",col="black")

x6G08_MD_578_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_578.txt")
x6G08_MD_578_proj<-project.pca(x6G08_MD_578_raw,pcadata)
points(x6G08_MD_578_proj[3],x6G08_MD_578_proj[2],pch=20)
text(x6G08_MD_578_proj[3],x6G08_MD_578_proj[2],pos=4,label="578",col="blue")

x6G08_MD_579_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_579.txt")
x6G08_MD_579_proj<-project.pca(x6G08_MD_579_raw,pcadata)
points(x6G08_MD_579_proj[3],x6G08_MD_579_proj[2],pch=20)
text(x6G08_MD_579_proj[3],x6G08_MD_579_proj[2],pos=4,label="579",col="blue")

x6G08_MD_580_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_580.txt")
x6G08_MD_580_proj<-project.pca(x6G08_MD_580_raw,pcadata)
points(x6G08_MD_580_proj[3],x6G08_MD_580_proj[2],pch=20)
text(x6G08_MD_580_proj[3],x6G08_MD_580_proj[2],pos=4,label="580",col="green")

x6G08_MD_581_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_581.txt")
x6G08_MD_581_proj<-project.pca(x6G08_MD_581_raw,pcadata)
points(x6G08_MD_581_proj[3],x6G08_MD_581_proj[2],pch=20)
text(x6G08_MD_581_proj[3],x6G08_MD_581_proj[2],pos=4,label="581",col="green")

x6G08_MD_582_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_582.txt")
x6G08_MD_582_proj<-project.pca(x6G08_MD_582_raw,pcadata)
points(x6G08_MD_582_proj[3],x6G08_MD_582_proj[2],pch=20)
text(x6G08_MD_582_proj[3],x6G08_MD_582_proj[2],pos=4,label="582",col="green")

x6G08_MD_583_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_583.txt")
x6G08_MD_583_proj<-project.pca(x6G08_MD_583_raw,pcadata)
points(x6G08_MD_583_proj[3],x6G08_MD_583_proj[2],pch=20)
text(x6G08_MD_583_proj[3],x6G08_MD_583_proj[2],pos=4,label="583",col="blue")

x6G08_MD_584_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_584.txt")
x6G08_MD_584_proj<-project.pca(x6G08_MD_584_raw,pcadata)
points(x6G08_MD_584_proj[3],x6G08_MD_584_proj[2],pch=20)
text(x6G08_MD_584_proj[3],x6G08_MD_584_proj[2],pos=4,label="584",col="red")

x6G08_MD_585_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_585.txt")
x6G08_MD_585_proj<-project.pca(x6G08_MD_585_raw,pcadata)
points(x6G08_MD_585_proj[3],x6G08_MD_585_proj[2],pch=20)
text(x6G08_MD_585_proj[3],x6G08_MD_585_proj[2],pos=4,label="585",col="blue")

x6G08_MD_586_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_586.txt")
x6G08_MD_586_proj<-project.pca(x6G08_MD_586_raw,pcadata)
points(x6G08_MD_586_proj[3],x6G08_MD_586_proj[2],pch=20)
text(x6G08_MD_586_proj[3],x6G08_MD_586_proj[2],pos=4,label="586",col="blue")

x6G08_MD_587_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_587.txt")
x6G08_MD_587_proj<-project.pca(x6G08_MD_587_raw,pcadata)
points(x6G08_MD_587_proj[3],x6G08_MD_587_proj[2],pch=20)
text(x6G08_MD_587_proj[3],x6G08_MD_587_proj[2],pos=4,label="587",col="black")

x6G08_MD_588_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_588.txt")
x6G08_MD_588_proj<-project.pca(x6G08_MD_588_raw,pcadata)
points(x6G08_MD_588_proj[3],x6G08_MD_588_proj[2],pch=20)
text(x6G08_MD_588_proj[3],x6G08_MD_588_proj[2],pos=4,label="588",col="green")

x6G08_MD_589_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_589.txt")
x6G08_MD_589_proj<-project.pca(x6G08_MD_589_raw,pcadata)
points(x6G08_MD_589_proj[3],x6G08_MD_589_proj[2],pch=20)
text(x6G08_MD_589_proj[3],x6G08_MD_589_proj[2],pos=4,label="589",col="blue")

x6G08_MD_590_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_590.txt")
x6G08_MD_590_proj<-project.pca(x6G08_MD_590_raw,pcadata)
points(x6G08_MD_590_proj[3],x6G08_MD_590_proj[2],pch=20)
text(x6G08_MD_590_proj[3],x6G08_MD_590_proj[2],pos=4,label="590",col="black")

x6G08_MD_591_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_591.txt")
x6G08_MD_591_proj<-project.pca(x6G08_MD_591_raw,pcadata)
points(x6G08_MD_591_proj[3],x6G08_MD_591_proj[2],pch=20)
text(x6G08_MD_591_proj[3],x6G08_MD_591_proj[2],pos=4,label="591",col="blue")

x6G08_MD_592_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_592.txt")
x6G08_MD_592_proj<-project.pca(x6G08_MD_592_raw,pcadata)
points(x6G08_MD_592_proj[3],x6G08_MD_592_proj[2],pch=20)
text(x6G08_MD_592_proj[3],x6G08_MD_592_proj[2],pos=4,label="592",col="green")

x6G08_MD_593_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_593.txt")
x6G08_MD_593_proj<-project.pca(x6G08_MD_593_raw,pcadata)
points(x6G08_MD_593_proj[3],x6G08_MD_593_proj[2],pch=20)
text(x6G08_MD_593_proj[3],x6G08_MD_593_proj[2],pos=4,label="593",col="blue")

x6G08_MD_594_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_594.txt")
x6G08_MD_594_proj<-project.pca(x6G08_MD_594_raw,pcadata)
points(x6G08_MD_594_proj[3],x6G08_MD_594_proj[2],pch=20)
text(x6G08_MD_594_proj[3],x6G08_MD_594_proj[2],pos=4,label="594",col="blue")

x6G08_MD_595_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_595.txt")
x6G08_MD_595_proj<-project.pca(x6G08_MD_595_raw,pcadata)
points(x6G08_MD_595_proj[3],x6G08_MD_595_proj[2],pch=20)
text(x6G08_MD_595_proj[3],x6G08_MD_595_proj[2],pos=4,label="595",col="black")

x6G08_MD_596_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_596.txt")
x6G08_MD_596_proj<-project.pca(x6G08_MD_596_raw,pcadata)
points(x6G08_MD_596_proj[3],x6G08_MD_596_proj[2],pch=20)
text(x6G08_MD_596_proj[3],x6G08_MD_596_proj[2],pos=4,label="596",col="black")

x6G08_MD_597_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_597.txt")
x6G08_MD_597_proj<-project.pca(x6G08_MD_597_raw,pcadata)
points(x6G08_MD_597_proj[3],x6G08_MD_597_proj[2],pch=20)
text(x6G08_MD_597_proj[3],x6G08_MD_597_proj[2],pos=4,label="597",col="blue")

x6G08_MD_598_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_598.txt")
x6G08_MD_598_proj<-project.pca(x6G08_MD_598_raw,pcadata)
points(x6G08_MD_598_proj[3],x6G08_MD_598_proj[2],pch=20)
text(x6G08_MD_598_proj[3],x6G08_MD_598_proj[2],pos=4,label="598",col="blue")

x6G08_MD_599_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_599.txt")
x6G08_MD_599_proj<-project.pca(x6G08_MD_599_raw,pcadata)
points(x6G08_MD_599_proj[3],x6G08_MD_599_proj[2],pch=20)
text(x6G08_MD_599_proj[3],x6G08_MD_599_proj[2],pos=4,label="599",col="green")

x6G08_MD_600_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_600.txt")
x6G08_MD_600_proj<-project.pca(x6G08_MD_600_raw,pcadata)
points(x6G08_MD_600_proj[3],x6G08_MD_600_proj[2],pch=20)
text(x6G08_MD_600_proj[3],x6G08_MD_600_proj[2],pos=4,label="600",col="green")

x6G08_MD_601_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_601.txt")
x6G08_MD_601_proj<-project.pca(x6G08_MD_601_raw,pcadata)
points(x6G08_MD_601_proj[3],x6G08_MD_601_proj[2],pch=20)
text(x6G08_MD_601_proj[3],x6G08_MD_601_proj[2],pos=4,label="601",col="green")

x6G08_MD_602_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_602.txt")
x6G08_MD_602_proj<-project.pca(x6G08_MD_602_raw,pcadata)
points(x6G08_MD_602_proj[3],x6G08_MD_602_proj[2],pch=20)
text(x6G08_MD_602_proj[3],x6G08_MD_602_proj[2],pos=4,label="602",col="green")

x6G08_MD_603_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_603.txt")
x6G08_MD_603_proj<-project.pca(x6G08_MD_603_raw,pcadata)
points(x6G08_MD_603_proj[3],x6G08_MD_603_proj[2],pch=20)
text(x6G08_MD_603_proj[3],x6G08_MD_603_proj[2],pos=4,label="603",col="blue")

x6G08_MD_604_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_604.txt")
x6G08_MD_604_proj<-project.pca(x6G08_MD_604_raw,pcadata)
points(x6G08_MD_604_proj[3],x6G08_MD_604_proj[2],pch=20)
text(x6G08_MD_604_proj[3],x6G08_MD_604_proj[2],pos=4,label="604",col="black")

x6G08_MD_605_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_605.txt")
x6G08_MD_605_proj<-project.pca(x6G08_MD_605_raw,pcadata)
points(x6G08_MD_605_proj[3],x6G08_MD_605_proj[2],pch=20)
text(x6G08_MD_605_proj[3],x6G08_MD_605_proj[2],pos=4,label="605",col="blue")

x6G08_MD_606_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_606.txt")
x6G08_MD_606_proj<-project.pca(x6G08_MD_606_raw,pcadata)
points(x6G08_MD_606_proj[3],x6G08_MD_606_proj[2],pch=20)
text(x6G08_MD_606_proj[3],x6G08_MD_606_proj[2],pos=4,label="606",col="red")

x6G08_MD_607_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_607.txt")
x6G08_MD_607_proj<-project.pca(x6G08_MD_607_raw,pcadata)
points(x6G08_MD_607_proj[3],x6G08_MD_607_proj[2],pch=20)
text(x6G08_MD_607_proj[3],x6G08_MD_607_proj[2],pos=4,label="607",col="red")

x6G08_MD_608_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_608.txt")
x6G08_MD_608_proj<-project.pca(x6G08_MD_608_raw,pcadata)
points(x6G08_MD_608_proj[3],x6G08_MD_608_proj[2],pch=20)
text(x6G08_MD_608_proj[3],x6G08_MD_608_proj[2],pos=4,label="608",col="red")

x6G08_MD_609_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_609.txt")
x6G08_MD_609_proj<-project.pca(x6G08_MD_609_raw,pcadata)
points(x6G08_MD_609_proj[3],x6G08_MD_609_proj[2],pch=20)
text(x6G08_MD_609_proj[3],x6G08_MD_609_proj[2],pos=4,label="609",col="blue")

x6G08_MD_610_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_610.txt")
x6G08_MD_610_proj<-project.pca(x6G08_MD_610_raw,pcadata)
points(x6G08_MD_610_proj[3],x6G08_MD_610_proj[2],pch=20)
text(x6G08_MD_610_proj[3],x6G08_MD_610_proj[2],pos=4,label="610",col="black")

x6G08_MD_611_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_611.txt")
x6G08_MD_611_proj<-project.pca(x6G08_MD_611_raw,pcadata)
points(x6G08_MD_611_proj[3],x6G08_MD_611_proj[2],pch=20)
text(x6G08_MD_611_proj[3],x6G08_MD_611_proj[2],pos=4,label="611",col="blue")

x6G08_MD_612_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_612.txt")
x6G08_MD_612_proj<-project.pca(x6G08_MD_612_raw,pcadata)
points(x6G08_MD_612_proj[3],x6G08_MD_612_proj[2],pch=20)
text(x6G08_MD_612_proj[3],x6G08_MD_612_proj[2],pos=4,label="612",col="green")

x6G08_MD_613_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_613.txt")
x6G08_MD_613_proj<-project.pca(x6G08_MD_613_raw,pcadata)
points(x6G08_MD_613_proj[3],x6G08_MD_613_proj[2],pch=20)
text(x6G08_MD_613_proj[3],x6G08_MD_613_proj[2],pos=4,label="613",col="red")

x6G08_MD_614_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_614.txt")
x6G08_MD_614_proj<-project.pca(x6G08_MD_614_raw,pcadata)
points(x6G08_MD_614_proj[3],x6G08_MD_614_proj[2],pch=20)
text(x6G08_MD_614_proj[3],x6G08_MD_614_proj[2],pos=4,label="614",col="blue")

x6G08_MD_615_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_615.txt")
x6G08_MD_615_proj<-project.pca(x6G08_MD_615_raw,pcadata)
points(x6G08_MD_615_proj[3],x6G08_MD_615_proj[2],pch=20)
text(x6G08_MD_615_proj[3],x6G08_MD_615_proj[2],pos=4,label="615",col="blue")

x6G08_MD_616_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_616.txt")
x6G08_MD_616_proj<-project.pca(x6G08_MD_616_raw,pcadata)
points(x6G08_MD_616_proj[3],x6G08_MD_616_proj[2],pch=20)
text(x6G08_MD_616_proj[3],x6G08_MD_616_proj[2],pos=4,label="616",col="blue")

x6G08_MD_617_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_617.txt")
x6G08_MD_617_proj<-project.pca(x6G08_MD_617_raw,pcadata)
points(x6G08_MD_617_proj[3],x6G08_MD_617_proj[2],pch=20)
text(x6G08_MD_617_proj[3],x6G08_MD_617_proj[2],pos=4,label="617",col="green")

x6G08_MD_618_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_618.txt")
x6G08_MD_618_proj<-project.pca(x6G08_MD_618_raw,pcadata)
points(x6G08_MD_618_proj[3],x6G08_MD_618_proj[2],pch=20)
text(x6G08_MD_618_proj[3],x6G08_MD_618_proj[2],pos=4,label="618",col="green")

x6G08_MD_619_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_619.txt")
x6G08_MD_619_proj<-project.pca(x6G08_MD_619_raw,pcadata)
points(x6G08_MD_619_proj[3],x6G08_MD_619_proj[2],pch=20)
text(x6G08_MD_619_proj[3],x6G08_MD_619_proj[2],pos=4,label="619",col="blue")

x6G08_MD_620_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_620.txt")
x6G08_MD_620_proj<-project.pca(x6G08_MD_620_raw,pcadata)
points(x6G08_MD_620_proj[3],x6G08_MD_620_proj[2],pch=20)
text(x6G08_MD_620_proj[3],x6G08_MD_620_proj[2],pos=4,label="620",col="green")

x6G08_MD_621_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_621.txt")
x6G08_MD_621_proj<-project.pca(x6G08_MD_621_raw,pcadata)
points(x6G08_MD_621_proj[3],x6G08_MD_621_proj[2],pch=20)
text(x6G08_MD_621_proj[3],x6G08_MD_621_proj[2],pos=4,label="621",col="blue")

x6G08_MD_622_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_622.txt")
x6G08_MD_622_proj<-project.pca(x6G08_MD_622_raw,pcadata)
points(x6G08_MD_622_proj[3],x6G08_MD_622_proj[2],pch=20)
text(x6G08_MD_622_proj[3],x6G08_MD_622_proj[2],pos=4,label="622",col="blue")

x6G08_MD_623_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_623.txt")
x6G08_MD_623_proj<-project.pca(x6G08_MD_623_raw,pcadata)
points(x6G08_MD_623_proj[3],x6G08_MD_623_proj[2],pch=20)
text(x6G08_MD_623_proj[3],x6G08_MD_623_proj[2],pos=4,label="623",col="blue")

x6G08_MD_624_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_624.txt")
x6G08_MD_624_proj<-project.pca(x6G08_MD_624_raw,pcadata)
points(x6G08_MD_624_proj[3],x6G08_MD_624_proj[2],pch=20)
text(x6G08_MD_624_proj[3],x6G08_MD_624_proj[2],pos=4,label="624",col="blue")

x6G08_MD_625_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_625.txt")
x6G08_MD_625_proj<-project.pca(x6G08_MD_625_raw,pcadata)
points(x6G08_MD_625_proj[3],x6G08_MD_625_proj[2],pch=20)
text(x6G08_MD_625_proj[3],x6G08_MD_625_proj[2],pos=4,label="625",col="blue")

x6G08_MD_626_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_626.txt")
x6G08_MD_626_proj<-project.pca(x6G08_MD_626_raw,pcadata)
points(x6G08_MD_626_proj[3],x6G08_MD_626_proj[2],pch=20)
text(x6G08_MD_626_proj[3],x6G08_MD_626_proj[2],pos=4,label="626",col="blue")

x6G08_MD_627_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_627.txt")
x6G08_MD_627_proj<-project.pca(x6G08_MD_627_raw,pcadata)
points(x6G08_MD_627_proj[3],x6G08_MD_627_proj[2],pch=20)
text(x6G08_MD_627_proj[3],x6G08_MD_627_proj[2],pos=4,label="627",col="blue")

x6G08_MD_628_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_628.txt")
x6G08_MD_628_proj<-project.pca(x6G08_MD_628_raw,pcadata)
points(x6G08_MD_628_proj[3],x6G08_MD_628_proj[2],pch=20)
text(x6G08_MD_628_proj[3],x6G08_MD_628_proj[2],pos=4,label="628",col="red")

x6G08_MD_629_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_629.txt")
x6G08_MD_629_proj<-project.pca(x6G08_MD_629_raw,pcadata)
points(x6G08_MD_629_proj[3],x6G08_MD_629_proj[2],pch=20)
text(x6G08_MD_629_proj[3],x6G08_MD_629_proj[2],pos=4,label="629",col="red")

x6G08_MD_630_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_630.txt")
x6G08_MD_630_proj<-project.pca(x6G08_MD_630_raw,pcadata)
points(x6G08_MD_630_proj[3],x6G08_MD_630_proj[2],pch=20)
text(x6G08_MD_630_proj[3],x6G08_MD_630_proj[2],pos=4,label="630",col="blue")

x6G08_MD_631_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_631.txt")
x6G08_MD_631_proj<-project.pca(x6G08_MD_631_raw,pcadata)
points(x6G08_MD_631_proj[3],x6G08_MD_631_proj[2],pch=20)
text(x6G08_MD_631_proj[3],x6G08_MD_631_proj[2],pos=4,label="631",col="green")

x6G08_MD_632_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_632.txt")
x6G08_MD_632_proj<-project.pca(x6G08_MD_632_raw,pcadata)
points(x6G08_MD_632_proj[3],x6G08_MD_632_proj[2],pch=20)
text(x6G08_MD_632_proj[3],x6G08_MD_632_proj[2],pos=4,label="632",col="green")

x6G08_MD_633_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_633.txt")
x6G08_MD_633_proj<-project.pca(x6G08_MD_633_raw,pcadata)
points(x6G08_MD_633_proj[3],x6G08_MD_633_proj[2],pch=20)
text(x6G08_MD_633_proj[3],x6G08_MD_633_proj[2],pos=4,label="633",col="red")

x6G08_MD_634_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_634.txt")
x6G08_MD_634_proj<-project.pca(x6G08_MD_634_raw,pcadata)
points(x6G08_MD_634_proj[3],x6G08_MD_634_proj[2],pch=20)
text(x6G08_MD_634_proj[3],x6G08_MD_634_proj[2],pos=4,label="634",col="green")

x6G08_MD_635_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_635.txt")
x6G08_MD_635_proj<-project.pca(x6G08_MD_635_raw,pcadata)
points(x6G08_MD_635_proj[3],x6G08_MD_635_proj[2],pch=20)
text(x6G08_MD_635_proj[3],x6G08_MD_635_proj[2],pos=4,label="635",col="blue")

x6G08_MD_636_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_636.txt")
x6G08_MD_636_proj<-project.pca(x6G08_MD_636_raw,pcadata)
points(x6G08_MD_636_proj[3],x6G08_MD_636_proj[2],pch=20)
text(x6G08_MD_636_proj[3],x6G08_MD_636_proj[2],pos=4,label="636",col="black")

x6G08_MD_637_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_637.txt")
x6G08_MD_637_proj<-project.pca(x6G08_MD_637_raw,pcadata)
points(x6G08_MD_637_proj[3],x6G08_MD_637_proj[2],pch=20)
text(x6G08_MD_637_proj[3],x6G08_MD_637_proj[2],pos=4,label="637",col="blue")

x6G08_MD_638_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_638.txt")
x6G08_MD_638_proj<-project.pca(x6G08_MD_638_raw,pcadata)
points(x6G08_MD_638_proj[3],x6G08_MD_638_proj[2],pch=20)
text(x6G08_MD_638_proj[3],x6G08_MD_638_proj[2],pos=4,label="638",col="blue")

x6G08_MD_639_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_639.txt")
x6G08_MD_639_proj<-project.pca(x6G08_MD_639_raw,pcadata)
points(x6G08_MD_639_proj[3],x6G08_MD_639_proj[2],pch=20)
text(x6G08_MD_639_proj[3],x6G08_MD_639_proj[2],pos=4,label="639",col="blue")

x6G08_MD_640_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_640.txt")
x6G08_MD_640_proj<-project.pca(x6G08_MD_640_raw,pcadata)
points(x6G08_MD_640_proj[3],x6G08_MD_640_proj[2],pch=20)
text(x6G08_MD_640_proj[3],x6G08_MD_640_proj[2],pos=4,label="640",col="green")

x6G08_MD_641_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_641.txt")
x6G08_MD_641_proj<-project.pca(x6G08_MD_641_raw,pcadata)
points(x6G08_MD_641_proj[3],x6G08_MD_641_proj[2],pch=20)
text(x6G08_MD_641_proj[3],x6G08_MD_641_proj[2],pos=4,label="641",col="green")

x6G08_MD_642_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_642.txt")
x6G08_MD_642_proj<-project.pca(x6G08_MD_642_raw,pcadata)
points(x6G08_MD_642_proj[3],x6G08_MD_642_proj[2],pch=20)
text(x6G08_MD_642_proj[3],x6G08_MD_642_proj[2],pos=4,label="642",col="blue")

x6G08_MD_643_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_643.txt")
x6G08_MD_643_proj<-project.pca(x6G08_MD_643_raw,pcadata)
points(x6G08_MD_643_proj[3],x6G08_MD_643_proj[2],pch=20)
text(x6G08_MD_643_proj[3],x6G08_MD_643_proj[2],pos=4,label="643",col="blue")

x6G08_MD_644_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_644.txt")
x6G08_MD_644_proj<-project.pca(x6G08_MD_644_raw,pcadata)
points(x6G08_MD_644_proj[3],x6G08_MD_644_proj[2],pch=20)
text(x6G08_MD_644_proj[3],x6G08_MD_644_proj[2],pos=4,label="644",col="green")

x6G08_MD_645_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_645.txt")
x6G08_MD_645_proj<-project.pca(x6G08_MD_645_raw,pcadata)
points(x6G08_MD_645_proj[3],x6G08_MD_645_proj[2],pch=20)
text(x6G08_MD_645_proj[3],x6G08_MD_645_proj[2],pos=4,label="645",col="blue")

x6G08_MD_646_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_646.txt")
x6G08_MD_646_proj<-project.pca(x6G08_MD_646_raw,pcadata)
points(x6G08_MD_646_proj[3],x6G08_MD_646_proj[2],pch=20)
text(x6G08_MD_646_proj[3],x6G08_MD_646_proj[2],pos=4,label="646",col="blue")

x6G08_MD_647_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_647.txt")
x6G08_MD_647_proj<-project.pca(x6G08_MD_647_raw,pcadata)
points(x6G08_MD_647_proj[3],x6G08_MD_647_proj[2],pch=20)
text(x6G08_MD_647_proj[3],x6G08_MD_647_proj[2],pos=4,label="647",col="red")

x6G08_MD_648_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_648.txt")
x6G08_MD_648_proj<-project.pca(x6G08_MD_648_raw,pcadata)
points(x6G08_MD_648_proj[3],x6G08_MD_648_proj[2],pch=20)
text(x6G08_MD_648_proj[3],x6G08_MD_648_proj[2],pos=4,label="648",col="red")

x6G08_MD_649_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_649.txt")
x6G08_MD_649_proj<-project.pca(x6G08_MD_649_raw,pcadata)
points(x6G08_MD_649_proj[3],x6G08_MD_649_proj[2],pch=20)
text(x6G08_MD_649_proj[3],x6G08_MD_649_proj[2],pos=4,label="649",col="red")

x6G08_MD_650_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_650.txt")
x6G08_MD_650_proj<-project.pca(x6G08_MD_650_raw,pcadata)
points(x6G08_MD_650_proj[3],x6G08_MD_650_proj[2],pch=20)
text(x6G08_MD_650_proj[3],x6G08_MD_650_proj[2],pos=4,label="650",col="red")

x6G08_MD_651_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_651.txt")
x6G08_MD_651_proj<-project.pca(x6G08_MD_651_raw,pcadata)
points(x6G08_MD_651_proj[3],x6G08_MD_651_proj[2],pch=20)
text(x6G08_MD_651_proj[3],x6G08_MD_651_proj[2],pos=4,label="651",col="red")

x6G08_MD_652_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_652.txt")
x6G08_MD_652_proj<-project.pca(x6G08_MD_652_raw,pcadata)
points(x6G08_MD_652_proj[3],x6G08_MD_652_proj[2],pch=20)
text(x6G08_MD_652_proj[3],x6G08_MD_652_proj[2],pos=4,label="652",col="red")

x6G08_MD_653_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_653.txt")
x6G08_MD_653_proj<-project.pca(x6G08_MD_653_raw,pcadata)
points(x6G08_MD_653_proj[3],x6G08_MD_653_proj[2],pch=20)
text(x6G08_MD_653_proj[3],x6G08_MD_653_proj[2],pos=4,label="653",col="red")

x6G08_MD_654_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_654.txt")
x6G08_MD_654_proj<-project.pca(x6G08_MD_654_raw,pcadata)
points(x6G08_MD_654_proj[3],x6G08_MD_654_proj[2],pch=20)
text(x6G08_MD_654_proj[3],x6G08_MD_654_proj[2],pos=4,label="654",col="red")

x6G08_MD_655_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_655.txt")
x6G08_MD_655_proj<-project.pca(x6G08_MD_655_raw,pcadata)
points(x6G08_MD_655_proj[3],x6G08_MD_655_proj[2],pch=20)
text(x6G08_MD_655_proj[3],x6G08_MD_655_proj[2],pos=4,label="655",col="green")

x6G08_MD_656_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_656.txt")
x6G08_MD_656_proj<-project.pca(x6G08_MD_656_raw,pcadata)
points(x6G08_MD_656_proj[3],x6G08_MD_656_proj[2],pch=20)
text(x6G08_MD_656_proj[3],x6G08_MD_656_proj[2],pos=4,label="656",col="green")

x6G08_MD_657_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_657.txt")
x6G08_MD_657_proj<-project.pca(x6G08_MD_657_raw,pcadata)
points(x6G08_MD_657_proj[3],x6G08_MD_657_proj[2],pch=20)
text(x6G08_MD_657_proj[3],x6G08_MD_657_proj[2],pos=4,label="657",col="blue")

x6G08_MD_658_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_658.txt")
x6G08_MD_658_proj<-project.pca(x6G08_MD_658_raw,pcadata)
points(x6G08_MD_658_proj[3],x6G08_MD_658_proj[2],pch=20)
text(x6G08_MD_658_proj[3],x6G08_MD_658_proj[2],pos=4,label="658",col="green")

x6G08_MD_659_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_659.txt")
x6G08_MD_659_proj<-project.pca(x6G08_MD_659_raw,pcadata)
points(x6G08_MD_659_proj[3],x6G08_MD_659_proj[2],pch=20)
text(x6G08_MD_659_proj[3],x6G08_MD_659_proj[2],pos=4,label="659",col="green")

x6G08_MD_660_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_660.txt")
x6G08_MD_660_proj<-project.pca(x6G08_MD_660_raw,pcadata)
points(x6G08_MD_660_proj[3],x6G08_MD_660_proj[2],pch=20)
text(x6G08_MD_660_proj[3],x6G08_MD_660_proj[2],pos=4,label="660",col="red")

x6G08_MD_661_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_661.txt")
x6G08_MD_661_proj<-project.pca(x6G08_MD_661_raw,pcadata)
points(x6G08_MD_661_proj[3],x6G08_MD_661_proj[2],pch=20)
text(x6G08_MD_661_proj[3],x6G08_MD_661_proj[2],pos=4,label="661",col="blue")

x6G08_MD_662_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_662.txt")
x6G08_MD_662_proj<-project.pca(x6G08_MD_662_raw,pcadata)
points(x6G08_MD_662_proj[3],x6G08_MD_662_proj[2],pch=20)
text(x6G08_MD_662_proj[3],x6G08_MD_662_proj[2],pos=4,label="662",col="blue")

x6G08_MD_663_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_663.txt")
x6G08_MD_663_proj<-project.pca(x6G08_MD_663_raw,pcadata)
points(x6G08_MD_663_proj[3],x6G08_MD_663_proj[2],pch=20)
text(x6G08_MD_663_proj[3],x6G08_MD_663_proj[2],pos=4,label="663",col="green")

x6G08_MD_664_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_664.txt")
x6G08_MD_664_proj<-project.pca(x6G08_MD_664_raw,pcadata)
points(x6G08_MD_664_proj[3],x6G08_MD_664_proj[2],pch=20)
text(x6G08_MD_664_proj[3],x6G08_MD_664_proj[2],pos=4,label="664",col="red")

x6G08_MD_665_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_665.txt")
x6G08_MD_665_proj<-project.pca(x6G08_MD_665_raw,pcadata)
points(x6G08_MD_665_proj[3],x6G08_MD_665_proj[2],pch=20)
text(x6G08_MD_665_proj[3],x6G08_MD_665_proj[2],pos=4,label="665",col="red")

x6G08_MD_666_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_666.txt")
x6G08_MD_666_proj<-project.pca(x6G08_MD_666_raw,pcadata)
points(x6G08_MD_666_proj[3],x6G08_MD_666_proj[2],pch=20)
text(x6G08_MD_666_proj[3],x6G08_MD_666_proj[2],pos=4,label="666",col="red")

x6G08_MD_667_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_667.txt")
x6G08_MD_667_proj<-project.pca(x6G08_MD_667_raw,pcadata)
points(x6G08_MD_667_proj[3],x6G08_MD_667_proj[2],pch=20)
text(x6G08_MD_667_proj[3],x6G08_MD_667_proj[2],pos=4,label="667",col="red")

x6G08_MD_668_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_668.txt")
x6G08_MD_668_proj<-project.pca(x6G08_MD_668_raw,pcadata)
points(x6G08_MD_668_proj[3],x6G08_MD_668_proj[2],pch=20)
text(x6G08_MD_668_proj[3],x6G08_MD_668_proj[2],pos=4,label="668",col="red")

x6G08_MD_669_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_669.txt")
x6G08_MD_669_proj<-project.pca(x6G08_MD_669_raw,pcadata)
points(x6G08_MD_669_proj[3],x6G08_MD_669_proj[2],pch=20)
text(x6G08_MD_669_proj[3],x6G08_MD_669_proj[2],pos=4,label="669",col="red")

x6G08_MD_670_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_670.txt")
x6G08_MD_670_proj<-project.pca(x6G08_MD_670_raw,pcadata)
points(x6G08_MD_670_proj[3],x6G08_MD_670_proj[2],pch=20)
text(x6G08_MD_670_proj[3],x6G08_MD_670_proj[2],pos=4,label="670",col="red")

x6G08_MD_671_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_671.txt")
x6G08_MD_671_proj<-project.pca(x6G08_MD_671_raw,pcadata)
points(x6G08_MD_671_proj[3],x6G08_MD_671_proj[2],pch=20)
text(x6G08_MD_671_proj[3],x6G08_MD_671_proj[2],pos=4,label="671",col="red")

x6G08_MD_672_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_672.txt")
x6G08_MD_672_proj<-project.pca(x6G08_MD_672_raw,pcadata)
points(x6G08_MD_672_proj[3],x6G08_MD_672_proj[2],pch=20)
text(x6G08_MD_672_proj[3],x6G08_MD_672_proj[2],pos=4,label="672",col="red")

x6G08_MD_673_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_673.txt")
x6G08_MD_673_proj<-project.pca(x6G08_MD_673_raw,pcadata)
points(x6G08_MD_673_proj[3],x6G08_MD_673_proj[2],pch=20)
text(x6G08_MD_673_proj[3],x6G08_MD_673_proj[2],pos=4,label="673",col="red")

x6G08_MD_674_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_674.txt")
x6G08_MD_674_proj<-project.pca(x6G08_MD_674_raw,pcadata)
points(x6G08_MD_674_proj[3],x6G08_MD_674_proj[2],pch=20)
text(x6G08_MD_674_proj[3],x6G08_MD_674_proj[2],pos=4,label="674",col="red")

x6G08_MD_675_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_675.txt")
x6G08_MD_675_proj<-project.pca(x6G08_MD_675_raw,pcadata)
points(x6G08_MD_675_proj[3],x6G08_MD_675_proj[2],pch=20)
text(x6G08_MD_675_proj[3],x6G08_MD_675_proj[2],pos=4,label="675",col="red")

x6G08_MD_676_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_676.txt")
x6G08_MD_676_proj<-project.pca(x6G08_MD_676_raw,pcadata)
points(x6G08_MD_676_proj[3],x6G08_MD_676_proj[2],pch=20)
text(x6G08_MD_676_proj[3],x6G08_MD_676_proj[2],pos=4,label="676",col="red")

x6G08_MD_677_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_677.txt")
x6G08_MD_677_proj<-project.pca(x6G08_MD_677_raw,pcadata)
points(x6G08_MD_677_proj[3],x6G08_MD_677_proj[2],pch=20)
text(x6G08_MD_677_proj[3],x6G08_MD_677_proj[2],pos=4,label="677",col="red")

x6G08_MD_678_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_678.txt")
x6G08_MD_678_proj<-project.pca(x6G08_MD_678_raw,pcadata)
points(x6G08_MD_678_proj[3],x6G08_MD_678_proj[2],pch=20)
text(x6G08_MD_678_proj[3],x6G08_MD_678_proj[2],pos=4,label="678",col="red")

x6G08_MD_679_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_679.txt")
x6G08_MD_679_proj<-project.pca(x6G08_MD_679_raw,pcadata)
points(x6G08_MD_679_proj[3],x6G08_MD_679_proj[2],pch=20)
text(x6G08_MD_679_proj[3],x6G08_MD_679_proj[2],pos=4,label="679",col="red")

x6G08_MD_680_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_680.txt")
x6G08_MD_680_proj<-project.pca(x6G08_MD_680_raw,pcadata)
points(x6G08_MD_680_proj[3],x6G08_MD_680_proj[2],pch=20)
text(x6G08_MD_680_proj[3],x6G08_MD_680_proj[2],pos=4,label="680",col="red")

x6G08_MD_681_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_681.txt")
x6G08_MD_681_proj<-project.pca(x6G08_MD_681_raw,pcadata)
points(x6G08_MD_681_proj[3],x6G08_MD_681_proj[2],pch=20)
text(x6G08_MD_681_proj[3],x6G08_MD_681_proj[2],pos=4,label="681",col="red")

x6G08_MD_682_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_682.txt")
x6G08_MD_682_proj<-project.pca(x6G08_MD_682_raw,pcadata)
points(x6G08_MD_682_proj[3],x6G08_MD_682_proj[2],pch=20)
text(x6G08_MD_682_proj[3],x6G08_MD_682_proj[2],pos=4,label="682",col="red")

x6G08_MD_683_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_683.txt")
x6G08_MD_683_proj<-project.pca(x6G08_MD_683_raw,pcadata)
points(x6G08_MD_683_proj[3],x6G08_MD_683_proj[2],pch=20)
text(x6G08_MD_683_proj[3],x6G08_MD_683_proj[2],pos=4,label="683",col="red")

x6G08_MD_684_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_684.txt")
x6G08_MD_684_proj<-project.pca(x6G08_MD_684_raw,pcadata)
points(x6G08_MD_684_proj[3],x6G08_MD_684_proj[2],pch=20)
text(x6G08_MD_684_proj[3],x6G08_MD_684_proj[2],pos=4,label="684",col="red")

x6G08_MD_685_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_685.txt")
x6G08_MD_685_proj<-project.pca(x6G08_MD_685_raw,pcadata)
points(x6G08_MD_685_proj[3],x6G08_MD_685_proj[2],pch=20)
text(x6G08_MD_685_proj[3],x6G08_MD_685_proj[2],pos=4,label="685",col="red")

x6G08_MD_686_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_686.txt")
x6G08_MD_686_proj<-project.pca(x6G08_MD_686_raw,pcadata)
points(x6G08_MD_686_proj[3],x6G08_MD_686_proj[2],pch=20)
text(x6G08_MD_686_proj[3],x6G08_MD_686_proj[2],pos=4,label="686",col="red")

x6G08_MD_687_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_687.txt")
x6G08_MD_687_proj<-project.pca(x6G08_MD_687_raw,pcadata)
points(x6G08_MD_687_proj[3],x6G08_MD_687_proj[2],pch=20)
text(x6G08_MD_687_proj[3],x6G08_MD_687_proj[2],pos=4,label="687",col="red")

x6G08_MD_688_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_688.txt")
x6G08_MD_688_proj<-project.pca(x6G08_MD_688_raw,pcadata)
points(x6G08_MD_688_proj[3],x6G08_MD_688_proj[2],pch=20)
text(x6G08_MD_688_proj[3],x6G08_MD_688_proj[2],pos=4,label="688",col="red")

x6G08_MD_689_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_689.txt")
x6G08_MD_689_proj<-project.pca(x6G08_MD_689_raw,pcadata)
points(x6G08_MD_689_proj[3],x6G08_MD_689_proj[2],pch=20)
text(x6G08_MD_689_proj[3],x6G08_MD_689_proj[2],pos=4,label="689",col="red")

x6G08_MD_690_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_690.txt")
x6G08_MD_690_proj<-project.pca(x6G08_MD_690_raw,pcadata)
points(x6G08_MD_690_proj[3],x6G08_MD_690_proj[2],pch=20)
text(x6G08_MD_690_proj[3],x6G08_MD_690_proj[2],pos=4,label="690",col="red")

x6G08_MD_691_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_691.txt")
x6G08_MD_691_proj<-project.pca(x6G08_MD_691_raw,pcadata)
points(x6G08_MD_691_proj[3],x6G08_MD_691_proj[2],pch=20)
text(x6G08_MD_691_proj[3],x6G08_MD_691_proj[2],pos=4,label="691",col="red")

x6G08_MD_692_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_692.txt")
x6G08_MD_692_proj<-project.pca(x6G08_MD_692_raw,pcadata)
points(x6G08_MD_692_proj[3],x6G08_MD_692_proj[2],pch=20)
text(x6G08_MD_692_proj[3],x6G08_MD_692_proj[2],pos=4,label="692",col="green")

x6G08_MD_693_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_693.txt")
x6G08_MD_693_proj<-project.pca(x6G08_MD_693_raw,pcadata)
points(x6G08_MD_693_proj[3],x6G08_MD_693_proj[2],pch=20)
text(x6G08_MD_693_proj[3],x6G08_MD_693_proj[2],pos=4,label="693",col="red")

x6G08_MD_694_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_694.txt")
x6G08_MD_694_proj<-project.pca(x6G08_MD_694_raw,pcadata)
points(x6G08_MD_694_proj[3],x6G08_MD_694_proj[2],pch=20)
text(x6G08_MD_694_proj[3],x6G08_MD_694_proj[2],pos=4,label="694",col="red")

x6G08_MD_695_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_695.txt")
x6G08_MD_695_proj<-project.pca(x6G08_MD_695_raw,pcadata)
points(x6G08_MD_695_proj[3],x6G08_MD_695_proj[2],pch=20)
text(x6G08_MD_695_proj[3],x6G08_MD_695_proj[2],pos=4,label="695",col="red")

x6G08_MD_696_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_696.txt")
x6G08_MD_696_proj<-project.pca(x6G08_MD_696_raw,pcadata)
points(x6G08_MD_696_proj[3],x6G08_MD_696_proj[2],pch=20)
text(x6G08_MD_696_proj[3],x6G08_MD_696_proj[2],pos=4,label="696",col="red")

x6G08_MD_697_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_697.txt")
x6G08_MD_697_proj<-project.pca(x6G08_MD_697_raw,pcadata)
points(x6G08_MD_697_proj[3],x6G08_MD_697_proj[2],pch=20)
text(x6G08_MD_697_proj[3],x6G08_MD_697_proj[2],pos=4,label="697",col="red")

x6G08_MD_698_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_698.txt")
x6G08_MD_698_proj<-project.pca(x6G08_MD_698_raw,pcadata)
points(x6G08_MD_698_proj[3],x6G08_MD_698_proj[2],pch=20)
text(x6G08_MD_698_proj[3],x6G08_MD_698_proj[2],pos=4,label="698",col="red")

x6G08_MD_699_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_699.txt")
x6G08_MD_699_proj<-project.pca(x6G08_MD_699_raw,pcadata)
points(x6G08_MD_699_proj[3],x6G08_MD_699_proj[2],pch=20)
text(x6G08_MD_699_proj[3],x6G08_MD_699_proj[2],pos=4,label="699",col="red")

x6G08_MD_700_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_700.txt")
x6G08_MD_700_proj<-project.pca(x6G08_MD_700_raw,pcadata)
points(x6G08_MD_700_proj[3],x6G08_MD_700_proj[2],pch=20)
text(x6G08_MD_700_proj[3],x6G08_MD_700_proj[2],pos=4,label="700",col="red")

x6G08_MD_701_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_701.txt")
x6G08_MD_701_proj<-project.pca(x6G08_MD_701_raw,pcadata)
points(x6G08_MD_701_proj[3],x6G08_MD_701_proj[2],pch=20)
text(x6G08_MD_701_proj[3],x6G08_MD_701_proj[2],pos=4,label="701",col="red")

x6G08_MD_702_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_702.txt")
x6G08_MD_702_proj<-project.pca(x6G08_MD_702_raw,pcadata)
points(x6G08_MD_702_proj[3],x6G08_MD_702_proj[2],pch=20)
text(x6G08_MD_702_proj[3],x6G08_MD_702_proj[2],pos=4,label="702",col="green")

x6G08_MD_703_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_703.txt")
x6G08_MD_703_proj<-project.pca(x6G08_MD_703_raw,pcadata)
points(x6G08_MD_703_proj[3],x6G08_MD_703_proj[2],pch=20)
text(x6G08_MD_703_proj[3],x6G08_MD_703_proj[2],pos=4,label="703",col="red")

x6G08_MD_704_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_704.txt")
x6G08_MD_704_proj<-project.pca(x6G08_MD_704_raw,pcadata)
points(x6G08_MD_704_proj[3],x6G08_MD_704_proj[2],pch=20)
text(x6G08_MD_704_proj[3],x6G08_MD_704_proj[2],pos=4,label="704",col="red")

x6G08_MD_705_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_705.txt")
x6G08_MD_705_proj<-project.pca(x6G08_MD_705_raw,pcadata)
points(x6G08_MD_705_proj[3],x6G08_MD_705_proj[2],pch=20)
text(x6G08_MD_705_proj[3],x6G08_MD_705_proj[2],pos=4,label="705",col="red")

x6G08_MD_706_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_706.txt")
x6G08_MD_706_proj<-project.pca(x6G08_MD_706_raw,pcadata)
points(x6G08_MD_706_proj[3],x6G08_MD_706_proj[2],pch=20)
text(x6G08_MD_706_proj[3],x6G08_MD_706_proj[2],pos=4,label="706",col="red")

x6G08_MD_707_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_707.txt")
x6G08_MD_707_proj<-project.pca(x6G08_MD_707_raw,pcadata)
points(x6G08_MD_707_proj[3],x6G08_MD_707_proj[2],pch=20)
text(x6G08_MD_707_proj[3],x6G08_MD_707_proj[2],pos=4,label="707",col="red")

x6G08_MD_708_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_708.txt")
x6G08_MD_708_proj<-project.pca(x6G08_MD_708_raw,pcadata)
points(x6G08_MD_708_proj[3],x6G08_MD_708_proj[2],pch=20)
text(x6G08_MD_708_proj[3],x6G08_MD_708_proj[2],pos=4,label="708",col="red")

x6G08_MD_709_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_709.txt")
x6G08_MD_709_proj<-project.pca(x6G08_MD_709_raw,pcadata)
points(x6G08_MD_709_proj[3],x6G08_MD_709_proj[2],pch=20)
text(x6G08_MD_709_proj[3],x6G08_MD_709_proj[2],pos=4,label="709",col="red")

x6G08_MD_710_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_710.txt")
x6G08_MD_710_proj<-project.pca(x6G08_MD_710_raw,pcadata)
points(x6G08_MD_710_proj[3],x6G08_MD_710_proj[2],pch=20)
text(x6G08_MD_710_proj[3],x6G08_MD_710_proj[2],pos=4,label="710",col="red")

x6G08_MD_711_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_711.txt")
x6G08_MD_711_proj<-project.pca(x6G08_MD_711_raw,pcadata)
points(x6G08_MD_711_proj[3],x6G08_MD_711_proj[2],pch=20)
text(x6G08_MD_711_proj[3],x6G08_MD_711_proj[2],pos=4,label="711",col="red")

x6G08_MD_712_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_712.txt")
x6G08_MD_712_proj<-project.pca(x6G08_MD_712_raw,pcadata)
points(x6G08_MD_712_proj[3],x6G08_MD_712_proj[2],pch=20)
text(x6G08_MD_712_proj[3],x6G08_MD_712_proj[2],pos=4,label="712",col="red")

x6G08_MD_713_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_713.txt")
x6G08_MD_713_proj<-project.pca(x6G08_MD_713_raw,pcadata)
points(x6G08_MD_713_proj[3],x6G08_MD_713_proj[2],pch=20)
text(x6G08_MD_713_proj[3],x6G08_MD_713_proj[2],pos=4,label="713",col="red")

x6G08_MD_714_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_714.txt")
x6G08_MD_714_proj<-project.pca(x6G08_MD_714_raw,pcadata)
points(x6G08_MD_714_proj[3],x6G08_MD_714_proj[2],pch=20)
text(x6G08_MD_714_proj[3],x6G08_MD_714_proj[2],pos=4,label="714",col="red")

x6G08_MD_715_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_715.txt")
x6G08_MD_715_proj<-project.pca(x6G08_MD_715_raw,pcadata)
points(x6G08_MD_715_proj[3],x6G08_MD_715_proj[2],pch=20)
text(x6G08_MD_715_proj[3],x6G08_MD_715_proj[2],pos=4,label="715",col="red")

x6G08_MD_716_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_716.txt")
x6G08_MD_716_proj<-project.pca(x6G08_MD_716_raw,pcadata)
points(x6G08_MD_716_proj[3],x6G08_MD_716_proj[2],pch=20)
text(x6G08_MD_716_proj[3],x6G08_MD_716_proj[2],pos=4,label="716",col="red")

x6G08_MD_717_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_717.txt")
x6G08_MD_717_proj<-project.pca(x6G08_MD_717_raw,pcadata)
points(x6G08_MD_717_proj[3],x6G08_MD_717_proj[2],pch=20)
text(x6G08_MD_717_proj[3],x6G08_MD_717_proj[2],pos=4,label="717",col="red")

x6G08_MD_718_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_718.txt")
x6G08_MD_718_proj<-project.pca(x6G08_MD_718_raw,pcadata)
points(x6G08_MD_718_proj[3],x6G08_MD_718_proj[2],pch=20)
text(x6G08_MD_718_proj[3],x6G08_MD_718_proj[2],pos=4,label="718",col="red")

x6G08_MD_719_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_719.txt")
x6G08_MD_719_proj<-project.pca(x6G08_MD_719_raw,pcadata)
points(x6G08_MD_719_proj[3],x6G08_MD_719_proj[2],pch=20)
text(x6G08_MD_719_proj[3],x6G08_MD_719_proj[2],pos=4,label="719",col="red")

x6G08_MD_720_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_720.txt")
x6G08_MD_720_proj<-project.pca(x6G08_MD_720_raw,pcadata)
points(x6G08_MD_720_proj[3],x6G08_MD_720_proj[2],pch=20)
text(x6G08_MD_720_proj[3],x6G08_MD_720_proj[2],pos=4,label="720",col="red")

x6G08_MD_721_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_721.txt")
x6G08_MD_721_proj<-project.pca(x6G08_MD_721_raw,pcadata)
points(x6G08_MD_721_proj[3],x6G08_MD_721_proj[2],pch=20)
text(x6G08_MD_721_proj[3],x6G08_MD_721_proj[2],pos=4,label="721",col="red")

x6G08_MD_722_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_722.txt")
x6G08_MD_722_proj<-project.pca(x6G08_MD_722_raw,pcadata)
points(x6G08_MD_722_proj[3],x6G08_MD_722_proj[2],pch=20)
text(x6G08_MD_722_proj[3],x6G08_MD_722_proj[2],pos=4,label="722",col="red")

x6G08_MD_723_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_723.txt")
x6G08_MD_723_proj<-project.pca(x6G08_MD_723_raw,pcadata)
points(x6G08_MD_723_proj[3],x6G08_MD_723_proj[2],pch=20)
text(x6G08_MD_723_proj[3],x6G08_MD_723_proj[2],pos=4,label="723",col="red")

x6G08_MD_724_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_724.txt")
x6G08_MD_724_proj<-project.pca(x6G08_MD_724_raw,pcadata)
points(x6G08_MD_724_proj[3],x6G08_MD_724_proj[2],pch=20)
text(x6G08_MD_724_proj[3],x6G08_MD_724_proj[2],pos=4,label="724",col="red")

x6G08_MD_725_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_725.txt")
x6G08_MD_725_proj<-project.pca(x6G08_MD_725_raw,pcadata)
points(x6G08_MD_725_proj[3],x6G08_MD_725_proj[2],pch=20)
text(x6G08_MD_725_proj[3],x6G08_MD_725_proj[2],pos=4,label="725",col="red")

x6G08_MD_726_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_726.txt")
x6G08_MD_726_proj<-project.pca(x6G08_MD_726_raw,pcadata)
points(x6G08_MD_726_proj[3],x6G08_MD_726_proj[2],pch=20)
text(x6G08_MD_726_proj[3],x6G08_MD_726_proj[2],pos=4,label="726",col="red")

x6G08_MD_727_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_727.txt")
x6G08_MD_727_proj<-project.pca(x6G08_MD_727_raw,pcadata)
points(x6G08_MD_727_proj[3],x6G08_MD_727_proj[2],pch=20)
text(x6G08_MD_727_proj[3],x6G08_MD_727_proj[2],pos=4,label="727",col="red")

x6G08_MD_728_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_728.txt")
x6G08_MD_728_proj<-project.pca(x6G08_MD_728_raw,pcadata)
points(x6G08_MD_728_proj[3],x6G08_MD_728_proj[2],pch=20)
text(x6G08_MD_728_proj[3],x6G08_MD_728_proj[2],pos=4,label="728",col="green")

x6G08_MD_729_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_729.txt")
x6G08_MD_729_proj<-project.pca(x6G08_MD_729_raw,pcadata)
points(x6G08_MD_729_proj[3],x6G08_MD_729_proj[2],pch=20)
text(x6G08_MD_729_proj[3],x6G08_MD_729_proj[2],pos=4,label="729",col="red")

x6G08_MD_730_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_730.txt")
x6G08_MD_730_proj<-project.pca(x6G08_MD_730_raw,pcadata)
points(x6G08_MD_730_proj[3],x6G08_MD_730_proj[2],pch=20)
text(x6G08_MD_730_proj[3],x6G08_MD_730_proj[2],pos=4,label="730",col="red")

x6G08_MD_731_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_731.txt")
x6G08_MD_731_proj<-project.pca(x6G08_MD_731_raw,pcadata)
points(x6G08_MD_731_proj[3],x6G08_MD_731_proj[2],pch=20)
text(x6G08_MD_731_proj[3],x6G08_MD_731_proj[2],pos=4,label="731",col="red")

x6G08_MD_732_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_732.txt")
x6G08_MD_732_proj<-project.pca(x6G08_MD_732_raw,pcadata)
points(x6G08_MD_732_proj[3],x6G08_MD_732_proj[2],pch=20)
text(x6G08_MD_732_proj[3],x6G08_MD_732_proj[2],pos=4,label="732",col="red")

x6G08_MD_733_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_733.txt")
x6G08_MD_733_proj<-project.pca(x6G08_MD_733_raw,pcadata)
points(x6G08_MD_733_proj[3],x6G08_MD_733_proj[2],pch=20)
text(x6G08_MD_733_proj[3],x6G08_MD_733_proj[2],pos=4,label="733",col="red")

x6G08_MD_734_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_734.txt")
x6G08_MD_734_proj<-project.pca(x6G08_MD_734_raw,pcadata)
points(x6G08_MD_734_proj[3],x6G08_MD_734_proj[2],pch=20)
text(x6G08_MD_734_proj[3],x6G08_MD_734_proj[2],pos=4,label="734",col="red")

x6G08_MD_735_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_735.txt")
x6G08_MD_735_proj<-project.pca(x6G08_MD_735_raw,pcadata)
points(x6G08_MD_735_proj[3],x6G08_MD_735_proj[2],pch=20)
text(x6G08_MD_735_proj[3],x6G08_MD_735_proj[2],pos=4,label="735",col="red")

x6G08_MD_736_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_736.txt")
x6G08_MD_736_proj<-project.pca(x6G08_MD_736_raw,pcadata)
points(x6G08_MD_736_proj[3],x6G08_MD_736_proj[2],pch=20)
text(x6G08_MD_736_proj[3],x6G08_MD_736_proj[2],pos=4,label="736",col="red")

x6G08_MD_737_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_737.txt")
x6G08_MD_737_proj<-project.pca(x6G08_MD_737_raw,pcadata)
points(x6G08_MD_737_proj[3],x6G08_MD_737_proj[2],pch=20)
text(x6G08_MD_737_proj[3],x6G08_MD_737_proj[2],pos=4,label="737",col="red")

x6G08_MD_738_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_738.txt")
x6G08_MD_738_proj<-project.pca(x6G08_MD_738_raw,pcadata)
points(x6G08_MD_738_proj[3],x6G08_MD_738_proj[2],pch=20)
text(x6G08_MD_738_proj[3],x6G08_MD_738_proj[2],pos=4,label="738",col="red")

x6G08_MD_739_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_739.txt")
x6G08_MD_739_proj<-project.pca(x6G08_MD_739_raw,pcadata)
points(x6G08_MD_739_proj[3],x6G08_MD_739_proj[2],pch=20)
text(x6G08_MD_739_proj[3],x6G08_MD_739_proj[2],pos=4,label="739",col="red")

x6G08_MD_740_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_740.txt")
x6G08_MD_740_proj<-project.pca(x6G08_MD_740_raw,pcadata)
points(x6G08_MD_740_proj[3],x6G08_MD_740_proj[2],pch=20)
text(x6G08_MD_740_proj[3],x6G08_MD_740_proj[2],pos=4,label="740",col="red")

x6G08_MD_741_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_741.txt")
x6G08_MD_741_proj<-project.pca(x6G08_MD_741_raw,pcadata)
points(x6G08_MD_741_proj[3],x6G08_MD_741_proj[2],pch=20)
text(x6G08_MD_741_proj[3],x6G08_MD_741_proj[2],pos=4,label="741",col="red")

x6G08_MD_742_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_742.txt")
x6G08_MD_742_proj<-project.pca(x6G08_MD_742_raw,pcadata)
points(x6G08_MD_742_proj[3],x6G08_MD_742_proj[2],pch=20)
text(x6G08_MD_742_proj[3],x6G08_MD_742_proj[2],pos=4,label="742",col="red")

x6G08_MD_743_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_743.txt")
x6G08_MD_743_proj<-project.pca(x6G08_MD_743_raw,pcadata)
points(x6G08_MD_743_proj[3],x6G08_MD_743_proj[2],pch=20)
text(x6G08_MD_743_proj[3],x6G08_MD_743_proj[2],pos=4,label="743",col="red")

x6G08_MD_744_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_744.txt")
x6G08_MD_744_proj<-project.pca(x6G08_MD_744_raw,pcadata)
points(x6G08_MD_744_proj[3],x6G08_MD_744_proj[2],pch=20)
text(x6G08_MD_744_proj[3],x6G08_MD_744_proj[2],pos=4,label="744",col="red")

x6G08_MD_745_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_745.txt")
x6G08_MD_745_proj<-project.pca(x6G08_MD_745_raw,pcadata)
points(x6G08_MD_745_proj[3],x6G08_MD_745_proj[2],pch=20)
text(x6G08_MD_745_proj[3],x6G08_MD_745_proj[2],pos=4,label="745",col="red")

x6G08_MD_746_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_746.txt")
x6G08_MD_746_proj<-project.pca(x6G08_MD_746_raw,pcadata)
points(x6G08_MD_746_proj[3],x6G08_MD_746_proj[2],pch=20)
text(x6G08_MD_746_proj[3],x6G08_MD_746_proj[2],pos=4,label="746",col="red")

x6G08_MD_747_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_747.txt")
x6G08_MD_747_proj<-project.pca(x6G08_MD_747_raw,pcadata)
points(x6G08_MD_747_proj[3],x6G08_MD_747_proj[2],pch=20)
text(x6G08_MD_747_proj[3],x6G08_MD_747_proj[2],pos=4,label="747",col="red")

x6G08_MD_748_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_748.txt")
x6G08_MD_748_proj<-project.pca(x6G08_MD_748_raw,pcadata)
points(x6G08_MD_748_proj[3],x6G08_MD_748_proj[2],pch=20)
text(x6G08_MD_748_proj[3],x6G08_MD_748_proj[2],pos=4,label="748",col="blue")

x6G08_MD_749_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_749.txt")
x6G08_MD_749_proj<-project.pca(x6G08_MD_749_raw,pcadata)
points(x6G08_MD_749_proj[3],x6G08_MD_749_proj[2],pch=20)
text(x6G08_MD_749_proj[3],x6G08_MD_749_proj[2],pos=4,label="749",col="red")

x6G08_MD_750_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_750.txt")
x6G08_MD_750_proj<-project.pca(x6G08_MD_750_raw,pcadata)
points(x6G08_MD_750_proj[3],x6G08_MD_750_proj[2],pch=20)
text(x6G08_MD_750_proj[3],x6G08_MD_750_proj[2],pos=4,label="750",col="red")

x6G08_MD_751_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_751.txt")
x6G08_MD_751_proj<-project.pca(x6G08_MD_751_raw,pcadata)
points(x6G08_MD_751_proj[3],x6G08_MD_751_proj[2],pch=20)
text(x6G08_MD_751_proj[3],x6G08_MD_751_proj[2],pos=4,label="751",col="red")

x6G08_MD_752_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_752.txt")
x6G08_MD_752_proj<-project.pca(x6G08_MD_752_raw,pcadata)
points(x6G08_MD_752_proj[3],x6G08_MD_752_proj[2],pch=20)
text(x6G08_MD_752_proj[3],x6G08_MD_752_proj[2],pos=4,label="752",col="red")

x6G08_MD_753_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_753.txt")
x6G08_MD_753_proj<-project.pca(x6G08_MD_753_raw,pcadata)
points(x6G08_MD_753_proj[3],x6G08_MD_753_proj[2],pch=20)
text(x6G08_MD_753_proj[3],x6G08_MD_753_proj[2],pos=4,label="753",col="red")

x6G08_MD_754_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_754.txt")
x6G08_MD_754_proj<-project.pca(x6G08_MD_754_raw,pcadata)
points(x6G08_MD_754_proj[3],x6G08_MD_754_proj[2],pch=20)
text(x6G08_MD_754_proj[3],x6G08_MD_754_proj[2],pos=4,label="754",col="green")

x6G08_MD_755_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_755.txt")
x6G08_MD_755_proj<-project.pca(x6G08_MD_755_raw,pcadata)
points(x6G08_MD_755_proj[3],x6G08_MD_755_proj[2],pch=20)
text(x6G08_MD_755_proj[3],x6G08_MD_755_proj[2],pos=4,label="755",col="red")

x6G08_MD_756_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_756.txt")
x6G08_MD_756_proj<-project.pca(x6G08_MD_756_raw,pcadata)
points(x6G08_MD_756_proj[3],x6G08_MD_756_proj[2],pch=20)
text(x6G08_MD_756_proj[3],x6G08_MD_756_proj[2],pos=4,label="756",col="red")

x6G08_MD_757_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_757.txt")
x6G08_MD_757_proj<-project.pca(x6G08_MD_757_raw,pcadata)
points(x6G08_MD_757_proj[3],x6G08_MD_757_proj[2],pch=20)
text(x6G08_MD_757_proj[3],x6G08_MD_757_proj[2],pos=4,label="757",col="red")

x6G08_MD_758_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_758.txt")
x6G08_MD_758_proj<-project.pca(x6G08_MD_758_raw,pcadata)
points(x6G08_MD_758_proj[3],x6G08_MD_758_proj[2],pch=20)
text(x6G08_MD_758_proj[3],x6G08_MD_758_proj[2],pos=4,label="758",col="red")

x6G08_MD_759_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_759.txt")
x6G08_MD_759_proj<-project.pca(x6G08_MD_759_raw,pcadata)
points(x6G08_MD_759_proj[3],x6G08_MD_759_proj[2],pch=20)
text(x6G08_MD_759_proj[3],x6G08_MD_759_proj[2],pos=4,label="759",col="red")

x6G08_MD_760_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_760.txt")
x6G08_MD_760_proj<-project.pca(x6G08_MD_760_raw,pcadata)
points(x6G08_MD_760_proj[3],x6G08_MD_760_proj[2],pch=20)
text(x6G08_MD_760_proj[3],x6G08_MD_760_proj[2],pos=4,label="760",col="red")

x6G08_MD_761_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_761.txt")
x6G08_MD_761_proj<-project.pca(x6G08_MD_761_raw,pcadata)
points(x6G08_MD_761_proj[3],x6G08_MD_761_proj[2],pch=20)
text(x6G08_MD_761_proj[3],x6G08_MD_761_proj[2],pos=4,label="761",col="red")

x6G08_MD_762_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_762.txt")
x6G08_MD_762_proj<-project.pca(x6G08_MD_762_raw,pcadata)
points(x6G08_MD_762_proj[3],x6G08_MD_762_proj[2],pch=20)
text(x6G08_MD_762_proj[3],x6G08_MD_762_proj[2],pos=4,label="762",col="red")

x6G08_MD_763_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_763.txt")
x6G08_MD_763_proj<-project.pca(x6G08_MD_763_raw,pcadata)
points(x6G08_MD_763_proj[3],x6G08_MD_763_proj[2],pch=20)
text(x6G08_MD_763_proj[3],x6G08_MD_763_proj[2],pos=4,label="763",col="red")

x6G08_MD_764_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_764.txt")
x6G08_MD_764_proj<-project.pca(x6G08_MD_764_raw,pcadata)
points(x6G08_MD_764_proj[3],x6G08_MD_764_proj[2],pch=20)
text(x6G08_MD_764_proj[3],x6G08_MD_764_proj[2],pos=4,label="764",col="green")

x6G08_MD_765_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_765.txt")
x6G08_MD_765_proj<-project.pca(x6G08_MD_765_raw,pcadata)
points(x6G08_MD_765_proj[3],x6G08_MD_765_proj[2],pch=20)
text(x6G08_MD_765_proj[3],x6G08_MD_765_proj[2],pos=4,label="765",col="red")

x6G08_MD_766_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_766.txt")
x6G08_MD_766_proj<-project.pca(x6G08_MD_766_raw,pcadata)
points(x6G08_MD_766_proj[3],x6G08_MD_766_proj[2],pch=20)
text(x6G08_MD_766_proj[3],x6G08_MD_766_proj[2],pos=4,label="766",col="red")

x6G08_MD_767_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_767.txt")
x6G08_MD_767_proj<-project.pca(x6G08_MD_767_raw,pcadata)
points(x6G08_MD_767_proj[3],x6G08_MD_767_proj[2],pch=20)
text(x6G08_MD_767_proj[3],x6G08_MD_767_proj[2],pos=4,label="767",col="red")

x6G08_MD_768_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_768.txt")
x6G08_MD_768_proj<-project.pca(x6G08_MD_768_raw,pcadata)
points(x6G08_MD_768_proj[3],x6G08_MD_768_proj[2],pch=20)
text(x6G08_MD_768_proj[3],x6G08_MD_768_proj[2],pos=4,label="768",col="red")

x6G08_MD_769_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_769.txt")
x6G08_MD_769_proj<-project.pca(x6G08_MD_769_raw,pcadata)
points(x6G08_MD_769_proj[3],x6G08_MD_769_proj[2],pch=20)
text(x6G08_MD_769_proj[3],x6G08_MD_769_proj[2],pos=4,label="769",col="red")

x6G08_MD_770_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_770.txt")
x6G08_MD_770_proj<-project.pca(x6G08_MD_770_raw,pcadata)
points(x6G08_MD_770_proj[3],x6G08_MD_770_proj[2],pch=20)
text(x6G08_MD_770_proj[3],x6G08_MD_770_proj[2],pos=4,label="770",col="red")

x6G08_MD_771_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_771.txt")
x6G08_MD_771_proj<-project.pca(x6G08_MD_771_raw,pcadata)
points(x6G08_MD_771_proj[3],x6G08_MD_771_proj[2],pch=20)
text(x6G08_MD_771_proj[3],x6G08_MD_771_proj[2],pos=4,label="771",col="red")

x6G08_MD_772_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_772.txt")
x6G08_MD_772_proj<-project.pca(x6G08_MD_772_raw,pcadata)
points(x6G08_MD_772_proj[3],x6G08_MD_772_proj[2],pch=20)
text(x6G08_MD_772_proj[3],x6G08_MD_772_proj[2],pos=4,label="772",col="red")

x6G08_MD_773_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_773.txt")
x6G08_MD_773_proj<-project.pca(x6G08_MD_773_raw,pcadata)
points(x6G08_MD_773_proj[3],x6G08_MD_773_proj[2],pch=20)
text(x6G08_MD_773_proj[3],x6G08_MD_773_proj[2],pos=4,label="773",col="red")

x6G08_MD_774_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_774.txt")
x6G08_MD_774_proj<-project.pca(x6G08_MD_774_raw,pcadata)
points(x6G08_MD_774_proj[3],x6G08_MD_774_proj[2],pch=20)
text(x6G08_MD_774_proj[3],x6G08_MD_774_proj[2],pos=4,label="774",col="red")

x6G08_MD_775_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_775.txt")
x6G08_MD_775_proj<-project.pca(x6G08_MD_775_raw,pcadata)
points(x6G08_MD_775_proj[3],x6G08_MD_775_proj[2],pch=20)
text(x6G08_MD_775_proj[3],x6G08_MD_775_proj[2],pos=4,label="775",col="red")

x6G08_MD_776_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_776.txt")
x6G08_MD_776_proj<-project.pca(x6G08_MD_776_raw,pcadata)
points(x6G08_MD_776_proj[3],x6G08_MD_776_proj[2],pch=20)
text(x6G08_MD_776_proj[3],x6G08_MD_776_proj[2],pos=4,label="776",col="red")

x6G08_MD_777_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_777.txt")
x6G08_MD_777_proj<-project.pca(x6G08_MD_777_raw,pcadata)
points(x6G08_MD_777_proj[3],x6G08_MD_777_proj[2],pch=20)
text(x6G08_MD_777_proj[3],x6G08_MD_777_proj[2],pos=4,label="777",col="green")

x6G08_MD_778_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_778.txt")
x6G08_MD_778_proj<-project.pca(x6G08_MD_778_raw,pcadata)
points(x6G08_MD_778_proj[3],x6G08_MD_778_proj[2],pch=20)
text(x6G08_MD_778_proj[3],x6G08_MD_778_proj[2],pos=4,label="778",col="green")

x6G08_MD_779_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_779.txt")
x6G08_MD_779_proj<-project.pca(x6G08_MD_779_raw,pcadata)
points(x6G08_MD_779_proj[3],x6G08_MD_779_proj[2],pch=20)
text(x6G08_MD_779_proj[3],x6G08_MD_779_proj[2],pos=4,label="779",col="green")

x6G08_MD_780_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_780.txt")
x6G08_MD_780_proj<-project.pca(x6G08_MD_780_raw,pcadata)
points(x6G08_MD_780_proj[3],x6G08_MD_780_proj[2],pch=20)
text(x6G08_MD_780_proj[3],x6G08_MD_780_proj[2],pos=4,label="780",col="red")

x6G08_MD_781_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_781.txt")
x6G08_MD_781_proj<-project.pca(x6G08_MD_781_raw,pcadata)
points(x6G08_MD_781_proj[3],x6G08_MD_781_proj[2],pch=20)
text(x6G08_MD_781_proj[3],x6G08_MD_781_proj[2],pos=4,label="781",col="red")

x6G08_MD_782_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_782.txt")
x6G08_MD_782_proj<-project.pca(x6G08_MD_782_raw,pcadata)
points(x6G08_MD_782_proj[3],x6G08_MD_782_proj[2],pch=20)
text(x6G08_MD_782_proj[3],x6G08_MD_782_proj[2],pos=4,label="782",col="red")

x6G08_MD_783_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_783.txt")
x6G08_MD_783_proj<-project.pca(x6G08_MD_783_raw,pcadata)
points(x6G08_MD_783_proj[3],x6G08_MD_783_proj[2],pch=20)
text(x6G08_MD_783_proj[3],x6G08_MD_783_proj[2],pos=4,label="783",col="green")

x6G08_MD_784_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_784.txt")
x6G08_MD_784_proj<-project.pca(x6G08_MD_784_raw,pcadata)
points(x6G08_MD_784_proj[3],x6G08_MD_784_proj[2],pch=20)
text(x6G08_MD_784_proj[3],x6G08_MD_784_proj[2],pos=4,label="784",col="red")

x6G08_MD_785_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_785.txt")
x6G08_MD_785_proj<-project.pca(x6G08_MD_785_raw,pcadata)
points(x6G08_MD_785_proj[3],x6G08_MD_785_proj[2],pch=20)
text(x6G08_MD_785_proj[3],x6G08_MD_785_proj[2],pos=4,label="785",col="red")

x6G08_MD_786_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_786.txt")
x6G08_MD_786_proj<-project.pca(x6G08_MD_786_raw,pcadata)
points(x6G08_MD_786_proj[3],x6G08_MD_786_proj[2],pch=20)
text(x6G08_MD_786_proj[3],x6G08_MD_786_proj[2],pos=4,label="786",col="red")

x6G08_MD_787_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_787.txt")
x6G08_MD_787_proj<-project.pca(x6G08_MD_787_raw,pcadata)
points(x6G08_MD_787_proj[3],x6G08_MD_787_proj[2],pch=20)
text(x6G08_MD_787_proj[3],x6G08_MD_787_proj[2],pos=4,label="787",col="red")

x6G08_MD_788_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_788.txt")
x6G08_MD_788_proj<-project.pca(x6G08_MD_788_raw,pcadata)
points(x6G08_MD_788_proj[3],x6G08_MD_788_proj[2],pch=20)
text(x6G08_MD_788_proj[3],x6G08_MD_788_proj[2],pos=4,label="788",col="red")

x6G08_MD_789_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_789.txt")
x6G08_MD_789_proj<-project.pca(x6G08_MD_789_raw,pcadata)
points(x6G08_MD_789_proj[3],x6G08_MD_789_proj[2],pch=20)
text(x6G08_MD_789_proj[3],x6G08_MD_789_proj[2],pos=4,label="789",col="red")

x6G08_MD_790_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_790.txt")
x6G08_MD_790_proj<-project.pca(x6G08_MD_790_raw,pcadata)
points(x6G08_MD_790_proj[3],x6G08_MD_790_proj[2],pch=20)
text(x6G08_MD_790_proj[3],x6G08_MD_790_proj[2],pos=4,label="790",col="red")

x6G08_MD_791_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_791.txt")
x6G08_MD_791_proj<-project.pca(x6G08_MD_791_raw,pcadata)
points(x6G08_MD_791_proj[3],x6G08_MD_791_proj[2],pch=20)
text(x6G08_MD_791_proj[3],x6G08_MD_791_proj[2],pos=4,label="791",col="red")

x6G08_MD_792_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_792.txt")
x6G08_MD_792_proj<-project.pca(x6G08_MD_792_raw,pcadata)
points(x6G08_MD_792_proj[3],x6G08_MD_792_proj[2],pch=20)
text(x6G08_MD_792_proj[3],x6G08_MD_792_proj[2],pos=4,label="792",col="red")

x6G08_MD_793_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_793.txt")
x6G08_MD_793_proj<-project.pca(x6G08_MD_793_raw,pcadata)
points(x6G08_MD_793_proj[3],x6G08_MD_793_proj[2],pch=20)
text(x6G08_MD_793_proj[3],x6G08_MD_793_proj[2],pos=4,label="793",col="red")

x6G08_MD_794_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_794.txt")
x6G08_MD_794_proj<-project.pca(x6G08_MD_794_raw,pcadata)
points(x6G08_MD_794_proj[3],x6G08_MD_794_proj[2],pch=20)
text(x6G08_MD_794_proj[3],x6G08_MD_794_proj[2],pos=4,label="794",col="red")

x6G08_MD_795_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_795.txt")
x6G08_MD_795_proj<-project.pca(x6G08_MD_795_raw,pcadata)
points(x6G08_MD_795_proj[3],x6G08_MD_795_proj[2],pch=20)
text(x6G08_MD_795_proj[3],x6G08_MD_795_proj[2],pos=4,label="795",col="red")

x6G08_MD_796_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_796.txt")
x6G08_MD_796_proj<-project.pca(x6G08_MD_796_raw,pcadata)
points(x6G08_MD_796_proj[3],x6G08_MD_796_proj[2],pch=20)
text(x6G08_MD_796_proj[3],x6G08_MD_796_proj[2],pos=4,label="796",col="red")

x6G08_MD_797_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_797.txt")
x6G08_MD_797_proj<-project.pca(x6G08_MD_797_raw,pcadata)
points(x6G08_MD_797_proj[3],x6G08_MD_797_proj[2],pch=20)
text(x6G08_MD_797_proj[3],x6G08_MD_797_proj[2],pos=4,label="797",col="red")

x6G08_MD_798_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_798.txt")
x6G08_MD_798_proj<-project.pca(x6G08_MD_798_raw,pcadata)
points(x6G08_MD_798_proj[3],x6G08_MD_798_proj[2],pch=20)
text(x6G08_MD_798_proj[3],x6G08_MD_798_proj[2],pos=4,label="798",col="red")

x6G08_MD_799_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_799.txt")
x6G08_MD_799_proj<-project.pca(x6G08_MD_799_raw,pcadata)
points(x6G08_MD_799_proj[3],x6G08_MD_799_proj[2],pch=20)
text(x6G08_MD_799_proj[3],x6G08_MD_799_proj[2],pos=4,label="799",col="red")

x6G08_MD_800_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_800.txt")
x6G08_MD_800_proj<-project.pca(x6G08_MD_800_raw,pcadata)
points(x6G08_MD_800_proj[3],x6G08_MD_800_proj[2],pch=20)
text(x6G08_MD_800_proj[3],x6G08_MD_800_proj[2],pos=4,label="800",col="red")

x6G08_MD_801_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_801.txt")
x6G08_MD_801_proj<-project.pca(x6G08_MD_801_raw,pcadata)
points(x6G08_MD_801_proj[3],x6G08_MD_801_proj[2],pch=20)
text(x6G08_MD_801_proj[3],x6G08_MD_801_proj[2],pos=4,label="801",col="red")

x6G08_MD_802_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_802.txt")
x6G08_MD_802_proj<-project.pca(x6G08_MD_802_raw,pcadata)
points(x6G08_MD_802_proj[3],x6G08_MD_802_proj[2],pch=20)
text(x6G08_MD_802_proj[3],x6G08_MD_802_proj[2],pos=4,label="802",col="red")

x6G08_MD_803_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_803.txt")
x6G08_MD_803_proj<-project.pca(x6G08_MD_803_raw,pcadata)
points(x6G08_MD_803_proj[3],x6G08_MD_803_proj[2],pch=20)
text(x6G08_MD_803_proj[3],x6G08_MD_803_proj[2],pos=4,label="803",col="red")

x6G08_MD_804_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_804.txt")
x6G08_MD_804_proj<-project.pca(x6G08_MD_804_raw,pcadata)
points(x6G08_MD_804_proj[3],x6G08_MD_804_proj[2],pch=20)
text(x6G08_MD_804_proj[3],x6G08_MD_804_proj[2],pos=4,label="804",col="red")

x6G08_MD_805_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_805.txt")
x6G08_MD_805_proj<-project.pca(x6G08_MD_805_raw,pcadata)
points(x6G08_MD_805_proj[3],x6G08_MD_805_proj[2],pch=20)
text(x6G08_MD_805_proj[3],x6G08_MD_805_proj[2],pos=4,label="805",col="red")

x6G08_MD_806_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_806.txt")
x6G08_MD_806_proj<-project.pca(x6G08_MD_806_raw,pcadata)
points(x6G08_MD_806_proj[3],x6G08_MD_806_proj[2],pch=20)
text(x6G08_MD_806_proj[3],x6G08_MD_806_proj[2],pos=4,label="806",col="red")

x6G08_MD_807_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_807.txt")
x6G08_MD_807_proj<-project.pca(x6G08_MD_807_raw,pcadata)
points(x6G08_MD_807_proj[3],x6G08_MD_807_proj[2],pch=20)
text(x6G08_MD_807_proj[3],x6G08_MD_807_proj[2],pos=4,label="807",col="red")

x6G08_MD_808_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_808.txt")
x6G08_MD_808_proj<-project.pca(x6G08_MD_808_raw,pcadata)
points(x6G08_MD_808_proj[3],x6G08_MD_808_proj[2],pch=20)
text(x6G08_MD_808_proj[3],x6G08_MD_808_proj[2],pos=4,label="808",col="red")

x6G08_MD_809_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_809.txt")
x6G08_MD_809_proj<-project.pca(x6G08_MD_809_raw,pcadata)
points(x6G08_MD_809_proj[3],x6G08_MD_809_proj[2],pch=20)
text(x6G08_MD_809_proj[3],x6G08_MD_809_proj[2],pos=4,label="809",col="red")

x6G08_MD_810_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_810.txt")
x6G08_MD_810_proj<-project.pca(x6G08_MD_810_raw,pcadata)
points(x6G08_MD_810_proj[3],x6G08_MD_810_proj[2],pch=20)
text(x6G08_MD_810_proj[3],x6G08_MD_810_proj[2],pos=4,label="810",col="red")

x6G08_MD_811_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_811.txt")
x6G08_MD_811_proj<-project.pca(x6G08_MD_811_raw,pcadata)
points(x6G08_MD_811_proj[3],x6G08_MD_811_proj[2],pch=20)
text(x6G08_MD_811_proj[3],x6G08_MD_811_proj[2],pos=4,label="811",col="red")

x6G08_MD_812_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_812.txt")
x6G08_MD_812_proj<-project.pca(x6G08_MD_812_raw,pcadata)
points(x6G08_MD_812_proj[3],x6G08_MD_812_proj[2],pch=20)
text(x6G08_MD_812_proj[3],x6G08_MD_812_proj[2],pos=4,label="812",col="red")

x6G08_MD_813_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_813.txt")
x6G08_MD_813_proj<-project.pca(x6G08_MD_813_raw,pcadata)
points(x6G08_MD_813_proj[3],x6G08_MD_813_proj[2],pch=20)
text(x6G08_MD_813_proj[3],x6G08_MD_813_proj[2],pos=4,label="813",col="red")

x6G08_MD_814_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_814.txt")
x6G08_MD_814_proj<-project.pca(x6G08_MD_814_raw,pcadata)
points(x6G08_MD_814_proj[3],x6G08_MD_814_proj[2],pch=20)
text(x6G08_MD_814_proj[3],x6G08_MD_814_proj[2],pos=4,label="814",col="red")

x6G08_MD_815_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_815.txt")
x6G08_MD_815_proj<-project.pca(x6G08_MD_815_raw,pcadata)
points(x6G08_MD_815_proj[3],x6G08_MD_815_proj[2],pch=20)
text(x6G08_MD_815_proj[3],x6G08_MD_815_proj[2],pos=4,label="815",col="red")

x6G08_MD_816_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_816.txt")
x6G08_MD_816_proj<-project.pca(x6G08_MD_816_raw,pcadata)
points(x6G08_MD_816_proj[3],x6G08_MD_816_proj[2],pch=20)
text(x6G08_MD_816_proj[3],x6G08_MD_816_proj[2],pos=4,label="816",col="blue")

x6G08_MD_817_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_817.txt")
x6G08_MD_817_proj<-project.pca(x6G08_MD_817_raw,pcadata)
points(x6G08_MD_817_proj[3],x6G08_MD_817_proj[2],pch=20)
text(x6G08_MD_817_proj[3],x6G08_MD_817_proj[2],pos=4,label="817",col="red")

x6G08_MD_818_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_818.txt")
x6G08_MD_818_proj<-project.pca(x6G08_MD_818_raw,pcadata)
points(x6G08_MD_818_proj[3],x6G08_MD_818_proj[2],pch=20)
text(x6G08_MD_818_proj[3],x6G08_MD_818_proj[2],pos=4,label="818",col="blue")

x6G08_MD_819_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_819.txt")
x6G08_MD_819_proj<-project.pca(x6G08_MD_819_raw,pcadata)
points(x6G08_MD_819_proj[3],x6G08_MD_819_proj[2],pch=20)
text(x6G08_MD_819_proj[3],x6G08_MD_819_proj[2],pos=4,label="819",col="blue")

x6G08_MD_820_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_820.txt")
x6G08_MD_820_proj<-project.pca(x6G08_MD_820_raw,pcadata)
points(x6G08_MD_820_proj[3],x6G08_MD_820_proj[2],pch=20)
text(x6G08_MD_820_proj[3],x6G08_MD_820_proj[2],pos=4,label="820",col="blue")

x6G08_MD_821_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_821.txt")
x6G08_MD_821_proj<-project.pca(x6G08_MD_821_raw,pcadata)
points(x6G08_MD_821_proj[3],x6G08_MD_821_proj[2],pch=20)
text(x6G08_MD_821_proj[3],x6G08_MD_821_proj[2],pos=4,label="821",col="red")

x6G08_MD_822_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_822.txt")
x6G08_MD_822_proj<-project.pca(x6G08_MD_822_raw,pcadata)
points(x6G08_MD_822_proj[3],x6G08_MD_822_proj[2],pch=20)
text(x6G08_MD_822_proj[3],x6G08_MD_822_proj[2],pos=4,label="822",col="red")

x6G08_MD_823_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_823.txt")
x6G08_MD_823_proj<-project.pca(x6G08_MD_823_raw,pcadata)
points(x6G08_MD_823_proj[3],x6G08_MD_823_proj[2],pch=20)
text(x6G08_MD_823_proj[3],x6G08_MD_823_proj[2],pos=4,label="823",col="red")

x6G08_MD_824_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_824.txt")
x6G08_MD_824_proj<-project.pca(x6G08_MD_824_raw,pcadata)
points(x6G08_MD_824_proj[3],x6G08_MD_824_proj[2],pch=20)
text(x6G08_MD_824_proj[3],x6G08_MD_824_proj[2],pos=4,label="824",col="red")

x6G08_MD_825_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_825.txt")
x6G08_MD_825_proj<-project.pca(x6G08_MD_825_raw,pcadata)
points(x6G08_MD_825_proj[3],x6G08_MD_825_proj[2],pch=20)
text(x6G08_MD_825_proj[3],x6G08_MD_825_proj[2],pos=4,label="825",col="red")

x6G08_MD_826_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_826.txt")
x6G08_MD_826_proj<-project.pca(x6G08_MD_826_raw,pcadata)
points(x6G08_MD_826_proj[3],x6G08_MD_826_proj[2],pch=20)
text(x6G08_MD_826_proj[3],x6G08_MD_826_proj[2],pos=4,label="826",col="red")

x6G08_MD_827_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_827.txt")
x6G08_MD_827_proj<-project.pca(x6G08_MD_827_raw,pcadata)
points(x6G08_MD_827_proj[3],x6G08_MD_827_proj[2],pch=20)
text(x6G08_MD_827_proj[3],x6G08_MD_827_proj[2],pos=4,label="827",col="red")

x6G08_MD_828_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_828.txt")
x6G08_MD_828_proj<-project.pca(x6G08_MD_828_raw,pcadata)
points(x6G08_MD_828_proj[3],x6G08_MD_828_proj[2],pch=20)
text(x6G08_MD_828_proj[3],x6G08_MD_828_proj[2],pos=4,label="828",col="red")

x6G08_MD_829_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_829.txt")
x6G08_MD_829_proj<-project.pca(x6G08_MD_829_raw,pcadata)
points(x6G08_MD_829_proj[3],x6G08_MD_829_proj[2],pch=20)
text(x6G08_MD_829_proj[3],x6G08_MD_829_proj[2],pos=4,label="829",col="red")

x6G08_MD_830_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_830.txt")
x6G08_MD_830_proj<-project.pca(x6G08_MD_830_raw,pcadata)
points(x6G08_MD_830_proj[3],x6G08_MD_830_proj[2],pch=20)
text(x6G08_MD_830_proj[3],x6G08_MD_830_proj[2],pos=4,label="830",col="red")

x6G08_MD_831_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_831.txt")
x6G08_MD_831_proj<-project.pca(x6G08_MD_831_raw,pcadata)
points(x6G08_MD_831_proj[3],x6G08_MD_831_proj[2],pch=20)
text(x6G08_MD_831_proj[3],x6G08_MD_831_proj[2],pos=4,label="831",col="red")

x6G08_MD_832_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_832.txt")
x6G08_MD_832_proj<-project.pca(x6G08_MD_832_raw,pcadata)
points(x6G08_MD_832_proj[3],x6G08_MD_832_proj[2],pch=20)
text(x6G08_MD_832_proj[3],x6G08_MD_832_proj[2],pos=4,label="832",col="red")

x6G08_MD_833_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_833.txt")
x6G08_MD_833_proj<-project.pca(x6G08_MD_833_raw,pcadata)
points(x6G08_MD_833_proj[3],x6G08_MD_833_proj[2],pch=20)
text(x6G08_MD_833_proj[3],x6G08_MD_833_proj[2],pos=4,label="833",col="red")

x6G08_MD_834_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_834.txt")
x6G08_MD_834_proj<-project.pca(x6G08_MD_834_raw,pcadata)
points(x6G08_MD_834_proj[3],x6G08_MD_834_proj[2],pch=20)
text(x6G08_MD_834_proj[3],x6G08_MD_834_proj[2],pos=4,label="834",col="red")

x6G08_MD_835_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_835.txt")
x6G08_MD_835_proj<-project.pca(x6G08_MD_835_raw,pcadata)
points(x6G08_MD_835_proj[3],x6G08_MD_835_proj[2],pch=20)
text(x6G08_MD_835_proj[3],x6G08_MD_835_proj[2],pos=4,label="835",col="red")

x6G08_MD_836_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_836.txt")
x6G08_MD_836_proj<-project.pca(x6G08_MD_836_raw,pcadata)
points(x6G08_MD_836_proj[3],x6G08_MD_836_proj[2],pch=20)
text(x6G08_MD_836_proj[3],x6G08_MD_836_proj[2],pos=4,label="836",col="red")

x6G08_MD_837_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_837.txt")
x6G08_MD_837_proj<-project.pca(x6G08_MD_837_raw,pcadata)
points(x6G08_MD_837_proj[3],x6G08_MD_837_proj[2],pch=20)
text(x6G08_MD_837_proj[3],x6G08_MD_837_proj[2],pos=4,label="837",col="red")

x6G08_MD_838_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_838.txt")
x6G08_MD_838_proj<-project.pca(x6G08_MD_838_raw,pcadata)
points(x6G08_MD_838_proj[3],x6G08_MD_838_proj[2],pch=20)
text(x6G08_MD_838_proj[3],x6G08_MD_838_proj[2],pos=4,label="838",col="red")

x6G08_MD_839_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_839.txt")
x6G08_MD_839_proj<-project.pca(x6G08_MD_839_raw,pcadata)
points(x6G08_MD_839_proj[3],x6G08_MD_839_proj[2],pch=20)
text(x6G08_MD_839_proj[3],x6G08_MD_839_proj[2],pos=4,label="839",col="red")

x6G08_MD_840_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_840.txt")
x6G08_MD_840_proj<-project.pca(x6G08_MD_840_raw,pcadata)
points(x6G08_MD_840_proj[3],x6G08_MD_840_proj[2],pch=20)
text(x6G08_MD_840_proj[3],x6G08_MD_840_proj[2],pos=4,label="840",col="red")

x6G08_MD_841_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_841.txt")
x6G08_MD_841_proj<-project.pca(x6G08_MD_841_raw,pcadata)
points(x6G08_MD_841_proj[3],x6G08_MD_841_proj[2],pch=20)
text(x6G08_MD_841_proj[3],x6G08_MD_841_proj[2],pos=4,label="841",col="red")

x6G08_MD_842_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_842.txt")
x6G08_MD_842_proj<-project.pca(x6G08_MD_842_raw,pcadata)
points(x6G08_MD_842_proj[3],x6G08_MD_842_proj[2],pch=20)
text(x6G08_MD_842_proj[3],x6G08_MD_842_proj[2],pos=4,label="842",col="green")

x6G08_MD_843_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_843.txt")
x6G08_MD_843_proj<-project.pca(x6G08_MD_843_raw,pcadata)
points(x6G08_MD_843_proj[3],x6G08_MD_843_proj[2],pch=20)
text(x6G08_MD_843_proj[3],x6G08_MD_843_proj[2],pos=4,label="843",col="red")

x6G08_MD_844_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_844.txt")
x6G08_MD_844_proj<-project.pca(x6G08_MD_844_raw,pcadata)
points(x6G08_MD_844_proj[3],x6G08_MD_844_proj[2],pch=20)
text(x6G08_MD_844_proj[3],x6G08_MD_844_proj[2],pos=4,label="844",col="red")

x6G08_MD_845_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_845.txt")
x6G08_MD_845_proj<-project.pca(x6G08_MD_845_raw,pcadata)
points(x6G08_MD_845_proj[3],x6G08_MD_845_proj[2],pch=20)
text(x6G08_MD_845_proj[3],x6G08_MD_845_proj[2],pos=4,label="845",col="red")

x6G08_MD_846_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_846.txt")
x6G08_MD_846_proj<-project.pca(x6G08_MD_846_raw,pcadata)
points(x6G08_MD_846_proj[3],x6G08_MD_846_proj[2],pch=20)
text(x6G08_MD_846_proj[3],x6G08_MD_846_proj[2],pos=4,label="846",col="red")

x6G08_MD_847_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_847.txt")
x6G08_MD_847_proj<-project.pca(x6G08_MD_847_raw,pcadata)
points(x6G08_MD_847_proj[3],x6G08_MD_847_proj[2],pch=20)
text(x6G08_MD_847_proj[3],x6G08_MD_847_proj[2],pos=4,label="847",col="red")

x6G08_MD_848_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_848.txt")
x6G08_MD_848_proj<-project.pca(x6G08_MD_848_raw,pcadata)
points(x6G08_MD_848_proj[3],x6G08_MD_848_proj[2],pch=20)
text(x6G08_MD_848_proj[3],x6G08_MD_848_proj[2],pos=4,label="848",col="red")

x6G08_MD_849_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_849.txt")
x6G08_MD_849_proj<-project.pca(x6G08_MD_849_raw,pcadata)
points(x6G08_MD_849_proj[3],x6G08_MD_849_proj[2],pch=20)
text(x6G08_MD_849_proj[3],x6G08_MD_849_proj[2],pos=4,label="849",col="red")

x6G08_MD_850_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_850.txt")
x6G08_MD_850_proj<-project.pca(x6G08_MD_850_raw,pcadata)
points(x6G08_MD_850_proj[3],x6G08_MD_850_proj[2],pch=20)
text(x6G08_MD_850_proj[3],x6G08_MD_850_proj[2],pos=4,label="850",col="red")

x6G08_MD_851_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_851.txt")
x6G08_MD_851_proj<-project.pca(x6G08_MD_851_raw,pcadata)
points(x6G08_MD_851_proj[3],x6G08_MD_851_proj[2],pch=20)
text(x6G08_MD_851_proj[3],x6G08_MD_851_proj[2],pos=4,label="851",col="red")

x6G08_MD_852_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_852.txt")
x6G08_MD_852_proj<-project.pca(x6G08_MD_852_raw,pcadata)
points(x6G08_MD_852_proj[3],x6G08_MD_852_proj[2],pch=20)
text(x6G08_MD_852_proj[3],x6G08_MD_852_proj[2],pos=4,label="852",col="red")

x6G08_MD_853_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_853.txt")
x6G08_MD_853_proj<-project.pca(x6G08_MD_853_raw,pcadata)
points(x6G08_MD_853_proj[3],x6G08_MD_853_proj[2],pch=20)
text(x6G08_MD_853_proj[3],x6G08_MD_853_proj[2],pos=4,label="853",col="red")

x6G08_MD_854_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_854.txt")
x6G08_MD_854_proj<-project.pca(x6G08_MD_854_raw,pcadata)
points(x6G08_MD_854_proj[3],x6G08_MD_854_proj[2],pch=20)
text(x6G08_MD_854_proj[3],x6G08_MD_854_proj[2],pos=4,label="854",col="red")

x6G08_MD_855_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_855.txt")
x6G08_MD_855_proj<-project.pca(x6G08_MD_855_raw,pcadata)
points(x6G08_MD_855_proj[3],x6G08_MD_855_proj[2],pch=20)
text(x6G08_MD_855_proj[3],x6G08_MD_855_proj[2],pos=4,label="855",col="green")

x6G08_MD_856_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_856.txt")
x6G08_MD_856_proj<-project.pca(x6G08_MD_856_raw,pcadata)
points(x6G08_MD_856_proj[3],x6G08_MD_856_proj[2],pch=20)
text(x6G08_MD_856_proj[3],x6G08_MD_856_proj[2],pos=4,label="856",col="red")

x6G08_MD_857_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_857.txt")
x6G08_MD_857_proj<-project.pca(x6G08_MD_857_raw,pcadata)
points(x6G08_MD_857_proj[3],x6G08_MD_857_proj[2],pch=20)
text(x6G08_MD_857_proj[3],x6G08_MD_857_proj[2],pos=4,label="857",col="red")

x6G08_MD_858_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_858.txt")
x6G08_MD_858_proj<-project.pca(x6G08_MD_858_raw,pcadata)
points(x6G08_MD_858_proj[3],x6G08_MD_858_proj[2],pch=20)
text(x6G08_MD_858_proj[3],x6G08_MD_858_proj[2],pos=4,label="858",col="red")

x6G08_MD_859_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_859.txt")
x6G08_MD_859_proj<-project.pca(x6G08_MD_859_raw,pcadata)
points(x6G08_MD_859_proj[3],x6G08_MD_859_proj[2],pch=20)
text(x6G08_MD_859_proj[3],x6G08_MD_859_proj[2],pos=4,label="859",col="red")

x6G08_MD_860_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_860.txt")
x6G08_MD_860_proj<-project.pca(x6G08_MD_860_raw,pcadata)
points(x6G08_MD_860_proj[3],x6G08_MD_860_proj[2],pch=20)
text(x6G08_MD_860_proj[3],x6G08_MD_860_proj[2],pos=4,label="860",col="red")

x6G08_MD_861_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_861.txt")
x6G08_MD_861_proj<-project.pca(x6G08_MD_861_raw,pcadata)
points(x6G08_MD_861_proj[3],x6G08_MD_861_proj[2],pch=20)
text(x6G08_MD_861_proj[3],x6G08_MD_861_proj[2],pos=4,label="861",col="red")

x6G08_MD_862_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_862.txt")
x6G08_MD_862_proj<-project.pca(x6G08_MD_862_raw,pcadata)
points(x6G08_MD_862_proj[3],x6G08_MD_862_proj[2],pch=20)
text(x6G08_MD_862_proj[3],x6G08_MD_862_proj[2],pos=4,label="862",col="red")

x6G08_MD_863_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_863.txt")
x6G08_MD_863_proj<-project.pca(x6G08_MD_863_raw,pcadata)
points(x6G08_MD_863_proj[3],x6G08_MD_863_proj[2],pch=20)
text(x6G08_MD_863_proj[3],x6G08_MD_863_proj[2],pos=4,label="863",col="red")

x6G08_MD_864_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_864.txt")
x6G08_MD_864_proj<-project.pca(x6G08_MD_864_raw,pcadata)
points(x6G08_MD_864_proj[3],x6G08_MD_864_proj[2],pch=20)
text(x6G08_MD_864_proj[3],x6G08_MD_864_proj[2],pos=4,label="864",col="red")

x6G08_MD_865_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_865.txt")
x6G08_MD_865_proj<-project.pca(x6G08_MD_865_raw,pcadata)
points(x6G08_MD_865_proj[3],x6G08_MD_865_proj[2],pch=20)
text(x6G08_MD_865_proj[3],x6G08_MD_865_proj[2],pos=4,label="865",col="red")

x6G08_MD_866_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_866.txt")
x6G08_MD_866_proj<-project.pca(x6G08_MD_866_raw,pcadata)
points(x6G08_MD_866_proj[3],x6G08_MD_866_proj[2],pch=20)
text(x6G08_MD_866_proj[3],x6G08_MD_866_proj[2],pos=4,label="866",col="red")

x6G08_MD_867_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_867.txt")
x6G08_MD_867_proj<-project.pca(x6G08_MD_867_raw,pcadata)
points(x6G08_MD_867_proj[3],x6G08_MD_867_proj[2],pch=20)
text(x6G08_MD_867_proj[3],x6G08_MD_867_proj[2],pos=4,label="867",col="red")

x6G08_MD_868_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_868.txt")
x6G08_MD_868_proj<-project.pca(x6G08_MD_868_raw,pcadata)
points(x6G08_MD_868_proj[3],x6G08_MD_868_proj[2],pch=20)
text(x6G08_MD_868_proj[3],x6G08_MD_868_proj[2],pos=4,label="868",col="red")

x6G08_MD_869_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_869.txt")
x6G08_MD_869_proj<-project.pca(x6G08_MD_869_raw,pcadata)
points(x6G08_MD_869_proj[3],x6G08_MD_869_proj[2],pch=20)
text(x6G08_MD_869_proj[3],x6G08_MD_869_proj[2],pos=4,label="869",col="red")

x6G08_MD_870_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_870.txt")
x6G08_MD_870_proj<-project.pca(x6G08_MD_870_raw,pcadata)
points(x6G08_MD_870_proj[3],x6G08_MD_870_proj[2],pch=20)
text(x6G08_MD_870_proj[3],x6G08_MD_870_proj[2],pos=4,label="870",col="red")

x6G08_MD_871_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_871.txt")
x6G08_MD_871_proj<-project.pca(x6G08_MD_871_raw,pcadata)
points(x6G08_MD_871_proj[3],x6G08_MD_871_proj[2],pch=20)
text(x6G08_MD_871_proj[3],x6G08_MD_871_proj[2],pos=4,label="871",col="red")

x6G08_MD_872_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_872.txt")
x6G08_MD_872_proj<-project.pca(x6G08_MD_872_raw,pcadata)
points(x6G08_MD_872_proj[3],x6G08_MD_872_proj[2],pch=20)
text(x6G08_MD_872_proj[3],x6G08_MD_872_proj[2],pos=4,label="872",col="red")

x6G08_MD_873_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_873.txt")
x6G08_MD_873_proj<-project.pca(x6G08_MD_873_raw,pcadata)
points(x6G08_MD_873_proj[3],x6G08_MD_873_proj[2],pch=20)
text(x6G08_MD_873_proj[3],x6G08_MD_873_proj[2],pos=4,label="873",col="red")

x6G08_MD_874_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_874.txt")
x6G08_MD_874_proj<-project.pca(x6G08_MD_874_raw,pcadata)
points(x6G08_MD_874_proj[3],x6G08_MD_874_proj[2],pch=20)
text(x6G08_MD_874_proj[3],x6G08_MD_874_proj[2],pos=4,label="874",col="red")

x6G08_MD_875_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_875.txt")
x6G08_MD_875_proj<-project.pca(x6G08_MD_875_raw,pcadata)
points(x6G08_MD_875_proj[3],x6G08_MD_875_proj[2],pch=20)
text(x6G08_MD_875_proj[3],x6G08_MD_875_proj[2],pos=4,label="875",col="red")

x6G08_MD_876_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_876.txt")
x6G08_MD_876_proj<-project.pca(x6G08_MD_876_raw,pcadata)
points(x6G08_MD_876_proj[3],x6G08_MD_876_proj[2],pch=20)
text(x6G08_MD_876_proj[3],x6G08_MD_876_proj[2],pos=4,label="876",col="red")

x6G08_MD_877_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_877.txt")
x6G08_MD_877_proj<-project.pca(x6G08_MD_877_raw,pcadata)
points(x6G08_MD_877_proj[3],x6G08_MD_877_proj[2],pch=20)
text(x6G08_MD_877_proj[3],x6G08_MD_877_proj[2],pos=4,label="877",col="red")

x6G08_MD_878_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_878.txt")
x6G08_MD_878_proj<-project.pca(x6G08_MD_878_raw,pcadata)
points(x6G08_MD_878_proj[3],x6G08_MD_878_proj[2],pch=20)
text(x6G08_MD_878_proj[3],x6G08_MD_878_proj[2],pos=4,label="878",col="red")

x6G08_MD_879_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_879.txt")
x6G08_MD_879_proj<-project.pca(x6G08_MD_879_raw,pcadata)
points(x6G08_MD_879_proj[3],x6G08_MD_879_proj[2],pch=20)
text(x6G08_MD_879_proj[3],x6G08_MD_879_proj[2],pos=4,label="879",col="red")

x6G08_MD_880_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_880.txt")
x6G08_MD_880_proj<-project.pca(x6G08_MD_880_raw,pcadata)
points(x6G08_MD_880_proj[3],x6G08_MD_880_proj[2],pch=20)
text(x6G08_MD_880_proj[3],x6G08_MD_880_proj[2],pos=4,label="880",col="red")

x6G08_MD_881_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_881.txt")
x6G08_MD_881_proj<-project.pca(x6G08_MD_881_raw,pcadata)
points(x6G08_MD_881_proj[3],x6G08_MD_881_proj[2],pch=20)
text(x6G08_MD_881_proj[3],x6G08_MD_881_proj[2],pos=4,label="881",col="red")

x6G08_MD_882_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_882.txt")
x6G08_MD_882_proj<-project.pca(x6G08_MD_882_raw,pcadata)
points(x6G08_MD_882_proj[3],x6G08_MD_882_proj[2],pch=20)
text(x6G08_MD_882_proj[3],x6G08_MD_882_proj[2],pos=4,label="882",col="red")

x6G08_MD_883_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_883.txt")
x6G08_MD_883_proj<-project.pca(x6G08_MD_883_raw,pcadata)
points(x6G08_MD_883_proj[3],x6G08_MD_883_proj[2],pch=20)
text(x6G08_MD_883_proj[3],x6G08_MD_883_proj[2],pos=4,label="883",col="red")

x6G08_MD_884_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_884.txt")
x6G08_MD_884_proj<-project.pca(x6G08_MD_884_raw,pcadata)
points(x6G08_MD_884_proj[3],x6G08_MD_884_proj[2],pch=20)
text(x6G08_MD_884_proj[3],x6G08_MD_884_proj[2],pos=4,label="884",col="red")

x6G08_MD_885_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_885.txt")
x6G08_MD_885_proj<-project.pca(x6G08_MD_885_raw,pcadata)
points(x6G08_MD_885_proj[3],x6G08_MD_885_proj[2],pch=20)
text(x6G08_MD_885_proj[3],x6G08_MD_885_proj[2],pos=4,label="885",col="red")

x6G08_MD_886_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_886.txt")
x6G08_MD_886_proj<-project.pca(x6G08_MD_886_raw,pcadata)
points(x6G08_MD_886_proj[3],x6G08_MD_886_proj[2],pch=20)
text(x6G08_MD_886_proj[3],x6G08_MD_886_proj[2],pos=4,label="886",col="red")

x6G08_MD_887_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_887.txt")
x6G08_MD_887_proj<-project.pca(x6G08_MD_887_raw,pcadata)
points(x6G08_MD_887_proj[3],x6G08_MD_887_proj[2],pch=20)
text(x6G08_MD_887_proj[3],x6G08_MD_887_proj[2],pos=4,label="887",col="red")

x6G08_MD_888_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_888.txt")
x6G08_MD_888_proj<-project.pca(x6G08_MD_888_raw,pcadata)
points(x6G08_MD_888_proj[3],x6G08_MD_888_proj[2],pch=20)
text(x6G08_MD_888_proj[3],x6G08_MD_888_proj[2],pos=4,label="888",col="red")

x6G08_MD_889_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_889.txt")
x6G08_MD_889_proj<-project.pca(x6G08_MD_889_raw,pcadata)
points(x6G08_MD_889_proj[3],x6G08_MD_889_proj[2],pch=20)
text(x6G08_MD_889_proj[3],x6G08_MD_889_proj[2],pos=4,label="889",col="red")

x6G08_MD_890_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_890.txt")
x6G08_MD_890_proj<-project.pca(x6G08_MD_890_raw,pcadata)
points(x6G08_MD_890_proj[3],x6G08_MD_890_proj[2],pch=20)
text(x6G08_MD_890_proj[3],x6G08_MD_890_proj[2],pos=4,label="890",col="red")

x6G08_MD_891_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_891.txt")
x6G08_MD_891_proj<-project.pca(x6G08_MD_891_raw,pcadata)
points(x6G08_MD_891_proj[3],x6G08_MD_891_proj[2],pch=20)
text(x6G08_MD_891_proj[3],x6G08_MD_891_proj[2],pos=4,label="891",col="red")

x6G08_MD_892_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_892.txt")
x6G08_MD_892_proj<-project.pca(x6G08_MD_892_raw,pcadata)
points(x6G08_MD_892_proj[3],x6G08_MD_892_proj[2],pch=20)
text(x6G08_MD_892_proj[3],x6G08_MD_892_proj[2],pos=4,label="892",col="red")

x6G08_MD_893_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_893.txt")
x6G08_MD_893_proj<-project.pca(x6G08_MD_893_raw,pcadata)
points(x6G08_MD_893_proj[3],x6G08_MD_893_proj[2],pch=20)
text(x6G08_MD_893_proj[3],x6G08_MD_893_proj[2],pos=4,label="893",col="red")

x6G08_MD_894_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_894.txt")
x6G08_MD_894_proj<-project.pca(x6G08_MD_894_raw,pcadata)
points(x6G08_MD_894_proj[3],x6G08_MD_894_proj[2],pch=20)
text(x6G08_MD_894_proj[3],x6G08_MD_894_proj[2],pos=4,label="894",col="red")

x6G08_MD_895_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_895.txt")
x6G08_MD_895_proj<-project.pca(x6G08_MD_895_raw,pcadata)
points(x6G08_MD_895_proj[3],x6G08_MD_895_proj[2],pch=20)
text(x6G08_MD_895_proj[3],x6G08_MD_895_proj[2],pos=4,label="895",col="green")

x6G08_MD_896_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_896.txt")
x6G08_MD_896_proj<-project.pca(x6G08_MD_896_raw,pcadata)
points(x6G08_MD_896_proj[3],x6G08_MD_896_proj[2],pch=20)
text(x6G08_MD_896_proj[3],x6G08_MD_896_proj[2],pos=4,label="896",col="red")

x6G08_MD_897_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_897.txt")
x6G08_MD_897_proj<-project.pca(x6G08_MD_897_raw,pcadata)
points(x6G08_MD_897_proj[3],x6G08_MD_897_proj[2],pch=20)
text(x6G08_MD_897_proj[3],x6G08_MD_897_proj[2],pos=4,label="897",col="green")

x6G08_MD_898_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_898.txt")
x6G08_MD_898_proj<-project.pca(x6G08_MD_898_raw,pcadata)
points(x6G08_MD_898_proj[3],x6G08_MD_898_proj[2],pch=20)
text(x6G08_MD_898_proj[3],x6G08_MD_898_proj[2],pos=4,label="898",col="green")

x6G08_MD_899_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_899.txt")
x6G08_MD_899_proj<-project.pca(x6G08_MD_899_raw,pcadata)
points(x6G08_MD_899_proj[3],x6G08_MD_899_proj[2],pch=20)
text(x6G08_MD_899_proj[3],x6G08_MD_899_proj[2],pos=4,label="899",col="red")

x6G08_MD_900_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_900.txt")
x6G08_MD_900_proj<-project.pca(x6G08_MD_900_raw,pcadata)
points(x6G08_MD_900_proj[3],x6G08_MD_900_proj[2],pch=20)
text(x6G08_MD_900_proj[3],x6G08_MD_900_proj[2],pos=4,label="900",col="red")

x6G08_MD_901_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_901.txt")
x6G08_MD_901_proj<-project.pca(x6G08_MD_901_raw,pcadata)
points(x6G08_MD_901_proj[3],x6G08_MD_901_proj[2],pch=20)
text(x6G08_MD_901_proj[3],x6G08_MD_901_proj[2],pos=4,label="901",col="red")

x6G08_MD_902_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_902.txt")
x6G08_MD_902_proj<-project.pca(x6G08_MD_902_raw,pcadata)
points(x6G08_MD_902_proj[3],x6G08_MD_902_proj[2],pch=20)
text(x6G08_MD_902_proj[3],x6G08_MD_902_proj[2],pos=4,label="902",col="red")

x6G08_MD_903_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_903.txt")
x6G08_MD_903_proj<-project.pca(x6G08_MD_903_raw,pcadata)
points(x6G08_MD_903_proj[3],x6G08_MD_903_proj[2],pch=20)
text(x6G08_MD_903_proj[3],x6G08_MD_903_proj[2],pos=4,label="903",col="red")

x6G08_MD_904_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_904.txt")
x6G08_MD_904_proj<-project.pca(x6G08_MD_904_raw,pcadata)
points(x6G08_MD_904_proj[3],x6G08_MD_904_proj[2],pch=20)
text(x6G08_MD_904_proj[3],x6G08_MD_904_proj[2],pos=4,label="904",col="red")

x6G08_MD_905_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_905.txt")
x6G08_MD_905_proj<-project.pca(x6G08_MD_905_raw,pcadata)
points(x6G08_MD_905_proj[3],x6G08_MD_905_proj[2],pch=20)
text(x6G08_MD_905_proj[3],x6G08_MD_905_proj[2],pos=4,label="905",col="red")

x6G08_MD_906_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_906.txt")
x6G08_MD_906_proj<-project.pca(x6G08_MD_906_raw,pcadata)
points(x6G08_MD_906_proj[3],x6G08_MD_906_proj[2],pch=20)
text(x6G08_MD_906_proj[3],x6G08_MD_906_proj[2],pos=4,label="906",col="red")

x6G08_MD_907_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_907.txt")
x6G08_MD_907_proj<-project.pca(x6G08_MD_907_raw,pcadata)
points(x6G08_MD_907_proj[3],x6G08_MD_907_proj[2],pch=20)
text(x6G08_MD_907_proj[3],x6G08_MD_907_proj[2],pos=4,label="907",col="red")

x6G08_MD_908_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_908.txt")
x6G08_MD_908_proj<-project.pca(x6G08_MD_908_raw,pcadata)
points(x6G08_MD_908_proj[3],x6G08_MD_908_proj[2],pch=20)
text(x6G08_MD_908_proj[3],x6G08_MD_908_proj[2],pos=4,label="908",col="red")

x6G08_MD_909_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_909.txt")
x6G08_MD_909_proj<-project.pca(x6G08_MD_909_raw,pcadata)
points(x6G08_MD_909_proj[3],x6G08_MD_909_proj[2],pch=20)
text(x6G08_MD_909_proj[3],x6G08_MD_909_proj[2],pos=4,label="909",col="red")

x6G08_MD_910_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_910.txt")
x6G08_MD_910_proj<-project.pca(x6G08_MD_910_raw,pcadata)
points(x6G08_MD_910_proj[3],x6G08_MD_910_proj[2],pch=20)
text(x6G08_MD_910_proj[3],x6G08_MD_910_proj[2],pos=4,label="910",col="red")

x6G08_MD_911_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_911.txt")
x6G08_MD_911_proj<-project.pca(x6G08_MD_911_raw,pcadata)
points(x6G08_MD_911_proj[3],x6G08_MD_911_proj[2],pch=20)
text(x6G08_MD_911_proj[3],x6G08_MD_911_proj[2],pos=4,label="911",col="red")

x6G08_MD_912_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_912.txt")
x6G08_MD_912_proj<-project.pca(x6G08_MD_912_raw,pcadata)
points(x6G08_MD_912_proj[3],x6G08_MD_912_proj[2],pch=20)
text(x6G08_MD_912_proj[3],x6G08_MD_912_proj[2],pos=4,label="912",col="red")

x6G08_MD_913_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_913.txt")
x6G08_MD_913_proj<-project.pca(x6G08_MD_913_raw,pcadata)
points(x6G08_MD_913_proj[3],x6G08_MD_913_proj[2],pch=20)
text(x6G08_MD_913_proj[3],x6G08_MD_913_proj[2],pos=4,label="913",col="green")

x6G08_MD_914_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_914.txt")
x6G08_MD_914_proj<-project.pca(x6G08_MD_914_raw,pcadata)
points(x6G08_MD_914_proj[3],x6G08_MD_914_proj[2],pch=20)
text(x6G08_MD_914_proj[3],x6G08_MD_914_proj[2],pos=4,label="914",col="red")

x6G08_MD_915_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_915.txt")
x6G08_MD_915_proj<-project.pca(x6G08_MD_915_raw,pcadata)
points(x6G08_MD_915_proj[3],x6G08_MD_915_proj[2],pch=20)
text(x6G08_MD_915_proj[3],x6G08_MD_915_proj[2],pos=4,label="915",col="red")

x6G08_MD_916_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_916.txt")
x6G08_MD_916_proj<-project.pca(x6G08_MD_916_raw,pcadata)
points(x6G08_MD_916_proj[3],x6G08_MD_916_proj[2],pch=20)
text(x6G08_MD_916_proj[3],x6G08_MD_916_proj[2],pos=4,label="916",col="red")

x6G08_MD_917_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_917.txt")
x6G08_MD_917_proj<-project.pca(x6G08_MD_917_raw,pcadata)
points(x6G08_MD_917_proj[3],x6G08_MD_917_proj[2],pch=20)
text(x6G08_MD_917_proj[3],x6G08_MD_917_proj[2],pos=4,label="917",col="red")

x6G08_MD_918_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_918.txt")
x6G08_MD_918_proj<-project.pca(x6G08_MD_918_raw,pcadata)
points(x6G08_MD_918_proj[3],x6G08_MD_918_proj[2],pch=20)
text(x6G08_MD_918_proj[3],x6G08_MD_918_proj[2],pos=4,label="918",col="red")

x6G08_MD_919_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_919.txt")
x6G08_MD_919_proj<-project.pca(x6G08_MD_919_raw,pcadata)
points(x6G08_MD_919_proj[3],x6G08_MD_919_proj[2],pch=20)
text(x6G08_MD_919_proj[3],x6G08_MD_919_proj[2],pos=4,label="919",col="red")

x6G08_MD_920_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_920.txt")
x6G08_MD_920_proj<-project.pca(x6G08_MD_920_raw,pcadata)
points(x6G08_MD_920_proj[3],x6G08_MD_920_proj[2],pch=20)
text(x6G08_MD_920_proj[3],x6G08_MD_920_proj[2],pos=4,label="920",col="red")

x6G08_MD_921_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_921.txt")
x6G08_MD_921_proj<-project.pca(x6G08_MD_921_raw,pcadata)
points(x6G08_MD_921_proj[3],x6G08_MD_921_proj[2],pch=20)
text(x6G08_MD_921_proj[3],x6G08_MD_921_proj[2],pos=4,label="921",col="red")

x6G08_MD_922_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_922.txt")
x6G08_MD_922_proj<-project.pca(x6G08_MD_922_raw,pcadata)
points(x6G08_MD_922_proj[3],x6G08_MD_922_proj[2],pch=20)
text(x6G08_MD_922_proj[3],x6G08_MD_922_proj[2],pos=4,label="922",col="red")

x6G08_MD_923_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_923.txt")
x6G08_MD_923_proj<-project.pca(x6G08_MD_923_raw,pcadata)
points(x6G08_MD_923_proj[3],x6G08_MD_923_proj[2],pch=20)
text(x6G08_MD_923_proj[3],x6G08_MD_923_proj[2],pos=4,label="923",col="red")

x6G08_MD_924_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_924.txt")
x6G08_MD_924_proj<-project.pca(x6G08_MD_924_raw,pcadata)
points(x6G08_MD_924_proj[3],x6G08_MD_924_proj[2],pch=20)
text(x6G08_MD_924_proj[3],x6G08_MD_924_proj[2],pos=4,label="924",col="red")

x6G08_MD_925_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_925.txt")
x6G08_MD_925_proj<-project.pca(x6G08_MD_925_raw,pcadata)
points(x6G08_MD_925_proj[3],x6G08_MD_925_proj[2],pch=20)
text(x6G08_MD_925_proj[3],x6G08_MD_925_proj[2],pos=4,label="925",col="red")

x6G08_MD_926_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_926.txt")
x6G08_MD_926_proj<-project.pca(x6G08_MD_926_raw,pcadata)
points(x6G08_MD_926_proj[3],x6G08_MD_926_proj[2],pch=20)
text(x6G08_MD_926_proj[3],x6G08_MD_926_proj[2],pos=4,label="926",col="red")

x6G08_MD_927_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_927.txt")
x6G08_MD_927_proj<-project.pca(x6G08_MD_927_raw,pcadata)
points(x6G08_MD_927_proj[3],x6G08_MD_927_proj[2],pch=20)
text(x6G08_MD_927_proj[3],x6G08_MD_927_proj[2],pos=4,label="927",col="red")

x6G08_MD_928_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_928.txt")
x6G08_MD_928_proj<-project.pca(x6G08_MD_928_raw,pcadata)
points(x6G08_MD_928_proj[3],x6G08_MD_928_proj[2],pch=20)
text(x6G08_MD_928_proj[3],x6G08_MD_928_proj[2],pos=4,label="928",col="red")

x6G08_MD_929_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_929.txt")
x6G08_MD_929_proj<-project.pca(x6G08_MD_929_raw,pcadata)
points(x6G08_MD_929_proj[3],x6G08_MD_929_proj[2],pch=20)
text(x6G08_MD_929_proj[3],x6G08_MD_929_proj[2],pos=4,label="929",col="red")

x6G08_MD_930_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_930.txt")
x6G08_MD_930_proj<-project.pca(x6G08_MD_930_raw,pcadata)
points(x6G08_MD_930_proj[3],x6G08_MD_930_proj[2],pch=20)
text(x6G08_MD_930_proj[3],x6G08_MD_930_proj[2],pos=4,label="930",col="red")

x6G08_MD_931_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_931.txt")
x6G08_MD_931_proj<-project.pca(x6G08_MD_931_raw,pcadata)
points(x6G08_MD_931_proj[3],x6G08_MD_931_proj[2],pch=20)
text(x6G08_MD_931_proj[3],x6G08_MD_931_proj[2],pos=4,label="931",col="red")

x6G08_MD_932_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_932.txt")
x6G08_MD_932_proj<-project.pca(x6G08_MD_932_raw,pcadata)
points(x6G08_MD_932_proj[3],x6G08_MD_932_proj[2],pch=20)
text(x6G08_MD_932_proj[3],x6G08_MD_932_proj[2],pos=4,label="932",col="red")

x6G08_MD_933_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_933.txt")
x6G08_MD_933_proj<-project.pca(x6G08_MD_933_raw,pcadata)
points(x6G08_MD_933_proj[3],x6G08_MD_933_proj[2],pch=20)
text(x6G08_MD_933_proj[3],x6G08_MD_933_proj[2],pos=4,label="933",col="red")

x6G08_MD_934_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_934.txt")
x6G08_MD_934_proj<-project.pca(x6G08_MD_934_raw,pcadata)
points(x6G08_MD_934_proj[3],x6G08_MD_934_proj[2],pch=20)
text(x6G08_MD_934_proj[3],x6G08_MD_934_proj[2],pos=4,label="934",col="red")

x6G08_MD_935_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_935.txt")
x6G08_MD_935_proj<-project.pca(x6G08_MD_935_raw,pcadata)
points(x6G08_MD_935_proj[3],x6G08_MD_935_proj[2],pch=20)
text(x6G08_MD_935_proj[3],x6G08_MD_935_proj[2],pos=4,label="935",col="red")

x6G08_MD_936_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_936.txt")
x6G08_MD_936_proj<-project.pca(x6G08_MD_936_raw,pcadata)
points(x6G08_MD_936_proj[3],x6G08_MD_936_proj[2],pch=20)
text(x6G08_MD_936_proj[3],x6G08_MD_936_proj[2],pos=4,label="936",col="red")

x6G08_MD_937_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_937.txt")
x6G08_MD_937_proj<-project.pca(x6G08_MD_937_raw,pcadata)
points(x6G08_MD_937_proj[3],x6G08_MD_937_proj[2],pch=20)
text(x6G08_MD_937_proj[3],x6G08_MD_937_proj[2],pos=4,label="937",col="green")

x6G08_MD_938_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_938.txt")
x6G08_MD_938_proj<-project.pca(x6G08_MD_938_raw,pcadata)
points(x6G08_MD_938_proj[3],x6G08_MD_938_proj[2],pch=20)
text(x6G08_MD_938_proj[3],x6G08_MD_938_proj[2],pos=4,label="938",col="green")

x6G08_MD_939_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_939.txt")
x6G08_MD_939_proj<-project.pca(x6G08_MD_939_raw,pcadata)
points(x6G08_MD_939_proj[3],x6G08_MD_939_proj[2],pch=20)
text(x6G08_MD_939_proj[3],x6G08_MD_939_proj[2],pos=4,label="939",col="red")

x6G08_MD_940_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_940.txt")
x6G08_MD_940_proj<-project.pca(x6G08_MD_940_raw,pcadata)
points(x6G08_MD_940_proj[3],x6G08_MD_940_proj[2],pch=20)
text(x6G08_MD_940_proj[3],x6G08_MD_940_proj[2],pos=4,label="940",col="red")

x6G08_MD_941_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_941.txt")
x6G08_MD_941_proj<-project.pca(x6G08_MD_941_raw,pcadata)
points(x6G08_MD_941_proj[3],x6G08_MD_941_proj[2],pch=20)
text(x6G08_MD_941_proj[3],x6G08_MD_941_proj[2],pos=4,label="941",col="red")

x6G08_MD_942_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_942.txt")
x6G08_MD_942_proj<-project.pca(x6G08_MD_942_raw,pcadata)
points(x6G08_MD_942_proj[3],x6G08_MD_942_proj[2],pch=20)
text(x6G08_MD_942_proj[3],x6G08_MD_942_proj[2],pos=4,label="942",col="red")

x6G08_MD_943_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_943.txt")
x6G08_MD_943_proj<-project.pca(x6G08_MD_943_raw,pcadata)
points(x6G08_MD_943_proj[3],x6G08_MD_943_proj[2],pch=20)
text(x6G08_MD_943_proj[3],x6G08_MD_943_proj[2],pos=4,label="943",col="red")

x6G08_MD_944_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_944.txt")
x6G08_MD_944_proj<-project.pca(x6G08_MD_944_raw,pcadata)
points(x6G08_MD_944_proj[3],x6G08_MD_944_proj[2],pch=20)
text(x6G08_MD_944_proj[3],x6G08_MD_944_proj[2],pos=4,label="944",col="red")

x6G08_MD_945_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_945.txt")
x6G08_MD_945_proj<-project.pca(x6G08_MD_945_raw,pcadata)
points(x6G08_MD_945_proj[3],x6G08_MD_945_proj[2],pch=20)
text(x6G08_MD_945_proj[3],x6G08_MD_945_proj[2],pos=4,label="945",col="red")

x6G08_MD_946_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_946.txt")
x6G08_MD_946_proj<-project.pca(x6G08_MD_946_raw,pcadata)
points(x6G08_MD_946_proj[3],x6G08_MD_946_proj[2],pch=20)
text(x6G08_MD_946_proj[3],x6G08_MD_946_proj[2],pos=4,label="946",col="red")

x6G08_MD_947_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_947.txt")
x6G08_MD_947_proj<-project.pca(x6G08_MD_947_raw,pcadata)
points(x6G08_MD_947_proj[3],x6G08_MD_947_proj[2],pch=20)
text(x6G08_MD_947_proj[3],x6G08_MD_947_proj[2],pos=4,label="947",col="red")

x6G08_MD_948_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_948.txt")
x6G08_MD_948_proj<-project.pca(x6G08_MD_948_raw,pcadata)
points(x6G08_MD_948_proj[3],x6G08_MD_948_proj[2],pch=20)
text(x6G08_MD_948_proj[3],x6G08_MD_948_proj[2],pos=4,label="948",col="red")

x6G08_MD_949_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_949.txt")
x6G08_MD_949_proj<-project.pca(x6G08_MD_949_raw,pcadata)
points(x6G08_MD_949_proj[3],x6G08_MD_949_proj[2],pch=20)
text(x6G08_MD_949_proj[3],x6G08_MD_949_proj[2],pos=4,label="949",col="red")

x6G08_MD_950_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_950.txt")
x6G08_MD_950_proj<-project.pca(x6G08_MD_950_raw,pcadata)
points(x6G08_MD_950_proj[3],x6G08_MD_950_proj[2],pch=20)
text(x6G08_MD_950_proj[3],x6G08_MD_950_proj[2],pos=4,label="950",col="red")

x6G08_MD_951_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_951.txt")
x6G08_MD_951_proj<-project.pca(x6G08_MD_951_raw,pcadata)
points(x6G08_MD_951_proj[3],x6G08_MD_951_proj[2],pch=20)
text(x6G08_MD_951_proj[3],x6G08_MD_951_proj[2],pos=4,label="951",col="red")

x6G08_MD_952_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_952.txt")
x6G08_MD_952_proj<-project.pca(x6G08_MD_952_raw,pcadata)
points(x6G08_MD_952_proj[3],x6G08_MD_952_proj[2],pch=20)
text(x6G08_MD_952_proj[3],x6G08_MD_952_proj[2],pos=4,label="952",col="red")

x6G08_MD_953_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_953.txt")
x6G08_MD_953_proj<-project.pca(x6G08_MD_953_raw,pcadata)
points(x6G08_MD_953_proj[3],x6G08_MD_953_proj[2],pch=20)
text(x6G08_MD_953_proj[3],x6G08_MD_953_proj[2],pos=4,label="953",col="red")

x6G08_MD_954_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_954.txt")
x6G08_MD_954_proj<-project.pca(x6G08_MD_954_raw,pcadata)
points(x6G08_MD_954_proj[3],x6G08_MD_954_proj[2],pch=20)
text(x6G08_MD_954_proj[3],x6G08_MD_954_proj[2],pos=4,label="954",col="red")

x6G08_MD_955_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_955.txt")
x6G08_MD_955_proj<-project.pca(x6G08_MD_955_raw,pcadata)
points(x6G08_MD_955_proj[3],x6G08_MD_955_proj[2],pch=20)
text(x6G08_MD_955_proj[3],x6G08_MD_955_proj[2],pos=4,label="955",col="red")

x6G08_MD_956_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_956.txt")
x6G08_MD_956_proj<-project.pca(x6G08_MD_956_raw,pcadata)
points(x6G08_MD_956_proj[3],x6G08_MD_956_proj[2],pch=20)
text(x6G08_MD_956_proj[3],x6G08_MD_956_proj[2],pos=4,label="956",col="red")

x6G08_MD_957_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_957.txt")
x6G08_MD_957_proj<-project.pca(x6G08_MD_957_raw,pcadata)
points(x6G08_MD_957_proj[3],x6G08_MD_957_proj[2],pch=20)
text(x6G08_MD_957_proj[3],x6G08_MD_957_proj[2],pos=4,label="957",col="red")

x6G08_MD_958_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_958.txt")
x6G08_MD_958_proj<-project.pca(x6G08_MD_958_raw,pcadata)
points(x6G08_MD_958_proj[3],x6G08_MD_958_proj[2],pch=20)
text(x6G08_MD_958_proj[3],x6G08_MD_958_proj[2],pos=4,label="958",col="red")

x6G08_MD_959_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_959.txt")
x6G08_MD_959_proj<-project.pca(x6G08_MD_959_raw,pcadata)
points(x6G08_MD_959_proj[3],x6G08_MD_959_proj[2],pch=20)
text(x6G08_MD_959_proj[3],x6G08_MD_959_proj[2],pos=4,label="959",col="red")

x6G08_MD_960_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_960.txt")
x6G08_MD_960_proj<-project.pca(x6G08_MD_960_raw,pcadata)
points(x6G08_MD_960_proj[3],x6G08_MD_960_proj[2],pch=20)
text(x6G08_MD_960_proj[3],x6G08_MD_960_proj[2],pos=4,label="960",col="red")

x6G08_MD_961_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_961.txt")
x6G08_MD_961_proj<-project.pca(x6G08_MD_961_raw,pcadata)
points(x6G08_MD_961_proj[3],x6G08_MD_961_proj[2],pch=20)
text(x6G08_MD_961_proj[3],x6G08_MD_961_proj[2],pos=4,label="961",col="red")

x6G08_MD_962_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_962.txt")
x6G08_MD_962_proj<-project.pca(x6G08_MD_962_raw,pcadata)
points(x6G08_MD_962_proj[3],x6G08_MD_962_proj[2],pch=20)
text(x6G08_MD_962_proj[3],x6G08_MD_962_proj[2],pos=4,label="962",col="red")

x6G08_MD_963_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_963.txt")
x6G08_MD_963_proj<-project.pca(x6G08_MD_963_raw,pcadata)
points(x6G08_MD_963_proj[3],x6G08_MD_963_proj[2],pch=20)
text(x6G08_MD_963_proj[3],x6G08_MD_963_proj[2],pos=4,label="963",col="red")

x6G08_MD_964_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_964.txt")
x6G08_MD_964_proj<-project.pca(x6G08_MD_964_raw,pcadata)
points(x6G08_MD_964_proj[3],x6G08_MD_964_proj[2],pch=20)
text(x6G08_MD_964_proj[3],x6G08_MD_964_proj[2],pos=4,label="964",col="red")

x6G08_MD_965_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_965.txt")
x6G08_MD_965_proj<-project.pca(x6G08_MD_965_raw,pcadata)
points(x6G08_MD_965_proj[3],x6G08_MD_965_proj[2],pch=20)
text(x6G08_MD_965_proj[3],x6G08_MD_965_proj[2],pos=4,label="965",col="red")

x6G08_MD_966_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_966.txt")
x6G08_MD_966_proj<-project.pca(x6G08_MD_966_raw,pcadata)
points(x6G08_MD_966_proj[3],x6G08_MD_966_proj[2],pch=20)
text(x6G08_MD_966_proj[3],x6G08_MD_966_proj[2],pos=4,label="966",col="red")

x6G08_MD_967_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_967.txt")
x6G08_MD_967_proj<-project.pca(x6G08_MD_967_raw,pcadata)
points(x6G08_MD_967_proj[3],x6G08_MD_967_proj[2],pch=20)
text(x6G08_MD_967_proj[3],x6G08_MD_967_proj[2],pos=4,label="967",col="red")

x6G08_MD_968_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_968.txt")
x6G08_MD_968_proj<-project.pca(x6G08_MD_968_raw,pcadata)
points(x6G08_MD_968_proj[3],x6G08_MD_968_proj[2],pch=20)
text(x6G08_MD_968_proj[3],x6G08_MD_968_proj[2],pos=4,label="968",col="red")

x6G08_MD_969_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_969.txt")
x6G08_MD_969_proj<-project.pca(x6G08_MD_969_raw,pcadata)
points(x6G08_MD_969_proj[3],x6G08_MD_969_proj[2],pch=20)
text(x6G08_MD_969_proj[3],x6G08_MD_969_proj[2],pos=4,label="969",col="red")

x6G08_MD_970_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_970.txt")
x6G08_MD_970_proj<-project.pca(x6G08_MD_970_raw,pcadata)
points(x6G08_MD_970_proj[3],x6G08_MD_970_proj[2],pch=20)
text(x6G08_MD_970_proj[3],x6G08_MD_970_proj[2],pos=4,label="970",col="red")

x6G08_MD_971_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_971.txt")
x6G08_MD_971_proj<-project.pca(x6G08_MD_971_raw,pcadata)
points(x6G08_MD_971_proj[3],x6G08_MD_971_proj[2],pch=20)
text(x6G08_MD_971_proj[3],x6G08_MD_971_proj[2],pos=4,label="971",col="red")

x6G08_MD_972_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_972.txt")
x6G08_MD_972_proj<-project.pca(x6G08_MD_972_raw,pcadata)
points(x6G08_MD_972_proj[3],x6G08_MD_972_proj[2],pch=20)
text(x6G08_MD_972_proj[3],x6G08_MD_972_proj[2],pos=4,label="972",col="red")

x6G08_MD_973_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_973.txt")
x6G08_MD_973_proj<-project.pca(x6G08_MD_973_raw,pcadata)
points(x6G08_MD_973_proj[3],x6G08_MD_973_proj[2],pch=20)
text(x6G08_MD_973_proj[3],x6G08_MD_973_proj[2],pos=4,label="973",col="red")

x6G08_MD_974_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_974.txt")
x6G08_MD_974_proj<-project.pca(x6G08_MD_974_raw,pcadata)
points(x6G08_MD_974_proj[3],x6G08_MD_974_proj[2],pch=20)
text(x6G08_MD_974_proj[3],x6G08_MD_974_proj[2],pos=4,label="974",col="red")

x6G08_MD_975_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_975.txt")
x6G08_MD_975_proj<-project.pca(x6G08_MD_975_raw,pcadata)
points(x6G08_MD_975_proj[3],x6G08_MD_975_proj[2],pch=20)
text(x6G08_MD_975_proj[3],x6G08_MD_975_proj[2],pos=4,label="975",col="red")

x6G08_MD_976_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_976.txt")
x6G08_MD_976_proj<-project.pca(x6G08_MD_976_raw,pcadata)
points(x6G08_MD_976_proj[3],x6G08_MD_976_proj[2],pch=20)
text(x6G08_MD_976_proj[3],x6G08_MD_976_proj[2],pos=4,label="976",col="red")

x6G08_MD_977_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_977.txt")
x6G08_MD_977_proj<-project.pca(x6G08_MD_977_raw,pcadata)
points(x6G08_MD_977_proj[3],x6G08_MD_977_proj[2],pch=20)
text(x6G08_MD_977_proj[3],x6G08_MD_977_proj[2],pos=4,label="977",col="red")

x6G08_MD_978_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_978.txt")
x6G08_MD_978_proj<-project.pca(x6G08_MD_978_raw,pcadata)
points(x6G08_MD_978_proj[3],x6G08_MD_978_proj[2],pch=20)
text(x6G08_MD_978_proj[3],x6G08_MD_978_proj[2],pos=4,label="978",col="red")

x6G08_MD_979_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_979.txt")
x6G08_MD_979_proj<-project.pca(x6G08_MD_979_raw,pcadata)
points(x6G08_MD_979_proj[3],x6G08_MD_979_proj[2],pch=20)
text(x6G08_MD_979_proj[3],x6G08_MD_979_proj[2],pos=4,label="979",col="red")

x6G08_MD_980_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_980.txt")
x6G08_MD_980_proj<-project.pca(x6G08_MD_980_raw,pcadata)
points(x6G08_MD_980_proj[3],x6G08_MD_980_proj[2],pch=20)
text(x6G08_MD_980_proj[3],x6G08_MD_980_proj[2],pos=4,label="980",col="green")

x6G08_MD_981_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_981.txt")
x6G08_MD_981_proj<-project.pca(x6G08_MD_981_raw,pcadata)
points(x6G08_MD_981_proj[3],x6G08_MD_981_proj[2],pch=20)
text(x6G08_MD_981_proj[3],x6G08_MD_981_proj[2],pos=4,label="981",col="red")

x6G08_MD_982_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_982.txt")
x6G08_MD_982_proj<-project.pca(x6G08_MD_982_raw,pcadata)
points(x6G08_MD_982_proj[3],x6G08_MD_982_proj[2],pch=20)
text(x6G08_MD_982_proj[3],x6G08_MD_982_proj[2],pos=4,label="982",col="red")

x6G08_MD_983_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_983.txt")
x6G08_MD_983_proj<-project.pca(x6G08_MD_983_raw,pcadata)
points(x6G08_MD_983_proj[3],x6G08_MD_983_proj[2],pch=20)
text(x6G08_MD_983_proj[3],x6G08_MD_983_proj[2],pos=4,label="983",col="green")

x6G08_MD_984_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_984.txt")
x6G08_MD_984_proj<-project.pca(x6G08_MD_984_raw,pcadata)
points(x6G08_MD_984_proj[3],x6G08_MD_984_proj[2],pch=20)
text(x6G08_MD_984_proj[3],x6G08_MD_984_proj[2],pos=4,label="984",col="green")

x6G08_MD_985_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_985.txt")
x6G08_MD_985_proj<-project.pca(x6G08_MD_985_raw,pcadata)
points(x6G08_MD_985_proj[3],x6G08_MD_985_proj[2],pch=20)
text(x6G08_MD_985_proj[3],x6G08_MD_985_proj[2],pos=4,label="985",col="red")

x6G08_MD_986_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_986.txt")
x6G08_MD_986_proj<-project.pca(x6G08_MD_986_raw,pcadata)
points(x6G08_MD_986_proj[3],x6G08_MD_986_proj[2],pch=20)
text(x6G08_MD_986_proj[3],x6G08_MD_986_proj[2],pos=4,label="986",col="red")

x6G08_MD_987_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_987.txt")
x6G08_MD_987_proj<-project.pca(x6G08_MD_987_raw,pcadata)
points(x6G08_MD_987_proj[3],x6G08_MD_987_proj[2],pch=20)
text(x6G08_MD_987_proj[3],x6G08_MD_987_proj[2],pos=4,label="987",col="green")

x6G08_MD_988_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_988.txt")
x6G08_MD_988_proj<-project.pca(x6G08_MD_988_raw,pcadata)
points(x6G08_MD_988_proj[3],x6G08_MD_988_proj[2],pch=20)
text(x6G08_MD_988_proj[3],x6G08_MD_988_proj[2],pos=4,label="988",col="green")

x6G08_MD_989_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_989.txt")
x6G08_MD_989_proj<-project.pca(x6G08_MD_989_raw,pcadata)
points(x6G08_MD_989_proj[3],x6G08_MD_989_proj[2],pch=20)
text(x6G08_MD_989_proj[3],x6G08_MD_989_proj[2],pos=4,label="989",col="green")

x6G08_MD_990_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_990.txt")
x6G08_MD_990_proj<-project.pca(x6G08_MD_990_raw,pcadata)
points(x6G08_MD_990_proj[3],x6G08_MD_990_proj[2],pch=20)
text(x6G08_MD_990_proj[3],x6G08_MD_990_proj[2],pos=4,label="990",col="red")

x6G08_MD_991_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_991.txt")
x6G08_MD_991_proj<-project.pca(x6G08_MD_991_raw,pcadata)
points(x6G08_MD_991_proj[3],x6G08_MD_991_proj[2],pch=20)
text(x6G08_MD_991_proj[3],x6G08_MD_991_proj[2],pos=4,label="991",col="red")

x6G08_MD_992_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_992.txt")
x6G08_MD_992_proj<-project.pca(x6G08_MD_992_raw,pcadata)
points(x6G08_MD_992_proj[3],x6G08_MD_992_proj[2],pch=20)
text(x6G08_MD_992_proj[3],x6G08_MD_992_proj[2],pos=4,label="992",col="green")

x6G08_MD_993_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_993.txt")
x6G08_MD_993_proj<-project.pca(x6G08_MD_993_raw,pcadata)
points(x6G08_MD_993_proj[3],x6G08_MD_993_proj[2],pch=20)
text(x6G08_MD_993_proj[3],x6G08_MD_993_proj[2],pos=4,label="993",col="green")

x6G08_MD_994_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_994.txt")
x6G08_MD_994_proj<-project.pca(x6G08_MD_994_raw,pcadata)
points(x6G08_MD_994_proj[3],x6G08_MD_994_proj[2],pch=20)
text(x6G08_MD_994_proj[3],x6G08_MD_994_proj[2],pos=4,label="994",col="red")

x6G08_MD_995_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_995.txt")
x6G08_MD_995_proj<-project.pca(x6G08_MD_995_raw,pcadata)
points(x6G08_MD_995_proj[3],x6G08_MD_995_proj[2],pch=20)
text(x6G08_MD_995_proj[3],x6G08_MD_995_proj[2],pos=4,label="995",col="green")

x6G08_MD_996_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_996.txt")
x6G08_MD_996_proj<-project.pca(x6G08_MD_996_raw,pcadata)
points(x6G08_MD_996_proj[3],x6G08_MD_996_proj[2],pch=20)
text(x6G08_MD_996_proj[3],x6G08_MD_996_proj[2],pos=4,label="996",col="green")

x6G08_MD_997_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_997.txt")
x6G08_MD_997_proj<-project.pca(x6G08_MD_997_raw,pcadata)
points(x6G08_MD_997_proj[3],x6G08_MD_997_proj[2],pch=20)
text(x6G08_MD_997_proj[3],x6G08_MD_997_proj[2],pos=4,label="997",col="green")

x6G08_MD_998_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_998.txt")
x6G08_MD_998_proj<-project.pca(x6G08_MD_998_raw,pcadata)
points(x6G08_MD_998_proj[3],x6G08_MD_998_proj[2],pch=20)
text(x6G08_MD_998_proj[3],x6G08_MD_998_proj[2],pos=4,label="998",col="blue")

x6G08_MD_999_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_999.txt")
x6G08_MD_999_proj<-project.pca(x6G08_MD_999_raw,pcadata)
points(x6G08_MD_999_proj[3],x6G08_MD_999_proj[2],pch=20)
text(x6G08_MD_999_proj[3],x6G08_MD_999_proj[2],pos=4,label="999",col="green")

x6G08_MD_1000_raw<-scan("/local/scratch/ejs1e14/MD_6G08/MD_PCA_ONLY/6G08_MD_1000.txt")
x6G08_MD_1000_proj<-project.pca(x6G08_MD_1000_raw,pcadata)
points(x6G08_MD_1000_proj[3],x6G08_MD_1000_proj[2],pch=20)
text(x6G08_MD_1000_proj[3],x6G08_MD_1000_proj[2],pos=4,label="1000",col="red")

q() 
