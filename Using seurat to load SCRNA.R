####installing DeuratDisk
install.packages("devtools")
devtools:: install_github("mojaveazure/seurat-disk")
install.packages("hdf5r")
library(hdf5r)
install.packages("Seurat")
#loading library
library(Seurat)
library(SeuratDisk)
####importing the hdf5  10x data
hdf5_data <- Read10X_h5(filename = "C:\\Users\\iyand\\Downloads\\scData\\20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix (2).h5", 
  use.names = TRUE, 
  unique.features = TRUE) ###feature matric barcode
####converting the feature matric barcode to seurat object
seurat_hdf5 <- CreateSeuratObject(counts = hdf5_data, min.cells = 3, min.features = 200)
#####Quality control
###mitochondria read
seurat_hdf5[["percentage.mt"]] <- PercentageFeatureSet(seurat_hdf5, pattern = "^MT-")
View(seurat_hdf5@meta.data) ####viewing my meta data
#####visualization using violine plot
VlnPlot(seurat_hdf5, features =  c("nFeature_RNA", "nCount_RNA", "percentage.mt"), ncol=3)
######
FeatureScatter(seurat_hdf5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
  geom_smooth(method ="lm")
##### filtring process by subsetting
seurat_hdf5 <- subset(seurat_hdf5, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percentage.mt <5)

######Normalization
seurat_hdf5 <- NormalizeData(seurat_hdf5)
##### identify 15 highly variable feature
seurat_hdf5<- FindVariableFeatures(seurat_hdf5, selection.method = "vst", nfeatures = 1000)
top15 <- head(VariableFeatures(seurat_hdf5), 15)

####pploting the variable
plot1 <- VariableFeaturePlot(seurat_hdf5)
LabelPoints(plot = plot1, points = top15, repel = TRUE)

####Data Scaling
seurat_hdf5 <- ScaleData(seurat_hdf5, features = rownames(seurat_hdf5))

####performing Linear dimensionality
seurat_hdf5 <- RunPCA(seurat_hdf5, features = VariableFeatures(object = seurat_hdf5))
#### dimentionality
ElbowPlot(seurat_hdf5)

######clustering
seurat_hdf5 <- FindNeighbors(seurat_hdf5, dims = 1:15)

seurat_hdf5 <- FindClusters(seurat_hdf5, resolution = 0.5)
####MARKER FOR CLUSTER
cluster_markers <- FindAllMarkers(seurat_hdf5, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)



