library(devtools)

barcodes <- read.csv("data/raw_feature_bc_matrix/barcodes.tsv", sep="\t", 
                     header=FALSE, col.names=c("Barcodes"))
features <- read.csv("data/raw_feature_bc_matrix/features.tsv", sep="\t", 
                     header=FALSE, col.names=c("Feature_ID", "Feature_name", 
                                               "Feature_type"))
library("Matrix")
counts <- readMM(file="data/raw_feature_bc_matrix/matrix.mtx")
# sce <- SingleCellExperiment::SingleCellExperiment(rowData=features, 
#                                                   colData=barcodes, 
#                                                   assays=c(counts=counts))

tissue.positions <- read.csv("data/spatial/tissue_positions_list.csv", 
                             header=FALSE, col.names=c("Barcodes", "in_tissue", 
                                                       "array_row", "array_col",
                                                       "pxl_col_in_fullres", "pxl_row_in_fullres"))


library("rjson")
scalefactors <- fromJSON(file="data/spatial/scalefactors_json.json")



# sce <- SingleCellExperiment::SingleCellExperiment(rowData=features, 
#                                                   colData=barcodes, 
#                                                   assays=c(counts=counts))

# 
# sum(barcodes[,1] %in% tissue.positions[,1])
# dim(barcodes)
# dim(tissue.positions)
# 
# sum(duplicated(tissue.positions[,c(3,4)]))
# 
# 
# sce <- SingleCellExperiment::SingleCellExperiment(rowData=features, colData=barcodes,
#                                                   assays=c(counts=counts))
# 

load_all()
ve <- VisiumExperiment(rowData=features, colData=barcodes, 
                           assays=c(counts=counts), 
                           spatialCoords=tissue.positions,
                           scaleFactors=scalefactors)
ve
spatialCoords(ve)
isInTissue(ve)


hres <- "data/spatial/tissue_hires_image.png"
# 
# png(filename=hres)
# par(bg = "transparent", usr = c(0, 51, 0, 451))
# 
# hexpolygon(20,20, dx = 5, density = 25, col = 2, lwd = 1.5)
# 
# 
# dev.off()
# 
 
library(png)
m <- readPNG(hres)
rimg <- as.raster(m)
library(lattice)
library(grid)
grid.raster(rimg, x=0, y=0, just="bottom", default.units = "native")
panel.fill(col=rgb(1,1,1,alpha=0.3))

hbin<-hexbin(20,20,xbins=50,IDs=TRUE)
mtrans<-hexTapply(hbin,dat$Points,sum,na.rm=TRUE)
cols <- rainbow(4)
grid.hexagons(hbin, style='lattice',
              ,minarea=0.1,maxarea=50,colorcut=c(0,.6,1),
              border=NA,
              pen=cols[mtrans+1])














