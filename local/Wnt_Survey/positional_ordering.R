library(URD)
library(gplots)
library(rstudioapi)

#set working directory to folder containing this script
setwd(dirname(getActiveDocumentContext()$path))

U.Endo <- readRDS("Hydra_URD_Endoderm.rds")

tentacle.spline <- geneSmoothFit(U.Endo, method = "spline", pseudotime = "pseudotime",
                                 cells = cellsInCluster(U.Endo, "segment", c("47", "46")),
                                 moving.window = 1, cells.per.window = 5, spar = 0.9, genes = U.Endo@count.data@Dimnames[[1]])


hypo.spline <- geneSmoothFit(U.Endo, method = "spline", pseudotime = "pseudotime",
                                 cells = cellsInCluster(U.Endo, "segment", c("47", "11")),
                                 moving.window = 1, cells.per.window = 5, spar = 0.9, genes = U.Endo@count.data@Dimnames[[1]])

pt.crop <- as.numeric(unlist(U.Endo@tree$segment.pseudotime.limits)[1])

foot.only.spline <- cropSmoothFit(tentacle.spline, pt.max = pt.crop)
hypostome.only.spline <- cropSmoothFit(hypo.spline, pt.min = pt.crop)
tentacle.only.spline <- cropSmoothFit(tentacle.spline, pt.min = pt.crop)

splines <- list(foot.only.spline, hypostome.only.spline, tentacle.only.spline)
names(splines) <- c("Foot/Body","Hypostome","Tentacle")

cols <- (scales::gradient_n_pal(RColorBrewer::brewer.pal(9, "YlOrRd")))(seq(0, 1, length.out = 50))

splines.reduced <- lapply(names(splines), function(n){
  s <- splines[[n]]
  l <- tolower(substr(n, 1, 1))
  colnames(s$mean.expression) <- paste0(1, as.character(round(as.numeric(colnames(s$mean.expression)),
                                                              digits = 2)))
  s$mean.expression <- matrixReduce(s$mean.expression)
  colnames(s$mean.smooth) <- paste0(1, as.character(round(as.numeric(colnames(s$mean.smooth)),
                                                          digits = 2)))
  s$mean.smooth <- matrixReduce(s$mean.smooth)
  return(s)
})

ms <- do.call("cbind", lapply(splines.reduced, function(s) s$mean.smooth))

ms.endo <- ms


U.Ecto <- readRDS("Hydra_URD_Ectoderm.rds")

basal.disc.cells <- cellsInCluster(U.Ecto, "res.1.5", "3")

body.column.and.battery.cells <- setdiff(cellsInCluster(U.Ecto, "segment", c("16", "26")), basal.disc.cells)

hypo.trajectory.cells <- c(
  cellsInCluster(U.Ecto, "segment", "25"),
  cellsInCluster(U.Ecto, "segment", "26")[order(U.Ecto@pseudotime[cellsInCluster(U.Ecto, "segment", "26"), "pseudotime"], decreasing = T)[1:250]] # Include the 250 'oldest' cells from before the branchpoint
) 

geneList <- U.Ecto@count.data@Dimnames[[1]]

# Calculate smoothed spline fits of the expression within each population
basal.disc.spline.single <- geneSmoothFit(U.Ecto, method="spline", 
                                          pseudotime="pseudotime", cells = basal.disc.cells, 
                                          genes = geneList, moving.window = 1, 
                                          cells.per.window = 5, spar=0.875) # Basal disc only

battery.spline.single <- geneSmoothFit(U.Ecto, method="spline", pseudotime="pseudotime", 
                                       cells = body.column.and.battery.cells, genes = geneList, 
                                       moving.window = 1, cells.per.window = 5, spar=0.875) # Body column and tentacles (which contain battery cells)

head.spline.single <- geneSmoothFit(U.Ecto, method="spline", pseudotime="pseudotime", 
                                    cells = hypo.trajectory.cells, genes = geneList, moving.window = 1, 
                                    cells.per.window = 5, spar=0.875) # Hypostome

## Combine pieces of spline fits into a single list for multi-plotting
# Identify the pseudotime of the branchpoint / breakpoint between basal disc and body column
pt.crop <- as.numeric(unlist(U.Ecto@tree$segment.pseudotime.limits)[1])
basal.crop <- max(U.Ecto@pseudotime[cellsInCluster(U.Ecto, "res.1.5", "3"), "pseudotime"])

# Crop according to the pseudotime of the branchpoint
body.only.spline.single <- cropSmoothFit(battery.spline.single, pt.max=pt.crop)
battery.only.spline.single <- cropSmoothFit(battery.spline.single, pt.min=pt.crop)
head.only.spline.single <- cropSmoothFit(head.spline.single, pt.min=pt.crop)

# Combine into a list; Names in the plots are determined by names of the smooth objects in the list
splines <- list(basal.disc.spline.single, body.only.spline.single, head.only.spline.single, battery.only.spline.single)
names(splines) <- c("Basal Disc", "Body Column", "Hypostome", "Tentacle")

splines.reduced <- lapply(names(splines), function(n) {
  s <- splines[[n]]
  l <- tolower(substr(n, 1, 2))
  colnames(s$mean.expression) <- paste0(l, as.character(round(as.numeric(colnames(s$mean.expression)), digits=2)))
  s$mean.expression <- matrixReduce(s$mean.expression)
  colnames(s$mean.smooth) <- paste0(l, as.character(round(as.numeric(colnames(s$mean.smooth)), digits=2)))
  s$mean.smooth <- matrixReduce(s$mean.smooth)
  return(s)
})

# Combine the four splines into a single table & rescale maximum across all regions
ms <- do.call("cbind", lapply(splines.reduced, function(s) s$mean.smooth))

ms.ecto <- ms

ms.ecto$ID <- rownames(ms.ecto)

ms.endo$ID <- rownames(ms.endo)

ms.both <- merge(ms.ecto, ms.endo, by = "ID")

ms.both$ID <- gsub("[|].*","",ms.both$ID)

rownames(ms.both) <- ms.both$ID

ms.both$ID <- NULL

ms.both <- as.matrix(ms.both)

save(ms.both, file = "ecto_endo_spline.RData")
