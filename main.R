# comm 04/11/2025

set.seed(1234)
libdir <- "renv/library/linux-arch-rolling/R-4.5/x86_64-pc-linux-gnu/"
setwd("~/Documents/communication/")
#
## set local directory for R packages
.libPaths(libdir)
.libPaths()
#

library("devtools")
devtools::install_github('linxihui/NNLM')
devtools::install_github('ZJUFanLab/SpaTalk')
library(SpaTalk)

file_path = "~/Documents/communication/data/"
files = list.files(file_path)

# while loop (forma burra mas logica)
f = list()
x = 0 
count = 1
while(length(files) > x) {
  f[[count]] = files[c(x+1,x+2)]
  x = x + 2
  count = count + 1
}

library(tibble)
for (i in f) {
  st_data = read.csv(paste0(file_path, i[1]))
  rownames(st_data) = st_data[,1]
  st_data[,1] = NULL
  head(tibble(st_data))
  
  st_meta = read.csv(paste0(file_path, i[2]))
  st_meta$X.1 = NULL
  colnames(st_meta) = c('x', 'y', 'celltype', 'spot')
  head(tibble(st_meta))
  
  # garantia, evita erros de leitura do csv "." ou "-" ex
  colnames(st_data) = st_meta$spot
  
  # revise gene symbols according to the NCNI gene symbols
  st_data <- rev_gene(data = as.matrix(st_data),
                      data_type = "count",
                      species = "Human",
                      geneinfo = geneinfo
  )
  
  st_data <- as.data.frame(as.matrix(st_data))
  
  st_meta$celltype <- factor(st_meta$celltype)
  st_meta$celltype = as.character(st_meta$celltype)
  
  # filtragem de celltypes com poucos spots
  cell_counts <- table(st_meta$celltype)
  valid_cells <- names(cell_counts[cell_counts >= 5])
  st_meta <- st_meta[st_meta$celltype %in% valid_cells, ]
  st_data <- st_data[, st_meta$spot]
  
  if (all(colnames(st_data) == st_meta$spot)) {
    print("spots do metadado e colunas de st_data batem")
  }
  
  # create SpaTalk data
  obj <- createSpaTalk(st_data = as.matrix(st_data),
                       st_meta = st_meta[, c(4, 1, 2)],
                       species = "Human",
                       if_st_is_sc = F,
                       spot_max_cell = 5,
                       celltype = st_meta$celltype
  )
  
  # calculo das vias paralelizando
  library(doParallel)
  registerDoParallel(cores = 4)
  
  obj <- find_lr_path(object = obj, 
                      lrpairs = lrpairs, 
                      pathways = pathways, 
                      if_doParallel = T
  )
  
  dec_cci_all(object = obj,
              if_doParallel = T)
  
  # celltypes <- c("6", "9")
  # 
  # celltypes <- as.character(c(0:13))
  # for (sender in celltypes) {
  #   for (receiver in celltypes) {
  #     message(paste0("Running ", sender, " → ", receiver, " ..."))
  #     
  #     res <- try(
  #       dec_cci(object = obj,
  #               celltype_sender = sender,
  #               celltype_receiver = receiver,
  #               if_doParallel = T,
  #               min_pairs = 0)
  #     )
  #     
  #     if (inherits(res, "try-error")) {
  #       message(paste0("⚠️ Failed for ", sender, " → ", receiver, " — skipping."))
  #     } else {
  #       obj <- res
  #     }
  #   }
  # }
  
  name = paste(strsplit(i[1], "_")[[1]][1], 
               strsplit(i[1], "_")[[1]][2],
               strsplit(i[1], "_")[[1]][3], sep = "_")
  
  save(obj, file = paste0("teste", name, ".RDS"))
}

# Infer all cell-cell communications
#obj <- dec_cci_all(object = obj, if_doParallel = F)



celltypes <- c("6", "9")

for (sender in celltypes) {
  for (receiver in celltypes) {
    message(paste0("Running ", sender, " → ", receiver, " ..."))
  }
}



