rm(list=ls())

setwd("~/Desktop/zaitlen_lab_desktop/")

counts = read.table("140287_bisulfite.txt", header=T)
counts$REFERENCE = toupper(counts$REFERENCE)

counts_only = counts[, 4:12]

count_bp <- function(bp, df) {
  print(bp)
  df_bp = subset(df, df$REFERENCE==bp)
  df_bp = df_bp[, 2:9]
  return(colSums(df_bp))
}

base_pairs = c("A", "C", "T", "G")
count_list = lapply(base_pairs, count_bp, df=counts_only)


count_df = data.frame(do.call(rbind, count_list))
count_df$bp = base_pairs
count_df = count_df[,c(9, 1, 2, 3, 4, 5, 6, 7, 8)]

counts_rev = count_df[, c(1:5)]
counts_fwd = count_df[, c(1, 6:9)]

rev = cbind(counts_rev[1], prop.table(as.matrix(counts_rev[-1]), margin = 1))
fwd = cbind(counts_fwd[1], prop.table(as.matrix(counts_fwd[-1]), margin = 1))

write.table(x = rev, "reverse.txt", row.names = F)
write.table(x = fwd, "forward.txt", row.names = F)
