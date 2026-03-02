rm(list = ls())
options(stringsAsFactors = F)

# settings ----------------------------------------------------------------

# arguments
argsVal <- commandArgs(trailingOnly = T)

wdir <- argsVal[1]

setwd(wdir)

df <- read.table("all.geno.txt")

### colnames ind, ALThomo, ALThet, REFhomo, missings
colnames(df) <- c("sample","ALThomo","ALThet","REFhomo","missing")

for (i in 1:nrow(df)) {
  
  if (df[i,"REFhomo"] < 10000000 | df[i,"missing"] > 1000000) {
    
    df[i,"decision"] <- "discard"
    
    } else {
      
      df[i,"decision"] <- "keep"
    } }

df_keep <- df[df$decision == "keep",]

df_discard <- df[df$decision == "discard",]

write.table(df_keep,file = "df_keep.txt",append = F,quote = F,sep = "\t",col.names = F,row.names = F)

write.table(df_discard,file = "df_discard.txt",append = F,quote = F,sep = "\t",col.names = F,row.names = F)
