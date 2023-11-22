library(dplyr)
library(data.table)

# Arg 1 is path of the new dataframe
# Arg 2 is path of the old dataframe (with ensembl_id as first col)
# Arg 3 is output of new one
args <- commandArgs(trailingOnly = TRUE)

# Read in the new dataframe
new_df <- fread(args[1], header = FALSE, sep = "\t", data.table = FALSE)

# Read in old dataframe
old_df <- fread(args[2], header = FALSE, sep = "\t", data.table = FALSE)

# Merge the two dataframes by cols V2 and V3, then replace first colun with V1 from the old dataframe
merged_df <- left_join(new_df, old_df, by = c("V1", "V3"), suffix=c('.delete', '')) %>%
        select(V1, V2, V3)



# Write out the merged dataframe to .tsv.gz file 
gz1 <- gzfile(args[3], "w")
write.table(merged_df, file = gz1, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
close(gz1)

