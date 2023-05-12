library(data.table)
df <-
  list.files(path = "/Users/waelg/OneDrive/Bureau/EPFL_4_2/Semester_project/data/Noordin Top CSV/CSV/", pattern = "*.csv") %>%
  map_df(~fread(.))
df

# setwd("C:/Users/waelg/OneDrive/Bureau/EPFL_4_2/Semester_project/")
# 
# list_csv_files <- list.files(path = "/Users/waelg/OneDrive/Bureau/EPFL_4_2/Semester_project/data/NoordinTopCSV/CSV/")
# df2 = do.call(rbind, lapply(list_csv_files, function(x) read.csv(x, stringsAsFactors = FALSE, sep = ",")))
# df2
