
library(tidyverse)
library(lubridate)

all_txt_files <- system('locate *', intern = TRUE)
basenames <- basename(all_txt_files)
file_info <- file.info(all_txt_files)

edit_dates <- date(file_info$atime)
known_edit_date <- ymd('2021/02/10')

all_txt_files[which(edit_dates == known_edit_date)]

