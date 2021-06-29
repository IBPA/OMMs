# Constants
MONO_COLUMNS <- c('Glucose', 'Galactose', 'Fructose', 'Xylose', 'Arabinose', 'Fucose', 'Rhamnose', 'GlcA', 'GalA', 'GlcNAc', 'GalNAc', 'Mannose', 'Allose', 'Ribose')

load_data <- function(data_file, metadata_file){
  # A) Load data
  df <- read.table(data_file, sep = ",", quote = '"', header = TRUE)
  df_data <- df[,MONO_COLUMNS]
  df_data[is.na(df_data)] <- 0
  df[,MONO_COLUMNS] <- df_data
  colnames(df)[4] <- "uid"
  
  # B) Load water content metadata
  df_meta <- read_excel(metadata_file)
  df_meta <- df_meta[,c("Unique Identifier Code","% water content")]
  colnames(df_meta) <- c("uid", "water")
  df_meta$uid <- as.numeric(str_remove(df_meta$uid, "^0+"))
  
  uids_invalid_water_content <- unlist(df_meta[is.na(df_meta$water) | df_meta$water<0 | df_meta$water>100, c("uid") ], use.names = FALSE)
  if (length(uids_invalid_water_content)>0){
    uids_str <- paste(uids_invalid_water_content, collapse = ",")
    warning(sprintf("Some foods have invalid water content (not between 0-100). Corresponding records with missing values will be removed and the water content of the rest will be clamped to be in range. The uids are: \n%s", uids_str))
  }
  
  df_meta$dry <- (100-df_meta$water)/100
  df_meta$dry <- pmin(pmax(df_meta$dry, 0.0), 1.0) # Clamp the values
  df_meta <- df_meta[!is.na(df_meta$dry),] # Remove missing value records
  
  df_meta$water <- NULL
  
  # C) Incorporate dry content % of each food into the corresponding glycan content vector
  df_comb <- merge(df, df_meta, by = "uid")
  df_comb[,MONO_COLUMNS] <- df_comb[,MONO_COLUMNS] * df_comb$dry
  #  ***** Skip this for now, due to issue in water % (not between 0-100)
  
  return(df_comb[,c("uid", MONO_COLUMNS)])
}

minmax_normalize <- function(df_data){
  scale_minmax <- function(x){
    (x- min(x)) /(max(x)-min(x))
  }
  df_data_scaled <- as.data.frame(lapply(df_data, scale_minmax))
  df_data_scaled[is.na(df_data_scaled)] <- 0
  
  return(df_data_scaled)
}

save_MMs <- function(df_food_proportions, filename) {
  # Write the non-zero proportions to the file for each meal.
  nz_proportions <- apply(df_food_proportions, 1, FUN = function(x){x[x!=0]})
  
  lines <- c()
  for (meal_idx in 1:length(nz_proportions)) {
    cur_meal_props <- nz_proportions[[meal_idx]]
    line_str <- paste(sprintf("%s:%.2f", names(cur_meal_props), unname(cur_meal_props)), collapse = ",")
    lines <- c(lines, line_str)
  }
  writeLines(lines, filename)
}

get_friendly_food_names <- function(uids, data_file){
  df <- read.table(data_file, sep = ",", quote = '"', header = TRUE)
  df_names <- df[df$Unique.Identifier.Code %in% uids, c("Unique.Identifier.Code", "Simple.name")]
  rownames(df_names) <- df_names$Unique.Identifier.Code
  return(unname(df_names[uids, c("Simple.name")]))
}

load_MMs <- function(filename) {
  all_lines <- str_split(readLines(filename), ",")
  fids <- unique(unlist(unname(data.frame(str_split(unlist(all_lines), ":"))[1,])))
  num_meals <- length(all_lines)
  
  MMs <- matrix(0, nrow = num_meals, ncol = length(fids), dimnames = list(1:num_meals, fids))
  for (i in 1:num_meals){
    line_cur <- all_lines[[i]]
    df_cur <- data.frame(str_split(line_cur, ":"))
    fids_cur <- unlist(unname(df_cur[1,]))
    proportions_cur <- unname(df_cur[2,])
    MMs[i,fids_cur] <- as.numeric(unlist(proportions_cur))
  }
  
  return(MMs)
}

get_information_content <- function(df_vectors) {
  # Remove const columns
  non_const_columns <- MONO_COLUMNS[apply(df_vectors, 2, var) != 0]
  df_data <- df_vectors[,non_const_columns]
  df_data <- scale(df_data, scale = FALSE)
  
  # Calculate Von Neumann entropy (i.e. entropy of eigen values of the covariance matrix)
  eigen_values <- abs(eigen(cov(df_data))$values)
  eigen_values <- eigen_values/sum(eigen_values)
  
  von_neumann_entropy <- -sum(unlist(lapply(eigen_values, FUN = function(x){x*log(x)})))
  
  return(von_neumann_entropy)
}
