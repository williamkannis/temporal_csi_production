#-------------------------------------------------------------------------------
#
#   Predictor preparation
#
#-------------------------------------------------------------------------------

# AUTHOR: William K. Annis

# CREATED: July 13, 2026

# DESCRIPTION: 


# House keeping  ---------------------------------------------------------------
rm(list = ls())

# Packages
library(dplyr)
devtools::load_all("~/Documents/work/R packages/secProd")

# directores
data_dir <- paste0(
  "~/Documents/Work/Everglades post-doc/",
  "Data analysis/Data cleaning/cleaned_data"
)
input_dir <- "input_data"
export_dir <- "prod_data"

# data
len_df <- readRDS(file.path(input_dir,"fslen_imputed_2026-07-09.rds"))
pisc_df <- readRDS(file.path(data_dir,"pisc_cleaned_2026-07-13.rds"))
samp_df <- readRDS(file.path(export_dir,"fs_sample_info_2026-08-12.rds"))
phy_df <- readRDS(file.path(data_dir,"phys_cleaned_2026-07-09.rds"))


# Sampling interval  -----------------------------------------------------------

samp_for <- samp_df %>% 
  left_join(
    phy_df %>% distinct(region,site,wateryear,cum),
    by = join_by(site,cum)
    ) %>% 
  mutate(
    interval_pred = round(interval / 30) * 30,
    interval_pred = case_when(
      interval_pred == 120 ~ 90,
      interval_pred == 150 ~ 180,
      T ~ interval_pred
    )
    )


# site-level predictors --------------------------------------------------------

phy_site <- phy_df %>% 
  select(
    region,
    site,
    wateryear,
    cum,
    depth_ave_30day,
    depth_ave_60day,
    depth_ave_90day,
    depth_ave_180day,
    dsldd,
    lastdaydry,
    peri_vol,
    plt_cov) %>% 
  
# Summarize data from throw level to site level
  group_by(region,site,wateryear,cum) %>% 
  summarize(
    across(everything(),~mean(.x,na.rm=T))
  ) %>% 
  mutate(across(everything(), ~ ifelse(is.nan(.x), NA, .x))) %>% 
  inner_join(samp_for) %>% 
  group_by(region,site) %>% 

# Calculate average predictor values during each sampling interval
  mutate(
    across(
      c(plt_cov,peri_vol,dsldd),
      ~ {
        nxt <- lead(.x,order_by = cum)
        ifelse(
          !is.na(.x) & !is.na(nxt),
          (.x + nxt) / 2,
          coalesce(.x, nxt)
        )
      },
      .names = "{.col}_int"
    ),

# Create column for average depth at end of sampling intera
    across(
      contains("depth"),
      ~lead(.x,order_by = cum),
      .names = "{.col}_lead"
      )
  ) %>% 
  ungroup() 

# Assign depth leads based on sampling interval duration
depth_names <- sapply(
  phy_site$interval_pred, 
  function(x) paste0("depth_ave_",x,"day_lead")
  )
phy_site$depth <- sapply(
  seq_len(nrow(phy_site)),
  function(i) {
    if (!depth_names[i] %in% names(phy_site)) return(NA)
    phy_site[[depth_names[i]]][i]
  }
  )

# remove intermediate variables
phy_site_final <- phy_site %>% 
  select(
    region,
    site,
    wateryear,
    cum,
    depth,
    contains("_int")
  )
summary(phy_site_final)

# TEMPORARY IMPUTE
set.seed(999)
phy_site_final <- group_impute(
  phy_site_final,
  plt_cov_int,
  by = c(region,site,wateryear)
  )
phy_site_final <- group_impute(
  phy_site_final,
  peri_vol_int,
  by = c(region,site,wateryear)
  )
phy_site_final <- phy_site_final %>% select(-impute)


# Site/Year-level predictors ---------------------------------------------------
phy_site_year <- phy_df %>% 
  select(
    region,
    site,
    wateryear,
    waterperiod,
    wet_sum_365day,
    depth_ave_365day,
    # dsldd,
    lastdaydry
    ) %>% 
  filter(waterperiod == 5) %>% 
  group_by(wateryear,region,site) %>% 
  summarize(
    across(everything(),~mean(.x,na.rm=T)),
    .groups = "drop"
  ) %>% 
  left_join(pisc_df) # %>%
  # inner_join(samp_for %>% distinct(region,site,wateryear)) 
summary(phy_site_year)



# TEMPORARY IMPUTE
set.seed(999)
phy_site_year_impute <- group_impute(
  phy_site_year,
  pisc_index,
  by = c(region)
) %>% 
  select(-impute)

# Create lag effects
phy_site_year_lag <- phy_site_year_impute %>% 
  group_by(region,site) %>% 
  mutate(
    across(
      .cols = c(wet_sum_365day,depth_ave_365day,lastdaydry,pisc_index),
      .fns = ~lag(.x,order_by = wateryear),
      .names = "{.col}_lag"
    ),
  ) %>% 
  ungroup() %>% 
  inner_join(samp_for %>% distinct(region,site,wateryear)) 

# Check for correlation
phy_site_year_lag %>% 
  select(-wateryear, -region, -site,  -waterperiod) %>% 
  cor(use = "complete.obs")

# Site-year PCA  ---------------------------------------------------------------

# Use full pisc and hydrology data sets to conduct PCA on biological and
# hydrology variables
# Prepare data input
pca_input <-phy_site_year_lag %>% 
  select(-wateryear, -region, -site,  -waterperiod) %>% 
  select(!contains("lag"))


# Run pca
pca_out <- vegan::rda(pca_input,scale = T)

# Examine results
summary(pca_out)
pca_out$CA$v
plot(pca_out, type = "n") 
points(pca_out, pch=19, display = "sites") 
text(pca_out, display = "species", col="blue") 

# Extract pcs that explain atleast 75% of variation and create data.frame
pca_result <-vegan::scores(pca_out,choices = c(1,2),display = "sites")
phy_site_year_pca <- cbind(phy_site_year_lag,pca_result)


# Region-year level predictors  ------------------------------------------------
phy_reg_year <- phy_site_year_pca %>% 
  group_by(region,wateryear) %>% 
  summarise(across(
    c(-site,-waterperiod),
    ~mean(.x,na.rm=T)
    ),
    .groups = "drop"
  )

# Check for correlation
phy_reg_year %>% 
  select(-wateryear, -region) %>% 
  cor(use = "complete.obs")


# Year-level predictors  -------------------------------------------------------
phy_year <- phy_site_year_pca %>% 
  group_by(wateryear) %>% 
  summarise(across(
    c(-site,-region,-waterperiod),
    ~mean(.x,na.rm=T)
    ),
    .groups = "drop"
    )

# Check for correlation
phy_year %>% 
  select(-wateryear) %>% 
  cor(use = "complete.obs")


# Export  ----------------------------------------------------------------------
saveRDS(
  phy_site_final,
  file.path(export_dir,paste0("phys_site_predictors_",Sys.Date(),".rds"))
  )
saveRDS(
  phy_site_year_pca, 
  file.path(export_dir,paste0("phys_siteyear_predictors_",Sys.Date(),".rds"))
  )
saveRDS(
  phy_reg_year,
  file.path(export_dir,paste0("phys_regionyear_predictors_",Sys.Date(),".rds"))
)
saveRDS(
  phy_year,
  file.path(export_dir,paste0("phys_year_predictors_",Sys.Date(),".rds"))
)

