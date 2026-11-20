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

# Directories
data_dir <- paste0(
  "~/Documents/Work/Everglades post-doc/",
  "Data analysis/Data cleaning/cleaned_data"
)
input_dir <- "input_data"
export_dir <- "prod_data"

# Data
len_df <- readRDS(file.path(input_dir,"fslen_imputed_2026-11-16.rds"))
pisc_df <- readRDS(file.path(data_dir,"pisc_cleaned_2026-11-03.rds"))
samp_df <- readRDS(file.path(export_dir,"fs_sample_info_2026-11-18.rds"))
phy_df <- readRDS(file.path(data_dir,"phys_cleaned_2026-11-03.rds"))


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

# Impute using median values for site and year
miss_meds <- phy_site_final %>% 
  group_by(region,site,wateryear) %>% 
  summarise(
    across(
      .cols = c(plt_cov_int,peri_vol_int),
      .fns = \(x) mean(x, na.rm = T),
      .names = "{.col}_median"
    )
  )

phy_site_impute <- phy_site_final %>% 
  left_join(miss_meds) %>% 
  mutate(
    across(
      .cols = c(plt_cov_int,peri_vol_int),
      .fns = \(x) case_when(
        is.na(x) ~ get(paste0(cur_column(), "_median")),
        T ~ x
      )
    )
  ) %>% 
  select(-contains("_median"))

summary(phy_site_impute)

# Check for correlation
phy_site_impute %>% 
  select(-wateryear, -region, -site) %>% 
  cor(use = "complete.obs")


# Site/Year-level predictors ---------------------------------------------------

# Filter annual predation data to only include analysis regions
pred_year <- pisc_df %>% 
  filter(region != "PHD")


# Summarize annual hydroperiod at site wateryear
hydr_year <- phy_df %>% 
  filter(
    waterperiod == 5,
    wateryear > 1996
    ) %>% 
  group_by(wateryear,region,site) %>% 
  summarize(
    across(wet_sum_365day,~mean(.x,na.rm=T)),
    .groups = "drop"
  )
  
# Combine and summarize at Region/wateryear
phy_reg_year <- hydr_year %>% 
  left_join(pred_year) %>%
  group_by(region,wateryear) %>% 
  summarise(across(
    -site,
    ~mean(.x,na.rm=T)
    ),
    .groups = "drop"
  )

summary(phy_reg_year)

# Check for correlation
phy_reg_year %>% 
  select(-wateryear, -region) %>% 
  cor(use = "complete.obs")


# Export  ----------------------------------------------------------------------
saveRDS(
  phy_site_impute,
  file.path(export_dir,paste0("phys_site_predictors_",Sys.Date(),".rds"))
  )

saveRDS(
  phy_reg_year,
  file.path(export_dir,paste0("phys_regionyear_predictors_",Sys.Date(),".rds"))
)


