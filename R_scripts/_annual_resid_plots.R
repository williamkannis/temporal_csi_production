
mod_dir <- "stan_outputs"
data_dir <- "hpc/data"

mod_files <- list.files(mod_dir)
data_files <- list.files(data_dir)

out <- readRDS(file.path(mod_dir,mod_files[4]))
data <- readRDS(file.path(data_dir,data_files[5]))

a <- out$summary("raw_residuals")

data$rg[data$st]
data$st

resid_df <- data.frame(
  year = data$yr,
  region = factor(data$rg[data$st]),
  resid = a$mean
) %>% 
  group_by(year,region) %>% 
  summarise(resid = median(resid))

library(ggplot2)
ggplot(resid_df,aes(x=year,y=resid,colour = region)) +
  geom_point()

abline(h=0)

resid_wide <- resid_df %>% 
  tidyr::pivot_wider(
    names_from = region,
    values_from = resid
    )

library(pheatmap)

pheatmap(resid_wide[,-1], 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         # color = hcl.colors(50, "Viridis"),
         border_color = "white")


corrplot::corrplot(as.matrix(resid_wide[,-1]),is.corr = F,method="square")
image(as.matrix(resid_wide[,-1]))




out_pb <- readRDS(file.path(mod_dir,mod_files[5]))
data_pb <- readRDS(file.path(data_dir,data_files[6]))

diag <- hurdle_diagnostics(
    fit = out_pb,
    y = data_pb$y,
    yr = data_pb$yr,
    st = data_pb$st,
    rg = data_pb$rg,
    data = NULL,
    ndraws = 1000,
    seed = 123,
    ppc_groups = NULL
)





