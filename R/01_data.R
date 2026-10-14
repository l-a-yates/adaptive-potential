## 01_data.R: load and check the data, save results/data/fitdat.rds,
## and write the design tables and descriptive plots.

library(ggplot2)

source(here::here("R", "00_functions.R"))
theme_set(theme_classic())

fitdat <- read.csv(here::here("data", "Data_Productivity.csv"))

fitdat$SET   <- as.factor(paste(fitdat$TREAT, fitdat$SET, sep = ":"))
fitdat$SIRE  <- as.factor(fitdat$SIRE)
fitdat$DAM   <- as.factor(paste0(fitdat$SIRE, fitdat$DAM))
fitdat$TREAT <- factor(fitdat$TREAT, levels = TREATS)

stopifnot(nrow(fitdat) == 6906)
saveRDS(fitdat, out("fitdat.rds"))


## ---- design -----------------------------------------------------------------
design <- fitdat %>%
  group_by(treat = TREAT) %>%
  summarise(n         = n(),
            n_sire    = n_distinct(SIRE),
            n_dam     = n_distinct(DAM),
            n_set     = n_distinct(SET),
            pct_zero  = round(100 * mean(TOTAL == 0), 1),
            mean_prod = round(mean(TOTAL), 1),
            max_prod  = max(TOTAL),
            .groups   = "drop") %>%
  mutate(temp = TEMP_OF[as.character(treat)],
         diet = DIET_LAB[DIET_OF[as.character(treat)]],
         lay_days = LAYDAYS[as.character(treat)],
         prod_per_day = round(mean_prod / lay_days, 2)) %>%
  relocate(treat, temp, diet)

print(as.data.frame(design), row.names = FALSE)
write.csv(design, out("tab_design.csv"), row.names = FALSE)

cat("\nN =", nrow(fitdat), "females,", nlevels(fitdat$SIRE), "sires,",
    nlevels(fitdat$DAM), "full-sib families,",
    nlevels(fitdat$SET), "treatment-by-set cells\n")
cat("overall zeros:", round(100 * mean(fitdat$TOTAL == 0), 1), "%   range:",
    min(fitdat$TOTAL), "-", max(fitdat$TOTAL), "\n")


## ---- sire x dam replication -------------------------------------------------
replication <- bind_rows(lapply(TREATS, function(tt) {
  x <- droplevels(fitdat[fitdat$TREAT == tt, ])
  tab <- as.data.frame.matrix(table(x$SIRE, sub("^[0-9]+", "", as.character(x$DAM))))
  tab$SIRE <- rownames(tab); tab$Treatment <- tt
  tab
}))
write.csv(replication, out("tab_replication.csv"), row.names = FALSE)

cat("\n=== females per non-empty treatment x sire x dam cell ===\n")
cnt <- with(fitdat, table(TREAT, SIRE, DAM))
print(table(cnt[cnt > 0]))


## ---- distribution of productivity -------------------------------------------
p_hist <- ggplot(fitdat, aes(x = TOTAL)) +
  geom_histogram(bins = 40, fill = "grey30") +
  facet_wrap(~ TREAT, ncol = 3, scales = "free_y") +
  labs(x = "Productivity (offspring per female)", y = "Count")
ggsave(out("fig_productivity_hist.pdf"), p_hist, height = 160, width = 180, units = "mm")


## ---- larval density ---------------------------------------------------------
p_dens <- ggplot(fitdat, aes(x = Offspring_Development, y = TOTAL)) +
  geom_point(alpha = 0.25, size = 0.6) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  facet_wrap(~ TREAT, ncol = 3, scales = "free") +
  labs(x = "Offspring emerging from development vial", y = "Productivity")
ggsave(out("fig_density_check.pdf"), p_dens, height = 200, width = 240, units = "mm")
