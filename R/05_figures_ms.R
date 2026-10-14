## 05_figures_ms.R: manuscript Figs 2-4

library(tidyr)
library(ggplot2)
library(ggdist)
library(patchwork)

source(here::here("R", "00_functions.R"))
theme_set(theme_classic(base_size = 10))

va  <- readRDS(out("va.rds"))
d   <- va$d
R_w <- cor_from_G(d$G_w)
rA  <- apply(R_w, 2:3, mean)

save_fig <- function(p, name, w, h) {
  ggsave(out(paste0(name, ".pdf")), p, width = w, height = h, units = "mm")
  ggsave(out(paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 600,
         device = ragg::agg_png)
}


## ---- Fig. 2 -----------------------------------------------------------------
## Three panels in one row: mean fitness, V_A(w), and its channel split.
## Posterior mean and 90% HPD of each column of a D x 6 matrix, as in the tables.
trt_summary <- function(M) bind_rows(lapply(TREATS, function(tt)
  cbind(temp = TEMP_OF[[tt]], diet = DIET_OF[[tt]], smry(M[, tt]))))
TRT_LAB <- c(T14DA = "14\nstd", T14DB = "14\nlow", T25DA = "25\nstd",
             T25DB = "25\nlow",  T28DA = "28\nstd", T28DB = "28\nlow")

## Diet by line type; the legend sits in (a) only.
p2a <- ggplot(trt_summary(log(va$Wbar)), aes(factor(temp), mean, group = diet)) +
  geom_line(aes(linetype = diet)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.2, stroke = 0.25) +
  scale_linetype_manual(values = c(DA = "solid", DB = "dashed"), labels = DIET_LAB) +
  ## extra room below for the legend
  scale_y_continuous(breaks = -3:2, expand = expansion(mult = c(0.18, 0.05))) +
  labs(title = bquote(bold("(a)") ~ "Mean fitness"), x = "Temperature (°C)",
       y = "log mean fitness\n(offspring per day of laying)", linetype = NULL) +
  theme(axis.title.y = element_text(size = 8.5, lineheight = 0.9),
        legend.position = c(0.68, 0.1), legend.background = element_blank(),
        legend.text = element_text(size = 8), legend.key.height = unit(3.5, "mm"))

pd2 <- position_dodge(width = 0.15)   # separate the overlapping intervals
p2b <- ggplot(trt_summary(d$VA_w), aes(factor(temp), mean, group = diet)) +
  geom_line(aes(linetype = diet), position = pd2) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.2, stroke = 0.25, position = pd2) +
  scale_linetype_manual(values = c(DA = "solid", DB = "dashed"), guide = "none") +
  labs(title = bquote(bold("(b)") ~ "Adaptive potential"), x = "Temperature (°C)",
       y = expression(V[A](w)))

p2c <- tibble(treat = rep(TREATS, 2), v = c(colMeans(d$count), colMeans(d$zi)),
              channel = factor(rep(c("count", "zero-inflation"), each = length(TREATS)),
                               c("zero-inflation", "count"))) %>%
  ggplot(aes(factor(TRT_LAB[treat], TRT_LAB), v, fill = channel)) +
  geom_col() +
  scale_fill_manual(values = c(count = "steelblue", `zero-inflation` = "coral")) +
  labs(title = bquote(bold("(c)") ~ "Variance decomposition"), x = NULL,
       y = expression(V[A](w)), fill = NULL) +
  theme(legend.position = c(0.72, 0.85), legend.text = element_text(size = 8),
        legend.key.size = unit(3.5, "mm"))

save_fig((p2a + p2b + p2c + plot_layout(nrow = 1)) &
           theme(plot.title = element_text(size = 10.5), plot.title.position = "plot"),
         "ms_fig2", 180, 70)


## ---- Fig. 3 -----------------------------------------------------------------
## Standardised sire breeding values (sire_bv, saved by 02) in two treatments,
## sires observed in both; y label on the first panel only, r_A as subtitle.
bv <- readRDS(FIT_PATH)$sire_bv
env_lab <- function(tt) paste0(TEMP_OF[[tt]], " °C\n", DIET_LAB[[DIET_OF[[tt]]]])

panel <- function(left, right, tag, title,
                  ylab = "Sire breeding value for relative fitness (standardised)") {
  bvs  <- function(tt) with(bv[bv$TREAT == tt & bv$observed, ], setNames(bv_mean, SIRE))
  l <- bvs(left); r <- bvs(right); keep <- intersect(names(l), names(r))
  dd <- tibble(SIRE   = rep(keep, each = 2),
               env    = rep(1:2, length(keep)),
               BV     = as.vector(rbind(as.numeric(scale(l[keep])), as.numeric(scale(r[keep])))),
               leftBV = rep(unname(l[keep]), each = 2))
  ## Jitter once in the data so lines and points share the same x.
  set.seed(123)
  dd$x <- dd$env + runif(nrow(dd), -0.075, 0.075)
  p <- ggplot(dd, aes(x, BV, group = SIRE, colour = leftBV)) +
    geom_line(alpha = 0.6, colour = "grey40") +
    geom_point(size = 1.3) +
    scale_x_continuous(breaks = 1:2, labels = c(env_lab(left), env_lab(right)),
                       expand = expansion(add = 0.4)) +
    scale_colour_gradientn(colours = c("#440154", "#21908c", "#fde725")) +
    theme(legend.position = "none") +
    labs(title = bquote(bold(.(tag)) ~ .(title)),
         subtitle = bquote(r[A] == .(sprintf("%.3f", rA[left, right]))),
         x = NULL, y = ylab)
  if (is.null(ylab)) p <- p + theme(axis.text.y = element_blank())
  p
}
p3s <- list(panel("T25DA", "T28DA", "(a)", "Temperature only"),
            panel("T25DA", "T28DB", "(b)", "Temperature × diet", ylab = NULL),
            panel("T14DB", "T28DB", "(c)", "Condition extremes", ylab = NULL))
## Shared y axis across the three panels.
ylim3 <- range(unlist(lapply(p3s, function(p) p$data$BV)))
save_fig(wrap_plots(p3s, nrow = 1) & scale_y_continuous(limits = ylim3), "ms_fig3", 180, 120)


## ---- Fig. 4 -----------------------------------------------------------------
## (a) diet effects and interactions on V_A(w), (b) margin contrasts on r_A.
DENS_COL <- "#4C72B0"

## Label each density with P(>0), just right of the density at fraction at of its
## height (slab geometry read back from the built plot). s: the summary table, with the
## y factor column y and, for a faceted plot, the facet factor column facet;
## pad: gap to the label; min_y: lowest label height above the row.
add_p0 <- function(p, s, y, facet = NULL, pad, min_y = 0, at = 0.5) {
  pos <- layer_data(p, 2) %>% group_by(PANEL, y) %>%
    summarise(peak = max(thickness), sc = scale[1],
              x = max(x[thickness >= peak * at]) + pad, .groups = "drop") %>%
    group_by(PANEL) %>% mutate(h = sc * peak / max(peak)) %>% ungroup() %>%   # per-panel scaling
    mutate(PANEL = as.integer(PANEL), y = as.integer(y))
  s <- mutate(s, PANEL = if (is.null(facet)) 1L else as.integer(.data[[facet]]),
              y = as.integer(.data[[y]]))
  p + geom_text(data = left_join(s, pos, by = c("PANEL", "y")),
                aes(x, y + pmax(h * at, min_y), label = sprintf("P(>0) = %.3f", p_gt0)),
                inherit.aes = FALSE, hjust = 0, size = 5.5 / .pt)
}

vc <- va_contrasts(d$VA_w) %>%
  mutate(quantity = sub("^(diet effect|interaction \\w+) ", "\\1\n", quantity),
         quantity = gsub("(\\d+)C", "\\1 °C", quantity),
         quantity = factor(quantity, rev(quantity)))
vc_sum <- summarise_draws_col(vc)

p4a <- ggplot(unnest(vc, draws), aes(draws, quantity)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  stat_slab(fill = DENS_COL, slab_alpha = 0.75) +
  geom_pointinterval(data = vc_sum, aes(x = mean, xmin = lower, xmax = upper)) +
  facet_wrap(~ " ") +   # placeholder strip title, matching (b)
  coord_cartesian(xlim = c(-1, 2.7)) +   # crops view only; thin upper tails beyond
  labs(x = expression("Difference in" ~ V[A](w)), y = NULL) +
  theme(strip.background = element_blank(), strip.text = element_text(hjust = 0.5))
p4a <- add_p0(p4a, vc_sum, "quantity", pad = 0.06, min_y = 0.36)   # clear the low, wide densities

lab_M <- c(Cold = "Cold margin\n(14 vs 25 °C)", Hot = "Hot margin\n(25 vs 28 °C)")
lab_P <- c(P1 = "P1: correlation\nacross\ntemperatures", P2 = "P2: correlation\nacross diets")
mc <- margin_contrasts(R_w) %>%
  mutate(margin = factor(lab_M[margin], lab_M), contrast = factor(lab_P[contrast], rev(lab_P)))
mc_sum <- summarise_draws_col(mc)

p4b <- ggplot(unnest(mc, draws), aes(draws, contrast)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  stat_slab(fill = DENS_COL, slab_alpha = 0.75, normalize = "panels", scale = 0.6) +
  geom_pointinterval(data = mc_sum, aes(x = mean, xmin = lower, xmax = upper)) +
  facet_wrap(~ margin, nrow = 1) +
  coord_cartesian(xlim = range(unlist(mc$draws)) + c(0, 0.35), expand = FALSE) +   # tight, plus room for labels
  scale_x_continuous(breaks = seq(-0.5, 1, 0.5)) +
  labs(x = expression("Difference in" ~ r[A]), y = NULL) +
  theme(strip.background = element_blank(), strip.text = element_text(hjust = 0.5),
        axis.text.y = element_text(hjust = 0.5))
p4b <- add_p0(p4b, mc_sum, "contrast", facet = "margin", pad = 0.04, at = 0.75)

save_fig(p4a + p4b + plot_layout(widths = c(1, 1)) +
           plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
           theme(plot.tag = element_text(face = "bold"), plot.tag.location = "plot"),
         "ms_fig4", 180, 85)

