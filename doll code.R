# packages

library(readr)
library(dplyr)
library(ape)
library(geosphere)
library(Matrix)
library(brms)
library(posterior)
library(tidyr)
library(ggplot2)
library(maps)
library(png)
library(grid)
library(gridExtra)
library(patchwork)


# data and tree

data_path <- "perfectly_matched_data.csv"
tree_path <- "perfectly_matched_tree.nex"



# presence for 0/1 binary
is_present <- function(x){
  x2 <- tolower(trimws(as.character(x)))
  x2 %in% c("1","true","t","yes","y","present","presence")
}


inv_logit <- function(x) 1 / (1 + exp(-x))


# align data to tree 

d0  <- read_csv(data_path, show_col_types = FALSE)
tr0 <- read.tree(tree_path)

# filter
d0 <- d0 %>% filter(tree_name %in% tr0$tip.label)

# prune phylogeny to match data
tr1 <- drop.tip(tr0, setdiff(tr0$tip.label, d0$tree_name))

# reorder to match tips
d1 <- d0 %>% slice(match(tr1$tip.label, tree_name))
stopifnot(all(d1$tree_name == tr1$tip.label))


# phylogenetic matrix
A <- vcv.phylo(tr1, corr = TRUE)
rownames(A) <- colnames(A) <- d1$tree_name
A_pd <- as.matrix(nearPD(A, corr = TRUE)$mat)   # force PD for brms

# spatial matrix
coords <- cbind(d1$Longitude, d1$Latitude)
D_km <- distm(coords, fun = distHaversine) / 1000

decay_km <- 500
G <- exp(-D_km / decay_km)
rownames(G) <- colnames(G) <- d1$tree_name
G <- (G + t(G)) / 2
diag(G) <- 1
G <- G + diag(1e-4, nrow(G))                    
G_pd <- as.matrix(nearPD(G, corr = TRUE)$mat) 


# filter by variables and standardise complexity
dM <- d1 %>%
  filter(
    !is.na(Play), !is.na(Magic),
    !is.na(Complexity),
    !is.na(Latitude), !is.na(Longitude)
  ) %>%
  mutate(
    Complexity_z = as.numeric((Complexity - mean(Complexity)) / (2*sd(Complexity))),
    phylo = factor(tree_name, levels = tree_name),
    space = factor(tree_name, levels = tree_name)
  )

# naming
idx <- match(dM$tree_name, rownames(A_pd))
A_M <- A_pd[idx, idx]
G_M <- G_pd[idx, idx]


# set priors
priors_main <- c(
  prior(normal(0,1), class="b"),
  prior(normal(0,1), class="Intercept"),
  prior(exponential(1), class="sd", group="phylo"),
  prior(exponential(1), class="sd", group="space")
)

# M0: no phylo/space
fit_magic_M0 <- brm(
  Magic ~ 1 + Complexity_z,
  data = dM,
  family = bernoulli(),
  prior = c(
    prior(normal(0,1), class="b"),
    prior(normal(0,1), class="Intercept")
  ),
  chains=4, cores=4, iter=4000,
  control=list(adapt_delta=0.97),
  seed=123
)

# M1: phylo only
fit_magic_M1 <- brm(
  Magic ~ 1 + Complexity_z + (1 | gr(phylo, cov = A_M)),
  data  = dM,
  data2 = list(A_M = A_M),
  family = bernoulli(),
  prior = c(
    prior(normal(0,1), class="b"),
    prior(normal(0,1), class="Intercept"),
    prior(exponential(1), class="sd", group="phylo")
  ),
  chains=4, cores=4, iter=4000,
  control=list(adapt_delta=0.97),
  seed=123
)

# M2: phylo + space
fit_magic_M2 <- brm(
  Magic ~ 1 + Complexity_z +
    (1 | gr(phylo, cov = A_M)) +
    (1 | gr(space, cov = G_M)),
  data  = dM,
  data2 = list(A_M = A_M, G_M = G_M),
  family = bernoulli(),
  prior = priors_main,
  chains=4, cores=4, iter=4000,
  control=list(adapt_delta=0.99, max_treedepth=12),
  seed=123
)

# play model
fit_play_M2 <- brm(
  Play ~ 1 + Complexity_z +
    (1 | gr(phylo, cov = A_M)) +
    (1 | gr(space, cov = G_M)),
  data  = dM,
  data2 = list(A_M = A_M, G_M = G_M),
  family = bernoulli(),
  prior = priors_main,
  chains=4, cores=4, iter=4000,
  control=list(adapt_delta=0.99, max_treedepth=12),
  seed=123
)

# model comparison for magic (M0 vs M1 vs M2) ---
loo_compare(loo(fit_magic_M0), loo(fit_magic_M1), loo(fit_magic_M2))


# maps

tint_png <- function(img, rgb, alpha_mult = 0.85) {
  if (dim(img)[3] == 3) {
    img <- array(
      c(img, array(1, dim(img)[1:2])),
      dim = c(dim(img)[1:2], 4)
    )
  }
  img[,,4] <- as.numeric(img[,,4]) * alpha_mult
  img[,,1] <- rgb[1]
  img[,,2] <- rgb[2]
  img[,,3] <- rgb[3]
  img
}

world_df <- map_data("world")

play_img  <- readPNG("play_doll.png")
magic_img <- readPNG("magic_doll.png")

play_col  <- c(0.75, 0.15, 0.15)
magic_col <- c(0.10, 0.20, 0.55)

play_icon  <- rasterGrob(tint_png(play_img,  play_col,  alpha_mult = 0.80), interpolate = TRUE)
magic_icon <- rasterGrob(tint_png(magic_img, magic_col, alpha_mult = 0.85), interpolate = TRUE)


dm <- dM %>%
  mutate(
    Play_present  = is_present(Play),
    Magic_present = is_present(Magic)
  )


add_icon_layer <- function(df, present_col, grob, w = 7, h = 7) {
  idx <- which(df[[present_col]])
  if (length(idx) == 0) return(list())
  lapply(idx, function(i) {
    annotation_custom(
      grob = grob,
      xmin = df$Longitude[i] - w/2,
      xmax = df$Longitude[i] + w/2,
      ymin = df$Latitude[i]  - h/2,
      ymax = df$Latitude[i]  + h/2
    )
  })
}

base_map <- ggplot() +
  geom_rect(aes(xmin=-180, xmax=180, ymin=-50, ymax=75),
            fill = "aliceblue", color = NA) +
  geom_polygon(
    data = world_df,
    aes(x = long, y = lat, group = group),
    fill = "grey97",
    color = "grey50",
    linewidth = 0.25
  ) +
  coord_quickmap(xlim = c(-180, 180), ylim = c(-50, 75), expand = FALSE) +
  theme_void() +
  theme(plot.margin = margin(10, 10, 10, 22))

p_play  <- base_map + add_icon_layer(dm, "Play_present",  play_icon,  w=7, h=7)
p_magic <- base_map + add_icon_layer(dm, "Magic_present", magic_icon, w=7, h=7)

p_maps <- (p_play / p_magic) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0.02, 0.98)
  )

p_maps

ggsave("play_magic_dolls_world_map.png", p_maps, width = 10, height = 7, dpi = 300, bg = "white")
ggsave("play_magic_dolls_world_map.pdf", p_maps, width = 10, height = 7, device = cairo_pdf, bg = "white")

# phylo plot

trait_map <- c(
  Toy        = "Play",
  Fertility  = "Fertility",
  Punishment = "Punishment",
  Spirits    = "Spirits",
  Protection = "Protection"
)

trait_labels <- names(trait_map)
trait_cols   <- unname(trait_map)

trait_colors <- c(
  Toy        = "red",
  Fertility  = "forestgreen",
  Punishment = "purple4",
  Spirits    = "darkblue",
  Protection = "orange3"
)

col_absent  <- "white"
col_missing <- "grey80"
col_outline <- "black"

fill_for <- function(status, trait_label){
  if (is.na(status) || status == "Missing") return(col_missing)
  if (status == "Absent") return(col_absent)
  trait_colors[[trait_label]]
}

# match data 
common <- intersect(tr1$tip.label, d1$tree_name)
tr2 <- keep.tip(tr1, common)
d2  <- d1[match(tr2$tip.label, d1$tree_name), ]

status_mat <- sapply(trait_cols, function(col){
  x <- d2[[col]]
  ifelse(is.na(x), "Missing",
         ifelse(is_present(x), "Present", "Absent"))
})
colnames(status_mat) <- trait_labels

draw_tree_plot <- function() {
  
  par(mar = c(3, 3, 9, 12), xpd = NA)
  
  plot.phylo(
    tr2,
    type = "phylogram",
    show.tip.label = TRUE,
    cex = 0.52,
    label.offset = 0,
    no.margin = TRUE,
    edge.width = 1
  )
  
  ntip  <- length(tr2$tip.label)
  y_tip <- 1:ntip
  
  lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  x_tip <- lp$xx[1:ntip]
  
  usr <- par("usr")
  x_range <- usr[2] - usr[1]
  K <- length(trait_labels)
  
  # matrix 
  matrix_block_width <- 0.045 * x_range
  pull_left          <- 0.040 * x_range
  gap_to_matrix      <- 0.003 * x_range
  
  x_matrix_right <- usr[2] - pull_left
  x_matrix_left  <- x_matrix_right - matrix_block_width
  x_cols         <- seq(x_matrix_left, x_matrix_right, length.out = K)
  
  x_line_end <- x_matrix_left - gap_to_matrix
  
  segments(
    x0 = x_tip, y0 = y_tip,
    x1 = x_line_end, y1 = y_tip,
    lty = 3, lwd = 0.8, col = "grey70"
  )
  
  # dots
  dot_cex <- 0.85
  
  for (j in seq_len(K)) {
    trait_label <- trait_labels[j]
    fills <- vapply(
      status_mat[, j],
      fill_for,
      character(1),
      trait_label = trait_label
    )
    
    points(
      x = rep(x_cols[j], ntip),
      y = y_tip,
      pch = 21,
      bg  = fills,
      col = col_outline,
      cex = dot_cex
    )
  }
  
  # legend
  legend(
    x = mean(usr[1:2]),
    y = usr[4] - 0.010 * diff(usr[3:4]),
    legend = trait_labels,
    horiz = TRUE,
    bty = "n",
    xjust = 0.5,
    pch = 21,
    pt.cex = 1.2,
    pt.bg = unname(trait_colors[trait_labels]),
    col = "black",
    x.intersp = 0.9
  )
}

png("tree_play_magic_baseR.png", width = 3600, height = 3200, res = 300)
draw_tree_plot()
dev.off()

pdf("tree_play_magic_baseR.pdf", width = 12, height = 3200/300, useDingbats = FALSE)
draw_tree_plot()
dev.off()


# Extract posteriors for figure
get_beta_draws <- function(fit, term, label){
  dr <- as_draws_df(fit)
  bname <- paste0("b_", term)
  if (!bname %in% names(dr)) stop("Missing parameter: ", bname)
  tibble(model = label, beta = dr[[bname]])
}

col_play_name  <- "red"
col_magic_name <- "darkblue"

draws_pm <- bind_rows(
  get_beta_draws(fit_magic_M2, "Complexity_z", "Magic"),
  get_beta_draws(fit_play_M2,  "Complexity_z", "Play")
) %>%
  mutate(model = factor(model, levels = c("Play","Magic")))

# Panel (a): Odds ratios per +2 SD Complexity (because of (x-mean)/(2*sd))
summ_or <- draws_pm %>%
  mutate(or = exp(beta)) %>%
  group_by(model) %>%
  summarise(
    mean = mean(or),
    lo95 = quantile(or, 0.025),
    hi95 = quantile(or, 0.975),
    lo50 = quantile(or, 0.25),
    hi50 = quantile(or, 0.75),
    .groups = "drop"
  )

p_or <- ggplot(summ_or, aes(y = model, x = mean, color = model)) +
  geom_vline(xintercept = 1, linewidth = 0.45) +
  geom_segment(aes(x = lo95, xend = hi95, yend = model), linewidth = 0.85) +
  geom_segment(aes(x = lo50, xend = hi50, yend = model), linewidth = 2.6) +
  geom_point(size = 3.0) +
  scale_x_log10() +
  scale_color_manual(values = c("Play" = col_play_name, "Magic" = col_magic_name)) +
  theme_classic(base_size = 12) +
  theme(axis.title.y = element_blank(), legend.position = "none") +
  labs(x = "Odds ratio per +2 SD Complexity (log scale)", title = "(a)")

# Panel (b): predicted probabilities across Complexity
grid_x <- tibble(Complexity_z = seq(-2.5, 2.5, length.out = 151))

lin_magic <- posterior_linpred(fit_magic_M2, newdata = grid_x, re_formula = NA, transform = FALSE)
lin_play  <- posterior_linpred(fit_play_M2,  newdata = grid_x, re_formula = NA, transform = FALSE)

p_magic <- inv_logit(lin_magic)
p_play  <- inv_logit(lin_play)

pred_summ <- bind_rows(
  tibble(
    model = "Magic",
    Complexity_z = grid_x$Complexity_z,
    p50 = apply(p_magic, 2, median),
    lo  = apply(p_magic, 2, quantile, probs = 0.025),
    hi  = apply(p_magic, 2, quantile, probs = 0.975)
  ),
  tibble(
    model = "Play",
    Complexity_z = grid_x$Complexity_z,
    p50 = apply(p_play, 2, median),
    lo  = apply(p_play, 2, quantile, probs = 0.025),
    hi  = apply(p_play, 2, quantile, probs = 0.975)
  )
) %>%
  mutate(model = factor(model, levels = c("Play","Magic")))

p_pred <- ggplot(pred_summ, aes(x = Complexity_z, y = p50, color = model, fill = model)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, linewidth = 0) +
  geom_line(linewidth = 1.25) +
  scale_color_manual(values = c("Play" = col_play_name, "Magic" = col_magic_name)) +
  scale_fill_manual(values  = c("Play" = col_play_name, "Magic" = col_magic_name)) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none") +
  labs(x = "Complexity (2 SD units)", y = "Predicted probability", title = "(b)")

final_fig <- arrangeGrob(p_or, p_pred, ncol = 1, heights = c(1, 1.25))
grid.newpage(); grid.draw(final_fig)

ggsave("play_magic_complexity_effects.pdf", final_fig, width = 6.8, height = 7.2, device = cairo_pdf)
ggsave("play_magic_complexity_effects.png", final_fig, width = 6.8, height = 7.2, dpi = 300, bg = "white")
