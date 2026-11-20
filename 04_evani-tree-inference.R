library(dplyr); library(ggplot2); library(vroom); library(latex2exp)
library(tidyr)
library(cowplot)

resdir <- "results"
figdir <- file.path(resdir, "figs")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
mc = c( "#9D4030", "#D0D55C", "#C3A97E", "#809D6F", "#769DA6", "#A69E33", "#734002", "#595622", "#8C873F",  "#474B71")

df <- vroom(file.path(resdir, "tree-eval.tsv"), na = "")
dft <- df %>%
  filter(builder=="NJ*" & study=="mutation") %>%
  mutate(nrate=cut(rate, c(5, 50, 300, 500), labels = c("5-50", "100-300", "400-500"), include.lowest = T)) %>%
  complete(nrate, method)
lvl=c("gdiff", "fastani", "mash", "dashing2", "skani")

p <- dft %>% mutate(method=factor(method, level=lvl)) %>%
  ggplot() +
  aes(y=method, x=nRF, fill=method) +
  facet_wrap(~nrate) +
  stat_summary(geom="bar", color="gray10", show.legend = F) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
  geom_text(
    data = dft %>%  complete(nrate, method) %>%  group_by(nrate, method) %>%  filter(all(is.na(nRF))),
    aes(y=method, x = 0.025, label = "x"),
    color = "black", 
    size = 6, 
    inherit.aes = FALSE,
  ) + 
  geom_text(
    # missing pairs as a percentage of the 15*14/2 = 105 pairs in the matrix
    data = dft %>%  complete(nrate, method) %>%  group_by(nrate, method) %>%  summarise(nk=100-sum(is.na(nRF))/n()*100, mp=mean(100*n_missing/(15*14/2))),
    aes(y=method, x = 0.45, label = paste(round(mp, 0), "%", sep = "")),
    # aes(x=method, y = 0.45, label = paste(round(nk, 2), "%"sep = "")),
    color = "black", 
    fontface="bold",
    size = 3.5, 
    inherit.aes = FALSE,
  ) + 
  theme_bw() +
  coord_cartesian(xlim=c(0, 0.5)) +
  scale_fill_manual(values = mc) +
  labs(y="", x="Normalized RF Error", title="EvANI benchmarking: species tree estimation using NJ*") +
  # theme(axis.text.x = element_text(angle=45, hjust=1)) +
  theme(title = element_text(size = 10)) +
  scale_y_discrete(name="") + theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )
p
ggsave("./E-evani_nrf_error.pdf", width=5, height = 2)
