require(dplyr)
require(tidyverse)
require(GenomicRanges)
library(MetBrewer)
require(ggstatsplot)
library(ggpubr)


plot_theme <- function () {
    theme_bw(base_size = 12, base_family = 'Helvetica') +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(hjust = 0.5),
              text = element_text(size=24, colour = 'black'),
              panel.border = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.ticks.x = element_line(colour = 'black'),
              axis.ticks.y = element_line(colour = 'black'),
              axis.text.x = element_text(colour = 'black'),
              axis.text.y = element_text(colour = 'black'))
}

plot_theme_facet <- function () {
    theme_bw(base_size = 12, base_family = 'Helvetica') +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(hjust = 0.5),
              text = element_text(size = 24, colour = 'black'),
              strip.background = element_blank(),
              panel.border = element_rect(colour="black"),
              axis.ticks.x = element_line(colour = 'black'),
              axis.ticks.y = element_line(colour = 'black'),
              axis.text.x = element_text(colour = 'black'),
              axis.text.y = element_text(colour = 'black'))
}
# change path to file depending on folder structure
all_pre_calculated = read_tsv('./dat/all_withMIN_071222.tsv', show_col_types = FALSE) %>%
    select(-X)

# change path to file depending on folder structure
LP1_pool1_locations = read_tsv('./dat/CAS LP1 bulk locations.tsv', show_col_types = FALSE)

LP1_pool1 =
    all_pre_calculated %>%
    filter(pool == 'LP1_0221') %>%
    inner_join(LP1_pool1_locations, by = 'tBC') %>%
    select(-contains('.'), -counts)

LP1_pool1_2 =
    all_pre_calculated %>%
    filter(pool == 'LP1_0621') %>%
    inner_join(LP1_pool1_locations, by = 'tBC') %>%
    select(-contains('.'), -counts)
# change path to file depending on folder structure
LP3_pool_locations = read_tsv('./dat/TRIP_LP3_mapped_locations_090722.tsv',
                              show_col_types = FALSE)
LP3_pool =
    all_pre_calculated %>%
    filter(pool == 'LP3') %>%
    distinct() %>%
    inner_join(LP3_pool_locations, by = 'tBC') %>%
    select(-contains('.'), -counts)
# change path to file depending on folder structure
LP1_locations = read.csv('./dat/mapped_locations_chromHMM_120321.tsv', sep = '\t')
LP1_latest_pool =
    all_pre_calculated %>%
    filter(pool == 'LP1_latest') %>%
    inner_join(LP1_locations, by = 'tBC') %>%
    select(-contains('.'), -counts)
inner_join(LP1_pool1, LP1_pool1_2, by = c('tBC', 'chr', 'location', 'annotation', 'strand')) %>%
    ggplot(aes(mean_z.x, mean_z.y)) +
    geom_point() +
    plot_theme()
LP1_pool1_unique = anti_join(LP1_pool1, LP1_pool1_2, by = 'tBC')
LP1_pools = bind_rows(LP1_pool1_2, LP1_pool1_unique)
all_pools = bind_rows(LP1_latest_pool, LP1_pools, LP3_pool) %>%
    mutate(id = paste(tBC, pool, sep = '_')) 
all_pools_gr = all_pools %>%
    mutate(start = location, end = location) %>%
    select(-location) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)
all_pools %>%
    select(id, tBC, chr, location, strand, annotation, mean, var, mean_z, var_z, twopower_MIN, pool) %>%
    write_tsv('dat/Supplementary Table 1.txt')
head(all_pools)
pdf("fig2_2.pdf")
all_pools %>%
    mutate(pool = case_when(
        pool == 'LP1_0221' ~ 'Pool 2',
        pool == 'LP1_0621' ~ 'Pool 3',
        pool == 'LP1_latest' ~ 'Pool 4',
        pool == 'LP3' ~ 'Pool 1'
    )) %>%
    ggplot(aes(log10(mean), log10(var))) +
    geom_point(aes(colour = pool)) +
    plot_theme() +
    scale_colour_manual(values = met.brewer('Hokusai1')[c(1, 3, 6, 7)]) +
    xlab('mean (log10)') +
    ylab('variance (log10)') +
    stat_cor(method = "pearson")
all_pools %>%
    mutate(pool = case_when(
        pool == 'LP1_0221' ~ 'Pool 2',
        pool == 'LP1_0621' ~ 'Pool 3',
        pool == 'LP1_latest' ~ 'Pool 4',
        pool == 'LP3' ~ 'Pool 1'
    )) %>%
    ggplot(aes(log10(mean), MIN/2)) +
    geom_point(aes(colour = pool)) +
    plot_theme() +
    scale_colour_manual(values = met.brewer('Hokusai1')[c(1, 3, 6, 7)]) +
    xlab('mean (log10)') +
    ylab('MIN (log10)') +
    stat_cor(method = "pearson")
dev.off()

pdf("SuppFigure_2C_2D.pdf")
all_pools %>%
    mutate(pool = case_when(
        pool == 'LP1_0221' ~ 'Pool 2',
        pool == 'LP1_0621' ~ 'Pool 3',
        pool == 'LP1_latest' ~ 'Pool 4',
        pool == 'LP3' ~ 'Pool 1'
    )) %>%
    ggplot(aes(log10(mean), log10(var/mean))) +
    geom_point(aes(colour = pool)) +
    plot_theme() +
    scale_colour_manual(values = met.brewer('Hokusai1')[c(1, 3, 6, 7)]) +
    xlab('mean (log10)') +
    ylab('fano factor (log10)') +
    stat_cor(method = "pearson")
all_pools %>%
    mutate(pool = case_when(
        pool == 'LP1_0221' ~ 'Pool 2',
        pool == 'LP1_0621' ~ 'Pool 3',
        pool == 'LP1_latest' ~ 'Pool 4',
        pool == 'LP3' ~ 'Pool 1'
    )) %>%
    ggplot(aes(log10(mean), log10(var/mean^2))) +
    geom_point(aes(colour = pool)) +
    plot_theme() +
    scale_colour_manual(values = met.brewer('Hokusai1')[c(1, 3, 6, 7)]) +
    xlab('mean (log10)') +
    ylab('CV2 (log10)') +
    stat_cor(method = "pearson")
dev.off()

