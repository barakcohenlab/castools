library(rtracklayer)
library(ggstatsplot)
library(circlize)
library(EnrichedHeatmap)
library(ggpubr)
library(MetBrewer)
library(GenomicRanges)
library(tidyverse)
library(ggplotHelper)

# Plot histone mods (Figure 3A/3B)
# change path to file depending on folder structure
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

all_pre_calculated = read_tsv('dat/all_withMIN_071222.tsv', show_col_types = FALSE) %>%
    select(-X)
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
    write_tsv('Supplementary Table 1.txt')
all_pools = bind_rows(LP1_latest_pool, LP1_pools, LP3_pool) %>%
    mutate(id = paste(tBC, pool, sep = '_')) 
head(all_pools)

all_pools_fitted_rank = all_pools %>%
        mutate(mean_rank = as.character(ntile(mean_z, 2)),
           var_rank = as.character(ntile(var, 2)),
           MIN_rank = as.character(ntile(twopower_MIN, 2)),
           cv2_rank = as.character(ntile(cv2_z, 2)),
           fano_rank = as.character(ntile(fano, 2))) %>%
    mutate(mean_rank = case_when(
        mean_rank == '1' ~ 'low mean',
        mean_rank == '2' ~ 'high mean'
    ), MIN_rank = case_when(
        MIN_rank == '1' ~ 'low MIN',
        MIN_rank == '2' ~ 'high MIN'
    ))

all_pools_fitted_rank_gr = all_pools_fitted_rank %>%
    mutate(start = location, end = location) %>%
    select(-location) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)

set.seed(498)
random_assign = sample(1:2, size = 939, replace = TRUE)

# change path to files depending on folder structure
histone_mods = fs::dir_ls(path = 'dat/Histone_modifications/', glob = '*_no_dup.bed.gz') %>%
  as.character() %>%
  purrr::set_names(function(x) {
    str_sub(fs::path_file(x), end = - str_length('_hg38_narrowpeaks_no_dup.bed.gz') - 1)
  })

narrowpeak_cols = c(signal_value = "numeric", 
                    p_value = "numeric",
                    q_value = "numeric", 
                    peak = "integer")

import_narrowpeak_file = function(file_path){
    import(file_path, format = 'BED', extraCols = narrowpeak_cols)
}

histone_mods_gr = histone_mods %>%
    map(~ import_narrowpeak_file(.x))

norm_to_mat = function(feature, gr){
    
    normalizeToMatrix(feature, 
                        gr, 
                        value_column = "score", 
                        extend = 5000, 
                        mean_mode = "w0", 
                        smooth = TRUE,
                        background = 0,
                        w = 50
    )
}

histone_mods_rank_normmat = histone_mods_gr %>%
    map(~ norm_to_mat(.x, all_pools_fitted_rank_gr))

plot_enriched_hm_comparison = function(current_mat, title, colour){
  
    col_fun = colorRamp2(quantile(current_mat, c(0, 0.99)), c("white", colour))
    
    EnrichedHeatmap(current_mat, 
                    column_title = title, 
                    col = col_fun,
                    name = title, 
                    top_annotation = HeatmapAnnotation(
                        lines = anno_enriched(
                          gp = gpar(
                            col = c('#348AA7', '#7A4419')
                            # ,lty = 1:2
                            ),
                            height = unit(2, 'in'),
                            axis_param = list(
                                side = 'left'
                            ))),
                    use_raster = TRUE)
}

selected_histones = c('H3K27ac', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K79me2', 'H3K9ac')

# change to MIN for MIN plots
lgd = Legend(at = c("low mean", "high mean"), title = "", 
    type = "lines", legend_gp = gpar(col = c('#348AA7', '#ECC30B')))

# , lty = 1:2
ht_list = NULL

    for (i in 1:length(histone_mods_rank_normmat)){
        current_name = names(histone_mods_rank_normmat)[i]
        if (current_name %in% selected_histones){
            ht_list = ht_list + 
            plot_enriched_hm_comparison(histone_mods_rank_normmat[[i]], current_name, 'blue')
        }
    }

# mean


add_anno_enriched = function(ht_list, name, ri, ci) {
    pushViewport(viewport(layout.pos.row = ri, layout.pos.col = ci))
    extract_anno_enriched(ht_list, name, newpage = FALSE)
    upViewport()
}


ht_list_1 =
    draw(ht_list, split = all_pools_fitted_rank$mean_rank, ht_gap = unit(0.7, "cm"), annotation_legend_list = list(lgd))
pdf('./fig3_A.pdf', width = 14, height = 3)
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 6)))
for (i in seq_along(selected_histones)){
    add_anno_enriched(ht_list_1, selected_histones[i], 1, i)
}
upViewport()
dev.off()

# MIN
lgd = Legend(at = c("low MIN", "high MIN"), title = "", 
    type = "lines", legend_gp = gpar(col = c('#348AA7', '#ECC30B')))

add_anno_enriched = function(ht_list, name, ri, ci) {
    pushViewport(viewport(layout.pos.row = ri, layout.pos.col = ci))
    extract_anno_enriched(ht_list, name, newpage = FALSE)
    upViewport()
}


ht_list_2 =
    draw(ht_list, split = all_pools_fitted_rank$MIN_rank, ht_gap = unit(0.7, "cm"), annotation_legend_list = list(lgd))
pdf('./fig3_B.pdf', width = 14, height = 3)

pushViewport(viewport(layout = grid.layout(nr = 1, nc = 6)))
for (i in seq_along(selected_histones)){
    add_anno_enriched(ht_list_2, selected_histones[i], 1, i)
}
upViewport()
dev.off()


calculate_average_histones = function(mod, df){

    as_tibble(histone_mods_rank_normmat[[mod]]) %>%
        mutate(average = rowMeans(., na.rm = TRUE)) %>%
        select(average) %>%
        bind_cols(df) %>%
        mutate(histone = mod, average = average + 1)
}

c('H3K27ac', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K79me2', 'H3K9ac') %>%
    map(~calculate_average_histones(.x, all_pools_fitted_rank)) %>%
    reduce(bind_rows) %>%
    ggplot(aes(mean_rank, log2(average))) +
    geom_violin(fill = 'gray') +
    # geom_boxplot(width = 0.1) +
    geom_jitter(width = 0.2, alpha = 0.2) +
    plot_theme_facet() +
    facet_wrap(~ histone, nrow = 2) +
    xlab('') +
    ylab('log2(histone signals)') +
    geom_signif(comparisons = list(c('low mean', 'high mean'))) +
    ylim(-1, 11)

ggsave('SuppFigure3A.pdf', width = 10, height = 6)

c('H3K27ac', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K79me2', 'H3K9ac') %>%
    map(~calculate_average_histones(.x, all_pools_fitted_rank)) %>%
    reduce(bind_rows) %>%
    ggplot(aes(MIN_rank, log2(average))) +
    geom_violin(fill = 'gray') +
    # geom_boxplot(width = 0.1) +
    geom_jitter(width = 0.2, alpha = 0.2) +
    plot_theme_facet() +
    facet_wrap(~ histone, nrow = 2) +
    xlab('') +
    ylab('log2(histone signals)') +
    geom_signif(comparisons = list(c('low MIN', 'high MIN'))) +
    ylim(-1, 11)

ggsave('SuppFigure3B.pdf', width = 10, height = 6)

# Permutate labels (Supplementary Figure 3C)

calculate_average_histones = function(mod, df){

    as_tibble(histone_mods_rank_normmat[[mod]]) %>%
        mutate(average = rowMeans(., na.rm = TRUE)) %>%
        select(average) %>%
        bind_cols(df) %>%
        mutate(histone = mod, average = average + 1)
}

set.seed(4309)
random_seeds = sample(1:5000, size = 1000)

calculate_mean_random_samples = function(current_seed){

    set.seed(current_seed)
    random_grouping = rep(1:2, each = 470)
    random_grouping_truncated = sample(random_grouping[1:939])

    df = all_pools_fitted_rank %>%
        mutate(random_rank = random_grouping_truncated) %>%
        mutate(random_rank = case_when(
            random_rank == 1 ~ 'group1',
            random_rank == 2 ~ 'group2'
        ))

    c('H3K27ac', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K79me2', 'H3K9ac') %>%
        map(~ calculate_average_histones(.x, df)) %>%
        reduce(bind_rows) %>%
        group_by(random_rank, histone) %>%
        summarise(mean_histones = mean(average), .groups = 'drop') %>%
        pivot_wider(id_cols = histone, names_from = random_rank, values_from = mean_histones) %>%
        mutate(diff = abs(group2-group1)) %>%
        mutate(id = current_seed)
}

all_random = random_seeds %>%
    map(~ calculate_mean_random_samples(.x)) %>%
    reduce(bind_rows)

mean_diffs = c('H3K27ac', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K79me2', 'H3K9ac') %>%
    map(~ calculate_average_histones(.x,  all_pools_fitted_rank)) %>%
    reduce(bind_rows) %>%
    group_by(mean_rank, histone) %>%
    summarise(mean_histones = mean(average)) %>%
    pivot_wider(id_cols = histone, names_from = mean_rank, values_from = mean_histones) %>%
    mutate(diff = abs(`low mean` - `high mean`))

MIN_diffs = c('H3K27ac', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K79me2', 'H3K9ac') %>%
    map(~ calculate_average_histones(.x,  all_pools_fitted_rank)) %>%
    reduce(bind_rows) %>%
    group_by(MIN_rank, histone) %>%
    summarise(mean_histones = mean(average)) %>%
    pivot_wider(id_cols = histone, names_from = MIN_rank, values_from = mean_histones) %>%
    mutate(diff = abs(`low MIN` - `high MIN`))

calculate_label_perm_p = function(current_mod, mean_MIN){

    if (mean_MIN == 'mean'){
        actual = filter(mean_diffs, histone == current_mod)$diff
    } else if (mean_MIN == 'MIN'){
        actual = filter(MIN_diffs, histone == current_mod)$diff
    }

    n_greater = all_random %>%
        filter(histone == current_mod & diff > actual) %>%
        nrow()

    return(n_greater/length(random_seeds))
}

selected_histones = c('H3K27ac', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K79me2', 'H3K9ac')

mean_perm_pvals = selected_histones %>%
    map(~ calculate_label_perm_p(.x, 'mean')) %>%
    set_names(selected_histones) %>%
    as.tibble() %>%
    pivot_longer(cols = colnames(.), names_to = 'histone', values_to = 'pvalue') %>%
    mutate(type = 'mean')

MIN_perm_pvals = selected_histones %>%
    map(~ calculate_label_perm_p(.x, 'MIN')) %>%
    set_names(selected_histones) %>%
    as.tibble() %>%
    pivot_longer(cols = colnames(.), names_to = 'histone', values_to = 'pvalue') %>%
    mutate(type = 'MIN')

bind_rows(mean_perm_pvals, MIN_perm_pvals) %>%
    ggplot(aes(histone, pvalue)) +
    geom_col(fill = 'lightgray', col = 'black') +
    facet_wrap(~ type) +
    plot_theme_facet() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(expand = c(0,0), limits = c(-0.01, 0.6))

ggsave('SuppFigure3C.pdf', width = 6, height = 5)
