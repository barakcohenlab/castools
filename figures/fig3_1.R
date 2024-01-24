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
