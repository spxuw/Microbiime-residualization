suppressPackageStartupMessages({
    library(readr); library(dplyr); library(stringr); library(purrr); 
    library(tidyr); library(ggplot2); library(patchwork); library(ggrepel)
})

auroc_csv <- "D:/桌面/Rstudio/大队列研究/age/auroc.csv"
out_dir   <- "D:/桌面/Rstudio/大队列研究/age/out_figs"

# ==========================================
# 老师指定的作图模板 (全局加黑加深，设定基准字号)
# ==========================================
custom_theme <- theme_bw(base_size = 12) +
    theme(
        text = element_text(color = "black"), # 强制全局基础文本为纯黑
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(linewidth = 0.2, color = "black"),
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.3, color = "black"),
        axis.line.x.top = element_blank(),
        axis.line.y.right = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 9, color = "black"),
        strip.background = element_blank(),
        axis.text.x = element_text(size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"),
        plot.title = element_text(color = "transparent", face = "bold", size = 13, hjust = 0.5),
        legend.position = "none" # 默认关闭，稍后在画图函数中强行开启
    )

# 1. 强制读取 Figure 1B 的基准清单
proj_list <- read_csv(file.path(out_dir, "Figure1_projects_final.csv"), show_col_types = FALSE)

# 2. 读取并清洗 auroc 数据
au <- read_csv(auroc_csv, show_col_types=FALSE)
col_ds   <- grep("dataset|file|data", names(au), value=TRUE, ignore.case=TRUE)[1]
col_bp   <- grep("BioProject|project|study|ID", names(au), value=TRUE, ignore.case=TRUE)[1]
col_null <- grep("auc_null.*mean|null_auc.*mean|pheno_auc_null_mean", names(au), value=TRUE, ignore.case=TRUE)[1]
col_rf   <- grep("full.*rf.*auc.*mean|pheno_auc_full_rf_mean", names(au), value=TRUE, ignore.case=TRUE)[1]
col_xgb  <- grep("full.*xgb.*auc.*mean|pheno_auc_full_xgb_mean", names(au), value=TRUE, ignore.case=TRUE)[1]

df_raw <- au |>
    mutate(
        study = as.character(.data[[col_bp]]),
        dataset_base = str_replace(as.character(.data[[col_ds]]), "\\.csv$", ""),
        source = case_when(
            dataset_base == "CRC" ~ "CRC", 
            str_detect(dataset_base, "^GMrepo_") ~ "GM", 
            str_detect(dataset_base, "^16S_") ~ "Res", 
            TRUE ~ NA_character_
        ),
        dataset_key = case_when(source == "GM" ~ "GMrepo", TRUE ~ dataset_base),
        study_label = if_else(source == "Res", paste0(study, " [", str_remove(dataset_key, "^16S_"), "]"), study)
    ) |> 
    inner_join(proj_list |> select(source, study_label, order), by=c("source","study_label")) |>
    # ============
# 剔除 auroc.csv 中同一 source 和 study_label 的重复行，确保每个队列只有一条数据！
distinct(source, study_label, .keep_all = TRUE) |> 
    # ==============================
mutate(
    auc_null = as.numeric(.data[[col_null]]),
    delta_rf  = as.numeric(.data[[col_rf]]) - auc_null,
    delta_xgb = as.numeric(.data[[col_xgb]]) - auc_null,
    color = case_when(source == "CRC" ~ "#C97DBB", source == "GM"  ~ "#E69F00", source == "Res" ~ "#5CB8E8")
) |> 
    group_by(source) |> 
    arrange(desc(order)) |> 
    mutate(
        point_id = row_number(),
        study_label = factor(study_label, levels = unique(study_label))
    ) |> 
    ungroup()

# 3. 相关性文字格式函数
format_p <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 0.001) return("< 0.001")
    sprintf("= %.3f", p)
}

# 4. 绘图核心函数
plot_scatter_final <- function(src, delta_col, title_text) {
    df_sub <- df_raw |> 
        filter(source == src) |> 
        rename(y_val = !!sym(delta_col)) |> 
        arrange(point_id) |> 
        mutate(study_label = factor(study_label, levels = unique(study_label)))
    
    threshold <- 0.010
    dist_mat <- as.matrix(dist(df_sub[, c("auc_null", "y_val")]))
    diag(dist_mat) <- Inf
    df_sub$is_overlap <- apply(dist_mat, 1, function(x) any(x < threshold, na.rm=TRUE))
    
    # Study-name legend/annotation removed; numbered point labels are retained.
    
    # Pearson correlation between age-only AUROC and microbiome-derived ΔAUROC
    ct <- cor.test(df_sub$auc_null, df_sub$y_val, method = "pearson")
    corr_label <- paste0("Pearson r = ", sprintf("%.2f", unname(ct$estimate)),
                         "\nP ", format_p(ct$p.value))
    
    ggplot(df_sub, aes(x = auc_null, y = y_val)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.5) +
        geom_smooth(method = "lm", se = FALSE, color = "grey35", linewidth = 0.6) +
        geom_point(aes(fill = study_label), color = "grey20", shape = 21, size = 5.2, stroke = 0.4, alpha = 0.9) +
        geom_text(data = filter(df_sub, !is_overlap), aes(label = point_id), 
                  size = 3.55, fontface = "bold", color = "black") +
        geom_text_repel(data = filter(df_sub, is_overlap), aes(label = point_id),
                        size = 3.55, fontface = "bold", color = "black",
                        point.padding = 0, box.padding = 0.5, 
                        min.segment.length = 0, segment.color = "grey30", segment.size = 0.4,
                        force = 4, seed = 42) +
        annotate("text",
                 x = Inf, y = Inf,
                 label = corr_label,
                 hjust = 1.05, vjust = 1.15,
                 size = 3.15, color = "black") +
        
        scale_fill_manual(values = setNames(rep(df_sub$color[1], nrow(df_sub)), df_sub$study_label), guide = "none") +
        
        # 【不挡点精髓】：X轴强制扩展 55% 留给右侧放图注，绝对不会挡住数据
        scale_x_continuous(expand = expansion(mult = c(0.1, 0.55))) +
        scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
        
        labs(title = title_text, x = "Age-only AUROC (Null)", y = expression(Delta*AUROC)) +
        custom_theme + # 应用老师模板
        theme(
            # Keep title layout space, but make the title text invisible.
            plot.title = element_text(face = "bold", size = 13, hjust = 0.5, color = "transparent"),
            plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "mm"),
            legend.position = "none"
        )
}

# 5. 保存结果 (应用导师指定的 scale 和宽高)
# RF
p_rf <- (plot_scatter_final("CRC", "delta_rf", "CRC") | 
             plot_scatter_final("GM",  "delta_rf", "GMrepo") | 
             plot_scatter_final("Res", "delta_rf", "ResMicroDB"))
ggsave(file.path(out_dir, "Figure3_Scatter_RF.pdf"), p_rf, width=12, height=4.5, scale=0.8)
ggsave(file.path(out_dir, "Figure3_Scatter_RF.png"), p_rf, width=12, height=4.5, dpi=600)

# XGBoost
p_xgb <- (plot_scatter_final("CRC", "delta_xgb", "CRC") | 
              plot_scatter_final("GM",  "delta_xgb", "GMrepo") | 
              plot_scatter_final("Res", "delta_xgb", "ResMicroDB"))
ggsave(file.path(out_dir, "Figure3_Scatter_XGB.pdf"), p_xgb, width=12, height=4.5, scale=0.8)
ggsave(file.path(out_dir, "Figure3_Scatter_XGB.png"), p_xgb, width=12, height=4.5, dpi=600)