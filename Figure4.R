suppressPackageStartupMessages({
    library(readr); library(dplyr); library(tidyr); library(tibble);
    library(randomForest); library(stringr); library(pROC); library(ggplot2); 
    library(ggrepel); library(ggpubr); library(patchwork)
})

# ==========================================
# 老师指定的作图模板 (全局加黑加深)
# ==========================================
custom_theme <- theme_bw() +
    theme(
        text = element_text(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(linewidth = 0.2, color = "black"),
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.3, color = "black"),
        axis.line.x.top = element_blank(),
        axis.line.y.right = element_blank(),
        strip.text = element_text(face = "plain", size = 8, color = "black"),
        strip.background = element_blank(),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.x = element_text(color = "black"),
        axis.title.y = element_text(color = "black"),
        plot.title = element_text(color = "black", face = "plain", size=10, hjust=0.5),
        legend.position = "none"
    )

out_dir <- "D:/桌面/Rstudio/大队列研究/age/out_figs"

# ==========================================
# Figure 4A: 绘制左侧流程图 (p_schematic)
# ==========================================
boxes <- data.frame(
    id = 1:7,
    x = c(50, 25, 75, 25, 75, 50, 50),
    y = c(90, 70, 70, 50, 50, 30, 10),
    w = c(45, 38, 38, 38, 38, 50, 45),
    h = c(12, 12, 12, 12, 12, 12, 12),
    
    label = c(
        "Discovery Cohort\n(6 Cohorts Iteratively)",
        "Model 1\nOriginal Abundance Matrix",
        "Model 2\nAge-Residualized Matrix",
        "Feature Set A\n(Original Top 5% & 10%)",
        "Feature Set B\n(Residualized Top 5% & 10%)",
        "Pooled Validation Cohort\n(5 Independent Cohorts)",
        "Diagnostic Evaluation\n(Cross-Cohort AUC Comparison)"
    ),
    
    fill = c("#F4FAFE", "#FFF5E6", "#E6F4FB", "#E69F00", "#56B4E9", "#F4FAFE", "#EAEAEA"),
    color = c("#333333", "#E69F00", "#56B4E9", "black", "black", "#333333", "#333333"),
    text_col = c("black", "black", "black", "white", "white", "black", "black"),
    font_face = c("plain", "plain", "plain", "plain", "plain", "plain", "plain")
)

venn_box <- data.frame(
    x = 50, y = 50, w = 15, h = 8,
    label = "Overlap\nAnalysis",
    fill = "white", color = "gray50", text_col = "black", font_face = "italic"
)

arrows <- data.frame(
    x = c(50, 50, 25, 75, 25, 75, 50),
    xend = c(25, 75, 25, 75, 50, 50, 50),
    y = c(84, 84, 64, 64, 44, 44, 24),
    yend = c(76, 76, 56, 56, 36, 36, 16)
)

p_schematic <- ggplot() +
    geom_segment(data = arrows, aes(x = x, y = y, xend = xend, yend = yend),
                 arrow = arrow(length = unit(0.3, "cm"), type = "closed"), 
                 linewidth = 1, color = "gray30") +
    geom_segment(aes(x = 35, y = 50, xend = 42, yend = 50), linetype = "dashed",
                 arrow = arrow(length = unit(0.2, "cm"), type = "closed"), color = "gray50", linewidth = 0.5) +
    geom_segment(aes(x = 65, y = 50, xend = 58, yend = 50), linetype = "dashed",
                 arrow = arrow(length = unit(0.2, "cm"), type = "closed"), color = "gray50", linewidth = 0.5) +
    geom_rect(data = boxes, aes(xmin = x - w/2, xmax = x + w/2, 
                                ymin = y - h/2, ymax = y + h/2, fill = fill, color = color),
              linewidth = 1) +
    geom_rect(data = venn_box, aes(xmin = x - w/2, xmax = x + w/2, 
                                   ymin = y - h/2, ymax = y + h/2, fill = fill, color = color),
              linewidth = 0.5, linetype = "dashed") +
    scale_fill_identity() +
    scale_color_identity() +
    # 【微调1】缩放了主流程图的文字大小（从 4.5 降到 3.5），稍微收紧行距（0.95），确保不超框
    geom_text(data = boxes, aes(x = x, y = y, label = label, color = text_col, fontface = font_face), 
              size = 3.5, lineheight = 0.95, family = "sans") +
    geom_text(data = venn_box, aes(x = x, y = y, label = label, color = text_col, fontface = font_face), 
              size = 3.0, lineheight = 0.95, family = "sans") +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
    theme_void() + 
    theme(plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"))


# ==============================================================================
# Figure 4B-G: 6队列轮转跨队列模型验证 (包含 50 次迭代)
# ==============================================================================
base_dir <- "D:/桌面/Rstudio/大队列研究/age"
meta_dir <- file.path(base_dir, "metadata")
prof_dir <- file.path(base_dir, "metaphlan4_profiles")

all_cohorts <- c(
    "This_study_cohort5__NSHII", "Public_study__YangJ_2020",
    "Public_study__LiuNN_2022", "This_study_cohort3__ONCOBIOME_IIGM_CZ", 
    "Public_study__GuptaA_2019", "Public_study__ThomasAM_2018b"
)

n_iterations <- 50

# 初始化一个列表，用来装右侧的 6 个箱线图
plot_list_bg <- list()

for (discovery_cohort in all_cohorts) {
    
    clean_name <- str_remove_all(discovery_cohort, "Public_study__|This_study_")
    message(sprintf("\n=======================================================\n🚀 正在启动流水线... 发现集切换为: %s\n=======================================================", clean_name))
    
    val_cohorts <- setdiff(all_cohorts, discovery_cohort)
    meta_file <- file.path(meta_dir, paste0(discovery_cohort, ".tsv"))
    profile_file <- file.path(prof_dir, paste0(discovery_cohort, ".tsv"))
    
    # --- [Phase 1] 发现集 (Discovery) 模型训练与特征提取 ---
    meta_raw <- read_tsv(meta_file, show_col_types = FALSE)
    s_col <- grep("(?i)^(Run|Sample|Sample_ID|SampleID|Name)$", colnames(meta_raw), value = TRUE)[1]
    d_col <- grep("(?i)^(Disease|Study\\.Group|Group)$", colnames(meta_raw), value = TRUE)[1]
    a_col <- grep("(?i)^Age$", colnames(meta_raw), value = TRUE)[1]
    
    meta <- meta_raw %>%
        select(all_of(c(s_col, d_col, a_col))) %>%
        rename(Run = !!sym(s_col), Disease = !!sym(d_col), Age = !!sym(a_col)) %>%
        mutate(Disease = case_when(grepl("(?i)crc|cancer|case|tumor|carcinoma", Disease) ~ "CRC", grepl("(?i)control|healthy|normal", Disease) ~ "healthy", TRUE ~ NA_character_)) %>% 
        drop_na(Disease, Age) %>% mutate(Disease = as.factor(Disease))
    
    prof <- read_tsv(profile_file, show_col_types = FALSE)
    prof_species <- prof %>% filter(grepl("s__", clade_name) & !grepl("t__", clade_name)) %>%
        mutate(clade_name = str_extract(clade_name, "s__.*")) %>% select(-any_of("NCBI_tax_id"))
    
    prof_t <- prof_species %>% pivot_longer(-clade_name, names_to = "Run", values_to = "Abundance") %>%
        pivot_wider(names_from = clade_name, values_from = Abundance, values_fill = 0)
    
    merged_data <- inner_join(meta, prof_t, by = "Run")
    pheno <- merged_data %>% select(Run, Disease, Age)
    features <- merged_data %>% select(-Run, -Disease, -Age)
    
    prev_threshold <- 0.10 * nrow(features)
    features_filtered <- features[, colSums(features > 0) >= prev_threshold]
    
    message(">> 计算年龄残差...")
    features_resid <- as.data.frame(matrix(NA, nrow = nrow(features_filtered), ncol = ncol(features_filtered)))
    rownames(features_resid) <- rownames(features_filtered)
    colnames(features_resid) <- colnames(features_filtered)
    for (taxa in colnames(features_filtered)) {
        model <- lm(features_filtered[[taxa]] ~ pheno$Age)
        features_resid[[taxa]] <- residuals(model)
    }
    gc() 
    
    set.seed(123) 
    message(">> 训练 全局Original模型 提取特征...")
    rf_orig_full <- randomForest(x = features_filtered, y = pheno$Disease, ntree = 500, importance = TRUE)
    imp_orig <- importance(rf_orig_full, type = 2)[, 1]
    imp_orig <- sort(imp_orig, decreasing = TRUE)
    
    message(">> 训练 全局Residualized模型 提取特征...")
    rf_resid_full <- randomForest(x = features_resid, y = pheno$Disease, ntree = 500, importance = TRUE)
    imp_resid <- importance(rf_resid_full, type = 2)[, 1]
    imp_resid <- sort(imp_resid, decreasing = TRUE)
    
    n_features <- length(imp_orig)
    top05_n <- max(1, round(n_features * 0.05))
    top10_n <- max(1, round(n_features * 0.10))
    
    features_list <- list(
        Original_Top05 = names(imp_orig)[1:top05_n], Original_Top10 = names(imp_orig)[1:top10_n],
        Resid_Top05    = names(imp_resid)[1:top05_n], Resid_Top10    = names(imp_resid)[1:top10_n],
        Imp_Orig_Full  = imp_orig, Imp_Resid_Full = imp_resid
    )
    
    cutoffs <- list(
        "Top05" = list(Orig = features_list$Original_Top05, Resid = features_list$Resid_Top05),
        "Top10" = list(Orig = features_list$Original_Top10, Resid = features_list$Resid_Top10)
    )
    
    auc_results_all <- data.frame() 
    
    for (cutoff_name in names(cutoffs)) {
        message(sprintf("\n--- 正在生成 %s 的结果 ---", cutoff_name))
        set_orig <- cutoffs[[cutoff_name]]$Orig
        set_resid <- cutoffs[[cutoff_name]]$Resid
        
        # [Phase 2 & 3: Venn图与哑铃图生成略过...]
        
        # ==========================================================
        # 🚨 [Phase 4] 50 次种子迭代 + 跨队列超级验证
        # ==========================================================
        message(sprintf("   >> 正在进行 50 次随机种子迭代验证 (Cutoff: %s) ...", cutoff_name))
        
        val_data_list <- list()
        all_needed_taxa <- unique(c(set_orig, set_resid))
        for (c_id in val_cohorts) {
            val_meta_raw <- read_tsv(file.path(meta_dir, paste0(c_id, ".tsv")), show_col_types = F)
            s_col <- grep("(?i)^(Run|Sample|Sample_ID|SampleID|Name)$", colnames(val_meta_raw), value = TRUE)[1]
            d_col <- grep("(?i)^(Disease|Study\\.Group|Group)$", colnames(val_meta_raw), value = TRUE)[1]
            val_meta <- val_meta_raw %>% select(all_of(c(s_col, d_col))) %>% rename(Run = !!sym(s_col), Disease = !!sym(d_col)) %>% mutate(Disease = case_when(grepl("(?i)crc|cancer|case|tumor|carcinoma", Disease) ~ "CRC", grepl("(?i)control|healthy|normal", Disease) ~ "healthy", TRUE ~ NA_character_)) %>% drop_na(Disease) %>% mutate(Disease = as.factor(Disease))
            val_prof <- read_tsv(file.path(prof_dir, paste0(c_id, ".tsv")), show_col_types = F) %>% filter(grepl("s__", clade_name) & !grepl("t__", clade_name)) %>% mutate(clade_name = str_extract(clade_name, "s__.*")) %>% select(-any_of("NCBI_tax_id")) %>% filter(clade_name %in% all_needed_taxa)
            val_prof_t <- val_prof %>% pivot_longer(-clade_name, names_to = "Run", values_to = "Abundance") %>% pivot_wider(names_from = clade_name, values_from = Abundance, values_fill = 0)
            
            merged_val <- inner_join(val_meta, val_prof_t, by = "Run")
            for(tx in set_orig) if(!tx %in% colnames(merged_val)) merged_val[[tx]] <- 0
            for(tx in set_resid) if(!tx %in% colnames(merged_val)) merged_val[[tx]] <- 0
            val_data_list[[c_id]] <- merged_val
        }
        
        # 50 次迭代核心循环
        for (seed_val in 1:n_iterations) {
            set.seed(seed_val)
            rf_orig_panel <- randomForest(x = features_filtered[, set_orig, drop = FALSE], y = pheno$Disease, ntree = 500)
            rf_resid_panel <- randomForest(x = features_resid[, set_resid, drop = FALSE], y = pheno$Disease, ntree = 500)
            
            super_pred_orig <- numeric(); super_pred_resid <- numeric(); super_labels <- character()
            
            for (c_id in val_cohorts) {
                val_data <- val_data_list[[c_id]]
                val_y <- val_data$Disease; val_x <- val_data %>% select(-Run, -Disease)
                pred_orig <- predict(rf_orig_panel, newdata = val_x[, set_orig, drop = FALSE], type = "prob")[, "CRC"]
                pred_resid <- predict(rf_resid_panel, newdata = val_x[, set_resid, drop = FALSE], type = "prob")[, "CRC"]
                super_pred_orig <- c(super_pred_orig, pred_orig); super_pred_resid <- c(super_pred_resid, pred_resid); super_labels <- c(super_labels, as.character(val_y))
            }
            
            auc_orig_val <- as.numeric(roc(super_labels, super_pred_orig, levels = c("healthy", "CRC"), direction = "<", quiet = T)$auc)
            auc_resid_val <- as.numeric(roc(super_labels, super_pred_resid, levels = c("healthy", "CRC"), direction = "<", quiet = T)$auc)
            
            auc_results_all <- rbind(auc_results_all, 
                                     data.frame(Iteration = seed_val, Cutoff = cutoff_name, Model = "Original Data", AUC = auc_orig_val),
                                     data.frame(Iteration = seed_val, Cutoff = cutoff_name, Model = "Age-Residualized Data", AUC = auc_resid_val)
            )
        }
    }
    
    # 保存数据
    auc_results_all$Model <- factor(auc_results_all$Model, levels = c("Original Data", "Age-Residualized Data"))
    auc_results_all$Cutoff <- factor(auc_results_all$Cutoff, levels = c("Top05", "Top10"), labels = c("Top 5% Features", "Top 10% Features"))
    
    write_csv(auc_results_all, file.path(out_dir, paste0("Figure6_SuperVal_AUC_Data_", clean_name, ".csv")))
    
    # --- [Phase 5] 绘制该 Discovery 队列专属的统计学箱线图 ---
    p_auc <- ggplot(auc_results_all, aes(x = Cutoff, y = AUC, fill = Model)) +
        geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85, color = "black", linewidth = 0.3) +
        geom_point(aes(color = Model), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.6), size = 2, alpha = 0.4) +
        scale_fill_manual(values = c("Original Data" = "#E69F00", "Age-Residualized Data" = "#56B4E9")) +
        scale_color_manual(values = c("Original Data" = "#CC8500", "Age-Residualized Data" = "#0072B2")) +
        
        stat_compare_means(aes(group = Model), paired = TRUE, method = "wilcox.test", 
                           label = "p.signif", size = 5, vjust = 0.5) +
        
        scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
        
        # 【微调2】简化 Y 轴坐标标题为 "Validation AUC"
        labs(x = NULL, y = "Validation AUC", 
             title = clean_name, fill = NULL, color = NULL) +
        custom_theme + # 套用老师的黑字模板
        theme(
            plot.title = element_text(face = "plain", hjust = 0.5, size = 11, margin = ggplot2::margin(b = 5)), 
            axis.text.x = element_text(face = "plain", size = 9, color = "black"), 
            legend.position = "top",
            legend.text = element_text(size = 11),
            legend.margin = ggplot2::margin(t=0, r=0, b=0, l=0, unit="mm"),
            plot.margin = ggplot2::margin(t = 10, r = 5, b = 5, l = 5, unit = "pt")
        ) +
        guides(
            fill = "none",
            color = guide_legend(override.aes = list(size = 4, alpha = 1))
        )
    
    # 将图片存入列表，不再单独保存 p_auc
    plot_list_bg[[clean_name]] <- p_auc
    message(sprintf("✅ %s 轮次执行完毕 (已完成 50 次迭代验证，数据已保存至 CSV)！", clean_name))
}


# ==============================================================================
# 终极拼图：A拉长一点放在最左边，B-G一共6个小图，分成3个为一列，共两列在右侧
# ==============================================================================
# 组装右侧 B-G (2列3行)，收集图注放置在顶部
p_right <- wrap_plots(plot_list_bg, ncol = 2) + 
    plot_layout(guides = "collect") & 
    theme(legend.position = "top")

# A 占 1 份宽，右侧 B-G 占 2 份宽。A 自然会被纵向拉长。
p_combined <- p_schematic + p_right + plot_layout(widths = c(1, 2))

ggsave(file.path(out_dir, "Figure4_Combined_Flow_and_Boxplots.pdf"), p_combined, width = 12, height = 8, scale = 0.8)
ggsave(file.path(out_dir, "Figure4_Combined_Flow_and_Boxplots.png"), p_combined, width = 12, height = 8, dpi = 600)

message("\n✅ Figure 4A (流程图) 和 Figure 4B-G (箱图) 成功合并！排版细节已彻底完善。")