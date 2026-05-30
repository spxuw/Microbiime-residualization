suppressPackageStartupMessages({
    library(readr); library(dplyr); library(stringr); library(purrr); library(tidyr); library(ggplot2); library(patchwork)
})

auroc_csv   <- "D:/桌面/Rstudio/大队列研究/age/auroc.csv"
out_dir     <- "D:/桌面/Rstudio/大队列研究/age/out_figs"

# ==========================================
# 老师指定的作图模板 (全局文字加黑加深)
# ==========================================
custom_theme <- theme_bw() +
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
        strip.text = element_text(face = "plain", size = 10.5, color = "black"),
        strip.background = element_blank(),
        axis.text.x = element_text(size = 9, color = "black", face = "plain"),
        axis.text.y = element_text(size = 9, color = "black", face = "plain"),
        axis.title.x = element_text(size = 9.5, color = "black", face = "plain"),
        axis.title.y = element_text(size = 9.5, color = "black", face = "plain"),
        plot.title = element_blank(),
        legend.position = "none" # 默认关闭，后续单独强行开启
    )

# =========================
# 1. 读取 1B 的权威清单，作为绝对基准
# =========================
proj_list <- read_csv(file.path(out_dir, "Figure1_projects_final.csv"), show_col_types = FALSE)

# =========================
# 2. 动态读取 auroc.csv 并提取列名
# =========================
au <- read_csv(auroc_csv, show_col_types=FALSE)

col_ds <- grep("dataset|file|data", names(au), value=TRUE, ignore.case=TRUE)[1]
col_bp <- grep("BioProject|project|study|ID", names(au), value=TRUE, ignore.case=TRUE)[1]

col_null <- grep("auc_null.*mean|null_auc.*mean|pheno_auc_null_mean", names(au), value=TRUE, ignore.case=TRUE)[1]
col_rf_full <- grep("full.*rf.*auc.*mean|pheno_auc_full_rf_mean", names(au), value=TRUE, ignore.case=TRUE)[1]
col_xgb_full <- grep("full.*xgb.*auc.*mean|pheno_auc_full_xgb_mean", names(au), value=TRUE, ignore.case=TRUE)[1]
col_rf_pperm <- grep("rf.*perm.*p|pheno.*delta.*rf.*perm.*p", names(au), value=TRUE, ignore.case=TRUE)[1]
col_xgb_pperm <- grep("xgb.*perm.*p|pheno.*delta.*xgb.*perm.*p", names(au), value=TRUE, ignore.case=TRUE)[1]

# =========================
# 3. 数据整理与【强行对齐 1B】
# =========================
df <- au |>
    mutate(
        dataset = as.character(.data[[col_ds]]), 
        study = as.character(.data[[col_bp]]), 
        dataset_base = str_replace(dataset, "\\.csv$", ""),
        source = case_when(
            dataset_base == "CRC" ~ "CRC", 
            str_detect(dataset_base, "^GMrepo_") ~ "GM", 
            str_detect(dataset_base, "^16S_") ~ "Res", 
            TRUE ~ NA_character_
        ),
        dataset_key = case_when(source == "GM" ~ "GMrepo", TRUE ~ dataset_base),
        study_label = if_else(source == "Res", paste0(study, " [", str_remove(dataset_key, "^16S_"), "]"), study)
    ) |> 
    right_join(proj_list |> dplyr::select(source, study_label, order), by=c("source","study_label")) |>
    mutate(color = case_when(source == "CRC" ~ "#C97DBB", source == "GM"  ~ "#E69F00", source == "Res" ~ "#5CB8E8"))

# =========================
# 4. 构建作图长数据
# =========================
make_long <- function(model){
    full_c <- if(model=="RF") col_rf_full else col_xgb_full
    p_c    <- if(model=="RF") col_rf_pperm else col_xgb_pperm
    
    d <- df |> transmute(
        source, study_label, color, order, 
        auc_null = if(is.na(col_null) || is.na(full_c)) NA_real_ else as.numeric(.data[[col_null]]), 
        auc_full = if(is.na(full_c)) NA_real_ else as.numeric(.data[[full_c]]), 
        p_perm   = if(is.na(p_c)) NA_real_ else as.numeric(.data[[p_c]])
    )
    
    list(
        auc = d |> pivot_longer(cols=c(auc_null, auc_full), names_to="part", values_to="auc", values_drop_na = FALSE) |> 
            mutate(part = recode(part, auc_null="Age-only (Null)", auc_full="Taxa+Age (Full)"), model = model),
        p   = d |> mutate(model=model, neglog10p = ifelse(!is.na(p_perm) & p_perm>0, -log10(p_perm), NA_real_))
    )
}

rf <- make_long("RF"); xgb <- make_long("XGBoost")
auc_all <- bind_rows(rf$auc, xgb$auc); p_all <- bind_rows(rf$p, xgb$p)

# 图中显示用名称：内部仍保留 CRC/GM/Res 作为 join/filter key
source_title <- c("CRC" = "CRC", "GM" = "GMrepo", "Res" = "ResMicroDB")

# =========================
# 5. 画图模块 (换装 custom_theme)
# =========================
plot_auc_one_source <- function(model_name, src){
    src_levels <- proj_list |> filter(source == src) |> arrange(order) |> pull(study_label)
    dd <- auc_all |> filter(model==model_name, source==src) |> mutate(study_label = factor(study_label, levels = src_levels))
    
    ggplot(dd, aes(x=study_label, y=auc, group=study_label)) +
        geom_line(aes(color=color), linewidth=0.6, alpha=0.9, show.legend=FALSE, na.rm=TRUE) +
        geom_point(aes(shape=part), size=2.6, color="grey20", fill="white", stroke=0.35, na.rm=TRUE) +
        scale_shape_manual(values=c("Age-only (Null)"=21, "Taxa+Age (Full)"=22), name=NULL) +
        scale_color_identity() + 
        scale_x_discrete(drop = FALSE) + 
        coord_flip(clip="off") + 
        labs(y="AUROC", x=NULL, title=NULL) + 
        custom_theme
}

plot_sig_one_source <- function(model_name, src){
    src_levels <- proj_list |> filter(source == src) |> arrange(order) |> pull(study_label)
    dd <- p_all |> filter(model==model_name, source==src) |> mutate(study_label = factor(study_label, levels = src_levels))
    
    ggplot(dd, aes(x=study_label, y=neglog10p)) +
        geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey60", linewidth=0.7) +
        geom_point(aes(color=color), size=2.6, alpha=0.9, show.legend=FALSE, na.rm=TRUE) +
        scale_color_identity() + 
        scale_x_discrete(drop = FALSE) + 
        coord_flip(clip="off") + 
        labs(y=expression(-log[10](p)), x=NULL, title=NULL) + 
        custom_theme
}

# =========================
# 6. 分别组装出图 (强制图注在对应的第一行底部)
# =========================

# --- 1. Random Forest 单独导出 ---
rf_auc_row <- (plot_auc_one_source("RF","CRC") | plot_auc_one_source("RF","GM") | plot_auc_one_source("RF","Res")) + 
    plot_layout(guides = "collect") & 
    theme(
        legend.position = "bottom",                     # 强行在此处把图注放到底部
        legend.title = element_blank(),
        legend.text = element_text(size=10.5, color="black", face="plain"),
        legend.key.size = unit(3.5, "mm")         # 略微压缩图注间距防占位
    )

rf_sig_row <- plot_sig_one_source("RF","CRC") | plot_sig_one_source("RF","GM") | plot_sig_one_source("RF","Res")

fig2_rf <- rf_auc_row / rf_sig_row

# 严格执行你要求的导出参数
ggsave(file.path(out_dir, "Figure2_RandomForest_Only.pdf"), fig2_rf, width=11, height=7, scale=0.7, dpi = 600)
ggsave(file.path(out_dir, "Figure2_RandomForest_Only.png"), fig2_rf, width=12, height=5.6, dpi=600)


# --- 2. XGBoost 单独导出 ---
xgb_auc_row <- (plot_auc_one_source("XGBoost","CRC") | plot_auc_one_source("XGBoost","GM") | plot_auc_one_source("XGBoost","Res")) + 
    plot_layout(guides = "collect") & 
    theme(
        legend.position = "bottom",                     # 强行在此处把图注放到底部
        legend.title = element_blank(),
        legend.text = element_text(size=10.5, color="black", face="plain"),
        legend.key.size = unit(3.5, "mm"),
        legend.margin = ggplot2::margin(1,1,1,1, "mm")           # 略微压缩图注间距防占位
    )

xgb_sig_row <- plot_sig_one_source("XGBoost","CRC") | plot_sig_one_source("XGBoost","GM") | plot_sig_one_source("XGBoost","Res")

fig2_xgb <- xgb_auc_row / xgb_sig_row

# 严格执行你要求的导出参数
ggsave(file.path(out_dir, "Figure2_XGBoost_Only.pdf"), fig2_xgb, width=12, height=5.6, scale=0.8)
ggsave(file.path(out_dir, "Figure2_XGBoost_Only.png"), fig2_xgb, width=12, height=5.6, dpi=600)

message("Figure 2 绘制完毕！模板、字体和保存参数均已更新，内部逻辑原封未动。")