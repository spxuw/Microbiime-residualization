suppressPackageStartupMessages({
    library(readr); library(dplyr); library(stringr); library(purrr)
    library(tidyr); library(ggplot2); library(patchwork)
})

# ==============================================================================
# 老师指定的导出模板 (【修改1】：加入 color="black" 强制全图文字变黑变深)
# ==============================================================================
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
        strip.text = element_text(face = "plain", size = 9, color = "black"),
        strip.background = element_blank(),
        axis.text.x = element_text(size = 9, color = "black", face = "plain"),
        axis.text.y = element_text(size = 9, color = "black", face = "plain"),
        axis.title.x = element_text(size = 9.5, color = "black", face = "plain"), 
        axis.title.y = element_text(size = 9.5, color = "black", face = "plain"),
        plot.title = element_blank(),
        legend.position = "none"
    )

# =========================
# Figure 1B: 以 auroc.csv 结果为基准过滤 Project
# =========================
summary_dir <- "D:/桌面/Rstudio/大队列研究/age/summary"
auroc_csv   <- "D:/桌面/Rstudio/大队列研究/age/auroc.csv"
out_dir     <- "D:/桌面/Rstudio/大队列研究/age/out_figs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

N_SHOW <- 13
PANEL_W <- 5.4
PANEL_H <- 0.12 * N_SHOW + 2.45
BAR_THICKNESS <- 0.52
YTEXT_SIZE <- 10
XTEXT_SIZE <- 11
LEG_TEXT_SIZE <- 9

COL_BLUE   <- "#5CB8E8"  # ResMicroDB
COL_ORANGE <- "#E69F00"  # GMrepo
COL_PURPLE <- "#C97DBB"  # CRC

# ---- 1. 读取 AUROC 文件构建“白名单” ----
get_colname <- function(df, cands){ hit <- cands[cands %in% names(df)]; hit[1] }

au <- read_csv(auroc_csv, show_col_types=FALSE)
col_ds <- get_colname(au, c("dataset","Dataset","file","data"))
col_bp <- get_colname(au, c("BioProject","bioproject","project","study","Study","ID"))

valid_au_studies <- au |>
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
    filter(!is.na(source), !is.na(study)) |>
    distinct(source, study_label)

# ---- 2. 读 Summary 并 Harmonize ----
detect_delim <- function(file){
    line <- readLines(file, n = 1, warn = FALSE)
    counts <- c(comma=str_count(line, ","), semi=str_count(line, ";"), tab=str_count(line, "\t"))
    nm <- names(which.max(counts))
    if (nm=="tab") return("\t")
    if (nm=="semi") return(";")
    ","
}
read_summary_any <- function(file){
    delim <- detect_delim(file)
    h <- read_delim(file, delim=delim, col_names=FALSE, n_max=1, show_col_types=FALSE)
    r <- read_delim(file, delim=delim, col_names=FALSE, n_max=1, skip=1, show_col_types=FALSE)
    if (ncol(r) == ncol(h) + 1) {
        header <- as.character(unlist(h[1, ]))
        coln <- c("row_id", header)
        dat <- read_delim(file, delim=delim, skip=1, col_names=coln, show_col_types=FALSE) |> select(-row_id)
    } else {
        dat <- read_delim(file, delim=delim, show_col_types=FALSE)
        if ("...1" %in% names(dat)) dat <- dat |> select(-...1)
    }
    dat
}
parse_n <- function(x){ suppressWarnings(as.numeric(str_replace_all(as.character(x), "[^0-9.+-]", ""))) }

is_control_label <- function(x){
    xl <- str_to_lower(str_squish(as.character(x)))
    xl %in% c("control", "health study participation") | str_detect(xl, "(^|\\b)(healthy|normal|health)(\\b|$)") | str_detect(xl, "\\bhc\\b")
}
harmonize <- function(x){
    xl <- str_to_lower(str_trim(x))
    if (is_control_label(x)) return("Control")
    if (str_detect(xl, "adenoma")) return("Adenoma")
    if (str_detect(xl, "carcinoma|colorectal|crc|cancer")) return("CRC")
    if (str_detect(xl, "covid")) return("COVID-19")
    if (str_detect(xl, "\\basthma\\b")) return("Asthma")
    if (str_detect(xl, "chronic rhinosinusitis with nasal polyps")) return("CRSwNP")
    if (str_detect(xl, "chronic rhinosinusitis without nasal polyps")) return("CRSsNP")
    if (str_detect(xl, "respiratory syncytial virus|\\brsv\\b")) return("RSV infection")
    if (str_detect(xl, "allergic rhinitis")) return("Allergic rhinitis")
    if (str_detect(xl, "pneumonia")) return("Pneumonia")
    if (str_detect(xl, "\\bhiv\\b")) return("HIV")
    if (str_detect(xl, "interstitial lung disease|\\bild\\b")) return("Interstitial lung disease")
    if (str_detect(xl, "crohn")) return("Crohn disease")
    if (str_detect(xl, "ulcerative colitis")) return("Ulcerative colitis")
    if (str_detect(xl, "multiple sclerosis")) return("Multiple sclerosis")
    if (str_detect(xl, "parkinson")) return("Parkinson's disease")
    if (str_detect(xl, "diabetes mellitus, type 1|type 1 diabetes|\\bt1d\\b")) return("Type 1 diabetes")
    if (str_detect(xl, "metabolic syndrome")) return("Metabolic syndrome")
    if (str_detect(xl, "obesity")) return("Obesity")
    str_to_sentence(str_replace_all(str_trim(x), "\\s+", " "))
}

sum_files <- list.files(summary_dir, pattern="_summary\\.csv$", full.names=TRUE, recursive=TRUE)
dat <- map_dfr(sum_files, function(f){
    ds <- str_replace(basename(f), "_summary\\.csv$", "")
    x  <- read_summary_any(f)
    if (all(c("BioProject","Phenotype") %in% names(x))) { x <- x |> rename(study = BioProject, phenotype_raw = Phenotype)
    } else if (all(c("Study.name","Disease") %in% names(x))) { x <- x |> rename(study = `Study.name`, phenotype_raw = Disease)
    } else if (all(c("project_id","phenotype") %in% names(x))) { x <- x |> rename(study = project_id, phenotype_raw = phenotype) }
    count_col <- intersect(c("n","N","count","Count","samples","Samples"), names(x))[1]
    tibble(dataset_base = ds, study = as.character(x$study), phenotype_raw = as.character(x$phenotype_raw), n = parse_n(x[[count_col]]))
}) |>
    mutate(
        source = case_when(str_detect(dataset_base, "^16S_") ~ "Res", dataset_base %in% c("GMrepo_genus","GMrepo_species") ~ "GM", dataset_base == "CRC" ~ "CRC", TRUE ~ "OtherSource"),
        dataset_key = case_when(source == "GM" ~ "GMrepo", TRUE ~ dataset_base),
        study_label = if_else(source == "Res", paste0(study, " [", str_remove(dataset_key, "^16S_"), "]"), study),
        phenotype = map_chr(phenotype_raw, harmonize)
    ) |> filter(source %in% c("CRC","GM","Res"), !is.na(n), n > 0) |>
    group_by(source, study_label, phenotype) |> summarise(n = max(n, na.rm=TRUE), .groups="drop")

# ---- 3. 双重过滤：必须在白名单内，且包含 Control 和 Case ----
dat2 <- dat |> 
    inner_join(valid_au_studies, by = c("source", "study_label")) |> 
    group_by(source, study_label) |>
    mutate(has_control = any(phenotype=="Control"), has_case = any(phenotype!="Control")) |>
    ungroup() |>
    filter(has_control, has_case)

# ---- 4. Plot ----
lighten <- function(col, amount = 0.28){ rgbv <- grDevices::col2rgb(col)/255; grDevices::rgb(rgbv[1,] + (1-rgbv[1,]) * amount, rgbv[2,] + (1-rgbv[2,]) * amount, rgbv[3,] + (1-rgbv[3,]) * amount) }
darken <- function(col, amount = 0.35){ rgbv <- grDevices::col2rgb(col)/255; grDevices::rgb(rgbv[1,] * (1-amount), rgbv[2,] * (1-amount), rgbv[3,] * (1-amount)) }
make_shades <- function(base_col, n){ grDevices::colorRampPalette(c(lighten(base_col,0.28), base_col, darken(base_col,0.35)))(max(n,1)) }

theme_compact <- function(){
    theme_classic(base_size = 13) + theme(
        plot.title = element_blank(), axis.text.y = element_text(size=YTEXT_SIZE, face="plain", color="black"), axis.text.x = element_text(size=XTEXT_SIZE, face="plain", color="black"),
        legend.position = c(0.985, 0.045), legend.justification = c(1, 0), legend.direction = "vertical", legend.box = "vertical",
        legend.title = element_blank(), legend.text  = element_text(size=LEG_TEXT_SIZE, face="plain", color="black"),
        legend.key.size = grid::unit(4.8, "mm"), legend.spacing.y = grid::unit(0.8, "mm"),
        legend.background = element_rect(fill = grDevices::adjustcolor("white", alpha.f = 0.82), colour = NA), plot.margin = ggplot2::margin(3, 10, 5, 3, unit = "mm")
    )
}

build_bystudy_plot <- function(df_source, common_list, title_text, base_fill, top_n=8, label_map){
    top_groups <- df_source |> filter(phenotype %in% common_list) |> group_by(phenotype) |> summarise(N=sum(n)) |> arrange(desc(N)) |> slice_head(n=top_n) |> pull(phenotype)
    df_plot <- df_source |> mutate(group = case_when(phenotype == "Control" ~ "Control", phenotype %in% top_groups ~ phenotype, TRUE ~ "Other")) |>
        group_by(study_label, group) |> summarise(N = sum(n), .groups="drop")
    keep_studies <- df_plot |> group_by(study_label) |> summarise(Total=sum(N)) |> arrange(desc(Total)) |> slice_head(n=N_SHOW) |> pull(study_label)
    df_plot <- df_plot |> filter(study_label %in% keep_studies)
    ord <- df_plot |> group_by(study_label) |> summarise(Total=sum(N)) |> arrange(Total) |> pull(study_label)
    
    present_groups <- df_plot |> group_by(group) |> summarise(T=sum(N)) |> arrange(desc(T)) |> pull(group)
    present_groups <- c(intersect("Control", present_groups), setdiff(present_groups, c("Control","Other")), intersect("Other", present_groups))
    df_plot <- df_plot |> mutate(study_label = factor(study_label, levels=ord), group = factor(group, levels=present_groups))
    
    diseases <- setdiff(present_groups, c("Control","Other"))
    dis_cols <- make_shades(base_fill, length(diseases)); names(dis_cols) <- diseases
    pal <- c("Control"="#BDBDBD", dis_cols, "Other"="#E5E5E5")
    
    lab_vec <- label_map[present_groups]; lab_vec[is.na(lab_vec)] <- present_groups[is.na(lab_vec)]
    p <- ggplot(df_plot, aes(x=study_label, y=N, fill=group)) + geom_col(width=BAR_THICKNESS) + coord_flip(clip="off") +
        scale_fill_manual(values=pal, drop=FALSE, labels=lab_vec) + scale_y_continuous(labels=scales::comma, expand = expansion(mult = c(0, 0.14))) +
        labs(x=NULL, y="Sample size (per study)", title=NULL) + theme_compact() + guides(fill = guide_legend(ncol=1, byrow=TRUE))
    list(plot=p, ord=ord, pal=pal)
}

res_crc <- build_bystudy_plot(dat2 |> filter(source=="CRC"), c("Adenoma","CRC"), "CRC cohorts", COL_PURPLE, 2, c("Control"="Control", "Adenoma"="Adenoma", "CRC"="CRC", "Other"="Other"))
res_gm <- build_bystudy_plot(dat2 |> filter(source=="GM"), c("Crohn disease","Ulcerative colitis","Multiple sclerosis","Type 1 diabetes","Parkinson's disease","Metabolic syndrome","Asthma","HIV","Obesity"), "GMrepo metagenome", COL_ORANGE, 8, c("Control"="Control","Crohn disease"="Crohn","Ulcerative colitis"="UC","Multiple sclerosis"="MS","Type 1 diabetes"="T1D","Parkinson's disease"="Parkinson","Metabolic syndrome"="MetS","Asthma"="Asthma","HIV"="HIV","Obesity"="Obesity","Other"="Other"))
res_res <- build_bystudy_plot(dat2 |> filter(source=="Res"), c("COVID-19","Asthma","CRSwNP","CRSsNP","RSV infection","Allergic rhinitis","Pneumonia","HIV","Interstitial lung disease"), "ResMicroDB 16S (airway)", COL_BLUE, 8, c("Control"="Control","COVID-19"="COVID-19","Asthma"="Asthma","CRSwNP"="CRSwNP","CRSsNP"="CRSsNP","RSV infection"="RSV","Allergic rhinitis"="AR","Pneumonia"="Pneumonia","HIV"="HIV","Interstitial lung disease"="ILD","Other"="Other"))

p_crc <- res_crc$plot; p_gm <- res_gm$plot; p_res <- res_res$plot

if(!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)

# =======================================================
# 1B 导出部分 (【修改2】：压缩图注大小，防止挡住柱子)
# =======================================================
p_all_B <- (p_crc | p_gm | p_res) + 
    plot_layout(guides="keep") & 
    custom_theme & 
    theme(
        legend.position = c(0.985, 0.045), 
        legend.justification = c(1, 0),
        legend.title = element_blank(),                                    # 去掉图注标题，省空间
        legend.key.size = unit(3.5, "mm"),                                   # 图注色块略增大
        legend.text = element_text(size = 8.2, color = "black", face = "plain"),           # 把图注字号压缩，并加黑
        legend.spacing.y = unit(0.5, "mm"),                                # 缩小行距
        legend.margin = ggplot2::margin(1, 1, 1, 1, "mm"),                          # 去掉图注内部的多余边距
        legend.background = element_rect(fill = alpha("white", 0.7), colour = NA) # 半透明背景防止彻底挡图
    )

ggsave(file.path(out_dir, "Figure1B_3panels_AUROC_Filtered.pdf"), p_all_B, width=12, height=3.2, scale=0.8)
ggsave(file.path(out_dir, "Figure1B_3panels_AUROC_Filtered.png"), p_all_B, width=12, height=3.2, dpi=600)

# Export lists for 1C and 2
proj_list <- bind_rows(
    tibble(source="CRC", study_label=res_crc$ord) |> mutate(order=row_number()),
    tibble(source="GM",  study_label=res_gm$ord)  |> mutate(order=row_number()),
    tibble(source="Res", study_label=res_res$ord) |> mutate(order=row_number())
)
write_csv(proj_list, file.path(out_dir, "Figure1_projects_final.csv"))

pal_all <- bind_rows(
    tibble(source="CRC", group=names(res_crc$pal), color_hex=unname(res_crc$pal)),
    tibble(source="GM",  group=names(res_gm$pal),  color_hex=unname(res_gm$pal)),
    tibble(source="Res", group=names(res_res$pal), color_hex=unname(res_res$pal))
)
write_csv(pal_all, file.path(out_dir, "Figure1_palette_all.csv"))

message("Done. Figure 1B now ONLY includes projects that successfully exported AUROC results.")



# =========================
# Figure 1C: Figure1C_age_boxplot_FINAL
# =========================

suppressPackageStartupMessages({
    library(readr); library(dplyr); library(stringr); library(purrr); library(tidyr); library(ggplot2)
})

# =========================
# Paths & settings
# =========================
summary_dir <- "D:/桌面/Rstudio/大队列研究/age/summary"
out_dir     <- "D:/桌面/Rstudio/大队列研究/age/out_figs"

N_SHOW <- 13
PANEL_W <- 5.4
PANEL_H <- 0.20 * N_SHOW + 2.0
WHISKER_LW <- 0.70; BOX_OUTLINE_LW <- 0.25; MEDIAN_LW <- 0.35; DODGE_W <- 0.72; BOX_W <- 0.52; CAP_W <- 0.18

# =========================
# 1. 读取 1B 的基准清单与颜色
# =========================
proj_list <- read_csv(file.path(out_dir, "Figure1_projects_final.csv"), show_col_types = FALSE)
pal_df <- read_csv(file.path(out_dir, "Figure1_palette_all.csv"), show_col_types = FALSE)
global_pal <- setNames(pal_df$color_hex, pal_df$group)

# =========================
# 2. 文件读取与解析
# =========================
detect_delim <- function(file){
    line <- readLines(file, n = 1, warn = FALSE)
    counts <- c(comma=str_count(line, ","), semi=str_count(line, ";"), tab=str_count(line, "\t"))
    nm <- names(which.max(counts))
    if (nm=="tab") return("\t")
    if (nm=="semi") return(";")
    ","
}

read_summary_any <- function(file){
    delim <- detect_delim(file)
    h <- read_delim(file, delim=delim, col_names=FALSE, n_max=1, show_col_types=FALSE)
    r <- read_delim(file, delim=delim, col_names=FALSE, n_max=1, skip=1, show_col_types=FALSE)
    if (ncol(r) == ncol(h) + 1) {
        header <- as.character(unlist(h[1, ]))
        coln <- c("row_id", header)
        dat <- read_delim(file, delim=delim, skip=1, col_names=coln, show_col_types=FALSE) |> select(-row_id)
    } else {
        dat <- read_delim(file, delim=delim, show_col_types=FALSE)
        if ("...1" %in% names(dat)) dat <- dat |> select(-...1)
    }
    dat
}

parse_num <- function(x){
    suppressWarnings(as.numeric(str_replace_all(as.character(x), "[^0-9.\\-+eE]", "")))
}

harmonize <- function(x){
    xl <- str_to_lower(str_trim(x))
    if (xl %in% c("control", "health study participation") | str_detect(xl, "(^|\\b)(healthy|normal|health)(\\b|$)")) return("Control")
    if (str_detect(xl, "adenoma")) return("Adenoma")
    if (str_detect(xl, "carcinoma|colorectal|crc|cancer")) return("CRC")
    if (str_detect(xl, "covid")) return("COVID-19")
    if (str_detect(xl, "\\basthma\\b")) return("Asthma")
    if (str_detect(xl, "chronic rhinosinusitis with nasal polyps")) return("CRSwNP")
    if (str_detect(xl, "chronic rhinosinusitis without nasal polyps")) return("CRSsNP")
    if (str_detect(xl, "respiratory syncytial virus|\\brsv\\b")) return("RSV infection")
    if (str_detect(xl, "allergic rhinitis")) return("Allergic rhinitis")
    if (str_detect(xl, "pneumonia")) return("Pneumonia")
    if (str_detect(xl, "\\bhiv\\b")) return("HIV")
    if (str_detect(xl, "interstitial lung disease|\\bild\\b")) return("Interstitial lung disease")
    if (str_detect(xl, "crohn")) return("Crohn disease")
    if (str_detect(xl, "ulcerative colitis")) return("Ulcerative colitis")
    if (str_detect(xl, "multiple sclerosis")) return("Multiple sclerosis")
    if (str_detect(xl, "parkinson")) return("Parkinson's disease")
    if (str_detect(xl, "diabetes mellitus, type 1|type 1 diabetes|\\bt1d\\b")) return("Type 1 diabetes")
    if (str_detect(xl, "metabolic syndrome")) return("Metabolic syndrome")
    if (str_detect(xl, "obesity")) return("Obesity")
    str_to_sentence(str_replace_all(str_trim(x), "\\s+", " "))
}

sum_files <- list.files(summary_dir, pattern="_summary\\.csv$", full.names=TRUE, recursive=TRUE)

all_long <- map_dfr(sum_files, function(f){
    x <- read_summary_any(f)
    ds_base <- str_replace(basename(f), "_summary\\.csv$", "")
    sc <- names(x)[names(x) %in% c("BioProject","project_id","Study.name")][1]
    gc <- names(x)[names(x) %in% c("Phenotype","phenotype","Disease")][1]
    nc <- names(x)[names(x) %in% c("n","N","count","Samples","samples")][1]
    
    tibble(
        dataset_base = ds_base,
        study = as.character(x[[sc]]),
        group_raw = as.character(x[[gc]]),
        n = parse_num(x[[nc]]),
        mean = if("mean" %in% names(x)) parse_num(x[["mean"]]) else NA_real_,
        sd = if("sd" %in% names(x)) parse_num(x[["sd"]]) else NA_real_,
        median = if("median" %in% names(x)) parse_num(x[["median"]]) else NA_real_,
        q25 = if("q25" %in% names(x)) parse_num(x[["q25"]]) else NA_real_,
        q75 = if("q75" %in% names(x)) parse_num(x[["q75"]]) else NA_real_
    )
}) |> 
    mutate(
        source = case_when(
            str_detect(dataset_base, "^16S_") ~ "Res",
            dataset_base %in% c("GMrepo_genus","GMrepo_species") ~ "GM",
            dataset_base == "CRC" ~ "CRC",
            TRUE ~ "Other"
        ),
        dataset_key = case_when(source == "GM" ~ "GMrepo", TRUE ~ dataset_base),
        study_label = if_else(source == "Res", paste0(study, " [", str_remove(dataset_key, "^16S_"), "]"), study),
        group = map_chr(group_raw, harmonize),
        group = if_else(group %in% names(global_pal), group, "Other"),
        mean = ifelse(is.na(mean) & !is.na(median), median, mean),
        sd   = ifelse(is.na(sd) & !is.na(q25) & !is.na(q75), (q75 - q25)/1.349, sd)
    ) |> filter(source %in% c("CRC","GM","Res"), !is.na(n), n > 1, !is.na(mean), !is.na(sd), sd >= 0)

# =========================
# 3. 计算箱线图统计量并对接 1B
# =========================
pool_mean_sd <- function(n, m, s){
    N <- sum(n); mean <- sum(n*m)/N; var <- (sum((n-1)*(s^2)) + sum(n*(m-mean)^2))/(N-1)
    list(n=N, mean=mean, sd=sqrt(var))
}

box_df <- all_long |>
    group_by(source, study_label, group_raw, group) |>
    slice_max(n, n=1, with_ties=FALSE) |>
    ungroup() |>
    group_by(source, study_label, group) |>
    summarise(p = list(pool_mean_sd(n, mean, sd)), n=p[[1]]$n, mean=p[[1]]$mean, sd=p[[1]]$sd, .groups="drop") |>
    rowwise() |>
    mutate(
        q25 = max(mean - 0.67448975 * sd, 0),
        upper = max(mean + 0.67448975 * sd, 0),
        iqr = upper - q25,
        ymin = max(q25 - 1.5 * iqr, 0),
        ymax = max(upper + 1.5 * iqr, 0)
    ) |> ungroup()

box_df <- box_df |>
    inner_join(proj_list |> select(source, study_label, order), by=c("source", "study_label"))

# =========================
# 4. 计算星号显著性
# =========================
welch_t_p <- function(m1, s1, n1, m2, s2, n2){
    if (any(is.na(c(m1,s1,n1,m2,s2,n2))) || n1<=1 || n2<=1 || s1<0 || s2<0) return(NA_real_)
    se <- sqrt(s1^2/n1 + s2^2/n2)
    if (se == 0) return(NA_real_)
    t  <- (m1 - m2) / se
    df <- (s1^2/n1 + s2^2/n2)^2 / ((s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1))
    2 * pt(-abs(t), df=df)
}

calc_study_star <- function(df_one){
    ctrl <- df_one |> filter(group=="Control") |> slice_max(order_by=n, n=1, with_ties=FALSE)
    if (nrow(ctrl)==0) return(tibble(p=NA_real_, star=""))
    others <- df_one |> filter(group!="Control", !is.na(mean), !is.na(sd), !is.na(n), n > 1, sd >= 0)
    if (nrow(others)==0) return(tibble(p=NA_real_, star=""))
    
    pvals <- mapply(
        FUN = function(m2,s2,n2) welch_t_p(ctrl$mean[1], ctrl$sd[1], ctrl$n[1], m2, s2, n2),
        m2 = others$mean, s2 = others$sd, n2 = others$n
    )
    pvals <- pvals[is.finite(pvals)]
    if (length(pvals)==0) return(tibble(p=NA_real_, star=""))
    pminv <- min(pvals, na.rm=TRUE)
    star_val <- if(pminv < 0.001) "***" else if(pminv < 0.01) "**" else if(pminv < 0.05) "*" else ""
    tibble(p=pminv, star=star_val)
}

p_tbl <- box_df |>
    group_by(source, study_label) |>
    group_split() |>
    map_dfr(function(df_one){
        tibble(source=df_one$source[1], study_label=df_one$study_label[1]) |> bind_cols(calc_study_star(df_one))
    })

ypos_tbl <- box_df |>
    group_by(source, study_label) |>
    summarise(y_star = max(ymax, na.rm=TRUE) + 1.0, .groups="drop") |>
    left_join(p_tbl, by=c("source","study_label")) |>
    filter(star != "")

# =========================
# 5. 画图 
# =========================
plot_panel_1C <- function(src, title_text){
    
    src_levels <- proj_list |> filter(source == src) |> arrange(order) |> pull(study_label)
    
    df <- box_df |> filter(source==src) |> mutate(study_label = factor(study_label, levels = src_levels))
    st <- ypos_tbl |> filter(source==src) |> mutate(study_label = factor(study_label, levels = src_levels))
    
    legend_labels <- function(x){
        abbr_map <- c(
            "Crohn disease" = "Crohn",
            "Ulcerative colitis" = "UC",
            "Multiple sclerosis" = "MS",
            "Type 1 diabetes" = "T1D",
            "Parkinson's disease" = "Parkinson",
            "Metabolic syndrome" = "MetS",
            "Allergic rhinitis" = "AR",
            "Interstitial lung disease" = "ILD",
            "RSV infection" = "RSV"
        )
        
        ifelse(x %in% names(abbr_map), abbr_map[x], x)
    }
    
    ggplot(df, aes(x=study_label, fill=group)) +
        geom_errorbar(aes(ymin=ymin, ymax=ymax), position=position_dodge(width=DODGE_W), width=CAP_W, linewidth=WHISKER_LW, color="grey25") +
        geom_boxplot(aes(ymin=q25, lower=q25, middle=mean, upper=upper, ymax=upper), stat="identity", position=position_dodge(width=DODGE_W), width=BOX_W, linewidth=BOX_OUTLINE_LW, colour="grey30", staplewidth=0, outlier.shape=NA) +
        geom_crossbar(aes(y=mean, ymin=mean, ymax=mean), stat="identity", position=position_dodge(width=DODGE_W), width=BOX_W*0.85, linewidth=MEDIAN_LW, colour="grey30") +
        geom_text(data=st, aes(x=study_label, y=y_star, label=star), inherit.aes=FALSE, size=4.3, fontface="plain", color="black") +
        scale_fill_manual(values=global_pal, labels=legend_labels) +
        coord_flip(clip="off") +
        labs(title=NULL, x=NULL, y="Age (years)", fill=NULL) +
        theme_classic(base_size = 13) +
        theme(
            plot.title = element_blank(),
            axis.text.y = element_text(size=10, face="plain", color="black"),
            axis.text.x = element_text(size=11, face="plain", color="black"),
            axis.title.y = element_text(size=9.5, face="plain", color="black"),
            legend.position = "right",
            legend.justification = "bottom",
            legend.direction = "vertical",
            legend.box = "vertical",
            legend.key.size = grid::unit(4.8, "mm"),
            plot.margin=ggplot2::margin(6,5,3,3, unit="mm")
        ) +
        guides(fill = guide_legend(ncol=1, byrow=TRUE))
}

p_crc_C <- plot_panel_1C("CRC", "CRC cohorts")
p_gm_C  <- plot_panel_1C("GM", "GMrepo metagenome")
p_res_C <- plot_panel_1C("Res", "ResMicroDB 16S (airway)")

if(!requireNamespace("patchwork", quietly=TRUE)) install.packages("patchwork")
library(patchwork)

# =======================================================
# 1C 导出部分 (【修改3】：同样压缩图注大小)
# =======================================================
p_all_C <- (p_crc_C | p_gm_C | p_res_C) + 
    patchwork::plot_layout(guides="keep") & 
    custom_theme & 
    theme(
        legend.position = "right", 
        legend.justification = "bottom",
        legend.title = element_blank(),                                    # 去掉图注标题，省空间
        legend.key.size = unit(3.5, "mm"),                                   # 图注色块略增大
        legend.text = element_text(size = 8.2, color = "black", face = "plain"),           # 把图注字号压缩，并加黑
        legend.spacing.y = unit(0.5, "mm"),                                # 缩小行距
        legend.margin = ggplot2::margin(0, 0, 0, 0, "mm")                           # 去掉图注多余边距
    )

ggsave(file.path(out_dir, "Figure1C_box_3panels_AUROC_Filtered.pdf"), p_all_C, width=12, height=3.2, scale=0.8)
ggsave(file.path(out_dir, "Figure1C_box_3panels_AUROC_Filtered.png"), p_all_C, width=12, height=3.2, dpi=600)

message("Success. Everything saved.")