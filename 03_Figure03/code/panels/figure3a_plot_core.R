library(tibble)
library(RColorBrewer)
library(dynutils)
library(stringr)
library(Hmisc)
library(plyr)
library(data.table)
library(dplyr)
library(tidyr)
library(scales)
library(pheatmap)
library(ggimage)
library(cowplot)

####------color-----------

method_palette = c(
  "velocyto" = "#D62728FF",
  "scVelo(sto)" = '#cab2d6',"scVelo(dyn)" = '#6a3d9a',
  "VeloAE" = '#33a02c' ,
  "dynamo" = '#FFB547FF',
  "VeloVAE" = '#b2df8a',
  "UniTVelo" = '#bcbd22',
  "DeepVelo(VAE)" = '#1f78b4',
  "cellDancer" = '#e377c2',
  "veloVI" = '#ff9896',
  "LatentVelo" = "#17BECFFF",
  "DeepVelo(GCN)" = '#a6cee3',
  "STT" = '#ff7f00',
  "MultiVelo"  ="#843C39FF",
  'PhyloVelo'="#E7CB94FF",
  'TFvelo' ="#7B4173FF",
  # Additional methods and colors
  "biVI" = "#2CA02CFF",
  # "velocycle" = "#E377C2FF",   # Deep pink
  "DeepKINET" = "#9467BDFF",   # Violet blue
  "TopicVelo" = "#1F77B4FF" ,   # Deep blue
  "SDEvelo" = "#8C564BFF",     # Brown
  "cell2fate" = "#F7DC6FFF",   # Bright yellow
  "TopoVelo" =   "#E377C2FF"
)

palettes = list('Metrics' = 'Blues',
                'Score' = 'Blues')

dataset_type_labels <- c(
  "bifurcating_converging" = "bifurcating converging",
  "bifurcating_cycle" = "bifurcating cycle",
  "bifurcating_loop" = "bifurcating loop",
  "bifurcating" = "bifurcating",
  "binary_tree" = "binary tree",
  "branching" = "branching",
  "consecutive_bifurcating" = "consecutive bifurcating",
  "converging" = "converging",
  "cycle" = "cycle",
  "cycle_simple" = "cycle (simple)",
  "disconnected" = "disconnected",
  "linear" = "linear",
  "linear_simple" = "linear (simple)",
  "trifurcating" = "trifurcating"
)

dataset_type_palette <- c(
  "bifurcating_converging" = "#E95C68",
  "bifurcating_cycle" = "#D7862B",
  "bifurcating_loop" = "#B89D27",
  "bifurcating" = "#84AC2F",
  "binary_tree" = "#51B848",
  "branching" = "#56C26A",
  "consecutive_bifurcating" = "#57C49A",
  "converging" = "#58C4C7",
  "cycle" = "#58B3E4",
  "cycle_simple" = "#5B92DA",
  "disconnected" = "#8E84CC",
  "linear" = "#B46DCC",
  "linear_simple" = "#D45DBC",
  "trifurcating" = "#EF4BAE"
)

dataset_seed_palette <- c(
  "seed1" = "#E2DEE9",
  "seed2" = "#B9B3D9",
  "seed3" = "#6D63B4"
)

####-----function---------------------
# add defaults for optional values
add_column_if_missing <- function(df, ...) {
  column_values <- list(...)
  for (column_name in names(column_values)) {
    default_val <- rep(column_values[[column_name]], nrow(df))

    if (column_name %in% colnames(df)) {
      df[[column_name]] <- ifelse(is.na(df[[column_name]]), default_val, df[[column_name]])
    } else {
      df[[column_name]] <- default_val
    }
  }
  df
}

check_method_name <- function(metrics){
  metrics <- metrics %>%
    filter(Method != "Dynamo_deterministic") %>%
    mutate(Method = gsub("Velocyto", "velocyto", Method,ignore.case = T))  %>%
    mutate(Method = gsub("scvelo_dynamical" , "scVelo(dyn)", Method,ignore.case = T))  %>%
    mutate(Method = gsub("scvelo_stochastic|scvelo$" , "scVelo(sto)", Method,ignore.case = T))  %>%
    mutate(Method = gsub("veloAE" , "VeloAE", Method,ignore.case = T))  %>%
    # mutate(Method = gsub("Dynamo_deterministic", "dynamo", Method))  %>%
    mutate(Method = gsub("Dynamo_stochastic|Dynamo", "dynamo", Method,ignore.case = T))  %>%
    mutate(Method = gsub("veloVAE", "VeloVAE", Method,ignore.case = T))  %>%
    mutate(Method = gsub("UniTvelo|UniteVelo", "UniTVelo", Method,ignore.case = T))  %>%
    mutate(Method = gsub("DeepVelo_SA|DeepveloSA", "DeepVelo(VAE)", Method,ignore.case = T)) %>%
    mutate(Method = gsub("Celldancer", "cellDancer", Method,ignore.case = T)) %>%
    mutate(Method = gsub("Latentvelo", "LatentVelo", Method,ignore.case = T)) %>%
    mutate(Method = gsub("DeepVelo_GB|DeepveloGB", "DeepVelo(GCN)", Method,ignore.case = T))%>%
    mutate(Method = gsub("stt", "STT", Method,ignore.case = T)) %>%
    mutate(Method = gsub("multivelo", "MultiVelo", Method,ignore.case = T))%>%
    mutate(Method = gsub("phylovelo", "PhyloVelo", Method,ignore.case = T)) %>%
    mutate(Method = gsub("tfvelo", "TFvelo", Method,ignore.case = T))

  metrics$Method <- factor(metrics$Method,levels = names(method_palette))
  metrics <- metrics[order(metrics$Method),]

  return(metrics)
}

normalize_by_col <- function(x){
  # ignore na
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)

  # ignore 0
  if (x_max == x_min) {
    return(rep(0, length(x)))
  } else {
    return ((x - x_min) / (x_max - x_min))
  }
}

parse_dataset_annotations <- function(dataset_ids) {
  dataset_type <- sub("_seed[0-9]+$", "", dataset_ids)
  dataset_seed <- sub("^.*_(seed[0-9]+)$", "\\1", dataset_ids)

  data.frame(
    dataset = dataset_ids,
    dataset_type = dataset_type,
    dataset_type_label = unname(dataset_type_labels[dataset_type]),
    dataset_seed = dataset_seed,
    stringsAsFactors = FALSE
  )
}


####---------plot_summmary-------------------
plot_summary = function(datapath,figpath,title_text = NULL,panel_label = NULL){

  ####------read data-----------------------
  metrics <- fread(datapath)

  ####--------process data------------------------
  colnames(metrics)[1] <- 'Method'

  metrics <- metrics %>%
    mutate(Method= gsub("^[0-9]+_", "", Method))

  metrics <- check_method_name(metrics)

  metrics$Method <- factor(metrics$Method,levels = names(method_palette))

  metrics <- metrics[order(metrics$Method)]

  # add median & mean
  metrics$Median <- apply(metrics[,2:43], 1, median, na.rm = TRUE)
  metrics$Mean <- apply(metrics[,2:43], 1, mean, na.rm = TRUE)

  metrics <- as.data.frame(metrics)
  # normalize
  # metrics[,2:43] <- apply(metrics[,2:43],2,normalize_by_col)
  # metrics[,44:45] <- apply(metrics[,44:45],2,normalize_by_col)

  # normalize by group
  # x_min <- -1
  # x_max <-1
  x_min <- min(metrics[,2:43], na.rm = TRUE)-0.01
  x_max <- max(metrics[,2:43], na.rm = TRUE)
  metrics[,2:43] <- (metrics[,2:43] - x_min) / (x_max - x_min)

  # mean and median
  x_min <- min(metrics[,44:45], na.rm = TRUE)-0.01
  x_max <- max(metrics[,44:45], na.rm = TRUE)
  metrics[,44:45] <- (metrics[,44:45] - x_min) / (x_max - x_min)

  # # drop NA rows
  # metrics %>%
  #   filter(!if_all(.fns = is.na))

  # x_min <- -1
  # x_max <- 1
  # metrics[,2:45] <- (metrics[,2:45] - x_min) / (x_max - x_min)

  ####----- basic info ------------
  row_height <- 1.1
  row_space <- .1
  row_bigspace <- .5
  col_width <- 1.1
  col_space <- .2
  col_bigspace <- .5
  segment_data <- NULL

  data =as.data.frame(metrics)
  dataset_cols <- colnames(metrics)[2:43]
  dataset_annotations <- parse_dataset_annotations(dataset_cols)

  ####----- row info and pos--------
  row_info <- as.data.frame(metrics$Method)
  # colnames(row_info) <- 'Method'
  colnames(row_info) <- 'id'
  row_pos <-
    row_info %>%
    mutate(
      group_i = row_number(),
      row_i = row_number(),
      colour_background = group_i %% 2 == 1,
      # do_spacing = c(FALSE, diff(as.integer(factor(group))) != 0),
      do_spacing = FALSE,
      ysep = ifelse(do_spacing, row_height + 2 * row_space, row_space),
      y = - (row_i * row_height + cumsum(ysep)),
      ymin = y - row_height / 2,
      ymax = y + row_height / 2
    )

  ####----- col info and pos--------
  column_info <- data.frame(id = colnames(metrics))

  column_info$group = "Metrics"
  column_info[column_info$id =='Method','group'] = 'Text'
  column_info[column_info$id %in% c('Median','Mean'),'group'] = 'Score'


  column_info$geom <- 'circle'
  column_info[column_info$id =='Method','geom'] = 'text'
  column_info[column_info$id %in% c('Median','Mean'),'geom'] = 'bar'

  column_info$width <- 1
  column_info[column_info$id =='Method','width'] = 6
  column_info[column_info$id %in% c('Median','Mean'),'width'] = 3

  column_info$overlay <- FALSE

  column_pos <-
    column_info %>%
    mutate(
      do_spacing = c(FALSE, diff(as.integer(factor(group))) != 0),
      xsep = case_when(
        overlay ~ c(0, -head(width, -1)),
        do_spacing ~ col_bigspace,
        TRUE ~ col_space
      ),
      xwidth = case_when(
        overlay & width < 0 ~ width - xsep,
        overlay ~ -xsep,
        TRUE ~ width
      ),
      xmax = cumsum(xwidth + xsep),
      xmin = xmax - xwidth,
      x = xmin + xwidth / 2
    )

  dataset_annotations$xmin <- column_pos$xmin[match(dataset_annotations$dataset, column_info$id)]
  dataset_annotations$xmax <- column_pos$xmax[match(dataset_annotations$dataset, column_info$id)]
  dataset_annotations$x <- column_pos$x[match(dataset_annotations$dataset, column_info$id)]

  ####-----circle data------------
  # gather circle data
  ind_circle <- which(column_info$geom == "circle")

  dat_mat <- as.matrix(data[, ind_circle])
  col_palette <- data.frame(metric = colnames(dat_mat),
                            group = column_info[match(colnames(dat_mat), column_info$id), "group"])

  # col_palette$name_palette <- 'RdPu'
  col_palette$name_palette <- lapply(col_palette$group, function(x) palettes[[as.character(x)]])

  circle_data <- data.frame(label = unlist(lapply(colnames(dat_mat),
                                                  function(x) rep(x, nrow(dat_mat)))),
                            x0 = unlist(lapply(column_pos$x[ind_circle],
                                               function(x) rep(x, nrow(dat_mat)))),
                            y0 = rep(row_pos$y, ncol(dat_mat)),
                            r = row_height/2*as.vector(sqrt(dat_mat))
  )

  # scale by col
  # for(l in unique(circle_data$label)){
  #   ind_l <- which(circle_data$label == l)
  #   circle_data[ind_l, "r"] <- rescale(circle_data[ind_l, "r"], to = c(0.05, 0.55), from = range(circle_data[ind_l, "r"], na.rm = T))
  # }

  colors <- NULL


  for(i in 1:ncol(dat_mat)){
    palette <- colorRampPalette(rev(brewer.pal(9, col_palette$name_palette[[i]])))(nrow(data)-sum(is.na(dat_mat[,i])))
    colors <- c(colors, palette[rank(dat_mat[,i], ties.method = "average", na.last = "keep")])
  }

  circle_data$colors <- colors

  ####-------bar/rect data-----------
  ind_bar <- which(column_info$geom == "bar")
  dat_mat <- as.matrix(data[, ind_bar])
  # colnames(dat_mat) <- 'Median'

  col_palette <- data.frame(metric = colnames(dat_mat),
                            group = column_info[match(colnames(dat_mat), column_info$id), "group"])


  col_palette$name_palette <- lapply(col_palette$group, function(x) palettes[[as.character(x)]])
  # col_palette$name_palette <- 'BuPu'

  rect_data <- data.frame(label = unlist(lapply(colnames(dat_mat),
                                                function(x) rep(x, nrow(dat_mat)))),
                          method = rep(row_info$id, ncol(dat_mat)),
                          value = as.vector(dat_mat),
                          xmin = unlist(lapply(column_pos[ind_bar, "xmin"],
                                               function(x) rep(x, nrow(dat_mat)))),
                          xmax = unlist(lapply(column_pos[ind_bar, "xmax"],
                                               function(x) rep(x, nrow(dat_mat)))),
                          ymin = rep(row_pos$ymin, ncol(dat_mat)),
                          ymax = rep(row_pos$ymax, ncol(dat_mat)),
                          xwidth = unlist(lapply(column_pos[ind_bar, "xwidth"],
                                                 function(x) rep(x, nrow(dat_mat))))
  )

  rect_data <- rect_data %>%
    add_column_if_missing(hjust = 0) %>%
    mutate(
      xmin = xmin + (1 - value) * xwidth * hjust,
      xmax = xmax - (1 - value) * xwidth * (1 - hjust)
    )

  colors <- NULL

  for(i in 1:ncol(dat_mat)){
    palette <- colorRampPalette(rev(brewer.pal(9, col_palette$name_palette[[i]])))(nrow(data)-sum(is.na(dat_mat[,i])))
    colors <- c(colors, palette[rank(dat_mat[,i], ties.method = "average", na.last = "keep")])
  }

  rect_data$colors <- colors

  ####------top annotations-------
  type_ymin <- 2.8
  type_ymax <- 3.7
  seed_ymin <- 1.95
  seed_ymax <- 2.65

  dataset_type_blocks <- dataset_annotations %>%
    mutate(block = cumsum(c(TRUE, diff(match(dataset_type, names(dataset_type_labels))) != 0))) %>%
    group_by(block, dataset_type, dataset_type_label) %>%
    summarise(
      xmin = min(xmin),
      xmax = max(xmax),
      .groups = "drop"
    )

  top_rect_data <- bind_rows(
    dataset_type_blocks %>%
      transmute(
        label = "dataset_type",
        method = dataset_type,
        value = 1,
        xmin = xmin,
        xmax = xmax,
        ymin = type_ymin,
        ymax = type_ymax,
        xwidth = xmax - xmin,
        colors = unname(dataset_type_palette[dataset_type]),
        alpha = 1,
        border = TRUE,
        border_colour = "grey30"
      ),
    dataset_annotations %>%
      transmute(
        label = "dataset_seed",
        method = dataset_seed,
        value = 1,
        xmin = xmin,
        xmax = xmax,
        ymin = seed_ymin,
        ymax = seed_ymax,
        xwidth = xmax - xmin,
        colors = unname(dataset_seed_palette[dataset_seed]),
        alpha = 1,
        border = TRUE,
        border_colour = "grey30"
      )
  )

  rect_data <- bind_rows(rect_data, top_rect_data)

  ####-------text data------------
  ind_text <- which(column_info$geom == "text")
  dat_mat <- as.matrix(data[, ind_text])
  colnames(dat_mat) <- 'Method'

  text_data <- data.frame(label_value = as.vector(dat_mat),
                          group = rep(colnames(dat_mat), each = nrow(dat_mat)),
                          xmin = unlist(lapply(column_pos[ind_text, "xmin"],
                                               function(x) rep(x, nrow(dat_mat)))),
                          xmax = unlist(lapply(column_pos[ind_text, "xmax"],
                                               function(x) rep(x, nrow(dat_mat)))),
                          ymin = rep(row_pos$ymin, ncol(dat_mat)),
                          ymax = rep(row_pos$ymax, ncol(dat_mat)),
                          size = 4, fontface = "plain", stringsAsFactors = F)

  #####-------- add top3 ranking for each bar column-------------
  # if(usability || atac_best){
  cols_bar <- unique(rect_data$label)
  cols_bar <- as.character(cols_bar[!is.na(cols_bar)])
  for(c in cols_bar){
    rect_tmp <- rect_data[rect_data$label == c,]
    rect_tmp <- add_column(rect_tmp, "label_value" = as.character(rank(-rect_tmp$value, ties.method = "min")))
    rect_tmp <- rect_tmp[rect_tmp$label_value %in% c("1", "2", "3","4","5","6"), c("label_value", "xmin", "xmax", "ymin", "ymax")]
    rect_tmp <- add_column(rect_tmp, "size" = 2.5, .after = "ymax")
    rect_tmp <- add_column(rect_tmp, "colors" = "black", .after = "size")
    rect_tmp <- add_column(rect_tmp, "fontface" = "plain", .after = "colors")
    rect_tmp <- add_column(rect_tmp, "group" = "top3", .after = "fontface")
    text_data <- bind_rows(text_data, rect_tmp)
  }
  # }


  ####-----add column names-------------
  # ADD COLUMN NAMES
  df <- column_pos %>% filter(id != "Method") %>% filter(id != "Ranking")

  if (nrow(df) > 0) {
    segment_data <- segment_data %>% bind_rows(
      df %>% transmute(x = x, xend = x, y = 0.05, yend = 0.22, size = .5)
    )
    text_data <-
      bind_rows(
        text_data,
        df %>% transmute(
          xmin = x, xmax = x, ymin = 0.28, ymax = 0.52,
          angle = 30, vjust = 0, hjust = 0,
          label_value = id,
          size = 2.4
        )
      )
  }

  ####-----title and dataset annotations------
  if (is.null(title_text)) {
    title_text <- "Cosine similarity of high-dim RNA velocity on cell level"
  }

  annotation_text <- bind_rows(
    data.frame(
      xmin = min(dataset_annotations$xmin) - 6.2,
      xmax = min(dataset_annotations$xmin) - 0.5,
      ymin = type_ymin,
      ymax = type_ymax,
      angle = 0,
      vjust = 0.5,
      hjust = 1,
      label_value = "Datasets type",
      size = 4,
      fontface = "plain",
      group = "annotation",
      stringsAsFactors = FALSE
    ),
    data.frame(
      xmin = min(dataset_annotations$xmin) - 6.2,
      xmax = min(dataset_annotations$xmin) - 0.5,
      ymin = seed_ymin,
      ymax = seed_ymax,
      angle = 0,
      vjust = 0.5,
      hjust = 1,
      label_value = "Datasets seed",
      size = 4,
      fontface = "plain",
      group = "annotation",
      stringsAsFactors = FALSE
    ),
    data.frame(
      xmin = min(dataset_annotations$xmin),
      xmax = max(dataset_annotations$xmax),
      ymin = type_ymax + 1.15,
      ymax = type_ymax + 1.95,
      angle = 0,
      vjust = 0.5,
      hjust = 0.5,
      label_value = title_text,
      size = 5.5,
      fontface = "plain",
      group = "annotation",
      stringsAsFactors = FALSE
    )
  )

  if (!is.null(panel_label)) {
    annotation_text <- bind_rows(
      annotation_text,
      data.frame(
        xmin = min(dataset_annotations$xmin) - 9.2,
        xmax = min(dataset_annotations$xmin) - 8.4,
        ymin = type_ymax + 1.05,
        ymax = type_ymax + 2.05,
        angle = 0,
        vjust = 0.5,
        hjust = 0,
        label_value = panel_label,
        size = 8,
        fontface = "bold",
        group = "annotation",
        stringsAsFactors = FALSE
      )
    )
  }

  text_data <- bind_rows(text_data, annotation_text)


  ####----------row annotation---------------
  # # GENERATE ROW ANNOTATION
  # if (plot_row_annotation) {
  #   row_annotation <-
  #     row_pos %>%
  #     select(group, ymin, ymax) %>%
  #     group_by(group) %>%
  #     summarise(
  #       ymin = min(ymin),
  #       ymax = max(ymax),
  #       y = (ymin + ymax) / 2
  #     ) %>%
  #     ungroup() %>%
  #     mutate(xmin = -.5, xmax = 5) %>%
  #     filter(!is.na(group), group != "")
  #
  #   text_data <- text_data %>% bind_rows(
  #     row_annotation %>%
  #       transmute(xmin, xmax, ymin = ymax + row_space, label_value = group %>% gsub("\n", " ", .),
  #                 hjust = 0, vjust = .5, fontface = "bold", size = 4) %>%
  #       mutate(ymax = ymin + row_height)
  #   )
  # }
  #
  # ####----------image data---------------
  # # gather image data
  # ind_img <- which(column_info$geom == "image")
  # if(length(ind_img) > 0){
  #   dat_mat <- as.matrix(data[, ind_img])
  #
  #   image_data <- data.frame(x = unlist(lapply(column_pos$x[ind_img],
  #                                              function(x) rep(x, nrow(dat_mat)))),
  #                            y = rep(row_pos$y, ncol(dat_mat)),
  #                            image = mapvalues(dat_mat, from = c("graph", "embed", "gene"),
  #                                              to = c("./img/graph.png", "./img/embedding.png", "./img/matrix.png")),
  #                            stringsAsFactors = FALSE
  #   )
  #
  # }
  #
  # suppressWarnings({
  #   minimum_x <- min(column_pos$xmin, segment_data$x, segment_data$xend,
  #                    text_data$xmin, na.rm = TRUE)
  #   maximum_x <- max(column_pos$xmax, segment_data$x, segment_data$xend,
  #                    text_data$xmax, na.rm = TRUE)
  #   minimum_y <- min(row_pos$ymin, segment_data$y, segment_data$yend,
  #                    text_data$ymin, na.rm = TRUE)
  #   maximum_y <- max(row_pos$ymax, segment_data$y, segment_data$yend,
  #                    text_data$ymax, na.rm = TRUE)
  # })

  ####---- legend---------
  # CREATE LEGEND for circle scores

  minimum_x <- min(column_pos$xmin, segment_data$x, segment_data$xend,
                   text_data$xmin, na.rm = TRUE)
  maximum_x <- max(column_pos$xmax, segment_data$x, segment_data$xend,
                   text_data$xmax, na.rm = TRUE)
  minimum_y <- min(row_pos$ymin, segment_data$y, segment_data$yend,
                   text_data$ymin, na.rm = TRUE)
  maximum_y <- max(row_pos$ymax, segment_data$y, segment_data$yend,
                   text_data$ymax, na.rm = TRUE)

  ####-----bottom mean bars-------
  column_means <- colMeans(data[, dataset_cols], na.rm = TRUE)
  mean_palette <- colorRampPalette(rev(brewer.pal(9, palettes[["Metrics"]])))(
    length(column_means) - sum(is.na(column_means))
  )

  mean_rect_data <- dataset_annotations %>%
    mutate(
      value = column_means[dataset],
      ymax = minimum_y - 0.35,
      ymin = ymax - value * 1.8,
      xwidth = xmax - xmin,
      colors = mean_palette[rank(value, ties.method = "average", na.last = "keep")],
      alpha = 1,
      border = TRUE,
      border_colour = "grey30",
      label = "column_mean",
      method = dataset
    ) %>%
    select(label, method, value, xmin, xmax, ymin, ymax, xwidth, colors, alpha, border, border_colour)

  rect_data <- bind_rows(rect_data, mean_rect_data)
  text_data <- bind_rows(
    text_data,
    data.frame(
      xmin = min(dataset_annotations$xmin) - 3.8,
      xmax = min(dataset_annotations$xmin) - 0.5,
      ymin = min(mean_rect_data$ymin),
      ymax = max(mean_rect_data$ymax),
      angle = 0,
      vjust = 0.5,
      hjust = 1,
      label_value = "Mean",
      size = 4,
      fontface = "plain",
      group = "annotation",
      stringsAsFactors = FALSE
    )
  )

  minimum_y <- min(minimum_y, min(mean_rect_data$ymin, na.rm = TRUE))

  x_min_output <- minimum_x + 8
  x_min_scaling <- minimum_x + 25
  x_min_ranking <- minimum_x + 34
  x_min_score <-  minimum_x + 43

  leg_max_y <- minimum_y - 1.2

  # # Create legend for Output
  # leg_min_x <- x_min_output
  # output_title_data <- data.frame(xmin = leg_min_x,
  #                                 xmax = leg_min_x+ 2,
  #                                 ymin = leg_max_y - 1,
  #                                 ymax = leg_max_y,
  #                                 label_value = "Output",
  #                                 hjust = 0, vjust = 0,
  #                                 fontface = "bold",
  #                                 size = 3)
  #
  # output_img <- data.frame(x = leg_min_x+0.5,
  #                          y = c(leg_max_y-2, leg_max_y-3.2,leg_max_y-4.4),
  #                          image = c("./img/matrix.png", "./img/embedding.png", "./img/graph.png")
  # )
  # if(atac || atac_best){
  #   output_text <- data.frame(xmin = leg_min_x+1.5,
  #                             xmax = leg_min_x+3,
  #                             ymin = c(leg_max_y-2.2, leg_max_y-3.4,leg_max_y-4.6),
  #                             ymax = c(leg_max_y-2.2, leg_max_y-3.4,leg_max_y-4.6),
  #                             label_value = c("feature", "embed", "graph"),
  #                             hjust = 0, vjust = 0,
  #                             fontface = "plain",
  #                             size = 3)
  # } else{
  #   output_text <- data.frame(xmin = leg_min_x+1.5,
  #                             xmax = leg_min_x+3,
  #                             ymin = c(leg_max_y-2.2, leg_max_y-3.4,leg_max_y-4.6),
  #                             ymax = c(leg_max_y-2.2, leg_max_y-3.4,leg_max_y-4.6),
  #                             label_value = c("gene", "embed", "graph"),
  #                             hjust = 0, vjust = 0,
  #                             fontface = "plain",
  #                             size = 3)
  # }
  #
  # text_data <- bind_rows(text_data, output_text, output_title_data)
  # image_data <- bind_rows(image_data, output_img)
  #
  # # Create legend for scaling
  # if(!atac && !atac_best){
  #   leg_min_x <- x_min_scaling
  #   scaling_title_data <- data.frame(xmin = leg_min_x,
  #                                    xmax = leg_min_x+ 2,
  #                                    ymin = leg_max_y - 1,
  #                                    ymax = leg_max_y,
  #                                    label_value = "Scaling",
  #                                    hjust = 0, vjust = 0,
  #                                    fontface = "bold",
  #                                    size = 3)
  #
  #   scaling_text <- data.frame(xmin = c(leg_min_x, leg_min_x+1),
  #                              xmax = c(leg_min_x+0.5, leg_min_x+3),
  #                              ymin = c(rep(leg_max_y-2,2), rep(leg_max_y-3,2)),
  #                              ymax = c(rep(leg_max_y-1,2), rep(leg_max_y-2,2)),
  #                              label_value = c("+", ": scaled", "-", ": unscaled"),
  #                              hjust = 0, vjust = 0,
  #                              fontface = c("bold","plain", "bold", "plain"),
  #                              size = c(5,3,5,3))
  #
  #   text_data <- bind_rows(text_data, scaling_title_data, scaling_text)
  # }

  ####-----dataset type legend------
  type_keys <- names(dataset_type_labels)
  type_labels <- unname(dataset_type_labels)
  type_legend <- data.frame(
    type_key = type_keys,
    type_label = type_labels,
    col_id = rep(c(0, 1), each = ceiling(length(type_keys) / 2))[seq_along(type_keys)],
    row_id = rep(seq_len(ceiling(length(type_keys) / 2)), 2)[seq_along(type_keys)],
    stringsAsFactors = FALSE
  )

  type_legend$xmin <- x_min_output + 2 + type_legend$col_id * 10
  type_legend$xmax <- type_legend$xmin + 0.8
  type_legend$ymax <- leg_max_y - 1.4 - (type_legend$row_id - 1) * 0.9
  type_legend$ymin <- type_legend$ymax - 0.8

  rect_data <- bind_rows(
    rect_data,
    type_legend %>%
      transmute(
        label = "legend_type",
        method = type_key,
        value = 1,
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        xwidth = xmax - xmin,
        colors = unname(dataset_type_palette[type_key]),
        alpha = 1,
        border = TRUE,
        border_colour = "grey30"
      )
  )

  text_data <- bind_rows(
    text_data,
    data.frame(
      xmin = x_min_output,
      xmax = x_min_output + 1,
      ymin = mean(range(type_legend$ymin, type_legend$ymax)) - 0.4,
      ymax = mean(range(type_legend$ymin, type_legend$ymax)) + 0.4,
      angle = 90,
      vjust = 0.5,
      hjust = 0.5,
      label_value = "Datasets type",
      size = 3.5,
      fontface = "bold",
      group = "annotation",
      stringsAsFactors = FALSE
    ),
    type_legend %>%
      transmute(
        xmin = xmax + 0.25,
        xmax = xmax + 6.8,
        ymin = ymin,
        ymax = ymax,
        angle = 0,
        vjust = 0.5,
        hjust = 0,
        label_value = type_label,
        size = 3.2,
        fontface = "plain",
        group = "annotation"
      )
  )

  ####-----dataset seed legend------
  seed_keys <- names(dataset_seed_palette)
  seed_legend <- data.frame(
    seed_key = seed_keys,
    row_id = seq_along(seed_keys),
    stringsAsFactors = FALSE
  )

  seed_legend$xmin <- x_min_scaling + 2
  seed_legend$xmax <- seed_legend$xmin + 0.8
  seed_legend$ymax <- leg_max_y - 1.4 - (seed_legend$row_id - 1) * 0.9
  seed_legend$ymin <- seed_legend$ymax - 0.8

  rect_data <- bind_rows(
    rect_data,
    seed_legend %>%
      transmute(
        label = "legend_seed",
        method = seed_key,
        value = 1,
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        xwidth = xmax - xmin,
        colors = unname(dataset_seed_palette[seed_key]),
        alpha = 1,
        border = TRUE,
        border_colour = "grey30"
      )
  )

  text_data <- bind_rows(
    text_data,
    data.frame(
      xmin = x_min_scaling,
      xmax = x_min_scaling + 1,
      ymin = mean(range(seed_legend$ymin, seed_legend$ymax)) - 0.3,
      ymax = mean(range(seed_legend$ymin, seed_legend$ymax)) + 0.3,
      angle = 90,
      vjust = 0.5,
      hjust = 0.5,
      label_value = "Datasets seed",
      size = 3.5,
      fontface = "bold",
      group = "annotation",
      stringsAsFactors = FALSE
    ),
    seed_legend %>%
      transmute(
        xmin = xmax + 0.25,
        xmax = xmax + 3.5,
        ymin = ymin,
        ymax = ymax,
        angle = 0,
        vjust = 0.5,
        hjust = 0,
        label_value = seed_key,
        size = 3.2,
        fontface = "plain",
        group = "annotation"
      )
  )

  # CREATE LEGEND for ranking colors

  ####-----rank legend-----------
  # CREATE LEGEND for ranking colors
  leg_min_x <- x_min_ranking
  rank_groups <- as.character(column_info[column_info$geom == "bar", "group"])
  # rank_groups <- as.character(column_info[column_info$geom == "circle", "group"])
  # rank_groups <- as.character(column_info[column_info$geom %in% c("circle",'bar'), "group"])

  rank_minimum_x <- list("Metrics" = leg_min_x,
                         "Score" = leg_min_x+1)
  leg_max_x <- leg_min_x+2


  rank_title_data <- data.frame(xmin = leg_min_x,
                                xmax = leg_min_x+ 2,
                                ymin = leg_max_y - 1,
                                ymax = leg_max_y,
                                label_value = "Ranking",
                                hjust = 0, vjust = 0,
                                fontface = "bold")


  for(rg in rank_groups){
    rank_palette <- colorRampPalette(rev(brewer.pal(9, palettes[[rg]])))(5)

    rank_data <- data.frame(xmin = rank_minimum_x[[rg]],
                            xmax = rank_minimum_x[[rg]] + .8,
                            ymin = seq(leg_max_y-4, leg_max_y - 2, by = .5),
                            ymax = seq(leg_max_y-3.5, leg_max_y -1.5, by = .5),
                            border = TRUE,
                            colors = rank_palette
    )
    rect_data <- bind_rows(rect_data, rank_data)

  }

  # create arrow for ranking
  arrow_data <- data.frame(x = leg_max_x + 1.5,
                           xend = leg_max_x +1.5,
                           y = leg_max_y-4,
                           yend = leg_max_y -1.5)


  # add text next to the arrow
  arrow_text <- data.frame(xmin = leg_max_x +2,
                           xmax = leg_max_x +2.5,
                           ymin = c(leg_max_y-2, leg_max_y-4),
                           ymax = c(leg_max_y-1.5, leg_max_y-3.5 ),
                           label_value = c("1", as.character(nrow(data))),
                           hjust = 0, vjust = 0, size = 2.5)


  text_data <- bind_rows(text_data, rank_title_data, arrow_text)

  ####---------circle legend--------
  cir_minimum_x <- x_min_score

  cir_legend_size <- 1
  cir_legend_space <- .1

  cir_legend_dat <-
    data.frame(
      value = seq(0, 1, by = .2),
      r = row_height/2*seq(0, 1, by = .2)
    )

  # cir_legend_dat$r <- rescale(cir_legend_dat$r, to = c(0.05, 0.55), from = range(cir_legend_dat$r, na.rm = T))

  x0 <- vector("integer", nrow(cir_legend_dat))
  for(i in 1:length(x0)){
    if(i == 1){
      x0[i] <- cir_minimum_x + cir_legend_space + cir_legend_dat$r[i]
    }
    else {
      x0[i] <- x0[i-1] + cir_legend_dat$r[i-1] + cir_legend_space + cir_legend_dat$r[i]
    }
  }

  cir_legend_dat$x0 <- x0
  cir_legend_min_y <- leg_max_y-4
  cir_legend_dat$y0 <- cir_legend_min_y + 1 + cir_legend_dat$r

  cir_legend_dat$colors <- NULL
  cir_maximum_x <- max(cir_legend_dat$x0)

  cir_title_data <- data_frame(xmin = cir_minimum_x,
                               xmax = cir_maximum_x,
                               ymin = leg_max_y -1,
                               ymax = leg_max_y,
                               label_value = "Score",
                               hjust = 0, vjust = 0, fontface = "bold")

  cir_value_data <- data.frame(xmin = cir_legend_dat$x0 - cir_legend_dat$r,
                               xmax = cir_legend_dat$x0 + cir_legend_dat$r,
                               ymin = cir_legend_min_y,
                               ymax = cir_legend_min_y +3,
                               hjust = .5, vjust = 0, size = 2.5,
                               label_value = ifelse(cir_legend_dat$value %in% c(0, 1),
                                                    paste0(cir_legend_dat$value*100, "%"), ""))

  circle_data <- bind_rows(circle_data, cir_legend_dat)
  text_data <- bind_rows(text_data, cir_title_data, cir_value_data)

  legend_box <- data.frame(
    label = "legend_box",
    method = "legend_box",
    value = 1,
    xmin = x_min_output - 1.3,
    xmax = cir_maximum_x + 4.8,
    ymin = leg_max_y - 7.2,
    ymax = leg_max_y + 0.6,
    xwidth = (cir_maximum_x + 4.8) - (x_min_output - 1.3),
    colors = NA,
    alpha = 1,
    border = TRUE,
    border_colour = "grey50",
    stringsAsFactors = FALSE
  )

  rect_data <- bind_rows(rect_data, legend_box)

  minimum_y <- min(minimum_y, min(text_data$ymin, na.rm = TRUE))

  ####----GGPLOT ------------
  g <-
    ggplot() +
    coord_equal(expand = FALSE) +
    scale_alpha_identity() +
    scale_colour_identity() +
    scale_fill_identity() +
    scale_size_identity() +
    scale_linetype_identity() +
    cowplot::theme_nothing()

  # PLOT ROW BACKGROUNDS
  df <- row_pos %>% filter(colour_background)

  if (nrow(df) > 0) {
    g <- g +
      geom_rect(aes(xmin = min(column_pos$xmin)-.25,
                    xmax = max(column_pos$xmax)+.25,
                    ymin = ymin - (row_space / 2),
                    ymax = ymax + (row_space / 2)),
                df, fill = "#DDDDDD")
  }
  g

  ####------ PLOT CIRCLES ------------
  if (length(ind_circle) > 0) {
    g <- g +
      ggforce::geom_circle(
        aes(x0 = x0, y0 = y0, fill= colors, r = r),
        circle_data,
        size=.25)
  }
  g

  ####--- PLOT RECTANGLES ------------
  if (nrow(rect_data) > 0) {
    # add defaults for optional values
    rect_data <- rect_data %>%
      add_column_if_missing(alpha = 1, border = TRUE, border_colour = "black") %>%
      mutate(border_colour = ifelse(border, border_colour, NA))

    g <- g +
      geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin,ymax = ymax,
                    fill = colors, colour = border_colour, alpha = alpha),
                rect_data,
                size = .25)
  }


  ####------ PLOT TEXT------------------

  text_data <- text_data %>%
    add_column_if_missing(
      hjust = .5,
      vjust = .5,
      size = 3,
      fontface = "plain",
      colors = "black",
      lineheight = 1,
      angle = 0
    ) %>%
    mutate(
      angle2 = angle / 360 * 2 * pi,
      cosa = cos(angle2) %>% round(2),
      sina = sin(angle2) %>% round(2),
      alphax = ifelse(cosa < 0, 1 - hjust, hjust) * abs(cosa) + ifelse(sina > 0, 1 - vjust, vjust) * abs(sina),
      alphay = ifelse(sina < 0, 1 - hjust, hjust) * abs(sina) + ifelse(cosa < 0, 1 - vjust, vjust) * abs(cosa),
      x = (1 - alphax) * xmin + alphax * xmax,
      y = (1 - alphay) * ymin + alphay * ymax
    ) %>%
    filter(label_value != "")

  # Set fontface for legend bold
  text_data[text_data$label_value == "Ranking", "fontface"] <- "bold"

  # Set fontface for ranking numbers bold
  # if(usability || atac_best){
  #   text_data[1:nrow(data), "fontface"] <- "bold"
  # }
  # subset text_data to left-aligned rows

  text_data_left <- text_data[which(text_data$group == "Method" | text_data$group == "top3"), ]
  text_data <- text_data[-which(text_data$group == "Method" | text_data$group == "top3"), ]

  g <- g + geom_text(aes(x = x, y = y, label = label_value,
                         colour = colors, hjust = hjust, vjust = vjust,
                         size = size, fontface = fontface, angle = angle),
                     data = text_data)
  g

  text_data_left[text_data_left$group == "Method", "x"] <- text_data_left[text_data_left$group == "Method", "x"] - 3

  # if(usability || atac_best){
  #   text_data_left[text_data_left$group == "top3", "x"] <- text_data_left[text_data_left$group == "top3", "xmin"] + .3
  #   text_data_left[text_data_left$group == "Method", "x"] <- text_data_left[text_data_left$group == "Method", "x"] + .5
  # }

  g <- g + geom_text(aes(x = x, y = y,
                         label = label_value, colour = colors, hjust = "left",
                         vjust = vjust, size = size, fontface = fontface, angle = angle),
                     data = text_data_left)
  g


  g <- g +
    geom_text(aes(x = x, y = y, label = label_value, colour = colors,
                  hjust = "left", vjust = vjust, size = size,
                  fontface = fontface, angle = angle),
              data = text_data_left)
  g


  #### -------- plot segments-----------------
  # PLOT SEGMENTS
  if (nrow(segment_data) > 0) {
    # add defaults for optional values
    segment_data <- segment_data %>% add_column_if_missing(size = .5, colour = "black", linetype = "solid")

    g <- g + geom_segment(aes(x = x, xend = xend, y = y, yend = yend, size = size, colour = colour, linetype = linetype), segment_data)
  }

  # PLOT ARROW RANKING
  if (nrow(arrow_data) > 0) {
    # add defaults for optional values
    arrow_data <- arrow_data %>% add_column_if_missing(size = .5, colour = "black", linetype = "solid")

    g <- g + geom_segment(aes(x = x, xend = xend, y = y, yend = yend, size = size, colour = colour, linetype = linetype), arrow_data, arrow = arrow(length = unit(0.1, "cm")), lineend = "round", linejoin = "bevel")
  }

  # # PLOT IMAGES
  # if(length(ind_img) > 0){
  #   for(r in 1:nrow(image_data)){
  #     g <- g + cowplot::draw_image(image = image_data$image[r], x = image_data[r, "x"]-.5, y = image_data[r, "y"]-.5)
  #   }
  #
  # }


  ####------  add  size------------
  # reserve a bit more room for text that wants to go outside the frame

  minimum_x <- minimum_x - 2
  maximum_x <- maximum_x + 5
  minimum_y <- minimum_y - 2
  maximum_y <- maximum_y + 4

  g$width <- maximum_x - minimum_x
  g$height <- maximum_y - minimum_y

  g <- g +
    expand_limits(x = c(minimum_x, maximum_x), y = c(minimum_y, maximum_y))

  ggsave(figpath,
         g,
         device = "pdf", width = 420, height = 300, units = "mm")

}
