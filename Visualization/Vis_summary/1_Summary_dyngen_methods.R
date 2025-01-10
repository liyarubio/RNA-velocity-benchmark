library(tibble)
library(RColorBrewer)
library(dynutils)
library(stringr)
library(Hmisc)
library(plyr)
library(data.table)
library(dplyr)
library(scales)
library(ggimage)
library(cowplot)

# source("knit_table.R") # Please put knit_table.R in your working dir
# setwd('/media/liyaru/LYR/Diff_change/2_benchmark/25_multiomic_benchmark/scib-reproducibility-main/visualization')

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
  'TFvelo' ="#7B4173FF"
)

palettes = list('Metrics' = 'Blues',
                'Score' = 'Blues')

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


####---------plot_summmary-------------------
plot_summary = function(datapath,outpath,figpath){
  
  ####------read data-----------------------
  # metrics <- fread('/media/liyaru/LYR/Benchmark/dyngen/eval_1016/cosine_sim_velocity.csv')
  
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
  
  fwrite(metrics,outpath)
  
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
      df %>% transmute(x = x, xend = x, y = -.3, yend = -.1, size = .5)
    )
    text_data <-
      bind_rows(
        text_data,
        df %>% transmute(
          xmin = x, xmax = x, ymin = 0, ymax = -0.5,
          angle = 30, vjust = 0, hjust = 0,
          label_value = id, 
          size = 3
        )
      )
  }
  
  
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
  
  x_min_output <- minimum_x+0.5
  x_min_scaling <- minimum_x + 5.5
  x_min_ranking <- minimum_x + 10.5
  x_min_score <-  minimum_x + 17
  
  leg_max_y <- minimum_y - .5
  
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

####--------plot---------------
# gene level
datapath = '/media/liyaru/LYR/Benchmark/dyngen/eval_data_1030/eval_data/gene_level/cosine_sim_velocity.T.csv' # gene level
figpath = '/media/liyaru/LYR/Benchmark/Fig/2.dyngen_cos_v_summary_normalize_by_group_gene_level.pdf'
outpath = '/media/liyaru/LYR/Benchmark/dyngen/eval_data_1030/eval_data/gene_level/cosine_sim_velocity.T_reorder.csv'
plot_summary(datapath,outpath,figpath)

# cell level high-dim
datapath = '/media/liyaru/LYR/Benchmark/dyngen/eval_data_1030/eval_data/cell_level/cosine_sim_velocity.csv' # gene level
figpath = '/media/liyaru/LYR/Benchmark/Fig/2.dyngen_cos_v_summary_normalize_by_group_cell_level.pdf'
outpath = '/media/liyaru/LYR/Benchmark/dyngen/eval_data_1030/eval_data/cell_level/cosine_sim_velocity_reorder.csv'
plot_summary(datapath,outpath,figpath)

# cell level low-dim
datapath = '/media/liyaru/LYR/Benchmark/dyngen/eval_data_1030/eval_data/cell_level/cosine_sim_velocity_pca2.csv' # gene level
figpath = '/media/liyaru/LYR/Benchmark/Fig/2.dyngen_cos_v_summary_normalize_by_group_cell_level_low_dim.pdf'
outpath = '/media/liyaru/LYR/Benchmark/dyngen/eval_data_1030/eval_data/cell_level/cosine_sim_velocity_pca2_reorder.csv'
plot_summary(datapath,outpath,figpath)


####--------pheatmap-----------------
datapath = '/media/liyaru/LYR/Benchmark/dyngen/eval_data_1030/eval_data/cell_level/cosine_sim_velocity.csv'

metrics <- fread(datapath)

colnames(metrics)[1] <- 'Method'

metrics <- metrics %>%
  mutate(Method= gsub("^[0-9]+_", "", Method))

metrics <- check_method_name(metrics)


metrics$Method <- factor(metrics$Method,levels = names(method_palette))

metrics <- metrics[order(metrics$Method)]

metrics <- as.data.frame(metrics)

metrics <- column_to_rownames(metrics,'Method')

annotation_col =data.frame(datasets = colnames(metrics))
annotation_col <- separate(annotation_col,col = 'datasets',into = c('Type','Seed'),sep = '_seed',remove = F)
annotation_col <- column_to_rownames(annotation_col,var = 'datasets')

type_color = hue_pal()(14)
names(type_color) <- unique(annotation_col$Type)
seed_color = brewer.pal(3,'Purples')
names(seed_color) <- unique(annotation_col$Seed)
ann_colors = list(
  Type = type_color,
  Seed = seed_color
)

p<- pheatmap(metrics,
         cluster_cols = F,cluster_rows = F,
         annotation_col = annotation_col,
         annotation_colors = ann_colors)

pdf('/media/liyaru/LYR/Benchmark/Fig/13.dyngen_cell_level_velocity.pdf',width = 10,height = 5)
p
dev.off()
