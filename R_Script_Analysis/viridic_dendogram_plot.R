# Script adapted and customized from viridic R plotting script

#install.packages("systemfonts")
#install.packages("dplyr")
#library(systemfonts)
#library(dplyr)
#anc_reg_path <- system_fonts() %>% filter(family == "Avenir Next Condensed", style == "Regular") %>% pull(path)
#anc_ital_path <- system_fonts() %>% filter(family == "Avenir Next Condensed", style == "Italic") %>% pull(path)
#register_font(name = "Avenir Next Condensed", plain = anc_reg_path, italic = anc_ital_path)

#install.packages("dendextend")
#install.packages("extrafont")
#install.packages("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
library(grid)
library(dendextend)

setwd("~/Downloads/04_VIRIDIC_out")

# Path to your VIRIDIC output directory
viridic_out <- "~/Downloads/04_VIRIDIC_out"

font_family = "Avenir Next Condensed"

# Load RDS files
sim_matrix_DF <- readRDS(paste0(viridic_out, "/sim_MA.RDS"))
qsLenFrac_DF <- readRDS(paste0(viridic_out, "/fractqslen_MA.RDS"))
qFracAlig_DF <- readRDS(paste0(viridic_out, "/qFractAlign_MA.RDS"))
sFracAlig_DF <- readRDS(paste0(viridic_out, "/sFractAlign_MA.RDS"))
len_DF <- readRDS(paste0(viridic_out, "/len_DF.RDS"))

# Convert to matrix
sim_matrix_MA <- as.matrix(sim_matrix_DF)

# ---- 1. Create dendrogram ----
# Use hierarchical clustering on similarity matrix
dend_row <- as.dendrogram(hclust(dist(sim_matrix_MA)))
dend_col <- as.dendrogram(hclust(dist(t(sim_matrix_MA))))

# ---- 2. Custom color settings ----
# Blues for similarity
sim_colors <- c("#FFFFFF", "#E6F2FF", "#CCE5FF", "#99CCFF", "#66B2FF", "#3399FF", "#0080FF", "#4D94FF", "#0066CC", "#004C99")
sim_breaks <- c(0, 4.9999, 9.9999, 19.9999, 39.9999, 49.9999, 59.9999, 79.9999, 89.9999, 100)

# Color for genome length ratio
len_colors <- c("#e60000", "#ff4d4d", "#ff9999", "#ffb3b3", "#ffffff")
len_breaks <- c(0.6, 0.7, 0.8, 0.9, 1)

# Orange gradient for aligned genome fraction
alig_colors <- c("#b37700", "#e69500", "#ffad33", "#f0c27b", "#ffd699", "#ffffff")
alig_breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)

# Create color functions
col_sim_fun <- colorRamp2(breaks = sim_breaks, colors = sim_colors)
col_sqLenfrac_fun <- colorRamp2(breaks = len_breaks, colors = len_colors)
col_Aligfrac_fun <- colorRamp2(breaks = alig_breaks, colors = alig_colors)

# Function to create heatmap with specified font
create_heatmap <- function(font_family) {
  # ---- 3. Heatmap Annotation with Dendrogram ----
  ha <- HeatmapAnnotation(
    'Genome length' = anno_barplot(len_DF$qlen, 
                                   which = "column", 
                                   border = FALSE, 
                                   baseline = 0,
                                   gp = gpar(fill = "gray80"),
                                   height = unit(2, "cm")),   # Height of the barplot
    which = "col", 
    show_annotation_name = TRUE, 
    annotation_name_rot = 0,  # Make the name horizontal
    annotation_name_gp = gpar(fontfamily = font_family),
    annotation_height = unit(2, "cm"),  # Total height of annotation including label
    gap = unit(1, "cm"),  # Gap between annotation and heatmap
    annotation_name_side = "right",  # Position label on right side
    annotation_name_offset = unit(0.1, "cm")  # Distance of label from annotation
  )
  
  # ---- 4. Combine dendrogram with heatmap ----
  ht <- Heatmap(
    sim_matrix_MA, 
    name = "Intergenomic similarity", 
    col = col_sim_fun,
    top_annotation = ha,
    
    show_heatmap_legend = FALSE, 
    
    # Apply font to row and column names
    row_names_gp = gpar(fontsize = 10, fontfamily = font_family),
    column_names_gp = gpar(fontsize = 10, fontfamily = font_family),
    
    cluster_rows = dend_row,    # Use dendrogram
    cluster_columns = dend_col, # Use dendrogram
    
    cell_fun = function(j, i, x, y, width, height, fill) {
      if(i > j) {
        # Lower triangle - fractions
        grid.rect(x = x, y = y, width = width,  height = height, 
                  gp = gpar(col = "gray80", fill = NA))
        
        # qAligFract
        grid.rect(x = x, y = y+height*0.33, width = width, height = height*0.33, 
                  gp = gpar(fill = col_Aligfrac_fun(qFracAlig_DF[i,j]), col = NA))
        
        # Display aligned fraction number
        grid.text(sprintf("%.1f", qFracAlig_DF[i, j]), x, y = y+height*0.33, 
                  gp = gpar(fontfamily = font_family, fontsize = 7))
        
        # qsLenFract
        grid.rect(x = x, y = y, width = width, height = height*0.33, 
                  gp = gpar(fill = col_sqLenfrac_fun(qsLenFrac_DF[i,j]), col = NA))
        
        # Display length ratio number
        grid.text(sprintf("%.1f", qsLenFrac_DF[i, j]), x, y, 
                  gp = gpar(fontfamily = font_family, fontsize = 7))
        
        # sAligFract
        grid.rect(x = x, y = y-height*0.33, width = width, height = height*0.33, 
                  gp = gpar(fill = col_Aligfrac_fun(sFracAlig_DF[i,j]), col = NA))
        
        # Display aligned fraction number
        grid.text(sprintf("%.1f", sFracAlig_DF[i, j]), x, y = y-height*0.33, 
                  gp = gpar(fontfamily = font_family, fontsize = 7))
      } else {
        # Upper triangle - similarity
        grid.rect(x = x, y = y, width = width, height = height, 
                  gp = gpar(fill = col_sim_fun(sim_matrix_MA[i,j]), col = "gray80"))
        
        # Display similarity numbers
        grid.text(sprintf("%.1f", sim_matrix_MA[i, j]), x, y, 
                  gp = gpar(fontfamily = font_family, fontsize = 10))
      }
    },
    
    show_row_names = TRUE, 
    show_column_names = TRUE
  )
  
  return(ht)
}

# Function to create legends with specified font
create_legends <- function(font_family) {
  # ---- 5. Legends ----
  lgd_Alig <- Legend(col_fun = col_Aligfrac_fun, title = "Aligned genome fraction", 
                     at = c(0, 0.25, 0.5, 0.75, 1), direction = "horizontal",
                     title_gp = gpar(fontfamily = font_family),
                     labels_gp = gpar(fontfamily = font_family))
  
  lgd_LenFrac <- Legend(col_fun = col_sqLenfrac_fun, title = "Genome length ratio", 
                        at = c(0.6, 0.7, 0.8, 0.9, 1), direction = "horizontal",
                        title_gp = gpar(fontfamily = font_family),
                        labels_gp = gpar(fontfamily = font_family))
  
  lgd_Similarity <- Legend(col_fun = col_sim_fun, title = "Intergenomic similarity",
                           direction = "horizontal", 
                           title_gp = gpar(fontfamily = font_family),
                           labels_gp = gpar(fontfamily = font_family))
  
  # Create separate legend combinations
  separate_legends <- list(lgd_Alig, lgd_LenFrac)
  
  # Combine the legends horizontally with spacing
  combined_legends <- packLegend(
    lgd_Alig,
    lgd_LenFrac,
    lgd_Similarity,
    direction = "horizontal",
    gap = unit(1, "cm") 
  )
  
  return(list(separate = separate_legends, combined = combined_legends))
}

# ---- 6. Save the plot ----
# Calculate appropriate dimensions
num_char <- max(nchar(colnames(sim_matrix_DF)))
ht_width <- max(7, 0.5 * ncol(sim_matrix_DF))
ht_height <- max(7, 0.5 * nrow(sim_matrix_DF))

# ---- Create PNG with Avenir Next Condensed Bold font ----
# Create heatmap and legends
ht_png <- create_heatmap(font_family)
legends_png <- create_legends(font_family)

# Save as PNG with legends at top

par(mar = c(5, 4, 4, 6)) 

png("viridic_with_dendrogram_top.png", width = ht_width*200, height = ht_height*350, res = 300)
draw(ht_png, 
     annotation_legend_list = legends_png$separate, 
     annotation_legend_side = "top", 
     heatmap_legend_side = "top", 
     merge_legends = FALSE)
dev.off()

png("viridic_with_dendrogram_aligned.png", width = ht_width*300, height = ht_height*350, res = 300)
draw(ht_png, 
     annotation_legend_list = combined_legends, 
     annotation_legend_side = "bottom", 
     heatmap_legend_side = "bottom",
     merge_legends = TRUE)

dev.off()

#######
# For the combined legends version
png("viridic_with_dendrogram_side_combined.png", width = ht_width*1200, height = ht_height*800, res = 1000)

# Create stacked vertical legends
stacked_legends <- packLegend(
  lgd_Alig,
  lgd_LenFrac,
  lgd_Similarity,
  direction = "vertical",  # Stack vertically instead of horizontally
  gap = unit(0.5, "cm")   # Gap between legends
)

draw(ht_png, 
     annotation_legend_list = stacked_legends, 
     annotation_legend_side = "right",  # Put on right side
     heatmap_legend_side = "right",    # Put on right side
     merge_legends = TRUE,
     padding = unit(c(2, 2, 2, 15), "mm"))  # Extra padding on right for legends
dev.off()
