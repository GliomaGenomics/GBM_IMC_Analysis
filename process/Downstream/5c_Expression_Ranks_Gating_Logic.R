labels = list()
labels$undefined = rownames(gating_info$marker_expression)

# Endothelial ----

labels$Endothelial = identify_cells(gating_info$maker_ranks[,c("SMA","CD31")] >= 15)
labels = update_undefined(labels)

# Immune ----

labels$`NK cell` = identify_cells(
    gating_info$maker_ranks[labels$undefined,c("CD45","NKP46")] >= 14
    )

labels = update_undefined(labels)


labels$`T cell` = identify_cells(
    gating_info$maker_ranks[labels$undefined ,c("CD45","CD3")] >= 11
    )

labels$Macrophage = identify_cells(
    gating_info$maker_ranks[labels$undefined, c("CD45","IBA1")] >= 15
    )

labels = compare_vectors(labels, "T cell", "Macrophage", return_common = F)

labels = update_undefined(labels)

check_distinct_vectors(labels)


labels$Microglia = identify_cells(
    gating_info$maker_ranks[labels$undefined, c("CD45","IBA1", "TMEM119")] >= 12
    )


labels = update_undefined(labels)
check_distinct_vectors(labels)


# Stroma ----

labels$Neuron = identify_cells(
    gating_info$maker_ranks[labels$undefined, c("NeuN_FOX3"), drop = FALSE] >= 18
)


labels$neuron_filt_1 = identify_cells(
    gating_info$maker_ranks[labels$Neuron, c("NKP46"), drop = FALSE] >= 16
)

labels$neuron_filt_2 = identify_cells(
    gating_info$maker_ranks[labels$Neuron, c("CD3"), drop = FALSE] >= 16
)

labels = combine_labels(labels, 
                        c("neuron_filt_1", "neuron_filt_2"), 
                        combined_name = "neuron_filt")


labels$Neuron = setdiff(labels$Neuron, labels$neuron_filt)

labels$neuron_filt = NULL

labels = update_undefined(labels)
check_distinct_vectors(labels)



labels$Astrocyte = identify_cells(
    gating_info$maker_ranks[labels$undefined, c("GFAP"), drop = FALSE] >= 17
)


labels$Oligodendrocyte = identify_cells(
    gating_info$maker_ranks[labels$undefined, c("MOG"), drop = FALSE] >= 18
)


# labels = compare_vectors(labels, "Astrocyte", "Neuron", return_common = F)
# labels = compare_vectors(labels, "Oligodendrocyte", "Neuron", return_common = F)
labels = compare_vectors(labels, "Astrocyte", "Oligodendrocyte", return_common = F)


labels = update_undefined(labels)
check_distinct_vectors(labels)

# Neoplastic ----

labels$AC = identify_cells(
    gating_info$maker_ranks[labels$undefined,
                            c("SLC1A3_EAAT1", "HOPX")] >= 17
)



labels$MES1 = identify_cells(
    gating_info$maker_ranks[labels$undefined, 
                            c("SOD2","CHI3L1")] >= 17
)


labels$MES2 = identify_cells(
    gating_info$maker_ranks[labels$undefined, 
                            c("ANEXIN_A2","ANXA1")] >= 17
)


labels = compare_vectors(labels, "MES1", "MES2", return_common = T)

labels = combine_labels(labels, 
               combine_elements = c("MES1","MES2","MES1_MES2"),
               combined_name = "MES")



labels$NPC = identify_cells(
    gating_info$maker_ranks[labels$undefined, 
                            c("DLL3", "BCAN")] >= 15
)


labels$OPC = identify_cells(
    gating_info$maker_ranks[labels$undefined, 
                            c("SCD5", "OLIG1")] >= 15
)


unique_combs = combn(c("AC","MES","NPC","OPC"), 2, simplify = F)

for (i in seq_along(unique_combs)) {
    
 labels =  compare_vectors(labels, 
                           vector1 = unique_combs[[i]][[1]],
                           vector2 =  unique_combs[[i]][[2]],
                           return_common = FALSE)
    
} 

labels = update_undefined(labels)
check_distinct_vectors(labels)




labels$AC_filt1 = identify_cells(
    gating_info$maker_ranks[labels$AC, c("CD3"), drop = F] >= 17
)

labels$AC_filt2 = identify_cells(
    gating_info$maker_ranks[labels$AC, c("NKP46"), drop = F] >= 17
)

labels = combine_labels(labels, c("AC_filt1", "AC_filt2"), "AC_filt")

labels$AC = setdiff(labels$AC, labels$AC_filt)
labels$AC_filt = NULL


labels = update_undefined(labels)
check_distinct_vectors(labels)

# Cell States ----

labels = list()

labels$high_JARID2_C_Terminus = identify_cells(
    gating_info$maker_ranks[,c("JARID2_C_Terminus"), drop=F] >= 7
    )
labels$low_JARID2_C_Terminus = identify_cells(
    gating_info$maker_ranks[,c("JARID2_C_Terminus"), drop=F] <= 3
)



labels$high_JARID2_N_Terminus = identify_cells(
    gating_info$maker_ranks[,c("JARID2_N_Terminus"), drop=F] >= 7
)
labels$low_JARID2_N_Terminus = identify_cells(
    gating_info$maker_ranks[,c("JARID2_N_Terminus"), drop=F] <= 3
)



labels$high_proliferation = identify_cells(
    gating_info$maker_ranks[,c("Ki67"), drop=F] >= 5
)



labels$high_hypoxia = identify_cells(
    gating_info$maker_ranks[,c("HIF1A"), drop=F] >= 5
)




labels$transcriptionally_respressive = identify_cells(
    gating_info$maker_ranks[,c("EZH2"), drop=F] >= 7
)


labels$EMT = identify_cells(
    gating_info$maker_ranks[,c("SNAI1"), drop=F] >= 7
)


labels$proliferating_stem_cell = identify_cells(
    gating_info$maker_ranks[,c("SOX2"), drop=F] >= 7
)


labels$quiescent_stem_cell_1 = identify_cells(
    gating_info$maker_ranks[,c("TNC"), drop=F] >= 7
)

labels$quiescent_stem_cell_2 = identify_cells(
    gating_info$maker_ranks[,c("TGFBeta"), drop=F] >= 7
)


labels = compare_vectors(labels, 
                         "quiescent_stem_cell_1", 
                         "quiescent_stem_cell_2", 
                         return_common = F)

labels = combine_labels(labels, 
                        combine_elements = c("quiescent_stem_cell_1","quiescent_stem_cell_2"),
                        combined_name = "quiescent_stem_cell")


spe@metadata$cell_states = labels


# End----
