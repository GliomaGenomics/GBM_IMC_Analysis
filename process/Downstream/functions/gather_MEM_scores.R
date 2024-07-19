gather_MEM_scores = function(MEM_score_matrix, threshold = 3){
    
    MEM_scores = t(MEM_score_matrix)
    
    scores = list(
        Up = vector("character", ncol(MEM_scores)),
        Down = vector("character", ncol(MEM_scores))
    )
    
    for (cluster in 1:ncol(MEM_scores)) {
        
        cluster_score = sort(round(MEM_scores[,cluster]))
        
        cluster_score = cluster_score[which(abs(cluster_score) > threshold)]
        
        up = sort(cluster_score[which(cluster_score > 0)], decreasing = T)
        down = cluster_score[which(cluster_score < 0)]
        
        scores[["Up"]][[cluster]] = ifelse(length(up) >= 1, 
                                           paste0(paste(names(up), up, sep = "+"), 
                                                  collapse = ", "), 
                                           "none")
        
        scores[["Down"]][[cluster]] = ifelse(length(down) >= 1, 
                                             paste0(paste(names(down), down, sep = ""), 
                                                    collapse = ", "), 
                                             "none")
    }
    
    cbind(cluster = colnames(MEM_scores), bind_rows(scores))
    
}