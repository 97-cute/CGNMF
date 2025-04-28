
library(proxy)
k_clusters <- 7   # 聚类数量

kmeansFunc<-function (data, k) 
{
  set.seed(12345678)
  if (nrow(data) < 3e+05) {
    numStart = 100
  } else {
    numStart = 1
  }
  cl <- suppressWarnings(try(kmeans(data, k, nstart = numStart, 
                                    iter.max = 100), silent = TRUE))
  
  # 返回聚类标签和聚类质心
  return(list(kmeans = cl$cluster, centers = cl$centers))
}
objective_value <- function(B, P,coord_df,tau,k1,k2,numK) {
  # 计算欧式距离
  library(stats)
  dist_matrix <- stats::dist(coord_df, method = "euclidean")
  sim_euclidean <- 1 / (0.001 + as.matrix(dist_matrix))
  # 计算余弦相似度的函数
  cosine_similarity <- function(x, y) {
    sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
  }
  #求导函数的一部分
  custom_function1 <- function(x, y) {
    # 计算二范数
    norm_x <- sqrt(sum(x^2))
    norm_y <- sqrt(sum(y^2))
    numerator <- sum(y) * (norm_x^2) 
    denominator <- norm_y * (norm_x^3)
    return(numerator / denominator)
  }
  custom_function2 <- function(x, y) {
    # 计算二范数
    norm_x <- sqrt(sum(x^2))
    norm_y <- sqrt(sum(y^2))
    numerator <-sum(x %*% t(y)) * x
    denominator <- norm_y * (norm_x^3)
    return(numerator / denominator)
  }
  # 计算对比学习损失
  N <- nrow(P)  # 样本数量
  L_contrast <- 0
  # 遍历每个样本 Pi
  Q <- matrix(NA, nrow=N, ncol=numK)
  W <- matrix(NA, nrow=N, ncol=numK)
  sim <- (1-as.matrix(dist(P, method = "cosine")))
  rownames(sim) <- rownames(P)
  colnames(sim) <- rownames(P)
  for (i in 1:N) {
    candidate_positive_samples <- list()
    similar_indices_P <- order(sim[i, ], decreasing = TRUE)[-1][1:k1]
    similar_indices_A <- order(sim_euclidean[i, ], decreasing=TRUE)[-1][1:k1]
    
    intersection <- intersect(similar_indices_P, similar_indices_A)
    
    if (length(intersection) > 0) {
      intersection <- intersection
    } else {
      intersection <- similar_indices_A
    }
    candidate_positive_samples<-P[intersection,]
    #cluster_index <- clusters$cluster[i]
    #centroid <- centroids[cluster_index, ]
    euclidean_distances <- c()
    # 计算当前样本与质心的距离
    if(length(intersection)==1){
      euclidean_distances <- exp(sqrt(sum((candidate_positive_samples - P[i,]) ^ 2)))
    }else{
      for (j in 1:length(intersection)) {
        distance <- exp(sqrt(sum((candidate_positive_samples[j, ] - P[i,]) ^ 2)))
        euclidean_distances[j] <- distance
      }
    }
    
    w= exp(1/euclidean_distances)/sum(exp(1/euclidean_distances))
    # 计算加权正样本
    weighted_pos_samples <- w%*%candidate_positive_samples
    
    # 计算负样本集
    
    neg_indices_P <- order(sim[i, ], decreasing=FALSE)[1:k2]
    neg_indices_A <- order(sim_euclidean[i, ], decreasing=FALSE)[1:k2]
    
    neg_intersection <- intersect(neg_indices_P, neg_indices_A)
    
    if (length(neg_intersection) > 0) {
      neg_intersection <- neg_intersection
    } else {
      neg_intersection <- neg_indices_P
    }
    neg_samples<-P[neg_intersection,]
    
    # 计算正样本相似度
    cos_sim_positive <- cosine_similarity(P[i, ], weighted_pos_samples)  
    
    # 计算负样本的累积相似度与负样本求导部分
    sum_neg_sim <- 0
    sum_neg_dao1 <- 0
    sum_neg_dao2 <- 0
    if(length(neg_intersection)==1){
      cos_sim_negative <- cosine_similarity(P[i, ], neg_samples)
      sum_neg_sim <- exp(cos_sim_negative / tau)
      sum_neg_dao1 <-  exp(cos_sim_negative / tau)*custom_function1(P[i, ], neg_samples)
      sum_neg_dao2 <-exp(cos_sim_negative / tau)*custom_function2(P[i, ], neg_samples)
    }else{
      for (j in 1:length(neg_intersection)) {
        cos_sim_negative <- cosine_similarity(P[i, ], neg_samples[j, ])
        sum_neg_sim <- sum_neg_sim + exp(cos_sim_negative / tau)
        sum_neg_dao1 <- sum_neg_dao1 +exp(cos_sim_negative / tau)*custom_function1(P[i, ], neg_samples[j, ])
        sum_neg_dao2 <- sum_neg_dao2 +exp(cos_sim_negative / tau)*custom_function2(P[i, ], neg_samples[j, ])
      }
    }
    # 计算对比损失
    loss_positive <- exp(cos_sim_positive / tau)
    loss <- -log(loss_positive / (loss_positive + sum_neg_sim))
    #计算对比损失求导
    Q[i, ]<-(1/(N*tau))*(
      custom_function2(P[i, ], weighted_pos_samples)
      +(1/(loss_positive + sum_neg_sim))*(loss_positive*custom_function1(P[i, ], weighted_pos_samples)+sum_neg_sim*sum_neg_dao1)
    )
    W[i, ]<-(1/(N*tau))*(
      custom_function1(P[i, ], weighted_pos_samples)
      +(1/(loss_positive + sum_neg_sim))*(loss_positive*custom_function2(P[i, ], weighted_pos_samples)+sum_neg_sim*sum_neg_dao2)
    )
  }
  return(list(Q=Q,W=W))
} 




updateB_P<-function(expr_matrix,zero_threshold,numK,ineibor,k1,k2, max_iter,epsilon,tau,a1,a2){
  expr_matrix <- expr_matrix[rowMeans(expr_matrix == 0) < zero_threshold, ]
  rownames(expr_matrix) <- paste0("Gene", 1:nrow(expr_matrix))
  rawData <- list("dataset1" = expr_matrix)
  library(rliger)
  ifnb_liger <- createLiger(rawData, removeMissing = F, verbose = F)
  #class(ifnb_liger)
  ifnb_liger <- rliger::normalize(ifnb_liger)
  #class(ifnb_liger)
  # 选择变量基因
  ifnb_liger <- selectGenes(ifnb_liger)
  # 运行 scaleNotCenter
  ifnb_liger <- scaleNotCenter(ifnb_liger)
  Y<-as.matrix(ifnb_liger@datasets[["dataset1"]]@scaleData)
  
  # 初始化 B 矩阵
  set.seed(123)
  n_rows <- nrow(Y)
  n_cols <- ncol(Y)
  
  #numK=20
  B <- matrix(runif(n_rows * numK ), nrow=n_rows, ncol=numK )
  P <- matrix(runif(n_cols *numK), nrow=n_cols, ncol=numK )
  library(RANN)
  library(Matrix)
  createA <- function(locationList,ineibor) {
    nSlices = length(locationList)
    AList = list()
    for (islice in 1:nSlices) {
      location = as.data.frame(locationList[[islice]])
      norm_cords = location[, c("x", "y")]
      rownames(norm_cords) <- rownames(location)
      #ineibor = 15
      near_data = nn2(norm_cords[, 1:2], k = ineibor)
      neibors = near_data$nn.idx
      neibors = neibors[, -1]
      Nmat = Matrix(0, nrow = nrow(neibors), ncol = nrow(neibors), sparse = TRUE)
      for (icol in 1:ncol(neibors)) {
        edges = data.frame(i = 1:nrow(neibors), j = neibors[, icol])
        adjacency = sparseMatrix(i = as.integer(edges$i),
                                 j = as.integer(edges$j), x = 1, dims = rep(nrow(neibors), 2), use.last.ij = TRUE)
        Nmat = Nmat + adjacency
      }
      Nmat = Nmat * t(Nmat)
      rownames(Nmat) = colnames(Nmat) = rownames(norm_cords)
      AList[[islice]] = Nmat
    }
    return(AList)
  }
  # cat(paste0("### 图构建开始! ...\n"))
  colnames(coord_df) <- c("x", "y")
  coord_df <- as.data.frame(coord_df)
  locationList <- list(coord_df)
  AList <- createA(locationList,ineibor)  #欧式距离创建邻接矩阵 A
  A <- as.matrix(AList[[1]])# 转换为矩阵形式
  D <- diag(rowSums(A))# 度矩阵 D
  L <- D - A# 拉普拉斯矩阵 L
  
  #loss_values <- c()  # 存储每次迭代的损失值
  objective<-objective_value(B, P,coord_df,tau,k1,k2,numK)
  Q<-objective$Q
  W<-objective$W
  for (iter in 1:max_iter) {
    # 保存前一次的W和H
    B_old <- B
    P_old <- P
    
    # 更新H矩阵
    B <- (Y%*%P)*B/(B%*%t(P)%*%P)  
    B[B < 0] <- 0  # 保证非负性
    # 更新W矩阵
    P <- (t(Y)%*%B+a1*(A%*%P)+(1/2)*a2*W)*P/(P%*%t(B)%*%B+a1*(D%*%P)+(1/2)*a2*Q)
    P[P < 0] <- 0  # 保证非负性
    # 计算损失值
    objective <- objective_value(B, P,coord_df,tau,k1,k2,numK)
    Q<-objective$Q
    W<-objective$W
    # 计算参数变化
    P_change <- max(abs(P - P_old))
    
    # 收敛条件检查
    if (P_change < epsilon) {
      # cat("收敛于迭代", iter, "次", "\n")
      break
    }
    
    # 打印当前损失值（可选）
    # if (iter %% 50 == 0) {
    #   cat("迭代", iter, "次", "\n")
    # }
  }
  return(P = P)
}

library(SeuratDisk)
seurat_object <- LoadH5Seurat("E:/STARmap_20180505_BY3_1k.h5seurat")#数据导入路径
library(Seurat)
seurat_object = UpdateSeuratObject(seurat_object)
TC<-seurat_object@meta.data
expr_matrix<-seurat_object@assays[["RNA"]]@counts
true_labels <- TC$label
coord_matrix<-seurat_object@reductions[["spatial"]]@cell.embeddings
# 将其转换为数据框以便进行操作
colnames(coord_matrix) <- c("x", "y")
coord_df <- as.data.frame(coord_matrix)
coord_df$x <- as.numeric(coord_df$x)
coord_df$y <- as.numeric(coord_df$y)
str(coord_df)
coord_df$x <- coord_df$x - min(coord_df$x) 
coord_df$y <- coord_df$y - min(coord_df$y) 
scaleFactor <- max(coord_df$x, coord_df$y)
coord_df$x <- coord_df$x / scaleFactor
coord_df$y <- coord_df$y / scaleFactor
coord_df <- as.matrix(coord_df[, c("x", "y")])

library(flexclust)  # 如果kmeansFunc或后续分析需要此库
library(fpc)        # 用于计算ARI

run_model <- function(zero_threshold, numK, ineibor, k1, k2, tau, a1, a2) {
  result <- updateB_P(expr_matrix,zero_threshold=0.9,numK=20,ineibor=15,k1=10,k2=60,max_iter = 500, epsilon = 0.05, tau = 0.3, a1 = 6, a2 = 5)
  s1 <- kmeansFunc(result, k_clusters)
  J1 <- s1[["kmeans"]]  # 提取聚类结果
  ari <- randIndex(table(J1, true_labels))
  # cat("a1 =", a1, ", a2 =", a2, ", ARI =", ari, "\n")
  return(ari)
}

run_model(zero_threshold=0.9,numK=20,ineibor=15,k1=10,k2=60, tau = 0.3, a1 = 6, a2 = 5)



