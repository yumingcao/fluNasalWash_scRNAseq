#' define infected cells using zero inflated negative binomial model and SVM 
#'
#' @param input the expression set with raw UMI counts
#' @param pd the phenotype data
#' @param virus_count column name of pd that contains virus count per cell
#' @param UMIsum column name of pd that contains total umi counts per cell
#' @param celltype column name of pd that contains cell type information per cell
#' @param pcnames column prefix for PC
#' @param output_dir output directory to save zinb model rds file
#' @param output_prefix output prefix for zinb model rds file
#' @export
#' 


define_infection <- function(input, 
                             pd = NULL, 
                             virus_count, 
                             UMIsum = NULL, 
                             celltype = NULL, 
                             pcnames, 
                             output_dir = "zinb_model/", 
                             output_prefix) {
  if (is.null(pd)) {
    pd = pData(input)
  }
  if (is.null(UMIsum)) {
    input$UMIsum <- colSums(exprs(input))
  }
  colnames(pd)[which(colnames(pd) == virus_count)] <- "virus_count"
  colnames(pd)[which(colnames(pd) == UMIsum)] <- "UMIsum"
  
  # perform ZINB
  message('Performing ZINB modeling on viral counts...');flush.console()
  
  model_zinb = hurdle((virus_count) ~ log(UMIsum) | log(UMIsum), 
                      data = pd, 
                      dist ="negbin")
  outname <- paste(output_prefix, "zinb_model.rds", sep = "\t")
  mainDir <- getwd()
  dir.create(file.path(mainDir, output_dir), showWarnings = F)
  saveRDS(model_zinb, file = paste(output_dir, outname, sep = ""))
  
  pd$predict_infection_zinb <- predict(model_zinb, pd, type = "response")
  pd$pzero <- round(predprob(model_zinb, newdata = pd)[,1])
  phat = predprob(model_zinb, newdata = pd)[,1]
  pd$pvals = 1 - pzanegbin(q = pd$virus_count, 
                         size = 1/model_zinb$theta, 
                         munb = pd$predict_infection_zinb, 
                         pobs0 = phat)
  pd$padj_infection_zinb=p.adjust(pd$pvals, method = "fdr")
  
  input$pred_infection <- "no"
  input$pred_infection[pd$padj_infection_zinb <= 0.05] <- "yes"
  message('ZINB modeling on viral counts Done');flush.console()
  
  # perform SVM 
  # if cell types are provided, perform per cell type
  message('Performing SVM classification...');flush.console()
  if (!is.null(celltype)) {
    for (i in unique(pd[,celltype])) {
      cell_subset = input[,which(pd[,celltype] == i)]
      bystander_train_idx <- which(cell_subset$pred_infection == "no" & 
                                     cell_subset$virus_count == 0 )
      infected_train_idx <- which(cell_subset$pred_infection == "yes" & 
                                    cell_subset$virus_count > 0 )
      input$Pass_Gate <- NA
      cell_subset$Pass_Gate[bystander_train_idx] <- "Bystander"
      cell_subset$Pass_Gate[infected_train_idx] <- "Infected"
      
      cell_subset <- flow_svm(cell_subset, pcnames = pcnames) #flow_svm() from SignallingSingleCell
    }
    input$Pass_Gate[match(colnames(cell_subset), colnames(input))] <- cell_subset$Pass_Gate
    input$SVM_Classify[match(colnames(cell_subset), colnames(input))] <- cell_subset$SVM_Classify
    
  } else {
    
    bystander_train_idx <- which(input$pred_infection == "no" & 
                                   input$virus_count == 0 )
    infected_train_idx <- which(input$pred_infection == "yes" & 
                                  input$virus_count > 0 )
    input$Pass_Gate <- NA
    input$Pass_Gate[bystander_train_idx] <- "Bystander"
    input$Pass_Gate[infected_train_idx] <- "Infected"
    
    input <- flow_svm(input, pcnames = pcnames)
  }
  
  return(input)
}