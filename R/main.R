#' @title Encode DNA Sequences
#' @description Convert DNA sequences into categorical nucleotide
#' features for machine learning prediction.
#' @param dna_strings A character vector containing DNA sequences
#' of equal length.
#' @return A data frame containing nucleotide positions as factors.
#' @examples
#' m6APrediction:::dna_encoding(c("GGACA", "ATCGA"))
dna_encoding <- function(dna_strings){
  nn <- nchar( dna_strings[1] )
  seq_m <- matrix( unlist( strsplit(dna_strings, "") ), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}
#' @title Multiple Sample Prediction for m6A Sites
#' @description Predict m6A probabilities and statuses for multiple
#' samples using a trained random forest model.
#' @param ml_fit A trained random forest model.
#' @param feature_df A data frame containing seven required columns:
#'   gc_content, RNA_type, RNA_region, exon_length,
#'   distance_to_junction, evolutionary_conservation, and DNA_5mer.
#' @param positive_threshold A numeric threshold for classifying
#' samples as Positive or Negative. Default is 0.5.
#' @return A data frame containing the original features and two
#' additional columns: predicted_m6A_prob and predicted_m6A_status.
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                              package = "m6APrediction"))
#' example_df <- read.csv(system.file("extdata",
#'                                    "m6A_input_example.csv",
#'                                    package = "m6APrediction"))
#' prediction_multiple(ml_fit, example_df, positive_threshold = 0.6)
#' @export
#' @import randomForest
#' @importFrom stats predict
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) #Check errors if incorrect column names of input data.frame
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA","lincRNA","lncRNA","pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region,levels=c("CDS","intron","3'UTR", "5'UTR"))
  encoded_df <- dna_encoding(feature_df$DNA_5mer)
  prediction_df <- cbind(feature_df, encoded_df)
  predicted_prob <- predict(ml_fit, newdata = prediction_df, type = "prob")
  feature_df$predicted_m6A_prob <- predicted_prob[, "Positive"]
  feature_df$predicted_m6A_status <- ifelse(
    feature_df$predicted_m6A_prob > positive_threshold,
    "Positive",
    "Negative" )
  return(feature_df)
}
#' @title Single Sample Prediction for m6A Sites
#' @description Predict the m6A probability and status for a single
#' sample using a trained random forest model.
#' @param ml_fit A trained random forest model.
#' @param gc_content A numeric value representing GC content.
#' @param RNA_type A character string specifying the RNA type:
#'   "mRNA", "lincRNA", "lncRNA", or "pseudogene".
#' @param RNA_region A character string specifying the RNA region:
#'   "CDS", "intron", "3'UTR", or "5'UTR".
#' @param exon_length A numeric value representing exon length.
#' @param distance_to_junction A numeric distance to a junction.
#' @param evolutionary_conservation A numeric conservation score.
#' @param DNA_5mer A character string representing a DNA 5-mer.
#' @param positive_threshold A numeric threshold for classification.
#' Default is 0.5.
#' @return A named character vector containing predicted_m6A_prob
#' and predicted_m6A_status.
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                              package = "m6APrediction"))
#' prediction_single(ml_fit, gc_content = 0.6,
#'                   RNA_type = "mRNA", RNA_region = "CDS",
#'                   exon_length = 12, distance_to_junction = 5,
#'                   evolutionary_conservation = 0.8,
#'                   DNA_5mer = "GGACA",
#'                   positive_threshold = 0.5)
#' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  single_df <- data.frame(gc_content = gc_content, RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
                          RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),exon_length = exon_length,distance_to_junction = distance_to_junction,evolutionary_conservation = evolutionary_conservation,DNA_5mer = DNA_5mer)
  prediction_result <- prediction_multiple(ml_fit, single_df, positive_threshold)
  returned_vector <- c(predicted_m6A_prob = prediction_result$predicted_m6A_prob[1], predicted_m6A_status = prediction_result$predicted_m6A_status[1])
  return(returned_vector)
}
