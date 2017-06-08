# To boxplot the recall, precision, accuracy, f-score and matthews correlation coefficient of the various predictors against the combined positive dataset and the negative dataset
# note that it is possible for a predictor to have zero true positive and zero false positive even if it made predictions, if those predicitions were not overlapping with either positive or negative datasets.
# in these cases, precision will be 1 and MCC will be 0.

boxplot_gipred <- function(result2, pal2, pred_labels) {
  result_mcc <- subset(result2, metrics=="Matthews Correlation Coefficient")
  result_nonmcc <- subset(result2, metrics!="Matthews Correlation Coefficient")
  
  mcc_plot <- ggplot(result_mcc, aes(x=predictor, y=value, fill=predictor)) + geom_boxplot(width=0.5) +
    theme_light() +
    theme(legend.position="bottom", legend.direction="horizontal", legend.title=element_blank(), axis.ticks=element_blank(), axis.text.x = element_blank()) +
    labs(x="",y="") + guides(fill=guide_legend(nrow=3, byrow=TRUE)) +
    scale_fill_manual(values=pal2, labels=pred_labels) +
    facet_wrap(~metrics)
  nonmcc_plots <- ggplot(result_nonmcc, aes(x=predictor, y=value, fill=predictor)) + geom_boxplot(width=0.5) +
    theme_light() +
    theme(legend.position="bottom", legend.direction="horizontal", legend.title=element_blank(), axis.ticks=element_blank(), axis.text.x = element_blank()) +
    labs(x="",y="") + guides(fill=guide_legend(nrow=3, byrow=TRUE)) +
    scale_fill_manual(values=pal2) +
    facet_wrap(~metrics) + theme(legend.position="none")
  leg <- get_legend(mcc_plot)
  mcc_plot <- mcc_plot + theme(legend.position="none")
  
  pdf("boxplot_accuracy.pdf", height = 6, width = 10)
  grid.arrange(mcc_plot, nonmcc_plots, leg, layout_matrix=cbind(c(1,1,3),c(2,2,3),c(2,2,3)), ncol=3, heights=c(1,1,.4), widths=c(3,2,2))
  dev.off()
}