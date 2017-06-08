dotplot_gipred <- function(result, pal, pred_labels) {
  result_mcc <- subset(result, metrics=="Matthews Correlation Coefficient")
  result_nonmcc <- subset(result, metrics!="Matthews Correlation Coefficient")
  
  mcc_plot <- ggplot(result_mcc, aes(x=predictor, y=value, color = predictor, fill = predictor)) + 
    geom_dotplot(binaxis='y', binwidth = 0.05, stackdir="center", dotsize = 0.5, method="dotdensity", stackratio = 0.8) + theme_light() +
    scale_color_manual(values=pal, labels=pred_labels) +
    scale_fill_manual(values=pal, labels=pred_labels) +
    theme(legend.position="bottom", legend.box="horizontal", legend.title=element_blank(), axis.ticks=element_blank(), axis.text.x = element_blank()) + 
    labs(x="", y="") + guides(fill=guide_legend(nrow=2, byrow = TRUE), color=guide_legend(nrow=2, byrow = TRUE), size = "none") +
    facet_grid(~metrics)
  nonmcc_plots <- ggplot(result_nonmcc, aes(x=predictor, y=value, color = predictor, fill = predictor)) + 
    geom_dotplot(binaxis='y', binwidth = 0.05, stackdir="center", dotsize = 0.5, method="dotdensity", stackratio = 0.8) + theme_light() +
    scale_color_manual(values=pal, labels=pred_labels) +
    scale_fill_manual(values=pal, labels=pred_labels) +
    theme(legend.position="bottom", legend.box="horizontal", legend.title=element_blank(), axis.ticks=element_blank(), axis.text.x = element_blank()) + 
    labs(x="", y="") + guides(fill=guide_legend(nrow=2, byrow = TRUE), color=guide_legend(nrow=2, byrow = TRUE), size = "none") +
    facet_wrap(~metrics) + theme(legend.position="none")
  leg <- get_legend(mcc_plot)
  mcc_plot <- mcc_plot + theme(legend.position="none")
  pdf("dotplot_accuracy_lit.pdf", height = 6, width = 11)
  grid.arrange(mcc_plot, nonmcc_plots, leg, layout_matrix=cbind(c(1,1,3),c(2,2,3),c(2,2,3)), ncol=3, heights=c(1,1,.4), widths=c(3,2,2))
  dev.off()
}