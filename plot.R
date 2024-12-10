setwd("~/breeder_equation_project/")
library(ggplot2)
library(cowplot)

data1 <- readRDS("estimate_accuracy/accuracy_table_long.rds")
data2 <- readRDS("estimate_accuracy2/accuracy_table_long.rds")
data5 <- readRDS("estimate_accuracy5/accuracy_table_long.rds")
data4 <- readRDS("estimate_accuracy4/accuracy_table_long.rds")

data1 <- cbind(data = "poplar_Rincent", 
               data1)
data2 <- cbind(data = "wheat_Rincent", 
               data2)
data5 <- cbind(data = "Winn", 
               trait = "GrainYield", 
               data5)
colnames(data4)[1] <- "jth_fold"
data4 <- cbind(data = "Krause", 
               trait = "GrainYield", 
               data4)

accuracy <- rbind(rbind(rbind(data1, data2), data5), data4)
accuracy[accuracy$accuracy > 1, ]
accuracy[accuracy$accuracy > 1, "accuracy"] <- 1
accuracy[accuracy$accuracy == 1, ]

accuracy_summary <- aggregate(. ~ data + trait + adjustment + `prediction method`, 
                              data=accuracy[, !colnames(accuracy) %in% c("type", "jth_fold")], 
                              FUN=function(x){c(mean=mean(x), se=sd(x)/sqrt(length(x)))})
saveRDS(accuracy_summary, "plot/accuracy_summary.rds")

accuracy_sum <- data.frame(accuracy_summary[[1]])
for (i in 2:length(accuracy_summary)){
  accuracy_sum = cbind(accuracy_sum, accuracy_summary[[i]])
}
colnames(accuracy_sum) <- c("data", "trait", "adjustment", "prediction_method", "accuracy.mean", 
                            "accuracy.se", "h2.mean", "h2.se")



p1 = ggplot(accuracy_sum[accuracy_sum$trait == "CIRC.ORL" & 
                           accuracy_sum$adjustment == "predictive ability", ], 
            aes(x=prediction_method)) +
  geom_point(aes(y=accuracy.mean, colour=prediction_method), size=0.75) + 
  geom_errorbar(aes(ymin=accuracy.mean-accuracy.se, ymax=accuracy.mean+accuracy.se, 
                    colour=prediction_method), width=0.1) + 
  facet_wrap(~adjustment) +
  ylim(0, 1) + 
  ylab(expression("r"["y"])) + 
  scale_colour_manual(values=c("blue", "gold3")) + 
  theme_minimal_grid(font_size=7) +
  theme(axis.title.x = element_blank(),
    #     axis.text.x = element_blank(),
    #     axis.ticks.x = element_blank(),
        legend.position="none") + 
  ggtitle(paste("trait: poplar circumference"))

p2 = ggplot(accuracy_sum[accuracy_sum$trait == "CIRC.ORL" & 
                           accuracy_sum$adjustment == "prediction accuracy", ],
            aes(x=prediction_method)) +
  geom_point(aes(y=accuracy.mean, colour=prediction_method), size=0.75) + 
  geom_errorbar(aes(ymin=accuracy.mean-accuracy.se, ymax=accuracy.mean+accuracy.se, 
                    colour=prediction_method), width=0.1) + 
  facet_wrap(~adjustment) +
  ylim(0, 1) + 
  ylab("r") + 
  scale_colour_manual(values=c("blue", "gold3")) + 
  theme_minimal_grid(font_size=7) +
  theme(axis.title.x = element_blank(),
    #     axis.text.x = element_blank(),
    #     axis.ticks.x = element_blank(),
        legend.position="none") 

# save_plot(paste("plot/", "accuracy.pdf", sep=""), 
#           plot_grid(p1, p2, nrow=1, align = "vh", axis="btlr"))

p3 = ggplot(accuracy_sum[accuracy_sum$trait == "IRR_GrainYield" & 
                           accuracy_sum$adjustment == "predictive ability", ], 
            aes(x=prediction_method)) +
  geom_point(aes(y=accuracy.mean, colour=prediction_method), size=0.75) + 
  geom_errorbar(aes(ymin=accuracy.mean-accuracy.se, ymax=accuracy.mean+accuracy.se, 
                    colour=prediction_method), width=0.1) + 
  facet_wrap(~adjustment) +
  ylim(0, 1) + 
  ylab(expression("r"["y"])) + 
  scale_colour_manual(values=c("blue", "gold3", "red")) + 
  theme_minimal_grid(font_size=7) +
  theme(axis.title.x = element_blank(),
        #     axis.text.x = element_blank(),
        #     axis.ticks.x = element_blank(),
        legend.position="none") + 
  ggtitle(paste("trait: wheat grain yield"))

p4 = ggplot(accuracy_sum[accuracy_sum$trait == "IRR_GrainYield" & 
                           accuracy_sum$adjustment == "prediction accuracy", ],
            aes(x=prediction_method)) +
  geom_point(aes(y=accuracy.mean, colour=prediction_method), size=0.75) + 
  geom_errorbar(aes(ymin=accuracy.mean-accuracy.se, ymax=accuracy.mean+accuracy.se, 
                    colour=prediction_method), width=0.1) + 
  facet_wrap(~adjustment) +
  ylim(0, 1) + 
  ylab("r") + 
  scale_colour_manual(values=c("blue", "gold3", "red")) + 
  theme_minimal_grid(font_size=7) +
  theme(axis.title.x = element_blank(),
        #     axis.text.x = element_blank(),
        #     axis.ticks.x = element_blank(),
        legend.position="none") 

# save_plot(paste("plot/","accuracy.pdf", sep=""), 
#           plot_grid(p3, p4, nrow=1, align = "vh", axis="btlr"))

p5 = ggplot(accuracy_sum[accuracy_sum$data == "Winn" & 
                           accuracy_sum$trait == "GrainYield" & 
                           accuracy_sum$adjustment == "predictive ability", ], 
            aes(x=prediction_method)) +
  geom_point(aes(y=accuracy.mean, colour=prediction_method), size=0.75) + 
  geom_errorbar(aes(ymin=accuracy.mean-accuracy.se, ymax=accuracy.mean+accuracy.se, 
                    colour=prediction_method), width=0.1) + 
  facet_wrap(~adjustment) +
  ylim(0, 1) + 
  ylab(expression("r"["y"])) + 
  scale_colour_manual(values=c("blue", "gold3", "red")) + 
  theme_minimal_grid(font_size=7) +
  theme(axis.title.x = element_blank(),
        #     axis.text.x = element_blank(),
        #     axis.ticks.x = element_blank(),
        legend.position="none") + 
  ggtitle(paste("trait: wheat grain yield"))

p6 = ggplot(accuracy_sum[accuracy_sum$data == "Winn" & 
                           accuracy_sum$trait == "GrainYield" & 
                           accuracy_sum$adjustment == "prediction accuracy", ],
            aes(x=prediction_method)) +
  geom_point(aes(y=accuracy.mean, colour=prediction_method), size=0.75) + 
  geom_errorbar(aes(ymin=accuracy.mean-accuracy.se, ymax=accuracy.mean+accuracy.se, 
                    colour=prediction_method), width=0.1) + 
  facet_wrap(~adjustment) +
  ylim(0, 1) + 
  ylab("r") + 
  scale_colour_manual(values=c("blue", "gold3", "red")) + 
  theme_minimal_grid(font_size=7) +
  theme(axis.title.x = element_blank(),
        #     axis.text.x = element_blank(),
        #     axis.ticks.x = element_blank(),
        legend.position="none") 

# save_plot(paste("plot/","accuracy.pdf", sep=""), 
#           plot_grid(p5, p6, nrow=1, align = "vh", axis="btlr"))

p7 = ggplot(accuracy_sum[accuracy_sum$data == "Krause" & 
                           accuracy_sum$trait == "GrainYield" & 
                           accuracy_sum$adjustment == "predictive ability", ], 
            aes(x=prediction_method)) +
  geom_point(aes(y=accuracy.mean, colour=prediction_method), size=0.75) + 
  geom_errorbar(aes(ymin=accuracy.mean-accuracy.se, ymax=accuracy.mean+accuracy.se, 
                    colour=prediction_method), width=0.1) + 
  facet_wrap(~adjustment) +
  ylim(0, 1) + 
  ylab(expression("r"["y"])) + 
  scale_colour_manual(values=c("blue", "gold3")) + 
  theme_minimal_grid(font_size=7) +
  theme(axis.title.x = element_blank(),
        #     axis.text.x = element_blank(),
        #     axis.ticks.x = element_blank(),
        legend.position="none") + 
  ggtitle(paste("trait: wheat grain yield"))

p8 = ggplot(accuracy_sum[accuracy_sum$data == "Krause" & 
                           accuracy_sum$trait == "GrainYield" & 
                           accuracy_sum$adjustment == "prediction accuracy", ],
            aes(x=prediction_method)) +
  geom_point(aes(y=accuracy.mean, colour=prediction_method), size=0.75) + 
  geom_errorbar(aes(ymin=accuracy.mean-accuracy.se, ymax=accuracy.mean+accuracy.se, 
                    colour=prediction_method), width=0.1) + 
  facet_wrap(~adjustment) +
  ylim(0, 1) + 
  ylab("r") + 
  scale_colour_manual(values=c("blue", "gold3")) + 
  theme_minimal_grid(font_size=7) +
  theme(axis.title.x = element_blank(),
        #     axis.text.x = element_blank(),
        #     axis.ticks.x = element_blank(),
        legend.position="none") 

save_plot(paste("plot/","accuracy.pdf", sep=""), 
          plot_grid(plot_grid(p1, p2, nrow=1, align = "vh", axis="btlr"), 
                    plot_grid(p3, p4, nrow=1, align = "vh", axis="btlr"), 
                    plot_grid(p5, p6, nrow=1, align = "vh", axis="btlr"), 
                    plot_grid(p7, p8, nrow=1, align = "vh", axis="btlr"), 
                    ncol=1, align="vh", axis="btlr", 
                    labels = "auto"),
          base_width = 5.5, base_height = 7)


