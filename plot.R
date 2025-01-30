setwd("~/breeder_equation_project/")
library(ggplot2)
library(cowplot)

data1 <- readRDS("estimate_accuracy/accuracy_table_long.rds")
data2 <- readRDS("estimate_accuracy2/accuracy_table_long.rds")
data5 <- readRDS("estimate_accuracy5/accuracy_table_long.rds")
data4 <- readRDS("estimate_accuracy4/accuracy_table_long.rds")

data1 <- cbind(`Case Study` = "cite{rincent2018phenomic}", 
               Crop = "Poplar", 
               data1)
data2 <- cbind(`Case Study` = "cite{rincent2018phenomic}", 
               Crop = "Wheat", 
               data2)
data5 <- cbind(`Case Study` = "cite{winn2023phenomic}", 
               Crop = "Wheat", 
               trait = "GrainYield", 
               data5)
# colnames(data4)[1] <- "jth_fold"
data4 <- cbind(`Case Study` = "cite{krause_hyperspectral_2019}", 
               Crop = "Wheat", 
               data4)
colnames(data4)[3] <- "trait"

accuracy <- rbind(rbind(rbind(data1, data2), data5), data4)
colnames(accuracy)[3] <- "Trait"
accuracy[accuracy$accuracy > 1, ]
accuracy[accuracy$accuracy > 1, "accuracy"] <- 1
accuracy[accuracy$accuracy == 1, ]



accuracy_summary <- aggregate(. ~ `Case Study` + Crop + Trait + adjustment + `prediction method`, 
                              data=accuracy[, !colnames(accuracy) %in% c("type", "jth_fold")], 
                              FUN=function(x){c(mean=mean(x), se=sd(x)/sqrt(length(x)))})
saveRDS(accuracy_summary, "plot/accuracy_summary.rds")
accuracy_summary <- readRDS("plot/accuracy_summary.rds")

accuracy_sum <- data.frame(accuracy_summary[[1]])
for (i in 2:length(accuracy_summary)){
  accuracy_sum = cbind(accuracy_sum, accuracy_summary[[i]])
}
colnames(accuracy_sum) <- c("Case Study", "Crop", "Trait", "adjustment", "prediction_method", "accuracy.mean", 
                            "accuracy.se", "h2.mean", "h2.se")



p1 = ggplot(accuracy_sum[accuracy_sum$Trait == "CIRC.ORL" & 
                           accuracy_sum$adjustment == "predictive ability", ], 
            aes(x=prediction_method)) +
  geom_point(aes(y=accuracy.mean, colour=prediction_method), size=0.75) + 
  geom_errorbar(aes(ymin=accuracy.mean-accuracy.se, ymax=accuracy.mean+accuracy.se, 
                    colour=prediction_method), width=0.1, linewidth=0.35) + 
  facet_wrap(~adjustment) +
  ylim(0, 1) + 
  ylab(expression("r"["y"])) + 
  scale_colour_manual(values=c("blue", "gold3")) + 
  theme_minimal_grid(font_size=7) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, vjust = 0.6),
    #     axis.ticks.x = element_blank(),
        legend.position="none") + 
  ggtitle(paste("Case 1: poplar circumference"))

p2 = ggplot(accuracy_sum[accuracy_sum$Trait == "CIRC.ORL" & 
                           accuracy_sum$adjustment == "prediction accuracy", ],
            aes(x=prediction_method)) +
  geom_point(aes(y=accuracy.mean, colour=prediction_method), size=0.75) + 
  geom_errorbar(aes(ymin=accuracy.mean-accuracy.se, ymax=accuracy.mean+accuracy.se, 
                    colour=prediction_method), width=0.1, linewidth=0.35) + 
  facet_wrap(~adjustment) +
  ylim(0, 1) + 
  ylab("r") + 
  scale_colour_manual(values=c("blue", "gold3")) + 
  theme_minimal_grid(font_size=7) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, vjust = 0.6),
    #     axis.ticks.x = element_blank(),
        legend.position="none") 

# save_plot(paste("plot/", "accuracy.pdf", sep=""), 
#           plot_grid(p1, p2, nrow=1, align = "vh", axis="btlr"))

p3 = ggplot(accuracy_sum[accuracy_sum$Trait == "IRR_GrainYield" & 
                           accuracy_sum$adjustment == "predictive ability", ], 
            aes(x=prediction_method)) +
  geom_point(aes(y=accuracy.mean, colour=prediction_method), size=0.75) + 
  geom_errorbar(aes(ymin=accuracy.mean-accuracy.se, ymax=accuracy.mean+accuracy.se, 
                    colour=prediction_method), width=0.1, linewidth=0.35) + 
  facet_wrap(~adjustment) +
  ylim(0, 1) + 
  ylab(expression("r"["y"])) + 
  scale_colour_manual(values=c("blue", "gold3", "red")) + 
  theme_minimal_grid(font_size=7) +
  theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 60, vjust = 0.6),
        #     axis.ticks.x = element_blank(),
        legend.position="none") + 
  ggtitle(paste("Case 1: wheat grain yield"))

p4 = ggplot(accuracy_sum[accuracy_sum$Trait == "IRR_GrainYield" & 
                           accuracy_sum$adjustment == "prediction accuracy", ],
            aes(x=prediction_method)) +
  geom_point(aes(y=accuracy.mean, colour=prediction_method), size=0.75) + 
  geom_errorbar(aes(ymin=accuracy.mean-accuracy.se, ymax=accuracy.mean+accuracy.se, 
                    colour=prediction_method), width=0.1, linewidth=0.35) + 
  facet_wrap(~adjustment) +
  ylim(0, 1) + 
  ylab("r") + 
  scale_colour_manual(values=c("blue", "gold3", "red")) + 
  theme_minimal_grid(font_size=7) +
  theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 60, vjust = 0.6),
        #     axis.ticks.x = element_blank(),
        legend.position="none") 

# save_plot(paste("plot/","accuracy.pdf", sep=""), 
#           plot_grid(p3, p4, nrow=1, align = "vh", axis="btlr"))

p5 = ggplot(accuracy_sum[accuracy_sum$`Case Study` == "cite{winn2023phenomic}" & 
                           accuracy_sum$Trait == "GrainYield" & 
                           accuracy_sum$adjustment == "predictive ability", ], 
            aes(x=prediction_method)) +
  geom_point(aes(y=accuracy.mean, colour=prediction_method), size=0.75) + 
  geom_errorbar(aes(ymin=accuracy.mean-accuracy.se, ymax=accuracy.mean+accuracy.se, 
                    colour=prediction_method), width=0.1, linewidth=0.35) + 
  facet_wrap(~adjustment) +
  ylim(0, 1) + 
  ylab(expression("r"["y"])) + 
  scale_colour_manual(values=c("blue", "gold3", "red")) + 
  theme_minimal_grid(font_size=7) +
  theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 60, vjust = 0.6),
        #     axis.ticks.x = element_blank(),
        legend.position="none") + 
  ggtitle(paste("Case 2: wheat grain yield"))

p6 = ggplot(accuracy_sum[accuracy_sum$`Case Study` == "cite{winn2023phenomic}" & 
                           accuracy_sum$Trait == "GrainYield" & 
                           accuracy_sum$adjustment == "prediction accuracy", ],
            aes(x=prediction_method)) +
  geom_point(aes(y=accuracy.mean, colour=prediction_method), size=0.75) + 
  geom_errorbar(aes(ymin=accuracy.mean-accuracy.se, ymax=accuracy.mean+accuracy.se, 
                    colour=prediction_method), width=0.1, linewidth=0.35) + 
  facet_wrap(~adjustment) +
  ylim(0, 1) + 
  ylab("r") + 
  scale_colour_manual(values=c("blue", "gold3", "red")) + 
  theme_minimal_grid(font_size=7) +
  theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 60, vjust = 0.6),
        #     axis.ticks.x = element_blank(),
        legend.position="none") 

# save_plot(paste("plot/","accuracy.pdf", sep=""), 
#           plot_grid(p5, p6, nrow=1, align = "vh", axis="btlr"))

p7 = ggplot(accuracy_sum[accuracy_sum$`Case Study` == "cite{krause_hyperspectral_2019}" & 
                           accuracy_sum$Trait == "2013-14_Moderate Drought" & 
                           accuracy_sum$adjustment == "predictive ability", ], 
            aes(x=prediction_method)) +
  geom_point(aes(y=accuracy.mean, colour=prediction_method), size=0.75) + 
  geom_errorbar(aes(ymin=accuracy.mean-accuracy.se, ymax=accuracy.mean+accuracy.se, 
                    colour=prediction_method), width=0.1, linewidth=0.35) + 
  facet_wrap(~adjustment) +
  ylim(0, 1) + 
  ylab(expression("r"["y"])) + 
  scale_colour_manual(values=c("blue", "gold3")) + 
  theme_minimal_grid(font_size=7) +
  theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 60, vjust = 0.6),
        #     axis.ticks.x = element_blank(),
        legend.position="none") + 
  ggtitle(paste("Case 3: wheat grain yield"))

p8 = ggplot(accuracy_sum[accuracy_sum$`Case Study` == "cite{krause_hyperspectral_2019}" & 
                           accuracy_sum$Trait == "2013-14_Moderate Drought" & 
                           accuracy_sum$adjustment == "prediction accuracy", ],
            aes(x=prediction_method)) +
  geom_point(aes(y=accuracy.mean, colour=prediction_method), size=0.75) + 
  geom_errorbar(aes(ymin=accuracy.mean-accuracy.se, ymax=accuracy.mean+accuracy.se, 
                    colour=prediction_method), width=0.1, linewidth=0.35) + 
  facet_wrap(~adjustment) +
  ylim(0, 1) + 
  ylab("r") + 
  scale_colour_manual(values=c("blue", "gold3")) + 
  theme_minimal_grid(font_size=7) +
  theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 60, vjust = 0.6),
        #     axis.ticks.x = element_blank(),
        legend.position="none") 

save_plot(paste("plot/","accuracy.pdf", sep=""), 
          plot_grid(plot_grid(p1, p2, p3, p4, nrow=1, align = "vh", axis="btlr", labels=c("a", NA, "b", NA)), 
                    # plot_grid(p3, p4, nrow=1, align = "vh", axis="btlr"), 
                    plot_grid(p5, p6, p7, p8, nrow=1, align = "vh", axis="btlr", labels=c("c", NA, "d", NA)), 
                    # plot_grid(p7, p8, nrow=1, align = "vh", axis="btlr"), 
                    ncol=1, align="vh", axis="btlr"
                    # labels = "auto"
                    ),
          base_width = 5.5, base_height = 3.5)



accuracy_sum[, 6:9] <- sapply(accuracy_sum[, 6:9], FUN=function(x){round(x, digits = 3)})
accuracy_sum$`Case Study` <- factor(accuracy_sum$`Case Study`, levels=as.factor(c(
                                                                  "cite{rincent2018phenomic}",
                                                                  "cite{winn2023phenomic}", 
                                                                  "cite{krause_hyperspectral_2019}")))
accuracy_sum$Trait <- factor(accuracy_sum$Trait, levels=as.factor(
  c("HT.ORL", "BF.ORL", "BF.SAV", "BS.ORL", "BS.SAV", "CIRC.ORL", "CIRC.SAV", "RUST.ORL", "DRY_GrainYield", 
    "DRY_HeadingDate", "IRR_GrainYield", "IRR_HeadingDate", "GrainYield", 
    "2013-14_Heat", "2013-14_Moderate Drought", "2013-14_Optimal Bed", "2013-14_Optimal Flat", "2013-14_Severe Drought", 
    "2014-15_Moderate Drought", "2014-15_Optimal Bed", "2014-15_Optimal Flat", "2014-15_Severe Drought", 
    "2015-16_Heat", "2015-16_Moderate Drought", "2015-16_Optimal Bed", "2015-16_Optimal Flat", "2015-16_Severe Drought", 
    "2016-17_Heat", "2016-17_Moderate Drought", "2016-17_Optimal Bed", "2016-17_Optimal Flat", "2016-17_Severe Drought")
))
accuracy_sum$adjustment <- factor(accuracy_sum$adjustment, 
                                  levels=as.factor(c("predictive ability", "prediction accuracy")))
accuracy_sum$prediction_method <- factor(accuracy_sum$prediction_method, levels=as.factor(
  c("genomic", "phenomic", "phenomic-grain", "phenomic-leaf", "phenomic-fixed", "phenomic-mixed")
))
accuracy_sum <- accuracy_sum[order(accuracy_sum$`Case Study`, 
                                   accuracy_sum$Trait, 
                                   accuracy_sum$prediction_method, 
                                   accuracy_sum$adjustment), ]

saveRDS(accuracy_sum, "plot/accuracy_sum.rds")
accuracy_sum <- readRDS("plot/accuracy_sum.rds")
accuracy_sum[, 1:5] <- sapply(accuracy_sum[, 1:5], FUN=function(x){as.character(x)})



accuracy_wide <- data.frame(`Case Study` = rep(NA, 69),
                            Crop = rep(NA, 69), 
                            Trait = rep(NA, 69), 
                            Prediction = rep(NA, 69), 
                            r = rep(NA, 69), 
                            r_y = rep(NA, 69), 
                            h2 = rep(NA, 69))

for (i in 1:69){
  accuracy_wide[i, 1:4] = accuracy_sum[(i-1)*2+1, c(1, 2, 3, 5)] 
  accuracy_wide[i, 5] = paste(accuracy_sum[(i-1)*2+1, 6], "(", accuracy_sum[(i-1)*2+1, 7], ")", sep="")
  accuracy_wide[i, 6] = paste(accuracy_sum[(i-1)*2+2, 6], "(", accuracy_sum[(i-1)*2+2, 7], ")", sep="")
  accuracy_wide[i, 7] = paste(accuracy_sum[(i-1)*2+1, 8], "(", accuracy_sum[(i-1)*2+1, 9], ")", sep="")
}

accuracy_wide$Trait[32:69] <- sub("-.*_", "", accuracy_wide$Trait[32:69])
accuracy_wide$Trait[32:69] <- sub("Moderate Drought", "ModDrt", accuracy_wide$Trait[32:69])
accuracy_wide$Trait[32:69] <- sub("Optimal Bed", "OptBed", accuracy_wide$Trait[32:69])
accuracy_wide$Trait[32:69] <- sub("Optimal Flat", "OptFlt", accuracy_wide$Trait[32:69])
accuracy_wide$Trait[32:69] <- sub("Severe Drought", "SevDrt", accuracy_wide$Trait[32:69])
accuracy_wide$Trait[32:69] <- paste(accuracy_wide$Trait[32:69], "_GrainYield", sep="")
accuracy_wide$Trait[1:31] <- c("Orl_Height", "Orl_Height", 
                         "Orl_BudFlush", "Orl_BudFlush", 
                         "Sav_BudFlush", "Sav_BudFlush", 
                         "Orl_BudSet", "Orl_BudSet", 
                         "Sav_BudSet", "Sav_BudSet", 
                         "Orl_Circumference", "Orl_Circumference", 
                         "Sav_Circumference", "Sav_Circumference", 
                         "Orl_Rust", "Orl_Rust", 
                         "Dry_GrainYield", "Dry_GrainYield", "Dry_GrainYield", 
                         "Dry_HeadingDate", "Dry_HeadingDate", "Dry_HeadingDate", 
                         "Irr_GrainYield", "Irr_GrainYield", "Irr_GrainYield", 
                         "Irr_HeadingDate", "Irr_HeadingDate", "Irr_HeadingDate", 
                         "GrainYield", "GrainYield", "GrainYield")
colnames(accuracy_wide)[1] <- "Case Study"
# accuracy_wide$`Case Study` <- sub("\\\\", "", accuracy_wide$`Case Study`)

saveRDS(accuracy_wide, "plot/accuracy_wide.rds")
accuracy_wide <- readRDS("plot/accuracy_wide.rds")
write.csv(accuracy_wide, "plot/accuracy_wide.csv", row.names = F, quote = F)





p1 = ggplot(accuracy_sum[accuracy_sum$trait == "CIRC.ORL" & 
                           accuracy_sum$adjustment == "predictive ability", ], 
            aes(x=prediction_method)) +
  geom_point(aes(y=accuracy.mean, colour=prediction_method), size=0.75) + 
  geom_errorbar(aes(ymin=accuracy.mean-accuracy.se, ymax=accuracy.mean+accuracy.se, 
                    colour=prediction_method), width=0.1) + 
  # facet_wrap(~adjustment) +
  ylim(0, 1) + 
  ylab("accuracy") + 
  scale_colour_manual(values=c("blue", "gold3")) + 
  theme_minimal_grid(font_size=9) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, vjust = 0.6),
    #     axis.ticks.x = element_blank(),
        legend.position="none") + 
  ggtitle(paste("Case 1: poplar circumference"))

# p2 = ggplot(accuracy_sum[accuracy_sum$trait == "CIRC.ORL" & 
#                            accuracy_sum$adjustment == "prediction accuracy", ],
#             aes(x=prediction_method)) +
#   geom_point(aes(y=accuracy.mean, colour=prediction_method), size=0.75) + 
#   geom_errorbar(aes(ymin=accuracy.mean-accuracy.se, ymax=accuracy.mean+accuracy.se, 
#                     colour=prediction_method), width=0.1) + 
#   facet_wrap(~adjustment) +
#   ylim(0, 1) + 
#   ylab("r") + 
#   scale_colour_manual(values=c("blue", "gold3")) + 
#   theme_minimal_grid(font_size=7) +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_text(angle = 60, vjust = 0.6),
#     #     axis.ticks.x = element_blank(),
#         legend.position="none") 

p3 = ggplot(accuracy_sum[accuracy_sum$trait == "IRR_GrainYield" & 
                           accuracy_sum$adjustment == "predictive ability" &
                           accuracy_sum$prediction_method != "phenomic-leaf", ], 
            aes(x=prediction_method)) +
  geom_point(aes(y=accuracy.mean, colour=prediction_method), size=0.75) + 
  geom_errorbar(aes(ymin=accuracy.mean-accuracy.se, ymax=accuracy.mean+accuracy.se, 
                    colour=prediction_method), width=0.1) + 
  # facet_wrap(~adjustment) +
  ylim(0, 1) + 
  # ylab("accuracy") + 
  scale_x_discrete(labels= c("genomic", "phenomic")) + 
  scale_colour_manual(values=c("blue", "gold3")) + 
  theme_minimal_grid(font_size=9) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, vjust = 0.6),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        #     axis.ticks.x = element_blank(),
        legend.position="none") + 
  ggtitle(paste("Case 2: wheat grain yield"))

save_plot(paste("plot/", "p1_p3.pdf", sep=""),
          plot_grid(p1, p3, nrow=1, align = "vh", axis="btlr", labels=c("a", "b")))



# p1 = ggplot(accuracy_sum[accuracy_sum$trait == "CIRC.ORL" & 
#                            accuracy_sum$adjustment == "predictive ability", ], 
#             aes(x=prediction_method)) +
#   geom_point(aes(y=accuracy.mean, colour=prediction_method), size=0.75) + 
#   geom_errorbar(aes(ymin=accuracy.mean-accuracy.se, ymax=accuracy.mean+accuracy.se, 
#                     colour=prediction_method), width=0.1) + 
#   facet_wrap(~adjustment) +
#   ylim(0, 1) + 
#   ylab(expression("r"["y"])) + 
#   scale_colour_manual(values=c("blue", "gold3")) + 
#   theme_minimal_grid(font_size=9) +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_text(angle = 60, vjust = 0.6),
#         #     axis.ticks.x = element_blank(),
#         legend.position="none") + 
#   ggtitle(paste("Case 1: poplar circumference"))
# 
# p2 = ggplot(accuracy_sum[accuracy_sum$trait == "CIRC.ORL" & 
#                            accuracy_sum$adjustment == "prediction accuracy", ],
#             aes(x=prediction_method)) +
#   geom_point(aes(y=accuracy.mean, colour=prediction_method), size=0.75) + 
#   geom_errorbar(aes(ymin=accuracy.mean-accuracy.se, ymax=accuracy.mean+accuracy.se, 
#                     colour=prediction_method), width=0.1) + 
#   facet_wrap(~adjustment) +
#   ylim(0, 1) + 
#   ylab("r") + 
#   scale_colour_manual(values=c("blue", "gold3")) + 
#   theme_minimal_grid(font_size=9) +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_text(angle = 60, vjust = 0.6),
#         #     axis.ticks.x = element_blank(),
#         legend.position="none") 
# 
# # save_plot(paste("plot/", "accuracy.pdf", sep=""), 
# #           plot_grid(p1, p2, nrow=1, align = "vh", axis="btlr"))
# 
# p3 = ggplot(accuracy_sum[accuracy_sum$trait == "IRR_GrainYield" & 
#                            accuracy_sum$adjustment == "predictive ability" &
#                            accuracy_sum$prediction_method != "phenomic-leaf", ], 
#             aes(x=prediction_method)) +
#   geom_point(aes(y=accuracy.mean, colour=prediction_method), size=0.75) + 
#   geom_errorbar(aes(ymin=accuracy.mean-accuracy.se, ymax=accuracy.mean+accuracy.se, 
#                     colour=prediction_method), width=0.1) + 
#   facet_wrap(~adjustment) +
#   ylim(0, 1) + 
#   ylab(expression("r"["y"])) + 
#   scale_x_discrete(labels= c("genomic", "phenomic")) + 
#   scale_colour_manual(values=c("blue", "gold3")) + 
#   theme_minimal_grid(font_size=9) +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_text(angle = 60, vjust = 0.6),
#         #     axis.ticks.x = element_blank(),
#         legend.position="none") + 
#   ggtitle(paste("Case 2: wheat grain yield"))
# 
# p4 = ggplot(accuracy_sum[accuracy_sum$trait == "IRR_GrainYield" & 
#                            accuracy_sum$adjustment == "prediction accuracy" &
#                            accuracy_sum$prediction_method != "phenomic-leaf", ],
#             aes(x=prediction_method)) +
#   geom_point(aes(y=accuracy.mean, colour=prediction_method), size=0.75) + 
#   geom_errorbar(aes(ymin=accuracy.mean-accuracy.se, ymax=accuracy.mean+accuracy.se, 
#                     colour=prediction_method), width=0.1) + 
#   facet_wrap(~adjustment) +
#   ylim(0, 1) + 
#   ylab("r") + 
#   scale_x_discrete(labels= c("genomic", "phenomic")) + 
#   scale_colour_manual(values=c("blue", "gold3")) + 
#   theme_minimal_grid(font_size=9) +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_text(angle = 60, vjust = 0.6),
#         #     axis.ticks.x = element_blank(),
#         legend.position="none") 
# 
# # save_plot(paste("plot/","accuracy.pdf", sep=""), 
# #           plot_grid(p3, p4, nrow=1, align = "vh", axis="btlr"))
# 
# save_plot(paste("plot/","p1_p2_p3_p4.pdf", sep=""), 
#           plot_grid(plot_grid(p1, p2, p3, p4, nrow=1, align = "vh", axis="btlr", 
#                               labels=c("a", NA, "b", NA))), 
#           base_asp = 3.236)



accuracy_sum <- readRDS("plot/accuracy_sum.rds")
accuracy_sum[, 1:5] <- sapply(accuracy_sum[, 1:5], FUN=function(x){as.character(x)})



dim(accuracy_sum)
# [1] 138   9
dim(accuracy_sum[accuracy_sum$`Case Study` == "cite{rincent2018phenomic}" &
                   accuracy_sum$Crop == "Poplar", ])
# [1] 32  9
dim(accuracy_sum[accuracy_sum$`Case Study` == "cite{rincent2018phenomic}" &
                   accuracy_sum$Crop == "Wheat", ])
# [1] 24  9
dim(accuracy_sum[accuracy_sum$`Case Study` == "cite{winn2023phenomic}", ])
# [1] 6 9
dim(accuracy_sum[accuracy_sum$`Case Study` == "cite{krause_hyperspectral_2019}", ])
# [1] 76  9
138-32-76
# 30

accuracy_diff <- data.frame(`Case Study` = rep(NA, 37), 
                            Crop = NA, 
                            Trait = NA, 
                            # adjustment = NA, 
                            Difference_Type = NA, 
                            Difference_ry = NA, 
                            Difference_r = NA)
for (i in 1:8){
  # print(paste("i: ", i, sep=""))
  # print((i-1)*4+1)
  accuracy_diff[i, 1:3] = accuracy_sum[(i-1)*4+1, 1:3]
  accuracy_diff[i, 4] = "Phenomic - Genomic"
  accuracy_diff[i, 5] = accuracy_sum[(i-1)*4+2, 6] - accuracy_sum[(i-1)*4+1, 6]
  accuracy_diff[i, 6] = accuracy_sum[(i-1)*4+4, 6] - accuracy_sum[(i-1)*4+3, 6]
}
for (i in seq(9, 15, 2)){
  # print(paste("i: ", i, sep=""))
  # print(((i-9)*6/2+1+32) : ((i-9)*6/2+6+32))
  accuracy_diff[i:(i+1), 1:3] = accuracy_sum[((i-9)*6/2+1+32) : ((i-9)*6/2+2+32), 1:3]
  accuracy_diff[i, 4] = "Phenomic_Grain - Genomic"
  accuracy_diff[i+1, 4] = "Phenomic_Leaf - Genomic"
  accuracy_diff[i, 5] = accuracy_sum[(i-9)*6/2+2+32, 6] - accuracy_sum[(i-9)*6/2+1+32, 6]
  accuracy_diff[i+1, 5] = accuracy_sum[(i-9)*6/2+3+32, 6] - accuracy_sum[(i-9)*6/2+1+32, 6]
  accuracy_diff[i, 6] = accuracy_sum[(i-9)*6/2+5+32, 6] - accuracy_sum[(i-9)*6/2+4+32, 6]
  accuracy_diff[i+1, 6] = accuracy_sum[(i-9)*6/2+6+32, 6] - accuracy_sum[(i-9)*6/2+4+32, 6]
}
for (i in 17){
  # print(paste("i: ", i, sep=""))
  # print(((i-17)*6/2+1+56) : ((i-17)*6/2+6+56))
  accuracy_diff[i:(i+1), 1:3] = accuracy_sum[((i-17)*6/2+1+56) : ((i-17)*6/2+2+56), 1:3]
  accuracy_diff[i, 4] = "Phenomic_Fixed - Genomic"
  accuracy_diff[i+1, 4] = "Phenomic_Random - Genomic"
  accuracy_diff[i, 5] = accuracy_sum[(i-17)*6/2+2+56, 6] - accuracy_sum[(i-17)*6/2+1+56, 6]
  accuracy_diff[i+1, 5] = accuracy_sum[(i-17)*6/2+3+56, 6] - accuracy_sum[(i-17)*6/2+1+56, 6]
  accuracy_diff[i, 6] = accuracy_sum[(i-17)*6/2+5+56, 6] - accuracy_sum[(i-17)*6/2+4+56, 6]
  accuracy_diff[i+1, 6] = accuracy_sum[(i-17)*6/2+6+56, 6] - accuracy_sum[(i-17)*6/2+4+56, 6]
}
for (i in 19:37){
  # print(paste("i: ", i, sep=""))
  # print((i-19)*4+1+62)
  accuracy_diff[i, 1:3] = accuracy_sum[(i-19)*4+1+62, 1:3]
  accuracy_diff[i, 4] = "Phenomic - Genomic"
  accuracy_diff[i, 5] = accuracy_sum[(i-19)*4+2+62, 6] - accuracy_sum[(i-19)*4+1+62, 6]
  accuracy_diff[i, 6] = accuracy_sum[(i-19)*4+4+62, 6] - accuracy_sum[(i-19)*4+3+62, 6]
}

saveRDS(accuracy_diff, "plot/accuracy_diff.rds")
accuracy_diff <- readRDS("plot/accuracy_diff.rds")



p9 <- ggplot(accuracy_diff, aes(x=Difference_ry, y=Difference_r)) + 
  geom_point(size=0.75, colour=prediction_method) + 
  geom_hline(yintercept=0) + 
  geom_vline(xintercept=0) + 
  xlim(-1, 1) + 
  ylim(-1, 1) + 
  xlab(expression(r["y, phenomic"]-r["y, genomic"])) + 
  ylab(expression(r[phenomic]-r[genomic])) + 
  theme_minimal_grid(font_size=7)
# save_plot(paste("plot/", "p9.pdf", sep=""),
#           plot_grid(p9))



accuracy_sum <- readRDS("plot/accuracy_sum.rds")
accuracy_sum[, 1:5] <- sapply(accuracy_sum[, 1:5], FUN=function(x){as.character(x)})

accuracy_diff2 <- data.frame(`Case Study` = rep(NA, 69),
                             Crop = rep(NA, 69), 
                             Trait = rep(NA, 69), 
                             Prediction = rep(NA, 69), 
                             Difference = NA)

for (i in 1:69){
  accuracy_diff2[i, 1:4] = accuracy_sum[(i-1)*2+1, c(1, 2, 3, 5)] 
  accuracy_diff2[i, 5] = accuracy_sum[(i-1)*2+2, 6] - accuracy_sum[(i-1)*2+1, 6]
  #   paste(accuracy_sum[(i-1)*2+1, 6], "(", accuracy_sum[(i-1)*2+1, 7], ")", sep="")
  # accuracy_wide[i, 6] = paste(accuracy_sum[(i-1)*2+2, 6], "(", accuracy_sum[(i-1)*2+2, 7], ")", sep="")
  # accuracy_wide[i, 7] = paste(accuracy_sum[(i-1)*2+1, 8], "(", accuracy_sum[(i-1)*2+1, 9], ")", sep="")
}

accuracy_diff2[accuracy_diff2$Prediction %in% c("phenomic-grain", "phenomic-leaf", 
                                                "phenomic-fixed", "phenomic-mixed"), 
               "Prediction"] <- "phenomic"

p10 <- ggplot(accuracy_diff2, aes(x=Prediction, y=Difference)) + 
  geom_point(size=0.75) + 
  geom_hline(yintercept=0) + 
  xlab("prediction method") + 
  ylab(expression(r-r[y])) + 
  theme_minimal_grid(font_size=7)
save_plot(paste("plot/", "p9_p10.pdf", sep=""),
          plot_grid(p9, p10, nrow=1, align = "vh", axis="btlr", labels=c("a", "b")), 
          base_width = 5.5, base_height = 1.75)







 

 
















