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
accuracy_summary <- readRDS("plot/accuracy_summary.rds")

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
        axis.text.x = element_text(angle = 60, vjust = 0.6),
    #     axis.ticks.x = element_blank(),
        legend.position="none") + 
  ggtitle(paste("Case 1: poplar circumference"))

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
        axis.text.x = element_text(angle = 60, vjust = 0.6),
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
            axis.text.x = element_text(angle = 60, vjust = 0.6),
        #     axis.ticks.x = element_blank(),
        legend.position="none") + 
  ggtitle(paste("Case 1: wheat grain yield"))

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
            axis.text.x = element_text(angle = 60, vjust = 0.6),
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
            axis.text.x = element_text(angle = 60, vjust = 0.6),
        #     axis.ticks.x = element_blank(),
        legend.position="none") + 
  ggtitle(paste("Case 2: wheat grain yield"))

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
            axis.text.x = element_text(angle = 60, vjust = 0.6),
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
            axis.text.x = element_text(angle = 60, vjust = 0.6),
        #     axis.ticks.x = element_blank(),
        legend.position="none") + 
  ggtitle(paste("Case 3: wheat grain yield"))

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



accuracy_sum[, 5:8] <- sapply(accuracy_sum[, 5:8], FUN=function(x){round(x, digits = 3)})
accuracy_sum$data <- factor(accuracy_sum$data, levels=as.factor(c("poplar_Rincent", "wheat_Rincent", "Winn", "Krause")))
accuracy_sum$trait <- factor(accuracy_sum$trait, levels=as.factor(
  c("HT.ORL", "BF.ORL", "BF.SAV", "BS.ORL", "BS.SAV", "CIRC.ORL", "CIRC.SAV", "RUST.ORL", "DRY_GrainYield", 
    "DRY_HeadingDate", "IRR_GrainYield", "IRR_HeadingDate", "GrainYield")
))
accuracy_sum$adjustment <- factor(accuracy_sum$adjustment, 
                                  levels=as.factor(c("predictive ability", "prediction accuracy")))
accuracy_sum$prediction_method <- factor(accuracy_sum$prediction_method, levels=as.factor(
  c("genomic", "phenomic", "phenomic-grain", "phenomic-leaf", "phenomic-fixed", "phenomic-mixed")
))
accuracy_sum <- accuracy_sum[order(accuracy_sum$data, 
                                   accuracy_sum$trait, 
                                   accuracy_sum$prediction_method, 
                                   accuracy_sum$adjustment), ]
accuracy_sum[, 1:4] <- sapply(accuracy_sum[, 1:4], FUN=function(x){as.character(x)})



accuracy_wide <- data.frame(Case_Study = rep(NA, 33),
                            Trait = rep(NA, 33), 
                            Prediction = rep(NA, 33), 
                            r = rep(NA, 33), 
                            r_y = rep(NA, 33), 
                            h2 = rep(NA, 33))

for (i in 1:33){
  accuracy_wide[i, 1:3] = accuracy_sum[(i-1)*2+1, c(1, 2, 4)] 
  accuracy_wide[i, 4] = paste(accuracy_sum[(i-1)*2+1, 5], "(", accuracy_sum[(i-1)*2+1, 6], ")", sep="")
  accuracy_wide[i, 5] = paste(accuracy_sum[(i-1)*2+2, 5], "(", accuracy_sum[(i-1)*2+2, 6], ")", sep="")
  accuracy_wide[i, 6] = paste(accuracy_sum[(i-1)*2+1, 7], "(", accuracy_sum[(i-1)*2+1, 8], ")", sep="")
}

accuracy_wide$Trait <- c("Orl_Height", "Orl_Height", 
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
                         "GrainYield", "GrainYield", "GrainYield", "GrainYield", "GrainYield")

saveRDS(accuracy_wide, "plot/accuracy_wide.rds")
accuracy_wide <- readRDS("plot/accuracy_wide.rds")

accuracy_wide
       Case_Study             Trait     Prediction            r          r_y             h2
1  poplar_Rincent        Orl_Height        genomic  0.565(0.03) 0.684(0.032)  0.687(0.043)
2  poplar_Rincent        Orl_Height       phenomic  0.656(0.01) 0.714(0.077)  0.687(0.043)
3  poplar_Rincent      Orl_BudFlush        genomic 0.722(0.019) 0.801(0.026)  0.813(0.012)
4  poplar_Rincent      Orl_BudFlush       phenomic 0.095(0.055) 0.017(0.021)  0.813(0.012)
5  poplar_Rincent      Sav_BudFlush        genomic  0.729(0.02) 0.848(0.041)  0.748(0.043)
6  poplar_Rincent      Sav_BudFlush       phenomic 0.333(0.039) 0.201(0.028)  0.748(0.043)
7  poplar_Rincent        Orl_BudSet        genomic 0.751(0.031)   0.86(0.05)  0.745(0.049)
8  poplar_Rincent        Orl_BudSet       phenomic  0.545(0.02)  0.55(0.028)  0.745(0.049)
9  poplar_Rincent        Sav_BudSet        genomic  0.608(0.04) 0.669(0.043)   0.828(0.03)
10 poplar_Rincent        Sav_BudSet       phenomic 0.337(0.017) 0.352(0.016)   0.828(0.03)
11 poplar_Rincent Orl_Circumference        genomic  0.61(0.018) 0.722(0.039) 0.726(0.048)
12 poplar_Rincent Orl_Circumference       phenomic 0.718(0.022) 0.754(0.021) 0.726(0.048)
13 poplar_Rincent Sav_Circumference        genomic 0.777(0.017) 0.987(0.013) 0.542(0.024)
14 poplar_Rincent Sav_Circumference       phenomic 0.809(0.014) 0.736(0.024) 0.542(0.024)
15 poplar_Rincent          Orl_Rust        genomic 0.622(0.015) 0.681(0.024) 0.838(0.035)
16 poplar_Rincent          Orl_Rust       phenomic 0.443(0.046) 0.452(0.044) 0.838(0.035)
17  wheat_Rincent    Dry_GrainYield        genomic 0.238(0.067) 0.407(0.121) 0.616(0.128)
18  wheat_Rincent    Dry_GrainYield phenomic-grain 0.543(0.047) 0.641(0.084) 0.616(0.128)
19  wheat_Rincent    Dry_GrainYield  phenomic-leaf 0.314(0.066) 0.334(0.077) 0.616(0.128)
20  wheat_Rincent   Dry_HeadingDate        genomic 0.487(0.041) 0.657(0.066) 0.616(0.135)
21  wheat_Rincent   Dry_HeadingDate phenomic-grain 0.806(0.033) 0.855(0.013) 0.616(0.135)
22  wheat_Rincent   Dry_HeadingDate  phenomic-leaf  0.848(0.02) 0.893(0.017) 0.616(0.135)
23  wheat_Rincent    Irr_GrainYield        genomic   0.4(0.059) 0.668(0.088)  0.58(0.159)
24  wheat_Rincent    Irr_GrainYield phenomic-grain 0.527(0.053)  0.65(0.071)  0.58(0.159)
25  wheat_Rincent    Irr_GrainYield  phenomic-leaf 0.276(0.083) 0.229(0.089)  0.58(0.159)
26  wheat_Rincent   Irr_HeadingDate        genomic 0.511(0.033) 0.662(0.056) 0.682(0.096)
27  wheat_Rincent   Irr_HeadingDate phenomic-grain 0.618(0.027)  0.74(0.035) 0.682(0.096)
28  wheat_Rincent   Irr_HeadingDate  phenomic-leaf   0.8(0.023) 0.862(0.027) 0.682(0.096)
29           Winn        GrainYield        genomic  0.307(0.01) 0.711(0.032) 0.189(0.014)
30           Winn        GrainYield phenomic-fixed 0.472(0.014) 0.757(0.192) 0.189(0.014)
31           Winn        GrainYield phenomic-mixed  0.47(0.015)  0.545(0.21) 0.189(0.014)
32         Krause        GrainYield        genomic 0.283(0.028) 0.619(0.052)  0.201(0.02)
33         Krause        GrainYield       phenomic  0.15(0.051) 0.166(0.056)  0.201(0.02)



























 

 
















