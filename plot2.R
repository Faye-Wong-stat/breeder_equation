setwd("~/breeder_equation_project/")
library(ggplot2)
library(cowplot)
library(foreach)



r_y <- seq(0, 1, 0.01)
h2 <- c(0.1, 0.4, 0.8)
cor_e <- c(0, 0.2, 0.4, 0.6, 0.8)

bias <- data.frame(r_y = rep(r_y, each=15), 
                   h2 = rep(rep(h2, each=5), 101), 
                   cor_e = rep(cor_e, 101*3))

bias$r_g <- bias$r_y / sqrt(bias$h2)
bias$r_p <- (bias$r_y - bias$cor_e*sqrt(1-bias$h2)) / sqrt(bias$h2)

bias$r_g[bias$r_g > 1] <- 1
bias$r_p[bias$r_p > 1] <- 1
bias$r_p[bias$r_p < -1] <- -1

bias$h2 <- paste("h^2 == ", bias$h2, sep="")
bias$cor_e <- factor(bias$cor_e)

p1 <- ggplot(bias, aes(r_y, r_g)) + 
  geom_line(color="blue") + 
  geom_abline(intercept=0, slope=1) + 
  facet_grid(~h2, labeller = label_parsed) + 
  coord_cartesian(ylim = c(0, 1)) + 
  xlab(expression(rho[y])) + 
  ylab(expression(rho)) + 
  # coord_fixed(ratio=1) + 
  theme_minimal_grid(font_size=7) + 
  theme(axis.title.x = element_blank()) 

p2 <- ggplot(bias, aes(r_y, r_p, color=cor_e)) + 
  geom_line() + 
  geom_abline(intercept=0, slope=1) + 
  facet_grid(~h2, labeller = label_parsed) + 
  # ylim(0, 1) + 
  coord_cartesian(ylim = c(0, 1)) +
  # coord_fixed(ratio=1) + 
  xlab(expression(rho[y])) + 
  ylab(expression(rho)) + 
  guides(color=guide_legend(title=expression(rho[e]))) + 
  theme_minimal_grid(font_size=7) + 
  theme(legend.position = "bottom", 
        legend.justification = "center")

legend2 <- get_legend(p2)

save_plot(paste("plot2/","bias.pdf", sep=""), 
          plot_grid(p1, p2 + theme(legend.position = "none"), 
                    legend2,
                    nrow=3, align = "h", axis="tb", labels=c("A", "B"), 
                    rel_heights = c(1, 1, 0.1)), 
          base_width = 5.6)



k = 20

n_gs = 10^seq(log10(100),log10(10000),length=10)
rel_costs = c(1,1/2,1/10,1/100)
res_i = foreach(n_g = n_gs,.combine = rbind) %do% {
  foreach(rel_cost = rel_costs,.combine = rbind) %do% {
    n_p = n_g/rel_cost
    i_p = dnorm(qnorm(1-k/n_p))*n_p/k
    i_g = dnorm(qnorm(1-k/n_g))*n_g/k
    data.frame(n_g=n_g,n_p=n_p,rel_cost=rel_cost,
               i_p=i_p,i_g=i_g,i_ratio = i_p/i_g)
  }
}

p = ggplot(res_i,aes(x=n_g,y=i_ratio)) + theme_minimal() +
  geom_line(aes(group = rel_cost,color=factor(n_p/n_g))) +
  xlab("n in genomic selection") + ylab(expression(i[p]/i[g])) +
  scale_x_log10() + expand_limits(y=.9) + labs(color = expression(n[p]/n[g])) +
  guides(color = guide_legend(nrow=1)) +
  theme_minimal_grid(font_size=7) + 
  theme(legend.position = 'bottom', 
        legend.justification = "center");p
save_plot(p,file = 'plot2/i_ratio.pdf',base_width = 2.8,base_height = 2.8)
