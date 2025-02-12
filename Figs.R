# library(ggplot2)
# library(cowplot)
# library(foreach)
# 
# 
# 
# ability = seq(0,1,length=100)#c(0,.4,.6,.8,.99)
# h2s = c(.1,.4,.8)#seq(0,1,length=100)
# cor_es = seq(0,.8,length=5)#c(-.1,0,.3,.5,.9)
# accuracy = foreach(h2 = h2s,.combine = rbind) %do% {
#   foreach(cor_e = cor_es,.combine = rbind) %do% {
#     va = h2
#     ve = 1-va
#     vy = 1
#     vhat = 1
#     cov_e = cor_e*sqrt(vhat*ve)
#     # ability = (cov_a+cov_e)/sqrt((va+ve)*vhat)
#     cov_a = ability*sqrt((va+ve)*vhat) - cov_e
#     # cov_a = accuracy*sqrt(vhat*va)
#     # cor(hat,y) = cov(hat,y)/sqrt(vy*vhat) = cov_v+cov_e
#     data.frame(h2=h2,cor_e=cor_e,
#                ability = ability,
#                accuracy_g = pmin(1,ability/sqrt(h2)),
#                accuracy_p = pmax(-1,pmin(1,cov_a/sqrt(va*vhat))),
#                accuracy_pb = pmax(-1,pmin(1,(ability-cor_e*sqrt(1-h2))/sqrt(h2))),
#                accuracy_p2 = pmax(-1,pmin(1,ability/sqrt(h2)-cor_e*sqrt(1-h2)/sqrt(h2)))
# 
#     )
#   }
# }
# 
# # ggplot(accuracy) +
# #   geom_hline(yintercept = ability) +
# #   geom_line(aes(x=cor_e,y=accuracy_p,group = ability,color=factor(ability))) +
# #   facet_wrap(~h2)
# #
# # ggplot(subset(accuracy,cor_e==0.5)) +
# #   geom_line(aes(x=h2,y=accuracy_g,group = ability,color=factor(ability)))
# #
# 
# (p1=ggplot(subset(accuracy,cor_e==0)) + xlim(c(0,1))+ylim(c(-1,1))+
#     geom_abline(slope=1,intercept=0)+
#     geom_line(aes(x=ability,y=accuracy_g),color='blue',linewidth=1.5) +
#     facet_wrap(~h2,labeller = label_both) + theme_minimal()+
#     xlab('Predictive ability') + ylab('Accuracy'))
# 
# accuracy$cor_e = factor(accuracy$cor_e)
# (p2=ggplot(accuracy) +
#     # geom_hline(yintercept = ability) +
#     geom_abline(slope=1,intercept=0)+
#     geom_line(aes(x=ability,y=accuracy_p,group = cor_e,color=cor_e),linewidth=1) +
#     facet_wrap(~h2,labeller = label_both) + theme_minimal() +
#     xlab('Predictive ability') + ylab('Accuracy'))
# 
# p3 = plot_grid(plot_grid(p1,p2+theme(legend.position = 'none'),ncol=1),
#           get_legend(p2+theme(legend.position = 'bottom')),ncol=1,rel_heights = c(.9,.1))
# 
# save_plot('bias_both.pdf',p3,base_width = 6.5,base_height = 4)
# 
# save_plot('bias_genomic.pdf',p1,base_height = 4)
# save_plot('bias_phenomic.pdf',p2+theme(legend.position = 'none'),base_height = 4)
# save_plot('bias_phenomic_legend.pdf',get_legend(p2),base_height = 4)
# 
# 
# 
# 
# # say we select 20 lines as parents regardless of population size
# k = 20
# 
# n_gs = 10^seq(log10(100),log10(10000),length=10)
# rel_costs = c(1,1/2,1/10,1/100)
# res_i = foreach(n_g = n_gs,.combine = rbind) %do% {
#   foreach(rel_cost = rel_costs,.combine = rbind) %do% {
#     n_p = n_g/rel_cost
#     i_p = dnorm(qnorm(1-k/n_p))*n_p/k
#     i_g = dnorm(qnorm(1-k/n_g))*n_g/k
#     data.frame(n_g=n_g,n_p=n_p,rel_cost=rel_cost,
#                i_p=i_p,i_g=i_g,i_ratio = i_p/i_g)
#   }
# }
# 
# p = ggplot(res_i,aes(x=n_g,y=i_ratio)) + theme_minimal() +
#   geom_line(aes(group = rel_cost,color=factor(n_p/n_g))) +
#   xlab('# lines Genomic') + ylab('Relative intensity') +
#   scale_x_log10() + expand_limits(y=.9) + labs(color = 'relative # lines') +
#   guides(color = guide_legend(nrow=2)) +
#   theme(legend.position = 'bottom');p
# save_plot(p,file = 'i_ratio.pdf',base_height = 3,base_asp=1)
# 
# 
# ### old code below
# 
# x = seq(-2,-2)
# 
# 
# y1 = rnorm(1e3,10,2)
# y2 = rnorm(3e3,10,2)
# y1s = sort(y1,decreasing = T)[1:50]
# y2s = sort(y2,decreasing = T)[1:50]
# 
# range(c(y1,y2))
# h1 = hist(y1,breaks = c(-Inf,seq(2,18,by=.5),Inf),plot=F)
# h2 = hist(y2,breaks = c(-Inf,seq(2,18,by=.5),Inf),plot=F)
# 
# plot(h1,ylim = c(0,max(h2$counts)),xlim=c(2,18),freq=T)
# hist(y1s,breaks = h1$breaks,freq=T,col=2,add=T)
# plot(h2,ylim = c(0,max(h2$counts)),xlim=c(2,18),freq=T)
# hist(y2s,breaks = h1$breaks,freq=T,col=2,add=T)
# plot(h1,breaks = h1$breaks,freq=T,col='grey70',add=T)
# hist(y1s,breaks = h1$breaks,freq=T,col=4,add=T)
# 
# x = seq(2,18,length=200)
# n1 = 200;n2 = 500
# mu=10;sd=2
# y1 = dnorm(x,mu,sd)*n1
# y2 = dnorm(x,mu,sd)*n2
# y3 = dnorm(x,mu,sd*1.5)*n2
# t1 = qnorm(1-20/200,mu,sd,lower.tail=T)
# t2 = qnorm(1-20/1000,mu,sd,lower.tail=T)
# t3 = qnorm(1-20/1000,mu,sd*1.5,lower.tail=T)
# xs1 = seq(t1,max(x),length=20)
# ys1 = dnorm(xs1,mu,sd)*n1
# # xs1 = c(xs1,xs1[length(xs1):1])
# # ys1 = c(ys1,rep(0,length(xs1)))
# xs2 = seq(t2,max(x),length=20)
# ys2 = dnorm(xs2,mu,sd)*n2
# xs3 = seq(t3,max(x),length=20)
# ys3 = dnorm(xs3,mu,sd*1.5)*n2
# # xs2 = c(xs2,xs2[length(xs2):1])
# 
# (p1 = ggplot() + theme_nothing() + xlab('Trait') + ylab('Frequency') +
#   geom_ribbon(data=data.frame(x = x,y=y2),aes(x=x,ymax = y,ymin=0),fill='grey90') +
#   geom_line(data=data.frame(x = x,y=y2),aes(x=x,y=y)) +
#   geom_ribbon(data=data.frame(x = xs2,y=ys2),aes(x=x,ymax = y,ymin=0),fill='salmon') +
#   geom_ribbon(data=data.frame(x = x,y=y1),aes(x=x,ymax = y,ymin=0),fill='grey70',alpha=.3) +
#   geom_line(data = data.frame(y=y1),aes(x=x,y=y)) +
#   geom_ribbon(data=data.frame(x = xs1,y=ys1),aes(x=x,ymax = y,ymin=0),fill='blue',alpha=.3))
# 
# (p3 = ggplot() + theme_nothing() + xlab('Trait') + ylab('Frequency') +
#     geom_ribbon(data=data.frame(x = x,y=y3),aes(x=x,ymax = y,ymin=0),fill='grey90') +
#     geom_line(data=data.frame(x = x,y=y3),aes(x=x,y=y)) +
#     geom_ribbon(data=data.frame(x = xs3,y=ys3),aes(x=x,ymax = y,ymin=0),fill='magenta') +
#   geom_ribbon(data=data.frame(x = x,y=y2),aes(x=x,ymax = y,ymin=0),fill='grey90',alpha=.3) +
#     geom_line(data=data.frame(x = x,y=y2),aes(x=x,y=y)) +
#     geom_ribbon(data=data.frame(x = xs2,y=ys2),aes(x=x,ymax = y,ymin=0),fill='salmon') +
#     geom_ribbon(data=data.frame(x = x,y=y1),aes(x=x,ymax = y,ymin=0),fill='grey70',alpha=.3) +
#     geom_line(data = data.frame(y=y1),aes(x=x,y=y)) +
#     geom_ribbon(data=data.frame(x = xs1,y=ys1),aes(x=x,ymax = y,ymin=0),fill='blue',alpha=.3))
# 
# 
# (p2 = ggplot() + theme_nothing() + xlab('Trait') + ylab('Frequency') + ylim(c(0,max(y2))) +
#     geom_ribbon(data=data.frame(x = x,y=y1),aes(x=x,ymax = y,ymin=0),fill='grey70',alpha=.6) +
#     geom_line(data = data.frame(y=y1),aes(x=x,y=y)) +
#     geom_ribbon(data=data.frame(x = xs1,y=ys1),aes(x=x,ymax = y,ymin=0),fill='blue',alpha=.6))
# 
# save_plot('small_normal.pdf',p2,base_height = 3)
# save_plot('large_normal.pdf',p1,base_height = 3)
# save_plot('large_normal_wide.pdf',p3,base_height = 3)
# 
# 
# n = 40
# sd=2
# BV = rnorm(n,mu,sd)
# Y_low = BV + rnorm(n,0,sd)
# Y_high = BV + rnorm(n,0,sd/3)
# (p1=ggplot(data.frame(x=Y_low,y=BV)) + ylim(range(c(Y_low,Y_high))) +
#     geom_point(aes(x=x,y=y)) + theme_nothing() + background_grid('none'))
# (p1=ggplot(data.frame(x=Y_low,y=BV)) + xlim(range(c(Y_low,Y_high))) +
#     geom_point(aes(x=x,y=y)) + geom_abline(slope = 1,intercept = 0) +
#     theme_minimal() + theme(axis.text = element_blank()) + background_grid('none') +
#     xlab('Estimated') + ylab('Breeding Value'))
# (p2=ggplot(data.frame(x=Y_low,y=BV)) + xlim(range(c(Y_low,Y_high))) +
#     geom_point(aes(x=x,y=y)) + geom_abline(slope = 1,intercept = 0) +
#     geom_point(data = data.frame(x = Y_high,y=BV),aes(x=x,y=y),color='salmon',size=1.5)+
#     theme_minimal() + theme(axis.text = element_blank()) + background_grid('none') +
#     xlab('Estimated') + ylab('Breeding Value'))
# cor(BV,Y_low)
# cor(BV,Y_high)
# save_plot('cor_low.pdf',p1,base_height = 2,base_asp = 1)
# save_plot('cor_high.pdf',p2,base_height = 2,base_asp = 1)
# ggplot(data.frame(x=Y_high,y=BV)) + geom_point(aes(x=x,y=y)) + theme_minimal() + background_grid('none')
# plot(Y_low,BV)
# 
# 
# 
# n = 80
# sd = 2
# se = sd
# GEBV = rnorm(n,mu,sd)
# BV = GEBV + rnorm(n,0,se)
# GV = BV + rnorm(n,0,sqrt(se^2/2))
# Y = GV + rnorm(n,0,sqrt(se^2/2))
# Y_hat = GV + rnorm(n,0,sqrt(se^2/2))
# d = data.frame(GEBV,BV,Y,Y_hat,GV,ID = 1:n)
# d_tall = tidyr::pivot_longer(d,cols=-c('ID','GEBV'))
# (p1=ggplot(d) + coord_equal(xlim=range(unlist(d[,c(1,4)])),ylim=range(unlist(d[,c(2,3)]))) + xlim(range(unlist(d))) +
#     # geom_smooth(aes(x=GEBV,y=BV),method='lm',color='grey30',fullrange=T) +
#     geom_abline(slope=1,intercept=0) +
#     # geom_segment(aes(x=BV,y=BV,xend=GEBV,yend=BV)) +
#     geom_point(aes(x=GEBV,y=BV)) +
#     theme_minimal() + theme(axis.text = element_blank()) + background_grid('none') +
#     xlab('Predicted BV') + ylab('Breeding Value'))
# (p2=ggplot(d) + coord_equal(xlim=range(unlist(d[,c(1,4)])),ylim=range(unlist(d[,c(2,3)]))) + xlim(range(unlist(d))) +
#     # geom_smooth(aes(x=GEBV,y=Y),method='lm',color='grey30',fullrange=T) +
#     geom_abline(slope=1,intercept=0) +
#     # geom_point(aes(x=GEBV,y=BV),alpha = 0.5) +
#     # geom_segment(aes(x=BV,y=BV,xend=GEBV,yend=BV)) +
#     # geom_segment(aes(x=GEBV,y=BV,xend=GEBV,yend=Y)) +
#     geom_point(aes(x=GEBV,y=Y),color = 'blue') +
#     theme_minimal() + theme(axis.text = element_blank()) + background_grid('none') +
#     xlab('Predicted BV') + ylab('Observed Trait (Y)'))
# cor(d$GEBV,d$BV)
# cor(d$GEBV,d$Y)
# save_plot('cor_BV.pdf',p1,base_height = 2,base_asp = 1)
# save_plot('cor_Y.pdf',p2,base_height = 2,base_asp = 1)
# 
# (p3=ggplot(d) + coord_equal(xlim=range(unlist(d[,c(1,4)])),ylim=range(unlist(d[,c(2,3)]))) + xlim(range(unlist(d))) +
#     # geom_smooth(aes(x=GEBV,y=Y),method='lm',color='grey30',fullrange=T) +
#     geom_abline(slope=1,intercept=0) +
#     geom_point(aes(x=GEBV,y=BV),alpha = 0.3) +
#     geom_segment(aes(x=BV,y=BV,xend=GEBV,yend=BV)) +
#     geom_segment(aes(x=GEBV,y=BV,xend=GEBV,yend=Y)) +
#     geom_segment(aes(x=GEBV,y=Y,xend=Y_hat,yend=Y)) +
#     geom_point(aes(x=GEBV,y=Y),color = 'salmon',alpha=0.3) +
#     geom_point(aes(x=Y_hat,y=Y),color = 'blue') +
#     theme_minimal() + theme(axis.text = element_blank()) + background_grid('none') +
#     xlab('Predicted BV') + ylab('Y'))
# 
# 
# (p1=ggplot(d) + coord_equal(xlim=range(unlist(d[,c(1,4)])),ylim=range(unlist(d[,c(2,3)]))) + xlim(range(unlist(d))) +
#     geom_smooth(aes(x=GEBV,y=BV),method='lm',color='grey30',fullrange=T) +
#     geom_point(aes(x=GEBV,y=BV)) +
#     theme_minimal() + theme(axis.text = element_blank()) + background_grid('none') +
#     xlab('Predicted BV') + ylab('Breeding Value'))
# (p2=ggplot(d) + coord_equal(xlim=range(unlist(d[,c(1,4)])),ylim=range(unlist(d[,c(2,3)]))) + xlim(range(unlist(d))) +
#     geom_smooth(aes(x=GEBV,y=Y),method='lm',color='grey30',fullrange=T) +
#     geom_point(aes(x=GEBV,y=BV),alpha = 0.1) +
#     geom_segment(aes(x=GEBV,y=BV,xend=GEBV,yend=Y)) +
#     geom_point(aes(x=GEBV,y=Y),color = 'salmon') +
#     theme_minimal() + theme(axis.text = element_blank()) + background_grid('none') +
#     xlab('Predicted BV') + ylab('Y'))
# (p3=ggplot(d) + coord_cartesian(xlim=range(unlist(d[,c(1,4)])),ylim=range(unlist(d[,c(2,3)]))) + xlim(range(unlist(d))) +
#     # geom_smooth(aes(x=Y_hat,y=Y),method='lm',color='grey30',fullrange=T) +
#     geom_point(aes(x=GEBV,y=BV),alpha = 0.1) +
#     # geom_segment(aes(x=GEBV,y=BV,xend=GEBV,yend=Y-(GV-GEBV))) +
#     # geom_segment(aes(x=GEBV,y=Y-(GV-GEBV),xend=GV,yend=Y)) +
#     geom_segment(aes(x=Y_hat,y=Y,xend=Y_hat,yend=)) +
#     geom_segment(aes(x=Y_hat,y=Y,xend=Y_hat,yend=)) +
#     geom_point(aes(x=GV,y=Y),color = 'salmon') +
#     theme_minimal() + theme(axis.text = element_blank()) + background_grid('none') +
#     xlab('Predicted Y') + ylab('Y'))
# (p2=ggplot(data.frame(x=GEBV,y=Y)) + ylim(range(c(Y_low,Y_high))) +
#     geom_point(data.frame(x=GEBV,BV),aes(x=x,y=y)) +
#     geom_point(data.frame(x=GEBV,y=BV),aes(x=x,y=y)) +
#     theme_minimal() + theme(axis.text = element_blank()) + background_grid('none') +
#     xlab('Predicted') + ylab('Y'))
# 
# 
# Y_low = BV + rnorm(n,0,sd/1.5)
# Y_high = BV + rnorm(n,0,sd/3)
# (p1=ggplot(data.frame(x=Y_low,y=BV)) + ylim(range(c(Y_low,Y_high))) +
#     geom_point(aes(x=x,y=y)) + theme_nothing() + background_grid('none'))
# (p1=ggplot(data.frame(x=Y_low,y=BV)) + xlim(range(c(Y_low,Y_high))) +
#     geom_point(aes(x=x,y=y)) + geom_abline(slope = 1,intercept = 0) +
#     theme_minimal() + theme(axis.text = element_blank()) + background_grid('none') +
#     xlab('Estimated') + ylab('Breeding Value'))
# (p2=ggplot(data.frame(x=Y_low,y=BV)) + xlim(range(c(Y_low,Y_high))) +
#     geom_point(aes(x=x,y=y)) + geom_abline(slope = 1,intercept = 0) +
#     geom_point(data = data.frame(x = Y_high,y=BV),aes(x=x,y=y),color='salmon',size=1.5)+
#     theme_minimal() + theme(axis.text = element_blank()) + background_grid('none') +
#     xlab('Estimated') + ylab('Breeding Value'))
# cor(BV,Y_low)
# cor(BV,Y_high)
# save_plot('cor_low.pdf',p1,base_height = 2,base_asp = 1)
# save_plot('cor_high.pdf',p2,base_height = 2,base_asp = 1)
# 
# 
# 
# library(foreach)
# accuracy = seq(0,1,length=100)
# h2s = c(.1,.4,.8,.99)
# cor_es = c(-.1,0,.3,.5,.9)
# ability = foreach(h2 = h2s,.combine = rbind) %do% {
#   foreach(cor_e = cor_es,.combine = rbind) %do% {
#     va = h2
#     ve = 1-va
#     vy = 1
#     vhat = 1
#     cov_a = accuracy*sqrt(vhat*va)
#     cov_e = cor_e*sqrt(vhat*ve)
#     # cor(hat,y) = cov(hat,y)/sqrt(vy*vhat) = cov_v+cov_e
#     data.frame(h2=h2,cor_e=cor_e,
#       accuracy = accuracy,
#       ability_g = pmin(accuracy*sqrt(h2),1),
#       ability_p = pmax(-1,pmin(1,(cov_a+cov_e)/sqrt((va+ve)*vhat)))
#     )
#   }
# }
# 
# ggplot(subset(ability,cor_e==0.5),aes(x=ability_g,y = accuracy)) +
#   geom_abline(slope=1,intercept=0) +
#   geom_line(aes(group = h2,color=h2))
# 
# 
# ggplot(subset(ability,ability_p>0),aes(x=ability_p,y = accuracy)) +
#   geom_abline(slope=1,intercept=0) +
#   geom_line(aes(group = cor_e,color=factor(cor_e))) + facet_wrap(~h2)
# 
# ggplot(subset(ability,ability_p>0),aes(x=ability_p,y = accuracy)) +
#   geom_abline(slope=1,intercept=0) +
#   geom_line(aes(group = h2,color=h2)) + facet_wrap(~cor_e)
# 
# P = diag(1,3)
# P[1,1] = 1
# P[1,2] = .7
# P[1,3] = .5
# P[lower.tri(P)] = t(P[upper.tri(P)])
# P = diag(sqrt(c(1,.8,.2)))%*% P %*% diag(sqrt(c(1,.8,.2)))
# chol(P)
# L = cbind(c(1,0,0),c(0,1,1))
# ability$ability_p2 = sapply(1:nrow(ability),function(i) {
#   P = diag(1,3)
#   # P[1,1] = 1
#   P[1,2] = ability$accuracy[i]
#   P[1,3] = ability$cor_e[i]
#   P[lower.tri(P)] = t(P[upper.tri(P)])
#   P = diag(sqrt(c(1,ability$h2[i],1-ability$h2[i])))%*% P %*% diag(sqrt(c(1,ability$h2[i],1-ability$h2[i])))
#   c = t(L) %*% P %*% L
#   c[1,2]
# })
# 
# 
# 
# 
# wheat = readRDS('accuracy_table_long_wheat.rds')
# wheat = subset(wheat,type %in% c('g_acc_pear','g_acc_gcor','p_acc_pear_grain','p_acc_gcor_grain') & fold == 'five folds')
# wheat$accuracy = pmin(1,wheat$accuracy)
# means = aggregate(accuracy~trait+type+adjustment+prediction,wheat,FUN=mean)
# means$SE = aggregate(accuracy~trait+type+adjustment+prediction,wheat,FUN=function(x) sd(x)/sqrt(length(x)))$accuracy
# dry_yield = subset(means,trait == 'IRR_GrainYield')
# dry_yield$Data = ifelse(dry_yield$prediction == 'genomic','Genomic','Phenomic')
# dry_yield$Measure = ifelse(grepl('pear',dry_yield$type),'Ability','Accuracy')
# position = position_dodge(width=1)
# (p1=ggplot(dry_yield) + theme_minimal()+ylim(c(0,1))+
#   geom_errorbar(aes(x = Measure,ymin = accuracy-2*SE,ymax = accuracy+2*SE,group = Data),width=.5,position = position)+
#   geom_bar(aes(x=Measure,y=accuracy,group = Data,fill=Data),stat='identity',position = position))
# save_plot('Rincent.pdf',p1,base_height = 4)
# 
# 
# 
# wheat = readRDS('accuracy_table_long.rds')
# wheat = subset(wheat,type %in% c('g_acc_pear','g_acc_gcor','p_acc_pear','p_acc_gcor') & fold == 'five folds')
# wheat$accuracy = pmin(1,wheat$accuracy)
# means = aggregate(accuracy~trait+type+adjustment+prediction,wheat,FUN=mean)
# means$SE = aggregate(accuracy~trait+type+adjustment+prediction,wheat,FUN=function(x) sd(x)/sqrt(length(x)))$accuracy
# dry_yield = subset(means,trait == 'CIRC.SAV')
# dry_yield$Data = ifelse(dry_yield$prediction == 'genomic','Genomic','Phenomic')
# dry_yield$Measure = ifelse(grepl('pear',dry_yield$type),'Ability','Accuracy')
# position = position_dodge(width=1)
# (p1=ggplot(dry_yield) + theme_minimal()+ ylim(c(0,1))+
#     geom_errorbar(aes(x = Measure,ymin = accuracy-2*SE,ymax = accuracy+2*SE,group = Data),width=.5,position = position)+
#     geom_bar(aes(x=Measure,y=accuracy,group = Data,fill=Data),stat='identity',position = position))
# save_plot('Rincent_HT.ORL.pdf',p1,base_height = 4)
# save_plot('Rincent_CIRC.SAV.pdf',p1,base_height = 4)
# 
# 
# 
