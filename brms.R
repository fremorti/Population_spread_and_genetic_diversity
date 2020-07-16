library(car)
library(tidyverse)
library(ggplot2)
library(brms)
library(gridExtra)

#reading edge.csv
edge <- read.csv('edge.csv', header = 1, sep = ';')
edge <- gather(edge, key = 'day', value = 'edge', -environment, -line, -div, convert = 1)
edge$day <- as.numeric(substring(edge$day, 2, 3))
edge$id <- paste0(edge$line, edge$div, edge$environment)
edge$environment <- factor(edge$environment, levels = c('B', 'G', 'T'), labels = c('bean', 'gradient', 'tomato'))

#calculating coefficient of variance (CV)
agg <- aggregate(edge[5], list(edge$day, edge$environment, edge$div), sd)
agg1 <- aggregate(edge[5], list(edge$day, edge$environment, edge$div), mean)
colnames(agg) <- c('day', 'environment', 'div', 'sd')
agg$mean <- agg1[,'edge']
agg$CV <- agg$sd/agg$mean

#reading dens.csv
dens <- read.csv('dens.csv', header = 1, sep = ';')
d <- dens[7:20]
dens$edge <- apply(t(apply(d, 1, cumsum))>=dens$total, 1, function(x) {min(which(x))})
dens$id <- paste0(dens$line, dens$div)
dens$sep <- paste0(dens$line, dens$div, dens$day, dens$environment)
dens$environment <- factor(dens$environment, levels = c('B', 'G', 'T'), labels = c('bean', 'gradient', 'tomato'))



##########################
# mean population spread #
##########################

#model
eform1 <- brmsformula(edge  ~ 1 + environment*day*div + (1 + day|line))
get_prior(eform1, data = edge)
fite1 <- brm(data = edge, family = gaussian,
             formula = eform1,
             prior = c(prior(normal(0,4), class = Intercept),
                       prior(normal(0,4), class = b),
                       prior(cauchy(0,2), class = sd), 
                       prior(cauchy(0,2), class = sigma),
                       prior(lkj(2), class = cor)),
             iter = 5000, warmup = 2000, chains = 2, cores = 2)
fite1$fit



#extract posteriors
parame1 <- fite1 %>% posterior_samples() %>% mutate('slope Gradient' = .$'b_day:divsfl'+.$'b_environmentgradient:day:divsfl', 
                                                    'slope Tomato' = .$'b_day:divsfl'+.$'b_environmenttomato:day:divsfl',
                                                    'slope Bean' = .$'b_day:divsfl',
                                                    dG = .$b_divsfl + .$'b_environmentgradient:divsfl',
                                                    dT = .$b_divsfl + .$'b_environmenttomato:divsfl',
                                                    dB = .$b_divsfl)

#collect differences in slopes and intercepts from posterior
slope_differences <- parame1 %>% dplyr::select('slope Bean', 'slope Gradient','slope Tomato') %>% gather(key = environment, value = d)
int_differences <- parame1 %>% dplyr::select(dG, dT, dB) %>% gather(key = environment, value = d)

#plot differences in intercepts
(e1 <- ggplot(int_differences, aes(x = environment, y = d, color = environment))+
    geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, aes(fill = environment))+
    geom_hline(yintercept = 0, linetype = 'dashed', size = 1)+
    ggtitle('differences in estimated intercepts')+
    scale_x_discrete(labels = c('bean', 'gradient', 'tomato'))+
    scale_fill_brewer( palette = "Dark2")+
    scale_color_brewer( palette = "Dark2")+
    ylab('difference in intercept')+
    theme(legend.position = "none"))
#plot differences in slopes
(e2<- ggplot(slope_differences, aes(x = environment, y = d, color = environment))+
    geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, aes(fill = environment))+
    geom_hline(yintercept = 0, linetype = 'dashed', size = 1)+
    ggtitle('differences in estimated spread rate')+
    scale_fill_brewer( palette = "Dark2")+
    scale_x_discrete(labels = c('bean', 'gradient', 'tomato'))+
    scale_color_brewer( palette = "Dark2")+
    ylab('difference in slope')+
    theme(legend.position = "none"))

#calculate and plot model fit
nde <- tibble(day = rep(1:35, 6), environment = rep(c('bean', 'tomato', 'gradient'), 2, each = 35), div = rep(c('mix', 'sfl'), each = 105))
poste1 <- fitted(fite1, newdata = nde, re_formula = NA, probs = c(0.09, 0.91))%>%
  as_tibble()%>% bind_cols(nde)

(e3 <- ggplot(poste1, aes(x = day, y = Estimate, group = div))+
    geom_point(data = edge, aes(x = day, y = edge, color = div))+
    geom_line(aes(color = div), size = 1.25)+
    geom_ribbon(aes(ymin = Q9, ymax = Q91), fill = "black", alpha = 0.1) +
    ylab('mean spread in each environment')+
    labs(color = 'genetic diversity:')+
    ggtitle('population edge')+
    scale_color_discrete(labels = c("mix", "single female line"))+
    facet_grid(.~environment)+
    theme(legend.position=c(0.65,1.25), legend.direction = 'horizontal', plot.margin=unit(c(8.5,5.5, 5.5, 5.5),'pt')))

#assemble fig 2
grid.arrange(grobs = list(e3, e2), widths = c(1), layout_matrix = rbind(c(1), c(2)))
#assemble fig S1
grid.arrange(grobs = list(e3, e1, e2), widths = c(1), layout_matrix = rbind(c(1), c(2), c(3)))


################################################
# Coefficient of variance of population spread #
################################################

# this models both mean and std dev in population spread
eform2 <- brmsformula(edge  ~ 1 + environment*day*div, sigma ~ 1 + day*environment*div)
get_prior(eform2, data = edge)
fite2 <- brm(data = edge, family = gaussian,
             formula = eform2,
             prior = c(prior(normal(0,4), class = Intercept),
                       prior(normal(0,4), class = b)),
             iter = 5000, warmup = 2000, chains = 2, cores = 2)
fite2$fit

#sample posterior and calculated the estimated CV
poste2 <- fitted(fite2, newdata = nde, re_formula = NA, probs = c(0.09, 0.91))%>%
  as_tibble()%>% bind_cols(nde)
poste1_mean <- fitted(fite2, newdata = nde, re_formula = NA, summary = 0)%>%as_tibble()
poste1_sigma <- fitted(fite2, newdata = nde, re_formula = NA, dpar = 'sigma', summary = 0)%>% as_tibble()
poste1_CV <- (poste1_sigma/poste1_mean)
pCV_mean <- poste1_CV%>%as_tibble()%>%summarise_all(mean)%>%gather(key = 'row', value = 'mean')
pCV_Q9 <- poste1_CV%>%as_tibble()%>%summarise_all(quantile, probs = 0.09)%>%gather(key = 'row', value = 'Q9')
pCV_Q91 <- poste1_CV%>%as_tibble()%>%summarise_all(quantile, probs = 0.91)%>%gather(key = 'row', value = 'Q91')
poste1___ <- poste2%>%mutate(CV_mean = pCV_mean$mean, CV_Q9 = pCV_Q9$Q9, CV_Q91 = pCV_Q91$Q91)


#fig 3: plot CV
(splot <- ggplot(poste1___, aes(x = day, y = CV_mean, color = div))+
    geom_point(data = agg, aes(x = day, y = CV, color = div))+
    geom_line(size = 1.3)+
    geom_ribbon(aes(ymin = CV_Q9, ymax = CV_Q91), fill = "black", alpha = 0.15, linetype=0) +
    facet_grid(.~environment)+
    ggtitle('Variability in population spread')+
    labs(color = 'genetic diversity:')+
    scale_color_discrete(labels = c('mix', 'single female line'))+
    ylab('CV')+
    theme(legend.position = 'bottom'))


####################
# population sizes #
####################

#model population size regressed over time
fits1 <- brm(data = dens, family = 'negbinomial',
             total  ~ 1 + environment*week*div + (1+week|line),
             prior = c(prior(normal(0,4), class = Intercept),
                       prior(normal(0,4), class = b),
                       prior(cauchy(0,2), class = sd),
                       prior(lkj(2), class = cor),
                       prior(cauchy(0,2), class = shape)),
             iter = 5000, warmup = 2000, chains = 2, cores = 2)
fits1$fit


#sample posterior
params1 <- fits1 %>% posterior_samples() %>% mutate(slopeG = .$'b_week:divsfl'+.$'b_environmentgradient:week:divsfl', 
                                                    slopeT = .$'b_week:divsfl'+.$'b_environmenttomato:week:divsfl',
                                                    slopeB = .$'b_week:divsfl',
                                                    dG = .$b_divsfl + .$'b_environmentgradient:divsfl',
                                                    dT = .$b_divsfl + .$'b_environmenttomato:divsfl',
                                                    dB = .$b_divsfl,
                                                    sGm = .$b_week + .$'b_environmentgradient:week',
                                                    sTm = .$b_week + .$'b_environmenttomato:week',
                                                    sBm = .$b_week,
                                                    sGs = .$b_week + .$'b_environmentgradient:week'+ .$'b_week:divsfl'+.$'b_environmentgradient:week:divsfl', 
                                                    sTs = .$b_week + .$'b_environmenttomato:week' + .$'b_week:divsfl'+.$'b_environmenttomato:week:divsfl',
                                                    sBs = .$b_week + .$'b_week:divsfl',
                                                    iGm = .$b_Intercept + .$b_environmentgradient,
                                                    iTm = .$b_Intercept + .$b_environmenttomato,
                                                    iBm = .$b_Intercept,
                                                    iGs = .$b_Intercept + .$b_environmentgradient + .$b_divsfl + .$'b_environmentgradient:divsfl', 
                                                    iTs = .$b_Intercept + .$b_environmenttomato + .$b_divsfl + .$'b_environmenttomato:divsfl',
                                                    iBs = .$b_Intercept + .$b_divsfl,
                                                    ID = rownames(.))
#collect differences in intercepts and slopes
slope_differences <- params1 %>% select(slopeG, slopeT, slopeB) %>% gather(key = environment, value = d)
int_differences <- params1 %>% select(dG, dT, dB) %>% gather(key = environment, value = d)

#plot differences in intercepts
(s1 <- ggplot(int_differences, aes(x = environment, y = d, color = environment))+
    geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, aes(fill = environment))+
    geom_hline(yintercept = 0, linetype = 'dashed', size = 1)+
    scale_fill_brewer( palette = "Dark2")+
    scale_color_brewer( palette = "Dark2")+
    scale_x_discrete(labels = c('bean', 'gradient', 'tomato'))+
    ggtitle('estimated differences in initial population sizes')+
    ylab('difference in log intercept (sfl-mix)')+
    theme(legend.position = "none"))
#plot differences in slopes
(s2<- ggplot(slope_differences, aes(x = environment, y = d, color = environment))+
    geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, aes(fill = environment))+
    geom_hline(yintercept = 0, linetype = 'dashed', size = 1)+
    scale_fill_brewer( palette = "Dark2")+
    scale_color_brewer( palette = "Dark2")+
    scale_x_discrete(labels = c('bean', 'gradient', 'tomato'))+
    ggtitle('estimated differences in population increase')+
    ylab('difference in slope (sfl-mix)')+
    theme(legend.position = "none"))

#fig 4: plot model fits
nds <- tibble(week = rep(1:5, 6), environment = rep(c('bean', 'tomato', 'gradient'), 2, each = 5), div = rep(c('mix', 'sfl'), each = 15))
posts1 <- fitted(fits1, newdata = nds, re_formula = NA, probs = c(0.09, 0.91))%>%
  as_tibble()%>% bind_cols(nds)

(s3 <- ggplot(posts1, aes(x = week, y = Estimate, group = div))+
    geom_line(aes(color = div))+
    geom_ribbon(aes(ymin = Q9, ymax = Q91), fill = "black", alpha = 0.1) +
    geom_point(data = dens, aes(x = week, y = total, color = div))+
    ylab('total population size')+
    ggtitle('population size')+
    labs(color = 'genetic diversity:')+
    scale_color_discrete(labels = c('mix', 'single female line'))+
    facet_grid(.~environment)+
    theme(legend.position=c(0.7,1.25), legend.direction = 'horizontal', plot.margin=unit(c(8.5,5.5, 5.5, 5.5),'pt')))

#assemble fig S3.1
grid.arrange(grobs = list(s3, s1, s2), widths = c(1), layout_matrix = rbind(c(1),c(2), c(3)))

#collect slopes and intercepts per posterior sample
is_corr1 <- params1 %>% select(sGm, sTm, sBm, sGs, sTs, sBs, iGm, iTm, iBm, iGs, iTs, iBs, ID)%>%pivot_longer(cols = c(sGm, sTm, sBm, sGs, sTs, sBs, iGm, iTm, iBm, iGs, iTs, iBs), names_to = 'c', values_to = 'value')%>%
  mutate(si = str_sub(.$c, 1, 1), environment = str_sub(.$c, 2, 2), div = str_sub(.$c, 3, 3))%>%pivot_wider(id_cols = c(environment, div, ID), names_from = si, values_from = value)
is_dcorr1 <- params1 %>% select(dG, dT, dB, sG = slopeG, sT = slopeT, sB = slopeB, ID)%>%pivot_longer(cols = c(dG, dT, dB, sG, sT, sB), names_to = 'c', values_to = 'value')%>%
  mutate(si = str_sub(.$c, 1, 1), environment = str_sub(.$c, 2, 2))%>%pivot_wider(id_cols = c(environment, ID), names_from = si, values_from = value)

#plot correlations in slopes and intercepts
(s1_corr<- ggplot(is_corr1, aes(x = i, y = s, color = div))+
    geom_point(alpha = 0.05)+
    facet_grid(~environment)+
    ggtitle('estimated slope intercept correlation')+
    ylab('estimated slope')+
    xlab('estimated intercept'))

#fig S3.2: plot correlations in differences in slopes and differences in intercepts
(s1_corr<- ggplot(is_dcorr1, aes(x = d, y = s))+
    geom_point(alpha = 0.05, color = 'blue4')+
    facet_grid(~environment)+
    geom_hline(yintercept = 0)+
    geom_vline(xintercept = 0)+
    ggtitle('estimated slope difference-intercept difference correlation')+
    ylab('estimated slope difference')+
    xlab('estimated intercept difference'))

#model population size regressed over furthest occupied patch
fits3 <- brm(data = dens, family = 'Negbinomial',
             total  ~ 1 + environment*edge*div + (1 + edge|line),
             prior = c(prior(normal(0,4), class = Intercept),
                       prior(normal(0,4), class = b),
                       prior(cauchy(0,2), class = sd),
                       prior(lkj(2), class = cor),
                       prior(cauchy(0,2), class = shape)),
             iter = 5000, warmup = 2000, chains = 2, cores = 2)
fits3$fit

#sample posterior
params3 <- fits3 %>% posterior_samples() %>% mutate(slopeG = .$'b_edge:divsfl'+.$'b_environmentgradient:edge:divsfl', 
                                                    slopeT = .$'b_edge:divsfl'+.$'b_environmenttomato:edge:divsfl',
                                                    slopeB = .$'b_edge:divsfl',
                                                    dG = .$b_divsfl + .$'b_environmentgradient:divsfl',
                                                    dT = .$b_divsfl + .$'b_environmenttomato:divsfl',
                                                    dB = .$b_divsfl,
                                                    sGm = .$b_edge + .$'b_environmentgradient:edge',
                                                    sTm = .$b_edge + .$'b_environmenttomato:edge',
                                                    sBm = .$b_edge,
                                                    sGs = .$b_edge + .$'b_environmentgradient:edge'+ .$'b_edge:divsfl'+.$'b_environmentgradient:edge:divsfl', 
                                                    sTs = .$b_edge + .$'b_environmenttomato:edge' + .$'b_edge:divsfl'+.$'b_environmenttomato:edge:divsfl',
                                                    sBs = .$b_edge + .$'b_edge:divsfl',
                                                    iGm = .$b_Intercept + .$b_environmentgradient,
                                                    iTm = .$b_Intercept + .$b_environmenttomato,
                                                    iBm = .$b_Intercept,
                                                    iGs = .$b_Intercept + .$b_environmentgradient + .$b_divsfl + .$'b_environmentgradient:divsfl', 
                                                    iTs = .$b_Intercept + .$b_environmenttomato + .$b_divsfl + .$'b_environmenttomato:divsfl',
                                                    iBs = .$b_Intercept + .$b_divsfl,
                                                    ID = rownames(.))

slope_differences <- params3 %>% select(slopeG, slopeT, slopeB) %>% gather(key = environment, value = d)
int_differences <- params3 %>% select(dG, dT, dB) %>% gather(key = environment, value = d)

#collect slopes and intercepts per posterior sample
is_corr3 <- params3 %>% select(sGm, sTm, sBm, sGs, sTs, sBs, iGm, iTm, iBm, iGs, iTs, iBs, ID)%>%pivot_longer(cols = c(sGm, sTm, sBm, sGs, sTs, sBs, iGm, iTm, iBm, iGs, iTs, iBs), names_to = 'c', values_to = 'value')%>%
  mutate(si = str_sub(.$c, 1, 1), environment = str_sub(.$c, 2, 2), div = str_sub(.$c, 3, 3))%>%pivot_wider(id_cols = c(environment, div, ID), names_from = si, values_from = value)
is_dcorr3 <- params3 %>% select(dG, dT, dB, sG = slopeG, sT = slopeT, sB = slopeB, ID)%>%pivot_longer(cols = c(dG, dT, dB, sG, sT, sB), names_to = 'c', values_to = 'value')%>%
  mutate(si = str_sub(.$c, 1, 1), environment = str_sub(.$c, 2, 2))%>%pivot_wider(id_cols = c(environment, ID), names_from = si, values_from = value)

#plot differences in log intercepts
(s3_1 <- ggplot(int_differences, aes(x = environment, y = d, color = environment))+
    geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, aes(fill = environment))+
    geom_hline(yintercept = 0, linetype = 'dashed', size = 1)+
    scale_fill_brewer( palette = "Dark2")+
    scale_color_brewer( palette = "Dark2")+
    scale_x_discrete(labels = c('bean', 'gradient', 'tomato'))+
    ggtitle('differences in estimated log intercepts')+
    ylab('difference in intercept (sfl-mix)')+
    theme(legend.position = "none"))
#plot differences in log slopes
(s3_2<- ggplot(slope_differences, aes(x = environment, y = d, color = environment))+
    geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), alpha = 0.3, aes(fill = environment))+
    geom_hline(yintercept = 0, linetype = 'dashed', size = 1)+
    scale_fill_brewer( palette = "Dark2")+
    scale_color_brewer( palette = "Dark2")+
    scale_x_discrete(labels = c('bean', 'gradient', 'tomato'))+
    ggtitle('differences in estimated population density')+
    ylab('difference in slope (sfl-mix)')+
    theme(legend.position = "none"))

#plot correlations in slopes and intercepts
(s3_corr<- ggplot(is_corr3, aes(x = i, y = s, color = div))+
    geom_point(alpha = 0.05)+
    facet_grid(~environment)+
    ggtitle('estimated slope intercept correlation')+
    ylab('estimated slope')+
    xlab('estimated intercept'))
#fig S4.2: plot correlations in differences in slopes and differences in intercepts
(s3_dcorr<- ggplot(is_dcorr3, aes(x = d, y = s, color = 1))+
    geom_point(alpha = 0.05, color = 'blue4')+
    facet_grid(~environment)+
    geom_hline(yintercept = 0)+
    geom_vline(xintercept = 0)+
    ggtitle('estimated slope difference-intercept difference correlation')+
    ylab('estimated slope difference')+
    xlab('estimated intercept difference'))
    

#maximum edge for each combination of diversity and environment to construct new data to fit the model to, only in the range of edges per combination
maxes <- c(max(dens[dens$div == 'mix' & dens$environment == 'bean',]$edge),
           max(dens[dens$div == 'mix' & dens$environment == 'gradient',]$edge),
           max(dens[dens$div == 'mix' & dens$environment == 'tomato',]$edge),
           max(dens[dens$div == 'sfl' & dens$environment == 'bean',]$edge),
           max(dens[dens$div == 'sfl' & dens$environment == 'gradient',]$edge),
           max(dens[dens$div == 'sfl' & dens$environment == 'tomato',]$edge))
nds <- tibble(edge = c(1:maxes[1], 1:maxes[2], 1:maxes[3], 1:maxes[4], 1:maxes[5], 1:maxes[6]), 
              environment = c(rep('bean', maxes[1]), rep('gradient', maxes[2]), rep('tomato', maxes[3]),rep('bean', maxes[4]), rep('gradient', maxes[5]), rep('tomato', maxes[6])), 
              div = c(rep('mix', sum(maxes[1:3])), rep('sfl', sum(maxes[4:6]))))
#calculate the fit
posts3 <- fitted(fits3, newdata = nds, re_formula = NA, probs = c(0.09, 0.91))%>%
  as_tibble()%>% bind_cols(nds)

#fig 5: plot model fits
(s3_3 <- ggplot(posts3, aes(x = edge, y = Estimate, group = div))+
    geom_line(aes(color = div))+
    geom_ribbon(aes(ymin = Q9, ymax = Q91), fill = "black", alpha = 0.1) +
    geom_point(data = dens, aes(x = edge, y = total, color = div))+
    scale_color_discrete(labels = c('mix', 'single female line'))+
    xlab('occupied patches')+
    ylab('total population size')+
    ggtitle('effect range size on population size')+
    facet_grid(.~environment)+
    labs(color = 'genetic diversity:')+
    theme(legend.position='bottom'))
#assemble fig S4.1
grid.arrange(grobs = list(s3_3,s3_1, s3_2), widths = c(1), layout_matrix = rbind(c(1),c(2), c(3)))
