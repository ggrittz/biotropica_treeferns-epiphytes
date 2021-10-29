rm(list=ls())

library(ggplot2)
library(iNEXT)
library(RColorBrewer)
library(ggpubr)
library(tableone)
library(tidyverse)
library(data.table)
library(rstatix)

#Data necessary to estimate the empirical diversity profile 
df = read.csv('C:/Users/Master/FURB/André Luís de Gasper - Giesta/Biotropica_review 2/asymptotic_diversity_profile_estimated.csv', 
              sep = ';')
df = df %>% filter(Target == "Diversity")
df$Target <- gsub("Diversity", "Estimated diversity", df$Target)

div_profile = ggplot(df, aes(x=Order.q, y=Estimate, colour=Community)) +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2)) + 
  #scale_y_continuous(breaks = seq(0, 250, 25)) + 
  theme_bw() + 
  geom_point(aes(shape=Community), size = 1, data=df) +
  geom_line(aes(linetype=Target), size = 1.5, data = df) +
  geom_ribbon(aes(ymin=LCL, ymax=UCL,
                  fill=Community, colour=NULL), alpha=0.2) +
  labs(x="Order of q", y="Species diversity") +
  theme(legend.position = "bottom", 
        legend.title=element_blank(),
        text=element_text(size=16),
        legend.box = "vertical",
        plot.title = element_text(size = 18)) +
  guides(fill=guide_legend(nrow = 2, byrow = T)) +
  scale_fill_manual(values = c("#B34436", "#FF8E80", '#24B36C', '#66FFB4')) + 
  scale_color_manual(values = c("#B34436", "#FF8E80", '#24B36C', '#66FFB4')) +
  ggtitle("(a) Asymptotic diversity profile")

{
  svg("C:/Users/Master/FURB/André Luís de Gasper - Giesta/diversity_profile_new.svg", width = 6, height= 8)
  plot(div_profile)
  dev.off()
  }


  
#Abundance vectors of accidentals only and genuine epiphytes only by host species
acc_als=c(30,10,10,6,5,4,3,3,3,3,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
sum(acc_als)
acc_cya=c(40,23,14,13,6,5,4,4,4,4,4,3,3,3,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
sum(acc_cya)
true_als=c(253,77,38,20,9,4,3,3,1)
sum(true_als)
true_cya=c(51,51,31,28,18,3,2,1,1,1)
sum(true_cya)

mylist <- list('True epiphytes on A. setosa' = true_als, 
               'True epiphytes on C. phalerata' = true_cya, 
               'Acc. epiphytes on A. setosa' = acc_als, 
               'Acc. epiphytes on C. phalerata' = acc_cya)


#Running Hill Numbers with Cmax (2 * minimum abundance vector or max. reference sample size - whichever is larger, Chao et al. 2014)
res_mylist <-iNEXT(mylist, q = c(0,1,2), datatype = "abundance", nboot = 200, endpoint = 408)

#RAREFACTION AND EXTRAPOLATION
#Fortifying dataframe for ggplot workaround
df <- ggplot2::fortify(res_mylist, type = 1)
head(df)

df.point <- df[which(df$method=="observed"),]
df.line <- df[which(df$method!="observed"),]
df.line$method <- factor(df.line$method, 
                         c("interpolated", "extrapolated"),
                         c("interpolation", "extrapolation"))

div_plot <- ggplot(df, aes(x=x, y=y, colour=site)) +
  theme_bw() + 
  geom_point(aes(shape=site), size=1, data=df.point) +
  geom_line(aes(linetype=method), lwd=1.5, data=df.line) +
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,
                  fill=site, colour=NULL), alpha=0.2) +
  facet_wrap(~order) + 
  labs(x="Number of individuals", y="Species diversity") +
  theme(legend.position = "bottom", 
        legend.title=element_blank(),
        text=element_text(size=16),
        legend.box = "vertical",
        plot.title = element_text(size = 18)) +
  guides(fill=guide_legend(nrow = 2, byrow = T)) +
  scale_fill_manual(values = c("#B34436", "#FF8E80", '#24B36C', '#66FFB4')) + 
  scale_color_manual(values = c("#B34436", "#FF8E80", '#24B36C', '#66FFB4')) +
  ggtitle("(b) Non-asymptotic size-based rarefaction and extrapolation")

{
svg("C:/Users/Master/FURB/André Luís de Gasper - Giesta/rarext_plot.svg", width = 6, height= 8)
plot(div_plot)
dev.off()
}

#### MERGING FIGURES ####
figure_div <- ggpubr::ggarrange(div_profile, div_plot,
                            nrow = 2, ncol = 1, common.legend = F, legend = "bottom")
figure_div

{
  png("C:/Users/Master/FURB/André Luís de Gasper - Giesta/new_plots.png", width = 8, height= 10, units = "in", res = 300)
  par(mar=c(5, 1, 5, 1))
  plot(figure_div)
  dev.off()
}

##################################
#### QUASI-POISSON REGRESSION ####
##################################
#Loading and adjusting data (both genuine and accidental epiphytes)
acc <- read.csv('C:/Users/Master/FURB/André Luís de Gasper - Giesta/Biotropica_review/new_results/correlation_data_acc.csv', 
                header = TRUE, sep = ';')
names(acc) <- c("species", "humidity", "height", "cbh", "abundance", "richness")

#Obtaining radius
acc$radius = acc$cbh/3.141593
acc$radius = round(acc$radius, 2)

#Approximating to the area of a cilinder as A=2Ï€rh+2Ï€r^2
acc$area = (2*pi*acc$radius*acc$height)+(2*pi*(acc$radius^2))

#Area variables correlated?
cor(acc$int, acc$area)

#Interaction term to test later
acc$int <- acc$height*acc$humidity

#Standardizing variables
acc <- acc %>% mutate_at(c("humidity", "area", "int"), ~(scale(.) %>% as.vector))

#Since response variables are not normally distributed and data is overdispersed, 
#we use the quasi-Poisson family distribution to fit the data
#performance::check_overdispersion(acc_model)

#########################
#### ABUNDANCE MODEL ####
#########################
ab_model <- glm(formula = abundance ~ height, family = quasipoisson(link="log"), data = acc)
ab_model2 <- glm(formula = abundance ~ area + humidity + area*humidity, family = quasipoisson(link="log"), data = acc)

#Significant difference between models?
anova(ab_model, ab_model2, test = "Chisq")

#Regression parameters for abundance model
tableone::ShowRegTable(ab_model)
#Pseudo R-squared for quasi-Poisson model: 1 - (residual.deviance/null.deviance)
1 - (ab_model$deviance/ab_model$null.deviance) #0.08 partial-RÂ²

#Visualizing...

#Test for all variables separately
#New range area to predict
range(acc$int)
xrange <- seq(-1.4, 4, 0.25)

#Predicted model
yab <- predict(ab_model, list(int = xrange), type="response")
plot(acc$int, acc$abundance, pch = 16, xlab = "interaction", ylab = "ab")
lines(xrange, yab)


#########################
#### RICHNESS MODEL #####
#########################
rich_model <- glm(formula = richness ~ area + humidity + area*humidity, family = quasipoisson(link="log"), data = acc)
#Regression parameters for richness model
tableone::ShowRegTable(rich_model)
#Pseudo R-squared for quasi-Poisson model: 1 - (residual.deviance/null.deviance)
1 - (rich_model$deviance/rich_model$null.deviance) #0.08 partial-RÂ²

#Visualizing...

#Test for all variables separately
#New range area to predict
range(acc$int)
xrange <- seq(-1.4, 4, 0.25)

#Predicted model
yab <- predict(ab_model, list(int = xrange), type="response")
plot(acc$int, acc$abundance, pch = 16, xlab = "interaction", ylab = "ab")
lines(xrange, yab)



######################
### NMDS ANALYSIS ####
######################
library(vegan)
library(viridis)
library(ggplot2)

data <- read.csv('nmds_matrix_onlyacc.csv', header = TRUE, sep = ';')
spp <- data[, 2:ncol(data)]
group <- data$site

#Binarizing abundance data (Jaccard distance)
spp <- ifelse(spp > 0,1,0)
spp <- as.data.frame(spp)

#PERMANOVA to test if accidental composition differ between each host
adonis(formula = spp ~ group, permutations = 999, method = "jaccard")

#Running NMDS
nmds_results <- metaMDS(spp, k = 2, maxit = 999, trymax = 500, wscore = TRUE, distance = "jaccard")

#Site data from NMDS object
data_scores <- as.data.frame(scores(nmds_results))

# Now add the extra aquaticSiteType column
data_scores <- cbind(data_scores, group)
colnames(data_scores)[3] <- "Host"

#Creating the hull object
grp.a <- data_scores[data_scores$Host == "A. setosa", ][chull(data_scores[data_scores$Host == 
                                                                   "A. setosa", c("NMDS1", "NMDS2")]), ] #Hull for A. setosa
grp.b <- data_scores[data_scores$Host == "C. phalerata", ][chull(data_scores[data_scores$Host == 
                                                                   "C. phalerata", c("NMDS1", "NMDS2")]), ]  #Hull for C. phalerata

hull.data <- rbind(grp.a, grp.b) #combine grp.a and grp.b
hull.data

#Plotting
nmds_plot <- ggplot() +
  geom_point(data = data_scores, aes(x = NMDS1, y = NMDS2, 
                                     color = Host), size = 2) +
  geom_polygon(data = hull.data, aes(x = NMDS1, y = NMDS2, group = Host), alpha = 0.15) + 
  scale_color_manual(values = inferno(15)[c(3, 8)],
                     name = "Host") +
  coord_equal() + 
  annotate(geom = "label", x = 0, y = -1.5, size = 8,
           label = paste("Stress: ", round(nmds_results$stress, digits = 3))) +
  theme_bw() +
  theme(axis.text = element_text(size = 20), 
        legend.position = "bottom",
        text = element_text(size = 18))


png("C:/Users/ggrit/FURB/Andr? Lu?s de Gasper - Giesta/Biotropica_review/new_results/nmds.svg", width = 10, height= 8)
plot(nmds_plot)
  dev.off()

##################

#Stratification graph: the distribution of accidentals vs true epiphytes throughout strata

library(data.table)

data = read.csv('new_results/strat_graph_biotropica.csv', header = TRUE, sep = ';')
#data$class = factor(data$class)
setDT(data)
data[int%in% "A", int := "A (0-1m)"]
data[int%in% "B", int := "B (1-2m)"]
data[int%in% "C", int := "C (2-3m)"]

colorset = c("Non-Acc" = "yellowgreen", "Acc" = "darkorange")

abun_plot <- ggplot(data) +
  geom_bar(aes(x = class, y = abun, fill = class), stat = "identity") + facet_wrap(~int) + theme_bw() +
  scale_fill_manual(values=colorset) + guides(fill=F) + xlab("\nEcological category") +ylab("Abundance\n") +
  theme(axis.text = element_text(size = 18), text = element_text(size = 20), axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

#Converting abundance to richness
data$richness = data[, c("abun")]
data$richness[data$richness > 0] <- 1 #ifelse() works too

rich_plot <- ggplot(data) +
  geom_bar(aes(x = class, y = richness, fill = class), stat = "identity") + facet_wrap(~int) + theme_bw() +
  scale_fill_manual(values=colorset) + guides(fill=F) + xlab("\nEcological category") +ylab("Richness\n") +
  theme(axis.text = element_text(size = 18), text = element_text(size = 20), strip.background = element_blank(),
        strip.text.x = element_blank())

final_plot = ggarrange(abun_plot, rich_plot, ncol = 1, nrow = 2)

final_plot
ggsave('new_results/strat_graph.jpg', width = 10, height = 10)
dev.off()


##################################
#### HEIGHT-INTERVAL ANALYSIS ####
##################################
data = read.csv('strat_graph_biotropica.csv', header = TRUE, sep = ';')
#data$class = factor(data$class)
setDT(data)
data[int%in% "A", int := "Lower (0-1m)"]
data[int%in% "B", int := "Medium (1-2m)"]
data[int%in% "C", int := "Upper (2-3m)"]

data = subset(data, class == "Acc")

data = subset(data, abun > 0)

# Build the linear model
model  <- lm(abun ~ int, data = data)
# Create a QQ plot of residuals
ggqqplot(residuals(model))

res.kruskal <- data %>% kruskal_test(abun ~ int)
res.kruskal

data %>% kruskal_effsize(abun ~ int)

# Pairwise comparisons
pwc <- data %>% 
  dunn_test(abun ~ int, p.adjust.method = "bonferroni") 
pwc

# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "int")
ggboxplot(data, x = "int", y = "abun") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc)
  ) + ylab("Abundance") + xlab("Height class") + theme_bw()


ggsave("C:/Users/Master/FURB/André Luís de Gasper - Giesta/Biotropica_review/new_results/kruskal_newnames.png", height = 10, width = 10)

rm(list=ls())
