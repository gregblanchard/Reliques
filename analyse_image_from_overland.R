

library(caTools)

library(raster)

#### get beta div map ####


dirname <- "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/zones_test"
filename <- file.path(dirname, "BetaDiversity_BCdiss_PCO_10_Fullres")
BetaDiversity <- raster(filename)
# plot(BetaDiversity)

# filename <- file.path(dirname, "Shannon_10_Fullres")
# Shannon <- raster(filename)


#### get sampling points ####

points_spl_utm <- readRDS(
 "/home/thesardfou/Documents/projets/Reliques/signature_spectrale_fragmentation/maps/beta_div_spectrale_NC/data_frag/points_spl_utm_landscape_metrics_topo.rds")

#### get beta_div values for sampling points ####

# select points in same extant 
points_spl_utm_zone <- crop(points_spl_utm, extent(BetaDiversity))

beta_div_pts <- extract(BetaDiversity, points_spl_utm_zone)

points_spl_utm_zone$beta_div <- beta_div_pts


#### analyses ####
points_spl_utm_zone$log_edge <- log(points_spl_utm_zone$dist_edge+1)
points_spl_utm_zone$beta_div01 <- 
  (points_spl_utm_zone$beta_div - min(points_spl_utm_zone$beta_div, na.rm = T)) / (max(points_spl_utm_zone$beta_div, na.rm = T) - min(points_spl_utm_zone$beta_div, na.rm = T))

data <- as.data.frame(points_spl_utm_zone)
dim(data)
data <- data[!is.na(data$beta_div),]
dim(data)
colnames(data)

lm_edge <- lm(beta_div ~ dist_edge , data = data)
summary(lm_edge)
lm_edge <- lm(beta_div ~ log_edge , data = data)
summary(lm_edge)

sub_data <- data[sample(3000),]
lm_edge <- lm(beta_div ~ log_edge , data = sub_data)
summary(lm_edge)
hist(lm_edge$residuals)

shapiro.test(lm_edge$residuals)
lmtest::bptest(lm_edge)

hist(data$beta_div)
hist(data$dist_edge)
hist(data$log_edge)

plot(beta_div ~ log_edge , data = data)

hist((data$beta_div)^1.8)

####################################################
#### model selection ####
####################################################
library(MuMIn)

##### simple model selection for landscape buffers ####
msel_aic_uni = function(Y,X, empty = T){
  aic = c()
  for(i in 1:ncol(X)){
    lm = lm (Y~X[,i])
    aic = c(aic,AIC(lm))
    if(empty & i == ncol(X)) aic = c(aic, AIC(lm(Y~1)))
  }
  if(empty) tab = cbind(c(colnames(X), "intercept"), aic) else  tab = cbind(colnames(X), aic)
  return(list(pred = colnames(X)[order(aic)][1], tab = tab))
}  
  # habitat amount
  aic_sel_HA <- msel_aic_uni(data_for_mod$beta_div01, data_for_mod[,c( paste0(rep("total.area_", 4), c(250,500,1000,2000)))])
  aic_sel_HA
  
  # effective mesh size
  aic_sel_EMS <- msel_aic_uni(data_for_mod$beta_div01, data_for_mod[,c( paste0(rep("effective.mesh.size_", 4), c(250,500,1000,2000)))])
  aic_sel_EMS
  
  
colnames(data)
land_var <- paste0(rep(c("total.edge_", "total.area_", "effective.mesh.size_"), each = 4), c(250,500,1000,2000))
# from selection
land_var <- c(aic_sel_HA$pred, aic_sel_EMS$pred)

topo_var <- c("elevation","slope","aspect" ,"curvature","substrat")
sel_var <- c("log_edge",land_var, topo_var)
data_for_mod <- data[, c("beta_div01", sel_var)]

# no NA values
data_for_mod <- data_for_mod[complete.cases(data_for_mod),]

# no points < d m from edge
d =0
data_for_mod <- data_for_mod[data_for_mod$log_edge > log(d),]

# standardize variables 
data_for_mod <- as.data.frame(apply(data_for_mod, 2, function(x) if(!any(is.numeric(x))){
  scale(as.numeric(as.factor(x)))
}else{
  scale(is.numeric(x))
}))
# Variance Inflation Factors

glm1 <- glm(beta_div01~ . , data = data_for_mod, na.action = "na.fail")
car::vif(glm1)
library(parallel)
cl <- makeCluster(3)
clusterExport(cl, "data_for_mod")
dd <- pdredge(glm1,cluster=cl,rank = "AIC", beta = "partial.sd", extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
      F = s$fstatistic[[1]])
  })
  )
stopCluster(cl)

dd
# selection few best models based on delta AIC
models.list <- get.models(dd,subset =  delta < 2) 

# average models and coeficiants based on weigth
muminavg <- model.avg(models.list, revised.var = TRUE, fit = TRUE, beta = "partial.sd")
summary(muminavg)

# transform AIC-weigthed standardized coeficients
weighted_coef <- coef(muminavg)
sign_coef <- weighted_coef
# transform to relative influence 
relat_inf <- abs(weighted_coef)/sum(abs(weighted_coef))
relat_inf <- relat_inf[2:length(relat_inf)]

# get conf. intervals for conditional average (selected models) 
weighted_confint <- confint(muminavg)

#### plot ####
library(GGally)

df_coef <- data.frame(term = names(weighted_coef[2:length(weighted_coef)]),
                      estimate = weighted_coef[2:length(weighted_coef)],
                      conf.low = weighted_confint[2:length(weighted_coef),1],
                      conf.high = weighted_confint[2:length(weighted_coef),2])

ggcoef(df_coef, sort = "ascending", vline_linetype =  "dotted", mapping = aes(x = estimate, y = term, colour = term),
       errorbar_height = .1) + 
  theme_bw() + 
  theme(legend.position="none")

library(gridExtra)

coef_plot <- ggcoef(df_coef, sort = "ascending", vline_linetype =  "dotted", mapping = aes(x = estimate, y = term, colour = term),
                    errorbar_height = .1) + 
  theme_bw() + 
  theme(legend.position="none")


df <- data.frame(relat_inf = relat_inf,term = names(relat_inf))
df <- df[order(df_coef$estimate, decreasing = F),]

relat_inf_plot <- ggplot(data=df, aes(x=term, y=relat_inf, fill = term)) +
  geom_bar(stat="identity")+
  scale_colour_gradient2()+
  coord_flip()+
  ylim(0, max(df$relat_inf)*1.3)+
  scale_x_discrete(limits = df$term)+
  theme_classic() + 
  theme(legend.position="none")


grid.arrange(coef_plot, relat_inf_plot, ncol=2)


#### simple plot ####
plot(data_for_mod$beta_div01~data_for_mod$total.area_2000)

#### best model ####
models.best <- get.models(dd,subset =  delta == 0)[[1]] 
summary(models.best)
r.sq <- r.squaredLR(models.best)
r.sq
plot(models.best)

#### test with glm ####
glm_edge <- glm(beta_div01 ~ log_edge, data = sub_data, family = "inverse.gaussian")
summary(glm_edge)
plot(glm_edge)
library(glmmTMB)    
glm_edge <- glmmTMB(beta_div01 ~ log_edge + (1|substrat), data = sub_data, family=beta_family(link="logit"))
summary(glm_edge)

#### test with glmnet ####

library(glmnet)

colnames(data)

land_var <- paste0(rep(c("total.edge_", "total.area_", "effective.mesh.size_"), each = 4), c(250,500,1000,2000))

topo_var <- c("elevation","slope","aspect" ,"curvature","substrat")
sel_var <- c("log_edge",land_var, topo_var)
expl_var <- data[, sel_var]

# no NA values
keep_pts <- complete.cases(expl_var) & complete.cases(data$beta_div01) 

# transform to numeric matrix 
expl_var <- as.matrix(expl_var[keep_pts,])
expl_var <- apply(expl_var, 2, function(x) if(!any(is.numeric(x))){
  as.numeric(as.factor(x))
}else{
  is.numeric(x)
})


beta_div01 <- data$beta_div01[keep_pts]

fit = glmnet(expl_var, beta_div01)
fit
summary(fit)
plot(fit)
coef(fit,s=0.05)
# cross validation 
cvfit = cv.glmnet(expl_var, beta_div01)
plot(cvfit)
cvfit$lambda.min
coef(cvfit, s = "lambda.min")

small.lambda.index <- which(cvfit$lambda == cvfit$lambda.min)
small.lambda.index <- which(cvfit$lambda == cvfit$lambda.1se)

small.lambda.betas <- cvfit$glmnet.fit$beta[, small.lambda.index]
small.lambda.betas
