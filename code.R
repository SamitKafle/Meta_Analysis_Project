#load library
library(metafor)

#load data
allom=read.csv("allometric_effectsize14.csv")
tls=read.csv("tls_effectsize14.csv")

# Calculate effect sizes: log response ratio for TLS vs DEST
es_tls <- escalc(measure = "ROM",
                 m1i = TLS_mean, sd1i = TLS_SD, n1i = n,
                 m2i = DEST_mean, sd2i = DEST_SD, n2i = n,
                 data = tls)

# Calculate effect sizes: log response ratio for ALLOM vs DEST
es_allom <- escalc(measure = "ROM",
                   m1i = ALLOM_mean, sd1i = ALLOM_SD, n1i = n,
                   m2i = DEST_mean, sd2i = DEST_SD, n2i = n,
                   data = allom)
# Meta-analysis TLS vs DEST
res_tls <- rma(yi, vi, data = es_tls, method = "REML")
summary(res_tls)

# Meta-analysis ALLOM vs DEST
res_allom <- rma(yi, vi, data = es_allom, method = "REML")
summary(res_allom)


# Forest plots
forest(res_tls, slab = tls$Author, xlab = "Log Response Ratio (TLS vs DEST)", main = "")
forest(res_allom, slab = allom$Author, xlab = "Log Response Ratio (ALLOM vs DEST)", main = "")


# Funnel plots for publication bias
funnel(res_tls, main = "Funnel Plot: TLS vs DEST")
funnel(res_allom, main = "Funnel Plot: ALLOM vs DEST")
# Step 4: Add Egger's test p-value to the plot
pval <- formatC(reg_tls$pval, format = "f", digits = 3)  # Format p-value to 3 decimal places
text(x = 0.6, y = 0.4, labels = paste0("Egger's p = ", pval), adj = c(0,1), cex = 1)
# Egger's test for asymmetry
reg_tls=regtest(res_tls, model = "lm")
reg_allom=regtest(res_allom, model = "lm")
reg_tls
reg_allom


###############################################################################

# Moderator  (species type)
res_mod_allom <- rma(yi, vi, mods = ~ type+equation_type, data = es_allom,method = "REML")
summary(res_mod_allom)

# Moderator ( equation type)
res_mod_allom2 <- rma(yi, vi, mods = ~ equation_type-1, data = es_allom, method = "REML")
res_mod_allom2



# load orchard package
pacman::p_load(devtools, tidyverse, metafor, patchwork, R.rsp, orchaRd, emmeans, ape, phytools, flextable)

# we can create a table of results for moderator
res2 <- orchaRd::mod_results(res_mod_allom, mod = "type", group = "ID")
res2

p2 <- orchaRd::orchard_plot(res2,
                            mod = "type", group = "Author", xlab = "lnRR")##tree type as moderator
p2


p6_per <- orchaRd::orchard_plot(res_mod_allom, mod = "type", group = "Author", 
                                xlab = "Bias (%)", 
                                angle = 0, N = "n", g = FALSE, transfm = "percentr") #plot with log ratio transformed to percentage

p6_per


res3 <- orchaRd::mod_results(res_mod_allom2, mod = "equation_type", group = "Author")#equation type as moderator
res3

p3 <- orchaRd::orchard_plot(res3,
                            mod = "equation_type", group = "Author", xlab = "lnRR")
p3

p7_per <- orchaRd::orchard_plot(res_mod_allom2, mod = "equation_type", group = "Author", 
                                xlab = "Bias (%)", 
                                angle = 0, N = "n", g = FALSE, transfm = "percentr") 

p7_per




##### Moderator ( dbh)
res_mod_allom3 <- rma(yi, vi, mods = ~ mean_dbh, data = es_allom, method = "REML")
res_mod_allom3
res_mod_tls <-rma(yi, vi, mods = ~ mean_dbh, data = es_tls, method = "REML")
res_mod_tls

#,
##weights = "prop", by = "Environment"
lim_bubble <- orchaRd::mod_results(res_mod_allom3, mod = "mean_dbh", group = "ID")
orchaRd::bubble_plot(lim_bubble, group = "Article",  mod = "mean_dbh", xlab = "mean_dbh", legend.pos = "top.left")


lim_bubble <- orchaRd::mod_results(res_mod_tls, mod = "mean_dbh", group = "ID")
orchaRd::bubble_plot(lim_bubble, group = "Article",  mod = "mean_dbh", xlab = "mean_dbh", legend.pos = "top.left")

##leave one out test
loo_test_allom <- leave_one_out(res_allom, group = "Author")##for allometric
loo_test_tls <- leave_one_out(res_tls, group = "Author")##for tls
#plot
orchard_leave1out(leave1out = loo_test_allom,
                  xlab = "lnRR",
                  ylab = "Study left out",
                  trunk.size = 1.2,
                  branch.size = 1.5,
                  alpha = 0.08,
                  legend.pos = "top.out")

orchard_leave1out(leave1out = loo_test_tls,
                  xlab = "lnRR",
                  ylab = "Study left out",
                  trunk.size = 1.2,
                  branch.size = 1.5,
                  alpha = 0.08,
                  legend.pos = "top.out")

###################################################caterpillar plots
orchaRd::caterpillars(res_mod_allom, mod="1", xlab = "lnRR",group="ID")##allometric vs destructive
orchaRd::caterpillars(res_mod_tls, mod="1", xlab = "lnRR",group="ID")##tls vs destructive



###################3

