# This file contains a list which holds the default settings for each
# phenotype:
#
#   (a) tranformation function
#   (b) covariate names
#   (c) outlier removal function
#   (d) specified outliers
#
model.info <- list(

  ppi3 = list(
    transformation = function (x) logit10(project.onto.interval(x,1,99)/100),
    covariates = c("ppibox3","ppibox5"),
    outlier.function = function(x) x < (-1),
    outliers = c("53204","55129","55182","55325","55766","55762","50918",
                 "50988","51075","55332","53155","55591","53569","55222",
                 "55693","53191","55750","53671","55318","55005","53577",
                 "53480","55462","55630","55460","55541","53680","55061")),
  
  ppi6 = list(
    transformation = function (x) logit10(project.onto.interval(x,1,99)/100),
    covariates = c("ppibox3", "ppibox5"),
    outlier.function = function(x) x < (-1),
    outliers = c("53110","55129","55348","55766","55770","50916","50988",
                 "53524","53526","55363","53585","53213","55591","55564",
                 "57057","55262","53401","50999","55693","57034","55750",
                 "55068","53549","53116","55674","53421","55304","55603",
                 "55095","55166","55541","55061")),
  
  ppi12 = list(
    transformation   = function (x) logit10(project.onto.interval(x,1,99)/100),
    covariates       = c("ppibox3","ppibox5"),
    outlier.function = NULL,
    outliers         = c("53524","53261","53128","55182","53404","53401",
                         "55693","53562","55129","53113","55080","55095",
                         "55304","55603","55462","55541","50923"))
)
