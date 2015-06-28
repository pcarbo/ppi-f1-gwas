# This file contains a list which holds the default settings for each
# phenotype:
#
#   (a) tranformation function
#   (b) covariate names
#   (c) outlier removal function
#   (d) specified outliers
#
model.info <- list(
  
  # Different covariates between studies (using CAC covariates).
  ppi12 = list(
    transformation   = function (x) logit10(project.onto.interval(x,1,99)/100),
    covariates       = c("ppibox3","ppibox5"),
    outlier.function = NULL,
    outliers         = c("53524","53261","53128","55182","53404","53401",
                         "55693","53562","55129","53113","55080","55095",
                         "55304","55603","55462","55541","50923")),
)
