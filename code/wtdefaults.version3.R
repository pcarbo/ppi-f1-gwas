# This file contains a list which holds the default settings for each phenotype's:
#   (a) tranformation function
#   (b) covariate names
#   (c) outlier removal function
#   (d) specified outliers

model.info <- list(
#   
#   # Not analyzed separately in either
#   d50bw = list(
#     transformation = NULL,
#     covariates = c("sex", "littersize"),
#     outlier.function = NULL,
#     outliers = NULL),
#   
#   d100bw = list(
#     transformation = NULL,
#     covariates = c("sex", "littersize", "study"),
#     outlier.function = NULL,
#     outliers = NULL),
#   
#   weightgain = list(
#     transformation = NULL,
#     covariates = c("study"),
#     outlier.function = NULL,
#     outliers = c("53213")),
#   
#   d1totalactivity = list(
#     transformation = function (x) log10(x),
#     covariates = "study",
#     outlier.function = NULL,
#     outliers = c("50931", "55689")),
#   
#   d1centerdur = list(
#     transformation = function (x) logit10(project.onto.interval(x / 400, 0.001, 0.999)),
#     covariates = c("sex"),
#     outlier.function = function (x) x > 1.5,
#     outliers = c("50931", "53273", "55523", "55689")),
#   
#   fstdist = list(
#     transformation = function (x) log10(x),
#     covariates = c("agouti", "study"),
#     outlier.function = NULL,
#     outliers = NULL),
#   
#   # Not analyzed separately in either
#   highmobfreq = list(
#     transformation = NULL,
#     covariates = NULL,
#     outlier.function = NULL,
#     outliers = NULL),
#   
#   # Not analyzed separately in TCF
#   highmobdur = list(
#     transformation = function (x) log10(x + .001),
#     covariates = c("agouti", "d100bw"),
#     outlier.function = NULL,
#     outliers = NULL),
#   
#   # Not analyzed separately in either
#   mobfreq = list(
#     transformation = NULL,
#     covariates = NULL,
#     outlier.function = NULL,
#     outliers = NULL),
#   
#   # Not analyzed separately in either
#   mobdur = list(
#     transformation = NULL,
#     covariates = NULL,
#     outlier.function = NULL,
#     outliers = NULL),
#   
#   # Different transformation and covariates between studies
#   immobfreq = list(
#     transformation = NULL,
#     covariates = NULL,
#     outlier.function = NULL,
#     outliers = NULL),
#   
#   # Not analyzed separately in TCF
#   immobdur = list(
#     transformation = function (x) log10(pmax(0.01, x)),
#     covariates = c("agouti", "study"),
#     outlier.function = function (x) x < (-1.9),
#     outliers = c("55232", "55363", "55381", "55382", "55613")),
#   
#   ppi3 = list(
#     transformation = function (x) logit10(project.onto.interval(x, 1, 99) / 100),
#     covariates = c("ppibox3"),
#     outlier.function = function(x) x < (-1),
#     outliers = c("53204", "55129", "55182", "55325", "55766", "55762", "50918", "50988", "51075",
#                  "55332", "53155", "55591", "53569", "55222", "55693", "53191", "55750", "53671","55318",
#                  "55005", "53577", "53480", "55462", "55630", "55460", "55541", "53680", "55061")),
#   
#   ppi6 = list(
#     transformation = function (x) logit10(project.onto.interval(x, 1, 99) / 100),
#     covariates = c("ppibox3", "ppibox5"),
#     outlier.function = function(x) x < (-1),
#     outliers = c("53110", "55129", "55348", "55766", "55770", "50916", "50988", "53524", "53526",
#                  "55363", "53585", "53213", "55591", "55564", "57057", "55262", "53401", "50999",
#                  "55693", "57034", "55750", "55068", "53549", "53116", "55674", "53421", "55304",
#                  "55603", "55095", "55166", "55541", "55061")),
  
  # Different covariates between studies (using CAC covariates)
  ppi12 = list(
    transformation = function (x) logit10(project.onto.interval(x, 1, 99) / 100),
    covariates = c("ppibox3", "ppibox5"),
    outlier.function = NULL,
    outliers = c("53524", "53261", "53128", "55182", "53404", "53401", "55693", "53562", "55129",
                 "53113", "55080", "55095", "55304", "55603", "55462", "55541", "50923"))
#   
#   # Different covariates between studies (using CAC covariates)
#   ppiavg = list(
#     transformation = function (x) logit10(project.onto.interval(x, 1, 99) / 100),
#     covariates = c("ppibox3", "ppibox5"),
#     outlier.function = function(x) x < (-1),
#     outliers = c("50967", "53110", "55129", "55766")),
#   
#   # Different covariates between studies (using TCF covariates)
#   startle = list(
#     transformation = function (x) log10(x),
#     covariates = c("sex", "d50bw", "ppibox1", "ppibox2", "ppibox3", "ppibox4"),
#     outlier.function = function(x) x < (-1),
#     outliers = NULL),
#   
#   p120b4 = list(
#     transformation = function (x) log10(x),
#     covariates = c("p120b1", "sex", "d50bw", "ppibox1", "ppibox2", "ppibox3", "ppibox4"),
#     outlier.function = function (x) x < (-0.7),
#     outliers = c("53106", "55865", "57004")),
#   
#   # Different covariates used between studies (left NULL)
#   avcontextd2 = list(
#     transformation = function (x) logit10(project.onto.interval(x, 1, 99) / 100),
#     covariates = NULL,
#     outlier.function = NULL,
#     outliers = NULL),
#   
#   # Different covariates used between studies, using union of both sets
#   avaltcontextd3 = list(
#     transformation = function (x) logit10(project.onto.interval(x, 1, 99) / 100),
#     covariates = c("study"),
#     outlier.function = NULL,
#     outliers = c("50970", "53595", "52409")),
#   
#   # Different covariates used between studies, using union of both sets
#   avtoned3 = list(
#     transformation = function (x) logit10(project.onto.interval(x, 1, 99) / 100),
#     covariates = c("sex", "study"),
#     outlier.function = NULL,
#     outliers = NULL),
#   
#   # Not analyzed separately in either study
#   d1fecalboli = list(
#     transformation = function (x) qt.random.tie(x),
#     covariates = NULL,
#     outlier.function = NULL,
#     outliers = NULL),
#   
#   # Not analyzed separately in either study
#   d2fecalboli = list(
#     transformation = function (x) qt.random.tie(x),
#     covariates = NULL,
#     outlier.function = NULL,
#     outliers = NULL),
#   
#   # Not analyzed separately in either study
#   d3fecalboli = list(
#     transformation = function (x) qt.random.tie(x),
#     covariates = NULL,
#     outlier.function = NULL,
#     outliers = NULL)
#   
)
