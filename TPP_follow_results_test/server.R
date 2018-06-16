library(shiny)
library(shinydashboard)





######################################################################
#                                                                    #
#                           Load packages                            #
#                                                                    #
######################################################################

library(HDInterval) # needed to calcluate HDI credible intervals in the Bayesian parameter estimation robustness test
library(ggplot2) # for visualization
library(emdist) # to calcluate earth mover's distance (EMD)

######################################################################
#                                                                    #
#                            Functions                               #
#                                                                    #
######################################################################

######################################################################
#                  Bayes factor calculation functions                #
######################################################################


### Functions for Bayes factor caclulation using beta prior
# These functions are required to run the Bayes factor analysis 
# we thank Richard Morey for his help in developing these functions


fullAlt_beta = Vectorize(function(p, y, N, alpha, beta){
  exp(dbinom(y, N, p, log = TRUE) + dbeta(p, alpha, beta, log = TRUE)) 
},"p")

normalize_beta = function(alpha, beta, interval){
  diff(pbeta(interval, alpha, beta))
}

restrictedAlt_beta = function(p,y,N,y_prior,N_prior,interval){
  alpha = y_prior + 1
  beta = N_prior - y_prior + 1
  fullAlt_beta(p, y, N, alpha, beta) / normalize_beta(alpha, beta, interval) * (p>interval[1] & p<interval[2])
}

margLike_beta = function(y, N, y_prior, N_prior, interval){
  integrate(restrictedAlt_beta, interval[1], interval[2], 
            y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval)[[1]]
}

BF01_beta = Vectorize(function(y, N, y_prior, N_prior, interval, null_prob){
  dbinom(y, N, null_prob) / margLike_beta(y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval)
},"y")


######################################################################
#                       Other supporting functions                   #
######################################################################

### Function calculating the highest density interval using sampling
# We use hdi() from the library(HDInterval)
# this function is needed for the Bayesian parameter estimation robustness test

mode_HDI <- function(scale, density, crit_width = 0.95, n_samples = 1e5){
  samp <- sample(x = scale, size = n_samples, replace = TRUE, prob = density)
  hdi_result = hdi(samp, credMass=crit_width)
  result = c(scale[which(density == max(density))], # mode
             hdi_result[1], # lower bound
             hdi_result[2]) # upper bound
  
  # only needed for the names of the result
  Crit_lb = (1-crit_width)/2
  Crit_ub = crit_width + (1-crit_width)/2
  
  names(result) = c("mode", paste(Crit_lb*100, "%", sep = ""), paste(Crit_ub*100, "%", sep = ""))
  return(result)
}




######################################################################
#                                                                    #
#                          Data parameters                           #
#                                                                    #
######################################################################

### analysis parameters
# these are the analysis parameters currently specified in our protocol

# number of erotic trials performed by participants in the study (if no missing trials)
trial_size_per_participant = 18

# probability of success if M0 is true
M0_prob = 0.5

# interim analysis points (in total number of erotic trials performed)
when_to_check = c(30060, 37836, 45612, 53406, 61182, 68958, 76734, 84528, 92304, 100080)

# thresholds to infer support for M0 (high) or M1 (low)
Inference_threshold_BF_high = 25
Inference_threshold_BF_low = 1/Inference_threshold_BF_high

# this information is used both for calculating replication Bayes factor, and the Bayesian parameter estimation robustness test. 
# Here we use data from Bem's experiment 1, 828 successes within 1560 erotic trials
y_prior = 828 #number of successes in erotic trials in Bem's experiment 1
N_prior = 1560 # number of erotic trials in Bem's experiment 1


# smallest effect size of interest in the NHST equivalence test 
minimum_effect_threshold_NHST = 0.01
# p threshold for the NHST proportion test robustness test
Inference_threshold_robustness_NHST = 0.005

# in the Bayesian parameter estimation robustness test this will determine the region of practical 
#equivalence (ROPE) interval. The ROPE is interpreted similarly to SESOI, but not entireli the same. 
# See Kruschke, J. K., & Liddell, T. M. (2017). The Bayesian New Statistics: Hypothesis testing, 
# estimation, meta-analysis, and power analysis from a Bayesian perspective. 
# Psychonomic Bulletin & Review, 1-29. 
minimum_effect_threshold_Bayes_Par_Est = 0.006
# this threshold is used to set the HDI width to check against the ROPE in the Bayesian parameter 
# estimation robustness test, if ths parameter is set to 0.05 for example, it means that we would 
# expect that 95% of the probability mass would be within the ROPE to accept a hypothesis
Inference_threshold_robustness_Bayes_Par_Est = 0.05 












######################################################################
#                                                                    #
#                         Data management                            #
#                                                                    #
######################################################################

# get up to date date from github
pilot_data_pre <- read.csv("https://raw.githubusercontent.com/gy0p4k/transparent-psi-results/master/results.csv")

# only read pilot sessions
# pilot sessions are marked with session_type = "pilot" and they were all collected in the lab with the laboratory_ID_code = "lab_ELTE_01"
lab_ID <- "lab_ELTE_01"
pilot_data <- pilot_data_pre[pilot_data_pre[,"session_type"] == "pilot" & pilot_data_pre[,"laboratory_ID_code"] == "lab_ELTE_01", ]


# Number of participants tested in the pilot test
sample_size_participants = length(unique(pilot_data[, "participant_ID"]))


pilot_data[, "trial_number"] = as.numeric(pilot_data[, "trial_number"])

#extract data from erotic trials 
data_BF = pilot_data[!is.na(pilot_data[, "trial_number"]) & pilot_data[, "reward_type"] == "erotic", ]

data_BF[,"participant_ID"] = droplevels(data_BF[,"participant_ID"])


######################################################################
#                                                                    #
#                          Data analysis                             #
#                                                                    #
######################################################################


######################################################################
#                       Primary confirmatory test                    #
######################################################################

#================================================================#
#                Calculate number of successes                   #
#================================================================#

# number of successes and total N of trials
sides_match = data_BF[,"guessed_side"] == data_BF[,"target_side"]
successes = sum(sides_match)
total_N = nrow(data_BF)


#================================================================#
#        Calculating Bayes factors using different priors        #
#================================================================#

### Replication Bayes factor, with the Bem 2011 experiment 1 results providing the prior information

BF_replication <- BF01_beta(y = successes, N = total_N, y_prior = y_prior, N_prior = N_prior, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null


### Bayes factor with uniform prior
# using a non-informative flat prior distribution with alpha = 1 and beta = 1

BF_uniform <- BF01_beta(y = successes, N = total_N, y_prior = 0, N_prior = 0, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null


### Bayes factor with BUJ prior
# the BUJ prior is calculated from Bem's paper where the prior distribution is defined as a
# normal distribution with a mean at 0 and 90th percentele is at medium effect size d = 0.5 
# (we asume that this is one-tailed). Source: Bem, D. J., Utts, J., & Johnson, W. O. (2011). 
# Must psychologists change the way they analyze their data? Journal of Personality and Social Psychology, 101(4), 716-719.
# We simulate this in this binomial framework with a one-tailed beta distribution with alpha = 7 and beta = 7.
# This distribution has 90% of its probability mass under p = 0.712, which we determined 
# to be equivalent to d = 0.5 medium effect size. We used the formula to convert d to log odds ratio logodds = d*pi/sqrt(3), 
# found here: Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R. (2009). 
# Converting Among Effect Sizes. In Introduction to Meta-Analysis (pp. 45-49): John Wiley & Sons, Ltd.
# Then, log odds ratio vas converted to probability using the formula: p = exp(x)/(1+exp(x))
# The final equation: exp(d*pi/sqrt(3))/(1+exp(d*pi/sqrt(3)))

BF_BUJ <- BF01_beta(y = successes, N = total_N, y_prior = 6, N_prior = 12, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null


#================================================================#
#                    Main analysis inference                     #
#================================================================#

# determine inference (supported model) based on the Bayes factors calculated above  
if(Inference_threshold_BF_low >= max(c(BF_replication, BF_uniform, BF_BUJ))) {
  inference_BF = "M1"
  break} else if(Inference_threshold_BF_high <= min(c(BF_replication, BF_uniform, BF_BUJ))) {
    inference_BF = "M0"
    break} else {inference_BF = "Inconclusive"}


##################################################
#                 Robustness tests               #
##################################################

#================================================================#
#               Robustness test of BF results with NHST          #
#================================================================#

# robustness of BF results is tested with NHST proportion tests
# here we perform both an equality test and an equivalence test to draw statistical inference


# equality test
equality_test_p = prop.test(x = successes, n = total_N, p = M0_prob, alternative = "greater")$p.value

# equivalence test
equivalence_test_p = prop.test(x = successes, n = total_N,
                               p = M0_prob+minimum_effect_threshold_NHST, 
                               alternative = "less")$p.value

# descritives of the frequentist estimate
pbar = successes/total_N
SE = sqrt(pbar * (1-pbar)/total_N)
E = qnorm(.9975) * SE
proportion_995CI = round(pbar + c(-E , E), 3)



# making inference decision
if((Inference_threshold_robustness_NHST > equality_test_p) & 
   (Inference_threshold_robustness_NHST <= equivalence_test_p)){
  inference_robustness_NHST = "M1"} else if((Inference_threshold_robustness_NHST > equivalence_test_p) & 
                                            (Inference_threshold_robustness_NHST <= equality_test_p)){
    inference_robustness_NHST = "M0" 
  } else {inference_robustness_NHST = "Inconclusive"}


#=======================================================================#
#   Robustness test of BF results with Bayesian parameter estimation    #
#=======================================================================#

# robustness of BF results is tested by calculating HDI of the posteriro distribution and checking its relation to
# the region of practical equivalence (ROPE), promoted in Kruschke, J. K., & Liddell, T. M. 
# (2017). The Bayesian New Statistics: Hypothesis testing, estimation, meta-analysis, and power 
# analysis from a Bayesian perspective. Psychonomic Bulletin & Review, 1-29. 

# calculate posterior distribution using beta distribution
prior_alpha = y_prior + 1
prior_beta = N_prior-y_prior+1

posterior_alpha = prior_alpha + successes
posterior_beta = prior_beta + total_N - successes

scale = seq(0, 1, length = 1001)
posterior_density = dbeta(scale, posterior_alpha, posterior_beta)

# calculate HDI for the posterior distribution
# (here we calculate the upper and lower bound of the 90% of the probability mass
# because we use a one-tailed test. This means that the total probability mass below
# the upper bound of the 90% HDI will be 95%)
hdi_result = mode_HDI(scale = scale, density = posterior_density, crit_width = 1-Inference_threshold_robustness_Bayes_Par_Est*2, n_samples = 1e6)

# parameters for decision making
HDI_lb = hdi_result[2]
HDI_ub = hdi_result[3]

# making inference decision
ROPE = M0_prob+minimum_effect_threshold_Bayes_Par_Est
if(HDI_lb >= ROPE){inference_robustness_Bayes_Par_Est = "M1"
} else if(HDI_ub <= ROPE){inference_robustness_Bayes_Par_Est = "M0"
} else {inference_robustness_Bayes_Par_Est = "Inconclusive"}


#=======================================================================#
#      Determine final inference of all robustness tests combined       #
#=======================================================================#

# the main analysis inference is only robust if all robustness tests came to the same inference as the main test
inferences = c(inference_robustness_NHST, inference_robustness_Bayes_Par_Est)
inference_robustness = if(all(inferences == inferences[1])){inferences[1]} else {"mixed"}

Robust = if(inference_BF == inference_robustness){"robust"} else {"not robust"}


######################################################################
#                                                                    #
#                   Interpretation and visualization                 #
#                                                                    #
######################################################################


# visualize results
# ALL THREE BF values need to cross the threshold for the interpretation to be accepted

BF_results_for_plotting = cbind(as.data.frame(c(BF_replication, BF_uniform, BF_BUJ)), c("BF_replication", "BF_uniform", "BF_BUJ"))
names(BF_results_for_plotting) = c("Bayes_factor_01", "BF_type")

plot <- ggplot(BF_results_for_plotting, aes(y = Bayes_factor_01, x = BF_type))+
  geom_point()+
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=c(Inference_threshold_BF_high), ymax=c(Inf)), alpha = 0.2, fill=c("pink"))+
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=c(Inference_threshold_BF_low), ymax=c(Inference_threshold_BF_high)), alpha = 0.2, fill=c("grey80"))+
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=c(0), ymax=c(Inference_threshold_BF_low)), alpha = 0.2, fill=c("lightgreen"))+
  geom_point(size = 3.5, shape = 21, fill = "white")+
  scale_y_log10(limits = c(0.005,200), breaks=c(0.01, Inference_threshold_BF_low, 0.1, 0.33, 0, 3, 10, Inference_threshold_BF_high, 100))+
  geom_hline(yintercept = c(Inference_threshold_BF_low, Inference_threshold_BF_high), linetype = "dashed")+
  geom_text(aes(x=0.5, y=c(100, 1, 0.01), label=c("Supports M0", "Inconclusive", "Supports M1"), angle = 270))




##################################################
#       Report of confirmatory analysis          #
##################################################

# number of erotic trials performed (without the last participant, whose data might
# be cut in half because of the interim stopping point)
numtrials_without_last_part = nrow(data_BF[data_BF[,"participant_ID"] != data_BF[nrow(data_BF),"participant_ID"],])
# number of missing trials/data points (minus the last participant, who might be cut in half)
missing_data_points = (length(unique(data_BF[data_BF[,"participant_ID"] != data_BF[nrow(data_BF),"participant_ID"], "participant_ID"])) * trial_size_per_participant) - numtrials_without_last_part 


general_text = paste("The study stopped at ", total_N, " erotic trials gathered from a total of ",
                     sample_size_participants, " participants.", " There has been ", missing_data_points, " (", round(missing_data_points/total_N*100,2),"%) missing data points due to unfinished sessions.", 
                     " We observed a total of ", round(mean(sides_match), 4)*100, "% successful guesses within ",
                     total_N, " erotic trials (99.5% CI = ", proportion_995CI[1]*100, "%", ", ", proportion_995CI[2]*100, "%", 
                     "; posterior mode = ", hdi_result[1]*100, "%", ", posterior 90% HDI = ", HDI_lb*100, "%", ", ",
                     HDI_ub*100, "%", ").", sep = "")

robustness_text = if(Robust == "robust"){paste(" The results proved to be robust to different statistical approaches, increasing our confidence in our inference.")} else if(Robust == "not robust"){
  paste(" However, the results did not prove to be robust to different statistical approaches.")}





######################################################################
#                                                                    #
#                         Exploratory analysis                       #
#                                                                    #
######################################################################
# The existence of individual differences between participants will be evaluated 
# using a frequentist approach where we assume that it is possible that only a small 
# subgroup of the population we sample from is capable of predicting future events
# with better than chance accuracy.

# EXPLORATORY ANALYSIS RESULTS WILL NOT AFFECT THE CONCLUSIONS OF OUR STUDY

#=======================================================================#
#           Comparison of expected and observed distributions           #
#=======================================================================#

# calculate proportion of successful guesses for each participant in the observed data
data_BF_split = split(data_BF, f = data_BF[,"participant_ID"])
success_proportions_empirical = sapply(data_BF_split, function(x) mean(x[,"guessed_side"] == x[,"target_side"]))


# samples 1,000,000 participants from a population with H0 success rate
# this is used for the stochastic dominance test as the null model
# we call this the theoretical sample, because it approximates the theoretical null model
sim_null_participant_num = 1000000
success_proportions_theoretical <- rbinom(sim_null_participant_num, size = trial_size_per_participant, prob=M0_prob)/trial_size_per_participant

# determine possible values of success rates
possible_success_rates = 0
for(i in 1:trial_size_per_participant){
  possible_success_rates[i+1] = round(1/(trial_size_per_participant/i), 2)
}
possible_success_rates_char = as.character(possible_success_rates)
success_proportions_theoretical_char_rounded = as.character(round(success_proportions_theoretical, 2))


# Round the success rates of each participant to the nearest "possible success rate" value. 
# This is necessary to allow for missing trials for some participants, which may generate
# success rates that are different from the "possible success rate" for those completing all trials. 
success_proportions_empirical_rounded <- sapply(success_proportions_empirical, function(x) possible_success_rates[which.min(abs(possible_success_rates - x))])

success_proportions_empirical_char_rounded = as.character(round(success_proportions_empirical_rounded, 2))

success_rates_theoretical = NA
for(i in 1:length(possible_success_rates)){
  success_rates_theoretical[i] = sum(success_proportions_theoretical_char_rounded == possible_success_rates_char[i])
}
success_rates_theoretical_prop = matrix(success_rates_theoretical/sum(success_rates_theoretical))


success_rates_empirical = NA
for(i in 1:length(possible_success_rates)){
  success_rates_empirical[i] = sum(success_proportions_empirical_char_rounded == possible_success_rates_char[i])
}
success_rates_empirical_prop = matrix(success_rates_empirical/sum(success_rates_empirical))



# plot 
histogram_plot_data = as.data.frame(c(success_rates_theoretical_prop, success_rates_empirical_prop))
histogram_plot_data = cbind(histogram_plot_data, factor(c(rep("Expected if M0 is true", length(success_rates_theoretical_prop)), rep("Observed", length(success_rates_empirical_prop)))))
histogram_plot_data = cbind(histogram_plot_data, factor(rep(possible_success_rates_char, 2)))
names(histogram_plot_data) = c("proportion", "group", "success")

figure_1 =  ggplot(histogram_plot_data, aes(y = proportion, x = success, group = group))+
  geom_bar(aes(fill = group), alpha = 0.5, stat = "identity", position = "identity")+
  scale_fill_manual(values = c("darkgrey", "black")) +
  xlab("Successful guess rate") +
  ylab("Proportion") +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.line = element_line(colour = "black", size = 1.2),
        axis.text.x = element_text(angle = 90, color = "black", face = "bold", size = 12, margin = margin(1,0,0,0,"mm"), vjust = 0.5),
        axis.text.y = element_text(face = "bold", color = "black", size = 12),
        axis.title = element_text(size = 16))


# earth mover's distance
emd2d(success_rates_theoretical_prop,success_rates_empirical_prop)




















shinyServer(function(input, output){
  
#  sourceData <- reactive({
#    invalidateLater(1000000,session)
#    
#    functionThatGetsData()
#  })
  
  
  
  
  
  
  
  
  
  
  output$text_summary <- renderText({
    
    
    ######################################################################
    #                                                                    #
    #                         Data management                            #
    #                                                                    #
    ######################################################################
    
    # get up to date date from github
    pilot_data_pre <- read.csv("https://raw.githubusercontent.com/gy0p4k/transparent-psi-results/master/results.csv")
    
    # only read pilot sessions
    # pilot sessions are marked with session_type = "pilot" and they were all collected in the lab with the laboratory_ID_code = "lab_ELTE_01"
    lab_ID <- "lab_ELTE_01"
    pilot_data <- pilot_data_pre[pilot_data_pre[,"session_type"] == "pilot" & pilot_data_pre[,"laboratory_ID_code"] == "lab_ELTE_01", ]
    
    
    # Number of participants tested in the pilot test
    sample_size_participants = length(unique(pilot_data[, "participant_ID"]))
    
    
    pilot_data[, "trial_number"] = as.numeric(pilot_data[, "trial_number"])
    
    #extract data from erotic trials 
    data_BF = pilot_data[!is.na(pilot_data[, "trial_number"]) & pilot_data[, "reward_type"] == "erotic", ]
    
    data_BF[,"participant_ID"] = droplevels(data_BF[,"participant_ID"])
    
    
    ######################################################################
    #                                                                    #
    #                          Data analysis                             #
    #                                                                    #
    ######################################################################
    
    
    ######################################################################
    #                       Primary confirmatory test                    #
    ######################################################################
    
    #================================================================#
    #                Calculate number of successes                   #
    #================================================================#
    
    # number of successes and total N of trials
    sides_match = data_BF[,"guessed_side"] == data_BF[,"target_side"]
    successes = sum(sides_match)
    total_N = nrow(data_BF)
    
    
    #================================================================#
    #        Calculating Bayes factors using different priors        #
    #================================================================#
    
    ### Replication Bayes factor, with the Bem 2011 experiment 1 results providing the prior information
    
    BF_replication <- BF01_beta(y = successes, N = total_N, y_prior = y_prior, N_prior = N_prior, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null
    
    
    ### Bayes factor with uniform prior
    # using a non-informative flat prior distribution with alpha = 1 and beta = 1
    
    BF_uniform <- BF01_beta(y = successes, N = total_N, y_prior = 0, N_prior = 0, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null
    
    
    ### Bayes factor with BUJ prior
    # the BUJ prior is calculated from Bem's paper where the prior distribution is defined as a
    # normal distribution with a mean at 0 and 90th percentele is at medium effect size d = 0.5 
    # (we asume that this is one-tailed). Source: Bem, D. J., Utts, J., & Johnson, W. O. (2011). 
    # Must psychologists change the way they analyze their data? Journal of Personality and Social Psychology, 101(4), 716-719.
    # We simulate this in this binomial framework with a one-tailed beta distribution with alpha = 7 and beta = 7.
    # This distribution has 90% of its probability mass under p = 0.712, which we determined 
    # to be equivalent to d = 0.5 medium effect size. We used the formula to convert d to log odds ratio logodds = d*pi/sqrt(3), 
    # found here: Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R. (2009). 
    # Converting Among Effect Sizes. In Introduction to Meta-Analysis (pp. 45-49): John Wiley & Sons, Ltd.
    # Then, log odds ratio vas converted to probability using the formula: p = exp(x)/(1+exp(x))
    # The final equation: exp(d*pi/sqrt(3))/(1+exp(d*pi/sqrt(3)))
    
    BF_BUJ <- BF01_beta(y = successes, N = total_N, y_prior = 6, N_prior = 12, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null
    
    
    #================================================================#
    #                    Main analysis inference                     #
    #================================================================#
    
    # determine inference (supported model) based on the Bayes factors calculated above  
    if(Inference_threshold_BF_low >= max(c(BF_replication, BF_uniform, BF_BUJ))) {
      inference_BF = "M1"
      break} else if(Inference_threshold_BF_high <= min(c(BF_replication, BF_uniform, BF_BUJ))) {
        inference_BF = "M0"
        break} else {inference_BF = "Inconclusive"}
    
    
    ##################################################
    #                 Robustness tests               #
    ##################################################
    
    #================================================================#
    #               Robustness test of BF results with NHST          #
    #================================================================#
    
    # robustness of BF results is tested with NHST proportion tests
    # here we perform both an equality test and an equivalence test to draw statistical inference
    
    
    # equality test
    equality_test_p = prop.test(x = successes, n = total_N, p = M0_prob, alternative = "greater")$p.value
    
    # equivalence test
    equivalence_test_p = prop.test(x = successes, n = total_N,
                                   p = M0_prob+minimum_effect_threshold_NHST, 
                                   alternative = "less")$p.value
    
    # descritives of the frequentist estimate
    pbar = successes/total_N
    SE = sqrt(pbar * (1-pbar)/total_N)
    E = qnorm(.9975) * SE
    proportion_995CI = round(pbar + c(-E , E), 3)
    
    
    
    # making inference decision
    if((Inference_threshold_robustness_NHST > equality_test_p) & 
       (Inference_threshold_robustness_NHST <= equivalence_test_p)){
      inference_robustness_NHST = "M1"} else if((Inference_threshold_robustness_NHST > equivalence_test_p) & 
                                                (Inference_threshold_robustness_NHST <= equality_test_p)){
        inference_robustness_NHST = "M0" 
      } else {inference_robustness_NHST = "Inconclusive"}
    
    
    #=======================================================================#
    #   Robustness test of BF results with Bayesian parameter estimation    #
    #=======================================================================#
    
    # robustness of BF results is tested by calculating HDI of the posteriro distribution and checking its relation to
    # the region of practical equivalence (ROPE), promoted in Kruschke, J. K., & Liddell, T. M. 
    # (2017). The Bayesian New Statistics: Hypothesis testing, estimation, meta-analysis, and power 
    # analysis from a Bayesian perspective. Psychonomic Bulletin & Review, 1-29. 
    
    # calculate posterior distribution using beta distribution
    prior_alpha = y_prior + 1
    prior_beta = N_prior-y_prior+1
    
    posterior_alpha = prior_alpha + successes
    posterior_beta = prior_beta + total_N - successes
    
    scale = seq(0, 1, length = 1001)
    posterior_density = dbeta(scale, posterior_alpha, posterior_beta)
    
    # calculate HDI for the posterior distribution
    # (here we calculate the upper and lower bound of the 90% of the probability mass
    # because we use a one-tailed test. This means that the total probability mass below
    # the upper bound of the 90% HDI will be 95%)
    hdi_result = mode_HDI(scale = scale, density = posterior_density, crit_width = 1-Inference_threshold_robustness_Bayes_Par_Est*2, n_samples = 1e6)
    
    # parameters for decision making
    HDI_lb = hdi_result[2]
    HDI_ub = hdi_result[3]
    
    # making inference decision
    ROPE = M0_prob+minimum_effect_threshold_Bayes_Par_Est
    if(HDI_lb >= ROPE){inference_robustness_Bayes_Par_Est = "M1"
    } else if(HDI_ub <= ROPE){inference_robustness_Bayes_Par_Est = "M0"
    } else {inference_robustness_Bayes_Par_Est = "Inconclusive"}
    
    
    #=======================================================================#
    #      Determine final inference of all robustness tests combined       #
    #=======================================================================#
    
    # the main analysis inference is only robust if all robustness tests came to the same inference as the main test
    inferences = c(inference_robustness_NHST, inference_robustness_Bayes_Par_Est)
    inference_robustness = if(all(inferences == inferences[1])){inferences[1]} else {"mixed"}
    
    Robust = if(inference_BF == inference_robustness){"robust"} else {"not robust"}
    
    
    ######################################################################
    #                                                                    #
    #                   Interpretation and visualization                 #
    #                                                                    #
    ######################################################################
    
    
    # visualize results
    # ALL THREE BF values need to cross the threshold for the interpretation to be accepted
    
    BF_results_for_plotting = cbind(as.data.frame(c(BF_replication, BF_uniform, BF_BUJ)), c("BF_replication", "BF_uniform", "BF_BUJ"))
    names(BF_results_for_plotting) = c("Bayes_factor_01", "BF_type")
    
    plot <- ggplot(BF_results_for_plotting, aes(y = Bayes_factor_01, x = BF_type))+
      geom_point()+
      geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=c(Inference_threshold_BF_high), ymax=c(Inf)), alpha = 0.2, fill=c("pink"))+
      geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=c(Inference_threshold_BF_low), ymax=c(Inference_threshold_BF_high)), alpha = 0.2, fill=c("grey80"))+
      geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=c(0), ymax=c(Inference_threshold_BF_low)), alpha = 0.2, fill=c("lightgreen"))+
      geom_point(size = 3.5, shape = 21, fill = "white")+
      scale_y_log10(limits = c(0.005,200), breaks=c(0.01, Inference_threshold_BF_low, 0.1, 0.33, 0, 3, 10, Inference_threshold_BF_high, 100))+
      geom_hline(yintercept = c(Inference_threshold_BF_low, Inference_threshold_BF_high), linetype = "dashed")+
      geom_text(aes(x=0.5, y=c(100, 1, 0.01), label=c("Supports M0", "Inconclusive", "Supports M1"), angle = 270))
    
    
    
    
    ##################################################
    #       Report of confirmatory analysis          #
    ##################################################
    
    # number of erotic trials performed (without the last participant, whose data might
    # be cut in half because of the interim stopping point)
    numtrials_without_last_part = nrow(data_BF[data_BF[,"participant_ID"] != data_BF[nrow(data_BF),"participant_ID"],])
    # number of missing trials/data points (minus the last participant, who might be cut in half)
    missing_data_points = (length(unique(data_BF[data_BF[,"participant_ID"] != data_BF[nrow(data_BF),"participant_ID"], "participant_ID"])) * trial_size_per_participant) - numtrials_without_last_part 
    
    
    general_text = paste("The study stopped at ", total_N, " erotic trials gathered from a total of ",
                         sample_size_participants, " participants.", " There has been ", missing_data_points, " (", round(missing_data_points/total_N*100,2),"%) missing data points due to unfinished sessions.", 
                         " We observed a total of ", round(mean(sides_match), 4)*100, "% successful guesses within ",
                         total_N, " erotic trials (99.5% CI = ", proportion_995CI[1]*100, "%", ", ", proportion_995CI[2]*100, "%", 
                         "; posterior mode = ", hdi_result[1]*100, "%", ", posterior 90% HDI = ", HDI_lb*100, "%", ", ",
                         HDI_ub*100, "%", ").", sep = "")
    
    robustness_text = if(Robust == "robust"){paste(" The results proved to be robust to different statistical approaches, increasing our confidence in our inference.")} else if(Robust == "not robust"){
      paste(" However, the results did not prove to be robust to different statistical approaches.")}
    
    
    
    
    
    
    
    
    if(inference_BF == "M1"){
      paste(general_text, " Observing this success rate is ", round(1/max(c(BF_replication, BF_uniform, BF_BUJ)),0),
            " times more likely if humans can guess future randomly determined events than if they are guessing randomly. Taken at face value, the data provide strong evidence that the probability of successfully guessing later computer-generated random events is higher than chance level as previously reported by Bem (2011) and others (Bem, Tressoldi, Rabeyron, & Duggan, 2015).",
            robustness_text, sep = "")
    } else if(inference_BF == "M0"){
      paste(general_text, " Observing this success rate is ", round(min(c(BF_replication, BF_uniform, BF_BUJ)),0),
            " times more likely if humans are guessing randomly than if they can guess future randomly determined events. Taken at face value, the data provide strong evidence that the probability of successfully guessing later computer-generated random events is not higher than chance level as previously reported by Bem (2011) and others (Bem, Tressoldi, Rabeyron, & Duggan, 2015).",
            robustness_text, sep = "")
    } else if(inference_BF == "Inconclusive" & max(c(BF_replication, BF_uniform, BF_BUJ)) < 1){
      paste(general_text, " Observing this success rate is ", round(1/max(c(BF_replication, BF_uniform, BF_BUJ)),0),
            " times more likely if humans can guess future randomly determined events than if they are guessing randomly. However, this study outcome did not reach the pre-specified criteria of strong support for either model.",
            robustness_text, sep = "")
    } else if(inference_BF == "Inconclusive" & min(c(BF_replication, BF_uniform, BF_BUJ)) > 1){
      paste(general_text, " Observing this success rate is ", round(min(c(BF_replication, BF_uniform, BF_BUJ)),0),
            " times more likely if humans are guessing randomly than if they can guess future randomly determined events. However, this study outcome did not reach the pre-specified criteria of strong support for either model.",
            robustness_text, sep = "")
    } else {paste(general_text, "However, this study outcome did not reach the pre-specified criteria of strong support for either model.",
                  robustness_text, sep = "")}
    
  })
  
  
  
  
  
  
  
  
  
  
  
  
  output$plot_refresh <- renderPlot({
    invalidateLater(2000)
    
    dratta <- as.data.frame(cbind(c(1:10), c(5 + as.numeric(substr(as.numeric(Sys.time()), 10, 12)), 7:15)))
    
    plot(dratta)
  })

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  output$plot1a <- renderPlot({
    invalidateLater(2000)
    
    
    
    ######################################################################
    #                                                                    #
    #                         Data management                            #
    #                                                                    #
    ######################################################################
    
    # get up to date date from github
    pilot_data_pre <- read.csv("https://raw.githubusercontent.com/gy0p4k/transparent-psi-results/master/results.csv")
    
    # only read pilot sessions
    # pilot sessions are marked with session_type = "pilot" and they were all collected in the lab with the laboratory_ID_code = "lab_ELTE_01"
    lab_ID <- "lab_ELTE_01"
    pilot_data <- pilot_data_pre[pilot_data_pre[,"session_type"] == "pilot" & pilot_data_pre[,"laboratory_ID_code"] == "lab_ELTE_01", ]
    
    
    # Number of participants tested in the pilot test
    sample_size_participants = length(unique(pilot_data[, "participant_ID"]))
    
    
    pilot_data[, "trial_number"] = as.numeric(pilot_data[, "trial_number"])
    
    #extract data from erotic trials 
    data_BF = pilot_data[!is.na(pilot_data[, "trial_number"]) & pilot_data[, "reward_type"] == "erotic", ]
    
    data_BF[,"participant_ID"] = droplevels(data_BF[,"participant_ID"])
    
    
    ######################################################################
    #                                                                    #
    #                          Data analysis                             #
    #                                                                    #
    ######################################################################
    
    
    ######################################################################
    #                       Primary confirmatory test                    #
    ######################################################################
    
    #================================================================#
    #                Calculate number of successes                   #
    #================================================================#
    
    # number of successes and total N of trials
    sides_match = data_BF[,"guessed_side"] == data_BF[,"target_side"]
    successes = sum(sides_match)
    total_N = nrow(data_BF)
    
    
    #================================================================#
    #        Calculating Bayes factors using different priors        #
    #================================================================#
    
    ### Replication Bayes factor, with the Bem 2011 experiment 1 results providing the prior information
    
    BF_replication <- BF01_beta(y = successes, N = total_N, y_prior = y_prior, N_prior = N_prior, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null
    
    
    ### Bayes factor with uniform prior
    # using a non-informative flat prior distribution with alpha = 1 and beta = 1
    
    BF_uniform <- BF01_beta(y = successes, N = total_N, y_prior = 0, N_prior = 0, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null
    
    
    ### Bayes factor with BUJ prior
    # the BUJ prior is calculated from Bem's paper where the prior distribution is defined as a
    # normal distribution with a mean at 0 and 90th percentele is at medium effect size d = 0.5 
    # (we asume that this is one-tailed). Source: Bem, D. J., Utts, J., & Johnson, W. O. (2011). 
    # Must psychologists change the way they analyze their data? Journal of Personality and Social Psychology, 101(4), 716-719.
    # We simulate this in this binomial framework with a one-tailed beta distribution with alpha = 7 and beta = 7.
    # This distribution has 90% of its probability mass under p = 0.712, which we determined 
    # to be equivalent to d = 0.5 medium effect size. We used the formula to convert d to log odds ratio logodds = d*pi/sqrt(3), 
    # found here: Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R. (2009). 
    # Converting Among Effect Sizes. In Introduction to Meta-Analysis (pp. 45-49): John Wiley & Sons, Ltd.
    # Then, log odds ratio vas converted to probability using the formula: p = exp(x)/(1+exp(x))
    # The final equation: exp(d*pi/sqrt(3))/(1+exp(d*pi/sqrt(3)))
    
    BF_BUJ <- BF01_beta(y = successes, N = total_N, y_prior = 6, N_prior = 12, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null
    
    
    #================================================================#
    #                    Main analysis inference                     #
    #================================================================#
    
    # determine inference (supported model) based on the Bayes factors calculated above  
    if(Inference_threshold_BF_low >= max(c(BF_replication, BF_uniform, BF_BUJ))) {
      inference_BF = "M1"
      break} else if(Inference_threshold_BF_high <= min(c(BF_replication, BF_uniform, BF_BUJ))) {
        inference_BF = "M0"
        break} else {inference_BF = "Inconclusive"}
    
    
    ##################################################
    #                 Robustness tests               #
    ##################################################
    
    #================================================================#
    #               Robustness test of BF results with NHST          #
    #================================================================#
    
    # robustness of BF results is tested with NHST proportion tests
    # here we perform both an equality test and an equivalence test to draw statistical inference
    
    
    # equality test
    equality_test_p = prop.test(x = successes, n = total_N, p = M0_prob, alternative = "greater")$p.value
    
    # equivalence test
    equivalence_test_p = prop.test(x = successes, n = total_N,
                                   p = M0_prob+minimum_effect_threshold_NHST, 
                                   alternative = "less")$p.value
    
    # descritives of the frequentist estimate
    pbar = successes/total_N
    SE = sqrt(pbar * (1-pbar)/total_N)
    E = qnorm(.9975) * SE
    proportion_995CI = round(pbar + c(-E , E), 3)
    
    
    
    # making inference decision
    if((Inference_threshold_robustness_NHST > equality_test_p) & 
       (Inference_threshold_robustness_NHST <= equivalence_test_p)){
      inference_robustness_NHST = "M1"} else if((Inference_threshold_robustness_NHST > equivalence_test_p) & 
                                                (Inference_threshold_robustness_NHST <= equality_test_p)){
        inference_robustness_NHST = "M0" 
      } else {inference_robustness_NHST = "Inconclusive"}
    
    
    #=======================================================================#
    #   Robustness test of BF results with Bayesian parameter estimation    #
    #=======================================================================#
    
    # robustness of BF results is tested by calculating HDI of the posteriro distribution and checking its relation to
    # the region of practical equivalence (ROPE), promoted in Kruschke, J. K., & Liddell, T. M. 
    # (2017). The Bayesian New Statistics: Hypothesis testing, estimation, meta-analysis, and power 
    # analysis from a Bayesian perspective. Psychonomic Bulletin & Review, 1-29. 
    
    # calculate posterior distribution using beta distribution
    prior_alpha = y_prior + 1
    prior_beta = N_prior-y_prior+1
    
    posterior_alpha = prior_alpha + successes
    posterior_beta = prior_beta + total_N - successes
    
    scale = seq(0, 1, length = 1001)
    posterior_density = dbeta(scale, posterior_alpha, posterior_beta)
    
    # calculate HDI for the posterior distribution
    # (here we calculate the upper and lower bound of the 90% of the probability mass
    # because we use a one-tailed test. This means that the total probability mass below
    # the upper bound of the 90% HDI will be 95%)
    hdi_result = mode_HDI(scale = scale, density = posterior_density, crit_width = 1-Inference_threshold_robustness_Bayes_Par_Est*2, n_samples = 1e6)
    
    # parameters for decision making
    HDI_lb = hdi_result[2]
    HDI_ub = hdi_result[3]
    
    # making inference decision
    ROPE = M0_prob+minimum_effect_threshold_Bayes_Par_Est
    if(HDI_lb >= ROPE){inference_robustness_Bayes_Par_Est = "M1"
    } else if(HDI_ub <= ROPE){inference_robustness_Bayes_Par_Est = "M0"
    } else {inference_robustness_Bayes_Par_Est = "Inconclusive"}
    
    
    #=======================================================================#
    #      Determine final inference of all robustness tests combined       #
    #=======================================================================#
    
    # the main analysis inference is only robust if all robustness tests came to the same inference as the main test
    inferences = c(inference_robustness_NHST, inference_robustness_Bayes_Par_Est)
    inference_robustness = if(all(inferences == inferences[1])){inferences[1]} else {"mixed"}
    
    Robust = if(inference_BF == inference_robustness){"robust"} else {"not robust"}
    
    
    ######################################################################
    #                                                                    #
    #                   Interpretation and visualization                 #
    #                                                                    #
    ######################################################################
    
    
    # visualize results
    # ALL THREE BF values need to cross the threshold for the interpretation to be accepted
    
    BF_results_for_plotting = cbind(as.data.frame(c(BF_replication, BF_uniform, BF_BUJ)), c("BF_replication", "BF_uniform", "BF_BUJ"))
    names(BF_results_for_plotting) = c("Bayes_factor_01", "BF_type")
    
    plot <- ggplot(BF_results_for_plotting, aes(y = Bayes_factor_01, x = BF_type))+
      geom_point()+
      geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=c(Inference_threshold_BF_high), ymax=c(Inf)), alpha = 0.2, fill=c("pink"))+
      geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=c(Inference_threshold_BF_low), ymax=c(Inference_threshold_BF_high)), alpha = 0.2, fill=c("grey80"))+
      geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=c(0), ymax=c(Inference_threshold_BF_low)), alpha = 0.2, fill=c("lightgreen"))+
      geom_point(size = 3.5, shape = 21, fill = "white")+
      scale_y_log10(limits = c(0.005,200), breaks=c(0.01, Inference_threshold_BF_low, 0.1, 0.33, 0, 3, 10, Inference_threshold_BF_high, 100))+
      geom_hline(yintercept = c(Inference_threshold_BF_low, Inference_threshold_BF_high), linetype = "dashed")+
      geom_text(aes(x=0.5, y=c(100, 1, 0.01), label=c("Supports M0", "Inconclusive", "Supports M1"), angle = 270))
    
    
    
    
    ##################################################
    #       Report of confirmatory analysis          #
    ##################################################
    
    # number of erotic trials performed (without the last participant, whose data might
    # be cut in half because of the interim stopping point)
    numtrials_without_last_part = nrow(data_BF[data_BF[,"participant_ID"] != data_BF[nrow(data_BF),"participant_ID"],])
    # number of missing trials/data points (minus the last participant, who might be cut in half)
    missing_data_points = (length(unique(data_BF[data_BF[,"participant_ID"] != data_BF[nrow(data_BF),"participant_ID"], "participant_ID"])) * trial_size_per_participant) - numtrials_without_last_part 
    
    
    general_text = paste("The study stopped at ", total_N, " erotic trials gathered from a total of ",
                         sample_size_participants, " participants.", " There has been ", missing_data_points, " (", round(missing_data_points/total_N*100,2),"%) missing data points due to unfinished sessions.", 
                         " We observed a total of ", round(mean(sides_match), 4)*100, "% successful guesses within ",
                         total_N, " erotic trials (99.5% CI = ", proportion_995CI[1]*100, "%", ", ", proportion_995CI[2]*100, "%", 
                         "; posterior mode = ", hdi_result[1]*100, "%", ", posterior 90% HDI = ", HDI_lb*100, "%", ", ",
                         HDI_ub*100, "%", ").", sep = "")
    
    robustness_text = if(Robust == "robust"){paste(" The results proved to be robust to different statistical approaches, increasing our confidence in our inference.")} else if(Robust == "not robust"){
    paste(" However, the results did not prove to be robust to different statistical approaches.")}
    
    print(plot)
  })
  
  
  
  
  
  
  
  
  
  
  
  
  output$plot1b <- renderPlot({
    invalidateLater(2000)
    
    
    
    ######################################################################
    #                                                                    #
    #                         Data management                            #
    #                                                                    #
    ######################################################################
    
    # get up to date date from github
    pilot_data_pre <- read.csv("https://raw.githubusercontent.com/gy0p4k/transparent-psi-results/master/results.csv")
    
    # only read pilot sessions
    # pilot sessions are marked with session_type = "pilot" and they were all collected in the lab with the laboratory_ID_code = "lab_ELTE_01"
    lab_ID <- "lab_ELTE_01"
    pilot_data <- pilot_data_pre[pilot_data_pre[,"session_type"] == "pilot" & pilot_data_pre[,"laboratory_ID_code"] == "lab_ELTE_01", ]
    
    
    # Number of participants tested in the pilot test
    sample_size_participants = length(unique(pilot_data[, "participant_ID"]))
    
    
    pilot_data[, "trial_number"] = as.numeric(pilot_data[, "trial_number"])
    
    #extract data from erotic trials 
    data_BF = pilot_data[!is.na(pilot_data[, "trial_number"]) & pilot_data[, "reward_type"] == "erotic", ]
    
    data_BF[,"participant_ID"] = droplevels(data_BF[,"participant_ID"])
    
    
    ######################################################################
    #                                                                    #
    #                          Data analysis                             #
    #                                                                    #
    ######################################################################
    
    
    ######################################################################
    #                       Primary confirmatory test                    #
    ######################################################################
    
    #================================================================#
    #                Calculate number of successes                   #
    #================================================================#
    
    # number of successes and total N of trials
    sides_match = data_BF[,"guessed_side"] == data_BF[,"target_side"]
    successes = sum(sides_match)
    total_N = nrow(data_BF)
    
    
    #================================================================#
    #        Calculating Bayes factors using different priors        #
    #================================================================#
    
    ### Replication Bayes factor, with the Bem 2011 experiment 1 results providing the prior information
    
    BF_replication <- BF01_beta(y = successes, N = total_N, y_prior = y_prior, N_prior = N_prior, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null
    
    
    ### Bayes factor with uniform prior
    # using a non-informative flat prior distribution with alpha = 1 and beta = 1
    
    BF_uniform <- BF01_beta(y = successes, N = total_N, y_prior = 0, N_prior = 0, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null
    
    
    ### Bayes factor with BUJ prior
    # the BUJ prior is calculated from Bem's paper where the prior distribution is defined as a
    # normal distribution with a mean at 0 and 90th percentele is at medium effect size d = 0.5 
    # (we asume that this is one-tailed). Source: Bem, D. J., Utts, J., & Johnson, W. O. (2011). 
    # Must psychologists change the way they analyze their data? Journal of Personality and Social Psychology, 101(4), 716-719.
    # We simulate this in this binomial framework with a one-tailed beta distribution with alpha = 7 and beta = 7.
    # This distribution has 90% of its probability mass under p = 0.712, which we determined 
    # to be equivalent to d = 0.5 medium effect size. We used the formula to convert d to log odds ratio logodds = d*pi/sqrt(3), 
    # found here: Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R. (2009). 
    # Converting Among Effect Sizes. In Introduction to Meta-Analysis (pp. 45-49): John Wiley & Sons, Ltd.
    # Then, log odds ratio vas converted to probability using the formula: p = exp(x)/(1+exp(x))
    # The final equation: exp(d*pi/sqrt(3))/(1+exp(d*pi/sqrt(3)))
    
    BF_BUJ <- BF01_beta(y = successes, N = total_N, y_prior = 6, N_prior = 12, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null
    
    
    #================================================================#
    #                    Main analysis inference                     #
    #================================================================#
    
    # determine inference (supported model) based on the Bayes factors calculated above  
    if(Inference_threshold_BF_low >= max(c(BF_replication, BF_uniform, BF_BUJ))) {
      inference_BF = "M1"
      break} else if(Inference_threshold_BF_high <= min(c(BF_replication, BF_uniform, BF_BUJ))) {
        inference_BF = "M0"
        break} else {inference_BF = "Inconclusive"}
    
    
    ##################################################
    #                 Robustness tests               #
    ##################################################
    
    #================================================================#
    #               Robustness test of BF results with NHST          #
    #================================================================#
    
    # robustness of BF results is tested with NHST proportion tests
    # here we perform both an equality test and an equivalence test to draw statistical inference
    
    
    # equality test
    equality_test_p = prop.test(x = successes, n = total_N, p = M0_prob, alternative = "greater")$p.value
    
    # equivalence test
    equivalence_test_p = prop.test(x = successes, n = total_N,
                                   p = M0_prob+minimum_effect_threshold_NHST, 
                                   alternative = "less")$p.value
    
    # descritives of the frequentist estimate
    pbar = successes/total_N
    SE = sqrt(pbar * (1-pbar)/total_N)
    E = qnorm(.9975) * SE
    proportion_995CI = round(pbar + c(-E , E), 3)
    
    
    
    # making inference decision
    if((Inference_threshold_robustness_NHST > equality_test_p) & 
       (Inference_threshold_robustness_NHST <= equivalence_test_p)){
      inference_robustness_NHST = "M1"} else if((Inference_threshold_robustness_NHST > equivalence_test_p) & 
                                                (Inference_threshold_robustness_NHST <= equality_test_p)){
        inference_robustness_NHST = "M0" 
      } else {inference_robustness_NHST = "Inconclusive"}
    
    
    #=======================================================================#
    #   Robustness test of BF results with Bayesian parameter estimation    #
    #=======================================================================#
    
    # robustness of BF results is tested by calculating HDI of the posteriro distribution and checking its relation to
    # the region of practical equivalence (ROPE), promoted in Kruschke, J. K., & Liddell, T. M. 
    # (2017). The Bayesian New Statistics: Hypothesis testing, estimation, meta-analysis, and power 
    # analysis from a Bayesian perspective. Psychonomic Bulletin & Review, 1-29. 
    
    # calculate posterior distribution using beta distribution
    prior_alpha = y_prior + 1
    prior_beta = N_prior-y_prior+1
    
    posterior_alpha = prior_alpha + successes
    posterior_beta = prior_beta + total_N - successes
    
    scale = seq(0, 1, length = 1001)
    posterior_density = dbeta(scale, posterior_alpha, posterior_beta)
    
    # calculate HDI for the posterior distribution
    # (here we calculate the upper and lower bound of the 90% of the probability mass
    # because we use a one-tailed test. This means that the total probability mass below
    # the upper bound of the 90% HDI will be 95%)
    hdi_result = mode_HDI(scale = scale, density = posterior_density, crit_width = 1-Inference_threshold_robustness_Bayes_Par_Est*2, n_samples = 1e6)
    
    # parameters for decision making
    HDI_lb = hdi_result[2]
    HDI_ub = hdi_result[3]
    
    # making inference decision
    ROPE = M0_prob+minimum_effect_threshold_Bayes_Par_Est
    if(HDI_lb >= ROPE){inference_robustness_Bayes_Par_Est = "M1"
    } else if(HDI_ub <= ROPE){inference_robustness_Bayes_Par_Est = "M0"
    } else {inference_robustness_Bayes_Par_Est = "Inconclusive"}
    
    
    #=======================================================================#
    #      Determine final inference of all robustness tests combined       #
    #=======================================================================#
    
    # the main analysis inference is only robust if all robustness tests came to the same inference as the main test
    inferences = c(inference_robustness_NHST, inference_robustness_Bayes_Par_Est)
    inference_robustness = if(all(inferences == inferences[1])){inferences[1]} else {"mixed"}
    
    Robust = if(inference_BF == inference_robustness){"robust"} else {"not robust"}
    
    
    ######################################################################
    #                                                                    #
    #                   Interpretation and visualization                 #
    #                                                                    #
    ######################################################################
    
    
    # visualize results
    # ALL THREE BF values need to cross the threshold for the interpretation to be accepted
    
    BF_results_for_plotting = cbind(as.data.frame(c(BF_replication, BF_uniform, BF_BUJ)), c("BF_replication", "BF_uniform", "BF_BUJ"))
    names(BF_results_for_plotting) = c("Bayes_factor_01", "BF_type")
    
    plot <- ggplot(BF_results_for_plotting, aes(y = Bayes_factor_01, x = BF_type))+
      geom_point()+
      geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=c(Inference_threshold_BF_high), ymax=c(Inf)), alpha = 0.2, fill=c("pink"))+
      geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=c(Inference_threshold_BF_low), ymax=c(Inference_threshold_BF_high)), alpha = 0.2, fill=c("grey80"))+
      geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=c(0), ymax=c(Inference_threshold_BF_low)), alpha = 0.2, fill=c("lightgreen"))+
      geom_point(size = 3.5, shape = 21, fill = "white")+
      scale_y_log10(limits = c(0.005,200), breaks=c(0.01, Inference_threshold_BF_low, 0.1, 0.33, 0, 3, 10, Inference_threshold_BF_high, 100))+
      geom_hline(yintercept = c(Inference_threshold_BF_low, Inference_threshold_BF_high), linetype = "dashed")+
      geom_text(aes(x=0.5, y=c(100, 1, 0.01), label=c("Supports M0", "Inconclusive", "Supports M1"), angle = 270))
    
    
    
    
    ##################################################
    #       Report of confirmatory analysis          #
    ##################################################
    
    # number of erotic trials performed (without the last participant, whose data might
    # be cut in half because of the interim stopping point)
    numtrials_without_last_part = nrow(data_BF[data_BF[,"participant_ID"] != data_BF[nrow(data_BF),"participant_ID"],])
    # number of missing trials/data points (minus the last participant, who might be cut in half)
    missing_data_points = (length(unique(data_BF[data_BF[,"participant_ID"] != data_BF[nrow(data_BF),"participant_ID"], "participant_ID"])) * trial_size_per_participant) - numtrials_without_last_part 
    
    
    general_text = paste("The study stopped at ", total_N, " erotic trials gathered from a total of ",
                         sample_size_participants, " participants.", " There has been ", missing_data_points, " (", round(missing_data_points/total_N*100,2),"%) missing data points due to unfinished sessions.", 
                         " We observed a total of ", round(mean(sides_match), 4)*100, "% successful guesses within ",
                         total_N, " erotic trials (99.5% CI = ", proportion_995CI[1]*100, "%", ", ", proportion_995CI[2]*100, "%", 
                         "; posterior mode = ", hdi_result[1]*100, "%", ", posterior 90% HDI = ", HDI_lb*100, "%", ", ",
                         HDI_ub*100, "%", ").", sep = "")
    
    robustness_text = if(Robust == "robust"){paste(" The results proved to be robust to different statistical approaches, increasing our confidence in our inference.")} else if(Robust == "not robust"){
    paste(" However, the results did not prove to be robust to different statistical approaches.")}

    print(plot)
  })
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  output$plot2 <- renderPlot({
    invalidateLater(2000)
    
    ######################################################################
    #                                                                    #
    #                         Data management                            #
    #                                                                    #
    ######################################################################
    
    # get up to date date from github
    pilot_data_pre <- read.csv("https://raw.githubusercontent.com/gy0p4k/transparent-psi-results/master/results.csv")
    
    # only read pilot sessions
    # pilot sessions are marked with session_type = "pilot" and they were all collected in the lab with the laboratory_ID_code = "lab_ELTE_01"
    lab_ID <- "lab_ELTE_01"
    pilot_data <- pilot_data_pre[pilot_data_pre[,"session_type"] == "pilot" & pilot_data_pre[,"laboratory_ID_code"] == "lab_ELTE_01", ]
    
    
    # Number of participants tested in the pilot test
    sample_size_participants = length(unique(pilot_data[, "participant_ID"]))
    
    
    pilot_data[, "trial_number"] = as.numeric(pilot_data[, "trial_number"])
    
    #extract data from erotic trials 
    data_BF = pilot_data[!is.na(pilot_data[, "trial_number"]) & pilot_data[, "reward_type"] == "erotic", ]
    
    data_BF[,"participant_ID"] = droplevels(data_BF[,"participant_ID"])
    
    
    ######################################################################
    #                                                                    #
    #                          Data analysis                             #
    #                                                                    #
    ######################################################################
    
    
    ######################################################################
    #                       Primary confirmatory test                    #
    ######################################################################
    
    #================================================================#
    #                Calculate number of successes                   #
    #================================================================#
    
    # number of successes and total N of trials
    sides_match = data_BF[,"guessed_side"] == data_BF[,"target_side"]
    successes = sum(sides_match)
    total_N = nrow(data_BF)
    
    
    #================================================================#
    #        Calculating Bayes factors using different priors        #
    #================================================================#
    
    ### Replication Bayes factor, with the Bem 2011 experiment 1 results providing the prior information
    
    BF_replication <- BF01_beta(y = successes, N = total_N, y_prior = y_prior, N_prior = N_prior, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null
    
    
    ### Bayes factor with uniform prior
    # using a non-informative flat prior distribution with alpha = 1 and beta = 1
    
    BF_uniform <- BF01_beta(y = successes, N = total_N, y_prior = 0, N_prior = 0, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null
    
    
    ### Bayes factor with BUJ prior
    # the BUJ prior is calculated from Bem's paper where the prior distribution is defined as a
    # normal distribution with a mean at 0 and 90th percentele is at medium effect size d = 0.5 
    # (we asume that this is one-tailed). Source: Bem, D. J., Utts, J., & Johnson, W. O. (2011). 
    # Must psychologists change the way they analyze their data? Journal of Personality and Social Psychology, 101(4), 716-719.
    # We simulate this in this binomial framework with a one-tailed beta distribution with alpha = 7 and beta = 7.
    # This distribution has 90% of its probability mass under p = 0.712, which we determined 
    # to be equivalent to d = 0.5 medium effect size. We used the formula to convert d to log odds ratio logodds = d*pi/sqrt(3), 
    # found here: Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R. (2009). 
    # Converting Among Effect Sizes. In Introduction to Meta-Analysis (pp. 45-49): John Wiley & Sons, Ltd.
    # Then, log odds ratio vas converted to probability using the formula: p = exp(x)/(1+exp(x))
    # The final equation: exp(d*pi/sqrt(3))/(1+exp(d*pi/sqrt(3)))
    
    BF_BUJ <- BF01_beta(y = successes, N = total_N, y_prior = 6, N_prior = 12, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null
    
    
    #================================================================#
    #                    Main analysis inference                     #
    #================================================================#
    
    # determine inference (supported model) based on the Bayes factors calculated above  
    if(Inference_threshold_BF_low >= max(c(BF_replication, BF_uniform, BF_BUJ))) {
      inference_BF = "M1"
      break} else if(Inference_threshold_BF_high <= min(c(BF_replication, BF_uniform, BF_BUJ))) {
        inference_BF = "M0"
        break} else {inference_BF = "Inconclusive"}
    
    
    ##################################################
    #                 Robustness tests               #
    ##################################################
    
    #================================================================#
    #               Robustness test of BF results with NHST          #
    #================================================================#
    
    # robustness of BF results is tested with NHST proportion tests
    # here we perform both an equality test and an equivalence test to draw statistical inference
    
    
    # equality test
    equality_test_p = prop.test(x = successes, n = total_N, p = M0_prob, alternative = "greater")$p.value
    
    # equivalence test
    equivalence_test_p = prop.test(x = successes, n = total_N,
                                   p = M0_prob+minimum_effect_threshold_NHST, 
                                   alternative = "less")$p.value
    
    # descritives of the frequentist estimate
    pbar = successes/total_N
    SE = sqrt(pbar * (1-pbar)/total_N)
    E = qnorm(.9975) * SE
    proportion_995CI = round(pbar + c(-E , E), 3)
    
    
    
    # making inference decision
    if((Inference_threshold_robustness_NHST > equality_test_p) & 
       (Inference_threshold_robustness_NHST <= equivalence_test_p)){
      inference_robustness_NHST = "M1"} else if((Inference_threshold_robustness_NHST > equivalence_test_p) & 
                                                (Inference_threshold_robustness_NHST <= equality_test_p)){
        inference_robustness_NHST = "M0" 
      } else {inference_robustness_NHST = "Inconclusive"}
    
    
    #=======================================================================#
    #   Robustness test of BF results with Bayesian parameter estimation    #
    #=======================================================================#
    
    # robustness of BF results is tested by calculating HDI of the posteriro distribution and checking its relation to
    # the region of practical equivalence (ROPE), promoted in Kruschke, J. K., & Liddell, T. M. 
    # (2017). The Bayesian New Statistics: Hypothesis testing, estimation, meta-analysis, and power 
    # analysis from a Bayesian perspective. Psychonomic Bulletin & Review, 1-29. 
    
    # calculate posterior distribution using beta distribution
    prior_alpha = y_prior + 1
    prior_beta = N_prior-y_prior+1
    
    posterior_alpha = prior_alpha + successes
    posterior_beta = prior_beta + total_N - successes
    
    scale = seq(0, 1, length = 1001)
    posterior_density = dbeta(scale, posterior_alpha, posterior_beta)
    
    # calculate HDI for the posterior distribution
    # (here we calculate the upper and lower bound of the 90% of the probability mass
    # because we use a one-tailed test. This means that the total probability mass below
    # the upper bound of the 90% HDI will be 95%)
    hdi_result = mode_HDI(scale = scale, density = posterior_density, crit_width = 1-Inference_threshold_robustness_Bayes_Par_Est*2, n_samples = 1e6)
    
    # parameters for decision making
    HDI_lb = hdi_result[2]
    HDI_ub = hdi_result[3]
    
    # making inference decision
    ROPE = M0_prob+minimum_effect_threshold_Bayes_Par_Est
    if(HDI_lb >= ROPE){inference_robustness_Bayes_Par_Est = "M1"
    } else if(HDI_ub <= ROPE){inference_robustness_Bayes_Par_Est = "M0"
    } else {inference_robustness_Bayes_Par_Est = "Inconclusive"}
    
    
    #=======================================================================#
    #      Determine final inference of all robustness tests combined       #
    #=======================================================================#
    
    # the main analysis inference is only robust if all robustness tests came to the same inference as the main test
    inferences = c(inference_robustness_NHST, inference_robustness_Bayes_Par_Est)
    inference_robustness = if(all(inferences == inferences[1])){inferences[1]} else {"mixed"}
    
    Robust = if(inference_BF == inference_robustness){"robust"} else {"not robust"}
    
    
    ######################################################################
    #                                                                    #
    #                   Interpretation and visualization                 #
    #                                                                    #
    ######################################################################
    
    
    # visualize results
    # ALL THREE BF values need to cross the threshold for the interpretation to be accepted
    
    BF_results_for_plotting = cbind(as.data.frame(c(BF_replication, BF_uniform, BF_BUJ)), c("BF_replication", "BF_uniform", "BF_BUJ"))
    names(BF_results_for_plotting) = c("Bayes_factor_01", "BF_type")
    
    plot <- ggplot(BF_results_for_plotting, aes(y = Bayes_factor_01, x = BF_type))+
      geom_point()+
      geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=c(Inference_threshold_BF_high), ymax=c(Inf)), alpha = 0.2, fill=c("pink"))+
      geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=c(Inference_threshold_BF_low), ymax=c(Inference_threshold_BF_high)), alpha = 0.2, fill=c("grey80"))+
      geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=c(0), ymax=c(Inference_threshold_BF_low)), alpha = 0.2, fill=c("lightgreen"))+
      geom_point(size = 3.5, shape = 21, fill = "white")+
      scale_y_log10(limits = c(0.005,200), breaks=c(0.01, Inference_threshold_BF_low, 0.1, 0.33, 0, 3, 10, Inference_threshold_BF_high, 100))+
      geom_hline(yintercept = c(Inference_threshold_BF_low, Inference_threshold_BF_high), linetype = "dashed")+
      geom_text(aes(x=0.5, y=c(100, 1, 0.01), label=c("Supports M0", "Inconclusive", "Supports M1"), angle = 270))
    
    
    
    
    ##################################################
    #       Report of confirmatory analysis          #
    ##################################################
    
    # number of erotic trials performed (without the last participant, whose data might
    # be cut in half because of the interim stopping point)
    numtrials_without_last_part = nrow(data_BF[data_BF[,"participant_ID"] != data_BF[nrow(data_BF),"participant_ID"],])
    # number of missing trials/data points (minus the last participant, who might be cut in half)
    missing_data_points = (length(unique(data_BF[data_BF[,"participant_ID"] != data_BF[nrow(data_BF),"participant_ID"], "participant_ID"])) * trial_size_per_participant) - numtrials_without_last_part 
    
    
    general_text = paste("The study stopped at ", total_N, " erotic trials gathered from a total of ",
                         sample_size_participants, " participants.", " There has been ", missing_data_points, " (", round(missing_data_points/total_N*100,2),"%) missing data points due to unfinished sessions.", 
                         " We observed a total of ", round(mean(sides_match), 4)*100, "% successful guesses within ",
                         total_N, " erotic trials (99.5% CI = ", proportion_995CI[1]*100, "%", ", ", proportion_995CI[2]*100, "%", 
                         "; posterior mode = ", hdi_result[1]*100, "%", ", posterior 90% HDI = ", HDI_lb*100, "%", ", ",
                         HDI_ub*100, "%", ").", sep = "")
    
    robustness_text = if(Robust == "robust"){paste(" The results proved to be robust to different statistical approaches, increasing our confidence in our inference.")} else if(Robust == "not robust"){
      paste(" However, the results did not prove to be robust to different statistical approaches.")}
    
    
    
    # plot results of the Bayesian parameter estimation used as a robustness test
    plot(scale, posterior_density, type="l", lty=1, xlab="x value", xlim = c(0.45, 0.55),
         ylab="Density")
    abline(v=c(M0_prob, ROPE), lty = c(1, 2))
    
    density_table = as.data.frame(cbind(scale, posterior_density))
    
    height_lim_lb = density_table[density_table[, "scale"] == HDI_lb, "posterior_density"]
    height_lim_ub = density_table[density_table[, "scale"] == HDI_ub, "posterior_density"]
    
    clip(0,1,-10,height_lim_lb)
    abline(v=HDI_lb, lty = 3)
    clip(0,1,-10,height_lim_ub)
    abline(v=HDI_ub, lty = 3)
  })
  
  
  
  
  
  
  
  
  
  
  
  
  
  output$plot3 <- renderPlot({
    invalidateLater(2000)
    
    
    ######################################################################
    #                                                                    #
    #                         Data management                            #
    #                                                                    #
    ######################################################################
    
    # get up to date date from github
    pilot_data_pre <- read.csv("https://raw.githubusercontent.com/gy0p4k/transparent-psi-results/master/results.csv")
    
    # only read pilot sessions
    # pilot sessions are marked with session_type = "pilot" and they were all collected in the lab with the laboratory_ID_code = "lab_ELTE_01"
    lab_ID <- "lab_ELTE_01"
    pilot_data <- pilot_data_pre[pilot_data_pre[,"session_type"] == "pilot" & pilot_data_pre[,"laboratory_ID_code"] == "lab_ELTE_01", ]
    
    
    # Number of participants tested in the pilot test
    sample_size_participants = length(unique(pilot_data[, "participant_ID"]))
    
    
    pilot_data[, "trial_number"] = as.numeric(pilot_data[, "trial_number"])
    
    #extract data from erotic trials 
    data_BF = pilot_data[!is.na(pilot_data[, "trial_number"]) & pilot_data[, "reward_type"] == "erotic", ]
    
    data_BF[,"participant_ID"] = droplevels(data_BF[,"participant_ID"])
    
    
    ######################################################################
    #                                                                    #
    #                          Data analysis                             #
    #                                                                    #
    ######################################################################
    
    
    ######################################################################
    #                       Primary confirmatory test                    #
    ######################################################################
    
    #================================================================#
    #                Calculate number of successes                   #
    #================================================================#
    
    # number of successes and total N of trials
    sides_match = data_BF[,"guessed_side"] == data_BF[,"target_side"]
    successes = sum(sides_match)
    total_N = nrow(data_BF)
    
    
    #================================================================#
    #        Calculating Bayes factors using different priors        #
    #================================================================#
    
    ### Replication Bayes factor, with the Bem 2011 experiment 1 results providing the prior information
    
    BF_replication <- BF01_beta(y = successes, N = total_N, y_prior = y_prior, N_prior = N_prior, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null
    
    
    ### Bayes factor with uniform prior
    # using a non-informative flat prior distribution with alpha = 1 and beta = 1
    
    BF_uniform <- BF01_beta(y = successes, N = total_N, y_prior = 0, N_prior = 0, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null
    
    
    ### Bayes factor with BUJ prior
    # the BUJ prior is calculated from Bem's paper where the prior distribution is defined as a
    # normal distribution with a mean at 0 and 90th percentele is at medium effect size d = 0.5 
    # (we asume that this is one-tailed). Source: Bem, D. J., Utts, J., & Johnson, W. O. (2011). 
    # Must psychologists change the way they analyze their data? Journal of Personality and Social Psychology, 101(4), 716-719.
    # We simulate this in this binomial framework with a one-tailed beta distribution with alpha = 7 and beta = 7.
    # This distribution has 90% of its probability mass under p = 0.712, which we determined 
    # to be equivalent to d = 0.5 medium effect size. We used the formula to convert d to log odds ratio logodds = d*pi/sqrt(3), 
    # found here: Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R. (2009). 
    # Converting Among Effect Sizes. In Introduction to Meta-Analysis (pp. 45-49): John Wiley & Sons, Ltd.
    # Then, log odds ratio vas converted to probability using the formula: p = exp(x)/(1+exp(x))
    # The final equation: exp(d*pi/sqrt(3))/(1+exp(d*pi/sqrt(3)))
    
    BF_BUJ <- BF01_beta(y = successes, N = total_N, y_prior = 6, N_prior = 12, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null
    
    
    #================================================================#
    #                    Main analysis inference                     #
    #================================================================#
    
    # determine inference (supported model) based on the Bayes factors calculated above  
    if(Inference_threshold_BF_low >= max(c(BF_replication, BF_uniform, BF_BUJ))) {
      inference_BF = "M1"
      break} else if(Inference_threshold_BF_high <= min(c(BF_replication, BF_uniform, BF_BUJ))) {
        inference_BF = "M0"
        break} else {inference_BF = "Inconclusive"}
    
    
    ##################################################
    #                 Robustness tests               #
    ##################################################
    
    #================================================================#
    #               Robustness test of BF results with NHST          #
    #================================================================#
    
    # robustness of BF results is tested with NHST proportion tests
    # here we perform both an equality test and an equivalence test to draw statistical inference
    
    
    # equality test
    equality_test_p = prop.test(x = successes, n = total_N, p = M0_prob, alternative = "greater")$p.value
    
    # equivalence test
    equivalence_test_p = prop.test(x = successes, n = total_N,
                                   p = M0_prob+minimum_effect_threshold_NHST, 
                                   alternative = "less")$p.value
    
    # descritives of the frequentist estimate
    pbar = successes/total_N
    SE = sqrt(pbar * (1-pbar)/total_N)
    E = qnorm(.9975) * SE
    proportion_995CI = round(pbar + c(-E , E), 3)
    
    
    
    # making inference decision
    if((Inference_threshold_robustness_NHST > equality_test_p) & 
       (Inference_threshold_robustness_NHST <= equivalence_test_p)){
      inference_robustness_NHST = "M1"} else if((Inference_threshold_robustness_NHST > equivalence_test_p) & 
                                                (Inference_threshold_robustness_NHST <= equality_test_p)){
        inference_robustness_NHST = "M0" 
      } else {inference_robustness_NHST = "Inconclusive"}
    
    
    #=======================================================================#
    #   Robustness test of BF results with Bayesian parameter estimation    #
    #=======================================================================#
    
    # robustness of BF results is tested by calculating HDI of the posteriro distribution and checking its relation to
    # the region of practical equivalence (ROPE), promoted in Kruschke, J. K., & Liddell, T. M. 
    # (2017). The Bayesian New Statistics: Hypothesis testing, estimation, meta-analysis, and power 
    # analysis from a Bayesian perspective. Psychonomic Bulletin & Review, 1-29. 
    
    # calculate posterior distribution using beta distribution
    prior_alpha = y_prior + 1
    prior_beta = N_prior-y_prior+1
    
    posterior_alpha = prior_alpha + successes
    posterior_beta = prior_beta + total_N - successes
    
    scale = seq(0, 1, length = 1001)
    posterior_density = dbeta(scale, posterior_alpha, posterior_beta)
    
    # calculate HDI for the posterior distribution
    # (here we calculate the upper and lower bound of the 90% of the probability mass
    # because we use a one-tailed test. This means that the total probability mass below
    # the upper bound of the 90% HDI will be 95%)
    hdi_result = mode_HDI(scale = scale, density = posterior_density, crit_width = 1-Inference_threshold_robustness_Bayes_Par_Est*2, n_samples = 1e6)
    
    # parameters for decision making
    HDI_lb = hdi_result[2]
    HDI_ub = hdi_result[3]
    
    # making inference decision
    ROPE = M0_prob+minimum_effect_threshold_Bayes_Par_Est
    if(HDI_lb >= ROPE){inference_robustness_Bayes_Par_Est = "M1"
    } else if(HDI_ub <= ROPE){inference_robustness_Bayes_Par_Est = "M0"
    } else {inference_robustness_Bayes_Par_Est = "Inconclusive"}
    
    
    #=======================================================================#
    #      Determine final inference of all robustness tests combined       #
    #=======================================================================#
    
    # the main analysis inference is only robust if all robustness tests came to the same inference as the main test
    inferences = c(inference_robustness_NHST, inference_robustness_Bayes_Par_Est)
    inference_robustness = if(all(inferences == inferences[1])){inferences[1]} else {"mixed"}
    
    Robust = if(inference_BF == inference_robustness){"robust"} else {"not robust"}
    
    
    ######################################################################
    #                                                                    #
    #                   Interpretation and visualization                 #
    #                                                                    #
    ######################################################################
    
    
    # visualize results
    # ALL THREE BF values need to cross the threshold for the interpretation to be accepted
    
    BF_results_for_plotting = cbind(as.data.frame(c(BF_replication, BF_uniform, BF_BUJ)), c("BF_replication", "BF_uniform", "BF_BUJ"))
    names(BF_results_for_plotting) = c("Bayes_factor_01", "BF_type")
    
    plot <- ggplot(BF_results_for_plotting, aes(y = Bayes_factor_01, x = BF_type))+
      geom_point()+
      geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=c(Inference_threshold_BF_high), ymax=c(Inf)), alpha = 0.2, fill=c("pink"))+
      geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=c(Inference_threshold_BF_low), ymax=c(Inference_threshold_BF_high)), alpha = 0.2, fill=c("grey80"))+
      geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=c(0), ymax=c(Inference_threshold_BF_low)), alpha = 0.2, fill=c("lightgreen"))+
      geom_point(size = 3.5, shape = 21, fill = "white")+
      scale_y_log10(limits = c(0.005,200), breaks=c(0.01, Inference_threshold_BF_low, 0.1, 0.33, 0, 3, 10, Inference_threshold_BF_high, 100))+
      geom_hline(yintercept = c(Inference_threshold_BF_low, Inference_threshold_BF_high), linetype = "dashed")+
      geom_text(aes(x=0.5, y=c(100, 1, 0.01), label=c("Supports M0", "Inconclusive", "Supports M1"), angle = 270))
    
    
    
    
    ##################################################
    #       Report of confirmatory analysis          #
    ##################################################
    
    # number of erotic trials performed (without the last participant, whose data might
    # be cut in half because of the interim stopping point)
    numtrials_without_last_part = nrow(data_BF[data_BF[,"participant_ID"] != data_BF[nrow(data_BF),"participant_ID"],])
    # number of missing trials/data points (minus the last participant, who might be cut in half)
    missing_data_points = (length(unique(data_BF[data_BF[,"participant_ID"] != data_BF[nrow(data_BF),"participant_ID"], "participant_ID"])) * trial_size_per_participant) - numtrials_without_last_part 
    
    
    general_text = paste("The study stopped at ", total_N, " erotic trials gathered from a total of ",
                         sample_size_participants, " participants.", " There has been ", missing_data_points, " (", round(missing_data_points/total_N*100,2),"%) missing data points due to unfinished sessions.", 
                         " We observed a total of ", round(mean(sides_match), 4)*100, "% successful guesses within ",
                         total_N, " erotic trials (99.5% CI = ", proportion_995CI[1]*100, "%", ", ", proportion_995CI[2]*100, "%", 
                         "; posterior mode = ", hdi_result[1]*100, "%", ", posterior 90% HDI = ", HDI_lb*100, "%", ", ",
                         HDI_ub*100, "%", ").", sep = "")
    
    robustness_text = if(Robust == "robust"){paste(" The results proved to be robust to different statistical approaches, increasing our confidence in our inference.")} else if(Robust == "not robust"){
      paste(" However, the results did not prove to be robust to different statistical approaches.")}
    
    
    
    
    
    ######################################################################
    #                                                                    #
    #                         Exploratory analysis                       #
    #                                                                    #
    ######################################################################
    # The existence of individual differences between participants will be evaluated 
    # using a frequentist approach where we assume that it is possible that only a small 
    # subgroup of the population we sample from is capable of predicting future events
    # with better than chance accuracy.
    
    # EXPLORATORY ANALYSIS RESULTS WILL NOT AFFECT THE CONCLUSIONS OF OUR STUDY
    
    #=======================================================================#
    #           Comparison of expected and observed distributions           #
    #=======================================================================#
    
    # calculate proportion of successful guesses for each participant in the observed data
    data_BF_split = split(data_BF, f = data_BF[,"participant_ID"])
    success_proportions_empirical = sapply(data_BF_split, function(x) mean(x[,"guessed_side"] == x[,"target_side"]))
    
    
    # samples 1,000,000 participants from a population with H0 success rate
    # this is used for the stochastic dominance test as the null model
    # we call this the theoretical sample, because it approximates the theoretical null model
    sim_null_participant_num = 1000000
    success_proportions_theoretical <- rbinom(sim_null_participant_num, size = trial_size_per_participant, prob=M0_prob)/trial_size_per_participant
    
    # determine possible values of success rates
    possible_success_rates = 0
    for(i in 1:trial_size_per_participant){
      possible_success_rates[i+1] = round(1/(trial_size_per_participant/i), 2)
    }
    possible_success_rates_char = as.character(possible_success_rates)
    success_proportions_theoretical_char_rounded = as.character(round(success_proportions_theoretical, 2))
    
    
    # Round the success rates of each participant to the nearest "possible success rate" value. 
    # This is necessary to allow for missing trials for some participants, which may generate
    # success rates that are different from the "possible success rate" for those completing all trials. 
    success_proportions_empirical_rounded <- sapply(success_proportions_empirical, function(x) possible_success_rates[which.min(abs(possible_success_rates - x))])
    
    success_proportions_empirical_char_rounded = as.character(round(success_proportions_empirical_rounded, 2))
    
    success_rates_theoretical = NA
    for(i in 1:length(possible_success_rates)){
      success_rates_theoretical[i] = sum(success_proportions_theoretical_char_rounded == possible_success_rates_char[i])
    }
    success_rates_theoretical_prop = matrix(success_rates_theoretical/sum(success_rates_theoretical))
    
    
    success_rates_empirical = NA
    for(i in 1:length(possible_success_rates)){
      success_rates_empirical[i] = sum(success_proportions_empirical_char_rounded == possible_success_rates_char[i])
    }
    success_rates_empirical_prop = matrix(success_rates_empirical/sum(success_rates_empirical))
    
    
    
    # plot 
    histogram_plot_data = as.data.frame(c(success_rates_theoretical_prop, success_rates_empirical_prop))
    histogram_plot_data = cbind(histogram_plot_data, factor(c(rep("Expected if M0 is true", length(success_rates_theoretical_prop)), rep("Observed", length(success_rates_empirical_prop)))))
    histogram_plot_data = cbind(histogram_plot_data, factor(rep(possible_success_rates_char, 2)))
    names(histogram_plot_data) = c("proportion", "group", "success")
    
    figure_1 =  ggplot(histogram_plot_data, aes(y = proportion, x = success, group = group))+
      geom_bar(aes(fill = group), alpha = 0.5, stat = "identity", position = "identity")+
      scale_fill_manual(values = c("darkgrey", "black")) +
      xlab("Successful guess rate") +
      ylab("Proportion") +
      theme(panel.border = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "bottom",
            axis.line = element_line(colour = "black", size = 1.2),
            axis.text.x = element_text(angle = 90, color = "black", face = "bold", size = 12, margin = margin(1,0,0,0,"mm"), vjust = 0.5),
            axis.text.y = element_text(face = "bold", color = "black", size = 12),
            axis.title = element_text(size = 16))
    
    
    # earth mover's distance
    emd2d(success_rates_theoretical_prop,success_rates_empirical_prop)
    
    
    print(figure_1)
  })
  
})