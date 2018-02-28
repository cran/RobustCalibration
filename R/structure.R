
## rcalibration Class


setClass("rcalibration", 		
         representation( 
           p_x = "integer",          ## dimension of the inputs
           p_theta="integer",        ## dimention of the calibration parameters
           num_obs = "integer",          ## experimental observations number
           
           ## data
           input = "matrix",           ## the design of experiments, size n x p_x
           output = "vector",           ## the experimental observations, size n x 1
           X="matrix",                  ## mean basis for experiment, size n x q
           have_trend="logical",         ## have mean or no
           q="integer",                 ## number of mean basis for experiment
           R0="list",                    ##abs difference of each type of input
           kernel_type="character",       #####type of kernel to specify
           alpha="vector",              ####roughness parameter in the kernel, only useful for pow_exp
           theta_range= "matrix",        ## the range of calibration parameters
           tilde_lambda="numeric",        ## parameter in the S-GaSP
           S="integer",                   ## number of MCMC
           S_0="integer",                 ## number of burn-in
           prior_par="vector",            ## prior parameters for jointly robust prior
           output_weights="vector",       ##whether the output contains weights
           sd_proposal="vector",            ##standard deviation of the MCMC, if we have a discrepancy, 
                                         ##the size is p_theta+p_x+1; if we don't, the size is p_theta
           discrepancy_type="character",    ## no-discrepancy,  GaSP and  S-GaSP
           simul_type="integer",          ## 0 means use RobustGaSP pacakge to fit the computer model, 
                                         ## 1 means the simulator is defined by the user  
                                         ## 2 means the simulator for Kilauea Volcano by the ascending-mode image 
                                        ##  3 means the simulator for Kilauea Volcano by the decending-mode image
           post_sample="matrix",    ## posterior samples after burn-in
           post_value="vector",             ## posterior value after burn-in after burn-in
           accept_S="vector",                ## number of proposed samples of the calibation parameters are accepted
           count_boundary="numeric"          ## number of proposed samples is outside the boundary of calibration parameters
         ), 
         # validity = function(object) {
         #   if (object@num_obs <= object@p) {
         #     return("the number of experiments must be larger than the spatial dimension")
         #   }
         #   
         #   if (ncol(object@output) != 1) {
         #     return("the response must have one dimension")
         #   }
         #   
         #   if (!identical(nrow(object@input), nrow(object@output))) {
         #     return("the number of observations is not equal to the number of experiments")
         #   }
         #   
         #   TRUE
         # }
)

# 
if(!isGeneric("show")) {
  setGeneric(name = "show",
             def = function(object) standardGeneric("show")
  )
}

setMethod("show", "rcalibration",
          function(object){
            show.rcalibration(object)
          }
)
# 
if(!isGeneric("predict")) {
  setGeneric(name = "predict",
             def = function(object, ...) standardGeneric("predict")
  )
}

setMethod("predict", "rcalibration",
          definition=function(object, testing_input,X_testing=matrix(0,dim(testing_input)[1],1),
                              n_thinning=10, testing_output_weights=rep(1,dim(testing_input)[1]), 
                              interval_est=NULL,interval_data=F, emulator=NULL,math_model=NULL,...){
            predict.rcalibration(object = object, testing_input = testing_input, X_testing=X_testing,
                          n_thinning=n_thinning, testing_output_weights=testing_output_weights, 
                          interval_est=interval_est,interval_data=interval_data, emulator=emulator,math_model=math_model,...)
          }
)

# 
# 
# 
# if(!isGeneric("Sample")) {
#   setGeneric(name = "Sample",
#              def = function(object, ...) standardGeneric("Sample")
#   )
# }
# 
# setMethod("Sample", "rgasp",
#           definition=function(object, testing_input, num_sample=1,testing_trend=matrix(1,dim(testing_input)[1],1), ...) {
#             Sample.rgasp(object = object, testing_input = testing_input, num_sample=num_sample,
#                           testing_trend=testing_trend , ...)
#           }
# )


