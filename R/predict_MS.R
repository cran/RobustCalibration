
##multiple sources
predict_MS.rcalibration_MS<-function(object, testing_input, X_testing=as.list(rep(0,object@num_sources)),
                                 testing_output_weights=NULL, 
                                 n_thinning=10,
                                  interval_est=NULL,interval_data=rep(F,length(testing_input)),math_model=NULL,...){
  

  if(object@measurement_bias==T){
    if(length(testing_input)!=(object@num_sources+1) ){
      stop("please specified the testing input for the model discrepancy. \n")
    }
  }else{
    if(length(testing_input)!=object@num_sources){
      stop("The number of sources in the testing input should match the number of sources in the object. \n")
    }
  }
  
  for(i_source in 1:object@num_sources){
    testing_input[[i_source]]=as.matrix(testing_input[[i_source]])
    ##make it a numeric matrix
    testing_input[[i_source]]=matrix(as.numeric(testing_input[[i_source]]), ncol = ncol(testing_input[[i_source]]))
    
  }
  if(object@measurement_bias==T){
    testing_input[[i_source+1]]=as.matrix( testing_input[[i_source+1]])
    ##make it a numeric matrix
    testing_input[[i_source+1]]=matrix(as.numeric(testing_input[[i_source+1]]), ncol = ncol(testing_input[[i_source+1]]))
    
  }
  
  if(is.null(testing_output_weights)){
    testing_output_weights=list()
    for(i_source in 1: object@num_sources){
      testing_output_weights[[i_source]]=rep(1,dim(testing_input[[i_source]])[1])
    }
  }
  
  
  if(object@measurement_bias==F){
    predict_obj=predict_MS_no_measurement_bias(object,testing_input,X_testing,testing_output_weights,interval_est,interval_data,math_model,n_thinning)
    
  }else{
    predict_obj=predict_MS_with_measurement_bias(object,testing_input,X_testing,testing_output_weights,interval_est,interval_data,math_model,n_thinning)
    
  }
  return(predict_obj)
  
  
}

predict_MS_no_measurement_bias<-function(object, testing_input, X_testing=as.list(rep(0,object@num_sources)),
                                        testing_output_weights=NULL, 
                                         interval_est=NULL,interval_data=rep(F,length(testing_input)),math_model=NULL,n_thinning=10){
  predictobj <- new("predictobj.rcalibration_MS")
  
  emulator=as.list(rep(NA,object@num_sources))
  
  for(i_source in 1:object@num_sources){
    if(object@simul_type[[i_source]]==0){
      if(length(object@emulator_ppgasp[[i_source]]@p>0)){ ##ppgasp
        emulator=object@emulator_ppgasp[[i_source]]
      }else{
        emulator=object@emulator_rgasp[[i_source]]
        
      }
      
     }
   }
    
  
  for(i_source in 1:object@num_sources){
    
    if(object@discrepancy_type[i_source]=='no-discrepancy'){
      record_cm_pred=0
      record_cm_pred_no_mean=0
      c_prop=1/4
      
      for(i_S in (1: SS)*n_thinning ){
        #print(i_S)
        
        if(i_S==floor(SS*n_thinning*c_prop)){
          cat(c_prop*100, 'percent is completed \n')
          c_prop=c_prop+1/4
        }
        
        
        theta=object@post_theta[i_S,] ##shared parameter
        sigma_2_delta=object@post_individual_par[[i_source]][i_S,1] #the first one is the sigma_2
        
        mean_cm_test=mathematical_model_eval(testing_input[[i_source]],theta,object@simul_type[i_source],object@emulator[[i_source]],
                                             object@emulator_type[i_source],object@loc_index_emulator[[i_source]], math_model[[i_source]]);
        mean_cm_test_no_mean=mean_cm_test
        if(object@have_trend[i_source]){
          theta_m=as.matrix(object@post_individual_par[[i_source]][i_S,2:(1+object@q[i_source]) ])
          
          #object@post_sample[i_S,(object@p_theta+2):(object@p_theta+1+object@q)]
          mean_cm_test=mean_cm_test+X_testing[[i_source]]%*%theta_m
        }
        
        
        record_cm_pred=record_cm_pred+mean_cm_test
        record_cm_pred_no_mean=record_cm_pred_no_mean+mean_cm_test_no_mean
        
        if(!is.null(interval_est)){
          qnorm_all=qnorm(interval_est[[i_source]]);
          for(i_int in 1:length(interval_est[[i_source]]) ){
            record_interval[,i_int]=record_interval[,i_int]+mean_cm_test+qnorm_all[i_int]*sqrt(sigma_2_delta/testing_output_weights[[i_source]])
          }
        }
        
      }
      
      record_cm_pred=record_cm_pred/SS
      record_cm_pred_no_mean=record_cm_pred_no_mean/SS
      #quilt.plot(x=input_ascending[,1], y=input_ascending[,2], z=record_cm_pred,nrow = 64, ncol = 64,main='real')
      
      #output.list[[i_source]]=list()
      
      #output.list[[i_source]]$math_model_mean=record_cm_pred
      
      predictobj@math_model_mean[[i_source]]=record_cm_pred
      predictobj@math_model_mean_no_trend[[i_source]]=record_cm_pred_no_mean
      if(!is.null(interval_est)){
        record_interval=record_interval/SS
        #output.list[[i_source]]$interval=record_interval
        predictobj@interval[[i_source]]=record_interval
        #ans.list[[2]]=record_interval
      }
      cat('Source',i_source, 'is completed \n')
      
      # return(output.list)
      
    }else{
      
      if(!is.null(interval_est)){
        c_star_record=rep(0,dim(testing_input[[i_source]])[1])
      }
      
      N_testing=dim(testing_input[[i_source]])[1]
      r0=as.list(1:object@p_x[i_source])
      for(i in 1:object@p_x[i_source]){
        r0[[i]]=abs(outer(object@input[[i_source]][,i],testing_input[[i_source]][,i],'-'))
      }
      
      record_cm_pred=0
      record_cm_pred_no_mean=0
      record_pred_mean=0
      
      SS=floor(dim(object@post_theta)[1]/n_thinning)
      
      #c_star=rep(0, n_testing)
      c_prop=1/4
      
      
      
      for(i_S in (1: SS)*n_thinning ){
        #print(i_S)
        
        if(i_S==floor(SS*n_thinning*c_prop)){
          cat(c_prop*100, 'percent is completed \n')
          c_prop=c_prop+1/4
        }
        
        theta=object@post_theta[i_S,]
        beta_delta=exp(object@post_individual_par[[i_source]][i_S,1:object@p_x[i_source]])
        eta_delta=exp(object@post_individual_par[[i_source]][i_S,object@p_x[i_source]+1])
        sigma_2_delta=object@post_individual_par[[i_source]][i_S,object@p_x[i_source]+2]
        if(object@have_trend[i_source]){
          theta_m=as.matrix(object@post_individual_par[[i_source]][i_S,(object@p_x[i_source]+3):(object@p_x[i_source]+2+object@q[i_source])])
        }
        
        mean_cm=mathematical_model_eval(object@input[[i_source]],theta,object@simul_type[[i_source]],object@emulator[[i_source]],
                                        object@emulator_type[i_source],object@loc_index_emulator[[i_source]],math_model[[i_source]]);
        
        output_minus_cm=object@output[[i_source]]- mean_cm
        
        if(object@have_trend[i_source]){
          output_minus_cm=output_minus_cm-object@X[[i_source]]%*%theta_m
        }
        ##remember 
        if(object@discrepancy_type[i_source]=="GaSP"){
          L=Get_R_new( beta_delta,   eta_delta,  
                       object@R0[[i_source]], object@kernel_type[[i_source]],object@alpha[[i_source]],
                       1/object@output_weights[[i_source]])
        }else if(object@discrepancy_type[i_source]=="S-GaSP"){
          L=Get_R_z_new( beta_delta,   eta_delta,  object@lambda_z[[i_source]][i_S], ##allow sampling of lambda_z
                         object@R0[[i_source]], object@kernel_type[[i_source]],object@alpha[[i_source]],
                         1/object@output_weights[[i_source]])
        }
        
        R_inv_y=backsolve(t(L),forwardsolve(L,output_minus_cm ))
        
        if(object@discrepancy_type[i_source]=="S-GaSP"){
          R_inv_y=Update_R_inv_y(R_inv_y,object@R0[[i_source]],beta_delta,object@kernel_type[i_source],object@alpha[[i_source]],
                                 object@lambda_z[[i_source]][i_S],object@num_obs[i_source])
        }
        
        r=separable_kernel(r0,beta_delta,kernel_type=object@kernel_type[i_source],object@alpha[[i_source]])
        
        rt_R_inv_y=t(r)%*%R_inv_y
        
        mean_cm_test=mathematical_model_eval(testing_input[[i_source]],theta,object@simul_type[i_source],object@emulator[[i_source]],
                                             object@emulator_type[i_source],object@loc_index_emulator[[i_source]],math_model[[i_source]]);
        mean_cm_test_no_mean=mean_cm_test
        #f_M_testing=mean_cm_test
        
        if(object@have_trend[i_source]){
          mean_cm_test=mean_cm_test+X_testing[[i_source]]%*%theta_m
        }
        
        
        pred_mean=mean_cm_test+rt_R_inv_y
        
        record_cm_pred=record_cm_pred+mean_cm_test
        record_cm_pred_no_mean=record_cm_pred_no_mean+mean_cm_test_no_mean
        
        record_pred_mean=record_pred_mean+pred_mean
        
        
        
        
        if(!is.null(interval_est)){
          if(object@discrepancy_type[i_source]=="GaSP"){
            R_tilde_inv_r=(backsolve(t(L),forwardsolve(L, r )))
            
            for(i_testing in 1:N_testing){
              c_star_record[i_testing]=1-r[,i_testing]%*%R_tilde_inv_r[,i_testing]
            }
            var_gasp_f=sigma_2_delta*c_star_record
            if(interval_data[i_source]==T){
              var_gasp_f=var_gasp_f+sigma_2_delta*eta_delta/testing_output_weights[[i_source]]
            }
            
            qnorm_all=qnorm(interval_est[[i_source]]);
            for(i_int in 1:length(interval_est[[i_source]]) ){
              record_interval[,i_int]=record_interval[,i_int]+pred_mean+qnorm_all[i_int]*sqrt(var_gasp_f)
              
            }
            
          }else if(object@discrepancy_type[i_source]=="S-GaSP"){
            R=separable_kernel(object@R0[[i_source]], beta_delta, object@kernel_type[i_source],
                               object@alpha[[i_source]])
            R_tilde=R+1/object@lambda_z[[i_source]][i_S]*diag(object@num_obs[i_source]);  
            
            #system.time(Chol_Eigen(R_tilde))
            
            L_R_tilde=Chol_Eigen(R_tilde)
            
            I_minus_R_R_tilde=diag(object@num_obs[i_source])-t(backsolve(t(L_R_tilde),forwardsolve(L_R_tilde,R )))
            
            
            R_tilde_inv_r=(backsolve(t(L_R_tilde),forwardsolve(L_R_tilde, r )))
            
            r_z=I_minus_R_R_tilde%*%r
            
            R_z_tilde_inv_r_z=(backsolve(t(L),forwardsolve(L, r_z )))
            
            for(i_testing in 1:N_testing){
              c_star_record[i_testing]=1-t(r[,i_testing])%*%R_tilde_inv_r[,i_testing]-r_z[,i_testing]%*%R_z_tilde_inv_r_z[,i_testing]
            }
            
            #t(r[,i_testing])%*%solve(R_tilde)%*%r[,i_testing] 
            
            var_gasp_f=sigma_2_delta*c_star_record
            
            if(interval_data[i_source]==T){
              var_gasp_f=var_gasp_f+sigma_2_delta*eta_delta/testing_output_weights[[i_source]]
            }
            
            qnorm_all=qnorm(interval_est[[i_source]]);
            for(i_int in 1:length(interval_est[[i_source]]) ){
              record_interval[,i_int]=record_interval[,i_int]+pred_mean+qnorm_all[i_int]*sqrt(var_gasp_f)
              
            }
            
          }
        }  
      }
      
      
      
      #quilt.plot(x=input_ascending[,1], y=input_ascending[,2], z=record_cm_pred,nrow = 64, ncol = 64,main='real')
      
      #quilt.plot(x=input_ascending[,1], y=input_ascending[,2], z=record_pred_mean,nrow = 64, ncol = 64,main='real')
      
      record_cm_pred=record_cm_pred/SS
      record_cm_pred_no_mean=record_cm_pred_no_mean/SS
      record_pred_mean=record_pred_mean/SS
      
      predictobj@math_model_mean[[i_source]]=record_cm_pred
      predictobj@math_model_mean_no_trend[[i_source]]=record_cm_pred_no_mean
      predictobj@mean[[i_source]]=record_pred_mean
      
      #output.list[[i_source]]=list()
      
      #output.list[[i_source]]$math_model_mean=record_cm_pred
      #output.list[[i_source]]$mean=record_pred_mean
      
      if(!is.null(interval_est)){
        record_interval=record_interval/SS
        predictobj@interval[[i_source]]=record_interval
        
        #output.list[[i_source]]$interval=record_interval
      }
      
    }
    
  }
  
  return(predictobj)
  
  
}





predict_MS_with_measurement_bias<-function(object, testing_input, X_testing=as.list(rep(0,object@num_sources)),
                                       testing_output_weights=NULL, 
                                         interval_est=NULL,interval_data=rep(F,length(testing_input)),
                                       math_model=NULL, n_thinning=10){
  
    emulator=as.list(rep(NA,object@num_sources))
    
    for(i_source in 1:object@num_sources){
      if(object@simul_type[[i_source]]==0){
        if(length(object@emulator_ppgasp[[i_source]]@p>0)){ ##ppgasp
          emulator=object@emulator_ppgasp[[i_source]]
        }else{
          emulator=object@emulator_rgasp[[i_source]]
          
        }
        
      }
    }
  
    predictobj <- new("predictobj.rcalibration_MS")
  
    if(!is.null(interval_est)){
      #record_interval_delta=matrix(0,dim(testing_input[[object@num_sources+1]]),length(interval_est[[object@num_sources+1]]));
      c_star_record_delta=rep(0,dim(testing_input[[object@num_sources+1]])[1])
    }
    
    SS=dim(object@post_individual_par[[1]])[1]/n_thinning
    
    
    num_sources_1=object@num_sources+1
    num_obs=object@num_obs[1]  ##shared sample size
    
    ##first get the model output ready
    for(i_source in 1:object@num_sources){
      record_cm_pred=0
      record_cm_pred_no_mean=0
      
      N_testing=dim(testing_input[[i_source]])[1]
      
      
      for(i_S in (1: SS)*n_thinning ){
        theta=object@post_theta[i_S,]
        mean_cm_test=mathematical_model_eval(testing_input[[i_source]],theta,object@simul_type[i_source],object@emulator[[i_source]],
                                             object@emulator_type[i_source],object@loc_index_emulator[[i_source]],math_model[[i_source]]);
        mean_cm_test_no_mean=mean_cm_test
        #f_M_testing=mean_cm_test
        
        if(object@have_trend[i_source]){
          theta_m=as.matrix(object@post_individual_par[[i_source]][i_S,(object@p_x[i_source]+3):(object@p_x[i_source]+2+object@q[i_source])])
          mean_cm_test=mean_cm_test+X_testing[[i_source]]%*%theta_m
        }
        record_cm_pred=record_cm_pred+mean_cm_test
        record_cm_pred_no_mean=record_cm_pred_no_mean+mean_cm_test_no_mean
        
      }
      
      record_cm_pred=record_cm_pred/SS
      record_cm_pred_no_mean=record_cm_pred_no_mean/SS
      
      predictobj@math_model_mean[[i_source]]=record_cm_pred
      predictobj@math_model_mean_no_trend[[i_source]]=record_cm_pred_no_mean
      
    }
    
     #quilt.plot(x=(testing_input[[6]][,1]), y=(testing_input[[6]][,2]),
     #          z=predictobj@math_model_mean_no_trend[[3]],nrow = 50, ncol = 50,zlim=c(-0.03,0.05))
    
    cat('Complete the mathematical model prediction \n')
    
    ##model bias

    r0=list()
    for(i in 1:object@p_x[num_sources_1]){
      r0[[i]]=abs(outer(object@input[[num_sources_1]][,i],testing_input[[num_sources_1]][,i],'-'))
    }
    
    delta_mean=0
    var_delta=0;
    for(i_S in (1: SS) ){
      #print(i_S)
      par_cur=object@post_individual_par[[num_sources_1]][i_S,]
      sigma_2_delta=par_cur[object@p_x[num_sources_1]+1]
      
      ##even if it is S-GaSP, when there is no noise, it is the same as GaSP  in prediction
      L_delta=Get_R_new(exp(par_cur[1:object@p_x[num_sources_1]]),0,
                        object@R0[[num_sources_1]],object@kernel_type[num_sources_1],object@alpha[[num_sources_1]], 
                        rep(1,num_obs));
      Sigma_inv_delta=(backsolve(t(L_delta),forwardsolve(L_delta, object@post_delta[i_S,] )))
      
      r=separable_kernel(r0,exp(par_cur[1:object@p_x[num_sources_1]]),kernel_type=object@kernel_type[num_sources_1],object@alpha[[num_sources_1]])
      
      pred_mean= t(r)%*%Sigma_inv_delta
      
      delta_mean=delta_mean+ pred_mean
      
      if(!is.null(interval_est)){
          R_tilde_inv_r=(backsolve(t(L_delta),forwardsolve(L_delta, r )))
          
          for(i_testing in 1:N_testing){
            c_star_record_delta[i_testing]=1-r[,i_testing]%*%R_tilde_inv_r[,i_testing]
          }
          var_delta=var_delta+sigma_2_delta*c_star_record_delta

      }
    }
    

    
    delta_mean=delta_mean/SS

    var_delta=var_delta/SS

    predictobj@delta_mean=delta_mean
    
    #quilt.plot(x=(testing_input[[num_sources_1]][,1]), y=(testing_input[[num_sources_1]][,2]),
    #           z=delta_mean,nrow = 50, ncol = 50,zlim=c(-0.03,0.05))
    
      
    cat('Complete the model discrepancy \n')
    
    ##measurement bias 

    
    for(i_source in 1: object@num_sources){
      #measurement_bias_with_mean=0
      measurement_bias_no_mean=0
      
      for(i in 1:object@p_x[i_source]){
        r0[[i]]=abs(outer(object@input[[i_source]][,i],testing_input[[i_source]][,i],'-'))
      }
      
      N_testing=dim(testing_input[[i_source]])[1]
      
      if(!is.null(interval_est)){
        record_interval=matrix(0,dim(testing_input[[i_source]])[1],length(interval_est[[i_source]]));
        c_star_record=rep(0,dim(testing_input[[i_source]])[1])
        
      }
      
      var_meaurement=0;
      
      for(i_S in 1: SS ){
        #print(i_S)
        theta=object@post_theta[i_S,]
        beta_delta=exp(object@post_individual_par[[i_source]][i_S,1:object@p_x[i_source]])
        eta_delta=exp(object@post_individual_par[[i_source]][i_S,object@p_x[i_source]+1])
        sigma_2_delta=object@post_individual_par[[i_source]][i_S,object@p_x[i_source]+2]
        
        #theta=object@post_theta[i_S,]
        #beta_delta=exp(object@post_individual_par[[i_source]][i_S,1:object@p_x])
        #eta_delta=exp(object@post_individual_par[[i_source]][i_S,object@p_x+1])
        

        mean_cm=mathematical_model_eval(object@input[[i_source]],theta,object@simul_type[[i_source]],object@emulator[[i_source]],
                                        object@emulator_type[i_source],object@loc_index_emulator[[i_source]],math_model[[i_source]]);

        output_minus_cm=object@output[[i_source]]- mean_cm
        
        if(object@have_trend[i_source]){
          theta_m=as.matrix(object@post_individual_par[[i_source]][i_S,(object@p_x[i_source]+3):(object@p_x[i_source]+2+object@q[i_source])])
          output_minus_cm=output_minus_cm-object@X[[i_source]]%*%theta_m
        }
        
        output_minus_cm_minus_delta=output_minus_cm-object@post_delta[i_S,]
        

        if(object@discrepancy_type[i_source]=="GaSP"){
          L=Get_R_new( beta_delta,   eta_delta,  
                       object@R0[[i_source]], object@kernel_type[[i_source]],object@alpha[[i_source]],
                       1/object@output_weights[[i_source]])
        }else if(object@discrepancy_type[i_source]=="S-GaSP"){
          L=Get_R_z_new( beta_delta,   eta_delta,  object@lambda_z[[i_source]][i_S],
                         object@R0[[i_source]], object@kernel_type[[i_source]],object@alpha[[i_source]],
                         1/object@output_weights[[i_source]])
        }
        
        R_inv_y=backsolve(t(L),forwardsolve(L,output_minus_cm_minus_delta ))
        
        if(object@discrepancy_type[i_source]=="S-GaSP"){
          R_inv_y=Update_R_inv_y(R_inv_y,object@R0[[i_source]],beta_delta,object@kernel_type[i_source],object@alpha[[i_source]],
                                 object@lambda_z[[i_source]][i_S],object@num_obs[i_source])
        }
        
        r=separable_kernel(r0,beta_delta,kernel_type=object@kernel_type[i_source],object@alpha[[i_source]])
        
        rt_R_inv_y=t(r)%*%R_inv_y
        
        # quilt.plot(x=(testing_input[[num_sources_1]][,1]), y=(testing_input[[num_sources_1]][,2]),
        #            z=rt_R_inv_y,nrow = 100, ncol = 100,zlim=c(-0.03,0.05))
        
        measurement_bias_no_mean=measurement_bias_no_mean+rt_R_inv_y
        #if(object@have_trend[i_source]){
        #   measurement_bias_with_mean=measurement_bias_with_mean+rt_R_inv_y+X_testing[[i_source]]%*%theta_m
        #}
        
        
        
        if(!is.null(interval_est)){
          if(object@discrepancy_type[i_source]=="GaSP"){
            R_tilde_inv_r=(backsolve(t(L),forwardsolve(L, r )))
            
            for(i_testing in 1:N_testing){
              c_star_record[i_testing]=1-r[,i_testing]%*%R_tilde_inv_r[,i_testing]
            }
            var_meaurement=sigma_2_delta*c_star_record
            if(interval_data[i_source]==T){
              var_meaurement=var_meaurement+sigma_2_delta*eta_delta/testing_output_weights[[i_source]]
            }
            
          }else if(object@discrepancy_type[i_source]=="S-GaSP"){
            R=separable_kernel(object@R0[[i_source]], beta_delta, object@kernel_type[i_source],
                               object@alpha[[i_source]])
            R_tilde=R+1/object@lambda_z[[i_source]][i_S]*diag(object@num_obs[i_source]);  
            
            #system.time(Chol_Eigen(R_tilde))
            
            L_R_tilde=Chol_Eigen(R_tilde)
            
            I_minus_R_R_tilde=diag(object@num_obs[i_source])-t(backsolve(t(L_R_tilde),forwardsolve(L_R_tilde,R )))
            
            
            R_tilde_inv_r=(backsolve(t(L_R_tilde),forwardsolve(L_R_tilde, r )))
            
            r_z=I_minus_R_R_tilde%*%r
            
            R_z_tilde_inv_r_z=(backsolve(t(L),forwardsolve(L, r_z )))
            
            for(i_testing in 1:N_testing){
              c_star_record[i_testing]=1-t(r[,i_testing])%*%R_tilde_inv_r[,i_testing]-r_z[,i_testing]%*%R_z_tilde_inv_r_z[,i_testing]
            }
            
            #t(r[,i_testing])%*%solve(R_tilde)%*%r[,i_testing] 
            
            var_meaurement=sigma_2_delta*c_star_record
            
            if(interval_data[i_source]==T){
              var_meaurement=var_meaurement+sigma_2_delta*eta_delta/testing_output_weights[[i_source]]
            }
            

          }
          
        }
        
      }
      measurement_bias_no_mean=measurement_bias_no_mean/SS
      #measurement_bias_with_mean=measurement_bias_with_mean/SS
      var_meaurement=var_meaurement/SS
      predictobj@measurement_bias_mean[[i_source]]=measurement_bias_no_mean
      
      if(!is.null(interval_est)){
        qnorm_all=qnorm(interval_est[[i_source]]);
        for(i_int in 1:length(interval_est[[i_source]]) ){
          record_interval[,i_int]=record_interval[,i_int]+measurement_bias_no_mean+qnorm_all[i_int]*sqrt(var_meaurement+var_delta)
        }
        
        predictobj@interval[[i_source]]=record_interval
      }
      predictobj@mean[[i_source]]=predictobj@delta_mean+predictobj@measurement_bias_mean[[i_source]]+predictobj@math_model_mean[[i_source]]
      
      #quilt.plot(x=(testing_input[[3]][,1]), y=(testing_input[[3]][,2]),
      #          z= predictobj@mean[[3]],nrow = 50, ncol = 50)
      
    }
        
      
     return(predictobj)
        
}




# 2D lattice code, comment it for now 
# predict_separable_2dim_MS<-function(object, testing_input_separable,
#                                              X_testing=NULL,math_model=NULL,...){
# 
# 
#   if(object@measurement_bias==T){
#     if(length(testing_input_separable)!=(object@num_sources+1) ){
#       stop("please specified the testing input for the model discrepancy. \n")
#     }
#   }else{
#     if(length(testing_input_separable)!=object@num_sources){
#       stop("The number of sources in the testing input should match the number of sources in the object. \n")
#     }
#   }
# 
#   if(object@measurement_bias==F){
#     predict_obj=predict_separable_2dim_MS_no_measurement_bias(object,testing_input_separable,X_testing,math_model)
# 
#   }else{
#     predict_obj=predict_separable_2dim_MS_with_measurement_bias(object,testing_input_separable,X_testing,math_model)
# 
#   }
#   return(predict_obj)
# 
# }
# 
# 
# predict_separable_2dim_MS_no_measurement_bias<-function(object,testing_input_separable,X_testing,math_model){
#   predictobj <- new("predictobj.rcalibration_MS")
#   
#   
#   r_separable=list()
#   r0_separable=list()
#   for(i_source in 1:object@num_sources){
#     SS=floor(dim(object@post_theta)[1])
#     
#       #N_testing=dim(testing_input[[i_source]])[1]
#       #r0=as.list(1:object@p_x[i_source])
#       for(i in 1:object@p_x[i_source]){
#         r0_separable[[i]]=abs(outer(object@input[[i_source]][,i],testing_input_separable[[i_source]][[i]],'-'))
#       }
#       
#       record_cm_pred=0
#       record_cm_pred_no_mean=0
#       record_pred_mean=0
#       
#       SS=floor(dim(object@post_theta)[1])
#       
#       #c_star=rep(0, n_testing)
#       c_prop=1/4
#       
#       
#       
#       for(i_S in (1: SS) ){
#         #print(i_S)
#         
#         if(i_S==floor(SS*c_prop)){
#           cat(c_prop*100, 'percent is completed \n')
#           c_prop=c_prop+1/4
#         }
#         
#         theta=object@post_theta[i_S,]
#         beta_delta=exp(object@post_individual_par[[i_source]][i_S,1:object@p_x[i_source]])
#         eta_delta=exp(object@post_individual_par[[i_source]][i_S,object@p_x[i_source]+1])
#         sigma_2_delta=object@post_individual_par[[i_source]][i_S,object@p_x[i_source]+2]
#         if(object@have_trend[i_source]){
#           theta_m=as.matrix(object@post_individual_par[[i_source]][i_S,(object@p_x[i_source]+3):(object@p_x[i_source]+2+object@q[i_source])])
#         }
#         
#         mean_cm=mathematical_model_eval(object@input[[i_source]],theta,object@simul_type[[i_source]],object@emulator[[i_source]],math_model[[i_source]]);
#         
#         output_minus_cm=object@output[[i_source]]- mean_cm
#         
#         if(object@have_trend[i_source]){
#           output_minus_cm=output_minus_cm-object@X[[i_source]]%*%theta_m
#         }
#         ##remember 
#         if(object@discrepancy_type[i_source]=="GaSP"){
#           L=Get_R_new( beta_delta,   eta_delta,  
#                        object@R0[[i_source]], object@kernel_type[[i_source]],object@alpha[[i_source]],
#                        1/object@output_weights[[i_source]])
#         }else if(object@discrepancy_type[i_source]=="S-GaSP"){
#           L=Get_R_z_new( beta_delta,   eta_delta,  object@lambda_z[[i_source]][i_S],
#                          object@R0[[i_source]], object@kernel_type[[i_source]],object@alpha[[i_source]],
#                          1/object@output_weights[[i_source]])
#         }
#         
#         R_inv_y=backsolve(t(L),forwardsolve(L,output_minus_cm ))
#         
#         if(object@discrepancy_type[i_source]=="S-GaSP"){
#           R_inv_y=Update_R_inv_y(R_inv_y,object@R0[[i_source]],beta_delta,object@kernel_type[i_source],object@alpha[[i_source]],
#                                  object@lambda_z[[i_source]][i_S],object@num_obs[i_source])
#         }
#         
#         
#         for(i_x in 1:2){
#           if(object@kernel_type[i_source]=="matern_5_2"){
#             r_separable[[i_x]]=matern_5_2_funct(r0_separable[[i_x]],beta_delta[i_x])
#           }else if(object@kernel_type[i_source]=="matern_3_2"){
#             r_separable[[i_x]]=matern_3_2_funct(r0_separable[[i_x]],beta_delta[i_x])
#           }else if(object@kernel_type[i_source]=="pow_exp"){
#             r_separable[[i_x]]=pow_exp_funct(r0_separable[[i_x]],beta_delta[i_x],object@alpha[[i_source]][i_x])
#           }
#         }
#         
#         
#         r2_dot_R_inv_y=r_separable[[2]]*as.vector(R_inv_y)
#         
#         rt_R_inv_y=as.vector(t(r_separable[[1]])%*%r2_dot_R_inv_y)
#         
#         #r=separable_kernel(r0,beta_delta,kernel_type=object@kernel_type[i_source],object@alpha[[i_source]])
#         
#         #rt_R_inv_y=t(r)%*%R_inv_y
#         
#         mean_cm_test=mathematical_model_eval(testing_input[[i_source]],theta,object@simul_type[i_source],object@emulator[[i_source]],math_model[[i_source]]);
#         mean_cm_test_no_mean=mean_cm_test
#         #f_M_testing=mean_cm_test
#         
#         if(object@have_trend[i_source]){
#           mean_cm_test=mean_cm_test+X_testing[[i_source]]%*%theta_m
#         }
#         
#         
#         pred_mean=mean_cm_test+rt_R_inv_y
#         
#         record_cm_pred=record_cm_pred+mean_cm_test
#         record_cm_pred_no_mean=record_cm_pred_no_mean+mean_cm_test_no_mean
#         
#         record_pred_mean=record_pred_mean+pred_mean
#         
#         
#         
#       }
#       
#       
#       
#       #quilt.plot(x=input_ascending[,1], y=input_ascending[,2], z=record_cm_pred,nrow = 64, ncol = 64,main='real')
#       
#       #quilt.plot(x=input_ascending[,1], y=input_ascending[,2], z=record_pred_mean,nrow = 64, ncol = 64,main='real')
#       
#       record_cm_pred=record_cm_pred/SS
#       record_cm_pred_no_mean=record_cm_pred_no_mean/SS
#       record_pred_mean=record_pred_mean/SS
#       
#       predictobj@math_model_mean[[i_source]]=record_cm_pred
#       predictobj@math_model_mean_no_trend[[i_source]]=record_cm_pred_no_mean
#       predictobj@mean[[i_source]]=record_pred_mean
#       
#       #output.list[[i_source]]=list()
#       
#       #output.list[[i_source]]$math_model_mean=record_cm_pred
#       #output.list[[i_source]]$mean=record_pred_mean
#       
# 
#     
#   }
# 
#     
#   
#   return(predictobj)
#   
#   
# }
# 
# 
# predict_separable_2dim_MS_with_measurement_bias<-function(object,testing_input_separable,X_testing,math_model){
#   
#   predictobj <- new("predictobj.rcalibration_MS")
#   
# 
#   SS=dim(object@post_individual_par[[1]])[1]
#   
#   
#   num_sources_1=object@num_sources+1
#   num_obs=object@num_obs[1]  ##shared sample size
#   testing_input=list()
#   ##first get the model output ready
#   for(i_source in 1:object@num_sources){
#     #print(i_source)
#     record_cm_pred=0
#     record_cm_pred_no_mean=0
#     
# 
#     testing_input[[i_source]]=expand.grid(testing_input_separable[[i_source]][[1]],testing_input_separable[[i_source]][[2]])
#     testing_input[[i_source]]=as.matrix(testing_input[[i_source]])
#     
#     #N_testing=dim(testing_input[[i_source]][[1]])[1]
#     
#     for(i_S in (1: SS) ){
#       theta=object@post_theta[i_S,]
#       mean_cm_test=mathematical_model_eval(testing_input[[i_source]],theta,object@simul_type[i_source],object@emulator[[i_source]],math_model[[i_source]]);
#       mean_cm_test_no_mean=mean_cm_test
#       #f_M_testing=mean_cm_test
#       
#       if(object@have_trend[i_source]){
#         theta_m=as.matrix(object@post_individual_par[[i_source]][i_S,(object@p_x[i_source]+3):(object@p_x[i_source]+2+object@q[i_source])])
#         mean_cm_test=mean_cm_test+X_testing[[i_source]]%*%theta_m
#       }
#       record_cm_pred=record_cm_pred+mean_cm_test
#       record_cm_pred_no_mean=record_cm_pred_no_mean+mean_cm_test_no_mean
#       
#     }
#     
#     record_cm_pred=record_cm_pred/SS
#     record_cm_pred_no_mean=record_cm_pred_no_mean/SS
#     
#     predictobj@math_model_mean[[i_source]]=record_cm_pred
#     predictobj@math_model_mean_no_trend[[i_source]]=record_cm_pred_no_mean
#     
#   }
#   
#   #quilt.plot(x=(testing_input_separable[[6]][[1]]), y=(testing_input_separable[[6]][[2]]),
#   #           z=predictobj@math_model_mean_no_trend[[3]],nrow = 50, ncol = 50,zlim=c(-0.03,0.05))
#   
#   #image2D(t(matrix(predictobj@math_model_mean_no_trend[[3]],541,469)))
#   cat('Complete the mathematical model prediction \n')
#   
#   
#   ##model bias
#   
#   r0_separable=list()
#   for(i_x in 1:2){
#     r0_separable[[i_x]]=as.matrix(abs(outer(object@input[[num_sources_1]][,i_x],testing_input_separable[[num_sources_1]][[i_x]],'-')))
#   }
#   
#   delta_mean=0
#   #var_delta=0;
#   r_separable=list()
#   
#   for(i_S in (1: SS) ){
#     #print(i_S)
#     par_cur=object@post_individual_par[[num_sources_1]][i_S,]
#     beta_delta=exp(par_cur[1:object@p_x[num_sources_1]])
#     sigma_2_delta=par_cur[object@p_x[num_sources_1]+1]
#     
#     ##even if it is S-GaSP, when there is no noise, it is the same as GaSP  in prediction
#     L_delta=Get_R_new(beta_delta,0,
#                       object@R0[[num_sources_1]],object@kernel_type[num_sources_1],object@alpha[[num_sources_1]], 
#                       rep(1,num_obs));
#     Sigma_inv_delta=(backsolve(t(L_delta),forwardsolve(L_delta, object@post_delta[i_S,] )))
#     
#     ##only work for two dimensions
#     for(i_x in 1:2){
#       if(object@kernel_type[num_sources_1]=="matern_5_2"){
#         r_separable[[i_x]]=matern_5_2_funct(r0_separable[[i_x]],beta_delta[i_x])
#       }else if(object@kernel_type[num_sources_1]=="matern_3_2"){
#         r_separable[[i_x]]=matern_3_2_funct(r0_separable[[i_x]],beta_delta[i_x])
#       }else if(object@kernel_type[num_sources_1]=="pow_exp"){
#           r_separable[[i_x]]=pow_exp_funct(r0_separable[[i_x]],beta_delta[i_x],object@alpha[[num_sources_1]][i_x])
#       }
#     }
#     
# 
#     r2_dot_R_inv_y=r_separable[[2]]*as.vector(Sigma_inv_delta)
#     
#     
#     pred_mean=as.vector(t(r_separable[[1]])%*%r2_dot_R_inv_y)
#     
# 
#     #pred_mean= t(r)%*%Sigma_inv_delta
#     
#     delta_mean=delta_mean+ pred_mean
#     
#   }
#   
#   
#   
#   delta_mean=delta_mean/SS
#   
#   predictobj@delta_mean=delta_mean
#   
#   #image2D(t(matrix(delta_mean,541,469)))
#   
#   #quilt.plot(x=(testing_input[[num_sources_1]][,1]), y=(testing_input[[num_sources_1]][,2]),
#   #           z=delta_mean,nrow = 50, ncol = 50,zlim=c(-0.03,0.05))
#   
#   
#   cat('Complete the model discrepancy \n')
#   
#   
#   ##measurement bias 
#   
#   for(i_source in 1: object@num_sources){
#     #print(i_source)
#     #measurement_bias_with_mean=0
#     measurement_bias_no_mean=0
#     
#     for(i in 1:object@p_x[i_source]){
#       r0_separable[[i]]=abs(outer(object@input[[i_source]][,i],testing_input_separable[[i_source]][[i]],'-'))
#     }
#     
#     #N_testing=dim(testing_input[[i_source]])[1]
#     
#     for(i_S in 1: SS ){
#       #print(i_S)
#       theta=object@post_theta[i_S,]
#       beta_delta=exp(object@post_individual_par[[i_source]][i_S,1:object@p_x[i_source]])
#       eta_delta=exp(object@post_individual_par[[i_source]][i_S,object@p_x[i_source]+1])
#       #sigma_2_delta=object@post_individual_par[[i_source]][i_S,object@p_x[i_source]+2]
# 
#       mean_cm=mathematical_model_eval(object@input[[i_source]],theta,object@simul_type[[i_source]],object@emulator[[i_source]],math_model[[i_source]]);
#       
#       output_minus_cm=object@output[[i_source]]- mean_cm
#       
#       if(object@have_trend[i_source]){
#         theta_m=as.matrix(object@post_individual_par[[i_source]][i_S,(object@p_x[i_source]+3):(object@p_x[i_source]+2+object@q[i_source])])
#         output_minus_cm=output_minus_cm-object@X[[i_source]]%*%theta_m
#       }
#       
#       output_minus_cm_minus_delta=output_minus_cm-object@post_delta[i_S,]
#       
#       
#       if(object@discrepancy_type[i_source]=="GaSP"){
#         L=Get_R_new( beta_delta,   eta_delta,  
#                      object@R0[[i_source]], object@kernel_type[[i_source]],object@alpha[[i_source]],
#                      1/object@output_weights[[i_source]])
#       }else if(object@discrepancy_type[i_source]=="S-GaSP"){
#         L=Get_R_z_new( beta_delta,   eta_delta,  object@lambda_z[[i_source]][i_S],
#                        object@R0[[i_source]], object@kernel_type[[i_source]],object@alpha[[i_source]],
#                        1/object@output_weights[[i_source]])
#       }
#       
#       R_inv_y=backsolve(t(L),forwardsolve(L,output_minus_cm_minus_delta ))
#       
#       if(object@discrepancy_type[i_source]=="S-GaSP"){
#         R_inv_y=Update_R_inv_y(R_inv_y,object@R0[[i_source]],beta_delta,object@kernel_type[i_source],object@alpha[[i_source]],
#                                object@lambda_z[[i_source]][i_S],object@num_obs[i_source])
#       }
#       
#       for(i_x in 1:2){
#         if(object@kernel_type[i_source]=="matern_5_2"){
#           r_separable[[i_x]]=matern_5_2_funct(r0_separable[[i_x]],beta_delta[i_x])
#         }else if(object@kernel_type[i_source]=="matern_3_2"){
#           r_separable[[i_x]]=matern_3_2_funct(r0_separable[[i_x]],beta_delta[i_x])
#         }else if(object@kernel_type[i_source]=="pow_exp"){
#           r_separable[[i_x]]=pow_exp_funct(r0_separable[[i_x]],beta_delta[i_x],object@alpha[[i_source]][i_x])
#         }
#       }
#       
#       
#       r2_dot_R_inv_y=r_separable[[2]]*as.vector(R_inv_y)
#       
#       pred_mean=as.vector(t(r_separable[[1]])%*%r2_dot_R_inv_y)
#       
#       #r=separable_kernel(r0,beta_delta,kernel_type=object@kernel_type[i_source],object@alpha[[i_source]])
#       
#       #rt_R_inv_y=t(r)%*%R_inv_y
#       
#       # quilt.plot(x=(testing_input[[num_sources_1]][,1]), y=(testing_input[[num_sources_1]][,2]),
#       #            z=rt_R_inv_y,nrow = 100, ncol = 100,zlim=c(-0.03,0.05))
#       
#       #image2D(t(matrix(pred_mean,541,469)))
#       
#       measurement_bias_no_mean=measurement_bias_no_mean+pred_mean
#     
#     }
#     
#     measurement_bias_no_mean=measurement_bias_no_mean/SS
#     #measurement_bias_with_mean=measurement_bias_with_mean/SS
#     #var_meaurement=var_meaurement/SS
#     predictobj@measurement_bias_mean[[i_source]]=measurement_bias_no_mean
#     
#     predictobj@mean[[i_source]]=predictobj@delta_mean+predictobj@measurement_bias_mean[[i_source]]+predictobj@math_model_mean[[i_source]]
#   }
#   
#  # image2D((matrix(predictobj@mean[[1]],541,469))
#   
#   
#   
#   return(predictobj)
#   
# }



