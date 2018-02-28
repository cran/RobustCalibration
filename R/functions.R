

mathematical_model_eval<-function(input,theta,simul_type, emulator,math_model){
  #VectorXd cm_obs_cur;
  p_theta=length(theta);
  p_x=dim(input)[2];
  n_testing=dim(input)[1];
  
  if(simul_type==0){ ###emulator
    testing_input=cbind(input, t(matrix(theta,p_theta,n_testing)));
    cm_obs_cur=Sample.rgasp(emulator,testing_input);
  }else if(simul_type==1){
    testing_input=cbind(input, t(matrix(theta,p_theta,n_testing)));
    cm_obs_cur=math_model(input, theta);
  }else if( (simul_type==2) | (simul_type==3) ){
    cm_obs_cur= Mogihammer(input,theta,simul_type);
  }
  return(cm_obs_cur);
}


##tell whether it has discrepancy function or no
post_sample <- function(input, output, R0_list, kernel_type, p_theta, output_weights,
                        par_cur, tilde_lambda,prior_par,theta_range,S, X, have_trend, alpha,sd_proposal,
                        discrepancy_type, simul_type,emulator,math_model){ 
  #method_type=Input_list[[17]];
  if(discrepancy_type=='no-discrepancy'){
    ans=post_sample_no_discrepancy(input, output, R0_list,  p_theta, output_weights,
                                   par_cur, theta_range,S, X, have_trend, alpha,sd_proposal,
                                   discrepancy_type, simul_type,emulator,math_model);
  }else{
    ans=post_sample_with_discrepancy(input, output, R0_list, kernel_type, p_theta, output_weights,
                                     par_cur, tilde_lambda,prior_par,theta_range,S, X, have_trend, alpha,sd_proposal,
                                     discrepancy_type, simul_type,emulator,math_model);
  }
  return(ans);
  
}

##sample the posterior with discrepancy
post_sample_with_discrepancy<-function(input, output, R0_list, kernel_type, p_theta, output_weights,
                             par_cur, tilde_lambda,prior_par,theta_range,S, X, have_trend, alpha,sd_proposal,
                             discrepancy_type,simul_type,emulator,math_model){
  

   mylist=as.list(1:4);
  
  #Initialization
  #MatrixXd input=Input_list[0];
  num_obs=dim(input)[1];
  p_x=dim(input)[2];
  
  #MatrixXd output=Input_list[1];
  #List R0_list=Input_list[2];
  #String   kernel_type=Input_list[3];
  #int      p_theta=Input_list[4];
  #VectorXd area_size=Input_list[6];
  #VectorXd par_cur=Input_list[7]; //initial values
  #double   tilde_lambda=0;  //here if it is GaSP it means tilde_lambda=0
  #int method_type=Input_list[16];  //if it is GaSP, it is 1. if it is S-GaSP, it is 2.
  #cm_eval=simul_type; #this is the computer model
  #Input_list[17]; //this is the computer model
  
  #if(method_type==2){
  #  tilde_lambda=Input_list[8];
  #}
  CL_a_b=prior_par;
  
  CL=CL_a_b[1:p_x];
  a=CL_a_b[p_x+1];
  b=CL_a_b[p_x+2];
  
  #MatrixXd theta_range=Input_list[10];
  #int      S=Input_list[11];
  
  #VectorXd X=Input_list[12];
  #have_mean=have_trend;
  p_theta_m=0;
  
  if(have_trend){
    p_theta_m=dim(X)[2];
    #.cols();
  }
  
  #VectorXd alpha=  Input_list[14];
  #VectorXd sd_all=  Input_list[15];
   sd_theta=sd_proposal[1:p_theta];
   sd_log_beta=sd_proposal[(p_theta+1):(p_theta+p_x)];
   sd_log_eta=sd_proposal[p_theta+p_x+1];
  
  
  inv_output_weights=1/output_weights;
  

  record_par=matrix(0,S,p_theta+p_x+2+p_theta_m);
  record_post=rep(0,S)
  #MatrixXd::Zero(S,p_theta+p_x+2+p_theta_m);
  #VectorXd record_post=VectorXd::Zero(S);
  
  param=par_cur;
  
  
  #MatrixXd L;
  if(discrepancy_type=='GaSP'){
    L=Get_R_new(exp(param[(p_theta+1):(p_theta+p_x)]),exp(param[p_theta+p_x+1]),R0_list,kernel_type,alpha, inv_output_weights);
  }else{
    L=Get_R_z_new(exp(param[(p_theta+1):(p_theta+p_x)]),exp(param[p_theta+p_x+1]),tilde_lambda,R0_list,kernel_type,alpha, inv_output_weights );
  }
  
  #MatrixXd L_propose;
  
  accept_N_theta=0;
  accept_N_beta=0;
  count_dec_record=rep(0,S); # this is to count how many sample points are outside the boundary and get rejected
  
  post_cur=0;
  
  theta_cur=rep(0,p_theta);
  theta_sample=rep(0,p_theta);
  
  #bool decision_0;
  #bool decision;
  param_propose=par_cur;
  xi_sample=rep(0,p_x);
  
  log_eta_sample=0;
  #post_propose;
  r_ratio=0;


  cm_obs_cur=mathematical_model_eval(input,param[1:p_theta],simul_type,emulator,math_model);
  
  
  c_prop=1/4
  
  #start of the sampling
  
  for (i_S in 1:S){
    if(i_S==floor(S*c_prop)){
      #cat(post_cur,'\n')
      cat(c_prop*S, ' of ', S, ' posterior samples are drawn \n')
      c_prop=c_prop+1/4
    }
    #cat(post_cur,'\n')
    #cat(par_cur,'\n')
    #cat(accept_N_theta,'\n')
    
    par_cur[(p_theta+p_x+2):(p_theta+p_x+2+p_theta_m)]=Sample_sigma_2_theta_m(par_cur,L,output,p_theta,p_x,X,have_trend,cm_obs_cur);
    
   # par_cur.segment(p_theta+p_x+1,1+p_theta_m)=Sample_sigma_2_theta_m(par_cur,L,output,p_theta,p_x,X,have_trend,cm_obs_cur);
    
    post_cur=Log_marginal_post(par_cur,L,output,p_theta,p_x,X,have_trend,CL,a,b,cm_obs_cur);
    
    #//sample theta
    theta_cur=par_cur[1:p_theta];
    decision_0=F;
    
    
    # for(int i_theta=0; i_theta<p_theta; ++i_theta){
    #   theta_sample(i_theta)=theta_cur(i_theta)+sd_theta(i_theta)*(theta_range(i_theta,1)-theta_range(i_theta,0))*distribution_stan_norm(generator);
    #   if((theta_sample(i_theta)>theta_range(i_theta,1))|| (theta_sample(i_theta)<theta_range(i_theta,0)) ){
    #     decision_0=true;
    #     count_dec_record(i_S)=1;
    #     break;
    #   }
    #   
    # }
    
    
    for(i_theta in 1:p_theta){
      theta_sample[i_theta]=rnorm(1,mean=theta_cur[i_theta],sd= (sd_theta[i_theta]*(theta_range[i_theta,2]- theta_range[i_theta,1]) ) )  ##here maybe truncated
       #cat(theta_sample[i_theta],'\n')
       #cat(theta_range[i_theta,2],'\n')
       #cat(decision_0,'\n')
        if((theta_sample[i_theta]>theta_range[i_theta,2])| (theta_sample[i_theta]<theta_range[i_theta,1]) ){
          decision_0=T  ##reject directly
        count_dec_record[i_S]=1
        break;
      }
    }


    #//if decision_0==true, then stay at original place
    #//ow see whether we accept the sample
    if(decision_0==F){
      param_propose=par_cur;
      param_propose[1:p_theta]=theta_sample;
      

      cm_obs_propose=mathematical_model_eval(input,theta_sample,simul_type,emulator,math_model);
      
      post_propose=Log_marginal_post(param_propose,L,output,p_theta,p_x,X,have_trend,CL,a,b, cm_obs_propose);
      r_ratio=exp(post_propose-post_cur);
      decision=Accept_proposal(r_ratio);
      if(decision){

        par_cur=param_propose;
        post_cur=post_propose;
        cm_obs_cur=cm_obs_propose;
        accept_N_theta=accept_N_theta+1;
        
        #cat(par_cur,'\n')
        
      }
      
    }
    
    #//sample xi and log_eta 
    for(i_x in 1: p_x){
      xi_sample[i_x]=par_cur[p_theta+i_x]+sd_log_beta[i_x]*rnorm(1);
      #distribution_stan_norm(generator);
    }
    log_eta_sample=par_cur[p_theta+p_x+1]+sd_log_eta*rnorm(1); 
    #distribution_stan_norm(generator);
    
    param_propose=par_cur;
    
    param_propose[(p_theta+1):(p_theta+p_x)]=xi_sample;
    param_propose[p_theta+p_x+1]=log_eta_sample;
    
    param=param_propose;
    
    
    if(discrepancy_type=='GaSP'){
      L_propose=Get_R_new(exp(param[(p_theta+1):(p_theta+p_x)]),exp(param[p_theta+p_x+1]),R0_list,kernel_type,alpha, inv_output_weights);
    }else{
      L_propose=Get_R_z_new(exp(param[(p_theta+1):(p_theta+p_x)]),exp(param[p_theta+p_x+1]),tilde_lambda,R0_list,kernel_type,alpha, inv_output_weights );
    }
    
    post_propose=Log_marginal_post(param,L_propose,output,p_theta,p_x,X,have_trend,CL,a,b,cm_obs_cur);
    
    r_ratio=exp(post_propose-post_cur);
    
    decision=Accept_proposal(r_ratio);
    
    if(decision){
      par_cur=param_propose;
      post_cur=post_propose;
      L=L_propose;
      accept_N_beta=accept_N_beta+1;
    }
    
    record_par[i_S,]=par_cur;
    record_post[i_S]=post_cur;
    
    #record_par.block(i_S,0,1,p_theta+p_x+2+p_theta_m)=par_cur.transpose();
    #record_post(i_S)=post_cur;
  }
  
  cat('Done with posterior sampling \n')
  

  mylist[[1]]=record_par;
  mylist[[2]]=record_post;
  
  mylist[[3]]=c(accept_N_theta,accept_N_beta);
  
  mylist[[4]]=count_dec_record;
  
  cat(accept_N_theta, ' of ', S, ' proposed posterior samples of calibration parameters are accepted \n')
  cat(accept_N_beta, ' of ', S, ' proposed posterior samples of range and nugget parameters are accepted \n')
  
  return(mylist);
  
  
}






post_sample_no_discrepancy<-function(input, output, R0_list,  p_theta, output_weights,
                                par_cur, theta_range,S, X, have_trend, alpha,sd_proposal,
                                discrepancy_type, simul_type,emulator,math_model){
  
  
  mylist=as.list(1:4);
  
  #Initialization
  num_obs=dim(input)[1];
  p_x=dim(input)[2];

  p_theta_m=0;
  
  if(have_trend){
    p_theta_m=dim(X)[2];
  }
  
  #List mylist;
  sd_theta=sd_proposal[1:p_theta];
  
  # simul_type=  Input_list[17];
  

  
  inv_output_weights=1/output_weights;

  
  record_par=matrix(0,S,p_theta+1+p_theta_m);
  record_post=rep(0,S)

  accept_N_theta=0;
  count_dec_record=rep(0,S); # this is to count how many sample points are outside the boundary and get rejected
  
  
  #MatrixXd record_par=MatrixXd::Zero(S,p_theta+1+p_theta_m);
  #VectorXd record_post=VectorXd::Zero(S);
  
  param=par_cur;
  

  post_cur=0;
  
  
  theta_cur=rep(0,p_theta);
  theta_sample=rep(0,p_theta);
  
  #bool decision_0;
  #bool decision;
  param_propose=par_cur;
  r_ratio=0;

  
  
  

  #start of the sampling
  


  #VectorXd cm_obs_cur;
  #VectorXd cm_obs_propose;
  
  #VectorXd emulator_par= VectorXd::Zero(p_x+p_theta);
  #if(simul_type==0){
  #  Input_list[5]; //additional parameter for emulator
  #}
  
  #MatrixXd simul_input=  Input_list[18];
  #VectorXd simul_output=  Input_list[19];
  
  
  cm_obs_cur=mathematical_model_eval(input,param[1:p_theta],simul_type,emulator,math_model);
  
  
  #cm_obs_cur=mathematical_model_eval(param.head(p_theta),input,simul_type,emulator_par,simul_input,simul_output);
  
  c_prop=1/4
  for (i_S in 1:S){
    if(i_S==floor(S*c_prop)){
      #cat(post_cur,'\n')
      cat(c_prop*S, ' of ', S, ' posterior samples are drawn \n')
      c_prop=c_prop+1/4
    }
    #cat(post_cur,'\n')
    
    #cat(par_cur,'\n')
    par_cur[(p_theta+1):(p_theta+1+p_theta_m)]=Sample_sigma_2_theta_m_no_discrepancy(par_cur,output,p_theta,X,have_trend, inv_output_weights,cm_obs_cur);
    
    # par_cur.segment(p_theta+p_x+1,1+p_theta_m)=Sample_sigma_2_theta_m(par_cur,L,output,p_theta,p_x,X,have_trend,cm_obs_cur);
    
    post_cur=Log_marginal_post_no_discrepancy(par_cur,output,p_theta,X,have_trend, inv_output_weights,cm_obs_cur);
    

    #par_cur.segment(p_theta,1+p_theta_m)=Sample_sigma_2_theta_m_no_discrepancy(par_cur,output,p_theta,X,have_mean, inv_output_weights,cm_obs_cur);
    
    #post_cur=Log_marginal_post_no_discrepancy(par_cur,output,p_theta,X,have_mean, inv_output_weights,cm_obs_cur);
    
    #sample theta
    theta_cur=par_cur[1:p_theta];
    decision_0=F;
    
    for(i_theta in 1:p_theta){
      theta_sample[i_theta]=rnorm(1,mean=theta_cur[i_theta],sd=sd_theta[i_theta]*(theta_range[i_theta,2]- theta_range[i_theta,1]) )  ##here maybe truncated
      if((theta_sample[i_theta]>theta_range[i_theta,2])| (theta_sample[i_theta]<theta_range[i_theta,1]) ){
        decision_0=T  ##reject directly
        count_dec_record[i_S]=1
        break;
      }
    }
    
    #//if decision_0==true, then stay at original place
    #//ow see whether we accept the sample

    if(decision_0==F){
      param_propose=par_cur;
      param_propose[1:p_theta]=theta_sample;
      cm_obs_propose=mathematical_model_eval(input,theta_sample,simul_type,emulator,math_model);
      
      post_propose=Log_marginal_post_no_discrepancy(param_propose,output,p_theta,X,have_trend, inv_output_weights,cm_obs_propose);

      #  Log_marginal_post(param_propose,L,output,p_theta,p_x,X,have_trend,CL,a,b, cm_obs_propose);
      r_ratio=exp(post_propose-post_cur);
      decision=Accept_proposal(r_ratio);
      if(decision){
        par_cur=param_propose;
        post_cur=post_propose;
        cm_obs_cur=cm_obs_propose;
        accept_N_theta=accept_N_theta+1;
      }
      
    }
    
    # theta_cur=par_cur.head(p_theta);
    # decision_0=false;
    # 
    # 
    # for(int i_theta=0; i_theta<p_theta; ++i_theta){
    #   theta_sample(i_theta)=theta_cur(i_theta)+sd_theta(i_theta)*(theta_range(i_theta,1)-theta_range(i_theta,0))*distribution_stan_norm(generator);
    #   if((theta_sample(i_theta)>theta_range(i_theta,1))|| (theta_sample(i_theta)<theta_range(i_theta,0)) ){
    #     decision_0=true;
    #     count_dec_record(i_S)=1;
    #     break;
    #   }
    #   
    # }
    # 
    # 
    # 
    # //if decision_0==true, then stay at original place
    # //ow see whether we accept the sample
    # if(decision_0==false){
    #   param_propose=par_cur;
    #   param_propose.head(p_theta)=theta_sample;
    #   
    #   cm_obs_propose=Mathematical_model_eval(theta_sample,input,simul_type,emulator_par,simul_input,simul_output);
    #   
    #   
    #   post_propose=Log_marginal_post_no_discrepancy(param_propose,output,p_theta,X,have_mean, inv_output_weights,cm_obs_propose);
    #   r_ratio=exp(post_propose-post_cur);
    #   decision=Accept_proposal(r_ratio);
    #   if(decision){
    #     par_cur=param_propose;
    #     post_cur=post_propose;
    #     cm_obs_cur=cm_obs_propose;
    #     accept_N_theta=accept_N_theta+1;
    #   }
    #   
    # }
    
    
    record_par[i_S,]=par_cur;
    record_post[i_S]=post_cur;
    
    #record_par.block(i_S,0,1,p_theta+p_x+2+p_theta_m)=par_cur.transpose();
    #record_post(i_S)=post_cur;
  }
  
  cat('Done with posterior sampling \n')
  mylist[[1]]=record_par;
  mylist[[2]]=record_post;
  
  mylist[[3]]=c(accept_N_theta);
  
  mylist[[4]]=count_dec_record;
  
  cat(accept_N_theta, ' of ', S, ' proposed posterior samples samples of calibration parameters are accepted \n')
  
    
  return(mylist);
  
  
}



