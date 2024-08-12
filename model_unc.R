
#development and validation samples


eventid = rownames(partA[partA$DAY30 == 1, ])
noneventid = rownames(partA[partA$DAY30 == 0, ])


devsample = partA[c(sample(eventid, 70), sample(noneventid, 930)), 
                  c("DAY30", "AGE","SHO", "ANT", "DIA",
                    "HYP","HRT","TTR","SEX") ]    # 1000 patients



devsample2h = partA[c(sample(eventid, 14), sample(noneventid, 186)), 
                    c("DAY30", "AGE","SHO", "ANT", "DIA",
                      "HYP","HRT","TTR","SEX")]   #200 patients



devsample5k = partA[c(sample(eventid, 350), sample(noneventid, 4650)), 
                    c("DAY30", "AGE","SHO", "ANT", "DIA",
                      "HYP","HRT","TTR","SEX")]   #5000 patients



indv5k = sample(1:nrow(partB), 5000)
valsample5k = partB[indv5k, c("DAY30", "AGE","SHO", "ANT", "DIA",
                              "HYP","HRT","TTR","SEX")]



# random forest ####


rf_grid = expand.grid(ntrees = c(500, 1000), 
                      mtry = c(2,3,4),
                      minnodesize = c(1,3,5,10))


sampvar_rf = function(B, devs, vals, grid){  
  
  print(paste("Started running at",format(Sys.time(), "%H:%M:%S")))
  
    # N*M dataframe for predictions from OG models
    M = nrow(grid)
    N = nrow(vals) # N is for validation!
  
    preds = data.frame(matrix(NA, nrow = N, ncol = M))
    
    for (m in 1:M){
      model = ranger(DAY30 ~ ., data = devs, num.trees = grid[m, 1],
                     mtry = grid[m, 2],
                     min.node.size = grid[m, 3], probability = TRUE)
      preds[ ,m] = predict(model, vals )$predictions[ ,2]
      
      if(m == M){
        print(paste("Preds df is done at",format(Sys.time(), "%H:%M:%S")))
      }
    }
    
    preds$ID = rownames(vals)
    
    # M long list with N*B dataframe for BS predictions
    bs = list()
    
    for (m in 1:M) {
      bs[[m]] = data.frame(matrix(NA, nrow = N, ncol = B ))
    }
    
    for (b in 1:B){
      b_ind = sample(1:nrow(devs), nrow(devs), replace = TRUE)
      bootsample = devs[b_ind, ]
      
      for (m in 1:M){
        
        bootmod =  ranger(DAY30 ~ ., data = bootsample, num.trees = grid[m,1],
                          mtry = grid[m,2],
                          min.node.size = grid[m,3], probability = TRUE)
        bs[[m]][ ,b] = predict(bootmod, vals )$predictions[ ,2]
      }
      
      if(b == 50){
        paste("First 50 BS iterations are done at",
              format(Sys.time(), "%H:%M:%S")) %>% print()
      }
      
      if(b == 100){
        paste("First 100 BS iterations are done at",
              format(Sys.time(), "%H:%M:%S")) %>% print()
      }
      
      if(b == 150){
        paste("First 150 BS iterations are done at",
              format(Sys.time(), "%H:%M:%S")) %>% print()
      }
      
      if(b == B){
        paste("All BS iterations are done at",
              format(Sys.time(), "%H:%M:%S")) %>% print()
      }
    }
    
    #sampling variances for n patients across m models
    var = data.frame(matrix(NA, nrow = N, ncol = M))
    
    for(m in 1:M){
      var[ ,m] = apply(bs[[m]], MARGIN = 1, FUN = var)
    }
    
    
    
    #standard error for n patient across m models
    se = data.frame(matrix(NA, nrow = N, ncol = M))
    
    for(m in 1:M){
      se[ ,m] = apply(bs[[m]], MARGIN = 1, FUN = sd)
    }
    
    
    
 
    return(list(preds = preds, var = var, se = se))
}




taudist = function(selist, N = 5000, M = 24){ 
  #this function is applicable to all model classes
  
  tau2 = rep(NA, N)
  muhat = rep(NA, N)

  for (n in 1:N){
    
    #some knn models have 0 sampling variance, RMA cannot run on them so they need to be filtered
    nnz_cells = which(selist$var[n, ] != 0)
    
    if (length(nnz_cells) > 0){
    
    re_mod = rma(yi = as.numeric(as.vector(selist$preds[n, nnz_cells])), 
                 vi = as.numeric(as.vector(selist$var[n, nnz_cells])),
                 control=list(stepadj=0.2, maxiter=10000))
    tau2[n] = re_mod$tau2
    muhat[n] = re_mod$b
    
    }
    
    else{
      tau2[n] = 0
      muhat[n] = 0
    }

  }
  
  tau = sqrt(tau2)
  
  sigma_pp = apply(selist$se[ ,1:M], MARGIN = 1, FUN = mean) 
  sigma_avgm = mean(sigma_pp)
  
  return(list(tau = tau,
              meantau = mean(tau), sdtau = sd(tau), 
              sigma_pp = sigma_pp,
              sigma_avgm = sigma_avgm,
              maxtau = max(tau), maxsigma = max(sigma_pp),
              muhat = muhat))
  
}





tauplots = function(selist, taulist, N = 5000, M = 24, family){ 
  #this function is applicable to all model classes
  
  predsgg = selist$preds[selist$preds$ID %in%
                           c("36857","17553","6388","15077", 
                             "570","37660", "25409", "39861", 
                             "20809", "1783", "1095", "31058"), ]
  
  
  #checking the spread of preds
  gg1 = ggplot(data = melt(predsgg, 
                           id.vars = c("ID")))+
    geom_jitter(aes(x = reorder(ID, value, mean), y = value), 
                alpha = 0.5, width = 0.25)+
    ylab(paste("predictions from", M , family, "models"))+
    xlab("patient ID")
    #scale_y_continuous( limits = c(0,0.5))
  
  #dist of tau
  gg2 = ggplot(data = as.data.frame(cbind(taulist$tau, log(taulist$tau))))+
    geom_density(aes(x = V1))+
    xlab("tau")
  
  
  #dist of tau and sigma
  gg3 = ggplot(data = as.data.frame(cbind(taulist$tau, taulist$sigma_pp)))+
    geom_density(aes(x = V1, color = "model_tau"), linewidth = 1)+
    geom_density(aes(x = V2, color = "sampling_sigma" ), linewidth = 1)+
    scale_color_manual(values = c("model_tau" = "lightblue",
                                  "sampling_sigma" = "darkolivegreen3"))+
    xlab("")
    #scale_x_continuous(limits = c(0,0.05))
  
  #dist of tau and sigma log
  gg4 = ggplot(data = as.data.frame(cbind(taulist$tau, taulist$sigma_pp)))+
    geom_density(aes(x = log(V1), color = "model_tau"), linewidth = 1)+
    geom_density(aes(x = log(V2), color = "sampling_sigma" ), linewidth = 1)+
    scale_color_manual(values = c("model_tau" = "lightblue",
                                  "sampling_sigma" = "darkolivegreen3"))+
    xlab("")
  
  
  return(list(preds = gg1, tau = gg2, tausigma = gg3,
              tausigmalog = gg4))
}




##first run - dev2h ####

      rf_1 = sampvar_rf(B = 200, devs = devsample2h, vals = valsample5k, 
                        grid = rf_grid) 

      tau_rf_1 = taudist(rf_1)
      
      plots_rf_1 = tauplots(rf_1, tau_rf_1, family = "Random forest")



##second run - dev1000 ####

      
    rf_2 = sampvar_rf(B = 200, devs = devsample, vals = valsample5k,
                      grid = rf_grid)


    tau_rf_2 = taudist(rf_2)


    plots_rf_2 = tauplots(rf_2, tau_rf_2, family = "Random forest")


## third run - dev5k ####

    rf_3 = sampvar_rf(B = 200, devs = devsample5k, vals = valsample5k,
                      grid = rf_grid)
    
    tau_rf_3 = taudist(rf_3)
    
    plots_rf_3 = tauplots(rf_3, tau_rf_3, family = "Random forest")


    
#regression methods####
    
    
sampvar_reg = function(B, devs, vals){
  
  paste("Started running at",format(Sys.time(), "%H:%M:%S")) %>% print()
  
  mods = list()
  mods[[1]] = glm(DAY30 ~ ., data = devs, family = "binomial")
  mods[[2]] = step(mods[[1]], direction = "backward", trace = 0)
  mods[[3]] = step(mods[[1]], direction = "forward", trace = 0)
  mods[[4]] = glm(DAY30 ~ I(AGE > 65) + SHO + ANT + DIA + HYP + HRT + TTR + SEX,
                  data = devs, family = "binomial")
  
  formula = formula(mods[[1]])
  
  mods[[5]] = glm(update(formula, . ~ . + I(AGE^2)), 
                  data = devs, family = "binomial")
  mods[[6]] = glm(update(formula, . ~ . + I(AGE^0.5)),
                  data = devs, family = "binomial")
  mods[[7]] = glm(update(formula, . ~ . + I(log(AGE))),
                  data = devs, family = "binomial")
  
  mods[[8]] = glm(update(formula, . ~ . + AGE:SEX),
                  data = devs, family = "binomial")
  
  mods[[9]] = glm(update(formula, . ~ . + SHO:ANT),
                  data = devs, family = "binomial")
  
  mods[[10]] = glm(update(formula, . ~ . + SHO:HRT),
                   data = devs, family = "binomial")
  
  mods[[11]] = glm(update(formula, . ~ . + AGE:HYP),
                   data = devs, family = "binomial")
  
  mods[[12]] = glm(update(formula, . ~ . + SHO:HYP),
                   data = devs, family = "binomial")
  
  
  mods[[13]] = gam(DAY30 ~ bs(AGE, df = 3) + SHO + ANT + DIA + HYP + HRT + TTR + SEX,
                   data = devs, family = "binomial")
  mods[[14]] = gam(DAY30 ~ bs(AGE, df = 4) + SHO + ANT + DIA + HYP + HRT + TTR + SEX,
                   data = devs, family = "binomial")
  
  mods[[15]] = gam(DAY30 ~ ns(AGE, df = 3 ) + SHO + ANT + DIA + HYP + HRT + TTR + SEX,
                   data = devs, family = "binomial")
  mods[[16]] = gam(DAY30 ~ ns(AGE, df = 4 ) + SHO + ANT + DIA + HYP + HRT + TTR + SEX,
                   data = devs, family = "binomial")
  mods[[17]] = gam(DAY30 ~ ns(AGE, df = 5 ) + SHO + ANT + DIA + HYP + HRT + TTR + SEX,
                   data = devs, family = "binomial")
  
  mods[[18]] = gam(DAY30 ~ s(AGE, df = 4) + SHO + ANT + DIA + HYP + HRT + TTR + SEX,
                   data = devs, family = "binomial")
  
  
  
  
  
  x = model.matrix(DAY30 ~ ., devs)
  y = devs$DAY30
  
  mods[[19]] = cv.glmnet(x = x, y = y, family = "binomial",
                         alpha = 1 )
  
  mods[[20]] = cv.glmnet(x = x, y = y, family = "binomial",
                         alpha = 0 )
  
  mods[[21]] = cv.glmnet(x = x, y = y, family = "binomial",
                         alpha = 0.5)
  
  
      
  N = nrow(vals)
  M = length(mods) + 3
      
      # n*m dataframe for predictions from OG models
      preds = data.frame(matrix(NA, nrow = N, ncol = M))
      
      for (m in 1:18){
        preds[ ,m] = predict(mods[[m]], newdata = vals, 
                             type = "response")
      }
      for (m in 19:21){
        preds[ ,m] = predict(mods[[m]], newx = model.matrix(DAY30 ~., vals),
                             type = "response", s = "lambda.min")[ ,1]
      }
      for(m in 22:24){
        preds[ ,m] = predict(mods[[m-3]], newx = model.matrix(DAY30 ~., vals),
                             type = "response", s = "lambda.1se")[ ,1]
        
        if (m == 24){
          paste("Preds df is done at",
                format(Sys.time(), "%H:%M:%S")) %>% print()
        }
      }
      
      preds$ID = rownames(vals)
      
      
      
      # m long list with n*B dataframe for BS predictions
      bs = list()
      
      for (m in 1:M) {
        bs[[m]] = data.frame(matrix(NA, nrow = N, ncol = B ))
      }
      
      
      for (b in 1:B) {
        b_ind = sample(1:nrow(devs), nrow(devs), replace = TRUE)
        bootsample = devs[b_ind, ]
        for (m in 1:12){
          bootmod =  glm(formula(mods[[m]]), data = bootsample, 
                         family = "binomial")
          bs[[m]][ ,b] = predict(bootmod, vals, type = "response")
        }
        for(m in 13:18){
          bootmod = gam(formula(mods[[m]]), data = bootsample, 
                        family = "binomial")
          bs[[m]][ ,b] = predict(bootmod, vals, type = "response")
        }
        for (m in 19:21){
          bootmod = cv.glmnet(x = model.matrix(DAY30~., bootsample), 
                              y=bootsample$DAY30, family = "binomial",
                              alpha = mods[[m]]$call$alpha )
          bs[[m]][ ,b] = predict(bootmod, newx = model.matrix(DAY30 ~., vals),
                                 type = "response", s = "lambda.min")[ ,1]
          bs[[m+3]][ ,b] = predict(bootmod, newx = model.matrix(DAY30 ~., vals),
                                   type = "response", s = "lambda.1se")[ ,1]
        }
        
        
        if(b == 50){
          paste("First 50 BS iterations are done at",
                format(Sys.time(), "%H:%M:%S")) %>% print()
        }
        
        if(b == 100){
          paste("First 100 BS iterations are done at",
                format(Sys.time(), "%H:%M:%S")) %>% print()
        }
        
        if(b == 150){
          paste("First 150 BS iterations are done at",
                format(Sys.time(), "%H:%M:%S")) %>% print()
        }
        
        if(b == B){
          paste("All BS iterations are done at",
                format(Sys.time(), "%H:%M:%S")) %>% print()
        }
        
      }
      
      
      #sampling variances for n patients across M models
      var = data.frame(matrix(NA, nrow = N, ncol = M))
      
      for(m in 1:M){
        var[ ,m] = apply(bs[[m]], MARGIN = 1, FUN = var)
      }
      
    
      
      #standard error for n patient across m models
      se = data.frame(matrix(NA, nrow = N, ncol = M))
      
      for(m in 1:M){
        se[ ,m] = apply(bs[[m]], MARGIN = 1, FUN = sd)
      }
      
      
      return(list(preds = preds, var = var, se = se))
      
}



## LR first run - dev2h####

    lr_1 = sampvar_reg(B = 200, devs = devsample2h, vals = valsample5k)
    
    tau_lr_1 = taudist(lr_1)
    
    plots_lr_1 = tauplots(lr_1, tau_lr_1, family = "regression")
    

    
##LR second run - dev1000 ####
    
    lr_2 = sampvar_reg(B = 200, devs = devsample, vals = valsample5k)
    
    tau_lr_2 = taudist(lr_2)
    
    plots_lr_2 = tauplots(lr_2, tau_lr_2, family = "regression")
    


##LR third run - dev5k ####

    lr_3 = sampvar_reg(B = 200, devs = devsample5k, vals = valsample5k)
    
    tau_lr_3 = taudist(lr_3)
    
    plots_lr_3 = tauplots(lr_3, tau_lr_3, family = "regression")
    



#XGBoost ####
    
    
data4xgb = function(sample){
      numerics = lapply(sample, function(x) as.numeric(as.character(x))) %>% as.data.frame()
      
      xlist = list()
      xlist[["data"]] = as(as.matrix(numerics[ ,2:ncol(numerics)]), "dgCMatrix")
      xlist[["label"]] = numerics$DAY30
      
      return(xlist)
    }
    
    
devsample_xgb = data4xgb(devsample)
devsample2h_xgb = data4xgb(devsample2h)
devsample5k_xgb = data4xgb(devsample5k)


valsample5k_xgb = data4xgb(valsample5k)

    
xgb_grid = expand.grid(nrounds = c(500, 1000),
                           max.depth = c(3,4,5),
                           eta = c(0.1,0.01),
                           colsample_bytree = c(0.75, 1))
    
    
sampvar_xgb = function(B, devs, vals, grid){
  
  paste("Started running at",format(Sys.time(), "%H:%M:%S"))  %>% print()
  
  N = nrow(vals$data)
  M = nrow(grid)
      
      #n*m dataframe for preictions from OG models
      preds = data.frame(matrix(NA, nrow = N, ncol = M))
      
      for (m in 1:M){
        model = xgboost(data = devs$data, label = devs$label,
                        nrounds = grid[m,1],
                        max.depth = grid[m,2],
                        eta = grid[m,3],
                        colsample_bytree = grid[m,4],
                        objective = "binary:logistic",
                        verbose = 0)
        
        preds[ ,m] = predict(model, vals$data)
        if (m == 24){
          paste("Preds df is done at",format(Sys.time(), "%H:%M:%S"))  %>% print()
        }
      }
      
      preds$ID = rownames(valsample5k)#no need to update-validation needs to be on5000
      
      # m long list with n*B dataframe for BS predictions
      bs = list()
      
      for (m in 1:M) {
        bs[[m]] = data.frame(matrix(NA, nrow = N, ncol = B ))
      }
      
      
      for (b in 1:B){
        b_ind = sample(1:nrow(devs$data), nrow(devs$data), replace = TRUE)
        bootdata = devs$data[b_ind, ]
        bootlabel = devs$label[b_ind]
        for (m in 1:M){
          bootmod =  xgboost(data = bootdata, label = bootlabel,
                             nrounds = xgb_grid[m,1],
                             max.depth = xgb_grid[m,2],
                             eta = xgb_grid[m,3],
                             colsample_bytree = xgb_grid[m,4],
                             objective = "binary:logistic",
                             verbose = 0)
          bs[[m]][ ,b] = predict(bootmod, vals$data )
          
        }
        
        #keeping track of progress
        if(b == 50){
          paste("First 50 BS iterations are done at",
                format(Sys.time(), "%H:%M:%S")) %>% print()
        }
        
        if(b == 100){
          paste("First 100 BS iterations are done at",
                format(Sys.time(), "%H:%M:%S")) %>% print()
        }
        
        if(b == 150){
          paste("First 150 BS iterations are done at",
                format(Sys.time(), "%H:%M:%S"))  %>% print()
        }
        
        if(b == B){
          paste("All BS iterations are done at",
                format(Sys.time(), "%H:%M:%S"))  %>% print()
        }
        
      }
      
      #sampling variances for n patients across M models
      var = data.frame(matrix(NA, nrow = N, ncol = M))
      
      for(m in 1:M){
        var[ ,m] = apply(bs[[m]], MARGIN = 1, FUN = var)
      }
      
      #standard error for n patient across M models
      se = data.frame(matrix(NA, nrow = N, ncol = M))
      
      for(m in 1:M){
        se[ ,m] = apply(bs[[m]], MARGIN = 1, FUN = sd)
      }
      
      
      return(list(preds = preds, var = var, se = se))
    }



##first run  - dev200 ####

    
  xgb_1 = sampvar_xgb(B=200, devs = devsample2h_xgb,
                        vals = valsample5k_xgb, grid = xgb_grid)

    
  tau_xgb_1 = taudist(xgb_1)
    
  plots_xgb_1 = tauplots(xgb_1, tau_xgb_1, family = "XGB")


    
##second run  - dev1000 ####

  xgb_2 = sampvar_xgb(B=200, devsample_xgb, valsample5k_xgb, xgb_grid)

  tau_xgb_2 = taudist(xgb_2 )

  plots_xgb_2 = tauplots(xgb_2, tau_xgb_2, family = "XGB" )



##third run - dev5k ####

      xgb_3 = sampvar_xgb(B = 200, devsample5k_xgb, valsample5k_xgb, xgb_grid)


      tau_xgb_3 = taudist(xgb_3)
   
      plots_xgb_3 = tauplots(xgb_3, tau_xgb_3, N = 5000, M = 24)






#KNN ####

data4knn = function(sample){
  sample$AGE = zscore(sample$AGE)
  return(sample)
}

devsampleknn2h = data4knn(devsample2h)
devsampleknn = data4knn(devsample)
devsampleknn5k = data4knn(devsample5k)
valsampleknn5k = data4knn(valsample5k)



sampvar_knn = function(B, devs, vals){
  
  
  paste("Started running at",format(Sys.time(), "%H:%M:%S"))  %>% print()
  
  N = nrow(vals)
  M = 24
  k_vec = floor(seq(7, nrow(devs)/12, length.out = M ))
  
  #predictions from OG models
  
  preds = as.data.frame(matrix(data = NA, nrow = N, ncol = M))
  
  for(m in 1:M){
    mod = knn3(DAY30 ~ ., devs, k = k_vec[m], prob = TRUE)
    preds[ ,m] = predict(mod, vals)[ ,2]
    
    if (m == 24){
      paste("Preds df is done at",
            format(Sys.time(), "%H:%M:%S"))  %>% print()
    }
  }
  
  preds$ID = rownames(valsample5k)
  
  
  # m long list with n*B dataframe for BS predictions
  bs = list()
  
  for (m in 1:M) {
    bs[[m]] = data.frame(matrix(NA, nrow = N, ncol = B ))
  }
  
  
  for (b in 1:B){
    b_ind = sample(1:nrow(devs), nrow(devs), replace = TRUE)
    bootsample = devs[b_ind, ]
    
    for(m in 1:M){
      bootmod = knn3(DAY30 ~ ., bootsample, k = k_vec[m], prob = TRUE)
      bs[[m]][ ,b] = predict(bootmod, vals)[ ,2]
    }
    
    #keeping track of progress
    if(b == 50){
      paste("First 50 BS iterations are done at",
            format(Sys.time(), "%H:%M:%S"))  %>% print()
    }
    
    if(b == 100){
      paste("First 100 BS iterations are done at",
            format(Sys.time(), "%H:%M:%S"))  %>% print()
    }
    
    if(b == 150){
      paste("First 150 BS iterations are done at",
            format(Sys.time(), "%H:%M:%S"))  %>% print()
    }
    
    if(b == B){
      paste("All BS iterations are done at",
            format(Sys.time(), "%H:%M:%S"))  %>% print()
    }
  }
  
  #sampling variances for n patients across M models
  var = data.frame(matrix(NA, nrow = N, ncol = M))
  
  for(m in 1:M){
    var[ ,m] = apply(bs[[m]], MARGIN = 1, FUN = var)
  }
  
  #standard error for n patient across M models
  se = data.frame(matrix(NA, nrow = N, ncol = M))
  
  for(m in 1:M){
    se[ ,m] = apply(bs[[m]], MARGIN = 1, FUN = sd)
  }
  
  
  return(list(preds = preds, var = var, se = se, k_vec = k_vec))
  
}



##knn run 1 - dev2h ####

      knn_1 = sampvar_knn(B = 200, devsampleknn2h, valsampleknn5k)
      
      tau_knn_1 = taudist(knn_1)
      
      plots_knn_1 = tauplots(knn_1, tau_knn_1, family = "knn")
    

##knn run 2 -  dev1k ####

      knn_2 = sampvar_knn(B = 200, 
                          devsampleknn, valsampleknn5k)
      
      
      tau_knn_2 = taudist(knn_2)
      
      plots_knn_2 = tauplots(knn_2, tau_knn_2, family = "knn")
      

##knn run3 - dev5k #### 
      
      knn_3 = sampvar_knn(B = 200, 
                          devsampleknn5k, valsampleknn5k)

      
      tau_knn_3 =  taudist(knn_3)

      plots_knn_3 = tauplots(knn_3, tau_knn_3, family = "knn")
      
      
      
      
      
   #All methods together####
      
      
      #pooling all results
      pool = function(rf, lr, xgb, knn,  M ){
        preds = as.data.frame(cbind(rf$preds[ ,1:24], lr$preds[ ,1:24],
                                    xgb$preds[ ,1:24], knn$preds))
        colnames(preds)[1:M] = paste("X", 1:M, sep = "")
        var = as.data.frame(cbind(rf$var, lr$var, xgb$var, knn$var))
        se = as.data.frame(cbind(rf$se, lr$se, xgb$se, knn$se))
        return(list(preds = preds, var = var, se = se))
      }
      
      
      tauplotsall = function(selist, taulist, ndev){ 
        
        predsgg = selist$preds[selist$preds$ID %in%
                                 c("36857","17553","6388","15077", 
                                   "570","37660", "25409", "39861", 
                                   "20809", "1783", "1095", "31058"), ]
        
        predsgg = melt(predsgg, id.vars = "ID")
        
        predsgg$variable = substring(predsgg$variable, first = 2)
        
        for(i in 1:nrow(predsgg)){
          
          if(predsgg$variable[i] %in% 1:24){
            predsgg$variable[i] = "RF"
          }
          else if(predsgg$variable[i] %in% 25:48){
            predsgg$variable[i] = "LR"
          }
          else if (predsgg$variable[i] %in% 49:72){
            predsgg$variable[i] = "XGB"
          }
          
          else if(predsgg$variable[i] %in% 73:96){
            predsgg$variable[i] = "KNN"
          }
        }
        
        
        
        #checking the spread of preds
        gg1 = ggplot(data = predsgg) +
          geom_jitter(aes(x = reorder(ID, value, mean), y = value, color = variable), 
                      alpha = 0.5, width = 0.25)+
          ylab(paste("predictions from all models (n=", ndev ,")"))+
          xlab("patient ID")+
          labs(color = "Model \nclass")+
          scale_y_continuous(limits = c(0,1))
        
        
        #dist of tau and sigma
        
        distdf = as.data.frame(cbind(taulist$tau, taulist$sigma_pp))
        
        distdf$V1 = ifelse(distdf$V1 > 0, distdf$V1, NA)
        zerotau =  sum(is.na(distdf$V1))
        
        gg2 = ggplot(data = distdf)+
          geom_density(aes(x = V1, color = "tau"), linewidth = 1)+
          geom_density(aes(x = V2, color = "sigma" ), linewidth = 1)+
          scale_color_manual(values = c("tau" = "lightblue",
                                        "sigma" = "darkolivegreen3"))+
          xlab("")+
          ylab("density")+
          scale_x_continuous(limits = c(0,0.1))+
          scale_y_continuous(limits = c(0,122))+
          theme(legend.position = c(.95, .95),
            legend.justification = c("right", "top"),
            legend.box.just = "right")+
          annotate("label", x = 0.08, y = 80, 
                     label = paste('No. of zero tau values:', zerotau), 
                     color = "black", size = 4, 
                     label.padding = unit(0.5, "lines"))
         
        
        #pairs of tau and sigma
        pairsdf =  cbind.data.frame(taulist$tau, 
                                    taulist$sigma_pp, 
                                    rownames(valsample5k))[1:700,]
        colnames(pairsdf) = c("tau", "sigma", "ID")
        
        gg3 = ggplot(pairsdf)+
          geom_point(mapping = aes(x = ID,y = tau, color = "tau"))+
          geom_point(mapping = aes(x = reorder(ID, sigma), y = sigma, 
                                   color = "sigma"))+
          theme(axis.text.x=element_blank(), 
                axis.ticks.x=element_blank())+
          xlab("patients")+
          ylab("value of tau and sigma")+
          ylim(0,0.22)+
          scale_color_manual(values = c("tau" = "lightblue",
                                        "sigma" = "darkolivegreen3"))+
          theme(legend.position = c(.05, .95),
            legend.justification = c("left", "top"),
            legend.box.just = "left")
        
        return(list(preds = gg1, tausigma = gg2, tspairs = gg3))
      }
      
      

      
      # dev2h
      
      
      all_1 = pool(rf_1, lr_1, xgb_1, knn_1, M = 96)
      
      tau_all_1 = taudist(all_1, M = 96)
      
      plots_all_1 = tauplotsall(all_1, tau_all_1, ndev = 200)
      
      plots_all_1$tausigma
      
      
      
      #dev 1000
      
      all_2 = pool(rf_2, lr_2, xgb_2, knn_2, M = 96)
      
      tau_all_2 = taudist(all_2, M = 96)
      
      plots_all_2 = tauplotsall(all_2, tau_all_2, ndev = 1000)
      
      plots_all_2$tausigma
      ggsave("tausigma2.png")
      
      #dev5k
      
      all_3 = pool(rf_3, lr_3, xgb_3, knn_3, M = 96)
      
      tau_all_3 = taudist(all_3, M = 96)
      
      plots_all_3 = tauplotsall(all_3, tau_all_3, ndev = 5000)
      
      plots_all_3$tausigma
      ggsave("tausigma3.png")
      
      summary2h = data.frame("LR" = c(tau_lr_1$meantau, tau_lr_1$sigma_avgm),
                             "RF" = c(tau_rf_1$meantau, tau_rf_1$sigma_avgm),
                             "XGB" = c(tau_xgb_1$meantau, tau_xgb_1$sigma_avgm),
                             "KNN" = c(tau_knn_1$meantau, tau_knn_1$sigma_avgm),
                             "All" = c(tau_all_1$meantau, tau_all_1$sigma_avgm))
      
      
      summary1k = data.frame("LR" = c(tau_lr_2$meantau, tau_lr_2$sigma_avgm),
                             "RF" = c(tau_rf_2$meantau, tau_rf_2$sigma_avgm),
                             "XGB" = c(tau_xgb_2$meantau, tau_xgb_2$sigma_avgm),
                             "KNN" = c(tau_knn_2$meantau, tau_knn_2$sigma_avgm),
                             "All" = c(tau_all_2$meantau, tau_all_2$sigma_avgm))
      
      summary5k = data.frame("LR" = c(tau_lr_3$meantau, tau_lr_3$sigma_avgm),
                             "RF" = c(tau_rf_3$meantau, tau_rf_3$sigma_avgm),
                             "XGB" = c(tau_xgb_3$meantau, tau_xgb_3$sigma_avgm),
                             "KNN" = c(tau_knn_3$meantau, tau_knn_3$sigma_avgm),
                             "All" = c(tau_all_3$meantau, tau_all_3$sigma_avgm))
      
      
      tables <- list(summary2h = summary2h, summary1k = summary1k, summary5k = summary5k)
      
      
      for (i in seq_along(tables)) {
        tables[[i]] <- as.data.frame(lapply(tables[[i]], function(x) round(x, 5)))
        rownames(tables[[i]]) <- c("tau", "sigma")
      }
  
      
      #checking if total uncertainty decreases with sample size (tau^2 + sigma^2)
      
      totalvar_2h = sqrt(mean(tau_all_1$tau^2) + mean(tau_all_1$sigma_pp^2))
      totalvar_1k = sqrt(mean(tau_all_2$tau^2) + mean(tau_all_2$sigma_pp^2))
      totalvar_5k = sqrt(mean(tau_all_3$tau^2) + mean(tau_all_3$sigma_pp^2))
      
      
      
      

#Super Learner####
      
#creating custom learners

      
learners_glmnet = create.Learner("SL.glmnet", detailed_names = TRUE,
                                       tune = list(alpha = c(0,0.5,1)))
      
learners_xgb = create.Learner("SL.xgboost", detailed_names = TRUE,
                              tune = list(max.depth = c(3,5),
                                          eta = c(0.1,0.01),
                                          colsample_bytree = c(0.75, 1) ))
      
learners_rf = create.Learner("SL.ranger", detailed_names = TRUE,
                             params = list(probability = TRUE),
                             tune = list(mtry = c(2,4),
                                         minnodesize = c(1,3,10)))


SL_liby = c(learners_glmnet$names, learners_rf$names,
           learners_xgb$names)

  
SL_1k_x = SuperLearner(Y = as.numeric(devsample$DAY30)-1, X = devsample[ ,2:9],
                           family = binomial(),
                            SL.library = SL_liby,
                            method = "method.CC_nloglik",
                          cvControl = list(V = 10, stratifyCV = TRUE, shuffle = TRUE)
)

preds_sl1_x = predict(SL_1k_x, newdata = valsample5k[ ,2:9], onlySL = TRUE)$pred

#comparing estimated predictions by RE model and SL ensemble
    

    ggplot(data = data.frame(RE = tau_all_2$muhat, SL = preds_sl1_x))+
      geom_point(aes(x = RE, y = SL), alpha = 0.4)+
      scale_x_continuous( limits = c(0,0.6))+
      scale_y_continuous( limits = c(0,0.75))+
      geom_abline(slope = 1)+
      xlab("Predictions from random effects model")+
      ylab("Predictions from super learner ensemble")
    