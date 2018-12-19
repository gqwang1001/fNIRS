library(stringr)
library(zoo)
library(plotly)
library(ggplot2)
library(reshape2)


rm(list = ls())
source("nbmfpca.R")

########## read data
    dir_vft = "../../fNIRS/processed12_18_2015/VFT/Hb VFT all_Data"
    files = list.files(dir_vft, pattern = "Oxy.csv$")
    subjects = as.numeric(str_extract(files, "[0-9]+"))-10000
    nchannels = 52
    nsubj = length(subjects)
    t_course = 1:1251
    for(i in 1:nsubj){
      if(i == 1){
        temp = read.csv(file = paste0(dir_vft,"/",files[i]),skip = 40)
        dat.vft = t(as.matrix(temp[t_course,-c(1,54:58)]))
      }else{
        temp = read.csv(file = paste0(dir_vft,"/",files[i]),skip = 40)
        dat.vft = rbind(vft_dat, t(as.matrix(temp[t_course,-c(1,54:58)])))
      }
    }
    saveRDS(dat.vft, 'vft_data.rds')

########### detect noisy channels
    wid = 50 #window width
    
    roll.var = rowMeans(t(rollapply(t(dat.vft),width = wid,sd)))
    
    nsubj = 53; nchannels = 52
    hist(roll.var,breaks = seq(0,3,by =1e-4))
    threshold =0.05
    length(which(roll.var>threshold))
    roll.var.mat = fnirs.reshape(roll.var,nsubj,nchannels = nchannels)
    ind_large_vft = which(abs(roll.var.mat)>threshold,arr.ind = T)
    saveRDS(ind_large_vft,file = "rej_ind_vft.rds")
    
    
    # visulize noisy channels
    jet.colors <-
      colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                         "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    
    data_plot = as.data.frame(roll.var.mat)
    rownames(roll.var.mat) = 1:nsubj
    colnames(roll.var.mat) = 1:nchannels
    
    roll.temp = melt(roll.var.mat)
    roll.temp$Subjects = roll.temp$Var1
    roll.temp$Channels = roll.temp$Var2
    roll.temp$RunningSD= roll.temp$value
    ggplot(data = roll.temp) + geom_tile(aes(x=Subjects, y=Channels, fill=RunningSD))+
      scale_fill_gradientn(colours = jet.colors(10))+theme_bw()+scale_x_continuous(breaks = seq(0,50,by=5))


############## Run unbalaced MFPCA
    
    vft.fit = nbmfpca(rejection.index = ind_large_vft, dat = dat.vft, nsubj = 53, nchannels = 52)
    tn = 1251
    ## plot overall mean
    par(mfrow = c(1,1))
    plot((1:tn)/10,vft.fit$mu,type = "l", xlab = "time(sec)", ylab = '',main = "VFT Overall Mean")
    ## plot channel-specific mean
    
    ngrid= 3
    par(mfrow=c(ngrid,ngrid))
    for(k in 1:(ngrid^2))
      plot(vft.fit$eta[k,],type = "l", ylab = "", xlab = "time(sec)", main = paste("eta", k))
    
    par(mfrow = c(2,3))
    percent1 = (vft.fit$lambda1)/sum(vft.fit$lambda1)
    percent2 = (vft.fit$lambda2)/sum(vft.fit$lambda2)
    
    for(i in 1:3){
      plot((1:tn)/10,vft.fit$phi1[,i], main = paste0("L1 eigenfunc",i, " ",round(percent1[i]*100,2),"%"),ylab="",xlab="time (sec)",type = 'l')
    }
    
    for(i in 1:3){
      plot((1:tn)/10,vft.fit$phi2[,i], main = paste0("L2 eigenfunc",i," ",round(percent2[i]*100,2),"%"),ylab="",xlab="time (sec)",type = 'l')
    }   
    
    
    library(xtable)
    klevel = vft.fit$K1
    results.table = data.frame(phi1 = vft.fit$lambda1[1:klevel], phi2=vft.fit$lambda2[1:klevel])
    results.table10000 = data.frame(phi1 = vft.fit$lambda1[1:klevel]*1e4, phi2=vft.fit$lambda2[1:klevel]*1e4)
    xtable(results.table10000)
    
    totalvaribility = sum(vft.fit$lambda1+vft.fit$lambda2)
    cumu.table = results.table/totalvaribility
    cumu.table$cumulative = cumsum(cumu.table$phi1+cumu.table$phi2)
    cumu.table[klevel+1,1:2] = apply(cumu.table[,1:2],2,sum)
    cumu.table[klevel+1,3] = max(cumu.table[1:klevel,3])
    xtable(cumu.table)
    
    
    par(mfrow = c(1,1), cex=0.7)
    plot(1:tn/10,dat.result$mu,type = "l",ylim = c(-0.025,0.025),xlab = "time (sec)", ylab = "", main="VFT mean (oxy-Hb)")
    for(i in 1:5)  lines(1:tn/10, vft.fit$eta[,i]+vft.fit$mu,col = i+1)
    legend("bottomleft", legend = c("ovealll", paste("channel",1:5)), 
           col = 1:6, pch = 20)


###### plot PC scores against clinic information

    comb_data = readRDS('demoInfo.rds')
    
    comb_data$cam = as.factor(comb_data$cam)
    comb_data$gender = as.numeric(comb_data$gender)
    s1.level1 = vft.fit$scores1[,1]
    s1.level2 = vft.fit$scores1[,2]
    
    library(ggplot2)
    library(ggpubr)
    s1_age = ggplot(data = comb_data,aes(age,s1.level1))+geom_point()+geom_smooth()+labs(title = 'score 1 vs. age')+theme(axis.title.x=element_blank(),axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5))
    s1_sat = ggplot(data = comb_data,aes(sat,s1.level1))+geom_point()+geom_smooth()+labs(title = 'score 1 vs. SAT')+theme(axis.title.x=element_blank(),axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5))
    s1_cam = ggplot(data = comb_data,aes(factor(cam),s1.level1))+geom_boxplot()+labs(title = 'score 1 vs. CAM')+theme(axis.title.x=element_blank(),axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5))
    cam_gender = ggplot(data = comb_data,aes(x=factor(cam),fill=factor(gender)))+geom_bar()+labs(title = 'CAM vs. age')+theme(axis.title.x=element_blank(),axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5))
    s1_yoe = ggplot(data = comb_data,aes(factor(YoE),s1.level1))+geom_boxplot() +labs(title = 'score 1 vs. age')+theme(axis.title.x=element_blank(),axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5))
    s1_gender = ggplot(data = comb_data,aes(factor(gender),s1.level1))+geom_boxplot()+labs(title = 'score 1 vs. gender')+theme(axis.title.x=element_blank(),axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5))
    
    s1_DSR = ggplot(data = comb_data,aes(dsr,s1.level1))+geom_point()+geom_smooth()+labs(title = 'score 1 vs. DSR')+theme(axis.title.x=element_blank(),axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5))
    
    s1_VFT = ggplot(data = comb_data,aes(vft,s1.level1))+geom_point()+geom_smooth()+labs(title = 'score 1 vs. VFT')+theme(axis.title.x=element_blank(),axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5))
    
    s1_s2 = ggplot(data = comb_data,aes(s1.level2,s1.level1))+geom_point()+labs(title = 'score 1 vs. score 2')+theme(axis.title.x=element_blank(),axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5))
    
    s1_meld = ggplot(data = comb_data,aes(meld,s1.level1))+geom_point()+geom_smooth()+labs(title = 'score 1 vs. MELD')+theme(axis.title.x=element_blank(),axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5))
    
    s1_cci = ggplot(data = comb_data,aes(cci,s1.level1))+geom_point()+geom_smooth()+labs(title = 'score 1 vs. CCI')+theme(axis.title.x=element_blank(),axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5))
    s1_inr = ggplot(data = comb_data,aes(inr,s1.level1))+geom_point()+geom_smooth()+labs(title = 'score 1 vs. INR')+theme(axis.title.x=element_blank(),axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5))
    
    ggarrange(s1_s2, s1_cam, s1_gender, s1_age, s1_sat, s1_VFT, s1_meld, s1_cci, s1_inr,ncol = 3, nrow = 3)