#main analysis of the birth weight ##
# singleton  data 9477 sample ##
library(dplyr)
library(stringi)
library(stringr)
library(ggpubr)
library(tidyr)
library(lmerTest)
library(sjstats)
library(lsr)
library(mice)
library(Matching)
library(MatchIt)
library(mgcv)
library(simr)
library(R.matlab)
##  ##
###load the site info #
all_sitealsub<-readMat('E:/twin/6.9/site.mat')
all_sitealsub<-unlist(all_sitealsub)%>%as.matrix()
all_sitealsub<-apply(all_sitealsub, 1, function(x){
  x<-str_split(x,'e')%>%unlist()
  x<-as.numeric(x[2])
})
all_sitealsub<-as.matrix(all_sitealsub)
cov_id<-readMat('E:/twin/6.9/cov_id.mat')%>%unlist()%>%as.matrix()
rownames(all_sitealsub)<-cov_id
ininfo<-read.csv('./5.0/dhx01.csv',fill = TRUE)
ininfo<-ininfo[-c(1,2),]
untwin_id_data<-read.table('E:/twin/6.9/un_twin_id.txt')
single_data<-untwin_id_data$x
birth_weight<-ininfo[which(!is.na(ininfo$birth_weight_lbs)),]
birth_weight<-birth_weight[which(birth_weight$eventname=='baseline_year_1_arm_1'),]
birth_weight<-birth_weight[which(!duplicated(birth_weight$src_subject_id)),]
rownames(birth_weight)<-birth_weight$src_subject_id
single_weight<-birth_weight[single_data,]
###
lbs_single<-single_weight$birth_weight_lbs%>%as.numeric()
orz_single<-single_weight$birth_weight_oz%>%as.numeric()
orz_single[which(is.na(orz_single))]<-0
single_weight_oz<-lbs_single*16+orz_single
## class the birth weight by  different gestational week #
single_gestation<-ges_allsample[single_data,]$ges_week_0_all%>%as.numeric()
#gender info#
single_gender<-twins_info[single_data,]$sex%>%as.numeric()## 1 is female 0 is male 
#####
combined_data<-data.frame(id=single_data,weright=single_weight_oz,gestational_weeek=single_gestation,gender=single_gender)
##
combined_real<-na.omit(combined_data)%>%as.data.frame() ##8344 sample for singleton # 
combined_real<-combined_real
combined_real$weright<-combined_real$weright*0.02834952
###
intersect(cbcl_socre_base$subjectkey,combined_real$id)%>%length()
## original data ##
get_original_birth_d<-function(combined_real){
  ## transform the OZ to KG #
  rownames(combined_real)<-combined_real$id
  ######  ###
  low_birth<-combined_real[which(combined_real$weright<2.5),]%>%rownames()
  normal_birth<-combined_real[which(combined_real$weright>=2.5&combined_real$weright<4),]%>%rownames()
  over_birth_weight<-combined_real[which(combined_real$weright>4),]%>%rownames()
  ##### return the id
  re_list<-list()
  re_list[["LBW"]]<-low_birth
  re_list[["Normal"]]<-normal_birth
  re_list[["overweight"]]<-over_birth_weight
  return(re_list)
}
LBW_data<- group_info[["LBW"]]
Normal_data<-group_info[["Normal"]]
overweight_data<-group_info[["overweight"]]
# weight with gestation week #

LBW_ges<- group_info_gestation_wk[["LBW"]]
Normal_ges<-group_info_gestation_wk[["Normal"]]
overweight_ges<-group_info_gestation_wk[["overweight"]]
###
LBW_se_inter<-intersect(LBW_ges,LBW_data)
nor_se_inter<-intersect(Normal_ges,Normal_data)
over_se_inter<-intersect(overweight_data,overweight_ges)
###
intersect(LBW_ges,lbw)%>%length()
intersect(overweight_ges,over_data)%>%length()



###
### diabetes and hypertension ##
diabetes_data<-twins_info[which(twins_info$devhx_10i3_p!='999'&twins_info$devhx_10i3_p!=''),]
hypertension_data<-twins_info[which(twins_info$devhx_10j3_p!='999'&twins_info$devhx_10j3_p!=''),]
rownames(race_info)<-race_info$subjectkey
##
income_data_com$demo_comb_income_v2_l<-as.numeric(income_data_com$demo_comb_income_v2_l)
income_ran1<-income_data_com[which(income_data_com$demo_comb_income_v2_l<=5),]$subjectkey
income_ran2<-income_data_com[which(income_data_com$demo_comb_income_v2_l==6),]$subjectkey
income_ran3<-income_data_com[which(income_data_com$demo_comb_income_v2_l==7),]$subjectkey
income_ran4<-income_data_com[which(income_data_com$demo_comb_income_v2_l==8),]$subjectkey
income_ran5<-income_data_com[which(income_data_com$demo_comb_income_v2_l>=9),]$subjectkey
unknow<-setdiff(union(untwins_data,twins_data),income_data_com$subjectkey)
####
income_d<-c(income_ran1,income_ran2,income_ran3,income_ran4,income_ran5,unknow)
id_data<-c(rep(1,length(income_ran1)),rep(2,length(income_ran2)),
           rep(3,length(income_ran3)),rep(4,length(income_ran4)),rep(5,length(income_ran5)),
           rep(6,length(unknow)))
data_d<-data.frame(income_d,id_data)
rownames(data_d)<-data_d$income_d
### Regression model#
compute_T_value<-function(data,normal_data,abnormal_data,cog_cbcl,fill_pos,PSM_M){
  # normal_data
  # abnormal_data: LGA or SGA data 
  # cog_cbcl: string for "cbcl" or "cog" 
  # fill_pos: Bool variant (True for filling the NA data)
  # PSM_M: Bool variant(True for performing the PSM match)
  library(compute.es)

  single_data<-normal_data
  single_data<-intersect(single_data,rownames(all_sitealsub))
  #abnormal_data<-low_weight_id
  twins_data<-abnormal_data
  ###
  #sinle_data
  income_data_com_single<-income_data_com[single_data,]$demo_comb_income_v2_l%>%as.numeric()
  age_info_father_single<-age_info_father[single_data,]$devhx_4_p%>%as.numeric()
  age_info_mother_single<-age_info_mother[single_data,]$devhx_3_p%>%as.numeric()
  mother_smoke_single<-mother_smoke[single_data,]$devhx_9_tobacco%>%as.numeric()
  mother_drink_single<-mother_alcohol[single_data,]$devhx_9_alcohol%>%as.numeric()
  mother_marijuana_sinlge<-mother_marijuana[single_data,]$devhx_9_marijuana%>%as.numeric()
  mother_morphine_single<-mother_morphine[single_data,]$devhx_9_her_morph%>%as.numeric()
  diabetes_data_single<-diabetes_data[single_data,]$devhx_10i3_p%>%as.numeric()
  hypertension_single<-hypertension_data[single_data,]$devhx_10j3_p%>%as.numeric()
  race_single<-race_encoder[single_data,]
  ###
  ges_single<-ges_allsample[single_data,]$ges_week_0_all%>%as.numeric()
  BMI_single<-physics_BMI[single_data,]$BMI%>%as.numeric()
  Puberty_single<-all_puber_info[single_data,]$all_puber%>%as.numeric()
  educa_single<-ecucation_par_info[single_data,]$demo_prnt_ed_v2_l%>%as.numeric()
  sex_single<-twins_info[single_data,]$sex%>%as.numeric()
  age_single<-twins_info[single_data,]$interview_age%>%as.numeric()
  single_site<-all_sitealsub[single_data,]
  ##twins data ##
  income_data_com_twins<-income_data_com[twins_data,]$demo_comb_income_v2_l%>%as.numeric()
  age_info_father_twins<-age_info_father[twins_data,]$devhx_4_p%>%as.numeric()
  age_info_mother_twins<-age_info_mother[twins_data,]$devhx_3_p%>%as.numeric()
  mother_smoke_twins<-mother_smoke[twins_data,]$devhx_9_tobacco%>%as.numeric()
  mother_drink_twins<-mother_alcohol[twins_data,]$devhx_9_alcohol%>%as.numeric()
  mother_marijuana_twins<-mother_marijuana[twins_data,]$devhx_9_marijuana%>%as.numeric()
  mother_morphine_twins<-mother_morphine[twins_data,]$devhx_9_her_morph%>%as.numeric()
  ges_twins<-ges_allsample[twins_data,]$ges_week_0_all%>%as.numeric()
  BMI_twins<-physics_BMI[twins_data,]$BMI%>%as.numeric()
  Puberty_twins<-all_puber_info[twins_data,]$all_puber%>%as.numeric()
  educa_twins<-ecucation_par_info[twins_data,]$demo_prnt_ed_v2_l%>%as.numeric()
  sex_twins<-twins_info[twins_data,]$sex%>%as.numeric()
  age_twins<-twins_info[twins_data,]$interview_age%>%as.numeric()
  twin_site<-all_sitealsub[twins_data,]
  race_twin<-race_encoder[twins_data,]
  diabetes_data_twin<-diabetes_data[twins_data,]$devhx_10i3_p%>%as.numeric()
  hypertension_twin<-hypertension_data[twins_data,]$devhx_10j3_p%>%as.numeric()
  #####
  income_var<-c(income_data_com_single,income_data_com_twins)
  age_var_father<-c(age_info_father_single,age_info_father_twins)
  age_var_mother<-c(age_info_mother_single,age_info_mother_twins)
  mother_var_smo<-c(mother_smoke_single,mother_smoke_twins)
  mother_var_alcho<-c(mother_drink_single,mother_drink_twins)
  mother_marijuana_var<-c(mother_marijuana_sinlge,mother_marijuana_twins)
  mother_morphine_var<-c(mother_morphine_single,mother_morphine_twins)
  BMI_var<-c(BMI_single,BMI_twins)
  ges_var<-c(ges_single,ges_twins)
  Puberty_var<-c(Puberty_single,Puberty_twins)
  educa_var<-c(educa_single,educa_twins)
  sex_var<-c(sex_single,sex_twins)
  age_var<-c(age_single,age_twins)
  site_info<-c(single_site,twin_site)
  race_info<-c(race_single,race_twin)
  diabetes_info<-c(diabetes_data_single,diabetes_data_twin)
  hypertension_info<-c(hypertension_single,hypertension_twin)
  group_diff<-c(rep(0,length(single_data)),rep(1,length(twins_data)))
  ####
  #print(length(twins_data))
  name_cbcl<-colnames(data)
  sig_data<-vector()
  p_values_sig<-vector()
  sample_num<-vector()
  single_mean_d<-vector()
  twins_mean_d<-vector()
  t_value<-vector()
  cohen_d<-vector()
  power_ci_out<-vector()
  power_t_test<-vector()
  ## 11:89
  save_list<-list()
  ##18:106
  
  if(cog_cbcl=="CBCL"){
    start_n=10
    seq_dep=4
  }else if(cog_cbcl=="cog"){
    start_n=10 
    seq_dep=1
  }else if(cog_cbcl=="struct") {
    start_n=1
    seq_dep=1
  }
  co_variants<-data.frame(group_diff,income_var,sex_var,age_var,site_info,diabetes_info,
                          age_var_father,BMI_var,ges_var,Puberty_var,educa_var,race_info,hypertension_info,
                          age_var_mother,mother_morphine_var,mother_var_smo,mother_var_alcho)
  #### fill the NA value in the datafram##
  if(fill_pos==TRUE){
    fill_data<-mice(co_variants,m=3,method = "pmm",maxit =1,seed=1)
    completedata<-complete(fill_data)
    co_variants<-completedata
  }
  #
  for (i in seq(from=start_n,to=ncol(data),seq_dep)) {
    #i=14
    cbcl_socre_single<-data[single_data,i]%>%as.numeric()
    cbcl_socre_twins<-data[twins_data,i]%>%as.numeric()
    ###
    if(length(cbcl_socre_single[which(is.na(cbcl_socre_single))])<(2*(length(cbcl_socre_single))/3))
    {
      
      cbcl_var_cat<-c(cbcl_socre_single,cbcl_socre_twins)
      glm_dataf<-cbind(cbcl_var_cat,co_variants)
      glm_dataf<-na.omit(glm_dataf)
      ##
      
      if(PSM_M==TRUE){
        ## Uisng PSM achieve the sample match ##
        ps_model <- matchit(group_diff~income_var+hypertension_info+sex_var+ges_var+race_info+mother_morphine_var+
                              age_var+age_var_father+BMI_var+Puberty_var+diabetes_info+mother_var_smo+
                              educa_var,distance = 'glm',
                            ratio=2, data = glm_dataf, method = "nearest")
        matched_data <- match.data(ps_model)
        glm_dataf<-matched_data
      }
      
      ##
      cbcl_socre_single<-glm_dataf[which(glm_dataf$group_diff==0),1]
      cbcl_socre_twins<-glm_dataf[which(glm_dataf$group_diff==1),1]
      ##the mean and the std
      sinlle_mean<-mean(cbcl_socre_single[which(!is.na(cbcl_socre_single))])%>%round(digits = 3)
      sinlle_std<-sd(cbcl_socre_single[which(!is.na(cbcl_socre_single))])%>%round(digits = 3)
      sinlle_mean_co<-paste(as.character(sinlle_mean),"(",as.character(sinlle_std),")",sep = "")
      ##
      twins_mean<-mean(cbcl_socre_twins[which(!is.na(cbcl_socre_twins))])%>%round(digits = 3)
      twins_std<-sd(cbcl_socre_twins[which(!is.na(cbcl_socre_twins))])%>%round(digits = 3)
      twins_mean_co<-paste(as.character(twins_mean),"(",as.character(twins_std),")",sep = "")
      
      ### get the glm datafram ges_var 
      twins_data_us<-glm_dataf[which(glm_dataf$group_diff==1),]
      single_data_us<-glm_dataf[which(glm_dataf$group_diff==0),]
      
      lin.mod = lmer(cbcl_var_cat~group_diff+income_var+hypertension_info+sex_var+ges_var+
                       age_var+age_var_father+BMI_var+Puberty_var+diabetes_info+
                       educa_var+age_var_mother+mother_var_alcho+mother_var_smo+(1|site_info),data=glm_dataf)
      
      # print(kappa(glm_dataf[,c(1,2,3,4,5,6,7,8,12,9,11,14,17,15,18)],exact= TRUE))
      pVal<-summary(lin.mod)$coefficients[2,5]
      print(pVal)
      #%>%round(digits = 3)
      tVal<-summary(lin.mod)$coefficients[2,4]%>%round(digits = 3)
      
      if(nrow(twins_data_us)>3){
        # print(name_cbcl[i])
        cohen_D<-cohensD(cbcl_socre_single[which(!is.na(cbcl_socre_single))],cbcl_socre_twins[which(!is.na(cbcl_socre_twins))])%>%round(digits = 3)
        sample_num<-c(sample_num,nrow(glm_dataf))
        single_mean_d<-c(single_mean_d,sinlle_mean_co)
        twins_mean_d<-c(twins_mean_d,twins_mean_co)
        sig_data<-c(sig_data,i)
        # print(pVal)
        p_values_sig<-c(p_values_sig,pVal)%>%round(digits = 3)
        t_value<-c(t_value,tVal)%>%round(digits = 3)
        cohen_d<-c(cohen_d,cohen_D)
        ##power analysis 
        test_data<-powerSim(fit = lin.mod, test = fixed('group_diff', method = 't'), nsim=50)
        power_values_t_test<-summary(test_data)$mean
        print(test_data)
        power_t_test<-c(power_t_test,power_values_t_test)
        fixed_effect <- fixef(lin.mod)  # 
        random_effect <- VarCorr(lin.mod)  #
        
        sim_model <- makeLmer(cbcl_var_cat~group_diff+income_var+hypertension_info+sex_var+ges_var+
                                age_var+age_var_father+BMI_var+Puberty_var+diabetes_info+
                                educa_var+age_var_mother+mother_var_alcho+mother_var_smo+(1|site_info), fixef=fixed_effect, VarCorr=random_effect, sigma=2, data=glm_dataf)
        
        
        # 设置样本量
        
        # 进行功效分析
        power_result <- powerSim(sim_model, nsim=50,test = fcompare(cbcl_var_cat~group_diff))
        power_values<-summary(power_result)$mean
        ci_lower<-summary(power_result)$lower 
        ci_upper<-summary(power_result)$upper ##
        power_ci_all<-paste(as.character(power_values),'(',as.character(round(ci_lower,2)),',',as.character(as.character(round(ci_upper,2))),')',
                            sep = "")
        print(power_ci_all)
        power_ci_out<-c(power_ci_out,power_ci_all)
      }
    }
  }
  ## In this way, the cohen'd was compute by the mean and sd #
  redata<-data.frame(name_cbcl[sig_data],twins_mean_d,single_mean_d,p_values_sig,t_value,cohen_d,power_ci_out,power_t_test)
  colnames(redata)<-c('Terms','Abormal_mean','Normal_mean','P','T','cohen_d','power_ci','power_t_test')
  plot_mean_sd<-function(redata1){
    redata_plot<-redata1[,1:3]
    colnames(redata_plot)<-c('Item','Abormal','Normal')
    redata_plot<-gather(redata_plot,category,value = 'iterm_value',2:3)
    redata_plot$SE_value<-sapply(redata_plot$iterm_value,function(x){
      x<-str_split(x,'[(]')%>%unlist()
      x<-str_split(x[2],'[)]')%>%unlist()
      x<-x[1]%>%as.numeric()
    })
    redata_plot$iterm_value<-sapply(redata_plot$iterm_value,function(x){
      x<-str_split(x,'[(]')%>%unlist()
      x<-x[1]%>%as.numeric()
    })
    ###
    p<-ggbarplot(redata_plot,x= 'Item',y= 'iterm_value',
                 color = "category",fill = "category", palette = "npg",position = position_dodge(0.8))+
      theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1,size=8))
    #+geom_errorbar(aes(ymax=iterm_value+SE_value),width=.2, position=position_dodge(.6))
    return(p)
  }
  if(cog_cbcl=='cog'){
    rownames(redata)<-redata$Terms
    redata<- redata[cog_names,]
    p1<-plot_mean_sd(redata)
    print(p1)
  }else if(cog_cbcl=='CBCL'){
    redata1<-redata[which(redata$Terms!='cbcl_scr_syn_totprob_r'),]
    ote_data<-redata[which(redata$Terms=='cbcl_scr_syn_totprob_r'),1:3]
    p1<-plot_mean_sd(redata1)
    p2<-plot_mean_sd(ote_data)
    c<-ggarrange(p2,p1,labels = c("A", "B"),ncol = 2, nrow = 1)
    print(c)
  }else{
    p1<-plot_mean_sd(redata)
    print(p1)
  }
  ## plot the datafram ##
  p_adjust_data<-p_values_sig<-p.adjust(redata$P,method = "BH")%>%round(digits = 3)
  redata<-cbind(redata,p_adjust_data)
  print('___________________________________________')
  print(nrow(twins_data_us))
  print(nrow(single_data_us))
  return(redata)
}
##
cbcl_low_nor<-compute_T_value(cbcl_socre_base,Normal_data,LBW_data,'CBCL',FALSE,FALSE)
#
cbcl_over_nor<-compute_T_value(cbcl_socre_base,Normal_data,overweight_data,'CBCL',FALSE,FALSE)
##
cog_low_AGA<-compute_T_value(cognitive_baseline,Normal_data,LBW_data,'cog',FALSE,FALSE)

cog_over_AGA<-compute_T_value(cognitive_baseline,Normal_data,overweight_data,'cog',FALSE,FALSE)
###
write.csv(cbcl_low_nor,'cbcl_low_power.csv',row.names = FALSE,quote = FALSE)
write.csv(cbcl_over_nor,'cbcl_large_power.csv',row.names = FALSE,quote = FALSE)
write.csv(cog_low_AGA,"LBW_Normal_cog_power.csv",row.names = FALSE,quote = FALSE)
write.csv(cog_over_AGA,"over_Nromal_cog_power.csv",row.names = FALSE,quote = FALSE)
###graph_data test #
graph2ppt(file='5.29_end.pptx',width=10,height=8,append=TRUE)
#####

###
###
SGA_Thicknesse<-compute_T_value(thickness_all[,149:151],Normal_data,LBW_data,'struct',FALSE,FALSE)
SGA_volume<-compute_T_value(volume_all[,149:151],Normal_data,LBW_data,'struct',FALSE,FALSE)
SGA_area<-compute_T_value(area_all[,149:151],Normal_data,LBW_data,'struct',FALSE,FALSE)
###
LGA_Thicknesse<-compute_T_value(thickness_all[,149:151],Normal_data,overweight_data,'struct',FALSE,FALSE)
LGA_volume<-compute_T_value(volume_all[,149:151],Normal_data,overweight_data,'struct',FALSE,FALSE)      
LGA_area<-compute_T_value(area_all[,149:151],Normal_data,overweight_data,'struct',FALSE,FALSE)
### save——csv——files
all_files_combined<-rbind(SGA_Thicknesse,SGA_area,SGA_volume,LGA_Thicknesse,LGA_area,LGA_volume)
write.csv(all_files_combined,'all_files_struc_power.csv',row.names = FALSE,quote = FALSE)

LBW_sig_volume<-LBW_volume[which(LBW_volume$p_adjust_data<=0.05),]
LBW_sig_area<-LBW_area[which(LBW_area$p_adjust_data<=0.05),]
LBW_sig_thickness<-LBW_Thicknesse[which(LBW_Thicknesse$p_adjust_data<=0.05),]
###
write.csv(LBW_Thicknesse,"MBW_thickness_all.csv",row.names = FALSE)
write.csv(LBW_volume,"MBW_volume_all.csv",row.names = FALSE)
write.csv(LBW_area,"MBW_area_all.csv",row.names = FALSE)

