#Get the score of the genes kept after the binarization

stat_score=apply(genoP,2,function(x){
  ns=names(which(x==1))
  mean(prop_net[ns,1])
})

keep=!is.na(stat_score)
prop_net=prop_net[,keep]
stat_score=stat_score[keep]

for(k in 1:ncol(prop_net)){
  zscale=scale(c(mean(stat_score[k]),prop_net[,k]))
  if(k==1){
    stat_score_df=c(zscale[1],quantile(zscale,probs=c(0.7,0.8,0.9,0.95,0.99)))
  }else{
    stat_score_df=rbind(stat_score_df,c(zscale[1],quantile(zscale,probs=c(0.7,0.8,0.9,0.95,0.99))))
  }
}

rownames(stat_score_df)=colnames(prop_net)
colnames(stat_score_df)[1]="average_score"
write.csv(stat_score_df,file="stat_score_df.csv",quote = F,row.names = T)
