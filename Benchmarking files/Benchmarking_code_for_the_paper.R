
#Dependencies for benchmarking (apart from those needed for RescueSC)
library(cellhashR)
library(DropletUtils)
library(deMULTIplex)
library(data.table)

#Creating some functions for the future; it is needed to run them first before performing benchmarking
CompareRescueTag<-function(dataset, other_tool_tags, tool_name){
  compared<-c()
  for(i in 1:length(dataset$orig.ident)){
    if(other_tool_tags[i] == dataset$Final_tags[i] & other_tool_tags[i] != 'Undetermined'){
      compared<-c(compared, "same_tag")
    } else if((other_tool_tags[i] == dataset$Final_tags[i] & other_tool_tags[i] == 'Undetermined') | (other_tool_tags[i]=="Multiplet" & dataset$Final_tags[i]=="Undetermined")){
      compared<-c(compared, "same_untagged")
    } else if((other_tool_tags[i] == "Undetermined" & dataset$Final_tags[i] != "Undetermined") | (other_tool_tags[i] == "Multiplet" & dataset$Final_tags[i] != "Undetermined")){
      compared<-c(compared, "rescued")
    } else if(other_tool_tags[i] != "Undetermined" & other_tool_tags[i] != "Multiplet" & dataset$Final_tags[i] != other_tool_tags[i] & dataset$Final_tags[i] != "Undetermined"){
      compared<-c(compared, "different_tag")
    } else{
      compared<-c(compared, "untagged_after")
    }
  }
  x<-paste("Comparison_RescueTag_vs_", tool_name, sep="")
  dataset[[x]]<-compared
  return(dataset)
}

ComparisonScatterPlots<-function(dataset,tool_name){
  if(tool_name=="HTODemux"){
    same_all<-subset(dataset, subset=Comparison_RescueTag_vs_HTODemux=="same_tag" | Comparison_RescueTag_vs_HTODemux=="same_untagged")
    rescued<-subset(dataset, subset=Comparison_RescueTag_vs_HTODemux=="rescued")
    different_tag<-subset(dataset, subset=Comparison_RescueTag_vs_HTODemux=="different_tag")
    untagged_after<-subset(dataset, subset=Comparison_RescueTag_vs_HTODemux=="untagged_after")
    comparison_column_same_all<-same_all$Comparison_RescueTag_vs_HTODemux
    comparison_column_rescued <- rescued$Comparison_RescueTag_vs_HTODemux
    comparison_column_different_tag <- different_tag$Comparison_RescueTag_vs_HTODemux
    comparison_column_untagged_after <- untagged_after$Comparison_RescueTag_vs_HTODemux
    df_for_plot<-data.frame("First_best"=same_all$Postnorm_first_best, "Delta"=same_all$Postnorm_delta, "Comparison"=comparison_column_same_all)
    same_plot<-ggplot(df_for_plot,aes(x=First_best,y=Delta,col=Comparison))+geom_point()+scale_color_manual(values = c("#ffe7a8ff", "#ccccccff"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+scale_x_continuous(limits = c(0,1), name = "First_best")+scale_y_continuous(limits=c(0,1), name="Delta")
    df_for_plot<-data.frame("First_best"=rescued$Postnorm_first_best, "Delta"=rescued$Postnorm_delta, "Comparison"=comparison_column_rescued)
    rescued_plot<-ggplot(df_for_plot,aes(x=First_best,y=Delta,col=Comparison))+geom_point()+scale_color_manual(values = c("#ffb380ff"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+scale_x_continuous(limits=c(0,1))+scale_y_continuous(limits = c(0,1))
    df_for_plot<-data.frame("First_best"=different_tag$Postnorm_first_best, "Delta"=different_tag$Postnorm_delta, "Comparison"=comparison_column_different_tag)
    different_tag_plot<-ggplot(df_for_plot,aes(x=First_best,y=Delta,col=Comparison))+geom_point()+scale_color_manual(values = c("#a8cbffff"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+scale_x_continuous(limits = c(0,1),name="First_best")+scale_y_continuous(limits = c(0,1), name = "Delta")
    df_for_plot<-data.frame("First_best"=untagged_after$Postnorm_first_best, "Delta"=untagged_after$Postnorm_delta, "Comparison"=comparison_column_untagged_after)
    untagged_after_plot<-ggplot(df_for_plot,aes(x=First_best,y=Delta,col=Comparison))+geom_point()+scale_color_manual(values = c("#ffa8a8ff"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+scale_x_continuous(limits = c(0,1),name="First_best")+scale_y_continuous(limits = c(0,1), name = "Delta")
    same_plot+rescued_plot+different_tag_plot+untagged_after_plot
  }
  else if(tool_name=="BFF_cluster"){
    same_all<-subset(dataset, subset=Comparison_RescueTag_vs_BFF_cluster=="same_tag" | Comparison_RescueTag_vs_BFF_cluster=="same_untagged")
    rescued<-subset(dataset, subset=Comparison_RescueTag_vs_BFF_cluster=="rescued")
    different_tag<-subset(dataset, subset=Comparison_RescueTag_vs_BFF_cluster=="different_tag")
    untagged_after<-subset(dataset, subset=Comparison_RescueTag_vs_BFF_cluster=="untagged_after")
    comparison_column_same_all<-same_all$Comparison_RescueTag_vs_BFF_cluster
    comparison_column_rescued <- rescued$Comparison_RescueTag_vs_BFF_cluster
    comparison_column_different_tag <- different_tag$Comparison_RescueTag_vs_BFF_cluster
    comparison_column_untagged_after <- untagged_after$Comparison_RescueTag_vs_BFF_cluster
    df_for_plot<-data.frame("First_best"=same_all$Postnorm_first_best, "Delta"=same_all$Postnorm_delta, "Comparison"=comparison_column_same_all)
    same_plot<-ggplot(df_for_plot,aes(x=First_best,y=Delta,col=Comparison))+geom_point()+scale_color_manual(values = c("#ffe7a8ff", "#ccccccff"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+scale_x_continuous(limits = c(0,1), name = "First_best")+scale_y_continuous(limits=c(0,1), name="Delta")
    df_for_plot<-data.frame("First_best"=rescued$Postnorm_first_best, "Delta"=rescued$Postnorm_delta, "Comparison"=comparison_column_rescued)
    rescued_plot<-ggplot(df_for_plot,aes(x=First_best,y=Delta,col=Comparison))+geom_point()+scale_color_manual(values = c("#ffb380ff"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+scale_x_continuous(limits=c(0,1))+scale_y_continuous(limits = c(0,1))
    df_for_plot<-data.frame("First_best"=different_tag$Postnorm_first_best, "Delta"=different_tag$Postnorm_delta, "Comparison"=comparison_column_different_tag)
    different_tag_plot<-ggplot(df_for_plot,aes(x=First_best,y=Delta,col=Comparison))+geom_point()+scale_color_manual(values = c("#a8cbffff"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+scale_x_continuous(limits = c(0,1),name="First_best")+scale_y_continuous(limits = c(0,1), name = "Delta")
    df_for_plot<-data.frame("First_best"=untagged_after$Postnorm_first_best, "Delta"=untagged_after$Postnorm_delta, "Comparison"=comparison_column_untagged_after)
    untagged_after_plot<-ggplot(df_for_plot,aes(x=First_best,y=Delta,col=Comparison))+geom_point()+scale_color_manual(values = c("#ffa8a8ff"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+scale_x_continuous(limits = c(0,1),name="First_best")+scale_y_continuous(limits = c(0,1), name = "Delta")
    same_plot+rescued_plot+different_tag_plot+untagged_after_plot
  }
  else if(tool_name=="BFF_raw"){
    same_all<-subset(dataset, subset=Comparison_RescueTag_vs_BFF_raw=="same_tag" | Comparison_RescueTag_vs_BFF_raw=="same_untagged")
    rescued<-subset(dataset, subset=Comparison_RescueTag_vs_BFF_raw=="rescued")
    different_tag<-subset(dataset, subset=Comparison_RescueTag_vs_BFF_raw=="different_tag")
    untagged_after<-subset(dataset, subset=Comparison_RescueTag_vs_BFF_raw=="untagged_after")
    comparison_column_same_all<-same_all$Comparison_RescueTag_vs_BFF_raw
    comparison_column_rescued <- rescued$Comparison_RescueTag_vs_BFF_raw
    comparison_column_different_tag <- different_tag$Comparison_RescueTag_vs_BFF_raw
    comparison_column_untagged_after <- untagged_after$Comparison_RescueTag_vs_BFF_raw
    df_for_plot<-data.frame("First_best"=same_all$Postnorm_first_best, "Delta"=same_all$Postnorm_delta, "Comparison"=comparison_column_same_all)
    same_plot<-ggplot(df_for_plot,aes(x=First_best,y=Delta,col=Comparison))+geom_point()+scale_color_manual(values = c("#ffe7a8ff", "#ccccccff"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+scale_x_continuous(limits = c(0,1), name = "First_best")+scale_y_continuous(limits=c(0,1), name="Delta")
    df_for_plot<-data.frame("First_best"=rescued$Postnorm_first_best, "Delta"=rescued$Postnorm_delta, "Comparison"=comparison_column_rescued)
    rescued_plot<-ggplot(df_for_plot,aes(x=First_best,y=Delta,col=Comparison))+geom_point()+scale_color_manual(values = c("#ffb380ff"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+scale_x_continuous(limits=c(0,1))+scale_y_continuous(limits = c(0,1))
    df_for_plot<-data.frame("First_best"=different_tag$Postnorm_first_best, "Delta"=different_tag$Postnorm_delta, "Comparison"=comparison_column_different_tag)
    different_tag_plot<-ggplot(df_for_plot,aes(x=First_best,y=Delta,col=Comparison))+geom_point()+scale_color_manual(values = c("#a8cbffff"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+scale_x_continuous(limits = c(0,1),name="First_best")+scale_y_continuous(limits = c(0,1), name = "Delta")
    df_for_plot<-data.frame("First_best"=untagged_after$Postnorm_first_best, "Delta"=untagged_after$Postnorm_delta, "Comparison"=comparison_column_untagged_after)
    untagged_after_plot<-ggplot(df_for_plot,aes(x=First_best,y=Delta,col=Comparison))+geom_point()+scale_color_manual(values = c("#ffa8a8ff"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+scale_x_continuous(limits = c(0,1),name="First_best")+scale_y_continuous(limits = c(0,1), name = "Delta")
    same_plot+rescued_plot+different_tag_plot+untagged_after_plot
  }
  else if(tool_name=="hashedDrops"){
    same_all<-subset(dataset, subset=Comparison_RescueTag_vs_hashedDrops=="same_tag" | Comparison_RescueTag_vs_hashedDrops=="same_untagged")
    rescued<-subset(dataset, subset=Comparison_RescueTag_vs_hashedDrops=="rescued")
    comparison_column_same_all<-same_all$Comparison_RescueTag_vs_hashedDrops
    comparison_column_rescued <- rescued$Comparison_RescueTag_vs_hashedDrops
    df_for_plot<-data.frame("First_best"=same_all$Postnorm_first_best, "Delta"=same_all$Postnorm_delta, "Comparison"=comparison_column_same_all)
    same_plot<-ggplot(df_for_plot,aes(x=First_best,y=Delta,col=Comparison))+geom_point()+scale_color_manual(values = c("#ffe7a8ff", "#ccccccff"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+scale_x_continuous(limits = c(0,1), name = "First_best")+scale_y_continuous(limits=c(0,1), name="Delta")
    df_for_plot<-data.frame("First_best"=rescued$Postnorm_first_best, "Delta"=rescued$Postnorm_delta, "Comparison"=comparison_column_rescued)
    rescued_plot<-ggplot(df_for_plot,aes(x=First_best,y=Delta,col=Comparison))+geom_point()+scale_color_manual(values = c("#ffb380ff"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+scale_x_continuous(limits=c(0,1))+scale_y_continuous(limits = c(0,1))
    same_plot+rescued_plot
  }
  else if(tool_name=="deMULTIplex"){
    same_all<-subset(dataset, subset=Comparison_RescueTag_vs_deMULTIplex=="same_tag" | Comparison_RescueTag_vs_deMULTIplex=="same_untagged")
    rescued<-subset(dataset, subset=Comparison_RescueTag_vs_deMULTIplex=="rescued")
    untagged_after<-subset(dataset, subset=Comparison_RescueTag_vs_deMULTIplex=="untagged_after")
    comparison_column_same_all<-same_all$Comparison_RescueTag_vs_deMULTIplex
    comparison_column_rescued <- rescued$Comparison_RescueTag_vs_deMULTIplex
    comparison_column_untagged_after <- untagged_after$Comparison_RescueTag_vs_deMULTIplex
    df_for_plot<-data.frame("First_best"=same_all$Postnorm_first_best, "Delta"=same_all$Postnorm_delta, "Comparison"=comparison_column_same_all)
    same_plot<-ggplot(df_for_plot,aes(x=First_best,y=Delta,col=Comparison))+geom_point()+scale_color_manual(values = c("#ffe7a8ff", "#ccccccff"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+scale_x_continuous(limits = c(0,1), name = "First_best")+scale_y_continuous(limits=c(0,1), name="Delta")
    df_for_plot<-data.frame("First_best"=rescued$Postnorm_first_best, "Delta"=rescued$Postnorm_delta, "Comparison"=comparison_column_rescued)
    rescued_plot<-ggplot(df_for_plot,aes(x=First_best,y=Delta,col=Comparison))+geom_point()+scale_color_manual(values = c("#ffb380ff"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+scale_x_continuous(limits=c(0,1))+scale_y_continuous(limits = c(0,1))
    df_for_plot<-data.frame("First_best"=untagged_after$Postnorm_first_best, "Delta"=untagged_after$Postnorm_delta, "Comparison"=comparison_column_untagged_after)
    untagged_after_plot<-ggplot(df_for_plot,aes(x=First_best,y=Delta,col=Comparison))+geom_point()+scale_color_manual(values = c("#ffa8a8ff"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+scale_x_continuous(limits = c(0,1),name="First_best")+scale_y_continuous(limits = c(0,1), name = "Delta")
    same_plot+rescued_plot+untagged_after_plot
  }
}

ComparisonRidgePlotsBest<-function(dataset, tool_tags_column){
    same_tag<-subset(dataset, subset=tool_tags_column == Final_tags & tool_tags_column != 'Undetermined' )
    same_undetermined<-subset(dataset, subset=(tool_tags_column == Final_tags & tool_tags_column == 'Undetermined') | (tool_tags_column=="Multiplet" & Final_tags=="Undetermined"))
    rescued_from_untagging<-subset(dataset, subset=(tool_tags_column == "Undetermined" & Final_tags != "Undetermined") | (tool_tags_column == "Multiplet" & Final_tags != "Undetermined"))
    different_tag_assigned<-subset(dataset, subset=tool_tags_column != "Undetermined" & tool_tags_column != "Multiplet" & Final_tags != tool_tags_column & Final_tags != "Undetermined")
    untagged_after_rescoring<-subset(dataset, subset=tool_tags_column != "Undetermined" & tool_tags_column != "Multiplet" & Final_tags == "Undetermined")
    x<-c(same_tag$Postnorm_first_best)
    x<-c(x,same_undetermined$Postnorm_first_best)
    x<-c(x,rescued_from_untagging$Postnorm_first_best)
    x<-c(x,different_tag_assigned$Postnorm_first_best)
    x<-c(x,untagged_after_rescoring$Postnorm_first_best)
    y<-c(rep('5_Same_tag_assigned', length(same_tag$Postnorm_first_best)))
    y<-c(y, rep('4_Same_cells_untagged', length(same_undetermined$Postnorm_first_best)))
    y<-c(y, rep('3_Rescued_from_untagging', length(rescued_from_untagging$Postnorm_first_best)))
    y<-c(y, rep('2_Different_tag_assigned', length(different_tag_assigned$Postnorm_first_best)))
    y<-c(y, rep('1_Untagged_after_scoring', length(untagged_after_rescoring$Postnorm_first_best)))
    df<-data.frame('Normalized_first_best_values'=x, 'Categories'=y)
    ggplot(df, aes(x=Normalized_first_best_values, y=Categories, fill=Categories))+geom_density_ridges(scale=1.5)+scale_x_continuous(expand = c(0,0), name="Normalized_first_best_values", breaks = seq(0,1, by=0.2))+theme_ridges(font_size = 19)+theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+scale_fill_manual(values = c("#ffa8a8ff","#a8cbffff", "#ffb380ff", "#ccccccff", "#ffe7a8ff"))                                                                                                                                    
}

ComparisonRidgePlotsBest_hashedDrops<-function(dataset, tool_tags_column){
  same_tag<-subset(dataset, subset=tool_tags_column == Final_tags & tool_tags_column != 'Undetermined' )
  same_undetermined<-subset(dataset, subset=(tool_tags_column == Final_tags & tool_tags_column == 'Undetermined') | (tool_tags_column=="Multiplet" & Final_tags=="Undetermined"))
  rescued_from_untagging<-subset(dataset, subset=(tool_tags_column == "Undetermined" & Final_tags != "Undetermined") | (tool_tags_column == "Multiplet" & Final_tags != "Undetermined"))
  x<-c(same_tag$Postnorm_first_best)
  x<-c(x,same_undetermined$Postnorm_first_best)
  x<-c(x,rescued_from_untagging$Postnorm_first_best)
  y<-c(rep('5_Same_tag_assigned', length(same_tag$Postnorm_first_best)))
  y<-c(y, rep('4_Same_cells_untagged', length(same_undetermined$Postnorm_first_best)))
  y<-c(y, rep('3_Rescued_from_untagging', length(rescued_from_untagging$Postnorm_first_best)))
  df<-data.frame('Normalized_first_best_values'=x, 'Categories'=y)
  ggplot(df, aes(x=Normalized_first_best_values, y=Categories, fill=Categories))+geom_density_ridges(scale=1.5)+scale_x_continuous(expand = c(0,0), name="Normalized_first_best_values", breaks = seq(0,1, by=0.2))+theme_ridges(font_size = 19)+theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+scale_fill_manual(values = c("#ffb380ff", "#ccccccff", "#ffe7a8ff"))                                                                                                                                    
  
}
ComparisonRidgePlotsBest_deMULTIplex<-function(dataset, tool_tags_column){
  same_tag<-subset(dataset, subset=tool_tags_column == Final_tags & tool_tags_column != 'Undetermined' )
  same_undetermined<-subset(dataset, subset=(tool_tags_column == Final_tags & tool_tags_column == 'Undetermined') | (tool_tags_column=="Multiplet" & Final_tags=="Undetermined"))
  rescued_from_untagging<-subset(dataset, subset=(tool_tags_column == "Undetermined" & Final_tags != "Undetermined") | (tool_tags_column == "Multiplet" & Final_tags != "Undetermined"))
  untagged_after_rescoring<-subset(dataset, subset=tool_tags_column != "Undetermined" & tool_tags_column != "Multiplet" & Final_tags == "Undetermined")
  x<-c(same_tag$Postnorm_first_best)
  x<-c(x,same_undetermined$Postnorm_first_best)
  x<-c(x,rescued_from_untagging$Postnorm_first_best)
  x<-c(x,untagged_after_rescoring$Postnorm_first_best)
  y<-c(rep('5_Same_tag_assigned', length(same_tag$Postnorm_first_best)))
  y<-c(y, rep('4_Same_cells_untagged', length(same_undetermined$Postnorm_first_best)))
  y<-c(y, rep('3_Rescued_from_untagging', length(rescued_from_untagging$Postnorm_first_best)))
  y<-c(y, rep('1_Untagged_after_scoring', length(untagged_after_rescoring$Postnorm_first_best)))
  df<-data.frame('Normalized_first_best_values'=x, 'Categories'=y)
  ggplot(df, aes(x=Normalized_first_best_values, y=Categories, fill=Categories))+geom_density_ridges(scale=1.5)+scale_x_continuous(expand = c(0,0), name="Normalized_first_best_values", breaks = seq(0,1, by=0.2))+theme_ridges(font_size = 19)+theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+scale_fill_manual(values = c("#ffa8a8ff", "#ffb380ff", "#ccccccff", "#ffe7a8ff"))                                                                                                                                    
}
ComparisonRidgePlotsDelta<-function(dataset, tool_tags_column){
  same_tag<-subset(dataset, subset=tool_tags_column == Final_tags & tool_tags_column != 'Undetermined' )
  same_undetermined<-subset(dataset, subset=(tool_tags_column == Final_tags & tool_tags_column == 'Undetermined') | (tool_tags_column=="Multiplet" & Final_tags=="Undetermined"))
  rescued_from_untagging<-subset(dataset, subset=(tool_tags_column == "Undetermined" & Final_tags != "Undetermined") | (tool_tags_column == "Multiplet" & Final_tags != "Undetermined"))
  different_tag_assigned<-subset(dataset, subset=tool_tags_column != "Undetermined" & tool_tags_column != "Multiplet" & Final_tags != tool_tags_column & Final_tags != "Undetermined")
  untagged_after_rescoring<-subset(dataset, subset=tool_tags_column != "Undetermined" & tool_tags_column != "Multiplet" & Final_tags == "Undetermined")
  x<-c(same_tag$Postnorm_delta)
  x<-c(x,same_undetermined$Postnorm_delta)
  x<-c(x,rescued_from_untagging$Postnorm_delta)
  x<-c(x,different_tag_assigned$Postnorm_delta)
  x<-c(x,untagged_after_rescoring$Postnorm_delta)
  y<-c(rep('5_Same_tag_assigned', length(same_tag$Postnorm_delta)))
  y<-c(y, rep('4_Same_cells_untagged', length(same_undetermined$Postnorm_delta)))
  y<-c(y, rep('3_Rescued_from_untagging', length(rescued_from_untagging$Postnorm_delta)))
  y<-c(y, rep('2_Different_tag_assigned', length(different_tag_assigned$Postnorm_delta)))
  y<-c(y, rep('1_Untagged_after_scoring', length(untagged_after_rescoring$Postnorm_delta)))
  df<-data.frame('Normalized_delta_values'=x, 'Categories'=y)
  ggplot(df, aes(x=Normalized_delta_values, y=Categories, fill=Categories))+geom_density_ridges(scale=1.5)+scale_x_continuous(expand = c(0,0), name="Normalized_delta_values", breaks = seq(0,1, by=0.2))+theme_ridges(font_size = 19)+theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+scale_fill_manual(values = c("#ffa8a8ff","#a8cbffff", "#ffb380ff", "#ccccccff", "#ffe7a8ff"))                                                                                                                                    
  
}

ComparisonRidgePlotsDelta_hashedDrops<-function(dataset, tool_tags_column){
  same_tag<-subset(dataset, subset=tool_tags_column == Final_tags & tool_tags_column != 'Undetermined' )
  same_undetermined<-subset(dataset, subset=(tool_tags_column == Final_tags & tool_tags_column == 'Undetermined') | (tool_tags_column=="Multiplet" & Final_tags=="Undetermined"))
  rescued_from_untagging<-subset(dataset, subset=(tool_tags_column == "Undetermined" & Final_tags != "Undetermined") | (tool_tags_column == "Multiplet" & Final_tags != "Undetermined"))
  x<-c(same_tag$Postnorm_delta)
  x<-c(x,same_undetermined$Postnorm_delta)
  x<-c(x,rescued_from_untagging$Postnorm_delta)
  y<-c(rep('5_Same_tag_assigned', length(same_tag$Postnorm_delta)))
  y<-c(y, rep('4_Same_cells_untagged', length(same_undetermined$Postnorm_delta)))
  y<-c(y, rep('3_Rescued_from_untagging', length(rescued_from_untagging$Postnorm_delta)))
  df<-data.frame('Normalized_delta_values'=x, 'Categories'=y)
  ggplot(df, aes(x=Normalized_delta_values, y=Categories, fill=Categories))+geom_density_ridges(scale=1.5)+scale_x_continuous(expand = c(0,0), name="Normalized_delta_values", breaks = seq(0,1, by=0.2))+theme_ridges(font_size = 19)+theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+scale_fill_manual(values = c("#ffb380ff", "#ccccccff", "#ffe7a8ff"))                                                                                                                                    
}

ComparisonRidgePlotsDelta_deMULTIplex<-function(dataset, tool_tags_column){
  same_tag<-subset(dataset, subset=tool_tags_column == Final_tags & tool_tags_column != 'Undetermined' )
  same_undetermined<-subset(dataset, subset=(tool_tags_column == Final_tags & tool_tags_column == 'Undetermined') | (tool_tags_column=="Multiplet" & Final_tags=="Undetermined"))
  rescued_from_untagging<-subset(dataset, subset=(tool_tags_column == "Undetermined" & Final_tags != "Undetermined") | (tool_tags_column == "Multiplet" & Final_tags != "Undetermined"))
  untagged_after_rescoring<-subset(dataset, subset=tool_tags_column != "Undetermined" & tool_tags_column != "Multiplet" & Final_tags == "Undetermined")
  x<-c(same_tag$Postnorm_delta)
  x<-c(x,same_undetermined$Postnorm_delta)
  x<-c(x,rescued_from_untagging$Postnorm_delta)
  x<-c(x,untagged_after_rescoring$Postnorm_delta)
  y<-c(rep('5_Same_tag_assigned', length(same_tag$Postnorm_delta)))
  y<-c(y, rep('4_Same_cells_untagged', length(same_undetermined$Postnorm_delta)))
  y<-c(y, rep('3_Rescued_from_untagging', length(rescued_from_untagging$Postnorm_delta)))
  y<-c(y, rep('1_Untagged_after_scoring', length(untagged_after_rescoring$Postnorm_delta)))
  df<-data.frame('Normalized_delta_values'=x, 'Categories'=y)
  ggplot(df, aes(x=Normalized_delta_values, y=Categories, fill=Categories))+geom_density_ridges(scale=1.5)+scale_x_continuous(expand = c(0,0), name="Normalized_delta_values", breaks = seq(0,1, by=0.2))+theme_ridges(font_size = 19)+theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+scale_fill_manual(values = c("#ffa8a8ff", "#ffb380ff", "#ccccccff", "#ffe7a8ff"))                                                                                                                                    
  
}

ComparisonPercentages<-function(dataset, tool_tags_column){
  same_tag<-subset(dataset, subset=tool_tags_column == Final_tags & tool_tags_column != 'Undetermined' )
  same_undetermined<-subset(dataset, subset=(tool_tags_column == Final_tags & tool_tags_column == 'Undetermined') | (tool_tags_column=="Multiplet" & Final_tags=="Undetermined"))
  rescued_from_untagging<-subset(dataset, subset=(tool_tags_column == "Undetermined" & Final_tags != "Undetermined") | (tool_tags_column == "Multiplet" & Final_tags != "Undetermined"))
  different_tag_assigned<-subset(dataset, subset=tool_tags_column != "Undetermined" & tool_tags_column != "Multiplet" & Final_tags != tool_tags_column & Final_tags != "Undetermined")
  untagged_after_rescoring<-subset(dataset, subset=tool_tags_column != "Undetermined" & tool_tags_column != "Multiplet" & Final_tags == "Undetermined")
  print(paste0("Same tag: ",length(same_tag$orig.ident)/length(dataset$orig.ident)*100))
  print(paste0("Same no tag: ",length(same_undetermined$orig.ident)/length(dataset$orig.ident)*100))
  print(paste0("Rescued: ",length(rescued_from_untagging$orig.ident)/length(dataset$orig.ident)*100))
  print(paste0("Different tag: ", length(different_tag_assigned$orig.ident)/length(dataset$orig.ident)*100))
  print(paste0("Untagged: ",length(untagged_after_rescoring$orig.ident)/length(dataset$orig.ident)*100))
  print(paste0("Total: ", length(same_tag$orig.ident)/length(dataset$orig.ident)+length(same_undetermined$orig.ident)/length(dataset$orig.ident)+length(rescued_from_untagging$orig.ident)/length(dataset$orig.ident)+length(different_tag_assigned$orig.ident)/length(dataset$orig.ident)+length(untagged_after_rescoring$orig.ident)/length(dataset$orig.ident)))
  
}

ComparisonPercentages_hashedDrops<-function(dataset, tool_tags_column){
  same_tag<-subset(dataset, subset=tool_tags_column == Final_tags & tool_tags_column != 'Undetermined' )
  same_undetermined<-subset(dataset, subset=(tool_tags_column == Final_tags & tool_tags_column == 'Undetermined') | (tool_tags_column=="Multiplet" & Final_tags=="Undetermined"))
  rescued_from_untagging<-subset(dataset, subset=(tool_tags_column == "Undetermined" & Final_tags != "Undetermined") | (tool_tags_column == "Multiplet" & Final_tags != "Undetermined"))
  print(paste0("Same tag: ",length(same_tag$orig.ident)/length(dataset$orig.ident)*100))
  print(paste0("Same no tag: ",length(same_undetermined$orig.ident)/length(dataset$orig.ident)*100))
  print(paste0("Rescued: ",length(rescued_from_untagging$orig.ident)/length(dataset$orig.ident)*100))
  print(paste0("Total: ", length(same_tag$orig.ident)/length(dataset$orig.ident)+length(same_undetermined$orig.ident)/length(dataset$orig.ident)+length(rescued_from_untagging$orig.ident)/length(dataset$orig.ident)))
  
}

ComparisonPercentages_deMULTIplex<-function(dataset, tool_tags_column){
  same_tag<-subset(dataset, subset=tool_tags_column == Final_tags & tool_tags_column != 'Undetermined' )
  same_undetermined<-subset(dataset, subset=(tool_tags_column == Final_tags & tool_tags_column == 'Undetermined') | (tool_tags_column=="Multiplet" & Final_tags=="Undetermined"))
  rescued_from_untagging<-subset(dataset, subset=(tool_tags_column == "Undetermined" & Final_tags != "Undetermined") | (tool_tags_column == "Multiplet" & Final_tags != "Undetermined"))
  untagged_after_rescoring<-subset(dataset, subset=tool_tags_column != "Undetermined" & tool_tags_column != "Multiplet" & Final_tags == "Undetermined")
  print(paste0("Same tag: ",length(same_tag$orig.ident)/length(dataset$orig.ident)*100))
  print(paste0("Same no tag: ",length(same_undetermined$orig.ident)/length(dataset$orig.ident)*100))
  print(paste0("Rescued: ",length(rescued_from_untagging$orig.ident)/length(dataset$orig.ident)*100))
  print(paste0("Untagged: ",length(untagged_after_rescoring$orig.ident)/length(dataset$orig.ident)*100))
  print(paste0("Total: ", length(same_tag$orig.ident)/length(dataset$orig.ident)+length(same_undetermined$orig.ident)/length(dataset$orig.ident)+length(rescued_from_untagging$orig.ident)/length(dataset$orig.ident)+length(untagged_after_rescoring$orig.ident)/length(dataset$orig.ident)))
}

##############################################
#                                            #
#           RescueTag benchmarking           #
#                                            #
##############################################

###################
#                 #
#   HTO dataset   #
#                 #
###################

#HTODemux

#Creating an object with the data
HTO_raw_data<-Read10X(path_to_the_raw_data_folder) # type a path to the "Human_PSC_derived_cardiomyocytes" folder
HTO_dataset<-CreateSeuratObject(counts=HTO_raw_data$`Gene Expression`)
HTO_dataset[["HTO"]]<-CreateAssayObject(counts = HTO_raw_data$`Antibody Capture`)

#Preparing a sample tag expression table for the future usage
write.csv(t(HTO_dataset$`Antibody Capture`), file="HTO_sample_tag_table.csv")
HTO_table_HTO_dataset<-read.csv("HTO_sample_tag_table.csv", row.names = 1)
HTO_table_HTO_dataset<-t(HTO_table_HTO_dataset)

#Preparing the dataset and running HTODemux
HTO_dataset <- NormalizeData(HTO_dataset)
HTO_dataset <- FindVariableFeatures(HTO_dataset, selection.method = "mean.var.plot")
HTO_dataset <- ScaleData(HTO_dataset, features = VariableFeatures(HTO_dataset))
HTO_dataset <- NormalizeData(HTO_dataset, assay = "HTO", normalization.method = "CLR")
HTO_dataset <- HTODemux(HTO_dataset, assay = "HTO", positive.quantile = 0.99)

#Renaming the sample tags
HTODemux_tags<-c()
for(i in 1:length(HTO_dataset$orig.ident)){
  if(HTO_dataset$hash.ID[i]=="Hashtag1"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_1")
  }
  else if(HTO_dataset$hash.ID[i]=="Hashtag2"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_2")
  }
  else if(HTO_dataset$hash.ID[i]=="Hashtag3"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_3")
  }
  else if(HTO_dataset$hash.ID[i]=="Hashtag4"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_4")
  }
  else if(HTO_dataset$hash.ID[i]=="Hashtag5"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_5")
  }
  else if(HTO_dataset$hash.ID[i]=="Hashtag6"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_6")
  }
  else if(HTO_dataset$hash.ID[i]=="Hashtag7"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_7")
  }
  else if(HTO_dataset$hash.ID[i]=="Hashtag8"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_8")
  }
  else if(HTO_dataset$hash.ID[i]=="Hashtag9"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_9")
  }
  else if(HTO_dataset$hash.ID[i]=="Hashtag10"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_10")
  }
  else if(HTO_dataset$hash.ID[i]=="Hashtag12"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_11")
  }
  else if(HTO_dataset$hash.ID[i]=="Hashtag14"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_12")
  }
  else if(HTO_dataset$hash.ID[i]=="Hashtag15"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_13")
  }
  else if(HTO_dataset$hash.ID[i]=="Doublet"){
    HTODemux_tags<-c(HTODemux_tags, "Multiplet")
  }
  else{
    HTODemux_tags<-c(HTODemux_tags, "Undetermined")
  }
}
HTO_dataset[["HTODemux_tags"]]<-HTODemux_tags

#BFF

#Running BFF
calls <- GenerateCellHashingCalls(barcodeMatrix = HTO_table_HTO_dataset, methods = c('bff_raw', 'bff_cluster'))
write.csv(calls, file="HTO_dataset_BFF_calls.csv")
calls<-read.csv("HTO_BFF_calls.csv", row.names = "cellbarcodes")
calls$X<-c()
BFF_cluster_primary_tags<-calls$bff_cluster

#BFF cluster

#Renaming sample tags
BFF_cluster_tags<-c()
for(i in 1:length(HTO_dataset$orig.ident)){
  if(BFF_cluster_primary_tags[i]=="Hashtag1"){
    BFF_cluster_tags<-c(BFF_cluster_tags, "Tag_1")
  }
  else if(BFF_cluster_primary_tags[i]=="Hashtag2"){
    BFF_cluster_tags<-c(BFF_cluster_tags, "Tag_2")
  }
  else if(BFF_cluster_primary_tags[i]=="Hashtag3"){
    BFF_cluster_tags<-c(BFF_cluster_tags, "Tag_3")
  }
  else if(BFF_cluster_primary_tags[i]=="Hashtag4"){
    BFF_cluster_tags<-c(BFF_cluster_tags, "Tag_4")
  }
  else if(BFF_cluster_primary_tags[i]=="Hashtag5"){
    BFF_cluster_tags<-c(BFF_cluster_tags, "Tag_5")
  }
  else if(BFF_cluster_primary_tags[i]=="Hashtag6"){
    BFF_cluster_tags<-c(BFF_cluster_tags, "Tag_6")
  }
  else if(BFF_cluster_primary_tags[i]=="Hashtag7"){
    BFF_cluster_tags<-c(BFF_cluster_tags, "Tag_7")
  }
  else if(BFF_cluster_primary_tags[i]=="Hashtag8"){
    BFF_cluster_tags<-c(BFF_cluster_tags, "Tag_8")
  }
  else if(BFF_cluster_primary_tags[i]=="Hashtag9"){
    BFF_cluster_tags<-c(BFF_cluster_tags, "Tag_9")
  }
  else if(BFF_cluster_primary_tags[i]=="Hashtag10"){
    BFF_cluster_tags<-c(BFF_cluster_tags, "Tag_10")
  }
  else if(BFF_cluster_primary_tags[i]=="Hashtag12"){
    BFF_cluster_tags<-c(BFF_cluster_tags, "Tag_11")
  }
  else if(BFF_cluster_primary_tags[i]=="Hashtag14"){
    BFF_cluster_tags<-c(BFF_cluster_tags, "Tag_12")
  }
  else if(BFF_cluster_primary_tags[i]=="Hashtag15"){
    BFF_cluster_tags<-c(BFF_cluster_tags, "Tag_13")
  }
  else if(BFF_cluster_primary_tags[i]=="Doublet"){
    BFF_cluster_tags<-c(BFF_cluster_tags, "Multiplet")
  }
  else{
    BFF_cluster_tags<-c(BFF_cluster_tags, "Undetermined")
  }
}
HTO_dataset[["BFF_cluster_tags"]]<-BFF_cluster_tags

#BFF raw

BFF_raw_primary_tags<-calls$bff_raw

#Renaming sample tags
BFF_raw_tags<-c()
for(i in 1:length(HTO_dataset$orig.ident)){
  if(BFF_raw_primary_tags[i]=="Hashtag1"){
    BFF_raw_tags<-c(BFF_raw_tags, "Tag_1")
  }
  else if(BFF_raw_primary_tags[i]=="Hashtag2"){
    BFF_raw_tags<-c(BFF_raw_tags, "Tag_2")
  }
  else if(BFF_raw_primary_tags[i]=="Hashtag3"){
    BFF_raw_tags<-c(BFF_raw_tags, "Tag_3")
  }
  else if(BFF_raw_primary_tags[i]=="Hashtag4"){
    BFF_raw_tags<-c(BFF_raw_tags, "Tag_4")
  }
  else if(BFF_raw_primary_tags[i]=="Hashtag5"){
    BFF_raw_tags<-c(BFF_raw_tags, "Tag_5")
  }
  else if(BFF_raw_primary_tags[i]=="Hashtag6"){
    BFF_raw_tags<-c(BFF_raw_tags, "Tag_6")
  }
  else if(BFF_raw_primary_tags[i]=="Hashtag7"){
    BFF_raw_tags<-c(BFF_raw_tags, "Tag_7")
  }
  else if(BFF_raw_primary_tags[i]=="Hashtag8"){
    BFF_raw_tags<-c(BFF_raw_tags, "Tag_8")
  }
  else if(BFF_raw_primary_tags[i]=="Hashtag9"){
    BFF_raw_tags<-c(BFF_raw_tags, "Tag_9")
  }
  else if(BFF_raw_primary_tags[i]=="Hashtag10"){
    BFF_raw_tags<-c(BFF_raw_tags, "Tag_10")
  }
  else if(BFF_raw_primary_tags[i]=="Hashtag12"){
    BFF_raw_tags<-c(BFF_raw_tags, "Tag_11")
  }
  else if(BFF_raw_primary_tags[i]=="Hashtag14"){
    BFF_raw_tags<-c(BFF_raw_tags, "Tag_12")
  }
  else if(BFF_raw_primary_tags[i]=="Hashtag15"){
    BFF_raw_tags<-c(BFF_raw_tags, "Tag_13")
  }
  else if(BFF_raw_primary_tags[i]=="Doublet"){
    BFF_raw_tags<-c(BFF_raw_tags, "Multiplet")
  }
  else{
    BFF_raw_tags<-c(BFF_raw_tags, "Undetermined")
  }
}
HTO_dataset[["BFF_raw_tags"]]<-BFF_raw_tags

#hashedDrops

#Running hashedDrops
hashedDrops_demultiplexing<-hashedDrops(HTO_table_HTO_dataset)
#Creating a vector with hashedDrops tags
hashedDrops_tags<-c()
for(i in 1:length(hashedDrops_demultiplexing$Total)){
  if(hashedDrops_demultiplexing$Confident[i] == TRUE & hashedDrops_demultiplexing$Doublet[i] == FALSE){
    hashedDrops_tags<-c(hashedDrops_tags, hashedDrops_demultiplexing$Best[i])
  }
  else if(hashedDrops_demultiplexing$Confident[i] == FALSE & hashedDrops_demultiplexing$Doublet[i] == TRUE){
    hashedDrops_tags<-c(hashedDrops_tags, "Multiplet")
  }
  else if(hashedDrops_demultiplexing$Confident[i] == FALSE & hashedDrops_demultiplexing$Doublet[i] == FALSE){
    hashedDrops_tags<-c(hashedDrops_tags, "Undetermined")
  }
}
#Renaming the tags
final_hashedDrops_tags<-c()
for(i in 1:length(hashedDrops_tags)){
  if(hashedDrops_tags[i]=="1"){
    final_hashedDrops_tags<-c(final_hashedDrops_tags, "Tag_1")
  }
  else if(hashedDrops_tags[i]=="2"){
    final_hashedDrops_tags<-c(final_hashedDrops_tags, "Tag_2")
  }
  else if(hashedDrops_tags[i]=="3"){
    final_hashedDrops_tags<-c(final_hashedDrops_tags, "Tag_3")
  }
  else if(hashedDrops_tags[i]=="4"){
    final_hashedDrops_tags<-c(final_hashedDrops_tags, "Tag_4")
  }
  else if(hashedDrops_tags[i]=="5"){
    final_hashedDrops_tags<-c(final_hashedDrops_tags, "Tag_5")
  }
  else if(hashedDrops_tags[i]=="6"){
    final_hashedDrops_tags<-c(final_hashedDrops_tags, "Tag_6")
  }
  else if(hashedDrops_tags[i]=="7"){
    final_hashedDrops_tags<-c(final_hashedDrops_tags, "Tag_7")
  }
  else if(hashedDrops_tags[i]=="8"){
    final_hashedDrops_tags<-c(final_hashedDrops_tags, "Tag_8")
  }
  else if(hashedDrops_tags[i]=="9"){
    final_hashedDrops_tags<-c(final_hashedDrops_tags, "Tag_9")
  }
  else if(hashedDrops_tags[i]=="10"){
    final_hashedDrops_tags<-c(final_hashedDrops_tags, "Tag_10")
  }
  else if(hashedDrops_tags[i]=="11"){
    final_hashedDrops_tags<-c(final_hashedDrops_tags, "Tag_11")
  }
  else if(hashedDrops_tags[i]=="12"){
    final_hashedDrops_tags<-c(final_hashedDrops_tags, "Tag_12")
  }
  else if(hashedDrops_tags[i]=="13"){
    final_hashedDrops_tags<-c(final_hashedDrops_tags, "Tag_13")
  }
  else if(hashedDrops_tags[i]=="Multiplet"){
    final_hashedDrops_tags<-c(final_hashedDrops_tags, "Multiplet")
  }
  else if(hashedDrops_tags[i]=="Undetermined"){
    final_hashedDrops_tags<-c(final_hashedDrops_tags, "Undetermined")
  }
}
HTO_dataset[["hashedDrops_tags"]]<-final_hashedDrops_tags

#RescueTag
HTO_dataset <- ScoringTags("HTO_sample_tag_table.csv", 1, 13, 1, 13, HTO_dataset) #if the HTO_sample_tag_table.csv file is not in the working directory, write a path to it 
HTO_dataset <- PreQCFilter(HTO_dataset)
HTO_dataset <- NormToDepth(HTO_dataset)
df_for_plot<-data.frame('First_best'=HTO_dataset$NTD_first_best)
df_for_plot %>% ggplot(aes(x=First_best))+geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+expand_limits(x=c(0,1))
HTO_dataset <- FilterLowSeqDepth(HTO_dataset) # for the distribution write either "normal" or "bimodal"
HTO_dataset <- NormTagQCParams(HTO_dataset)
HTO_dataset[["Ratio"]]<-HTO_dataset$Postnorm_delta/HTO_dataset$Postnorm_first_best
HTO_dataset<-FinalTagging(HTO_dataset)

#Comparison of HTODemux and RescueTag
HTODemux_tags<-HTO_dataset$HTODemux_tags
HTO_dataset<-CompareRescueTag(HTO_dataset, HTODemux_tags, "HTODemux")

#Scatter plots for main dataset, RescueTag vs HTODemux
ComparisonScatterPlots(HTO_dataset, "HTODemux")

#Displaying ridge plots, main dataset, RescueTag vs HTODemux
ComparisonRidgePlotsBest(HTO_dataset, HTODemux_tags)
ComparisonRidgePlotsDelta(HTO_dataset, HTODemux_tags)

#Displaying percentages of cells from each category in comparison
ComparisonPercentages(HTO_dataset, HTODemux_tags)

#Comparison plots for main dataset, RescueTag vs BFF cluster
BFF_cluster_tags<-HTO_dataset$BFF_cluster_tags
HTO_dataset<-CompareRescueTag(HTO_dataset, BFF_cluster_tags, "BFF_cluster")
ComparisonScatterPlots(HTO_dataset, "BFF_cluster")
ComparisonRidgePlotsBest(HTO_dataset, BFF_cluster_tags)
ComparisonRidgePlotsDelta(HTO_dataset, BFF_cluster_tags)
ComparisonPercentages(HTO_dataset, BFF_cluster_tags)

#Comparison plots for main dataset, RescueTag vs BFF raw
BFF_raw_tags<-HTO_dataset$BFF_raw_tags
HTO_dataset<-CompareRescueTag(HTO_dataset, BFF_raw_tags, "BFF_raw")
ComparisonScatterPlots(HTO_dataset, "BFF_raw")
ComparisonRidgePlotsBest(HTO_dataset, BFF_raw_tags)
ComparisonRidgePlotsDelta(HTO_dataset, BFF_raw_tags)
ComparisonPercentages(HTO_dataset, BFF_raw_tags)

#Comparison plots for main dataset, RescueTag vs hashedDrops
hashedDrops_tags<-HTO_dataset$hashedDrops_tags
HTO_dataset<-CompareRescueTag(HTO_dataset, hashedDrops_tags, "hashedDrops")
ComparisonScatterPlots(HTO_dataset, "hashedDrops")
ComparisonRidgePlotsBest_hashedDrops(HTO_dataset, hashedDrops_tags)
ComparisonRidgePlotsDelta_hashedDrops(HTO_dataset, hashedDrops_tags)
ComparisonPercentages_hashedDrops(HTO_dataset, hashedDrops_tags)

###################
#                 #
#   CMO dataset   #
#                 #
###################

#HTODemux

#Creating an object with the dataset
CMO_raw_data<-Read10X_h5("SC3_v3_NextGem_DI_CellPlex_Nuclei_30K_Multiplex_count_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE) #if the H5 file is not in your working directory, type the whole path to it as a first argument of the function
CMO_dataset<-CreateSeuratObject(counts=CMO_raw_data$`Gene Expression`)
CMO_dataset[['CMO']] = CreateAssayObject(counts = CMO_raw_data$`Multiplexing Capture`)
cells <- fread("SC3_v3_NextGem_DI_CellPlex_Nuclei_30K_Multiplex_multiplexing_analysis_assignment_confidence_table.csv",select = c("Barcodes")) #if the csv file is not in your working directory, type the whole path to it as a first argument of the function
CMO_dataset <- subset(CMO_dataset, cells = cells$Barcodes)
DefaultAssay(CMO_dataset)<-"CMO"
CMO_dataset = subset(x = CMO_dataset, subset = nCount_CMO > 0)
#Preparing dataset for HTODemux
CMO_dataset <- NormalizeData(CMO_dataset, assay = "CMO", normalization.method = "CLR")
#Running HTODemux
CMO_dataset <- HTODemux(CMO_dataset, assay = "CMO", positive.quantile = 0.99)
write.csv(t(CMO_dataset@assays$CMO$counts), file='CMO_Brain_tag_counts.csv')
CMO_table_CMO_dataset<-read.csv("CMO_Brain_tag_counts.csv", row.names = 1)

#Renaming the sample tags
HTODemux_tags<-c()
for(i in 1:length(CMO_dataset$orig.ident)){
  if(CMO_dataset$hash.ID[i]=="CMO301"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_1")
  }
  else if(CMO_dataset$hash.ID[i]=="CMO302"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_2")
  }
  else if(CMO_dataset$hash.ID[i]=="CMO303"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_3")
  }
  else if(CMO_dataset$hash.ID[i]=="CMO304"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_4")
  }
  else if(CMO_dataset$hash.ID[i]=="CMO305"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_5")
  }
  else if(CMO_dataset$hash.ID[i]=="CMO306"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_6")
  }
  else if(CMO_dataset$hash.ID[i]=="CMO307"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_7")
  }
  else if(CMO_dataset$hash.ID[i]=="CMO308"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_8")
  }
  else if(CMO_dataset$hash.ID[i]=="CMO309"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_9")
  }
  else if(CMO_dataset$hash.ID[i]=="CMO310"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_10")
  }
  else if(CMO_dataset$hash.ID[i]=="CMO311"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_11")
  }
  else if(CMO_dataset$hash.ID[i]=="CMO312"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_12")
  }
  else if(CMO_dataset$hash.ID[i]=="Doublet"){
    HTODemux_tags<-c(HTODemux_tags, "Multiplet")
  }
  else{
    HTODemux_tags<-c(HTODemux_tags, "Undetermined")
  }
}
CMO_dataset[["HTODemux_tags"]]<-HTODemux_tags

#deMULTIplex

#Running deMULTIplex
bar.table <- CMO_table_CMO_dataset
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results1 <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results1$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "none") +
  geom_vline(xintercept=threshold.results1$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))
round1.calls <- classifyCells(bar.table, q=findQ(threshold.results1$res, threshold.results1$extrema))
neg.cells <- names(round1.calls)[which(round1.calls == "Negative")]
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

threshold.results2 <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results2$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "none") +
  geom_vline(xintercept=threshold.results2$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))
round2.calls <- classifyCells(bar.table, q=findQ(threshold.results2$res, threshold.results2$extrema))
neg.cells <- c(neg.cells, names(round2.calls)[which(round2.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

threshold.results3 <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results3$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "none") +
  geom_vline(xintercept=threshold.results1$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))
round3.calls <- classifyCells(bar.table, q=findQ(threshold.results3$res, threshold.results3$extrema))
neg.cells <- c(neg.cells, names(round3.calls)[which(round3.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

threshold.results4 <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results4$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "none") +
  geom_vline(xintercept=threshold.results1$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))
round4.calls <- classifyCells(bar.table, q=findQ(threshold.results4$res, threshold.results4$extrema))
neg.cells <- c(neg.cells, names(round4.calls)[which(round4.calls == "Negative")])

final.calls <- c(round4.calls, rep("Negative",length(neg.cells)))
names(final.calls) <- c(names(round4.calls),neg.cells)

#Ordering sample tags in a proper order
final_ordered_tags<-c()
for(i in 1:length(final.calls)){
  final_ordered_tags<-c(final_ordered_tags, final.calls[rownames(CMO_table_CMO_dataset[i,])])
}
names(final_ordered_tags)<-c()

#Renaming sample tags
deMULTIplex_tags<-c()
for(i in 1:length(CMO_dataset$orig.ident)){
  if(final_ordered_tags[i]=="CMO301"){
    deMULTIplex_tags<-c(deMULTIplex_tags, "Tag_1")
  }
  else if(final_ordered_tags[i]=="CMO302"){
    deMULTIplex_tags<-c(deMULTIplex_tags, "Tag_2")
  }
  else if(final_ordered_tags[i]=="CMO303"){
    deMULTIplex_tags<-c(deMULTIplex_tags, "Tag_3")
  }
  else if(final_ordered_tags[i]=="CMO304"){
    deMULTIplex_tags<-c(deMULTIplex_tags, "Tag_4")
  }
  else if(final_ordered_tags[i]=="CMO305"){
    deMULTIplex_tags<-c(deMULTIplex_tags, "Tag_5")
  }
  else if(final_ordered_tags[i]=="CMO306"){
    deMULTIplex_tags<-c(deMULTIplex_tags, "Tag_6")
  }
  else if(final_ordered_tags[i]=="CMO307"){
    deMULTIplex_tags<-c(deMULTIplex_tags, "Tag_7")
  }
  else if(final_ordered_tags[i]=="CMO308"){
    deMULTIplex_tags<-c(deMULTIplex_tags, "Tag_8")
  }
  else if(final_ordered_tags[i]=="CMO309"){
    deMULTIplex_tags<-c(deMULTIplex_tags, "Tag_9")
  }
  else if(final_ordered_tags[i]=="CMO310"){
    deMULTIplex_tags<-c(deMULTIplex_tags, "Tag_10")
  }
  else if(final_ordered_tags[i]=="CMO311"){
    deMULTIplex_tags<-c(deMULTIplex_tags, "Tag_11")
  }
  else if(final_ordered_tags[i]=="CMO312"){
    deMULTIplex_tags<-c(deMULTIplex_tags, "Tag_12")
  }
  else if(final_ordered_tags[i]=="Doublet"){
    deMULTIplex_tags<-c(deMULTIplex_tags, "Multiplet")
  }
  else{
    deMULTIplex_tags<-c(deMULTIplex_tags, "Undetermined")
  }
}
CMO_dataset[["deMULTIplex_tags"]]<-deMULTIplex_tags

# RescueTag

CMO_dataset <- ScoringTags("CMO_Brain_tag_counts.csv", 1, 12, 1, 12, CMO_dataset) #if the CMO_Brain_tag_counts.csv file is not in the working directory, write a path to it
CMO_dataset <- PreQCFilter(CMO_dataset)
CMO_dataset <- NormToDepth(CMO_dataset)
df_for_plot<-data.frame('First_best'=CMO_dataset$NTD_first_best)
df_for_plot %>% ggplot(aes(x=First_best))+geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+expand_limits(x=c(0,1))
CMO_dataset <- FilterLowSeqDepth(CMO_dataset, "bimodal") 
CMO_dataset <- NormTagQCParams(CMO_dataset)
CMO_dataset[["Ratio"]]<-CMO_dataset$Postnorm_delta/CMO_dataset$Postnorm_first_best
CMO_dataset<-FinalTagging(CMO_dataset)

#Comparison of HTODemux and RescueTag

HTODemux_tags<-CMO_dataset$HTODemux_tags
CMO_dataset<-CompareRescueTag(CMO_dataset, HTODemux_tags, "HTODemux")

#Scatter plots for CMO dataset, RescueTag vs HTODemux
ComparisonScatterPlots(CMO_dataset, "HTODemux")

#Displaying ridge plots, CMO dataset, RescueTag vs HTODemux
ComparisonRidgePlotsBest(CMO_dataset, HTODemux_tags)
ComparisonRidgePlotsDelta(CMO_dataset, HTODemux_tags)


#Comparison of deMULTIplex and RescueTag
deMULTIplex_tags<-CMO_dataset$deMULTIplex_tags
CMO_dataset<-CompareRescueTag(CMO_dataset, deMULTIplex_tags, "deMULTIplex")

#Scatter plots for CMO dataset, RescueTag vs deMULTIplex
ComparisonScatterPlots(CMO_dataset, "deMULTIplex")

#Displaying ridge plots, CMO dataset, RescueTag vs deMULTIplex
ComparisonRidgePlotsBest(CMO_dataset, deMULTIplex_tags)
ComparisonRidgePlotsDelta(CMO_dataset, deMULTIplex_tags)

####################
#                  #
#   Main dataset   #
#                  #
####################

# Reading the dataset and setting up the sample tag expression table
Main_dataset<-readRDS("Main_dataset_raw_seurat.rds") #if the file is not in the working directory, type the path to it as an argument to the function
HTO_table_main_dataset<-read.csv("Main_dataset_Sample_Tag_ReadsPerCell.csv", row.names = 1) #if the file is not in the working directory, type the path to it as an argument to the function
HTO_table_main_dataset<-HTO_table_main_dataset[,1:4]
HTO_table_main_dataset<-t(HTO_table_main_dataset)

# Preparing the dataset and running HTODemux
Main_dataset <- NormalizeData(Main_dataset)
Main_dataset <- FindVariableFeatures(Main_dataset, selection.method = "mean.var.plot")
Main_dataset <- ScaleData(Main_dataset, features = VariableFeatures(Main_dataset))
Main_dataset[["HTO"]]<- CreateAssayObject(counts = HTO_table_main_dataset)
Main_dataset <- NormalizeData(Main_dataset, assay = "HTO", normalization.method = "CLR")
Main_dataset <- HTODemux(Main_dataset, assay = "HTO", positive.quantile = 0.99)

#Renaming the sample tags
HTODemux_tags<-c()
for(i in 1:length(Main_dataset$orig.ident)){
  if(Main_dataset$hash.ID[i]=="SampleTag01-flex-Read-Count"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_1")
  }
  else if(Main_dataset$hash.ID[i]=="SampleTag02-flex-Read-Count"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_2")
  }
  else if(Main_dataset$hash.ID[i]=="SampleTag03-flex-Read-Count"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_3")
  }
  else if(Main_dataset$hash.ID[i]=="SampleTag04-flex-Read-Count"){
    HTODemux_tags<-c(HTODemux_tags, "Tag_4")
  }
  else if(Main_dataset$hash.ID[i]=="Doublet"){
    HTODemux_tags<-c(HTODemux_tags, "Multiplet")
  }
  else{
    HTODemux_tags<-c(HTODemux_tags, "Undetermined")
  }
}
Main_dataset[["HTODemux_tags"]]<-HTODemux_tags

#BFF

#Changing the names of cell barcodes

#Setting up the variables with the sample tag expression tables from the other datasets
cmo<-read.csv("CMO_Brain_tag_counts.csv")
hto<-read.csv("HTO_sample_tag_table.csv")
unique_barcodes<-cmo$X
additional_barcodes<-hto$X

#Selecting unique barcodes for the main dataset
for(i in 1:1000){
  x<-length(unique(str_detect(cmo$X, hto$X[i])))
  if(x == 1){
    unique_barcodes<-c(unique_barcodes, additional_barcodes[i])
  }
  if(length(unique_barcodes) == 32000){
    break
  }
}

#Creating a final version of the sample tag expression table for BFF algorithms
HTO_table_main_dataset<-read.csv("Exact-32000-Female-3-months-exercise-JAL_Sample_Tag_ReadsPerCell.csv")
HTO_table_main_dataset$cell<-unique_barcodes
write.csv(HTO_table_main_dataset, file="Main_dataset_sample_tag_table_renamed.csv")
HTO_table_main_dataset<-read.csv("Main_dataset_sample_tag_table_renamed.csv", row.names = "cell")
HTO_table_main_dataset$X<-c()
HTO_table_main_dataset<-HTO_table_main_dataset[,1:4]
expression_matrix<-HTO_table_main_dataset
HTO_table_main_dataset<-t(HTO_table_main_dataset)

#Running BFF
calls <- GenerateCellHashingCalls(barcodeMatrix = HTO_table_main_dataset, methods = c('bff_raw', 'bff_cluster'))
write.csv(calls, file="Main_dataset_BFF_calls.csv")
calls<-read.csv("Main_dataset_BFF_calls.csv", row.names = "cellbarcodes")
calls$X<-c()
#Ordering tags of the cells in a proper order and adding them to the dataset
#BFF cluster
final_ordered_tags<-c()
for(i in 1:32000){
  final_ordered_tags<-c(final_ordered_tags, calls[rownames(expression_matrix[i,]),][2])
}
final_vector_tags<-c()
for(i in 1:32000){
  final_vector_tags<-c(final_vector_tags, final_ordered_tags[[i]])
}

#Renaming the sample tags
BFF_cluster_tags<-c()
for(i in 1:length(Main_dataset$orig.ident)){
  if(final_vector_tags[i]=="SampleTag01_flex_Read_Count"){
    BFF_cluster_tags<-c(BFF_cluster_tags, "Tag_1")
  }
  else if(final_vector_tags[i]=="SampleTag02_flex_Read_Count"){
    BFF_cluster_tags<-c(BFF_cluster_tags, "Tag_2")
  }
  else if(final_vector_tags[i]=="SampleTag03_flex_Read_Count"){
    BFF_cluster_tags<-c(BFF_cluster_tags, "Tag_3")
  }
  else if(final_vector_tags[i]=="SampleTag04_flex_Read_Count"){
    BFF_cluster_tags<-c(BFF_cluster_tags, "Tag_4")
  }
  else if(final_vector_tags[i]=="Doublet"){
    BFF_cluster_tags<-c(BFF_cluster_tags, "Multiplet")
  }
  else{
    BFF_cluster_tags<-c(BFF_cluster_tags, "Undetermined")
  }
}
Main_dataset[["BFF_cluster_tags"]]<-BFF_cluster_tags

#BFF raw
#Ordering sample tags in a proper order
final_ordered_tags<-c()
for(i in 1:32000){
  final_ordered_tags<-c(final_ordered_tags, calls[rownames(expression_matrix[i,]),][1])
}
final_vector_tags<-c()
for(i in 1:32000){
  final_vector_tags<-c(final_vector_tags, final_ordered_tags[[i]])
}

#Renaming the sample tags
BFF_raw_tags<-c()
for(i in 1:length(Main_dataset$orig.ident)){
  if(final_vector_tags[i]=="SampleTag01_flex_Read_Count"){
     BFF_raw_tags<-c( BFF_raw_tags, "Tag_1")
  }
  else if(final_vector_tags[i]=="SampleTag02_flex_Read_Count"){
     BFF_raw_tags<-c( BFF_raw_tags, "Tag_2")
  }
  else if(final_vector_tags[i]=="SampleTag03_flex_Read_Count"){
     BFF_raw_tags<-c( BFF_raw_tags, "Tag_3")
  }
  else if(final_vector_tags[i]=="SampleTag04_flex_Read_Count"){
     BFF_raw_tags<-c( BFF_raw_tags, "Tag_4")
  }
  else if(final_vector_tags[i]=="Doublet"){
     BFF_raw_tags<-c( BFF_raw_tags, "Multiplet")
  }
  else{
     BFF_raw_tags<-c( BFF_raw_tags, "Undetermined")
  }
}
Main_dataset[["BFF_raw_tags"]]<-BFF_raw_tags

#hashedDrops

#Running hashedDrops
hashedDrops_demultiplexing<-hashedDrops(HTO_table_main_dataset)

#Creating a vector with hashedDrops tags
hashedDrops_tags<-c()
for(i in 1:length(hashedDrops_demultiplexing$Total)){
  if(hashedDrops_demultiplexing$Confident[i] == TRUE & hashedDrops_demultiplexing$Doublet[i] == FALSE){
    hashedDrops_tags<-c(hashedDrops_tags, hashedDrops_demultiplexing$Best[i])
  }
  else if(hashedDrops_demultiplexing$Confident[i] == FALSE & hashedDrops_demultiplexing$Doublet[i] == TRUE){
    hashedDrops_tags<-c(hashedDrops_tags, "Multiplet")
  }
  else if(hashedDrops_demultiplexing$Confident[i] == FALSE & hashedDrops_demultiplexing$Doublet[i] == FALSE){
    hashedDrops_tags<-c(hashedDrops_tags, "Undetermined")
  }
}

#Renaming sample tags
final_hashedDrops_tags<-c()
for(i in 1:length(hashedDrops_tags)){
  if(hashedDrops_tags[i]=="1"){
    final_hashedDrops_tags<-c(final_hashedDrops_tags, "Tag_1")
  }
  else if(hashedDrops_tags[i]=="2"){
    final_hashedDrops_tags<-c(final_hashedDrops_tags, "Tag_2")
  }
  else if(hashedDrops_tags[i]=="3"){
    final_hashedDrops_tags<-c(final_hashedDrops_tags, "Tag_3")
  }
  else if(hashedDrops_tags[i]=="4"){
    final_hashedDrops_tags<-c(final_hashedDrops_tags, "Tag_4")
  }
  else if(hashedDrops_tags[i]=="Multiplet"){
    final_hashedDrops_tags<-c(final_hashedDrops_tags, "Multiplet")
  }
  else if(hashedDrops_tags[i]=="Undetermined"){
    final_hashedDrops_tags<-c(final_hashedDrops_tags, "Undetermined")
  }
}

Main_dataset[["hashedDrops_tags"]]<-final_hashedDrops_tags

#RescueTag
Main_dataset <- ScoringTags("Exact-32000-Female-3-months-exercise-JAL_Sample_Tag_ReadsPerCell.csv", 1, 4, 1, 4, Main_dataset)
Main_dataset <- PreQCFilter(Main_dataset)
Main_dataset <- NormToDepth(Main_dataset)
df_for_plot<-data.frame('First_best'=Main_dataset$NTD_first_best)
df_for_plot %>% ggplot(aes(x=First_best))+geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+expand_limits(x=c(0,1))
Main_dataset <- FilterLowSeqDepth(Main_dataset) # for the distribution write either "normal" or "bimodal"
Main_dataset <- NormTagQCParams(Main_dataset)
Main_dataset[["Ratio"]]<-Main_dataset$Postnorm_delta/Main_dataset$Postnorm_first_best
Main_dataset<-FinalTagging(Main_dataset)

#Comparison of HTODemux and RescueTag
HTODemux_tags<-Main_dataset$HTODemux_tags
Main_dataset<-CompareRescueTag(Main_dataset, HTODemux_tags, "HTODemux")

#Scatter plots for main dataset, RescueTag vs HTODemux
ComparisonScatterPlots(Main_dataset, "HTODemux")

#Displaying ridge plots, main dataset, RescueTag vs HTODemux
ComparisonRidgePlotsBest(Main_dataset, HTODemux_tags)
ComparisonRidgePlotsDelta(Main_dataset, HTODemux_tags)
ComparisonPercentages(Main_dataset, HTODemux_tags)

#Comparison plots for main dataset, RescueTag vs BFF cluster
BFF_cluster_tags<-Main_dataset$BFF_cluster_tags
Main_dataset<-CompareRescueTag(Main_dataset, BFF_cluster_tags, "BFF_cluster")
ComparisonScatterPlots(Main_dataset, "BFF_cluster")
ComparisonRidgePlotsBest(Main_dataset, BFF_cluster_tags)
ComparisonRidgePlotsDelta(Main_dataset, BFF_cluster_tags)
ComparisonPercentages(Main_dataset, BFF_cluster_tags)

#Comparison plots for main dataset, RescueTag vs BFF raw
BFF_raw_tags<-Main_dataset$BFF_raw_tags
Main_dataset<-CompareRescueTag(Main_dataset, BFF_raw_tags, "BFF_raw")
ComparisonScatterPlots(Main_dataset, "BFF_raw")
ComparisonRidgePlotsBest(Main_dataset, BFF_raw_tags)
ComparisonRidgePlotsDelta(Main_dataset, BFF_raw_tags)
ComparisonPercentages(Main_dataset, BFF_raw_tags)

#Comparison plots for main dataset, RescueTag vs hashedDrops
hashedDrops_tags<-Main_dataset$hashedDrops_tags
Main_dataset<-CompareRescueTag(Main_dataset, hashedDrops_tags, "hashedDrops")
ComparisonScatterPlots(Main_dataset, "hashedDrops")
ComparisonRidgePlotsBest_hashedDrops(Main_dataset, hashedDrops_tags)
ComparisonRidgePlotsDelta_hashedDrops(Main_dataset, hashedDrops_tags)
ComparisonPercentages_hashedDrops(Main_dataset, hashedDrops_tags)

###########################################################
#                                                         #
#           RescueCluster vs Classical analysis           #
#                                                         #
###########################################################

Main_dataset<-subset(Main_dataset, subset=Final_tags != "Multiplet" | Final_tags != "Undetermined")
Main_dataset[["percent.mt"]] <- PercentageFeatureSet(Main_dataset, pattern="^mt") #calculate percent mt
Idents(Main_dataset) <- "orig.ident"
VlnPlot(Main_dataset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=1)

Main_dataset <- Main_dataset %>% NormalizeData() %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% ScaleData() %>% SCTransform(vars.to.regress = c("percent.mt"))
Main_dataset <- RunPCA(Main_dataset, assay = "SCT")
ElbowPlot(Main_dataset)
Main_dataset <- RunHarmony(Main_dataset, group.by.vars = c("Final_tags"), reduction = "pca", assay.use = "SCT", reduction.save = "harmony", dims.use=1:6)
Main_dataset <- RunUMAP(Main_dataset, reduction = "harmony", assay = "SCT", dims = 1:6) 
Main_dataset <- FindNeighbors(object = Main_dataset, reduction = "harmony", dims = 1:6) 
Main_dataset <- FindClusters(Main_dataset, resolution = 0.5)

RefinedDimPlot(Main_dataset, 30)
VlnPlot(Main_dataset, features = 'percent.mt', group.by = 'seurat_clusters') # displays MGPC distribution in every cluster
VlnPlot(Main_dataset, features = 'nFeature_RNA', group.by = 'seurat_clusters') # displays number of genes per cell distribution in every cluster
ClusterQC(Main_dataset, "seurat_clusters")

#Classical analysis

Classical_main_dataset<-readRDS("Female-3-months-exercise-JAL_Seurat.rds")
Classical_main_dataset<- PreQCFilter(Classical_main_dataset)
Classical_main_dataset[["percent.mt"]] <- PercentageFeatureSet(Classical_main_dataset, pattern="^mt") #calculate percent mt
VlnPlot(Classical_main_dataset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.1)
Classical_main_dataset <- subset(Classical_main_dataset, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 30)
Classical_main_dataset <- subset(Classical_main_dataset, subset = Sample_Name !="Multiplet")
Classical_main_dataset <- NormalizeData(Classical_main_dataset) 
Classical_main_dataset <- FindVariableFeatures(Classical_main_dataset, selection.method="vst", nfeatures=2000) 
Classical_main_dataset <- ScaleData(Classical_main_dataset) 
Classical_main_dataset <- RunPCA(Classical_main_dataset, features=VariableFeatures(object=Classical_main_dataset)) 
ElbowPlot(Classical_main_dataset)
Classical_main_dataset <- FindNeighbors(Classical_main_dataset, dims = 1:6)
Classical_main_dataset <- FindClusters(Classical_main_dataset, resolution = 0.5)
Classical_main_dataset <- RunUMAP(Classical_main_dataset, dims = 1:6)
DimPlot(Classical_main_dataset, reduction = "umap", label=T)+NoLegend()

###################################################################
#                                                                 #
#           Annotation of clusters with AutoClusterType           #
#                                                                 #
###################################################################

#RescueCluster strategy
querry <- GetAssayData(Main_dataset, slot = "data")
querryMatrix <- as.matrix(querry)
mca_result <- scMCA(scdata = querryMatrix, numbers_plot = 3)
Main_dataset[["MCA"]] <- mca_result$scMCA
Main_dataset <- AutoClusterType("Reference_annotations_scMCA.csv", mca_result, Main_dataset)
Idents(Main_dataset) <- "cell.ID"
RefinedDimPlot(Main_dataset, 30) 
ClusterQC(Main_dataset, "cell.ID")

#Classical strategy
querry <- GetAssayData(Classical_main_dataset, slot = "data")
querryMatrix <- as.matrix(querry)
mca_result <- scMCA(scdata = querryMatrix, numbers_plot = 3)
Classical_main_dataset[["MCA"]] <- mca_result$scMCA
Classical_main_dataset <- AutoClusterType("Reference_annotations_scMCA.csv", mca_result, Classical_main_dataset)
Idents(Classical_main_dataset) <- "cell.ID"
RefinedDimPlot(Classical_main_dataset, 30) 

##################################################
#                                                #
#           RescueCluster benchmarking           #
#                                                #
##################################################

E12_5_1<-Read10X("Mouse_heart_dataset/Raw_data/E12_5_1")
E12_5_2<-Read10X("Mouse_heart_dataset/Raw_data/E12_5_2")
E14_5_1<-Read10X("Mouse_heart_dataset/Raw_data/E14_5_1")
E14_5_2<-Read10X("Mouse_heart_dataset/Raw_data/E14_5_2")
E16_5_1<-Read10X("Mouse_heart_dataset/Raw_data/E16_5_1")
E16_5_2<-Read10X("Mouse_heart_dataset/Raw_data/E16_5_2")
p3_1_2<-Read10X("Mouse_heart_dataset/Raw_data/p3_1_2")
E12_5_1<-CreateSeuratObject(E12_5_1, project = "E12_5_1")
E12_5_2<-CreateSeuratObject(E12_5_2, project = "E12_5_2")
E14_5_1<-CreateSeuratObject(E14_5_1, project = "E14_5_1")
E14_5_2<-CreateSeuratObject(E14_5_2, project = "E14_5_2")
E16_5_1<-CreateSeuratObject(E16_5_1, project = "E16_5_1")
E16_5_2<-CreateSeuratObject(E16_5_2, project = "E16_5_2")
p3_1_2<-CreateSeuratObject(p3_1_2, project = "p3_1_2")
Mouse_heart_dataset<-merge(E12_5_1, y=c(E12_5_2,E14_5_1,E14_5_2,E16_5_1,E16_5_2,p3_1_2), add.cell.ids=c("E12_5_1", "E12_5_2", "E14_5_1", "E14_5_2", "E16_5_1", "E16_5_2", "p3_1_2"))
Mouse_heart_dataset[["Sample_Tags"]]<-Mouse_heart_dataset$orig.ident
Mouse_heart_dataset[["orig.ident"]]<-rep("heart",length(Mouse_heart_dataset$orig.ident))

Mouse_heart_dataset <- PreQCFilter(Mouse_heart_dataset)
Mouse_heart_dataset[["percent.mt"]] <- PercentageFeatureSet(Mouse_heart_dataset, pattern="^mt") 
Idents(Mouse_heart_dataset) <- "orig.ident"
VlnPlot(Mouse_heart_dataset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=1)

Mouse_heart_dataset <- Mouse_heart_dataset %>% NormalizeData() %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% ScaleData() %>% SCTransform(vars.to.regress = c("percent.mt"))
Mouse_heart_dataset <- RunPCA(Mouse_heart_dataset, assay = "SCT")
ElbowPlot(Mouse_heart_dataset)
Mouse_heart_dataset <- RunHarmony(Mouse_heart_dataset, group.by.vars = c("Sample_Tags"), reduction = "pca", assay.use = "SCT", reduction.save = "harmony", dims.use=1:7) 
Mouse_heart_dataset <- RunUMAP(Mouse_heart_dataset, reduction = "harmony", assay = "SCT", dims = 1:7) 
Mouse_heart_dataset <- FindNeighbors(object = Mouse_heart_dataset, reduction = "harmony", dims = 1:7) 
Mouse_heart_dataset <- FindClusters(Mouse_heart_dataset, resolution = 0.5)

RefinedDimPlot(Mouse_heart_dataset, 30)
VlnPlot(Mouse_heart_dataset, features = 'percent.mt', group.by = 'seurat_clusters') 
VlnPlot(Mouse_heart_dataset, features = 'nFeature_RNA', group.by = 'seurat_clusters') 
ClusterQC(Mouse_heart_dataset, "seurat_clusters")

querry <- GetAssayData(Mouse_heart_dataset, slot = "data")
querryMatrix <- as.matrix(querry)
mca_result <- scMCA(scdata = querryMatrix, numbers_plot = 3)
Mouse_heart_dataset[["MCA"]] <- mca_result$scMCA
Mouse_heart_dataset <- AutoClusterType("Reference_annotations_scMCA.csv", mca_result, Mouse_heart_dataset)
Idents(Mouse_heart_dataset) <- "cell.ID"
RefinedDimPlot(Mouse_heart_dataset, 30) 
ClusterQC(Mouse_heart_dataset, "cell.ID")

##############################################
#                                            #
#           DEGs of cardiomyocytes           #
#                                            #
##############################################

cardiomyocytes<-subset(Mouse_heart_dataset, subset=cell.ID == c("Cardiomyocyte_Neonatal-heart/Stromal_cell_Neonatal-heart", "Cardiomyocyte_Neonatal-heart_1", "Cardiomyocyte_Neonatal-heart_2", "Cardiomyocyte_Neonatal-heart_3", "Cardiomyocyte_Neonatal-heart_4", "Cardiomyocyte_Neonatal-heart_5", "Cardiomyocyte_Neonatal-heart_6", "Cardiomyocyte_Neonatal-heart_7", "Cardiomyocyte_Neonatal-heart_8", "Cardiomyocyte_Neonatal-heart_9", "Cardiomyocyte_Neonatal-heart_10", "Cardiomyocyte_Neonatal-heart_11","Cardiomyocyte_Neonatal-heart_12", "Cardiomyocyte_Neonatal-heart_13", "Cardiomyocyte_Neonatal-heart_14","Cardiomyocyte_Neonatal-heart_15", "Cardiomyocyte_Neonatal-heart_16"))
cardio_clust<-unique(cardiomyocytes$seurat_clusters)
for (i in cardio_clust) {
  test_cardio<-cardiomyocytes
  Idents(test_cardio)<-"seurat_clusters"
  if(i==2){
    test_cardio<-RenameIdents(test_cardio, "6"="other","5"="other","0"="other","13"="other", "10"="other","16"="other","1"="other","15"="other","12"="other","11"="other","4"="other","3"="other","18"="other","19"="other","21"="other","20"="other")
    test_cardio[["seurat_clusters"]]<-Idents(test_cardio)
  }
  else if (i==6){
    test_cardio<-RenameIdents(test_cardio, "2"="other","5"="other","0"="other","13"="other", "10"="other","16"="other","1"="other","15"="other","12"="other","11"="other","4"="other","3"="other","18"="other","19"="other","21"="other","20"="other")
    test_cardio[["seurat_clusters"]]<-Idents(test_cardio)
  }
  else if (i==5){
    test_cardio<-RenameIdents(test_cardio, "2"="other","6"="other","0"="other","13"="other", "10"="other","16"="other","1"="other","15"="other","12"="other","11"="other","4"="other","3"="other","18"="other","19"="other","21"="other","20"="other")
    test_cardio[["seurat_clusters"]]<-Idents(test_cardio)
  }
  else if (i==0){
    test_cardio<-RenameIdents(test_cardio, "2"="other","6"="other","5"="other","13"="other", "10"="other","16"="other","1"="other","15"="other","12"="other","11"="other","4"="other","3"="other","18"="other","19"="other","21"="other","20"="other")
    test_cardio[["seurat_clusters"]]<-Idents(test_cardio)
  }
  else if (i==13){
    test_cardio<-RenameIdents(test_cardio, "2"="other","6"="other","5"="other","0"="other", "10"="other","16"="other","1"="other","15"="other","12"="other","11"="other","4"="other","3"="other","18"="other","19"="other","21"="other","20"="other")
    test_cardio[["seurat_clusters"]]<-Idents(test_cardio)
  }
  else if (i==10){
    test_cardio<-RenameIdents(test_cardio, "2"="other","6"="other","5"="other","0"="other", "13"="other","16"="other","1"="other","15"="other","12"="other","11"="other","4"="other","3"="other","18"="other","19"="other","21"="other","20"="other")
    test_cardio[["seurat_clusters"]]<-Idents(test_cardio)
  }
  else if (i==16){
    test_cardio<-RenameIdents(test_cardio, "2"="other","6"="other","5"="other","0"="other","13"="other", "10"="other","1"="other","15"="other","12"="other","11"="other","4"="other","3"="other","18"="other","19"="other","21"="other","20"="other")
    test_cardio[["seurat_clusters"]]<-Idents(test_cardio)
  }
  else if (i==1){
    test_cardio<-RenameIdents(test_cardio, "2"="other","6"="other","5"="other","0"="other","13"="other", "10"="other","16"="other","15"="other","12"="other","11"="other","4"="other","3"="other","18"="other","19"="other","21"="other","20"="other")
    test_cardio[["seurat_clusters"]]<-Idents(test_cardio)
  }
  else if (i==15){
    test_cardio<-RenameIdents(test_cardio, "2"="other","6"="other","5"="other","0"="other","13"="other", "10"="other","16"="other","1"="other","12"="other","11"="other","4"="other","3"="other","18"="other","19"="other","21"="other","20"="other")
    test_cardio[["seurat_clusters"]]<-Idents(test_cardio)
  }
  else if (i==12){
    test_cardio<-RenameIdents(test_cardio, "2"="other","6"="other","5"="other","0"="other","13"="other", "10"="other","16"="other","1"="other","15"="other","11"="other","4"="other","3"="other","18"="other","19"="other","21"="other","20"="other")
    test_cardio[["seurat_clusters"]]<-Idents(test_cardio)
  }
  else if (i==11){
    test_cardio<-RenameIdents(test_cardio, "2"="other","6"="other","5"="other","0"="other","13"="other", "10"="other","16"="other","1"="other","15"="other","12"="other","4"="other","3"="other","18"="other","19"="other","21"="other","20"="other")
    test_cardio[["seurat_clusters"]]<-Idents(test_cardio)
  }
  else if (i==4){
    test_cardio<-RenameIdents(test_cardio, "2"="other","6"="other","5"="other","0"="other","13"="other", "10"="other","16"="other","1"="other","15"="other","12"="other","11"="other","3"="other","18"="other","19"="other","21"="other","20"="other")
    test_cardio[["seurat_clusters"]]<-Idents(test_cardio)
  }
  else if (i==3){
    test_cardio<-RenameIdents(test_cardio, "2"="other","6"="other","5"="other","0"="other","13"="other", "10"="other","16"="other","1"="other","15"="other","12"="other","11"="other","4"="other","18"="other","19"="other","21"="other","20"="other")
    test_cardio[["seurat_clusters"]]<-Idents(test_cardio)
  }
  else if (i==18){
    test_cardio<-RenameIdents(test_cardio, "2"="other","6"="other","5"="other","0"="other","13"="other", "10"="other","16"="other","1"="other","15"="other","12"="other","11"="other","4"="other","3"="other","19"="other","21"="other","20"="other")
    test_cardio[["seurat_clusters"]]<-Idents(test_cardio)
  }
  else if (i==19){
    test_cardio<-RenameIdents(test_cardio, "2"="other","6"="other","5"="other","0"="other","13"="other", "10"="other","16"="other","1"="other","15"="other","12"="other","11"="other","4"="other","3"="other","18"="other","21"="other","20"="other")
    test_cardio[["seurat_clusters"]]<-Idents(test_cardio)
  }
  else if (i==21){
    test_cardio<-RenameIdents(test_cardio, "2"="other","6"="other","5"="other","0"="other","13"="other", "10"="other","16"="other","1"="other","15"="other","12"="other","11"="other","4"="other","3"="other","18"="other","19"="other","20"="other")
    test_cardio[["seurat_clusters"]]<-Idents(test_cardio)
  }
  else if (i==20){
    test_cardio<-RenameIdents(test_cardio, "2"="other","6"="other","5"="other","0"="other","13"="other", "10"="other","16"="other","1"="other","15"="other","12"="other","11"="other","4"="other","3"="other","18"="other","19"="other","21"="other")
    test_cardio[["seurat_clusters"]]<-Idents(test_cardio)
  }
  if (length(unique(test_cardio$seurat_clusters)) == 2){
    Idents(test_cardio)<-'seurat_clusters'
    x<-FindMarkers(test_cardio, ident.1 = cardio_clust[which(cardio_clust==i)], ident.2 = "other")
    if (exists('x') & nrow(x)>0){
      write.csv(x, file = paste(cardio_clust[which(cardio_clust==i)], '_', "vs all", '_DEG.csv', sep=''))
    }
  }
}

#Exploration of top 30 DEGs for every cardiomyocyte cluster

heart_degs<-list.files(path="Mouse heart DEGs/")
cardio_clusters<-sort(unique(cardiomyocytes$seurat_clusters))
for(i in 1:length(heart_degs)){
  x<-read.csv(paste("Mouse heart DEGs/", heart_degs[i], sep=""))
  FeaturePlot(cardiomyocytes, features = x$X[1:12])
  ggsave(filename = paste(as.character(cardio_clusters[i]),"_vs_all_1-12.pdf"), height = 6.94, width = 9.57)
  FeaturePlot(cardiomyocytes, features = x$X[13:24])
  ggsave(filename = paste(as.character(cardio_clusters[i]),"_vs_all_13-24.pdf"), height = 6.94, width = 9.57)
  FeaturePlot(cardiomyocytes, features = x$X[25:30])
  ggsave(filename = paste(as.character(cardio_clusters[i]),"_vs_all_25-30.pdf"), height = 6.94, width = 9.57)
}

#FeaturePlots with interesting differentially expressed genes in cardiomyocyte clusters

#Clusters with mature cardiomyocytes
FeaturePlot(cardiomyocytes, features = c("Cox6a2", "Atp5g3", "Ckm", "Actn2", "Pln", "Fabp3","Cox8b","Atp2a2","Mb"))

#Cluster with fibroblast contamination
FeaturePlot(cardiomyocytes, features = c("Col3a1", "Mfap4", "Col1a2", "Lum", "Dpt", "Tcf21", "Dcn"))

#Cluster with endothelial cell contamination
FeaturePlot(cardiomyocytes, features = c("Igfbp4", "Tm4sf1", "Egfl7", "Ramp2", "Ecscr", "Eng", "Gngt2", "Calcrl"))

#Clusters with atrial cardiomyocytes
FeaturePlot(cardiomyocytes, features = c("Pam", "Myl7", "Gpx3", "Gja5", "Epha4", "Dkk3", "Stard10", "Kcnj3", "Nppa", "Bmp10"))

#Clusters with inflamed/unhealthy cardiomyocytes and/or some erythrocyte contamination
FeaturePlot(cardiomyocytes, features = c("Hba-a1","Snca", "Isg20", "Alas2"))
