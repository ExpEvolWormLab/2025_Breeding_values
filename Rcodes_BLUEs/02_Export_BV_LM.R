rm(list=ls())
gc()
NGM_env=new.env();NaCl_env=new.env()
load("RData_BLUEs/NGM_LMER_BLUEs.RData",envir=NGM_env)
load("RData_BLUEs/NaCl_LMER_BLUEs.RData",envir=NaCl_env)

final_BLUE=NULL
for(traitID in NGM_env$vect_P_traits){
  final_BLUE.tp <- data.frame(traitID=traitID,pop_label=NGM_env$all_traits_corrected$pop_label,BLUE=NGM_env$all_traits_corrected[,traitID])
  if(traitID=="area.F"){
    final_BLUE.tp$traitID="area"
  }
  final_BLUE=rbind(final_BLUE,final_BLUE.tp)
}

write.table(final_BLUE,file="BLUES_NGM.txt",sep='\t',quote=FALSE,row.names=FALSE)


boxplot(final_BLUE$BLUE~final_BLUE$traitID)
boxplot(final_BLUE$BLUE[final_BLUE$traitID=="area"]~tstrsplit(final_BLUE$pop_label[final_BLUE$traitID=="area"],"L")[[1]],las=2,xlab="")
boxplot(final_BLUE$BLUE[final_BLUE$traitID=="area"]~substring(final_BLUE$pop_label[final_BLUE$traitID=="area"],1,3),las=2,xlab="")
boxplot(final_BLUE$BLUE[final_BLUE$traitID=="FS"]~substring(final_BLUE$pop_label[final_BLUE$traitID=="FS"],1,3),las=2,xlab="")

### Export the breeding values for each trait - NaCl ###
final_BLUE=NULL
for(traitID in NaCl_env$vect_P_traits){
  final_BLUE.tp <- data.frame(traitID=traitID,pop_label=NaCl_env$all_traits_corrected$pop_label,BLUE=NaCl_env$all_traits_corrected[,traitID])
  if(traitID=="area.F"){
    final_BLUE.tp$traitID="area"
  }
  final_BLUE=rbind(final_BLUE,final_BLUE.tp)
}

write.table(final_BLUE,file="BLUES_NaCl.txt",sep='\t',quote=FALSE,row.names=FALSE)

boxplot(final_BLUE$BLUE~final_BLUE$traitID)
boxplot(final_BLUE$BLUE[final_BLUE$traitID=="area"]~tstrsplit(final_BLUE$pop_label[final_BLUE$traitID=="area"],"L")[[1]],las=2,xlab="")

