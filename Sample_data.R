tmp<-read.csv2('C:/tmp/Production.csv')
tmp$default<-((tmp$event_family+tmp$event_health+tmp$event_home+tmp$event_work)>0)*1
tmp$default
dim(tmp)
tmp2<-tmp[(tmp$cid/2-round(tmp$cid/2-0.01))==0,c(0:18,23)]
write.csv2(tmp2, 'C:/tmp/ProductionS.csv')

tmp<-read.csv2('C:/tmp/Ksiazka/abt_sam_beh_train.csv')
View(tmp[1:10,])
tmp2<-tmp[(tmp$cid/2-round(tmp$cid/2-0.01))==0,
          c(names(tmp)[grep('app_', names(tmp))],
            names(tmp)[grep('act_cus', names(tmp))],
            'default_cus12')]

tmp2$default_flag<-(tmp2$default_cus12==1)*1

write.csv2(tmp2, 'C:/tmp/Ksiazka/abt_sam_beh_trainS.csv')


# names(tmp)[grep('act_cus', names(tmp))]
