
library(dplyr)

if(F){
  insert_size <-insert_size <- read.table("/home/user/data/lit/project/6mA/mapping/PC10_LY_MCC-M1.bedpe.sorted.bed")
  inflection_width <- read.table("/home/user/data/lit/project/6mA/iNPS_output/PC10_LY_MCC-M1_Gathering.like_bed",skip = 8,sep = '\t')
  insert_size %>% mutate(insert_size=V3-V2)
  mean(insert_size$V3-insert_size$V2)
  mean(inflection_width$V5)
  hist(insert_size$V3-insert_size$V2)
  hist(inflection_width$V5)
}

##### method 1 #####

visual_method_1 <- function(binned_level_dir){
  
  binned_level <- read.table(binned_level_dir)
  binned_level$V1 %>% str_split_fixed(.,"\\.",2) %>% data.frame() %>% mutate(binned_level) %>% 
    set_names("sample","bin","sample.bin","mA","A","FreqSum") -> binned_level.input
  
  binned_level.input %>%
    ggplot(mapping = aes(x = factor(bin,levels = c("nucleosome.bin_1","nucleosome.bin_2","nucleosome.bin_3",
                                                   "nucleosome.bin_4","nucleosome.bin_5",
                                                   "linker.bin_1","linker.bin_2","linker.bin_3","linker.bin_4","linker.bin_5")), 
                         y = mA,color=sample,group=sample)) + 
    geom_line()+
    geom_point()+
    theme(axis.text.x = element_text(hjust = 1, angle = 45),
          axis.title.x = element_blank())+
    ylab("6mA")
}


# binned_level <- read.table("/home/user/data2/lit/project/6mA/binned.level.txt")
# binned_level <- read.table("/home/user/data2/lit/project/6mA/binned.level.list_1.txt")
# binned_level <- read.table("/home/user/data2/lit/project/6mA/binned.level.list_2.txt")
binned_level <- read.table("/home/user/data/lit/project/6mA/feature/danpos/method_1/binned.level.txt")
binned_level <- read.table("/home/user/data/lit/project/6mA/feature/iNPS/method_1/binned.level.txt")


binned_level$V1 %>% str_split_fixed(.,"\\.",2) %>% data.frame() %>% mutate(binned_level) %>% 
  set_names("sample","bin","sample.bin","mA","A","FreqSum") -> binned_level.input

binned_level.input %>%
  ggplot(mapping = aes(x = factor(bin,levels = c("nucleosome.bin_1","nucleosome.bin_2","nucleosome.bin_3",
                                                 "nucleosome.bin_4","nucleosome.bin_5",
                                                 "linker.bin_1","linker.bin_2","linker.bin_3","linker.bin_4","linker.bin_5")), 
                       y = mA,color=sample,group=sample)) + 
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(hjust = 1, angle = 45),
        axis.title.x = element_blank())+
  ylab("6mA")

p <- binned_level.input %>%
ggplot(mapping = aes(x = factor(bin,levels = c("nucleosome.bin_1","nucleosome.bin_2","nucleosome.bin_3",
                                               "nucleosome.bin_4","nucleosome.bin_5",
                                               "linker.bin_1","linker.bin_2","linker.bin_3","linker.bin_4","linker.bin_5")), 
                     y = mA/A,color=sample,group=sample)) + 
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(hjust = 1, angle = 45),
        axis.title.x = element_blank())+
        ylab("6mA/A")+
  theme(legend.key.size = unit(0.2, "pt"))


ggsave(p,filename = "/home/user/data/lit/project/6mA/feature/danpos/method_1/binned_level_6mA.png",width = 5,height = 3)


p <- binned_level.input %>%
  ggplot(mapping = aes(x = factor(bin,levels = c("nucleosome.bin_1","nucleosome.bin_2","nucleosome.bin_3",
                                                 "nucleosome.bin_4","nucleosome.bin_5",
                                                 "linker.bin_1","linker.bin_2","linker.bin_3","linker.bin_4","linker.bin_5")), 
                       y = FreqSum/A,color=sample,group=sample)) + 
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(hjust = 1, angle = 45),
        axis.title.x = element_blank())+
  ylab("FreqSum/A")+
  theme(legend.key.size = unit(0.2, "pt"))

ggsave(p,filename = "binned_level_FreqSum.png",width = 5,height = 3)
 

over_all_level <- read.table("/home/user/data2/lit/project/6mA/level.txt")

over_all_level$V1 %>% str_split_fixed(.,"\\.",2) %>% data.frame() %>% mutate(over_all_level) %>% 
  set_names("sample","bin","sample.bin","mA","A","FreqSum") -> over_all_level.input


over_all_level.input %>%
  ggplot(mapping = aes(x = factor(bin,levels = c("nucleosome","linker")), 
                       y = mA,color=sample,group=sample)) + 
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(hjust = 1, angle = 45),
        axis.title.x = element_blank())+
  ylab("6mA")

p <- over_all_level.input %>%
  ggplot(mapping = aes(x = factor(bin,levels = c("nucleosome","linker")), 
                       y = mA/A,color=sample,group=sample)) + 
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(hjust = 1, angle = 45),
        axis.title.x = element_blank())+
  ylab("6mA/A")+
  # theme(legend.title = element_text(size = 10), 
  #       legend.text = element_text(size = 10))+
  theme(legend.key.size = unit(0.2, "pt"))

ggsave(p,filename = "over_all_level_6mA.png",width = 5,height = 3)


p <- over_all_level.input %>%
  ggplot(mapping = aes(x = factor(bin,levels = c("nucleosome","linker")), 
                       y = FreqSum/A,color=sample,group=sample)) + 
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(hjust = 1, angle = 45),
        axis.title.x = element_blank())+
  ylab("FreqSum/A")+
  theme(legend.key.size = unit(0.2, "pt"))


ggsave(p,filename = "over_all_level_FreqSum.png",width = 5,height = 3)

##### method 1 #####



##### method 2 #####

binned_level_nucleosome <- read.table("/home/user/data2/lit/project/6mA/binned.level.ten_bins.txt")
binned_level_linker <- read.table("/home/user/data2/lit/project/6mA/binned.level.n.txt")
binned_level <- rbind(binned_level_nucleosome,binned_level_linker)

binned_level$V1 %>% str_split_fixed(.,"\\.",2) %>% data.frame() %>% mutate(binned_level) %>% 
  set_names("sample","bin","sample.bin","mA","A","FreqSum") -> binned_level.input

name_boom <- function(fixname,seq){
  name_all=NULL
  for (i in seq){
    name=paste0(fixname,i)
    name_all=c(name_all,name)
  }
  return(name_all)
}

binned_level.input %>%
  ggplot(mapping = aes(x = factor(bin,levels = c( name_boom("linker.bin_left_",seq(1,5)),
                                                           name_boom("nucleosome.bin_",seq(1,10)),
                                                      name_boom("linker.bin_right_",seq(1,5)))),
                       y = mA,color=sample,group=sample)) + 
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(hjust = 1, angle = 45),
        axis.title.x = element_blank())+
  ylab("6mA")

p <- binned_level.input %>%
  ggplot(mapping = aes(x = factor(bin,levels = c( name_boom("linker.bin_left_",seq(1,5)),
                                                  name_boom("nucleosome.bin_",seq(1,10)),
                                                  name_boom("linker.bin_right_",seq(1,5)))),
                       y = mA/A,color=sample,group=sample)) + 
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(hjust = 1, angle = 45),
        axis.title.x = element_blank())+
  ylab("6mA/A")+
  theme(legend.key.size = unit(0.2, "pt"))

ggsave(p,filename = "binned_level_6mA.png",width = 5,height = 3)


p <- binned_level.input %>%
  ggplot(mapping = aes(x = factor(bin,levels = c( name_boom("linker.bin_left_",seq(1,5)),
                                                  name_boom("nucleosome.bin_",seq(1,10)),
                                                  name_boom("linker.bin_right_",seq(1,5)))),
                       y = FreqSum/A,color=sample,group=sample)) + 
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(hjust = 1, angle = 45),
        axis.title.x = element_blank())+
  ylab("FreqSum/A")+
  theme(legend.key.size = unit(0.2, "pt"))

ggsave(p,filename = "binned_level_FreqSum.png",width = 5,height = 3)


over_all_level <- read.table("/home/user/data2/lit/project/6mA/level.txt")

over_all_level$V1 %>% str_split_fixed(.,"\\.",2) %>% data.frame() %>% mutate(over_all_level) %>% 
  set_names("sample","bin","sample.bin","mA","A","FreqSum") -> over_all_level.input


over_all_level.input %>%
  ggplot(mapping = aes(x = factor(bin,levels = c("nucleosome","linker")), 
                       y = mA,color=sample,group=sample)) + 
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(hjust = 1, angle = 45),
        axis.title.x = element_blank())+
  ylab("6mA")

p <- over_all_level.input %>%
  ggplot(mapping = aes(x = factor(bin,levels = c("nucleosome","linker")), 
                       y = mA/A,color=sample,group=sample)) + 
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(hjust = 1, angle = 45),
        axis.title.x = element_blank())+
  ylab("6mA/A")+
  # theme(legend.title = element_text(size = 10), 
  #       legend.text = element_text(size = 10))+
  theme(legend.key.size = unit(0.2, "pt"))

ggsave(p,filename = "over_all_level_6mA.png",width = 5,height = 3)


p <- over_all_level.input %>%
  ggplot(mapping = aes(x = factor(bin,levels = c("nucleosome","linker")), 
                       y = FreqSum/A,color=sample,group=sample)) + 
  geom_line()+
  geom_point()+
  theme(axis.text.x = element_text(hjust = 1, angle = 45),
        axis.title.x = element_blank())+
  ylab("FreqSum/A")+
  theme(legend.key.size = unit(0.2, "pt"))


ggsave(p,filename = "over_all_level_FreqSum.png",width = 5,height = 3)