library(ggplotify)
library(patchwork)
library(viridis)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

load("pooled_sim1kg_25ntimes.Rdata")
h=c("(10,10;10)","(30,30;30)","(60,60;60)","(10,60;60)")
delta=0.8
sim = 100
ntimes=25
metric = 'AUC'
exact=FALSE
n_method=25
plot_list = list()
j=1
titles = c('A: h=(10,10,10)','B: h=(30,30,30)','C: h=(60,60,60)','D: h=(10,60,60)')
ylabs = c("AUC", "", "AUC","")
mycolors <- c('gray67','lightskyblue1','lightskyblue2','lightskyblue3','mistyrose1','mistyrose2','mistyrose3','orange','orange1','orange2','orange3')
for(i in 1:4){
    DNN = read.csv(paste0("/Users/tiangu/OneDrive - Harvard University/Tian Gu-Shared/SynTL/JBI/DNN_h",h[i],"_4layers.csv"), header = F)
    DNN$AUC = ifelse(DNN$V1==0, NA, DNN$V1)
    
    SER = read.csv(paste0("/Users/tiangu/OneDrive - Harvard University/Tian Gu-Shared/SynTL/JBI/SER_h",h[i],".csv"), header = F)
    SER$AUC = ifelse(SER$V1==0, NA, SER$V1)
    
    name_data = paste0('h=', h[i],'_delta=', delta[j],'_exact=',exact, '_ntimes=', ntimes)
    dat = data.frame(pooled_out[[name_data]][[metric]][,c(1:4)], 
                     DNN$AUC, SER$AUC, 
                     pooled_out[[name_data]][[metric]][,5],
                     pooled_out[[name_data]][[metric]][,c(7:10)])
    colnames(dat) = c('Target-only','1','2','3','DNN','SER','Direct','M=0','M=1','M=2','M=3')
    dat <- dat[,c('Target-only','1','2','3','DNN','SER','Direct','M=0','M=1','M=2','M=3')]
    long <- melt(dat)
    long$group = c(rep('Target',sim), rep('Source',sim*3), rep('Others',sim*3), rep('COMMUTE', sim*4))
    long$group = factor(long$group, c('Target','Source','Others','COMMUTE'))
    #mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(ncol(dat))
    
    p = ggplot(long, aes(x=variable, y=value, fill=variable)) +
        scale_y_continuous(limits = c(0.5,0.9)) +
        geom_boxplot()+
        scale_fill_manual(values = mycolors) +
        ylab(ylabs[i]) + xlab('') +
        ggtitle(titles[i]) +
        theme(text = element_text(size=19),
              legend.key.size = unit(0.8, 'cm'),
              legend.text=element_text(size=14),
              axis.text.x=element_text(angle=15,hjust=1),
              legend.position = "none") +
        facet_grid(cols = vars(group), scales = "free", space = "free") #+
        #scale_fill_brewer()
    plot_list[[i]] = p
}
(plot_list[[1]] | plot_list[[2]]) / (plot_list[[3]] | plot_list[[4]])

ggsave(path = '/Users/tiangu/OneDrive - Harvard University/Tian Gu-Shared/SynTL/JBI', width = 15, height = 8, filename='sim1.png', dpi=200)






