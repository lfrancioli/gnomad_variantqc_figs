library('ggpubr')
library('scales')
library('grDevices')

theme_classic()
model_colors = c(
  RF="#1f77b4",
  VQSR="#ff7f0e"
)
model_names = c(
  RF="Random Forests",
  VQSR="Variant Quality Score Recalibration"
)

n_trios = c(
  exomes=4568,
  genomes=212
)

# True = snv, False = indel
cutoffs = c(
  'True'=90,
  'False'=82
)

plot_titles = c(
  'exomes True' = "Exomes SNVs",
  'exomes False' = "Exomes Indels",
  'genomes True' = "Genomes SNVs",
  'genomes False' = "Genomes Indels"
)

labels = c('a','b','c','d','e','f','g','h','i','j','k','l','m')

# Load data
autosomes = read_tsv('autosomes.tsv')
chr20 = read_tsv('chr20.tsv')
concordance = read_tsv('concordance.tsv')



# General format
format_supp_plot = function(p, title=NA, legend=F){
  p = p +
    scale_color_manual(name='', values=model_colors, labels=model_names) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
    p = p + ggtitle("")
  
  if(!legend){
    p = p + theme(legend.position = "None")
  }else{
    p = p + theme(legend.text=element_text(size=12))
  }
  
  return(p)
}

# Concordance plots format
format_concordance = function(p, same_scale){
  p =  p +
    scale_x_continuous(labels=percent_format(accuracy = 1))
  
  if(same_scale){
    p = p +
      scale_y_continuous(labels=percent_format(accuracy = 1),  limits = c(0.8, 1.0))
  } else {
    p = p +
      scale_y_continuous(labels=percent)
  }
  return(p)
}

# Get a row
get_plot_row = function(plots, xlabel, ylabel, labels, legend=F){
  plot = ggarrange(plotlist=plots, nrow=1, ncol=4, common.legend = legend,labels=labels)
  return( plot %>%  
            annotate_figure(
              left = text_grob(ylabel, rot = 90),
              bottom = text_grob(xlabel)
            )
  )
}

get_header_row = function(nrows_y_label){
  row = ggarrange(
    text_grob("Exomes SNVs", face='bold', just='center'),
    text_grob("Exomes Indels", face='bold', just='center'),
    text_grob("Genomes SNVs", face='bold', just='center'),
    text_grob("Genomes Indels", face='bold', just='center'),
    nrow = 1,
    ncol = 4
  )
  
  if(nrows_y_label>0){
    row = annotate_figure(row, left=strrep("\n",nrows_y_label-1))
  }
  return(row)
}

## PR Plot
pr_plot = function(same_scale = F){
  truth_samples = c('NA12878'='NA12878', 'syndip'='Synthetic diploid')
  rows = list(get_header_row(2))
  for(s in names(truth_samples)){
    row=list()
    for(x in c('exomes', 'genomes')){
      for(y in c('True',  'False')){
        plot_data = concordance %>%
          filter(data_type == x & 
                   snv == y & 
                   rank_name == "truth_sample_rank" & 
                   truth_sample ==  s)
        p = plot_data %>%
          ggplot(aes(recall, precision, col=model)) +
          geom_point() + 
          geom_vline(xintercept = plot_data %>% filter(model == 'RF' & bin==cutoffs[y]) %$% recall, linetype='dashed') +
          geom_hline(yintercept = plot_data %>% filter(model == 'RF' & bin==cutoffs[y]) %$% precision, linetype='dashed')
        p = format_supp_plot(p, legend = s==names(truth_samples[1]))
       row = c(row,list(format_concordance(p, same_scale)))
      }
    }
    rows = c(rows,
             list(get_plot_row(row,
                          paste(truth_samples[s],'recall'), 
                          paste(truth_samples[s],'\nprecision',sep=''),
                          labels=labels[(4*(length(rows)-1)+1):(4*(length(rows)-1)+4)], 
                          legend=s==names(truth_samples[1]))))
  }
  return(ggarrange(plotlist=rows, nrow=length(rows), ncol = 1, heights = c(0.1, 1.2, 1)))
}

# Rare variants metrics plots
rare_variants_metrics_plot = function(){

  #dnms
  dnms = autosomes %>%
    filter(rank_id == 'rank') %>%
    group_by(data_type, model, snv, bin) %>%
    summarise(n_de_novo=sum(n_de_novo)) %>%
    group_by(data_type, model, snv) %>%
    arrange(bin) %>%
    mutate(cum_dnm = cumsum(n_de_novo)/n_trios[data_type])
  
  dnm_plots=list()
  for(x in c('exomes', 'genomes')){
    for(y in c('True',  'False')){
      plot_data = dnms %>%
        filter(data_type == x & 
                 snv == y)
      p =  plot_data %>%
        ggplot(aes(bin, cum_dnm, col=model))  + 
        geom_point() + geom_vline(xintercept=cutoffs[y], linetype='dashed')
      dnm_plots = c(dnm_plots,
                    list(format_supp_plot(p, legend = T))
      )
    }
  }
  dnm_row = get_plot_row(dnm_plots, 
                         'Model percentile', 
                         bquote(atop('Number of' ~ italic('de novo') ~ 'calls', 'per child (cumulative)')), 
                         labels = labels[1:4],
                         legend = T)
  
  #trans singletons
  trans_singletons = chr20 %>%
    filter(rank_id == 'rank') %>%
    group_by(data_type, model, snv, bin) %>%
    summarise(n_trans_singletons=sum(n_trans_singletons)) %>%
    group_by(data_type, model, snv) %>%
    arrange(bin) %>%
    mutate(cum_trans_singletons = cumsum(n_trans_singletons))
  
  trans_singletons_plots=list()
  for(x in c('exomes', 'genomes')){
    for(y in c('True',  'False')){
      p = trans_singletons %>%
        filter(data_type == x & 
                 snv == y) %>%
        ggplot(aes(bin, cum_trans_singletons, col=model))  + 
        geom_point() + geom_vline(xintercept=cutoffs[y], linetype='dashed')
      trans_singletons_plots=  c(trans_singletons_plots,
                                 list(format_supp_plot(p))
      )
    }
  }
  trans_singletons_row = get_plot_row(trans_singletons_plots, 
                                      'Model percentile', 
                                      'Number of transmitted singletons\n(chromosome 20, cumulative)', 
                                      labels = labels[5:8],
                                      legend = F)
  
  #trans singletons
  validated_dnm = autosomes %>%
    filter(rank_id == 'rank') %>%
    group_by(data_type, model, snv, bin) %>%
    summarise(n_validated_de_novos=sum(n_validated_de_novos)) %>%
    group_by(data_type, model, snv) %>%
    arrange(bin) %>%
    mutate(cum_validated_de_novos = cumsum(n_validated_de_novos))
  
  validated_dnm_plots=list()
  for(x in c('exomes')){
    for(y in c('True',  'False')){
      p = validated_dnm %>%
        filter(data_type == x & 
                 snv == y) %>%
        ggplot(aes(bin, cum_validated_de_novos, col=model))  + 
        geom_point() + geom_vline(xintercept=cutoffs[y], linetype='dashed')
      validated_dnm_plots=  c(validated_dnm_plots,
                                 list(format_supp_plot(p))
      )
    }
  }
  validated_dnm_row = get_plot_row(validated_dnm_plots, 
                                      'Model percentile', 
                                      bquote(atop('Number of validated' ~ italic('de novo') ~ 'mutations','(cumulative)')), 
                                      labels = labels[9:10],
                                      legend = F)
  
  return(ggarrange(
    get_header_row(2),
    dnm_row, 
    trans_singletons_row, 
    validated_dnm_row, 
    nrow = 4, 
    ncol = 1,
    heights = c(0.1, 1.2, 1, 1)))
}
