library(ggrepel)
library(ggplot2)
library(glmnet)
require(doMC)
registerDoMC(cores=8)


file_loc = '/users/sweiss/downloads/'

#read data
ad_targets = read.csv(paste0(file_loc,'ad_targets.csv'))
ads = read.csv(paste0(file_loc,'ads.csv'))
targets = read.csv(paste0(file_loc,'targets.csv'))

head(ad_targets)
head(ads)
head(targets)

targets_1 = targets
targets_1[,'name'] = gsub("\\(", "", targets_1[,'name'])
targets_1[,'name'] = gsub("\\)", "", targets_1[,'name'])

targets_1[,'category'] = gsub("\\(", "", targets_1[,'category'])
targets_1[,'category'] = gsub("\\)", "", targets_1[,'category'])

targets_1[,'prefix'] = gsub("\\(", "", targets_1[,'prefix'])
targets_1[,'prefix'] = gsub("\\)", "", targets_1[,'prefix'])

targets_1[,'targets_2'] = mapply(function(x, y) gsub(paste0(x,':'), '', y),
  targets_1[,'prefix'],  targets_1[,'name'])

#merge targets with ads
ad_targets_merge = merge(ad_targets, targets_1, by.x = 'target_id', by.y = 'id', how = 'all.x')
ad_targets_merge[,'cat_targets'] = paste(ad_targets_merge[,'category'], ad_targets_merge[,'targets_2'], sep = '__')

agg_topics = aggregate(cat_targets ~ ad_id, FUN = paste, collapse = ";", data = ad_targets_merge)

#remove highlights
remove_highlights = c('placements__','language__','interests__','age__','location_living_in__','location__',
                      'people_who_match__','behaviors__','and_must_also_match__','Friends of people who are connected to',
                      'People who like ','and_must_also_match__',' ')

for(x in remove_highlights){agg_topics[,'cat_targets'] = gsub(x,'',agg_topics[,'cat_targets'])}

#find 500 most common targets
targets = unlist(strsplit(paste(agg_topics[,'cat_targets'], collapse = ';'), ';'))
top_x_targets = names(table(targets)[order(-table(targets))][1:501])
top_x_targets = top_x_targets[-which(top_x_targets %in% c('',' '))]

#create ad by target matrix 
targets_name_mat = data.frame(matrix(0,nrow(agg_topics), length(top_x_targets)))
colnames(targets_name_mat) = top_x_targets
for(x in top_x_targets){
  
  targets_name_mat[grep(x,agg_topics[,'cat_targets'], fixed = TRUE),x] = 1
}


targets_name_mat$id =agg_topics[,'ad_id']
merged_data = merge(targets_name_mat, ads, on.x = 'id', on.y = 'id')
merged_data = merged_data[order(merged_data[,'created']),]
merged_data = subset(merged_data, spend_currency == "RUB")
merged_data_x = merged_data[,2:501]
merged_data_y = merged_data[,502:511]

#cleanup dates and rubles
merged_data_y[,'new_spend_amount'] = merged_data_y[,'spend_amount']
merged_data_y[which(merged_data_y[,'spend_currency'] == 'USD'),'new_spend_amount'] = (subset(merged_data_y, spend_currency == 'USD')[,'spend_amount'])*60
merged_data_y[,'date'] = as.Date(as.character(substr(merged_data_y[,'created'], 1,10)), format = '%Y-%m-%d')
merged_data_y[,'date_ended'] = as.Date(as.character(substr(merged_data_y[,'ended'], 1,10)), format = '%Y-%m-%d')
merged_data_y[is.na(merged_data_y[,'date_ended']),'date_ended'] = merged_data_y[is.na(merged_data_y[,'date_ended']),'date']
merged_data_y[,'clicks_per_cost'] = log((merged_data_y['clicks']+1)/(1+merged_data_y[,'new_spend_amount']))
merged_data_y[,'impressions_per_cost'] = log((merged_data_y['impressions']+1)/(1+merged_data_y[,'new_spend_amount']))



get_regressions_over_time = function(x, y, start_date_x, window_x, window_y){
  end_date_x = start_date_x + window_x
  start_date_y = start_date_x + window_x + 1
  end_date_y = start_date_x + window_x + 1 + window_y

  to_keep_x = which(y[,'date']  > start_date_x & y[,'date_ended'] <= end_date_x)
  to_keep_y = which(y[,'date_ended']  > start_date_y & y[,'date'] <= end_date_y)

  print(paste('start_date_x',start_date_x))
  print(paste('end_date_x ',end_date_x))
  
  print(paste('start_date_y ',start_date_y))
  print(paste('end_date_y ',end_date_y))
  

  x_train = (x[to_keep_x, ])
  y_train = y[to_keep_x, ]
  
  diff_days = as.numeric(y_train[,'date_ended'] - y_train[,'date'])+1
  diff_days[is.na(diff_days)] = 1
  
  reg_costs = cv.glmnet(as.matrix(x_train), log(1+y_train[,'new_spend_amount']/diff_days), parallel = TRUE, alpha = 0)
  reg_clicks = cv.glmnet(as.matrix(x_train), log(1+y_train[,'clicks']/diff_days), parallel = TRUE, alpha = 0)
  reg_impressions = cv.glmnet(as.matrix(x_train), log(1+y_train[,'impressions']/diff_days), parallel = TRUE, alpha = 0)
  reg_clicks_cost = cv.glmnet(as.matrix(x_train), y_train[,'clicks_per_cost'], parallel = TRUE, alpha = 0)
  reg_impressions_cost = cv.glmnet(as.matrix(x_train), y_train[,'impressions_per_cost'], parallel = TRUE, alpha = 0)
  
  coefs_costs = as.matrix(coef(reg_costs, s = 'lambda.min'))
  coefs_clicks = as.matrix(coef(reg_clicks, s = 'lambda.min'))
  coefs_impressions = as.matrix(coef(reg_impressions, s = 'lambda.min'))
  coefs_clicks_cost = as.matrix(coef(reg_clicks_cost, s = 'lambda.min'))
  coefs_impressions_cost = as.matrix(coef(reg_impressions_cost, s = 'lambda.min'))
  
  counts_topics = c(0,colSums(x[to_keep_y, ]))
  counts_topics_current = c(0,colSums(x[to_keep_x, ]))
  all_coefs_counts = data.frame(coefs_costs = coefs_costs, coefs_clicks=coefs_clicks,
                                coefs_impressions=coefs_impressions, counts_topics,coefs_clicks_cost,
                                coefs_impressions_cost, start_date_x, 
                                counts_topics/sum(counts_topics), counts_topics/nrow(x_train),counts_topics_current)
  colnames(all_coefs_counts) = c('coefs_costs', 'coefs_clicks','coefs_impressions','counts_topics','coefs_clicks_cost',
                                 'coefs_impressions_cost','start_date_x', 'norm_counts_topics', 'prop_assigned','counts_topics_current')

  return(list(all_coefs_counts, nrow(x_train)))

}

library(GGally)
pairs_plot = ggpairs(log(1+merged_data_y[,c('clicks','impressions','spend_amount')]))+theme_minimal()+
  theme(text = element_text(size=15)) 
sort(table(merged_data_y[,'spend_amount']), decreasing = TRUE)[1:40]


most_common_topics = sort(colSums(merged_data_x[which(merged_data_y[,'date'] > as.Date("2016-01-01") & merged_data_y[,'date']< as.Date("2017-02-01")),which(colSums(merged_data_x)>10 )]), decreasing = TRUE)
most_common_topics = data.frame(names = names(most_common_topics), counts = most_common_topics)
most_common_topics[,'names']= factor(most_common_topics[,'names'], most_common_topics[order(-most_common_topics[,'counts']),'names'])
common_feature_names = ggplot(most_common_topics[1:15,], aes(x = names, y = counts))+geom_bar(stat = 'identity')+theme_minimal()+
  ggtitle('Counts of Top 15 Most Frequently Used Targets')

date_sequence = seq(as.Date("2016-01-01"), as.Date("2017-02-01"), "weeks")

reg_coefs = lapply( date_sequence,function(x) 
  get_regressions_over_time(merged_data_x[,which(colSums(merged_data_x)>10 )], merged_data_y, x,28,7))



reg_coefs_1 = lapply(reg_coefs, function(x) {temp = x[[1]]; temp$names = rownames(temp);return(temp)} )
reg_coefs_2 = do.call(rbind,reg_coefs_1)

together_data_1 =  subset(reg_coefs_2, names != '(Intercept)')

pairs_plot_coefs = ggpairs(together_data_1[,c('coefs_costs','coefs_clicks','coefs_impressions',
                                        'coefs_clicks_cost','coefs_impressions_cost','counts_topics')])+theme_minimal()+
  theme(text = element_text(size=15)) 

make_plots_over_time = function(dat){
  dat[,'names'] = rownames(dat)
  plt = ggplot(subset(dat, counts_topics>0 & coefs_costs!=0), aes(x = coefs_costs, y = coefs_clicks, size = counts_topics, label = names))+
    geom_point()+ geom_text_repel()+theme_minimal()+geom_abline(intercept = 0, slope = 1)+ggtitle(dat[1,'start_date_x'])+
    theme(text = element_text(size=15)) 
    
  return(plt)
}
example_11_4_plot = make_plots_over_time(reg_coefs_1[[45]])


reg = lm(scale(log(counts_topics+1)) ~ scale(coefs_costs)+scale(coefs_clicks)+scale(coefs_clicks_cost)+scale(log(counts_topics_current+1))+as.factor(start_date_x), data = together_data_1 )
summary(reg)

pander(reg)

library(randomForest)

random_sample = sample(nrow(subset(together_data_1, names != '(Intercept)')), 4000)

rf = randomForest(log(counts_topics+1) ~ coefs_costs + coefs_clicks+ counts_topics_current+start_date_x, data = subset(together_data_1, names != '(Intercept)')[random_sample,] )
rf_only_dates = randomForest(log(counts_topics+1) ~  start_date_x+counts_topics_current, data = subset(together_data_1, names != '(Intercept)')[random_sample,] )

library(pdp)
pd <- partial(rf, pred.var = c("coefs_costs", "coefs_clicks"), train = subset(together_data_1, names != '(Intercept)')[-random_sample,] )
# Default PDP

# Add contour lines and use a different color palette
rwb <- colorRampPalette(c("red", "white", "blue"))
pdp2 <- plotPartial(pd, contour = TRUE, col.regions = rwb)
plot(pdp2)


#most cost effective ads 
together_data_2 = aggregate(coefs_clicks_cost~names, data = together_data_1, FUN = mean)

together_data_2[,'names']= factor(together_data_2[,'names'], together_data_2[order(-together_data_2[,'coefs_clicks_cost']),'names'])
together_data_2 = together_data_2[order(together_data_2[,'coefs_clicks_cost'], decreasing = TRUE),]


coefs_clicks_plots = ggplot(together_data_2[1:10,], aes(x = names, y = coefs_clicks_cost))+geom_bar(stat = 'identity')+theme_minimal()+
    ggtitle('Most Efficient Targets')+xlab('Target Names')+ylab('Clicks / Costs Coefficients')  + theme(axis.text.x = element_text(angle = 90))

save(pdp2,common_feature_names,
     reg, example_11_4_plot,pairs_plot_coefs,pairs_plot, together_data_1,coefs_clicks_plots,
     file= '/Users/sweiss/src/IRA_Facebook_objectives/russian_objectives_plots.rdata')
