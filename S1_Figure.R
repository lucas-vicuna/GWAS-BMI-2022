library(data.table)
library(dplyr)
library(lubridate)
library(ggplot2)
library(nlme)
library(lme4)
library(lmerTest)

# setwd('C:/Users/ebarr/Downloads/')

# setwd('C:/Users/ebarr/Onedrive/Documentos/Genetica')
setwd('C:/Users/ebarr/Onedrive/Descargas')


file1 <- paste0("rfmix15.chr",22,".stats.phen.pel_ceu.txt")
df <- fread(file1, sep='\t', header=T, select = 1, col.names = 'CODE')
# df <- fread(file1, sep='\t', header=T)


df1 <- fread('BaseMarzo2020.csv') 
# df1 <- readr::read_csv('BaseMarzo2020.csv') 

colnames(df1) <- c("CODE","date_birth","Sex","vis_date","Weight","Height","BMI","BAZ")

df1 <- df1 %>% inner_join(df)

df2 <- df1 %>% 
  mutate(vis_date = mdy(vis_date),
         date_birth = mdy(date_birth),
         Age = as.numeric(difftime(vis_date, date_birth, units = c("days"))/365)
  ) %>% 
  select(CODE, vis_date, Sex, Weight, Height, BMI, Age) %>% 
  filter(!is.na(BMI)) %>% 
  group_by(CODE) %>% 
  arrange(vis_date) %>% 
  mutate(time = row_number()) %>% 
  ungroup()

df46 <- 
  bind_rows(
    # df2 %>% transmute(CODE, BMI, Age, Sex = 'all', tipo = 'observed'),
    # df2 %>% transmute(CODE, BMI, Age, Sex, tipo = 'observed'),
    df2 %>% transmute(CODE, BMI = BMI_pred8, Age, Sex = 'all', tipo = 'fitted'),
    df2 %>% transmute(CODE, BMI = BMI_pred8, Age, Sex, tipo = 'fitted') ) %>% 
  filter(between(Age,2,10)) %>% 
  filter(CODE %in% (M8 %>% filter(AR>1) %>% .$CODE)) %>% 
  mutate(sex2 = case_when(
    Sex=='Masculino' ~ 'Male',
    Sex=='Femenino' ~ 'Female',
    TRUE ~ 'all'))

df46$inter2 <- interaction(df46$sex2, df46$tipo)

p146 <- 
df46 %>% 
  mutate(inter=inter2) %>%  #= factor(inter2, levels = c('all.observed','all.fitted','Female.observed','Female.fitted','Male.observed','Male.fitted'))) %>% 
  ggplot() +
  aes(x = Age, y = BMI, color = inter) +
    # scale_linetype_manual('', values=c('solid','11','solid','11','solid','11'), labels = c('All observed','All fitted','Female observed','Female fitted','Male observed','Male fitted')) +
  geom_segment(aes(x = 2, y = 16.46, xend = 4.55, yend = 16.46), linetype='dashed',color="lightgray", size=0.5)+
  geom_segment(aes(x = 4.55, y = 15.3, xend = 4.55, yend = 16.46), linetype='dashed',color="lightgray", size=0.5)+
  geom_segment(aes(x = 2, y = 16.60, xend = 4.35, yend = 16.60), linetype='dashed',color="lightgray", size=0.5)+
  geom_segment(aes(x = 4.35, y = 15.3, xend = 4.35, yend = 16.60), linetype='dashed',color="lightgray", size=0.5)+
  geom_smooth(se=F, formula = y ~ poly(x,3), method = 'lm',size=0.5) +
  scale_color_manual('', values=c('#E7B800','#00AFBB','#FC4E07'), labels = c('All','Girls','Boys')) +
  labs(y = expression(BMI~(kg/m^2)), x = 'Age (years)') +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.key.size = unit(2, 'mm'),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 20)) +
  theme(plot.margin = unit(c(1,1,2,1), "lines")) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim=c(15.5,19.5)) +
  annotation_custom(texty,xmin=1.5,xmax=1.5,ymin=16.5,ymax=16.5) +
  annotation_custom(textx,xmin=4.6,xmax=4.6,ymin=15.15,ymax=15.15)

gt <- ggplotGrob(p146)
gt$layout$clip[gt$layout$name=="panel"] <- "off"

pdf('3curvas.pdf', width=7, height=7*611/731)
grid.draw(gt)
dev.off() 

p146 %>% class
ggsave('Curve by Sex - observed vs fitted values.png', scale=1,width=15)