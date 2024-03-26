library(dplyr)
library(openxlsx)
library(rstatix)
library("DescTools")

df = read.csv(
  "proteins-mics-genomes.csv",
  sep="\t", colClasses=c("genome_id"="character"), row.names="genome_id")

df[is.na(df)] <- 0

colnames(df)[which(names(df) == "X..measurement_value....ampicillin..")] <- "Amp"
colnames(df)[which(names(df) == "X..measurement_value....aztreonam..")] <- "Az"
colnames(df)[which(names(df) == "X..measurement_value....cefazolin..")] <- "CefZ"
colnames(df)[which(names(df) == "X..measurement_value....ceftazidime..")] <- "CefT"
colnames(df)[which(names(df) == "X..measurement_value....imipenem..")] <- "Ip"
colnames(df)[which(names(df) == "X..measurement_value....meropenem..")] <- "Mp"

colnames(df)[colnames(df) == "Est.PBP2"] <- "Est.PPBP"
colnames(df)[colnames(df) == "Gly"] <- "GlyII"

## MICs to categories
df$Amp <- factor(df$Amp,
                levels = c(0,2,4,8,16,32,64),
                labels = c("0","1","1","1","2","3","3"))
df$Az <- factor(df$Az,
                 levels = c(0,1,2,4,8,11,16,32,64),
                 labels = c("0","1","1","1","1","2","2",
                            "3","3"))
df$CefZ <- factor(df$CefZ,
                levels = c(0,1,2,4,8,16,32,64),
                labels = c("0","1","1","1","2","2",
                           "3","3"))
df$CefT <- factor(df$CefT,
                levels = c(0,0.5,1,2,4,6,8,16,32,64,128,256),
                labels = c("0","1","1","1","1","2","2","2",
                           "3","3","3","3"))
df$Ip <- factor(df$Ip,
                  levels = c(0,0.25,0.5,1,2,4,8,16,32,64,256),
                  labels = c("0","1","1","1","1","2","2","2",
                             "3","3","3"))
df$Mp <- factor(df$Mp,
                levels = c(0,0.015,0.03,0.06,0.12,0.125,0.25,0.5,
                           1,2,4,8,11,16,32,64,128),
                labels = c("0","1","1","1","1","1","1","1","1","1",
                           "2","2","2","2","3","3","3"))

columns_to_process <- c("A", "A2", "AHp", "AHL", "AKS", "B1.B2",
                        "B3", "BHp", "GlyII", "Human", "PqsE",
                        "VarG", "AmpH", "C",  "CHp",  "Est.PPBP",
                        "EstI",  "EstII", "PKS", "D", "PBP2")

prot_presen_ausenc <- df %>%
  mutate(across(all_of(columns_to_process), ~ ifelse(. > 0, 1, .)))

## For each antibiotic
ip <- subset(prot_presen_ausenc, Ip!=0)

ip <- ip %>% 
  dplyr::select(Specie, A, A2, B1.B2, B3, C, D, AHp, AHL, AKS, BHp, GlyII, Human, PqsE, VarG, AmpH, 
                CHp, Est.PPBP, EstI, EstII, PKS, PBP2, Ip)

result_df <- data.frame()

for (col_name in names(ip)) {
  freq_table <- table(ip[[col_name]])
  
  freq_df <- as.data.frame(freq_table)
  
  colnames(freq_df) <- c("Value", "Frequency")
  
  freq_df$Variable <- col_name
  
  freq_df$Percentage <- round(freq_df$Frequency / sum(freq_df$Frequency),3) * 100
  
  result_df <- rbind(result_df, freq_df)
}

## For each protein group
A <- table(ip$A, ip$Ip)
Af <- A[,2:4]

y <- fisher.test(Af)
par <- pairwise_fisher_test(Af)

r_2 <- Rev(Af[,c(1,2)])
r_3 <- Rev(Af[,c(1,3)])

rr_2 <- RelRisk(t(r_2), conf.level=0.95)
rr_3 <- RelRisk(t(r_3), conf.level=0.95)

pct_table <- epiDisplay::tabpct(ip$A, ip$Ip, graph = FALSE, decimal = 5)

names <- list('Count' = A, 'Row' = pct_table$table.row.percent, 
              'Col' = pct_table$table.column.percent, 'Fisher' = y$p.value,
              'Par' = par, 'RR_2' = rr_2, 'RR_3' = rr_3)
openxlsx::write.xlsx(names, file = 'a-ip.xlsx')
