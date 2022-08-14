**Identifying heterogeneous subgroups of systemic connective tissue
diseases by applying a joint dimension reduction and clustering approach
to immunomarkers : a retrospective study** <br/>- cluster analysis
================
**Chia-Wei Chang<sup>a,#</sup>, Hsin-Yao Wang<sup>b,#</sup>, Wei-Lin
Lo<sup>c</sup>, Wei-Ting Lin<sup>b</sup>, Jia-Ruei Yu<sup>b</sup>, Yi-Ju
Tseng<sup>a,d,\*</sup>** <br/> <sup>a</sup> Department of Computer
Science, National Yang Ming Chiao Tung University, Hsinchu, Taiwan <br/>
<sup>b</sup> Department of Laxboratory Medicine, Chang Gung Memorial
Hospital at Linkou, Taoyuan City, Taiwan <br/> <sup>c</sup> Department
of Rheumatology, Chang Gung Memorial Hospital at Keelung, Keelung City,
Taiwan <br/> <sup>d</sup> Computational Health Informatics Program,
Boston Children’s Hospital, Boston, MA, USA <br/> <sup>\#</sup> Chang
and Wang contribute equally to this work <br/> <sup>\*</sup>
Corresponding Author <br/>

------------------------------------------------------------------------

# Set up for environment

## Load packages from CRAN

``` r
.cran_pkgs <- c("data.table",
                "dplyr",
                "magrittr",
                "stringr",
                "lubridate",
                "scales",
                "tidyr",
                
                "purrr",
                "furrr",
                "devtools",
                "knitr",
                
                "modeest",
                "vcd",
                "psych",
                "DescTools",
                "clustrd",
                "yardstick",

                "ggplot2",
                "ggsci",
                "ggrepel",
                "ggforce",
                "ggpubr",
                "cowplot",
                
                "tableone",
                "flextable"
                
                )

if (any(!.cran_pkgs %in% installed.packages())) {
  install.packages(.cran_pkgs[!.cran_pkgs %in% installed.packages()],
                   lib = Sys.getenv("R_LIBS_USER"),
                   dependencies = TRUE)
}

sapply(.cran_pkgs,
       function(x) suppressPackageStartupMessages(require(x,character.only = TRUE)))
```

    ## data.table      dplyr   magrittr    stringr  lubridate     scales      tidyr 
    ##       TRUE       TRUE       TRUE       TRUE       TRUE       TRUE       TRUE 
    ##      purrr      furrr   devtools      knitr    modeest        vcd      psych 
    ##       TRUE       TRUE       TRUE       TRUE       TRUE       TRUE       TRUE 
    ##  DescTools    clustrd  yardstick    ggplot2      ggsci    ggrepel    ggforce 
    ##       TRUE       TRUE       TRUE       TRUE       TRUE       TRUE       TRUE 
    ##     ggpubr    cowplot   tableone  flextable 
    ##       TRUE       TRUE       TRUE       TRUE

## Load external functions

``` r
source("./External_R_Functions/ggRadar2.R")
source("./External_R_Functions/cramer'V matrix.R")
```

<br/>

# Import datasets

## Analysis-ready dataset[^1]

From
[Data_Wrangling.md](https://github.com/DHLab-TSENG/Heterogeneity-in-SCTD-paper/blob/main/Data_Wrangling.md)

``` r
monoCTD_dataset <- 
  readRDS("./Dataset/A5_monoCTD_ready_analysis_dataset.rds") %>%
  # Subset subjects first diagnosed between 2001 and 2016
  .[data.table::between(year(First_diagnosis_date),lower = 2001,upper = 2016),]
```

<br/>

## Query table for immunomarkers

``` r
simul_classify_item <- 
  list(c("72-134","M25-238"),                  # Rheumatoid factor (RF)
       c("72-245","72A245","72B245","72C245"), # Antinuclear antibody (ANA)
       c("72-247","M25-150"))                  # anti-Double strand DNA antibody (anti-dsDNA)

removed_items <- c("ESR","CRP","CRP-Emr")      # Remove CRP, CRP-Emr, and ESR items

exam_table <- 
  fread("./Dataset/exam_table.csv") %>% 
  .[,ITEM_LABSH1IT := paste(ITEM,LABSH1IT,sep = "_")] %>% 
  .[,.SD,.SDcols = c("LABIT","LABSH1IT","ITEM","ITEM_LABSH1IT","LABNMABV")] %>% 
  # Remove redundant items in simultaneously classify items 
  .[!LABSH1IT %in% c(map(simul_classify_item, ~ .x[-1],) %>% flatten_chr(.)),] %>% 
  # Remove CRP, CRP-Emr, and ESR items
  .[str_detect(ITEM,paste(removed_items,collapse = "|"),negate = TRUE),] %>% 
  .[,unique(.SD)]
```

<br/>

------------------------------------------------------------------------

# Data preparation for feature selection

## Transform datset into wide format by immunomarker

``` r
cluster_analysis_data <- 
  copy(monoCTD_dataset) %>% 
   .[,':=' (ITEM     = str_split(ITEM_LABSH1IT,"_",simplify = TRUE)[,1],
            LABSH1IT = str_split(ITEM_LABSH1IT,"_",simplify = TRUE)[,2])] %>% 
   .[str_detect(ITEM,paste(removed_items,collapse = "|"),negate = TRUE),] %>% 
   .[,Group := as.character(Group)] %>%
   .[,.SD,.SDcols = -c("First_diagnosis_date","LABSH1IT","ITEM_LABSH1IT","SCDATE")] %>% 
   .[,AGE := mlv(AGE,method = "mfv"),by = .(IDCODE)] %>% 
   dcast(.,... ~ ITEM,value.var = "Imp_Status",fun.aggregate = toString,fill = "-") %>% 
  ## Substitute complete with abbreviation
   setnames(.,exam_table$ITEM,exam_table$LABNMABV) %>%
  ## Preprocess for MCA-Kmeans
   .[,c(exam_table$LABNMABV) := map(.SD, ~ str_replace_all(.x,c("Normal" = "N_","Abnormal" = "AbN_"))),
     .SDcols = exam_table$LABNMABV] %>%
   .[,c(exam_table$LABNMABV) := map2(.SD,exam_table$LABNMABV, ~ paste0(.x,.y)),
     .SDcols = exam_table$LABNMABV] %>%
   .[,c(exam_table$LABNMABV) := map(.SD, ~ factor(.x)),
     .SDcols = exam_table$LABNMABV]
```

<br/>

## Subset for selected SCTDs

``` r
selected_disease_groups <- 
  c("Sjogren's syndrome","Rheumatoid arthritis","Systemic lupus erythematosus")

cluster_analysis_data <- 
  cluster_analysis_data[Group %in% selected_disease_groups,]
```

<br/>

## Remove variables with single unique value (no information or all Normal)

``` r
redundant_vars <- 
  copy(cluster_analysis_data) %>% 
  .[,map(.SD, ~ uniqueN(.x)),.SDcols = exam_table$LABNMABV] %>% 
  { names(.)[str_which(.,"1")] }
```

<br/>

# Feature selection

## Kaiser–Meyer–Olkin (KMO) test for the measure of adequacy (MSA)

``` r
MSA_cutoff  <- 0.6
MSAi_cutoff <- 0.5
redundant_iter_vars <- redundant_vars

repeat {
  
  # Construct the correlation coefficient matrix (remove redundant variables in advance)
  cramerV_iter_dataset <- 
    if (length(redundant_iter_vars)>0) { 
      copy(cluster_analysis_data) %>% 
        .[,.SD,.SDcols = -(redundant_iter_vars)] %>% 
        cramerV_matrix(dataset = .,
                       vars = str_subset(exam_table$LABNMABV,
                                         paste(redundant_iter_vars,collapse = "|"),
                                         negate = TRUE)) %>% 
        setNames(.,
                 str_subset(exam_table$LABNMABV,
                            paste(redundant_iter_vars,collapse = "|"),
                            negate = TRUE))
      } else { 
        
        cluster_analysis_data[,.SD] %>% 
        cramerV_matrix(dataset = .,
                       vars = exam_table$LABNMABV) %>% 
        setNames(.,exam_table$LABNMABV)
        
      }
  
  # Prepare the correlation coefficient matrix for the KMO test
  mod_cramerV_iter_dataset <- 
    copy(cramerV_iter_dataset) %>% 
    as.data.table(.,keep.rownames = "ITEM") %>% 
    melt.data.table(.,id.vars = "ITEM",variable.name = "ITEM_2") %>% 
    .[,ITEM := factor(ITEM,levels = levels(ITEM_2))] %>% 
    dcast.data.table(.,... ~ ITEM_2,value.var = "value") %>% 
    as.matrix(.,rownames = "ITEM")
  
  # KMO test 
  KMO_iter_results  <- KMO(mod_cramerV_iter_dataset)
  
  # iterate till overall MSA is over the specified cutoff value
  if (KMO_iter_results$MSA>=MSA_cutoff) break
  
  # Find the variable with the smallest MSAi value to remove
  min_MSAi_var <- KMO_iter_results$MSAi[KMO_iter_results$MSAi <= MSAi_cutoff] %>% 
                  which.min(.) %>% 
                  names(.)
  
  redundant_iter_vars <- append(redundant_iter_vars,min_MSAi_var)

}

print(KMO_iter_results)
```

    ## Kaiser-Meyer-Olkin factor adequacy
    ## Call: KMO(r = mod_cramerV_iter_dataset)
    ## Overall MSA =  0.63
    ## MSA for each item = 
    ##      DC IgG          RF   FLC Kappa  FLC Lambda          C3          C4 
    ##        0.75        0.76        0.50        0.50        0.76        0.77 
    ##      CRYOID   CRYOFIBRI         ANA  Anti-dsDNA   Anti-ICSA   Anti-BMZA 
    ##        0.55        0.51        0.78        0.80        0.71        0.66 
    ##    Anti-TPO   Anti-THYG AQP4-autoAb      GAD-Ab   Anti-TSHR Anti-RNP-Sm 
    ##        0.50        0.50        0.69        0.43        0.44        0.69 
    ##    Anti-SSA    Anti-SSB     Anti-Jo    Anti-Scl   Anti-CENP   Anti-RibP 
    ##        0.52        0.52        0.48        0.67        0.52        0.60 
    ##        ACAG        ACAM      B2GP1G 
    ##        0.75        0.52        0.56

<br/>

## Check whether the correlation matrix is an identity matrix by Bartlett’s test of Sphericity

``` r
Bart_results <- 
  cortest.bartlett(mod_cramerV_iter_dataset,
                   # n = sqrt(nrow(cluster_analysis_data)))
                   n = nrow(cluster_analysis_data)/100)

print(Bart_results$p.value)
```

    ## [1] 0.02690535

<br/>

# Data preparation for cluster analysis

## Specify variables to be removed

``` r
removed_variables <- 
  # remove variables with single unique value (all results are Normal)
  copy(cluster_analysis_data) %>% 
  .[,map(.SD, ~ n_distinct(.x)),.SDcols = exam_table$LABNMABV] %>% 
  { names(.)[str_which(.,"1")] } %>% 
  
  # remove variables not passing the KMO test for sampling of adequacy
  c(.,setdiff(exam_table$LABNMABV,names(KMO_iter_results$MSAi))) %>% 
  unique(.) %>%
  # modify variable names that have been changed due to previous steps
  str_replace_all(.,"\\.","\\-") %>% 
  str_replace_all(.,c("DC-IgG" = "DC IgG","FLC-Kappa" = "FLC Kappa","CH50-Assay" = "CH50 Assay"))
```

<br/>

## Remove the specified variables from the data set for cluster analysis

``` r
cluster_analysis_ready_data <-  
  copy(cluster_analysis_data) %>% 
  .[,.SD,.SDcols = -c(removed_variables)] %>% 
  .[,.SD,.SDcols = -c("IDCODE","Group","AGE","SEX")]
```

<br/>

# Cluster analysis

## MCA k-means

``` r
plan(multisession, workers = 5)

alphak_vector <- 0.5
seed_vector   <- 365 # for MAC OS: 365; on Linux OS, use 20211111 instead.
result_list_name <- paste("seed =",seed_vector)

monoCTD_MCAkmeans_tune <-
  future_map(seed_vector,
             ~ tuneclus(data = cluster_analysis_ready_data,
                        nclusrange = 3:6,
                        ndimrange = 2:4,
                        criterion = "asw",
                        method = "MCAK",
                        alphak = alphak_vector,
                        seed = .x)
             ) %>% 
  setNames(.,
           result_list_name)

plan(sequential)
```

<br/>

------------------------------------------------------------------------

# Data profile for demographic and immunomarker variables

## Prepare clustering result datasets to extract data for the profile

``` r
PC_data <-
  # Retrieve the best result object
  pluck(monoCTD_MCAkmeans_tune,result_list_name[1]) %>% 
  pluck(.,"clusobjbest") %>% 
  # Retrieve PC coordinate and cluster data
  { cbind(pluck(.,"obscoord"),
          pluck(.,"cluster") ) } %>% 
  as.data.table(.) %>% 
  # Set names
  { setnames(.,c(paste("PC",seq(1,ncol(.)-1),sep = " "),"Cluster")) } %>% 
  # Add Group column
  cbind(.,
        cluster_analysis_data) %>% 
  # Factorise Cluster column
  .[,':=' (Cluster = factor(Cluster,
                            levels = seq(min(Cluster),max(Cluster),by = 1)),
           Group   = factor(Group,
                            levels = selected_disease_groups))]


attr_data <- 
  # Retrieve best result object
  pluck(monoCTD_MCAkmeans_tune,result_list_name[1]) %>% 
  pluck(.,"clusobjbest") %>% 
  # Retrieve attribute data
  pluck(.,"attcoord") %>% 
  as.data.table(.,keep.rownames = TRUE) %>% 
  { setnames(.,c("Attr",paste("PC",seq(1,ncol(.)-1),sep = " "))) } %>% 
  # Modify values in the attribute column and create the status column for plotting purposes
  .[,Attr_mdf := str_sub(Attr,
                         str_locate(Attr,"\\.AbN|\\.N")[,"start"]) %>% 
                 str_replace(.,"\\.","")] %>% 
  .[,Status   := str_split(Attr_mdf,"_",simplify = TRUE)[,1] %>% 
                 factor(.,levels = c("N","AbN"))]
```

<br/>

## Data preprocess

``` r
auto_antibody_string <- unique(sapply(str_split(attr_data$Attr_mdf,"_"),`[`,2))

table_one_by_diseases <- 
  # data preps
  copy(PC_data) %>% 
  { .[,.SD,.SDcols = -c(str_subset(names(.),"^PC|^Cluster"))] } %>% 
  .[,SEX := factor(SEX,levels = c("M","F"),labels = c("Male","Female"))] %>% 
  .[,c(auto_antibody_string) := map(.SD, ~ sapply(str_split(.x,"_"),`[`,1)),.SDcols = auto_antibody_string] %>% 
  .[,c(auto_antibody_string) := map(.SD, ~ factor(.x,levels = c("N","AbN"))),.SDcols = auto_antibody_string] %>% 
  # create the table one
  CreateTableOne(data = .,
                 strata = "Group",
                 vars = c("AGE","SEX",auto_antibody_string),
                 factorVars = c("SEX",auto_antibody_string),
                 test = TRUE,
                 testExact = TRUE,
                 smd = FALSE,
                 addOverall = FALSE)

# p values calculation
p_value_table <- 
  map(list(table_one_by_diseases$ContTable,
           table_one_by_diseases$CatTable) %>% 
      compact(.),
      ~ attr(.x,"pValues") %>% 
        as.data.table(.,keep.rownames = TRUE) %>% 
        .[,-3] %>% 
        setnames(.,c("Variable","p value")) %>% 
        .[,"p value" := ifelse(`p value`<0.001,"<0.001",sprintf("%.3f",`p value`))]) %>% 
    rbindlist(.)

# statistics calculations for continuous variables by SCTD group
table_one_continuous_var <- 
  # extract desired statistics and merge them into one data.table by SCTD group
  map(table_one_by_diseases$ContTable,
      ~ .x[,c("n","mean","sd")] %>% 
        as.data.table(x = .,keep.rownames = "Statistics") %>% 
        setnames(.,".","value")) %>% 
    rbindlist(.,use.names = TRUE,idcol = "Group") %>% 
  # retrieve data and modify values for an easy-to-read purpose
  dcast.data.table(., ... ~ Statistics,value.var = "value") %>% 
  .[,':=' (n    = comma(n),
           mean = sprintf("%.1f",mean),
            sd   = sprintf("%.1f",sd))] %>% 
  .[,"mean (sd)" := paste0(sprintf("%3s ",mean),"(",sprintf("%4s",sd),")")] %>% 
  # transform dataset for the tabulation purpose
   melt.data.table(.,id.vars = "Group",variable.name = "Statistics",value.name = "value") %>% 
   dcast.data.table(., ... ~ Group,value.var = "value") %>% 
  .[Statistics %in% "mean (sd)",] %>% 
  .[,':=' (Ori_Variable = "AGE",
           Variable     = "Age")] %>% 
  merge(.,
        p_value_table,
        by.x = "Ori_Variable",
        by.y = "Variable")
  

# statistics calculations for categorical variables by SCTD group
table_one_categorical_var <- 
  # extract desired statistics and merge them into one data.table by SCTD group
  map(table_one_by_diseases$CatTable,
      ~ map(.x, ~ .x[,c("level","freq","percent")]) %>% 
        rbindlist(.,idcol = "Ori_Variable")
      ) %>% 
  rbindlist(.,idcol = "Group") %>% 
  # retrieve data and modify values for an easy-to-read purpose
  .[!level %in% "N",] %>% 
  .[,':=' (freq    = as.integer(freq) %>% comma(.,accuracy = 1),
           percent = sprintf("%.2f",percent) %>% str_replace_all(.,"0.00","0"))] %>% 
  .[,"n (%)" := sprintf("%3s (%s)",freq,percent)] %>% 
  .[Ori_Variable=="SEX",Variable := level] %>% 
  .[is.na(Variable), Variable := Ori_Variable] %>% 
  # transform dataset for the tabulation purpose
  .[,.SD,.SDcols = -c("level")] %>% 
  dcast.data.table(.,Variable + Ori_Variable ~ Group, drop = TRUE,value.var = "n (%)") %>% 
  .[,"Statistics" := "n (%)"] %>% 
  merge(.,
        p_value_table,
        by.x = "Ori_Variable",
        by.y = "Variable")
```

<br/>

## Generate a table for the data profile[^2]

``` r
# Create customised column title
table_lable <- 
  map(table_one_by_diseases$ContTable,
      ~ .x[,1] %>% as.data.table) %>% 
  rbindlist(.,idcol = "Group") %>% 
  setnames(.,c("Group","n")) %>% 
  .[,table_lable := paste0(Group,"\n","(n = ",comma(n),")") %>% 
                    str_replace_all(.,"Sjogren's syndrome","Sjogren's\nsyndrome") %>% 
                    str_replace_all(.,"Rheumatoid arthritis","Rheumatoid\narthritis") %>% 
                    str_replace_all(.,"Systemic lupus erythematosus","Systemic lupus\nerythematosus")] %>% 
  .[,table_lable]


# For substituting immunomarker abbreviations with complete names
var_short_names <- str_subset(table_one_categorical_var$Variable,paste(exam_table$LABNMABV,collapse = "|"))
var_long_names  <- exam_table[order(LABNMABV),][LABNMABV %in% var_short_names,ITEM]

# Table one template
table_one_template <- 
  map(list(table_one_continuous_var,
           table_one_categorical_var),
      ~ .x[is.na(Variable),Variable := Ori_Variable] %>% 
        .[,.SD,.SDcols = -c("Ori_Variable")]) %>% 
  rbindlist(.,use.names = TRUE) %>% 
  # Variable group assigning
    .[Variable %in% c("Age"),':=' (Variable = paste("Age",Statistics,sep = ", "),
                                   Group    = NA_character_)] %>% 
    .[Variable %in% c("Male","Female"),Group := paste("Sex",Statistics,sep = ", ")] %>% 
    .[Variable %in% c("C3","C4"), Group := paste("Complement(= abnormal)",Statistics,sep = ", ")] %>% 
    .[str_detect(Variable,"Age|Male|Female|C3|C4",negate = TRUE), Group := paste("Auto-antibodies(= positive)",Statistics,sep = ", ")] %>% 
    .[,Group := factor(Group,levels = c(NA_character_,
                                        "Sex, n (%)",
                                        "Auto-antibodies(= positive), n (%)",
                                        "Complement(= abnormal), n (%)"),
                       exclude = NULL)] %>% 
  # Substitute immunomarker abbreviations with complete descriptions
    setkey(Variable) %>% 
    .[var_short_names,Variable := var_long_names] %>% 
  # For an easy-to-read purpose
    .[str_detect(Variable,"Age",negate = TRUE),Variable := sprintf("\t%s",Variable)] %>% 
    .[,.SD,.SDcols = -c("Statistics")] %>% 
    setcolorder(.,c("Group","Variable",names(table_one_by_diseases$CatTable))) %>% 
    .[order(Group,na.last = FALSE),] %>% 
  as_grouped_data(., groups = "Group",columns = NULL) %>% 
  setnames(.,c(names(table_one_by_diseases$CatTable)),table_lable)


# Set output formats
set_flextable_defaults(font.family = "Calibri",
                       font.size = 10)

# Table one generation
table_one_flextable <- 
  flextable(table_one_template,
            col_keys = c("Variable",table_lable,"p value"),
            cwidth =c(1.45,1,1,1),
            theme_fun = 
            ) %>% 
    compose(.,
            i = ~!is.na(Group),
            j = "Variable",
            value = as_paragraph(as_chunk(Group))) %>% 
    align(.,j = c(table_lable,"p value"),align = "right",part = "all") %>% 
  autofit(.)
```

<br/>

------------------------------------------------------------------------

# Presentation for clustering results

## Biplot of the first and second principal components

``` r
variable_list <- 
  map(c(sum(str_detect(names(PC_data),"^PC")):2),
      ~ paste0("`","PC ",c(.x-1,.x),"`") %>% 
        c(.,"Status","Attr_mdf","Cluster","Group")) %>% 
  rev(.)

cluster_color_set <- 
  c(pal_nejm(alpha = 0.8)(n_distinct(PC_data$Cluster))[-6],"grey10")

cluster_biplots <- 
  map(variable_list,
      ~ ggplot() +
          geom_point(data = PC_data[,unique(.SD),.SDcols = -c("IDCODE","SEX","AGE")],
                     aes_string(x = .x[1],y = .x[2],color = .x[5],shape = .x[6]),
                     position = position_jitter(height = 0.001,width = 0.001),
                     size = 1.5,
                     alpha = 0.5) +
          geom_point(data = attr_data[Status %in% c("AbN"),],
                     aes_string(x = .x[1],y = .x[2]),
                     size = 2,
                     shape = "cross",
                     color = "black",
                     alpha = 0.3) +
          geom_text_repel(data = attr_data[Status %in% c("AbN"),] %>% 
                                 .[,Attr_mdf := str_replace_all(Attr_mdf,"AbN_","")],
                          aes_string(x = .x[1],y = .x[2],label = .x[4]),
                          size = rel(3),
                          alpha = 0.2,
                          min.segment.length = 5,
                          box.padding = .3,
                          max.overlaps = 10) +
          geom_vline(xintercept = 0,color = "black",alpha = 0.4) +
          geom_hline(yintercept = 0,color = "black",alpha = 0.4) +
          geom_mark_ellipse(data = PC_data[,unique(.SD),.SDcols = -c("IDCODE","SEX","AGE")],
                            aes_string(x = .x[1],y = .x[2],fill = .x[5]),color = NA,
                            alpha = 0.2) +
          scale_linetype_discrete(name = "Immunomaker Test Results") +
          scale_shape_discrete(name = "SCTD",labels = map_chr(selected_disease_groups,
                                                              ~ str_split(.x,"\\ ") %>%
                                                                map_chr(.,~ str_sub(.x,1,1) %>% toupper(.) %>% str_flatten(.)))) +
          scale_colour_manual(name = "Cluster",values = cluster_color_set) +
          scale_fill_manual(name = "Cluster",values = cluster_color_set) +
          scale_x_continuous(name = str_replace_all(.x[1],"\\`","")) +
          scale_y_continuous(name = str_replace_all(.x[2],"\\`","")) +
          guides(shape = guide_legend(override.aes = list(size = 2.5)),
                 color = guide_legend(override.aes = list(size = 2.5),nrow = 1)) +
          theme_bw() +
          theme(axis.title.x    = element_text(face = "bold",size = rel(1),margin = margin(.5,0,0,0,"cm")),
                axis.title.y    = element_text(face = "bold",size = rel(1),margin = margin(0,.5,0,0,"cm")),
                axis.text       = element_text(face = "bold",size = rel(.9)),
                legend.title    = element_text(face = "bold",size = rel(1)),
                legend.text     = element_text(size = rel(.9)),
                legend.position = "bottom",
                legend.box.margin = margin(0,2,0,0,"cm"))
      )
```

<br/>

<br/>

## Profile the number and proportion of cases in each cluster

``` r
cluster_size_string <- 
  copy(PC_data) %>% 
  .[,.N,by = "Cluster"] %>% 
  .[,strings := label_comma()(N)] %>% 
  .[order(Cluster),]

profile_by_cluster_data <- 
  copy(PC_data) %>% 
    .[,.N,by = c("Cluster","Group")] %>% 
    .[order(Cluster,Group),] %>%
    .[,Cluster_label := factor(Cluster,
                               levels = cluster_size_string$Cluster,
                               labels = cluster_size_string$strings)] %>% 
    .[,Cluster_size       := sum(N),by = c("Cluster")] %>% 
    .[,Cluster_proportion := sqrt(round(Cluster_size/sum(N),3))*0.95] %>% 
    .[,Group_proportion   := round(N/Cluster_size,3)] %>% 
    .[,Group_position   := cumsum(Group_proportion)-0.5*Group_proportion,by = c("Cluster")] %>% 
    .[,Group := factor(Group,
                       levels = selected_disease_groups,
                       labels = map_chr(selected_disease_groups,~ str_split(.x,"\\s") %>% 
                                                                  map_chr(., ~ str_sub(.x,1,1) %>% 
                                                                  str_to_upper(.) %>% 
                                                                  str_c(.,collapse = ""))))]

profile_by_cluster_holistic_data <- 
  merge(profile_by_cluster_data,
        { copy(profile_by_cluster_data) %>% 
          .[,unique(.SD),.SDcols = c("Cluster","Cluster_proportion")] %>% 
          .[,x_position := round(cumsum(Cluster_proportion/sum(Cluster_proportion))-2.2*Cluster_proportion/sum(Cluster_proportion),3)] %>%
          .[,.SD,.SDcols = -c("Cluster_proportion")] },
        by = "Cluster",
        all.x = TRUE
  ) %>% 
  .[,Cluster_label := factor(Cluster,
                             levels = sort(unique(cluster_size_string$Cluster)),
                             labels = cluster_size_string$strings)]


profile_by_cluster_integrated <- 
  ggplot(profile_by_cluster_holistic_data,aes(width = Cluster_proportion)) +
      geom_bar(aes(x = x_position,y = Group_proportion,fill = Group),
               stat = "identity",position = "stack",alpha = 0.8) +
      scale_x_continuous(name = "Cluster\n(n)",
                         breaks = unique(profile_by_cluster_holistic_data$x_position),
                         labels = unique(profile_by_cluster_holistic_data$Cluster),
                         expand = c(0.001,0.001)) +
      scale_y_continuous(name = "Percentage",
                         breaks = seq(0,1,by = 0.2),
                         labels = seq(0,100,by = 20),
                         expand = c(0,0),) +
      scale_fill_manual(name = "CTDs",values = pal_jama()(3)) +
      facet_grid(~ Cluster_label, space = "free_x", scales = "free_x",switch = "x") +
      theme_bw() +
      theme(legend.title = element_text(face = "bold",size = rel(1)),
            legend.text = element_text(face = "bold",size = rel(.9)),
            legend.position = "bottom",
            legend.key.size = unit(.8,"line"),
            axis.title.x = element_text(face = "bold",size = rel(1),margin = margin(.3,0,0,0,unit = "cm")),
            axis.title.y = element_text(face = "bold",size = rel(1),margin = margin(0,.3,0,0,unit = "cm")),
            axis.text = element_text(face = "bold",size = rel(.9)),
            strip.placement = "outside",
            strip.background = element_blank(),
            strip.text = element_text(face = "bold",size = rel(.6),vjust = 3.5),
            panel.border = element_blank()
            )
```

<br/>

<br/>

## Profile the number and proportion of case in each SCTD group

``` r
cluster_size_string <- 
  copy(PC_data) %>% 
  .[,.N,by = "Group"] %>% 
  .[order(-N),] %>% 
  .[,Group_abr := factor(Group,
                         levels = unique(Group),
                         labels = map_chr(Group,~ str_split(.x,"\\s") %>% 
                                                  map_chr(., ~ str_sub(.x,1,1) %>% 
                                                               str_to_upper(.) %>% 
                                                               str_c(.,collapse = ""))))] %>% 
  .[,strings := label_comma()(N)]

profile_by_group_data <- 
  copy(PC_data) %>% 
    .[,.N,by = c("Cluster","Group")] %>% 
    .[order(Group,-Cluster),] %>% 
    .[,Group_label := factor(Group,
                             levels = cluster_size_string$Group,
                             labels = cluster_size_string$strings)] %>%
    .[,Group_size := sum(N),by = "Group"] %>% 
    .[,Cluster_proportion_within_group := round(N/Group_size,4)] %>% 
    .[,Cluster_proportion_within_group_position := cumsum(Cluster_proportion_within_group)-0.5*Cluster_proportion_within_group,by = c("Group")] %>% 
    .[,Group := factor(Group,
                       levels = cluster_size_string$Group,
                       labels = cluster_size_string$Group_abr)]

profile_by_group_holistic_data <- 
  merge(profile_by_group_data,
        { copy(profile_by_group_data) %>% 
          .[,unique(.SD),.SDcols = c("Group","Group_size")] %>% 
          .[,Group_proportion := round(Group_size/sum(Group_size),3),] %>% 
          .[,x_position := round(cumsum(Group_proportion/sum(Group_proportion))-0.5*Group_proportion/sum(Group_proportion),3)] %>% 
          .[,.SD,.SDcols = -c("Group_size")] },
        by = "Group",
        all.x = TRUE
  )

profile_by_group_integrated <- 
  ggplot(profile_by_group_holistic_data,aes(width = Group_proportion)) +
      geom_bar(aes(x = x_position,y = Cluster_proportion_within_group,fill = Cluster),
               stat = "identity",position = "stack",alpha = 0.7) +
      scale_x_continuous(name = "CTDs\n(n)",
                         breaks = unique(profile_by_group_holistic_data$x_position),
                         labels = unique(profile_by_group_holistic_data$Group),
                         expand = c(0.001,0.001)) +
      scale_y_continuous(name = "Percentage",
                         breaks = seq(0,1,by = 0.2),
                         labels = seq(0,100,by = 20),
                         expand = c(0,0)) +
      scale_fill_manual(values = cluster_color_set) +
      facet_grid(~ Group_label, space = "free_x", scales = "free_x",switch = "x") +
      theme_bw() +
      theme(legend.title = element_text(face = "bold",size = rel(1)),
            legend.text = element_text(face = "bold",size = rel(.9)),
            legend.position = "bottom",
            legend.key.size = unit(.6,"line"),
            axis.title.y = element_text(face = "bold",size = rel(1),margin = margin(0,.3,0,0,unit = "cm")),
            axis.title.x = element_text(face = "bold",size = rel(1),margin = margin(0,0,0,0,unit = "cm")),
            axis.text = element_text(face = "bold",size = rel(.9)),
            strip.placement = "outside",
            strip.background = element_blank(),
            strip.text = element_text(face = "bold",size = rel(.6),vjust = 3.5),
            panel.border = element_blank()
            ) +
      guides(fill = guide_legend(nrow=1,byrow=TRUE))
```

<br/>

## Combine the plots to present in a single frame - Figure 2[^3]

``` r
concatenated_clustering_results_plot <- 
  ggarrange(cluster_biplots[[1]],
            ggarrange(profile_by_cluster_integrated,
                      profile_by_group_integrated,
                      nrow = 2,labels = c("B","C")),
            ncol = 2,
            labels = c("A"),
            widths = c(1.5,1))


ggexport(concatenated_clustering_results_plot,
         filename = paste0("./Cluster_Analysis_files/Cluster_profile/",
                           "clustering_results_and_cluster_content_profile",
                           ".jpeg"),
         width = 12000,height = 8000,res = 1000,verbose = FALSE)
```

# Immunomarker profile among clusters

## Profile the proportion of status in each immunomarker by cluster

``` r
exam_items <- attr_data[,str_replace_all(Attr_mdf,"AbN_|N_","") %>% unique(.)]

profile_status_by_item_data <- 
  copy(PC_data) %>% 
  # Subset data by column names
  .[,.SD,.SDcols = c("IDCODE","Cluster","Group",exam_items)] %>% 
  # Transform data set into long format for further data portraying
  melt(.,
       measure.vars = exam_items,
       variable.name = "Exam_item",
       value.name = "Status") %>%
  # Modify results values by removing Exam_item strings in values
  .[,Status := str_replace_all(Status,
                               paste(c("_",exam_items),collapse = "|"),
                               "")] %>%
  # Tabulation
  .[,.N,by = c("Cluster","Exam_item","Status")] %>%
  # Transform data set for tabulating every possible permutation
  dcast.data.table(.,... ~ Status,value.var = "N",fill = 0) %>% 
  melt.data.table(.,id.vars = c("Cluster","Exam_item"),variable.name = "Status",value.name = "N") %>% 
  # Calculation for proportions
  .[,proportion := round(N/sum(N),3),by = c("Cluster","Exam_item")] %>%
  # Subset data set for those rows of Abnormal Status
  .[Status %in% "AbN",] %>%
  # Data ordering
  .[order(Cluster,-proportion),]
```

<br/>

## Prepare plot data for immunomarker profile

``` r
radar_plot_item_order <- 
  split(profile_status_by_item_data,by = "Cluster") %>% 
  map(., ~ .x[order(-proportion),] %>%
           .[proportion > 0.5,as.character(Exam_item)] ) %>%
  flatten_chr(.) %>% 
  unique(.) %>% 
  # set the above variables in the front of other variables
  { c(.,str_subset(levels(profile_status_by_item_data$Exam_item),
                   paste(.,collapse = "|"),
                   negate = TRUE)) } %>% 
  str_subset(.,
             split(profile_status_by_item_data,by = "Exam_item") %>% 
              # remove exam items in all clusters are less than 10%
              map(., ~ .x[all(proportion<=0.1),unique(as.character(Exam_item))]) %>% 
              flatten_chr(.) %>% 
              paste(.,collapse = "|"),
             negate = TRUE)

customised_rose_plot_item_order <-
  c("ANA", "C3", "C4", "Anti-dsDNA",
    "CRYOID", "ACAG", "Anti-RNP-Sm", "Anti-SSB", "Anti-SSA", "RF")

profile_status_by_item_plot_data <-  
  copy(profile_status_by_item_data) %>% 
  .[,.SD,.SDcols = c("Cluster","Exam_item","proportion")] %>% 
  .[Exam_item %in% customised_rose_plot_item_order,] %>% 
  .[,Exam_item := factor(Exam_item,levels = customised_rose_plot_item_order)] %>%
  .[,Cluster := factor(Cluster,
                       levels = levels(Cluster),
                       labels = map_chr(levels(Cluster),
                                        ~ str_c("Cluster ",.x)))] %>% 
    .[order(Cluster,Exam_item),] %>% 
  dcast.data.table(.,... ~ Exam_item,value.var = "proportion")
```

<br/>

## Profile immunomarker status by cluster - overlapped

``` r
# A function for change the colour transparency of a given colour set
t_col <- 
  function(color, percent = percent, name = NULL) {
    ## Get RGB values for named colour
    rgb.val <- col2rgb(color)
  
    ## Make new colour using input colour as base and alpha set by transparency
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                 max = 255,
                 alpha = (100 - percent) * 255 / 100,
                 names = name)
    
    ## Save the colour
    invisible(t.col)
}

profile_status_by_item_overlapped_plot <-
  ggRadar2(data = profile_status_by_item_plot_data,
           aes(group = Cluster),
           alpha = 0.05,
           size = .1,
           rescale = FALSE) +
    scale_color_manual(values = map_chr(cluster_color_set, ~ t_col(.x,percent = 40)),
                       labels = c(profile_status_by_item_plot_data[,str_replace_all(Cluster,"Cluster\\s","")])) +
    scale_fill_manual(values = map_chr(cluster_color_set, ~ t_col(.x,percent = 40)),
                      labels = c(profile_status_by_item_plot_data[,str_replace_all(Cluster,"Cluster\\s","")])) +
    scale_y_continuous(breaks = seq(0,1,by = .25),
                       labels = paste0(seq(0,100,by = 25),"%"),
                       limits = seq(0,1)) +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.background = element_rect(fill = "white"),
          panel.grid = element_line(size = 0.5,linetype = 2),
          legend.position = "bottom",
          legend.title = element_text(face = "bold",size = rel(1.2)),
          legend.text = element_text(face = "bold",size = rel(1.1)),
          axis.text.x = element_text(face = "bold",size = rel(1.1),margin = margin(10,0,0,0)),
          axis.text.y = element_text(face = "bold",size = rel(1.1),margin = margin(0,10,0,0)),
          axis.ticks.y.left = element_blank(),
          plot.background = element_rect(fill = "white", color = "white")) +
    guides(color = guide_legend(nrow = 1,byrow = TRUE))
```

<br/>

<br/>

## Profile immunomarker status by cluster - separated

``` r
profile_status_by_item_roseplot <-
  copy(profile_status_by_item_plot_data) %>% 
  melt.data.table(data = .,
                  id.vars = "Cluster",
                  variable.name = "Biomarker",
                  variable.factor = TRUE,
                  value.name = "Proportion") %>% 
  { ggplot(data = .) +
      geom_bar(aes(x = Biomarker, y = Proportion, fill = Cluster), stat = "identity", colour = "grey20", size = .3) +
      scale_fill_manual(values = map_chr(cluster_color_set, ~ t_col(.x,percent = 40))) +
      scale_y_continuous(breaks = seq(0,1,by = .25),
                         labels = paste0(seq(0,100,by = 25),"%"),
                         limits = seq(0,1)) +
      theme_bw() +
      coord_polar(theta = "x", start = -pi/45, clip = "off") +
      theme(panel.border = element_blank(),
            panel.background = element_rect(fill = "white"),
            panel.grid = element_line(size = 0.5,linetype = 2),
            # panel.spacing = unit(.2, "lines"),
            legend.position = "none",
            axis.title = element_blank(),
            axis.text.x = element_text(face = "bold",size = rel(1),margin = margin(10,0,0,0)),
            axis.text.y = element_text(face = "bold",size = rel(1),margin = margin(0,10,0,0)),
            axis.ticks.y.left = element_blank(),
            strip.text = element_text(face = "bold",size = rel(1.2)),
            strip.background = element_blank(),
            plot.background = element_rect(fill = "white", color = "white")) +
      facet_wrap(Cluster ~ .,nrow = 2)
      
  }
```

<br/>

<br/>

## Calculation of sensitivity and specificity for assigning to clusters

### Data preps

``` r
# --------------------------------------------------------------------------------------------------------------
sen_spe_data_prep <- 
  copy(PC_data) %>% 
  .[,c(exam_items) := map(.SD,
                          ~ str_split(.x,"_",simplify = TRUE)[,1] %>% 
                            factor(.,
                                   levels = c("N","AbN"),
                                   labels = c("0","1"))),
    .SDcols = exam_items]

# --------------------------------------------------------------------------------------------------------------
sen_spe_data <- 
  PC_data[,unique(Cluster) %>% seq_along(.)] %>% 
  map(.,
      ~ copy(sen_spe_data_prep) %>% 
        .[,Truth := ifelse(Cluster==.x,1,0) %>% factor(.,levels = c(0,1))] %>% 
        .[,.SD,.SDcols = c("Truth",exam_items)] %>% 
        melt.data.table(., id.vars = "Truth", variable.factor = FALSE, value.name = "Test" ,value.factor = TRUE)
      )

# --------------------------------------------------------------------------------------------------------------
sens_result_data <- 
  map(sen_spe_data,
      ~ copy(.x) %>% 
        .[,sens_vec(truth = Truth, estimate = Test, event_level = "second"), by = "variable"] %>% 
        .[,.metric := "sens"] %>% 
        setnames(.,"V1",".estimate")
      )

spec_result_data <- 
  map(sen_spe_data,
      ~ .x[,spec_vec(truth = Truth, estimate = Test, event_level = "second"), by = "variable"] %>% 
        .[,.metric := "spec"] %>% 
        setnames(.,"V1",".estimate")
      )

sens_spec_result_data <- 
  map2(sens_result_data,
       spec_result_data,
       ~ rbind(.x,.y) %>% 
         setDT(.) %>% 
         .[,.metric := factor(.metric, levels = c("spec","sens"), labels = c("Specificity","Sensitivity"))] %>% 
         .[, cumulative_estimate := cumsum(.estimate), by = c("variable")]
         # dcast.data.table(., variable ~ .metric, value.var = ".estimate") %>% 
         # .[, sum_sens_spec := rowSums(.SD), .SDcols = c("sens","spec")]
       ) %>% 
    setNames(.,paste(PC_data[,unique(Cluster) %>% seq_along(.)]) ) %>%
    rbindlist(.,idcol = "Cluster") %>% 
  .[variable %in% customised_rose_plot_item_order,] %>% 
  .[,variable := factor(variable,levels = rev(customised_rose_plot_item_order))] %>% 
  .[,Cluster := factor(Cluster,levels = PC_data[,levels(Cluster) %>% rev(.)])] %>%
  # melt.data.table(.,id.vars = c("Cluster","variable"), variable.name = "metric", value.name = "estimate") %>% 
  .[order(Cluster, .metric),]
```

<br/>

## Visualise the results

``` r
sens_spec_result_plot <- 
  copy(sens_spec_result_data) %>% 
  { ggplot(.,
           aes(x = variable, y = cumulative_estimate, fill = Cluster)) +
      geom_col(data = .[.metric=="Specificity"], width = .65, position = position_dodge(width = .7), alpha = .5) +
      geom_col(data = .[.metric=="Sensitivity"], width = .65, position = position_dodge(width = .7), alpha = 1) +
      geom_tile(aes(x = variable, y = NA_integer_, alpha = .metric)) +
      geom_hline(yintercept = 1,   size = .5, color = "grey40", linetype = 2) +
      geom_hline(yintercept = 1.5, size = .5, color = "grey20", linetype = 2) +
      scale_x_discrete(name = "Immunomarkers") +
      scale_alpha_discrete(name = "Metric", breaks = c("Sensitivity","Specificity")) +
      scale_y_continuous(name = "Sensitivity + Specificity",expand = c(0,0),limits = c(0,2)) +
      scale_fill_manual(values = rev(cluster_color_set)) +
      theme_bw()  +
      guides(fill = guide_legend(reverse = TRUE)) +
    theme(axis.title.x = element_text(face = "bold",size = rel(1.2),margin = margin(20,0,0,0)),
          axis.title.y = element_text(face = "bold",size = rel(1.2),margin = margin(0,20,0,0)),
          axis.text = element_text(face = "bold",size = rel(1)),
          legend.title = element_text(face = "bold",size = rel(1.2)),
          legend.text = element_text(face = "bold",size = rel(1))) +
    coord_flip() }
```

<br/>

## Combine the plots to present in a single frame - Figure 3[^4]

``` r
concatenated_immunomarkers_profile_plot <- 
  ggarrange(ggarrange(profile_status_by_item_roseplot,
                      profile_status_by_item_overlapped_plot,
                      nrow = 1, 
                      ncol = 2,
                      widths = c(3,2.05),
                      labels = c("A","B")),
            sens_spec_result_plot,
            labels = c("","C"),
            nrow = 2,
            ncol = 1)

ggexport(concatenated_immunomarkers_profile_plot,
         filename = paste0("./Cluster_Analysis_files/Cluster_profile/",
                           "immunomarker_profile",
                           ".jpeg"),
         width = 20000,height = 16000,res = 1200,verbose = FALSE)
```

<br/>

Next section:
[Clinical_Implication_Analysis.md](https://github.com/DHLab-TSENG/Heterogeneity-in-SCTD-paper/blob/main/Clinical_Implication_Analysis.md)

[^1]: Data availability in this repository is restricted due to the
    regulation of the law on the protection of patients’ data

[^2]: Please refer to the manuscript for **Table 1**

[^3]: Please refer to the manuscript for **Figure 2**

[^4]: Please refer to the manuscript for **Figure 3**
