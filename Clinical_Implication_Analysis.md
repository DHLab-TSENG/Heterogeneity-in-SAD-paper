**Identifying heterogeneous subgroups of systemic autoimmune diseases by
applying a joint dimension reduction and clustering approach to
immunomarkers : a retrospective study** <br/>- clinical implication
analysis
================
**Chia-Wei Chang<sup>a,#</sup>, Hsin-Yao Wang<sup>b,#</sup>, Wan-Ying
Lin<sup>c</sup>, Yu-Chiang Wang<sup>d</sup>, Wei-Lin Lo<sup>e</sup>,
Wei-Ting Lin<sup>b</sup>, Jia-Ruei Yu<sup>b</sup>, Yi-Ju
Tseng<sup>a,d,\*</sup>** <br/> <sup>a</sup> Department of Computer
Science, National Yang Ming Chiao Tung University, Hsinchu, Taiwan <br/>
<sup>b</sup> Department of Laxboratory Medicine, Chang Gung Memorial
Hospital at Linkou, Taoyuan City, Taiwan <br/> <sup>c</sup> Syu Kang
Sport Clinic, Taipei, Taiwan <br/> <sup>d</sup> d Department of
Medicine, Brigham and Women’s Hospital, Boston, USA <br/> <sup>e</sup>
Department of Rheumatology, Chang Gung Memorial Hospital at Keelung,
Keelung City, Taiwan <br/> <sup>f</sup> Computational Health Informatics
Program, Boston Children’s Hospital, Boston, MA, USA <br/> <sup>\#</sup>
Chang and Wang contribute equally to this work <br/> <sup>\*</sup>
Corresponding Author <br/>

------------------------------------------------------------------------

# Set up for environment

## Load libraries

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

                "ggplot2",
                "ggiraphExtra",
                "ggsci",
                "ggrepel",
                "ggforce",
                "ggpubr",
                "cowplot",
                "pheatmap",
                "viridis",
                "ggplotify",
                "gtable",
                
                "tableone",
                "flextable",
                "webshot"
                
                )

if (any(!.cran_pkgs %in% installed.packages())) {
  install.packages(.cran_pkgs[!.cran_pkgs %in% installed.packages()],
                   lib = Sys.getenv("R_LIBS_USER"),
                   dependencies = TRUE)
}

sapply(.cran_pkgs,
       function(x) suppressPackageStartupMessages(require(x,character.only = TRUE)))
```

    ##   data.table        dplyr     magrittr      stringr    lubridate       scales 
    ##         TRUE         TRUE         TRUE         TRUE         TRUE         TRUE 
    ##        tidyr        purrr        furrr     devtools        knitr      modeest 
    ##         TRUE         TRUE         TRUE         TRUE         TRUE         TRUE 
    ##          vcd        psych    DescTools      clustrd      ggplot2 ggiraphExtra 
    ##         TRUE         TRUE         TRUE         TRUE         TRUE         TRUE 
    ##        ggsci      ggrepel      ggforce       ggpubr      cowplot     pheatmap 
    ##         TRUE         TRUE         TRUE         TRUE         TRUE         TRUE 
    ##      viridis    ggplotify       gtable     tableone    flextable      webshot 
    ##         TRUE         TRUE         TRUE         TRUE         TRUE         TRUE

<br/>

## Install packages from GitHub

``` r
.ghub_pkgs <- c("DHLab-TSENG/dxpr", "ddsjoberg/gtsummary","rstudio/gt")

if (any(!sapply(strsplit(.ghub_pkgs,"/"),`[`,2) %in% installed.packages())){
  library(devtools)
  install_github(.ghub_pkgs[!sapply(strsplit(.ghub_pkgs,"/"),`[`,2) %in% installed.packages()],
                 lib = Sys.getenv("R_LIBS_USER"),
                 dependencies = TRUE,
                 force = TRUE)
}

sapply(sapply(strsplit(.ghub_pkgs,"/"),`[`,2),
       function(x) suppressPackageStartupMessages(require(x,character.only = TRUE)))
```

    ##      dxpr gtsummary        gt 
    ##      TRUE      TRUE      TRUE

<br/>

# Import datasets[^1]

From
[Data_Wrangling.md](https://github.com/DHLab-TSENG/Heterogeneity-in-SAD-paper/blob/main/Data_Wrangling.md)

``` r
raw_diagnosis_data    <- readRDS("./Dataset/D1_1_concatenated_diagnosis_data(dt_query_dx_data).rds")
CTD_diagnosis_data    <- readRDS("./Dataset/A1_1_diagnosis_eachDate_data.rds")
CTD_diagnosis_summary <- readRDS("./Dataset/A1_2_diagnosis_Count_eachDx_data.rds")
```

From
[Cluster_Analysis.md](https://github.com/DHLab-TSENG/Heterogeneity-in-SAD-paper/blob/main/Cluster_Analysis.md)

``` r
PC_data   <- readRDS("./Dataset/PC_data.rds")
attr_data <- readRDS("./Dataset/attr_data.rds")
```

------------------------------------------------------------------------

# Preprocess data

## Import diagnosis database and its profiling datasets

``` r
dx_database <- 
  raw_diagnosis_data %>% 
  # transform column types
  .[,DSSSQNO := as.character(DSSSQNO)] %>% 
  # Extract study subjects
  .[IDCODE %in% PC_data$IDCODE]

# Extract study subjects
CustomGrepGroup_groupedDT<- 
  CTD_diagnosis_data[ID %in% PC_data$IDCODE]

# Extract study subjects
CustomGrepGroup_summarised_groupedDT <- 
  CTD_diagnosis_summary[ID %in% PC_data$IDCODE]
```

<br/>

## Subset subjects’ data by follow-up period thresholds

``` r
Followup_threshold <- 3 # at least for three years

CustomGrepGroup_summarised_groupedDT_sufficient_followup <- 
  CustomGrepGroup_summarised_groupedDT[,as.numeric(period/365.25) >= Followup_threshold] %>% 
  { CustomGrepGroup_summarised_groupedDT[.,] }
```

<br/>

## Profile number of remaining subjects after follow-up period filter by SAD group

``` r
Sufficient_followup_profile <- 
  copy(CustomGrepGroup_summarised_groupedDT_sufficient_followup) %>% 
    .[,Sufficient_N := uniqueN(ID),by = c("Group")] %>% 
    .[,unique(.SD),.SDcols = c("Group","Sufficient_N")] %>% 
    .[order(Group),] %>% 
    merge(.,
          PC_data[,.N,by = "Group"],
          by = "Group") %>% 
    .[,Sufficient_proportion := round(Sufficient_N/N*100,1)]

print(Sufficient_followup_profile)
```

    ##                           Group Sufficient_N    N Sufficient_proportion
    ## 1:         Rheumatoid arthritis         3808 5696                  66.9
    ## 2:           Sjogren's syndrome         2043 3541                  57.7
    ## 3: Systemic lupus erythematosus         1887 2686                  70.3

<br/>

------------------------------------------------------------------------

# Comorbidities identified by CCS

## Convert ICD into CCS, retaining data within the given time range

``` r
# Query CCS by Diagnosis ICD-code
correction_table <- 
  icdDxToCCS(dxDataFile = dx_database,
             idColName = IDCODE,
             icdColName = mod_DSSID,
             dateColName = IPDAT,
             isDescription = FALSE,
             icd10usingDate = "2016-01-01") %>% 
  pluck(.,"Error") %>% 
  .[Suggestion!="",unique(.SD),.SDcols = c("ICD","Suggestion")]

CCSLong_follow_window <- 
  merge(dx_database,
        correction_table,
        by.x = "mod_DSSID",
        by.y = "ICD",
        all.x = TRUE) %>% 
  # Replace original ICD-code with suggestion one (correct one)
  .[!is.na(Suggestion),mod_DSSID := Suggestion] %>% 
  # Query CCS with correct ICD-code
  icdDxToCCS(dxDataFile = .,
             idColName = IDCODE,
             icdColName = mod_DSSID,
             dateColName = IPDAT,
             isDescription = TRUE,
             icd10usingDate = "2016-01-01")

# Extract CCS data where the original ICD-codes within given time range
plan(multisession, workers = 5)

CCS_within_DateRange <- 
  future_pmap_dfr(list(CustomGrepGroup_summarised_groupedDT_sufficient_followup$ID,
                       CustomGrepGroup_summarised_groupedDT_sufficient_followup$firstCaseDate,
                       CustomGrepGroup_summarised_groupedDT_sufficient_followup$endCaseDate),
                  ~ CCSLong_follow_window$groupedDT %>% 
                    .[ID %in% ..1 & Date %between% c(..2,..3),])

plan(sequential)
```

<br/>

## Generate CCS reference database by transforming website data derived from HCUP

``` r
CCS_catagories <- 
  system("curl -s https://hcup-us.ahrq.gov/toolssoftware/ccs/AppendixASingleDX.txt |
          awk 'NF && /^[1-9]/ {print $0}' ",
         intern = TRUE)

CCS_bundle <- 
  system("curl -s https://hcup-us.ahrq.gov/toolssoftware/ccs/AppendixASingleDX.txt |
          awk 'NR>4 && !/^[1-9#]/ {print $0}' ",
         intern = TRUE) %>% 
  str_trim(.,"both") %>%
  data.table(CCS = .)

split_index <-   
  copy(CCS_bundle) %>% 
  .[,.(which(str_length(CCS)==0))] %>% 
  { .[,c(c(.[1,V1], 
           V1 - lag(V1))[-2], 
         .[,nrow(CCS_bundle)-tail(V1,1)]) ] } %>% 
  map2(seq(1:length(CCS_catagories)),
       .,
       ~ rep(.x,.y)) %>% 
  flatten(.) %>% 
  data.table(CCS_Group_Index = as.character(.))

CCS_reference <- 
  cbind(CCS_bundle,split_index) %>% 
  split(.,by = "CCS_Group_Index") %>% 
  map(., ~ .x[,head(.SD,-1)]) %>%
  setNames(.,CCS_catagories) %>% 
  map(.,
      ~ .[,str_flatten(CCS," ") %>% str_split(.," ")] ) %>% 
  rbindlist(.,idcol = "CCS_category") %>% 
  setnames(.,"V1","CCS_code")
```

<br/>

## Summarise unmatched ICDs and resolve by modifying ICD codes

``` r
ambiguous_ICD <- 
  copy(CCSLong_follow_window$groupedDT) %>% 
  .[is.na(CCS_CATEGORY_DESCRIPTION),] %>% 
  .[,':=' (Year = year(Date),
          `Befor_2016-01-01` = ifelse(ymd(Date)<ymd("2016-01-01"),"TRUE","FALSE") )] %>% 
  .[,.N,by = c("ICD","Befor_2016-01-01")] %>% 
  .[order(`Befor_2016-01-01`,-N),] %>% 
  # Look up for corresponding CCS and matched patterns using ICDs
  .[,':=' (prefix_0 = map_lgl(ICD, ~ CCS_reference[str_detect(CCS_code,paste0("^[0]+",.x,"$")),.N]>0),
           suffix_0 = map_lgl(ICD, ~ CCS_reference[str_detect(CCS_code,paste0("^",.x,"[0]+[0-9]*")),.N]>0),
           suffix_1 = map_lgl(ICD, ~ CCS_reference[str_detect(CCS_code,paste0("^",.x,"[1]+[0-9]*")),.N]>0),
           prefix_0_ptrn = map_chr(ICD, ~ CCS_reference[str_detect(CCS_code,paste0("^[0]+",.x,"$")),
                                                        str_c(CCS_code,collapse = " | ")]),
           suffix_0_ptrn = map_chr(ICD, ~ CCS_reference[str_detect(CCS_code,paste0("^",.x,"[0]+[0-9]*")),
                                                        str_c(CCS_code,collapse = " | ")]),
           suffix_1_ptrn = map_chr(ICD, ~ CCS_reference[str_detect(CCS_code,paste0("^",.x,"[1]+[0-9]*")),
                                                        str_c(CCS_code,collapse = " | ")])
           )]

correction_table_padding0or1 <- 
  copy(ambiguous_ICD) %>% 
  .[prefix_0==TRUE & str_detect(ICD,"^[8,9]"),
    replacement := map_chr(prefix_0_ptrn, ~ str_split(.x, "\\ \\| ") %>%
                                            map_chr(.,~ tail(.x,1)) )] %>%
  .[prefix_0==TRUE & str_detect(ICD,"^[^8,9]"),
    replacement := map_chr(prefix_0_ptrn, ~ str_split(.x, "\\ \\| ") %>% 
                                            map_chr(.,~ head(.x,1)) )] %>% 
  # for 345, 345, 535 to xxx00
  .[is.na(replacement) & suffix_0==TRUE & str_length(ICD)==3 &
    map_lgl(suffix_0_ptrn, ~ { str_split(.x,"\\ \\| ") %>%
                                map_chr(., ~ head(.x,1)) %>%
                                map(., ~ str_length(.x)) == 4 } ) &
    map_lgl(suffix_0_ptrn, ~ { str_split(.x,"\\ \\| ") %>% 
                                map_int(., ~ length(.x)) > 1 } ),
    replacement := map_chr(ICD, ~ str_pad(.x,5,side = "right",pad = "0") )] %>% 
  # for 946, 116 to xxx0
  .[is.na(replacement) & suffix_0==TRUE & str_length(ICD)==3 &
    map_lgl(suffix_0_ptrn, ~ { str_split(.x,"\\ \\| ") %>%
                                map_chr(., ~ head(.x,1)) %>%
                                map(., ~ str_length(.x)) == 4 } ) &
    map_lgl(suffix_0_ptrn, ~ { str_split(.x,"\\ \\| ") %>% 
                                map_int(., ~ length(.x)) == 1 } ),
    replacement := map_chr(ICD, ~ str_pad(.x,4,side = "right",pad = "0") )] %>% 
  # for 295, 532, 531, 493 and others to xxx00
  .[is.na(replacement) & suffix_0==TRUE,
    replacement := map_chr(suffix_0_ptrn, ~ { str_split(.x,"\\ \\|") %>% 
                                               map_chr(., ~ head(.x,1)) } )] %>% 
  .[is.na(replacement) & suffix_1==TRUE,
    replacement := map_chr(suffix_1_ptrn, ~ { str_split(.x,"\\ \\|") %>% 
                                               map_chr(., ~ head(.x,1)) } )] %>% 
  # remove unfixable data
  .[!is.na(replacement),]
```

<br/>

## Rerun conversion from ICD to CCS using modified (corrected) ICD data

``` r
CCSLong_follow_window_final <- 
  merge(dx_database,
        correction_table,
        by.x = "mod_DSSID",
        by.y = "ICD",
        all.x = TRUE) %>% 
  .[!is.na(Suggestion),mod_DSSID := Suggestion] %>% 
  merge(.,
        correction_table_padding0or1[,unique(.SD),.SDcols = c("ICD","replacement")],
        by.x = "mod_DSSID",
        by.y = "ICD",
        all.x = TRUE) %>% 
  .[!is.na(replacement),mod_DSSID := replacement] %>% 
  icdDxToCCS(dxDataFile = .,
             idColName = IDCODE,
             icdColName = mod_DSSID,
             dateColName = IPDAT,
             isDescription = TRUE,
             icd10usingDate = "2016-01-01")

# Extract CCS data where the original ICD-codes within given time range
plan(multisession, workers = 5)

CCS_within_DateRange <- 
  future_pmap_dfr(list(CustomGrepGroup_summarised_groupedDT_sufficient_followup$ID,
                       CustomGrepGroup_summarised_groupedDT_sufficient_followup$firstCaseDate,
                       CustomGrepGroup_summarised_groupedDT_sufficient_followup$endCaseDate),
                  ~ CCSLong_follow_window$groupedDT %>% 
                    .[ID %in% ..1 & Date %between% c(..2,..3),])

plan(sequential)
```

<br/>

## Merge clustering results with CCS data set and generate profile data set for clinical implications analysis

``` r
removed_CCS <- c("Residual codes; unclassified",
                 "Other connective tissue disease",
                 "External cause codes")

CCS_within_DateRange_Cluster_Group <-
  copy(CCS_within_DateRange) %>%
  merge(.,
        PC_data[,.SD,.SDcols = c("IDCODE","Group","Cluster")],
        by.x = "ID",
        by.y = "IDCODE",
        all.x = TRUE)

CCS_within_DateRange_Cluster_Group_profile_data <- 
  copy(CCS_within_DateRange_Cluster_Group) %>% 
  .[str_detect(CCS_CATEGORY_DESCRIPTION,paste(removed_CCS,collapse = "|"),negate = TRUE),] %>% 
  .[,Cluster := paste("Cluster",Cluster)] %>%
  .[,unique(.SD),.SDcols = c("ID","CCS_CATEGORY_DESCRIPTION","Cluster")] %>% 
  # Modify description strings for an easy-to-read purpose
  .[CCS_CATEGORY_DESCRIPTION %in% "Inflammation; infection of eye (except that caused by tuberculosis or sexually transmitteddisease)",
    CCS_CATEGORY_DESCRIPTION := "Inflammation; infection of eye"] %>% 
  .[CCS_CATEGORY_DESCRIPTION %in% "Pneumonia (except that caused by tuberculosis or sexually transmitted disease)",
    CCS_CATEGORY_DESCRIPTION := "Pneumonia"] %>% 
  .[CCS_CATEGORY_DESCRIPTION %in% "Hodgkin`s disease",
    CCS_CATEGORY_DESCRIPTION := "Hodgkin's disease"] %>% 
  .[CCS_CATEGORY_DESCRIPTION %in% "Non-Hodgkin`s lymphoma",
    CCS_CATEGORY_DESCRIPTION := "Non-Hodgkin's lymphoma"] %>% 
  .[CCS_CATEGORY_DESCRIPTION %in% "Parkinson`s disease",
    CCS_CATEGORY_DESCRIPTION := "Parkinson's disease"]
```

<br/>

## Positive proportions in each CCS by Cluster

``` r
CCS_retain_variables <- 
  copy(CCS_within_DateRange_Cluster_Group_profile_data) %>% 
  .[,Cluster_size := uniqueN(ID),by = c("Cluster")] %>% 
  .[!is.na(CCS_CATEGORY_DESCRIPTION),.N,by = c("Cluster","CCS_CATEGORY_DESCRIPTION","Cluster_size")] %>% 
  .[,prop := round(N/Cluster_size*100,2)] %>% 
  dcast.data.table(.,CCS_CATEGORY_DESCRIPTION ~ Cluster,value.var = "prop",fill = 0) %>% 
  # Threshold setting
  { .[rowSums(. < 25) < n_distinct(CCS_within_DateRange_Cluster_Group_profile_data$Cluster)] } %>%
  .[,CCS_CATEGORY_DESCRIPTION]


table_one_CCS_data <- 
  dcast.data.table(
    data = CCS_within_DateRange_Cluster_Group_profile_data[CCS_CATEGORY_DESCRIPTION %in% CCS_retain_variables,],
    formula = ID + Cluster ~ CCS_CATEGORY_DESCRIPTION,
    fun.aggregate = length,
    fill = "0")

table_one_CCS <- 
  CreateTableOne(data = table_one_CCS_data,
                 strata = "Cluster",
                 vars = CCS_retain_variables,
                 factorVars = CCS_retain_variables,
                 test = TRUE,
                 includeNA = FALSE)
```

<br/>

## Prepare plot data for displaying positive rate of each CCS category by cluster

``` r
comorbidity_CCS_proportion_byGroup_data <- 
    merge(table_one_CCS_data,
          PC_data[,unique(.SD),.SDcols = c("IDCODE","Group")],
          by.x = "ID",
          by.y = "IDCODE",
          all.x = TRUE) %>% 
    melt.data.table(data = .,
                    id.vars = c("ID","Group","Cluster"),
                    variable.name = "Comorbidity",
                    value.name = "Affected",
                    variable.factor = FALSE,
                    value.factor = TRUE) %>% 
    { merge(.[,n_distinct(ID),by = c("Cluster","Group")],
            .[,sum(Affected),by = c("Cluster","Group","Comorbidity")],
            by = c("Cluster","Group"),
            all.y = TRUE)
      } %>% 
    setnames(.,c("V1.x","V1.y"),c("Cluster_size","Number_of_affected")) %>% 
    # Frequency computations
    .[,Proportion_of_affected := round(Number_of_affected/Cluster_size,3)*100,
      by = c("Cluster","Group","Comorbidity")] %>% 
    # Subset by variable
    .[,.SD,.SDcols = c("Cluster","Group","Comorbidity","Proportion_of_affected")] %>% 
    # Transform dataset into wide format
    dcast.data.table(., ... ~ Cluster,value.var = "Proportion_of_affected") %>% 
    .[order(Group,-`Cluster 1`),] %>% 
    .[,Comorbidity := factor(Comorbidity,levels = unique(Comorbidity))] %>% 
    # Split by SAD group
    split(.,by = "Group",keep.by = FALSE) %>% 
    # Transform format of dataset
    map(., ~ as.matrix(.x,rownames = "Comorbidity"))
```

<br/>

------------------------------------------------------------------------

# Comorbidities identified by ICD-code

## Enumerate comorbilities and corresponding ICD-codes in regex forms

``` r
grepTable <- 
  data.table(grepIcd = c("^714[0,1,2,3.*,4,89,9]+|^M13[0,1.{1,},8{1,}]+",
                         "^7291|^M791.*",
                         "^6927.*|^69282|^L589",
                         "^7040[0,1,9]+|^L63.*",
                         "^5280[0,9]+|^528[2,9]+|^K12[0,1,30,39]+|^K137[0,9]+",
                         
                         "^6954|^L93.*",
                         "^4430|^I7300",
                         "^4438[1,9]+|^4439|^446.*|^447[6,8,9]+|^I79.*|^I7389|^I739|^M30[0,1,8]+|^M31[0-7]+|^I77[6,9]+|^I7789",
                         "^345[0.*,1.*,2,3,4.*,5.*,8.*,9.*]+|^G40[1.*,2.*,3.*,4.*,5.*,8.*,9.*]+",
                         "^420[0,9.*]+|^56789|^5679|^M3212|^K65[8,9]+|^I31[8,9]+",
                         
                         "^58[3,4,5,6]+|^996[73,81]+|^V420|^V451|^N18.*|^N19.*|^T828.*|^Z940|^Z992|^Z9115",
                         "^7142|^M0524.*",
                         "^7140|^M063[1,2,3,4,5,6,7,9]+",
                         "^7931[1,9]+|^M0630|^R91[1,8]+",
                         "^3571|^M055.*",
                         
                         "^51189|^5119|^J90|^J918",
                         "^4168|^I2721",
                         "^41519|^I26[09,93,94,99]+",
                         "^2646|^E506|^37515|^H0412[1-3,9]|M3501",
                         "^5277|^K117",
                         
                         "^7856|^R59.*",
                         "^5271|^K111",
                         "^708[1,8,9]+|^L50[1,8,9]+",
                         "^7872.*|^R13.*",
                         "^7194[0-9]*|^M255[0-7,9]*[1,2,9]*"
                         ),
             Group    = c("Arthritis",
                          "Myalgia",
                          "Photosensitivity",
                          "Alopecia",
                          "Oral ulcers",
                          
                          "Discoid lupus",
                          "Raynaud's phenomenon",
                          "Vasculitis",
                          "Seizures",
                          "Serositis",
                          
                          "Renal compromise",
                          "Digital vasculitis",
                          "Skin nodulosis",
                          "Pulmonary nodulosis",
                          "Peripheral nerve involvement",
                          
                          "Pleural effusion",
                          "PAH",
                          "Pulmonary embolism",
                          "Xerophthalmia",
                          "Xerostomia",
                          
                          "Lymphadenopathy",
                          "Parotid enlargement",
                          "Urticaria",
                          "Dysphagia",
                          "Arthralgia"
                          )
             )
```

<br/>

## Convert ICD-code into comorbidities and check that the regex patterns correctly extract the disease groups

``` r
CustomGrepComorbidity <-
  map(grepTable$Group,
      ~ icdDxToCustomGrep(dxDataFile  = dx_database,
                          idColName   = IDCODE,
                          icdColName  = mod_DSSID,
                          dateColName = IPDAT,
                          customGroupingTable = grepTable[Group %in% .x,])) %>% 
  { list(map(., ~ .x$groupedDT %>% .[!is.na(GrepedGroup),]) %>% rbindlist,
         map(., ~ .x$summarised_groupedDT) %>% rbindlist) } %>% 
  setNames(.,c("groupedDT","summarised_groupedDT")) %>% 
  map(.,~ setnames(.x,"GrepedGroup","Comorbidity",skip_absent = TRUE))


CustomGrepComorbidity$groupedDT[,unique(ICD),by = "Comorbidity"] %>% 
  setnames(.,c("Comorbidity","V1"),c("Comorbidity","ICD")) %>% 
  .[order(Comorbidity,ICD)]
```

<br/>

## Retain data within the given time range

``` r
plan(multisession, workers = 5)

ICD_within_DateRange <- 
  future_pmap_dfr(list(CustomGrepGroup_summarised_groupedDT_sufficient_followup$ID,
                       CustomGrepGroup_summarised_groupedDT_sufficient_followup$firstCaseDate,
                       CustomGrepGroup_summarised_groupedDT_sufficient_followup$endCaseDate),
                  ~ copy(CustomGrepComorbidity$groupedDT) %>% 
                    .[ID %in% ..1 & Date %between% c(..2,..3),])

plan(sequential)
```

<br/>

## Merge clustering results with ICD dataset and generate profile dataset for clinical implication analysis

``` r
ICD_within_DateRange_Cluster_Group <-
  merge(ICD_within_DateRange,
        PC_data[,.SD,.SDcols = c("IDCODE","Group","Cluster")],
        by.x = "ID",
        by.y = "IDCODE",
        all.x = TRUE)

ICD_within_DateRange_Cluster_Group_profile_data <- 
  copy(ICD_within_DateRange_Cluster_Group) %>% 
  .[,Cluster := paste("Cluster",Cluster)] %>%
  .[,unique(.SD),.SDcols = c("ID","Comorbidity","Cluster")]
```

<br/>

## Positive proportions in each comorbidity by Cluster

``` r
ICD_retain_variables <- 
  copy(ICD_within_DateRange_Cluster_Group_profile_data) %>% 
  .[,Cluster_size := uniqueN(ID),by = c("Cluster")] %>% 
  .[!is.na(Comorbidity),.N,by = c("Cluster","Comorbidity","Cluster_size")] %>% 
  .[,prop := round(N/Cluster_size*100,2)] %>% 
  dcast.data.table(.,Comorbidity ~ Cluster,value.var = "prop",fill = 0) %>% 
  .[,Comorbidity]

table_one_ICD_data <- 
  dcast.data.table(data = ICD_within_DateRange_Cluster_Group_profile_data,
                   ID + Cluster ~ Comorbidity,
                   fun.aggregate = length,
                   fill = "0")

table_one_ICD <-          
  CreateTableOne(data = table_one_ICD_data,
                 strata = "Cluster",
                 vars = ICD_retain_variables,
                 factorVars = ICD_retain_variables,
                 test = TRUE,
                 includeNA = FALSE)
```

<br/>

## Prepare plot data for displaying positive rate of each ICD category by cluster

``` r
comorbidity_ICD_proportion_byGroup_data <- 
  merge(table_one_ICD_data,
        PC_data[,unique(.SD),.SDcols = c("IDCODE","Group")],
        by.x = "ID",
        by.y = "IDCODE",
        all.x = TRUE) %>% 
  melt.data.table(data = .,
                  id.vars = c("ID","Group","Cluster"),
                  variable.name = "Comorbidity",
                  value.name = "Affected",
                  variable.factor = FALSE,
                  value.factor = TRUE) %>% 
  { merge(.[,n_distinct(ID),by = c("Cluster","Group")],
          .[,sum(Affected),by = c("Cluster","Group","Comorbidity")],
          by = c("Cluster","Group"),
          all.y = TRUE)
    } %>% 
  setnames(.,c("V1.x","V1.y"),c("Cluster_size","Number_of_affected")) %>% 
  # Frequency computations
  .[,Proportion_of_affected := round(Number_of_affected/Cluster_size,3)*100,
    by = c("Cluster","Group","Comorbidity")] %>% 
  # Subset by variable
  .[,.SD,.SDcols = c("Cluster","Group","Comorbidity","Proportion_of_affected")] %>% 
  # Transform dataset into wide format
  dcast.data.table(., ... ~ Cluster,value.var = "Proportion_of_affected") %>% 
  .[order(Group,-`Cluster 1`),] %>% 
  .[,Comorbidity := factor(Comorbidity,levels = unique(Comorbidity))] %>% 
  # Split by SAD group
  split(.,by = "Group",keep.by = FALSE) %>% 
  # Transform format of dataset
  map(., ~ as.matrix(.x,rownames = "Comorbidity"))
```

<br/>

------------------------------------------------------------------------

# Visualize positive rate of comorbidities by cluster with heatmaps

## CCS categories

``` r
paletteLength <- 20
mid_point <- 65
min_value <- 0
max_value <- 100
Color_list <- colorRampPalette(c("darkgreen", "white", "darkgoldenrod2", "firebrick3"))(paletteLength)
Break_list <- c(seq(min_value,mid_point,length.out = ceiling(paletteLength/2) + 1), 
                seq(mid_point, max_value, length.out = floor(paletteLength/2) + 1)[-1] )  
  
bold_row_labels <- map(comorbidity_CCS_proportion_byGroup_data,
                        ~ lapply(rownames(.x),function(x) bquote(bold(.(x))) ) )
bold_col_labels <- map(comorbidity_CCS_proportion_byGroup_data,
                        ~ lapply(colnames(.x),function(x) bquote(bold(.(x))) ) )

comorbidity_CCS_prevalence_byGroup_plot <-
  pmap(list(comorbidity_CCS_proportion_byGroup_data,
            bold_row_labels,
            bold_col_labels,
            names(comorbidity_CCS_proportion_byGroup_data),
            c(FALSE,FALSE,FALSE)),
      ~ pheatmap(mat = ..1,
                 color = Color_list,
                 cellwidth = 6,
                 cellheight = 6,
                 clustering_method = "average",
                 cluster_cols = FALSE,
                 breaks = Break_list,
                 fontsize = 6,
                 legend_breaks = seq(min_value,max_value,by = paletteLength/2),
                 legend_labels = seq(min_value,max_value,by = paletteLength/2),
                 labels_row = as.expression(..2),
                 labels_col = as.expression(..3),
                 main = ifelse(..4 == "Sjogren's syndrome",
                               paste(..4,
                                     str_c(rep(" ",28),collapse = "")),
                               ifelse(..4 == "Rheumatoid arthritis",
                                      paste(..4,
                                            str_c(rep(" ",28),collapse = "")),
                                      ifelse(..4 == "Systemic lupus erythematosus",
                                             paste(..4,
                                                   str_c(rep(" ",10), collapse = "")),
                                             "")
                                      )
                               ),
                 legend = ..5,
                 silent = TRUE
                 )
      )
```

<br/>

## ICD categories

``` r
paletteLength <- 20
mid_point <- 65
min_value <- 0
max_value <- 100
Color_list <- colorRampPalette(c("darkgreen", "white", "darkgoldenrod2", "firebrick3"))(paletteLength)
Break_list <- c(seq(min_value,mid_point,length.out = ceiling(paletteLength/2) + 1), 
                seq(mid_point, max_value, length.out = floor(paletteLength/2) + 1)[-1] )  
  
bold_row_labels <- map(comorbidity_ICD_proportion_byGroup_data,
                        ~ lapply(rownames(.x),function(x) bquote(bold(.(x))) ) )
bold_col_labels <- map(comorbidity_ICD_proportion_byGroup_data,
                        ~ lapply(colnames(.x),function(x) bquote(bold(.(x))) ) )

comorbidity_ICD_prevalence_byGroup_plot <-
  pmap(list(comorbidity_ICD_proportion_byGroup_data,
            bold_row_labels,
            bold_col_labels,
            names(comorbidity_ICD_proportion_byGroup_data),
            c(TRUE,FALSE,FALSE)),
      ~ pheatmap(mat = ..1,
                 color = Color_list,
                 cellwidth = 6,
                 cellheight = 6,
                 clustering_method = "average",
                 cluster_cols = FALSE,
                 breaks = Break_list,
                 fontsize = 6,
                 legend_breaks = seq(min_value,max_value,by = paletteLength/2),
                 legend_labels = seq(min_value,max_value,by = paletteLength/2),
                 labels_row = as.expression(..2),
                 labels_col = as.expression(..3),
                 main = ifelse(..4 == "Sjogren's syndrome",
                               paste(str_c(rep(" ",80),collapse = ""),
                                     "",
                                     str_c(rep(" ",30),collapse = ""),
                                     "Rate of occurrence"),
                               ifelse(..4 == "Rhematoid arthritis",
                                      paste("",
                                            str_c(rep(" ",28),collapse = "")),
                                      ifelse(..4 == "Systemic lupus erythematosus",
                                             paste("",
                                                   str_c(rep(" ",10), collapse = "")),
                                             "")
                                      )
                               ),
                 legend = ..5,
                 silent = TRUE
                 )
      )
```

<br/>

## Combine heatmaps of CCS and ICD category into a uniframe - Figure 4[^2]

``` r
concatenated_comorbidity_prevalence_byGroup_plot <- 
  cbind(rbind(comorbidity_CCS_prevalence_byGroup_plot[[1]]$gtable,
              comorbidity_CCS_prevalence_byGroup_plot[[2]]$gtable,
              comorbidity_CCS_prevalence_byGroup_plot[[3]]$gtable),
        rbind(comorbidity_ICD_prevalence_byGroup_plot[[1]]$gtable,
              comorbidity_ICD_prevalence_byGroup_plot[[2]]$gtable,
              comorbidity_ICD_prevalence_byGroup_plot[[3]]$gtable)
        ) %>% 
  as.ggplot(.) +
  ggtitle(paste(str_c(rep(" ",12),collapse = ""),
                "A. broad-spectrum screening",
                str_c(rep(" ",40),collapse = ""),
                "B. common comorbidities",
                collapse = "")) +
  theme(title = element_text(size = rel(.8), face = "bold"),
        plot.margin = margin(1.5,0,.5,0, "cm"))

ggexport(concatenated_comorbidity_prevalence_byGroup_plot,
         filename = paste0("./Clinical_Implication_Analysis_files/",
                           "concatenated_comorbidity_prevalence_byGroup_plot",
                           ".pdf"),
         # width = 8000,height = 13000,res = 1000,verbose = FALSE) # for jpeg
         width = 8, height = 13, verbose = FALSE) # for pdf

# system(
#   paste0("cp \\",
# 
#          "./Clinical_Implication_Analysis_files/",
#          "concatenated_comorbidity_prevalence_byGroup_plot",
#          ".pdf",
#          " \\",
#          "./Heterogeneity-in-SAD-paper_files/",
#          "Figure4_concatenated_comorbidity_prevalence_byGroup_plot",
#          ".pdf"
#          )
#   )
```

## Prepare datasets for profiling the proportion differences between clusters

``` r
library(gtsummary)
library(gt)

# CCS
comorbidity_CCS_proportion_test_data <- 
  merge(table_one_CCS_data,
        PC_data[,unique(.SD),.SDcols = c("IDCODE","Group")],
        by.x = "ID",
        by.y = "IDCODE",
        all.x = TRUE)

# ICD
comorbidity_ICD_proportion_test_data <- 
  merge(table_one_ICD_data,
          PC_data[,unique(.SD),.SDcols = c("IDCODE","Group")],
          by.x = "ID",
          by.y = "IDCODE",
          all.x = TRUE)

comorbidity_proportion_diff_test_data <- 
  map(list(comorbidity_CCS_proportion_test_data,
           comorbidity_ICD_proportion_test_data),
      ~ .x[, .SD, .SDcols = -c("ID")] %>% 
      .[, names(.) := map(.SD, ~ as.character(.x))]
      )
```

<br/>

## Conduct statistical test for the proportion differences in the clusters - between six groups

``` r
fisher_test_fun <- 
  function(data, variable, by, ...) {
    copy(data) %>%
    setDT() %>% 
    .[, .SD, .SDcols = c(variable, by)] %>% 
    .[, c(variable) := map(.SD,
                           ~ factor(.x, 
                                    levels = c("0","1"),
                                    labels = c("Unaffected", "Affected"))
                           ),
      .SDcols = c(variable)] %>% 
    setkeyv(., c(variable, by)) %>%
    .[CJ(get(variable), get(by), unique = TRUE), 
      .N, 
      by = .EACHI] %>% 
    dcast.data.table(.,
                     ... ~ Cluster,
                     value.var = "N",
                     fill = 0,
                     drop = c(FALSE,FALSE)
                     ) %>% 
    .[,RowSum := rowSums(.SD), .SDcols = -c(variable)] %>% 
    
    melt.data.table(., 
                    id.vars = c("RowSum",variable), 
                    variable.name = "by") %>% 
    .[RowSum == 0,
      value := fifelse(all(RowSum == 0) &
                        by == sample(head(unique(by),1),1),
                       1,
                       value)] %>% 
    dcast.data.table(., ... ~ by) %>% 
    .[order(get(variable)),] %>% 
    .[, .SD, .SDcols = -c("RowSum",variable)] %>% 
    
    fisher.test(.,
                hybrid = TRUE, 
                simulate.p.value = TRUE, 
                B = 2000) %>% 
    broom::tidy(.) %>% 
    .$`p.value`
    }


between_six_groups_tbl <-
  map(comorbidity_proportion_diff_test_data,
      ~ group_by(.x, Group) %>% 
        nest()
      ) %>% 
  map(.,
      ~ pluck(.x, "data") %>% 
        map(.,
            ~ tbl_summary(.x,
                          by = Cluster,
                          type = all_categorical() ~ "dichotomous",
                          missing = "no") %>% 
              add_stat(fns = all_categorical() ~ fisher_test_fun) %>% 
              modify_header(
                list(all_stat_cols() ~ "**{level}**  \n(n = {style_number(n)})",
                     add_stat_1 ~ "**p-value**"),
                     label ~ "**Features**"
                ) %>%
              modify_footnote(all_stat_cols() ~ NA)
            ) %>% 
        tbl_merge(.,
                  tab_spanner = str_c("**",c(unique(.x$Group)),"**")
                  ) %>% 
        as_gt() %>%
          tab_options(table.font.size = "11",
                      # table.font.names = "Times New Roman"
                      ) %>%
          text_transform(fn = function(x) { str_replace(x, "0 \\(0%\\)", "\\-") },
                         locations = cells_body(columns = all_stat_cols())
                         ) %>% 
          text_transform(fn = function(x) { str_replace(x, "\\(", "<br>\\(") },
                         locations = cells_body(columns = all_stat_cols())
                         ) %>% 
          fmt_number(columns = starts_with("add_stat"),
                     decimals = 3) %>% 
          text_transform(fn = function(x) { str_replace(x, "(0.000)", "<0.001") },
                         locations = cells_body(columns = c(add_stat_1_1, add_stat_1_2, add_stat_1_3))
                         ) %>% 
          text_transform(fn = function(x) { str_replace(x, "(1.000)", "1") },
                         locations = cells_body(columns = c(add_stat_1_1, add_stat_1_2, add_stat_1_3))
                         ) %>% 
          tab_style(locations = cells_body(columns = label),
                      style = cell_text(weight = "bold")) %>% 
          tab_style(locations = cells_body(columns = add_stat_1_1,rows = add_stat_1_1 < 0.05),
                    style = cell_text(weight = "bold")) %>% 
          tab_style(locations = cells_body(columns = add_stat_1_2,rows = add_stat_1_2 < 0.05),
                    style = cell_text(weight = "bold")) %>% 
          tab_style(locations = cells_body(columns = add_stat_1_3,rows = add_stat_1_3 < 0.05),
                    style = cell_text(weight = "bold")) %>% 
          tab_style(location = cells_body(columns = c(add_stat_1_1, add_stat_1_2, add_stat_1_3)),
                    style = cell_fill(color = "grey70", alpha = 0.2)) %>% 
          tab_style(location = cells_body(columns = c(stat_1_1, stat_3_1, stat_5_1, stat_1_2, stat_3_2, stat_5_2, stat_1_3, stat_3_3, stat_5_3)),
                    style = cell_fill(color = "grey90", alpha = 0.2))
      )
```

## Conduct statistical test for the proportions differences in the clusters - between each two groups

``` r
set.seed(2023)

comp_pair <-
  PC_data[, sort(unique(Cluster)) %>% str_c("Cluster ",.)] %>%
  combn(., m = 2) %>% 
  as.data.table(.) %>% 
  as.list(.)

between_each_two_groups_comp <- 
  expand.grid(Group = PC_data[,levels(Group) %>% sort(.)],
              Category = c("CCS", "ICD"),
              Comp_pair = comp_pair,
              stringsAsFactors = FALSE) %>% 
  setDT(.) %>% 
  .[, Data := map(comorbidity_proportion_diff_test_data,
                  ~ .x %>% 
                    split(., by = "Group")) %>% 
              list_c(.) %>% 
              map(.,
                  ~ .x[, .SD, .SDcols = -c("Group")]) %>% 
              rep(., 
                  choose(PC_data[,uniqueN(Cluster)],2)) ] %>% 
  # Data manipulation for the statistical test
   .[, Test_data := 
        map2(Data,
             Comp_pair,
            ~ .x[Cluster %in% c(.y), ] %>% 
              melt.data.table(id.vars = "Cluster",
                              variable.name = "Variable",
                              value.name = "Affected") %>% 
              .[,Affected := factor(Affected, 
                                    levels = c("0","1"),
                                    labels = c("Unaffected", "Affected"))] %>% 
              split(., by = "Variable") %>% 
              map(.,
                  ~ .x[CJ(Cluster, Affected, unique = TRUE), 
                       .N, 
                       by = .EACHI,
                       on = .(Cluster, Affected)] %>% 
                    dcast.data.table(.,
                                     ... ~ Cluster,
                                     value.var = "N",
                                     fill = 0,
                                     drop = c(FALSE,FALSE)
                                     ) %>% 
                    .[,RowSum := rowSums(.SD), .SDcols = -c("Affected")] %>% 
                    melt.data.table(., 
                                id.vars = c("RowSum","Affected"), 
                                variable.name = "Cluster") %>% 
                    .[RowSum == 0,
                      value := fifelse(all(RowSum == 0) &
                                        Cluster == sample(unique(Cluster),1),
                                       1,
                                       value)] %>% 
                    dcast.data.table(., ... ~ Cluster) %>% 
                    .[order(Affected),] %>% 
                    .[, .SD, .SDcols = -c("RowSum","Affected")]) )] %>% 
  .[, Test_result := 
        map(Test_data,
            ~ .x %>% 
              map(.,
                  ~ fisher.test(.x,
                                simulate.p.value = TRUE,
                                B = 2000) %>% 
                    broom::tidy(.) ) %>% 
              rbindlist(.,idcol = "variable") %>% 
              { data.table( Variable = .$variable,
                            `p.value` = .$`p.value` %>% round(., 3),
                            `lower.95.CI` = .$`conf.low` %>% round(., 3),
                            `upper.95.CI` = .$`conf.high` %>% round(., 3) ) } %>% 
              .[, significance := fifelse(`p.value` < 0.05, "***", "")] )] %>% 
  .[, Comp_pair := 
        map_chr(Comp_pair,
                ~ str_c(.x, collapse = " vs. ") %>% 
                  str_replace(., "vs\\.\\sCluster", "vs\\.") ) ] %>% 
  .[order(Group, Comp_pair, Category),]
```

[^1]: Data availability in this repository is restricted due to the
    regulation of the law on the protection of patients’ data

[^2]: Please refer to the manuscript for **Figure 4**.
