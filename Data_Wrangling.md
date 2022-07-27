**Identifying heterogeneous subgroups of systemic connective tissue
diseases by applying a joint dimension reduction and clustering approach
to immunomarkers : a retrospective study** <br/>- data wrangling
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
                
                "ggplot2",
                "ggiraphExtra",
                "ggsci",
                "cowplot"
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
    ##        tidyr        purrr        furrr     devtools        knitr      ggplot2 
    ##         TRUE         TRUE         TRUE         TRUE         TRUE         TRUE 
    ## ggiraphExtra        ggsci      cowplot 
    ##         TRUE         TRUE         TRUE

## Install packages from git-hub

``` r
.ghub_pkgs <- c("DHLab-TSENG/dxpr")

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

    ## dxpr 
    ## TRUE

<br/>

------------------------------------------------------------------------

# Import data sets

## Raw data

## Reference data for exams

<br/>

------------------------------------------------------------------------

# Data pre-process

## Diagnosis data pre-process

``` r
ICD10_implt_date <- "2016-01-01"
diag_table_name  <- c("RSASOPDAF","RSASERDAF","RSASICDSUM")

icd_conversion_preps <- 
  raw_dataset[diag_table_name] %>% 
  map(., ~ list(dxDataFile  = .x,
                icdColName  = "DSSID",
                dateColName = "IPDAT"))

icd_conversion <- list()

for (i in diag_table_name) {
  
  dxDataFile <- copy(icd_conversion_preps[[i]][["dxDataFile"]])
  
  setnames(dxDataFile,
           c(icd_conversion_preps[[i]][["icdColName"]],icd_conversion_preps[[i]][["dateColName"]]),
           c("icdColName","dateColName"))
  
  icd_conversion_element <- 
    icdDxDecimalToShort(dxDataFile  = dxDataFile %>% .[!is.na(icdColName),],
                        icdColName  = icdColName,
                        dateColName = dateColName,
                        icd10usingDate = ICD10_implt_date) 
  
  icd_conversion <- append(icd_conversion,list(icd_conversion_element))
  
  rm(icd_conversion_element)
  
  if (i == last(diag_table_name)) { 
    names(icd_conversion) <- diag_table_name
    rm(icd_conversion_preps)
    gc(verbose = FALSE)
    }
  
} 
```

<br/>

### Summarise Top20 Error ICD

``` r
map(icd_conversion, 
    ~ plotICDError(errorFile = .x$Error,
                   icdVersion = all,
                   wrongICDType = all,
                   others = FALSE,
                   topN = 20))
```

    ## $RSASOPDAF
    ## $RSASOPDAF$graph

![](Data_Wrangling_files/figure-gfm/Summarise%20Error%20ICD-1.png)<!-- -->

    ## 
    ## $RSASOPDAF$ICD
    ##      ICD  count CumCountPerc IcdVersionInFile    WrongType Suggestion
    ##  1:  585 352697          29%            ICD 9 Wrong format       5859
    ##  2: 5640 164987       42.56%            ICD 9 Wrong format      56409
    ##  3: 2794 147258       54.67%            ICD 9 Wrong format      27949
    ##  4: 2740  62423        59.8%            ICD 9 Wrong format           
    ##  5: 3000  46289       63.61%            ICD 9 Wrong format      30009
    ##  6: 2554  42252       67.08%            ICD 9 Wrong format           
    ##  7: 5210  34720       69.94%            ICD 9 Wrong format      52109
    ##  8: 7805  33250       72.67%            ICD 9 Wrong format      78059
    ##  9: M791  32608       75.35%           ICD 10 Wrong format           
    ## 10: 7330  32187          78%            ICD 9 Wrong format      73309
    ## 11: V048  31911       80.62%            ICD 9 Wrong format      V0489
    ## 12: 2500  31783       83.24%            ICD 9 Wrong format           
    ## 13: 5714  31724       85.84%            ICD 9 Wrong format      57149
    ## 14: 5997  31312       88.42%            ICD 9 Wrong format           
    ## 15: 5234  26682       90.61%            ICD 9 Wrong format           
    ## 16: 6009  26107       92.76%            ICD 9 Wrong format           
    ## 17: 5280  23012       94.65%            ICD 9 Wrong format      52809
    ## 18: 5333  22564       96.51%            ICD 9 Wrong format           
    ## 19: 5339  21517       98.28%            ICD 9 Wrong format           
    ## 20: 7148  20979         100%            ICD 9 Wrong format      71489
    ## 
    ## 
    ## $RSASERDAF
    ## $RSASERDAF$graph

![](Data_Wrangling_files/figure-gfm/Summarise%20Error%20ICD-2.png)<!-- -->

    ## 
    ## $RSASERDAF$ICD
    ##      ICD count CumCountPerc IcdVersionInFile    WrongType Suggestion
    ##  1:  389  9711       25.68%            ICD 9 Wrong format       3899
    ##  2:  585  8581       48.37%            ICD 9 Wrong format       5859
    ##  3: 7806  7916        69.3%            ICD 9 Wrong format           
    ##  4: 5640  1704       73.81%            ICD 9 Wrong format      56409
    ##  5: 7890  1550        77.9%            ICD 9 Wrong format      78909
    ##  6:   90   951       80.42%            ICD 9 Wrong format           
    ##  7: 2765   927       82.87%            ICD 9 Wrong format           
    ##  8: M791   651       84.59%           ICD 10 Wrong format           
    ##  9: 7809   641       86.29%            ICD 9 Wrong format      78099
    ## 10: 7895   588       87.84%            ICD 9 Wrong format      78959
    ## 11: 5350   585       89.39%            ICD 9 Wrong format           
    ## 12: 5997   550       90.84%            ICD 9 Wrong format           
    ## 13: 8734   529       92.24%            ICD 9 Wrong format      87349
    ## 14: 2740   523       93.62%            ICD 9 Wrong format           
    ## 15: 7863   476       94.88%            ICD 9 Wrong format      78639
    ## 16: 2500   410       95.97%            ICD 9 Wrong format           
    ## 17:  539   402       97.03%            ICD 9 Wrong format           
    ## 18:  388   390       98.06%            ICD 9 Wrong format       3889
    ## 19:  465   381       99.07%            ICD 9 Wrong format       4659
    ## 20: 3000   352         100%            ICD 9 Wrong format      30009
    ## 
    ## 
    ## $RSASICDSUM
    ## $RSASICDSUM$graph

![](Data_Wrangling_files/figure-gfm/Summarise%20Error%20ICD-3.png)<!-- -->

    ## 
    ## $RSASICDSUM$ICD
    ##      ICD count CumCountPerc IcdVersionInFile    WrongType Suggestion
    ##  1: V581  3438       10.97%            ICD 9 Wrong format           
    ##  2:  389  3196       21.17%            ICD 9 Wrong format       3899
    ##  3:  585  2823       30.18%            ICD 9 Wrong format       5859
    ##  4: V451  2543        38.3%            ICD 9 Wrong format           
    ##  5: 2740  2082       44.95%            ICD 9 Wrong format           
    ##  6: 6000  1913       51.05%            ICD 9 Wrong format           
    ##  7: 5640  1796       56.78%            ICD 9 Wrong format      56409
    ##  8: 7054  1540        61.7%            ICD 9 Wrong format           
    ##  9: 2554  1495       66.47%            ICD 9 Wrong format           
    ## 10: 2880  1446       71.09%            ICD 9 Wrong format      28809
    ## 11:  414  1273       75.15%            ICD 9 Wrong format       4149
    ## 12: 7806  1233       79.08%            ICD 9 Wrong format           
    ## 13: 2848   985       82.23%            ICD 9 Wrong format      28489
    ## 14: 7895   960       85.29%            ICD 9 Wrong format      78959
    ## 15: 7032   921       88.23%            ICD 9 Wrong format           
    ## 16: 2765   852       90.95%            ICD 9 Wrong format           
    ## 17: 2874   751       93.35%            ICD 9 Wrong format      28749
    ## 18: 2824   741       95.71%            ICD 9 Wrong format      28249
    ## 19: 2794   710       97.98%            ICD 9 Wrong format      27949
    ## 20: 7070   633         100%            ICD 9 Wrong format      70709

<br/>

### Assign standardised (corrected) ICD Code to a new column and subset diagnosis related datasets

``` r
substitute_index <-
  map(diag_table_name, 
      ~ raw_dataset[[.x]] %>%
        .[,!is.na(.SD),.SDcols = "DSSID"]) %>%
  setNames(.,diag_table_name)

raw_dataset_list_dx <-
  map2(.x = diag_table_name,
       .y = substitute_index,
       ~ copy(raw_dataset[[.x]]) %>%
         .[.y[,1],] %>%
         cbind(.,icd_conversion[[.x]]$ICD) %>%
         setnames(.,"ICD","mod_DSSID")) %>%
  setNames(.,diag_table_name)
```

<br/>

## Define diseases and their corresponding ICD codes (grep with short form)

``` r
disease_ICD_table <-
  data.table(grepIcd = c("^7140+.*|^M05[4,5,7,8,9]+.*|^M06[0,2,3,8,9]+.*",
                         "^7100|^M32[1,8]+.*",
                         "^7101|^M34[0,1,9]+.*",
                         "^7102|^M350.*",
                         
                         "^7109|^M359",
                         "^710[3,4]+.*|^M33[1,2,9]+.*",
                         "^7108|^M351",
                         "^7200|^M459",
                         "^28981|^D6861"
                         ),
             Group   = c("Rheumatoid arthritis",
                         "Systemic lupus erythematosus",
                         "Systemic sclerosis",
                         "Sjogren's syndrome",
                         
                         "Undifferentiated connective tissue disease",
                         "Dermatopolymyositis",
                         "Mixed connective tissue disease",
                         "Ankylosing spondylitis",
                         "Antiphospholipid syndrome"
                         )
             )
```

<br/>

### Identify disease groups

``` r
dx_data       <- rbindlist(raw_dataset_list_dx,use.name = TRUE,fill = TRUE)
disease_group <- disease_ICD_table$Group

CustomGrepGroup <-
  map(disease_group,
      ~ icdDxToCustomGrep(dxDataFile  = dx_data,
                          idColName   = IDCODE,
                          icdColName  = DSSID,
                          dateColName = IPDAT,
                          customGroupingTable = disease_ICD_table[Group %in% .x,])) %>% 
  { list(map(., ~ .x$groupedDT %>% .[!is.na(GrepedGroup),]) %>% rbindlist,
         map(., ~ .x$summarised_groupedDT) %>% rbindlist) } %>% 
  setNames(.,c("groupedDT","summarised_groupedDT")) %>% 
  map(.,~ setnames(.x,"GrepedGroup","Group",skip_absent = TRUE))
```

<br/>

### Identify eligible cases for each disease with grouped ICD code

``` r
eligible_case_list <- 
  map(disease_group,
      ~ selectCases(dxDataFile = dx_data,
                    idColName = IDCODE,
                    icdColName = DSSID,
                    dateColName = IPDAT,
                    icd10usingDate = ICD10_implt_date,
                    groupDataType = customGrepIcdGroup,
                    customGroupingTable = disease_ICD_table,
                    caseCondition = .x,
                    caseCount = 2,
                    periodRange = c(30,365),
                    isDescription = FALSE,
                    caseName = .x) %>% 
      .[selectedCase %in% .x,]) %>% 
  setNames(.,disease_ICD_table$Group)
```

<br/>

### Identify eligible monoCTD cases and exclude UCTD cases

``` r
eligible_mono_CTD_dx_dataset <- 
  pluck(CustomGrepGroup,"groupedDT") %>% 
  .[ID %in% map_df(eligible_case_list,~ .x[,.(ID)])$ID, ] %>% 
  .[,Group := factor(Group,levels = disease_ICD_table$Group)] %>% 
  split(.,by = "ID") %>% 
  map(., ~ .x[.x[,uniqueN(Group)==1] ]) %>% 
  map(., ~ .x[.x[,str_detect(Group,"Undifferentiated connective tissue disease",negate = TRUE)] ]) %>% 
  keep(., ~ nrow(.x)>0) %>% 
  rbindlist(.)
```

<br/>

## Pre-process for lab data - merge with DATE-related data

``` r
# Retrieve DATE-related data from LAB-INDEX files
DAT_related_table <-
  copy(raw_dataset[["RSASLABINDX"]]) %>%
  .[,unique(.SD),.SDcols = c("IDCODE","LABNO","CLTDAT","BKGDAT","RCVDAT","ODRDPT")]

lab_exam_data <- 
  # Subset lab exam data by IDs for selected lab exam and eligible monoCTD cases
  copy(raw_dataset[["RSASLABRSLT"]]) %>% 
  .[LABIT %in% unique(exam_table$LABIT),] %>% 
  .[IDCODE %in% eligible_mono_CTD_dx_dataset[,unique(ID)],] %>% 
  
  # Merge 採檢日(CLTDAT), 收件日(RCVDAT), and 登記日(BKGDAT) from 檢驗索引檔(RSASLABINDX), mainly to query CLTDAT
  merge(.,DAT_related_table,by = c("IDCODE","LABNO"),all.x = TRUE) %>% 

  # Create a variable as an indicator of sample collection date
  # Priority: 採檢日(CLTDAT)* ➜ 收件日(RCVDAT)* ➜ 登記日(BKGDAT)* ➜ 輸入日(IPDAT)* ➜ 發報日(RDAT) ➜ 驗證日(VRFDAT)
  .[, SCDATE := CLTDAT] %>% 
  .[is.na(CLTDAT), SCDATE := RCVDAT, by = c("IDCODE","LABNO")] %>% 
  .[is.na(RCVDAT), SCDATE := BKGDAT, by = c("IDCODE","LABNO")] %>% 
  .[is.na(BKGDAT), SCDATE := IPDAT, by = c("IDCODE","LABNO")] %>% 
  
  # Annotate for each exam item with exam description
  .[LABIT=="72-112" & LABSH1IT=="72-111" & str_detect(LABNMABV,"^CRP"),LABSH1IT := "72-112"] %>% # unusual occurrence
  .[LABIT=="72-112" & LABSH1IT=="72-111" & str_detect(LABNMABV,"^HS"),LABIT := "72-111"] %>%     # unusual occurrence
  .[LABIT=="72-260" & is.na(LABSH1IT),LABSH1IT := "72-260"] %>%                                  # unusual occurrence
  # .[,LABNMABV := NULL] %>%
  merge(.,
        unique(exam_table[,.SD,.SDcols = c("LABIT","ITEM","LABSH1IT","ITEM_LABSH1IT")]),
        by = c("LABIT","LABSH1IT"),all.x = TRUE) %>% 
  
  # This step will remove 72-111
  .[LABIT %in% unique(exam_table$LABIT) ]
```

<br/>

## Pre-process for rheumatology and immunology reports data - merge with DATE-related data

``` r
rheuma_immuno_grepTable <- 
  exam_table[str_detect(LABIT,"^M"),] %>% 
  { data.table(exam_code = .[,as.character(LABIT)],
               exam_item = .[,ITEM],
               exam_ptrn = .[,Patterns]) }

rheuma_immuno_report_data_list <- 
  map(1:nrow(rheuma_immuno_grepTable), 
      ~ copy(raw_dataset[["RSASAIRRPT"]]) %>% 
        # Subset by exam_code for 風濕免疫科報告檔(RSASAIRRPT)
        .[str_detect(ITEM,rheuma_immuno_grepTable[.x,exam_code]),] %>% 
        .[,LABRESUVAL := 
            # Prune out headings
            sapply(str_split(REPORT01,"\\-{3,}"),'[',2) %>% 
            # Extract the expected sentence for numerical results
            str_sub(., start = str_locate(.,rheuma_immuno_grepTable[.x,exam_ptrn])[,1],end = -1) %>% 
            # Customized steps
            str_replace_all(.,c("Histone:" = "",
                                "Anti-Ribosomal P Ab:" = "",
                                "Rheumatoid Factor:" = "",
                                ">=" = ">",
                                "\\.{2}" = "\\.")) %>% 
            # Extract numerical results: cut off title
            { sapply(str_split(.,"\\:"),'[',2) } %>% 
            # Extract numerical results: remove units
            { sapply(str_split(.,"U\\/ml|ml|AU\\/ml|IU\\/mL|EliA|參考值|MPL units\\/|GPL units\\/|MPL\\-|GPL\\-|MPL\\/|GPL\\/"),'[',1) } %>% 
            # Second customized steps
            str_replace_all(.,c("^[\\s]+\\.|\\/[\\s]?$" = "")) ]
      ) %>% 
  setNames(.,rheuma_immuno_grepTable$exam_item)

rheuma_immuno_report_data <- 
  # Subset by IDs for eligible UCTD and monoCTD cases
  map(rheuma_immuno_report_data_list,
      ~ copy(.x) %>% 
        .[IDCODE %in% eligible_mono_CTD_dx_dataset[,unique(ID)],]) %>% 
  
  # Subset for certain columns (ODATE: 開單日期, EDATE: 收件日期, RDATE: 報告日期)
  map(., 
      ~ .x[,.SD,.SDcols = c("IDCODE","LOC","ODATE","EDATE","RDATE","SOURCE","LABRESUVAL")]) %>% 
  
  # Combine all exam results into one dataset
  rbindlist(.,idcol = "ITEM") %>% 
  
  # Assign exam code according to ITEM
  merge(.,rheuma_immuno_grepTable[,.SD,.SDcols = c("exam_item","exam_code")],by.x = "ITEM",by.y = "exam_item") %>% 
  
  # Create a variable as an exam item indicator
  .[,LABSH1IT := str_replace(exam_code,"\\.|\\*","")] %>% 
  .[,ITEM_LABSH1IT := paste(ITEM,LABSH1IT,sep = "_")] %>% 
  
  # Create a variable as an indication of sample collection date
  # Priority: 開單日(ODATE) ➜ 收件日期(EDATE) ➜ 報告日期(RDATE)
  .[,SCDATE := ODATE] %>% 
  .[is.na(ODATE),SCDATE := EDATE] %>% 
  .[is.na(EDATE),SCDATE := RDATE] %>% 
  
  # Remove interim variables
  .[,.SD,.SDcols = -c("exam_code")] %>% 
  
  # Patterns reconition for classifying exam results
  .[str_detect(LABRESUVAL,"\\-") & str_detect(LABRESUVAL,"\\+"),LABRESUVAL := "Equivocal"] %>% 
  .[str_detect(LABRESUVAL,"\\-"),                               LABRESUVAL := "Normal"] %>% 
  .[str_detect(LABRESUVAL,"\\+|.*[Aa]bove.*|.*[Pp]ositive.*" ), LABRESUVAL := "Abnormal"] %>% 
  .[str_detect(LABRESUVAL,"([0-9*]\\.?[0-9*])",negate = TRUE) & 
    str_detect(LABRESUVAL,"Normal|Abnormal|Equivocal",negate = TRUE),LABRESUVAL := NA_character_] %>% 
  .[str_detect(LABRESUVAL,"[:alpha:]+") & 
    str_detect(LABRESUVAL,"Normal|Abnormal|Equivocal",negate = TRUE),LABRESUVAL := str_replace_all(LABRESUVAL,"[:alpha:]+","")] %>% 
  .[str_detect(LABRESUVAL,"\\\\|\\/"),                               LABRESUVAL := str_replace_all(LABRESUVAL,"\\\\|\\/","")] %>% 
  
  # (optional) Rearrange dataset
  .[order(ITEM_LABSH1IT,IDCODE,SCDATE),]
```

<br/>

## Combine lab and rheumatology and immunology report data

``` r
var_subset <- c("IDCODE","ITEM_LABSH1IT","SCDATE","LABRESUVAL")

lab_report_results <- 
  map(list(lab_exam_data,rheuma_immuno_report_data),
      ~ copy(.x) %>% 
        .[,.SD,.SDcols = var_subset] %>% 
        .[,LABRESUVAL := str_replace(LABRESUVAL,"NA",NA_character_)]) %>% 
  rbindlist(.) %>% 
  .[order(IDCODE,SCDATE,ITEM_LABSH1IT)]
```

<br/>

## Determining status of lab/report reulsts using specified cut-off points

### Create a result data container for the combined lab and rheumatology and immunology report data

``` r
complete_sex_bday_table <-
  pluck(raw_dataset,"RSASHOSP") %>% 
  copy(.) %>% 
  .[IDCODE %in% eligible_mono_CTD_dx_dataset[,unique(ID)],] %>% 
  .[,unique(.SD),.SDcols = c("IDCODE","SEX","BIRTHDAY")] %>% 
  .[,BIRTHDAY := ymd(BIRTHDAY)]

exam_result_dataset <- 
  # Step 1: Retrieve SEX and BIRTHDAY from 歸戶檔(RSAHOSP), then calculating AGE
  merge(copy(lab_report_results),
        copy(complete_sex_bday_table),
        by = "IDCODE",all.x = TRUE) %>% 
  .[,AGE := BIRTHDAY %--% SCDATE %>% as.period(.,unit = "years") %>% year(.)] %>% 

  # Step 2: Retrieve LABVALUN from lab_exam_data
  merge(.,
        { copy(lab_exam_data) %>% 
          .[,unique(.SD),.SDcols = c("IDCODE","SCDATE","ITEM_LABSH1IT","LABVALUN","IPDAT")] %>%
          .[order(IDCODE,ITEM_LABSH1IT,LABVALUN),head(.SD,1L),by = c("IDCODE","ITEM_LABSH1IT","IPDAT")] %>% 
          .[,.SD,.SDcols = c("IDCODE","SCDATE","ITEM_LABSH1IT","LABVALUN")] } ,
        by = c("IDCODE","SCDATE","ITEM_LABSH1IT"),all.x = TRUE) %>% 
  
  # Step 3: Retrieve ITEM and LABSH1IT by splitting strings of ITEM_LABSH1IT
  .[,':=' (ITEM     = str_split(ITEM_LABSH1IT,"_",n = 2,simplify = TRUE)[,1],
           LABSH1IT = str_split(ITEM_LABSH1IT,"_",n = 2,simplify = TRUE)[,2])] %>% 
  
  # Step 4: Transform SCDATE from strings to date format
  .[,SCDATE := ymd(SCDATE)]
```

<br/>

### Clean raw results of lab/report data and merge processed data with reference information

``` r
exam_result_dataset_ori <- exam_result_dataset
exam_item <- unique(exam_result_dataset$ITEM_LABSH1IT)
indicator_vars <- c("Ref_LABVALUN","Ref_SEX","Ref_AGE","Ref_Duration")
classify_vars  <- c("Ref_Neg_L_Bound","Ref_Neg_U_Bound","Ref_Equ_L_Bound","Ref_Equ_U_Bound",
                    "Ref_Pos_L_Bound","Ref_Pos_U_Bound","non_Dich_Category")

normal_pattern    <- "Normal|[Nn]egative|NEGATIVE|Nonreactive|\\(\\-\\)|^\\-$"
abnormal_pattern  <- "[Aa]bnormal|[Pp]ositive|POSITIVE|[Aa]ntibody|Ab$|[Aa][Bb]\\s*\\(\\+\\)|\\(\\+\\)|\\(\\+|\\(1\\+\\)|1\\+|\\+"
equivocal_pattern <- "[Ee]quivocal|EQUIVOCAL|[Bb]orderline|BORDERLINE|[Ww]eakly|WEAKLY|WEAK|螢光"

UNIT_colname <- "LABVALUN"
SEX_colname  <- "SEX"
AGE_colname  <- "AGE"
DATE_colname <- "SCDATE"

for (i in exam_item) {
  
  ### Subset lab exam data by exam entries (ITEM_LABSH1IT) and preprocess exam values to extract pure values
  exam_result_data_subset <- 
        copy(exam_result_dataset) %>% 
        
        # Subset dataset by an experiment item name and code
        .[ITEM_LABSH1IT %in% i,] %>% 
        
        .[,mdf_LABRESUVAL := LABRESUVAL] %>% 
        # remove all space in imputed values
        .[str_detect(mdf_LABRESUVAL,"\\s"),mdf_LABRESUVAL := str_replace_all(mdf_LABRESUVAL,"\\s","")] %>% 
        # replace strings containing positive and negative patterns with POSITIVE and NEGATIVE
        .[str_detect(mdf_LABRESUVAL,"[Pp]ositive|POSITIVE"),mdf_LABRESUVAL := "Positive"] %>% 
        .[str_detect(mdf_LABRESUVAL,"[Nn]egative|NEGATIVE"),mdf_LABRESUVAL := "Negative"] %>% 
  
        # remove all concatenated 'NA' strings
        .[str_detect(mdf_LABRESUVAL,"NA\\,|\\,NA") &
          str_detect(mdf_LABRESUVAL,paste0(normal_pattern,abnormal_pattern,equivocal_pattern,collapse = "|"),negate = TRUE),
          mdf_LABRESUVAL := str_replace_all(mdf_LABRESUVAL,"NA\\,|\\,NA","")] %>% 
        
        # assign Abnormality results
        .[str_detect(mdf_LABRESUVAL,normal_pattern),   final_LABRESUVAL := "Normal"] %>% 
        .[str_detect(mdf_LABRESUVAL,abnormal_pattern), final_LABRESUVAL := "Abnormal"] %>% 
        .[str_detect(mdf_LABRESUVAL,equivocal_pattern),final_LABRESUVAL := "Equivocal"] %>% 
        
        # pattern 1
        .[str_detect(mdf_LABRESUVAL,"[A-Za-z]") & str_detect(mdf_LABRESUVAL,"[0-9]",negate = TRUE) &
          str_detect(mdf_LABRESUVAL,paste0(normal_pattern,abnormal_pattern,equivocal_pattern,collapse = "|"),negate = TRUE),
          final_LABRESUVAL := "NA"] %>%
        
        # pattern 2
        .[str_detect(mdf_LABRESUVAL,"[A-Za-z]") & str_detect(mdf_LABRESUVAL,"[0-9]+|[0-9]+\\.[0-9]+") &
          str_detect(mdf_LABRESUVAL,"[0-9]+\\.[0-9]+\\E\\-[1-9]",negate = TRUE) &
          str_detect(mdf_LABRESUVAL,paste0(normal_pattern,abnormal_pattern,equivocal_pattern,collapse = "|"),negate = TRUE),
          final_LABRESUVAL := str_replace_all(mdf_LABRESUVAL,"[A-Za-z()]","")] %>%
    
        # pattern 3
        .[str_detect(mdf_LABRESUVAL,"[0-9]\\.$"),final_LABRESUVAL := "0" ] %>% 
        .[str_detect(mdf_LABRESUVAL,"\\.{2,}"),final_LABRESUVAL := str_replace_all(mdf_LABRESUVAL,"\\.+","\\.") ] %>% 
        .[str_detect(mdf_LABRESUVAL,"[0-9]\\.[0-9]+\\E\\-[1-9]"),final_LABRESUVAL := as.numeric(mdf_LABRESUVAL) %>% as.character(.)] %>% 
    
        # pattern 4
        .[str_detect(mdf_LABRESUVAL,"^\\.$|^\\-$|\\[\\ *空白\\ *\\]|^\\]$"),final_LABRESUVAL := "NA"] %>% 
        
        # pattern 5
        .[str_detect(mdf_LABRESUVAL,"[>≥]") & 
          str_detect(mdf_LABRESUVAL,"\\:",negate = TRUE) &
          str_detect(mdf_LABRESUVAL,paste0(normal_pattern,abnormal_pattern,equivocal_pattern,collapse = "|"),negate = TRUE),
          final_LABRESUVAL := str_replace_all(mdf_LABRESUVAL,"[>≥A-Za-z()]","") %>% { as.numeric(.)*2 } %>% as.character(.) %>% 
                              str_replace_all(.,"[<≤A-Za-z()]","") %>% { as.numeric(.)/2 } %>% as.character(.)] %>%
        .[str_detect(mdf_LABRESUVAL,"[<≤]") & 
          str_detect(mdf_LABRESUVAL,"\\:",negate = TRUE) &
          str_detect(mdf_LABRESUVAL,paste0(normal_pattern,abnormal_pattern,equivocal_pattern,collapse = "|"),negate = TRUE),
          final_LABRESUVAL := str_replace_all(mdf_LABRESUVAL,"[<≤A-Za-z()]","") %>% { as.numeric(.)/2 } %>% as.character(.) %>% 
                              str_replace_all(.,"[>≥A-Za-z()]","") %>% { as.numeric(.)*2 } %>% as.character(.)] %>%
        
        # pattern 6
         .[str_detect(mdf_LABRESUVAL,"[>≥]") &
           str_detect(mdf_LABRESUVAL,"\\,") &
           str_detect(mdf_LABRESUVAL,paste0(normal_pattern,abnormal_pattern,equivocal_pattern,collapse = "|"),negate = TRUE),
           final_LABRESUVAL := str_split(mdf_LABRESUVAL,",") %>%
                               map(., ~ ifelse(str_detect(.x,"[>≥A-Za-z()]"),{ str_replace_all(.x,"[>≥A-Za-z()]","") %>% as.numeric(.)*2 },.x)) %>%
                               map(., ~ ifelse(str_detect(.x,"[<≤A-Za-z()]"),{ str_replace_all(.x,"[<≤A-Za-z()]","") %>% as.numeric(.)/2 },.x)) %>% 
                               map_chr(.,~ str_c(.x,collapse = ","))] %>%
         .[str_detect(mdf_LABRESUVAL,"[<≤]") &
           str_detect(mdf_LABRESUVAL,"\\,") &
           str_detect(mdf_LABRESUVAL,paste0(normal_pattern,abnormal_pattern,equivocal_pattern,collapse = "|"),negate = TRUE),
           final_LABRESUVAL := str_split(mdf_LABRESUVAL,",") %>%
                               map(., ~ ifelse(str_detect(.x,"[<≤A-Za-z()]"),{ str_replace_all(.x,"[<≤A-Za-z()]","") %>% as.numeric(.)/2 },.x)) %>%
                               map(., ~ ifelse(str_detect(.x,"[>≥A-Za-z()]"),{ str_replace_all(.x,"[>≥A-Za-z()]","") %>% as.numeric(.)*2 },.x)) %>%
                               map_chr(.,~ str_c(.x,collapse = ","))] %>% 
        
        # ANA only
        .[,ANA_mdf_LABRESUVAL := "NA",] %>% 
        { if(str_detect(i,"ANA")) {
             .[,ANA_mdf_LABRESUVAL := { map(mdf_LABRESUVAL, ~ str_split(.x,","))  %>%
                                        map(., ~ map(.x, ~ str_extract(.x,":[0-9]*"))) %>%
                                        map(., ~ map(.x, ~ str_replace(.x,":",""))) %>%
                                        map(., ~ map(.x, ~ unique(.x))) %>%
                                        map(., ~ map(.x, ~ str_replace_na(.x))) %>% 
                                        map(., ~ map_chr(.x, ~ str_c(.x,collapse = ","))) %>%
                                        map(., ~ map(.x, ~ str_replace(.x,"NA\\,|\\,NA",""))) %>%
                                        unlist(.) }] %>% 
            .[!is.na(LABRESUVAL) & !LABRESUVAL %in% c("-","[ 空白 ]"),.SD]
          } else { .[,.SD] }
        }
    

  ### Subset exam table with specified exam item
  exam_table_subset <- exam_table[ITEM_LABSH1IT %in% i,] %>% 
                       { .[,m_index := seq(1,nrow(.)) %>% as.character(.)] }
  
  #### Conditional on differences in unit of exam results (LABVALUN)
  if (exam_table[ITEM_LABSH1IT %in% i,
                 map(.SD, ~ str_detect(.x,pattern = "-",negate = TRUE) %>% sum),
                 .SDcols = indicator_vars[1]] %>% 
      rowSums(.) > 0 ) {

  subset_condition <- paste0(UNIT_colname,"==","'",exam_table_subset[,Ref_LABVALUN],"'")

  #### Conditional on differences in SEX and AGE
  } else if (exam_table[ITEM_LABSH1IT %in% i,
                        map(.SD, ~ str_detect(.x,pattern = "-",negate = TRUE) %>% sum),
                        .SDcols = indicator_vars[2:3]] %>% 
             rowSums(.) > 0 ) {

  subset_condition <- paste(paste0(SEX_colname,"==","'",exam_table_subset[,Ref_SEX],"'"),
                            paste0(AGE_colname,exam_table_subset[,Ref_AGE]),
                            sep = " & ")
  
  ####  Conditional on differences in the year that the exam implemented (SCDATE)
  } else if (exam_table[ITEM_LABSH1IT %in% i,
                        map(.SD, ~ str_detect(.x,pattern = "-",negate = TRUE) %>% sum),
                        .SDcols = indicator_vars[4]] %>% 
             rowSums(.) > 0 ) {
  
  subset_condition <- ifelse(exam_table_subset[,str_detect(Ref_Duration,"-")],
                             paste0("between(",
                                    "year(",DATE_colname,")",
                                    ",",
                                    exam_table_subset[2,str_split(Ref_Duration,"-")][1,],
                                    ",",
                                    exam_table_subset[2,str_split(Ref_Duration,"-")][2,],
                                    ")"),
                             paste0("year(",DATE_colname,")",exam_table_subset[,Ref_Duration]))
  
  #### No condition
  } else {
    
    subset_condition <- NULL
    
  }
  
  ### Generate merging index between subset of lab exam data and exam table with subset conditions
  if(!is.null(subset_condition)) {
     merge_index <- exam_result_data_subset[,map_dfc(subset_condition, ~ { eval(parse(text=.x)) })] %>%
                    setNames(.,subset_condition) %>% setDT(.) %>% 
                    .[,str_subset(subset_condition,"N/A") := map(.SD, ~ replace_na(.x,TRUE)),.SDcols = str_subset(subset_condition,"N/A")] %>% 
                    pmap_dfr(., ~ data.table(index = which(c(...) == TRUE))) 
  } else {
    merge_index <- data.table(index = rep("1",nrow(exam_result_data_subset)))
    }
     
  #### Add an index column back to lab exam data subset
  if(all(str_detect(subset_condition,"SEX|AGE"),!is.null(subset_condition))) {
    exam_result_data_subset[!is.na(SEX) & !is.na(AGE),m_index := merge_index$index %>% as.character(.)]
  } else {
    exam_result_data_subset[,m_index := merge_index$index %>% as.character(.)] 
    }
  
  ### Merge reference values of exam results with index composition (partial)
  exam_result_data_subset[exam_table_subset,on = c("ITEM_LABSH1IT","m_index"),
                          ':=' (Ref_LABVALUN      = i.Ref_LABVALUN,
                                Ref_SEX           = i.Ref_SEX,
                                Ref_AGE           = i.Ref_AGE,
                                Ref_Duration      = i.Ref_Duration,
                                Ref_Neg_L_Bound   = i.Ref_Neg_L_Bound,
                                Ref_Neg_U_Bound   = i.Ref_Neg_U_Bound,
                                Ref_Equ_L_Bound   = i.Ref_Equ_L_Bound,
                                Ref_Equ_U_Bound   = i.Ref_Equ_U_Bound,
                                Ref_Pos_L_Bound   = i.Ref_Pos_L_Bound,
                                Ref_Pos_U_Bound   = i.Ref_Pos_U_Bound,
                                non_Dich_Category = i.non_Dich_Category)] %>% 
  .[is.na(final_LABRESUVAL),final_LABRESUVAL  := mdf_LABRESUVAL] %>%
  .[ANA_mdf_LABRESUVAL != "NA",final_LABRESUVAL := ANA_mdf_LABRESUVAL] %>% 
  .[,.SD,.SDcols = c("ITEM_LABSH1IT","IDCODE",DATE_colname,"final_LABRESUVAL","ANA_mdf_LABRESUVAL",classify_vars)] %>% 
  
  ### Merge reference values and exam results back to original dataset for all exam entries
  exam_result_dataset[.,on = c("ITEM_LABSH1IT","IDCODE",DATE_colname),
                      ':=' (final_LABRESUVAL   = i.final_LABRESUVAL,
                            ANA_mdf_LABRESUVAL = i.ANA_mdf_LABRESUVAL,
                            Ref_Neg_L_Bound    = i.Ref_Neg_L_Bound,
                            Ref_Neg_U_Bound    = i.Ref_Neg_U_Bound,
                            Ref_Equ_L_Bound    = i.Ref_Equ_L_Bound,
                            Ref_Equ_U_Bound    = i.Ref_Equ_U_Bound,
                            Ref_Pos_L_Bound    = i.Ref_Pos_L_Bound,
                            Ref_Pos_U_Bound    = i.Ref_Pos_U_Bound,
                            non_Dich_Category  = i.non_Dich_Category)] %>% 
  
  ### Transform "NA" into NA
  { set(.,which(.[["ANA_mdf_LABRESUVAL"]]=="NA"),"ANA_mdf_LABRESUVAL",NA) } %>% 
  { set(.,which(.[["final_LABRESUVAL"]]=="NA"),"final_LABRESUVAL",NA) }
  
}
```

<br/>

### Determine status of exam results using specified cut-off points

``` r
exam_result_status_dataset <- 
  copy(exam_result_dataset) %>% 
  
  # Assign status of results for those have determined
  .[final_LABRESUVAL %in% c("Normal","Abnormal","Equivocal"),Status := final_LABRESUVAL] %>% 
    
  # Create a new variable to determine status for those with multiple exam results 
  .[!is.na(final_LABRESUVAL) & is.na(Status),
    "multi_final_LABRESUVAL_list" := map(.SD, ~ as.list(str_split(.x,"\\,"))),.SDcols = "final_LABRESUVAL"] %>% 
  
  # Status classify ---------------------------------------------------------------------------------------------------
  ## with only negative upper bound
  .[!is.na(final_LABRESUVAL) & is.na(Status) &
    Ref_Neg_L_Bound == "-" & Ref_Neg_U_Bound != "-" & Ref_Pos_U_Bound == "-" & Ref_Pos_L_Bound == "-",
    inequation_LABRESUVAL := map(.SD, ~ map2(.x,Ref_Neg_U_Bound, ~ paste0(.x,.y))),.SDcols = "multi_final_LABRESUVAL_list"] %>%
  .[!is.na(final_LABRESUVAL) & is.na(Status) &
    Ref_Neg_L_Bound == "-" & Ref_Neg_U_Bound != "-" & Ref_Pos_U_Bound == "-" & Ref_Pos_L_Bound == "-",
    evaluation_LABRESUVAL := map(inequation_LABRESUVAL, ~ map_lgl(.x, ~ eval(parse(text = .x)))) %>% 
                                                                        imap(., ~ all(.x)) %>% 
                                                                        flatten_lgl] %>%
  .[!is.na(final_LABRESUVAL) & is.na(Status) &
    Ref_Neg_L_Bound == "-" & Ref_Neg_U_Bound != "-" & Ref_Pos_U_Bound == "-" & Ref_Pos_L_Bound == "-",
    Status := fifelse(evaluation_LABRESUVAL,"Normal","Abnormal")] %>% 

  ## with negative lower and upper bound
  .[!is.na(final_LABRESUVAL) & is.na(Status) &
    Ref_Neg_L_Bound != "-" & Ref_Neg_U_Bound != "-" & Ref_Pos_U_Bound == "-" & Ref_Pos_L_Bound == "-",
    inequation_LABRESUVAL := map(.SD, ~ map2(.x,Ref_Neg_L_Bound, ~ paste0(.x,.y,"&",.x)) %>%
                                        map2(.,Ref_Neg_U_Bound, ~ paste0(.x,.y))),.SDcols = "multi_final_LABRESUVAL_list"] %>%
  .[!is.na(final_LABRESUVAL) & is.na(Status) &
    Ref_Neg_L_Bound != "-" & Ref_Neg_U_Bound != "-" & Ref_Pos_U_Bound == "-" & Ref_Pos_L_Bound == "-",
    evaluation_LABRESUVAL := map(inequation_LABRESUVAL, ~ map_lgl(.x, ~ eval(parse(text = .x)))) %>% 
                                                                        imap(., ~ all(.x)) %>% 
                                                                        flatten_lgl] %>%
  .[!is.na(final_LABRESUVAL) & is.na(Status) &
    Ref_Neg_L_Bound != "-" & Ref_Neg_U_Bound != "-" & Ref_Pos_U_Bound == "-" & Ref_Pos_L_Bound == "-",
    Status := fifelse(evaluation_LABRESUVAL,"Normal","Abnormal")] %>% 
  
  # with non-dichotomous status & for abnormal
  .[non_Dich_Category == "Y" & !is.na(final_LABRESUVAL) & is.na(Status),
    inequation_LABRESUVAL := map(.SD, ~ map2(.x,Ref_Pos_L_Bound, ~ paste0(.x,.y))),.SDcols = "multi_final_LABRESUVAL_list"] %>%
  .[non_Dich_Category == "Y" & !is.na(final_LABRESUVAL) & is.na(Status),
    evaluation_LABRESUVAL := map(inequation_LABRESUVAL, ~ map_lgl(.x, ~ eval(parse(text = .x)))) %>%
                                                                        imap(., ~ any(.x)) %>%
                                                                        flatten_lgl] %>%
  .[non_Dich_Category == "Y" & !is.na(final_LABRESUVAL) & is.na(Status),
    Status := fifelse(evaluation_LABRESUVAL,"Abnormal",NA_character_)] %>% 
  
  # with non-dichotomous status & for equivocal
  .[non_Dich_Category == "Y" & !is.na(final_LABRESUVAL) & is.na(Status),
    inequation_LABRESUVAL := map(.SD, ~ map2(.x,Ref_Equ_L_Bound, ~ paste0(.x,.y,"&",.x)) %>%
                                        map2(. ,Ref_Equ_U_Bound, ~ paste0(.x,.y))),.SDcols = "multi_final_LABRESUVAL_list"] %>%
  .[non_Dich_Category == "Y" & !is.na(final_LABRESUVAL) & is.na(Status),
    evaluation_LABRESUVAL := map(inequation_LABRESUVAL, ~ map_lgl(.x, ~ eval(parse(text = .x)))) %>% 
                                                                        imap(., ~ any(.x)) %>% 
                                                                        flatten_lgl] %>%
  .[non_Dich_Category == "Y" & !is.na(final_LABRESUVAL) & is.na(Status),
    Status := fifelse(evaluation_LABRESUVAL,"Equivocal",NA_character_)] %>% 
  
  # with non-dichotomous status & for normal
  .[non_Dich_Category == "Y" & !is.na(final_LABRESUVAL) & is.na(Status),
    inequation_LABRESUVAL := map(.SD, ~ map2(.x,Ref_Neg_U_Bound, ~ paste0(.x,.y))),.SDcols = "multi_final_LABRESUVAL_list"] %>%
  .[non_Dich_Category == "Y" & !is.na(final_LABRESUVAL) & is.na(Status),
    evaluation_LABRESUVAL := map(inequation_LABRESUVAL, ~ map_lgl(.x, ~ eval(parse(text = .x)))) %>% 
                                                                        imap(., ~ all(.x)) %>% 
                                                                        flatten_lgl] %>%
  .[non_Dich_Category == "Y" & !is.na(final_LABRESUVAL) & is.na(Status),
    Status := fifelse(evaluation_LABRESUVAL,"Normal",NA_character_)]
```

<br/>

### Determine simultaneously the status for exam items with multiple test results

``` r
simul_classify_item <- 
  list(c("72-134","M25-238"),                  # Rheumatoid factor (RF)
       c("72-245","72A245","72B245","72C245"), # Antinuclear antibody (ANA)
       c("72-247","M25-150")                   # anti-Double strand DNA antibody (anti-dsDNA)
       )

walk(simul_classify_item, 
     ~ exam_result_status_dataset %>% 
       .[LABSH1IT %in% .x,
         Status := fifelse(all(is.na(Status)),NA_character_,
                           fifelse(any(Status=="Abnormal",na.rm = TRUE),"Abnormal","Normal")),
         by = .(IDCODE,SCDATE)])

# Remove all duplicated exam item, keep the main item only
dup_exam_item <- 
  map(simul_classify_item, ~ .x[-1],) %>% 
  flatten_chr(.)

dup_exam_removed_result_status_dataset <- 
  copy(exam_result_status_dataset) %>% 
  .[!LABSH1IT %in% dup_exam_item,]
```

<br/>

------------------------------------------------------------------------

## Combine within-seven-day exam results

### Create merging index

``` r
merge_within_boundary_data <- 
  # Calculate rolling difference in days and assign an column to indicate whether the row has been grouped
  function(data = data,
           Group_ID = Group_ID,
           boundary_less_than_equal_to = boundary_less_than_equal_to) {
    
    combine_interim <- 
        copy(data) %>% 
        .[,TIMEDIFF := min(SCDATE) %--% SCDATE %>% 
                       as.period(.,unit = "days") %>% 
                       as.numeric(.,unit = "days"),
          by = c(Group_ID)] %>% 
        .[TIMEDIFF >= as.numeric(boundary_less_than_equal_to), TIMEDIFF := -1]
    
    # Collect rows have been grouped and those have not
    combine_index   <- combine_interim[TIMEDIFF!=-1,]
    combine_interim <- combine_interim[TIMEDIFF==-1,]
    
    # Iterate rolling difference calculation until all rows are grouped
    while( any(combine_interim$TIMEDIFF==-1) ) {
      
      combine_interim %>% 
        .[,TIMEDIFF := min(SCDATE) %--% SCDATE %>% 
                       as.period(.,unit = "days") %>% 
                       as.numeric(.,unit = "days"),
          by = c(Group_ID)] %>% 
        .[TIMEDIFF >= as.numeric(boundary_less_than_equal_to), TIMEDIFF := -1]
      
      combine_index   <- rbind(combine_index,combine_interim[TIMEDIFF!=-1,])
      combine_interim <- combine_interim[TIMEDIFF==-1,]
      
      }
     
    # Generate group index
    combine_index <- 
      copy(combine_index) %>% 
        .[,INDEX := cumsum(ifelse(shift(TIMEDIFF,n = 1,type = "lag",fill = 0)-TIMEDIFF!=1 & 
                                   shift(IDCODE,n = 1,type = "lag",fill = 0)==IDCODE | 
                                  is.na(shift(TIMEDIFF,n = 1,type = "lag",fill = 0)-TIMEDIFF) & 
                                   shift(IDCODE,n = 1,type = "lag",fill = 0)==IDCODE,
                                  0,1))] %>% 
        .[,.SD,.SDcols = c(Group_ID,"SCDATE","TIMEDIFF","INDEX")]
    
    return(combine_index)
  }

merge_index <- 
  merge_within_boundary_data(data = dup_exam_removed_result_status_dataset,
                             Group_ID = "IDCODE", 
                             boundary_less_than_equal_to = 7) %>% 
  .[order(IDCODE,SCDATE),]
```

<br/>

### Merge index with the original exam data to prepare final lab and report data

``` r
exam_result_status_data <- 
  # Append grouping index to exam and report result data
  cbind(copy(dup_exam_removed_result_status_dataset),
        copy(merge_index) %>% .[,.SD,.SDcols = c("TIMEDIFF","INDEX")]) %>% 
  # Replace SCDATE with the earliest SCDATE within each group of INDEX 
  .[,mod_SCDATE := min(SCDATE),by = "INDEX"] %>% 
  # Subset by variables
  .[,.SD,.SDcols = c("IDCODE","SCDATE","mod_SCDATE","ITEM_LABSH1IT",
                     "Status","SEX","AGE")] %>% 
  # Transform data set into wide format for completing ITEM_LABSH1IT
  dcast(., ... ~ ITEM_LABSH1IT,fun.aggregate = toString,value.var = "Status",
        fill = NA_character_,drop = c(TRUE,FALSE)) %>%
  # Transform data set back to long format
  melt(.,id.vars = c("IDCODE","SCDATE","mod_SCDATE","SEX","AGE"),
       variable.name = "ITEM_LABSH1IT",value.name = "Status") %>%
  
  .[,Status := str_replace_all(Status,"NA,\\s*|,\\s*NA","")] %>%
  .[,Status := str_replace_all(Status,"NA",NA_character_)] %>%
  #  Assessing Status final results
  .[!is.na(Status),overall_Status := fifelse(str_detect(Status,"Abnormal"),"Abnormal","Normal")] %>% 
  
  # added step for removing merge-within-seven-day step
  .[,SCDATE := mod_SCDATE] %>% 
  .[,-c("mod_SCDATE")]
```

<br/>

------------------------------------------------------------------------

# Merge lab results with diagnosis data

## Calculate time intervals between the first exam and the first diagnosis using raw data

``` r
# Retrieve first diagnosis date 
monoCTD_diagnosis_query_table <- 
  copy(eligible_mono_CTD_dx_dataset) %>% 
  .[,unique(.SD)] %>% 
  .[order(ID,Date,Group),head(.SD,1L),by = c("ID","Group")] %>% 
  setnames(.,"Date","First_diagnosis_date")

# Merge exam result status data with monoCTD diagnosis data
monoCTD_dup_exam_removed_result_status_dataset <- 
  copy(exam_result_status_data) %>%
  merge(.,
        monoCTD_diagnosis_query_table,
        by.x = "IDCODE",by.y = "ID",
        all.x = TRUE) %>% 
  ## remove exam results belong to ever-UCTD cases 
  .[!is.na(Group)]

# --------------------------------------------------------------------------------------------------------------------- #
# Calculate time intervals between the first exam and the first diagnosis using raw data
monoCTD_first_diag_exam_span_raw_data <- 
  copy(monoCTD_dup_exam_removed_result_status_dataset) %>% 
  .[,first_diag_exam_span := First_diagnosis_date %--% SCDATE %>% 
                             as.period(.,unit = "days") %>% 
                             day(.),] %>% 
  .[,first_diag_exam_span_abs := abs(first_diag_exam_span),] %>% 
  .[,over_an_year := fifelse(abs(first_diag_exam_span)>=180,"Yes","No")] %>% 
  # extract results with the minimal gap for each case and exam item
  .[order(IDCODE,ITEM_LABSH1IT,first_diag_exam_span_abs),head(.SD,1L),by = c("IDCODE","ITEM_LABSH1IT")]
```

<br/>

## Calculate interval between first diagnosis and each exam records, and subset the data set by the specified intervals

``` r
monoCTD_first_diag_exam_span_calculation_data <-
 copy(monoCTD_dup_exam_removed_result_status_dataset) %>%
  # calculate intervals between first diagnosis date and sample collection date
  .[,first_diag_exam_span := First_diagnosis_date %--% SCDATE %>% 
                             as.period(.,unit = "days") %>% 
                             day(.),] %>% 
  .[,first_diag_exam_dist := first_diag_exam_span^2 + 
                             abs(first_diag_exam_span) + 
                             first_diag_exam_span] %>% 
  # order with prioritising those with testing results than those with shorter interval
  .[order(IDCODE,ITEM_LABSH1IT,overall_Status,first_diag_exam_dist),.SD,by = c("IDCODE","ITEM_LABSH1IT")]

interval_lower_limit        <- -30
expand_interval_lower_limit <- -90 # Sensitivity analysis c(-90,-60)
interval_upper_limit        <- 30

monoCTD_first_diag_exam_span_within_interval_limit_imputed_data <- 
  copy(monoCTD_first_diag_exam_span_calculation_data) %>% 
  
  # Step 0: order dataset by the variables to ensure the most recent exam results of each patients' exam items at the first row
  .[order(IDCODE,ITEM_LABSH1IT,overall_Status,first_diag_exam_dist,SCDATE),] %>% 
  
  # Step 1: identify exam results within ± 30 days of first diagnosis date, 
  #         then assign non-NA results as imputed status values
  .[!is.na(overall_Status) & 
      between(first_diag_exam_span,
              interval_lower_limit,interval_upper_limit),
    Imp_Status := overall_Status] %>% 
  
  # Step 2: identify exam results within the range of expand lower limit and upper limit of first diagnosis date, 
  #         then assign non-NA results as imputed status values to those are still NA after Step 1
  .[!is.na(overall_Status) & is.na(Imp_Status) & 
      between(first_diag_exam_span,
              expand_interval_lower_limit,interval_upper_limit),
    Imp_Status := overall_Status] %>% 
  
  # Step 3: remove all exam results out of specified range of intervals (this step would remove subjects without any results within the specified interval)
  .[between(first_diag_exam_span,
            expand_interval_lower_limit,interval_upper_limit),] %>% 
  
  # Step 4: subset data set for first rows of each patients' each exam items
  .[,head(.SD,1L),by = c("IDCODE","ITEM_LABSH1IT")] %>% 
  
  # Step 5: For the patients having none of the exam results in the specified interval (disease specific), assign Abnormal
  .[is.na(Imp_Status) & 
      ITEM_LABSH1IT %in% "Rheumatoid factor (RF)_72-134" & 
      Group %in% "Rheumatoid arthritis",
    Imp_Status := "Abnormal"] %>% 
  .[is.na(Imp_Status) & 
      ITEM_LABSH1IT %in% "Antinuclear antibody (ANA)_72-245" & 
      Group %in% c("Systemic lupus erythematosus","Sjogren's syndrome"),
    Imp_Status := "Abnormal"] %>% 
  # Step 6: For those exam results without test results in the nearest records to
  #         first diagnosed date, assign Normal
  .[is.na(Imp_Status),Imp_Status := "Normal"]
```

<br/>

# Prepare ready-for-analysis dataset

``` r
selected_variables <- 
  c("IDCODE","AGE","SEX","Group","First_diagnosis_date",
    "ITEM_LABSH1IT","SCDATE","Imp_Status")

monoCTD_ready_analysis_dataset <- 
  copy(monoCTD_first_diag_exam_span_within_interval_limit_imputed_data) %>% 
  .[,.SD,.SDcols = selected_variables]
```

<br/>

    ## Classes 'data.table' and 'data.frame':   551376 obs. of  7 variables:
    ##  $ AGE                 : num  64 64 64 64 64 64 64 64 64 64 ...
    ##  $ SEX                 : chr  "F" "F" "F" "F" ...
    ##  $ Group               : Factor w/ 9 levels "Rheumatoid arthritis",..: 4 4 4 4 4 4 4 4 4 4 ...
    ##  $ First_diagnosis_date: Date, format: "2015-05-13" "2015-05-13" ...
    ##  $ ITEM_LABSH1IT       : Factor w/ 36 levels "AQP4 autoantibody_72-292",..: 1 2 3 4 5 6 7 8 9 10 ...
    ##  $ SCDATE              : Date, format: "2015-05-07" "2015-05-07" ...
    ##  $ Imp_Status          : chr  "Normal" "Abnormal" "Normal" "Normal" ...
    ##  - attr(*, ".internal.selfref")=<externalptr>

<br/>

<br/>

![</br>**Figure
1**</br>](./Data_Wrangling_files/Figure%201.%20Eligibility%20and%20the%20number%20of%20study%20subjects.jpeg)
**Figure 1. Eligibility and the number of study subjects.** (RA,
rheumatoid arthritis; SLE, systemic lupus erythematosus; SS, Sjögren’s
syndrome)
