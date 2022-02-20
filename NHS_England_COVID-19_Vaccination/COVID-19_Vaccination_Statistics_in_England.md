# COVID-19 Vaccination Statistics in England

I-Hsuan Lin

University of Manchester

February 20, 2022

## 1. Introduction

This notebook shows how to use `readxl` package to retreive the *Monthly COVID-19 Vaccinations* Dataset from
[NHS England](https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-vaccinations/) and create various plots to show key statistics with `ggplot2`.

### About this dataset

The majority of the figures in this publication provide information on individuals who are eligible for vaccination (those who have an NHS number and are currently alive) who have been vaccinated for COVID-19 in England or who live in England but have been vaccinated for COVID-19 outside of England.

The publication also includes figures based on all vaccinations administered in England, even if individuals vaccinated are resident outside of England, do not have an NHS number or are no longer alive.

All data in the monthly publication covers vaccinations administered up to midnight of the last day of the previous month.

- **Period**: 8th December 2020 to 31st January 2022
- **Source**: National Immunisation Management System (NIMS)
- **Basis**: England
- **Published**: 10th February 2022

## 2. Loading required libraries


```R
library(readxl)
library(ggplot2)
library(scales)
library(tidyverse)
```

## 3. Set output parameters


```R
# Set width
options(width = 110)

# Set output image size
options(repr.plot.width = 12, repr.plot.height = 10, repr.plot.res = 150)
```

## 4. Retrieve datasets

- `n1` - ONS Population Estimates - Mid-year 2020
- `n2` - COVID-19 Vaccinations By NHS Region of Residence and Age Group
   - `n2a` - 1st dose
   - `n2b` - 2nd dose
   - `n2c` - Booster or 3rd dose
- `n3` - COVID-19 Cumulative Vaccinations by Date of Vaccination and NHS Region of Residence
- `n4` - NIMS Population Estimates by Ethnicity
- `n5` - COVID-19 Vaccinations By Ethnicity and Age Group


```R
domain <- "https://www.england.nhs.uk"
path <- "statistics/wp-content/uploads/sites/2/2022/02"
filename <- "COVID-19-monthly-announced-vaccinations-10-February-2022.xlsx"

# Download file
download.file(url = paste0(domain, "/", path, "/", filename), destfile = filename, method = "curl")
```


```R
# Use read_excel to import from Excel file
n1 <- read_excel(filename, sheet = "11. Pop estimates (ONS 2020)", range = cell_rows(c(12:22)), 
                 .name_repair =  ~ make.names(.x, unique = TRUE)) %>% filter_all(any_vars(!is.na(.)))

n2a <- read_excel(filename, sheet = "1. Region & Age", range = cell_rows(c(23:33)), 
                  .name_repair =  ~ make.names(.x, unique = TRUE)) %>% filter_all(any_vars(!is.na(.)))
n2b <- read_excel(filename, sheet = "1. Region & Age", range = cell_rows(c(35:45)), 
                  .name_repair =  ~ make.names(.x, unique = TRUE)) %>% filter_all(any_vars(!is.na(.)))
n2c <- read_excel(filename, sheet = "1. Region & Age", range = cell_rows(c(47:57)), 
                  .name_repair =  ~ make.names(.x, unique = TRUE)) %>% filter_all(any_vars(!is.na(.)))

n3 <- read_excel(filename, sheet = "9. Vaccination Date", range = cell_rows(c(12:1000)), 
                 .name_repair =  ~ make.names(.x, unique = TRUE)) %>% filter_all(any_vars(!is.na(.)))

n4 <- read_excel(filename, sheet = "12. Pop estimates (NIMS)", range = cell_rows(c(11:31)), 
                 .name_repair =  ~ make.names(.x, unique = TRUE)) %>% filter_all(any_vars(!is.na(.)))

n5 <- read_excel(filename, sheet = "4. Ethnicity, Region & Age", range = cell_rows(c(11:32)), 
                 .name_repair =  ~ make.names(.x, unique = TRUE)) %>% filter_all(any_vars(!is.na(.)))
```

### ONS Population Estimates - Mid-year 2020


```R
n1
```

### COVID-19 Vaccinations By NHS Region of Residence and Age Group


```R
n2a
n2b
n2c
```

### COVID-19 Cumulative Vaccinations by Date of Vaccination and NHS Region of Residence


```R
head(n3)
```

### NIMS Population Estimates by Ethnicity


```R
n4
```

### COVID-19 Vaccinations By Ethnicity and Age Group


```R
n5
```

## 5. Set descriptions


```R
# '12-15','16-17','18-24','25-29','30-34','35-39','40-44','45-49','50-54','55-59','60-64','65-69',
# '70-74','75-79','80+'
ageGroup1 <- n1[1, 3:17] %>% as.character()
names(ageGroup1) <- paste0("G", 1:length(ageGroup1))
ageGroup1 <- factor(ageGroup1, levels = ageGroup1)
ageGroup1

# '18-24','25-29','30-34','35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-79','80+'
ageGroup2 <- ageGroup1[!ageGroup1 %in% c("12-15","16-17")]
ageGroup2 <- droplevels(ageGroup2)
ageGroup2

regionName <- n1$NHS.Region.name[3:9]
regionCode <- n1$NHS.Region.code[3:9]
names(regionName) <- regionCode
regionName <- as.factor(regionName)
regionName

ethnicity <- n4[3:nrow(n4), 1] %>% pull()
names(ethnicity) <- paste0("E", seq_along(ethnicity))
ethnicity <- factor(ethnicity, levels = ethnicity)
ethnicity
```

## 6. Manipulate `data.frame`

### ONS Population Estimates - Mid-year 2020 (tibble `n1`)


```R
# Wide format
population <- as.data.frame(n1[n1$NHS.Region.code %in% regionCode,]) %>% select(-NHS.Region.code) # drop 1st column
# Drop last 4 columns: 50+, 18+, 16+, 12+
population <- population[,1:(length(ageGroup1)+1)]
colnames(population) <- c("Region", names(ageGroup1))
population$Region <- factor(population$Region, levels = levels(regionName))
population <- mutate_at(population, vars(-Region), "as.integer")
population

# Convert to long format
population <- pivot_longer(population, !Region, names_to = "AgeGroup", values_to = "Total")
population$AgeGroup <- as.factor(population$AgeGroup)
population$AgeGroup <- factor(population$AgeGroup, levels = gtools::mixedsort(levels(population$AgeGroup)))
levels(population$AgeGroup) <- ageGroup1
head(population)
```

### COVID-19 Vaccinations By NHS Region of Residence and Age Group (tibble `n2a`, `n2b`, `n2c`)

#### First dose (region & age group)


```R
# Wide format
ra1 <- as.data.frame(n2a[n2a$NHS.Region.code %in% regionCode,]) %>% select(-NHS.Region.code) # drop 1st column
colnames(ra1) <- c("Region", names(ageGroup1))
ra1$Region <- factor(ra1$Region, levels = levels(regionName))
ra1 <- mutate_at(ra1, vars(-Region), "as.integer")
ra1

# Convert to long format
ra1 <- pivot_longer(ra1, !Region, names_to = "AgeGroup", values_to = "Count")
ra1$AgeGroup <- as.factor(ra1$AgeGroup)
ra1$AgeGroup <- factor(ra1$AgeGroup, levels = gtools::mixedsort(levels(ra1$AgeGroup)))
levels(ra1$AgeGroup) <- ageGroup1
ra1$Dose <- "First Dose"
head(ra1)
```

#### Second dose (region & age group)


```R
# Wide format
ra2 <- as.data.frame(n2b[n2b$NHS.Region.code %in% regionCode,]) %>% select(-NHS.Region.code) # drop 1st column
colnames(ra2) <- c("Region", names(ageGroup1))
ra2$Region <- factor(ra2$Region, levels = levels(regionName))
ra2 <- mutate_at(ra2, vars(-Region), "as.integer")
ra2

# Convert to long format
ra2 <- pivot_longer(ra2, !Region, names_to = "AgeGroup", values_to = "Count")
ra2$AgeGroup <- as.factor(ra2$AgeGroup)
ra2$AgeGroup <- factor(ra2$AgeGroup, levels = gtools::mixedsort(levels(ra2$AgeGroup)))
levels(ra2$AgeGroup) <- ageGroup1
ra2$Dose <- "Second Dose"
head(ra2)
```

#### Booster or 3rd dose (region & age group)


```R
# Wide format
ra3 <- as.data.frame(n2c[n2c$NHS.Region.code %in% regionCode,]) %>% select(-NHS.Region.code) # drop 1st column
ra3 <- ra3[,-2] # Drop "Under 18"
colnames(ra3) <- c("Region", names(ageGroup2))
ra3$Region <- factor(ra3$Region, levels = levels(regionName))
ra3 <- mutate_at(ra3, vars(-Region), "as.integer")
ra3

# Convert to long format
ra3 <- pivot_longer(ra3, !Region, names_to = "AgeGroup", values_to = "Count")
ra3$AgeGroup <- as.factor(ra3$AgeGroup)
ra3$AgeGroup <- factor(ra3$AgeGroup, levels = gtools::mixedsort(levels(ra3$AgeGroup)))
levels(ra3$AgeGroup) <- ageGroup2
ra3$Dose <- "Booster/Third Dose"
head(ra3)
```

#### Combine 1st, 2nd, 3rd/booster doses (region & age group)


```R
dat1 <- rbind(ra1, ra2, ra3)
dat1 <- merge(population, dat1) %>% arrange(Region, order(gtools::mixedorder(AgeGroup)))
dat1$Fraction <- dat1$Count / dat1$Total
dat1$Dose <- as.factor(dat1$Dose)
dat1$Dose <- factor(dat1$Dose, levels = c("First Dose","Second Dose","Booster/Third Dose"))
head(dat1)
```

### COVID-19 Cumulative Vaccinations by Date of Vaccination and NHS Region of Residence (tibble `n3`)


```R
dr <- as.data.frame(n3) %>% rename(Date = 1)
# Remove columns and rows with all NAs
dr <- dr[rowSums(is.na(dr)) != ncol(dr), colSums(is.na(dr)) != nrow(dr)]
dr <- dr[grep("^44", dr$Date, invert = FALSE), ]
head(dr)
```

#### First dose (date & region)


```R
# Wide format
idx <- grep("1st.dose", colnames(dr)) # 1st.dose
dr1 <- dr[, c(1, idx:(length(regionName)+(idx-1)))]
colnames(dr1) <- c("Date", as.character(regionName))
dr1$Date <- as.Date(as.integer(dr1$Date), origin = "1899-12-30") # Change data format
dr1 <- mutate_at(dr1, vars(-Date), "as.integer")
head(dr1)

# Convert to long format
dr1 <- pivot_longer(dr1, !Date, names_to = "Region", values_to = "Count")
dr1$Region <- factor(dr1$Region, levels = levels(regionName))
dr1$Dose <- "First Dose"
head(dr1)
```

#### Second dose (date & region)


```R
# Wide format
idx <- grep("2nd.dose", colnames(dr)) # 2nd.dose
dr2 <- dr[, c(1, idx:(length(regionName)+(idx-1)))]
colnames(dr2) <- c("Date", as.character(regionName))
dr2$Date <- as.Date(as.integer(dr2$Date), origin = "1899-12-30") # Change data format
dr2 <- mutate_at(dr2, vars(-Date), "as.integer")
head(dr2)

# Convert to long format
dr2 <- pivot_longer(dr2, !Date, names_to = "Region", values_to = "Count")
dr2$Region <- factor(dr2$Region, levels = levels(regionName))
dr2$Dose <- "Second Dose"
head(dr2)
```

#### Booster or 3rd dose (date & region)


```R
# Wide format
idx <- grep("Booster.or.3rd.dose", colnames(dr)) # Booster.or.3rd.dose
dr3 <- dr[, c(1, idx:(length(regionName)+(idx-1)))]
colnames(dr3) <- c("Date", as.character(regionName))
dr3$Date <- as.Date(as.integer(dr3$Date), origin = "1899-12-30") # Change data format
dr3 <- mutate_at(dr3, vars(-Date), "as.integer")
head(dr3)

# Convert to long format
dr3 <- pivot_longer(dr3, !Date, names_to = "Region", values_to = "Count")
dr3$Region <- factor(dr3$Region, levels = levels(regionName))
dr3$Dose <- "Booster/Third Dose"
head(dr3)
```

#### Combine 1st, 2nd, 3rd/booster doses (date & region)


```R
dat2 <- rbind(dr1, dr2, dr3)
totalPop <- population %>% group_by(Region) %>% summarise(Total = sum(Total))
dat2 <- merge(totalPop, dat2) %>% arrange(Region, Date)
dat2$Fraction <- dat2$Count / dat2$Total
dat2$Dose <- factor(dat2$Dose, levels = c("First Dose","Second Dose","Booster/Third Dose"))
head(dat2)
```

### NIMS Population Estimates by Ethnicity (tibble `n4`)


```R
# Wide format
epop <- as.data.frame(n4) %>% rename("Ethnicity" = 1)
epop <- epop[epop$Ethnicity %in% ethnicity,]
# Remove columns and rows with all NAs
epop <- epop[rowSums(is.na(epop)) != ncol(epop), colSums(is.na(epop)) != nrow(epop)]
# Drop last 4 columns: 50+, 18+, 16+, 12+
epop <- epop[,1:(length(ageGroup1)+1)]
colnames(epop) <- c("Ethnicity", names(ageGroup1))
epop <- mutate_at(epop, vars(-Ethnicity), "as.integer")
head(epop)

# Convert to long format
epop <- pivot_longer(epop, !Ethnicity, names_to = "AgeGroup", values_to = "Total")
epop$AgeGroup <- as.factor(epop$AgeGroup)
epop$AgeGroup <- factor(epop$AgeGroup, levels = gtools::mixedsort(levels(epop$AgeGroup)))
levels(epop$AgeGroup) <- ageGroup1
head(epop)
```

### COVID-19 Vaccinations By Ethnicity and Age Group (tibble `n5`)


```R
ea <- as.data.frame(n5) %>% rename(Ethnicity = 1)
ea <- ea[ea$Ethnicity %in% ethnicity,]
# Remove columns and rows with all NAs
ea <- ea[rowSums(is.na(ea)) != ncol(ea), colSums(is.na(ea)) != nrow(ea)]
head(ea)
```

#### First dose (ethnicity & age group)


```R
# Wide format
idx <- grep("1st.dose", colnames(ea)) # 1st.dose
ea1 <- ea[, c(1, idx:(length(ageGroup1)+(idx-1)))]
colnames(ea1) <- c("Ethnicity", names(ageGroup1))
ea1 <- mutate_at(ea1, vars(-Ethnicity), "as.integer")
head(ea1)

# Convert to long format
ea1 <- pivot_longer(ea1, !Ethnicity, names_to = "AgeGroup", values_to = "Count")
ea1$AgeGroup <- as.factor(ea1$AgeGroup)
ea1$AgeGroup <- factor(ea1$AgeGroup, levels = gtools::mixedsort(levels(ea1$AgeGroup)))
levels(ea1$AgeGroup) <- ageGroup1
ea1$Dose <- "First Dose"
head(ea1)
```

#### Second dose (ethnicity & age group)


```R
# Wide format
idx <- grep("2nd.dose", colnames(ea)) # 2nd.dose
ea2 <- ea[, c(1, idx:(length(ageGroup1)+(idx-1)))]
colnames(ea2) <- c("Ethnicity", names(ageGroup1))
ea2 <- mutate_at(ea2, vars(-Ethnicity), "as.integer")
head(ea2)

# Convert to long format
ea2 <- pivot_longer(ea2, !Ethnicity, names_to = "AgeGroup", values_to = "Count")
ea2$AgeGroup <- as.factor(ea2$AgeGroup)
ea2$AgeGroup <- factor(ea2$AgeGroup, levels = gtools::mixedsort(levels(ea2$AgeGroup)))
levels(ea2$AgeGroup) <- ageGroup1
ea2$Dose <- "Second Dose"
head(ea2)
```

#### Booster or 3rd dose (ethnicity & age group)


```R
# Wide format
idx <- grep("Booster.or.3rd.dose", colnames(ea)) # Booster.or.3rd.dose
ea3 <- ea[, c(1, idx:(length(ageGroup1)+(idx-2)))]
ea3 <- ea3[,-2] # Drop "Under 18"
colnames(ea3) <- c("Ethnicity", names(ageGroup2))
ea3 <- mutate_at(ea3, vars(-Ethnicity), "as.integer")
head(ea3)

# Convert to long format
ea3 <- pivot_longer(ea3, !Ethnicity, names_to = "AgeGroup", values_to = "Count")
ea3$AgeGroup <- as.factor(ea3$AgeGroup)
ea3$AgeGroup <- factor(ea3$AgeGroup, levels = gtools::mixedsort(levels(ea3$AgeGroup)))
levels(ea3$AgeGroup) <- ageGroup2
ea3$Dose <- "Booster/Third Dose"
head(ea3)
```

#### Combine 1st, 2nd, 3rd/booster doses (ethnicity & age group)


```R
head(epop)
head(ea1)
head(ea2)
head(ea3)
```


```R
dat3 <- rbind(ea1, ea2, ea3)
dat3 <- merge(epop, dat3) %>% arrange(Ethnicity, order(gtools::mixedorder(AgeGroup)))
dat3$Fraction <- dat3$Count / dat3$Total
dat3$AgeGroup <- as.factor(dat3$AgeGroup)
dat3$AgeGroup <- factor(dat3$AgeGroup, levels = gtools::mixedsort(levels(dat3$AgeGroup)))
levels(dat3$AgeGroup) <- ageGroup1
dat3$Dose <- as.factor(dat3$Dose)
dat3$Dose <- factor(dat3$Dose, levels = c("First Dose","Second Dose","Booster/Third Dose"))
head(dat3)
```

## 7. Create plots

### ONS population estimates in England


```R
population %>% ggplot(aes(AgeGroup, Total)) + geom_col(fill = "black") + 
    facet_wrap(~ Region, ncol = 4) + theme_classic(base_size = 16) +
    scale_y_continuous(labels = unit_format(unit = "K", scale = 1e-3)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = "ONS population in England", x = "Age group")
```

### Cumulative COVID-19 doses by age group


```R
total <- dat1 %>% summarise(sum(Count)) %>% pull()
total <- label_comma(accuracy = 1)(total)

dat1 %>% group_by(AgeGroup, Dose) %>% summarise(Doses = sum(Count)) %>%
    ggplot(aes(AgeGroup, Doses)) + geom_col(fill = "black") + theme_classic(base_size = 16) +
    scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) + 
    labs(title = "Cumulative Doses from All Regions in England", 
         subtitle = paste("Total:", total), x = "Age group")
```

### Cumulative COVID-19 doses by NHS region of residence and age group


```R
ggplot(dat1, aes(AgeGroup, Count)) + geom_col(fill = "black", width = 0.8) + 
    facet_wrap(~ Region, ncol = 4) + theme_classic(base_size = 16) +
    scale_y_continuous(labels = unit_format(unit = "K", scale = 1e-3)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = "Total Doses in England", x = "Age group", y = "Doses")
```

### Total COVID-19 doses by NHS region of residence and age group, coloured by 1st and 2nd dose


```R
ggplot(dat1, aes(AgeGroup, Count, fill = Dose)) + 
    geom_col(position = position_dodge2(width = 0.9, padding = 0.1, preserve = "single")) +
    facet_wrap(~ Region, ncol = 3) + theme_classic(base_size = 16) +
    scale_y_continuous(labels = unit_format(unit = "K", scale = 1e-3)) + 
    scale_fill_manual(values = c("cyan","red","blue"), labels = c("1st","2nd","3rd/Booster")) + 
    theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = "Total First, Second, Third/Booster Doses in England", x = "Age group", y = "Doses")
```

### Percentage of estimated population receiving COVID-19 vaccination, coloured by 1st and 2nd dose


```R
ggplot(dat1, aes(AgeGroup, Fraction, fill = Dose)) + 
    geom_col(position = position_dodge2(width = 0.9, padding = 0.1, preserve = "single")) + 
    facet_wrap(~ Region, ncol = 3) + theme_classic(base_size = 16) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) + geom_hline(yintercept = 0.75, linetype = 2) +
    scale_fill_manual(values = c("cyan","red","blue"), labels = c("1st","2nd","3rd/Booster")) + 
    theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = "Total First, Second, Third/Booster Doses in England", x = "Age group", y = "% Population")
```

### NIMS population estimates by ethnicity in England


```R
epop %>% ggplot(aes(AgeGroup, Total)) + geom_col(fill = "black") + 
    facet_wrap(~ Ethnicity, ncol = 3, scales = "free_y") + theme_classic(base_size = 12) +
    scale_y_continuous(labels = unit_format(unit = "K", scale = 1e-3)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = "ONS population in England", x = "Age group")
```

### Percentage of estimated population by ethnicity receiving COVID-19 vaccination, coloured by 1st and 2nd dose


```R
ggplot(dat3, aes(AgeGroup, Fraction, fill = Dose)) + 
    geom_col(position = position_dodge2(width = 0.9, padding = 0.1, preserve = "single")) + 
    facet_wrap(~ Ethnicity, ncol = 3) + theme_classic(base_size = 12) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) + geom_hline(yintercept = 0.75, linetype = 2) +
    scale_fill_manual(values = c("cyan","red","blue"), labels = c("1st","2nd","3rd/Booster")) + 
    theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = "Total First, Second, Third/Booster Doses in England", x = "Age group", y = "% Population")
```

### COVID-19 cumulative vaccinations by date of vaccination


```R
dat2 %>% group_by(Date, Dose) %>% summarise(Count = sum(Count), .groups = "drop") %>%
    ggplot(aes(Date, Count, fill = Dose)) + geom_area(position = position_dodge(width = 0)) + 
    theme_classic(base_size = 16) +
    scale_x_date(date_labels = "%d %b %y", date_breaks = "2 week") +
    scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
    scale_fill_manual(values = c("cyan","red","blue"), labels = c("1st","2nd","3rd/Booster")) + 
    theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = "Total First, Second, Third/Booster Doses in England", x = "Date", y = "Doses")
```

### COVID-19 cumulative vaccinations by date of vaccination and NHS region of residence


```R
ggplot(dat2, aes(Date, Count, color = Region)) + geom_line(size = 1) + 
    facet_wrap(~ Dose, nrow = 1) + theme_classic(base_size = 16) +
    scale_x_date(date_labels = "%d %b %y", date_breaks = "4 week") +
    scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
    scale_color_brewer(palette = "Dark2") + guides(color = guide_legend(override.aes = list(size = 5))) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
          panel.spacing = unit(1, "lines")) +
    labs(title = "Cumulative vaccinations by date", x = "Date", y = "Doses")
```

### Percentage of estimated population receiving COVID-19 vaccination, coloured by region


```R
ggplot(dat2, aes(Date, Fraction, color = Region)) + geom_line(size = 1) + 
    facet_wrap(~ Dose, nrow = 1) + theme_classic(base_size = 16) +
    scale_x_date(date_labels = "%d %b %y", date_breaks = "4 week") +
    scale_y_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
    scale_color_brewer(palette = "Dark2") + guides(color = guide_legend(override.aes = list(size = 5))) +
    geom_hline(yintercept = seq(0.25, 0.75, 0.25), linetype = 2) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
          panel.spacing = unit(1, "lines")) +
    labs(title = "Cumulative vaccinations by date", x = "Date", y = "% Population")
```

## 8. Session info


```R
sessionInfo()
```
