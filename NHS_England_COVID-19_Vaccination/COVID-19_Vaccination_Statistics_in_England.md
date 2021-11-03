# COVID-19 Vaccination Statistics in England

I-Hsuan Lin

University of Manchester

June 20, 2021

## 1. Introduction

This notebook shows how to use `readxl` package to retreive the *Monthly COVID-19 Vaccinations* Dataset from
[NHS England](https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-vaccinations/) and create various plots to show key statistics with `ggplot2`.

### About this dataset

This publication shows the number of vaccinations given to people in England who are eligible for vaccination (individuals who have an NHS number and are currently alive in the resident population). This data covers vaccinations administered up to midnight of the last day of each month.

- **Period**: 8th December 2020 to 31st May 2021
- **Source**: National Immunisation Management Service (NIMS)
- **Published**: 10th June 2021

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
options(repr.plot.width = 12, repr.plot.height = 8, repr.plot.res = 150)
```

## 4. Retrieve datasets

- `n1` - Population Estimates (ONS)
- `n2` - Vaccinations by Region and Age
- `n3` - Vaccinations by NHS Region and Vaccination Date


```R
domain <- "https://www.england.nhs.uk"
path <- "statistics/wp-content/uploads/sites/2/2021/06"
filename <- "COVID-19-monthly-announced-vaccinations-10-June-2021.xlsx"

# Download file
download.file(url = paste0(domain, "/", path, "/", filename), destfile = filename, method = "curl")
```


```R
# Use read_excel to import from Excel file
n1 <- read_excel(filename, sheet = "Population estimates (ONS)", range = cell_rows(c(12:22)), 
                 .name_repair = "minimal")
n2 <- read_excel(filename, sheet = "Region & Age", range = cell_rows(c(11:21)), .name_repair = "minimal")
n3 <- read_excel(filename, sheet = "Vaccination Date", range = cell_rows(c(11:1000)), .name_repair = "minimal")
```


```R
n1 # Population estimates (ONS)
```


```R
n2 # Region & Age
```


```R
head(n3) # Vaccination Date
```

## 5. Set descriptions


```R
# 'Under 16','16-29','30-34','35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-79','80+'
agegroups1 <- n1[1, 3:15] %>% as.character()
names(agegroups1) <- paste0("G", 1:length(agegroups1))
agegroups1 <- factor(agegroups1, levels = agegroups1)
agegroups1

# '16-29','30-34','35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-79','80+'
agegroups2 <- agegroups1[-1]
names(agegroups2) <- paste0("G", 1:length(agegroups2))
agegroups2 <- factor(agegroups2, levels = agegroups2)
agegroups2

regions <- n1[4:10,1] %>% pull()
names(regions) <- paste0("R", 1:length(regions))
regions <- as.factor(regions)
regions
```

## 6. Manipulate `data.frame`

### ONS population mapped to NHS Region (tibble `n1`)


```R
population <- as.data.frame(n1)
# Remove columns and rows with all NAs
population <- population[rowSums(is.na(population)) != ncol(population), 
                         colSums(is.na(population)) != nrow(population)]

# Drop last column with 16+ statistics
population <- population[-(1:2), -ncol(population)]
colnames(population) <- c("Region", names(agegroups1))
population <- mutate_at(population, vars(2:ncol(population)), "as.integer")
rownames(population) <- population$Region
population <- population[regions,]
population

# Drop first column with Under 16 statistics
population2 <- data.frame(Region = population$Region, population[,3:ncol(population)])
colnames(population2) <- c("Region", names(agegroups2))
population2
```

### COVID-19 vaccinations by NHS region of residence and age group (tibble `n2`)


```R
d <- as.data.frame(n2)
# Remove columns and rows with all NAs
d <- d[rowSums(is.na(d)) != ncol(d), colSums(is.na(d)) != nrow(d)]
d
```

#### First dose (region & age group)


```R
# Wide format
d1 <- d[-(1:2), 1:(length(agegroups2)+1)]
colnames(d1) <- c("Region", names(agegroups2))
d1$Region <- stringr::str_to_title(d1$Region)
d1$Region <- factor(d1$Region, levels = levels(regions))
d1

# Convert to long format
d1 <- pivot_longer(d1, !Region, names_to = "AgeGroup", values_to = "Count")
d1$Dose <- "First Dose"
head(d1)
```

#### Second dose (region & age group)


```R
# Wide format
d2 <- d[-(1:2), 14:(length(agegroups2)+13)]
colnames(d2) <- names(agegroups2)
d2 <- cbind(data.frame(Region = factor(levels(d1$Region))), d2)
d2

# Convert to long format
d2 <- pivot_longer(d2, !Region, names_to = "AgeGroup", values_to = "Count")
d2$Dose <- "Second Dose"
head(d2)
```

#### Combine 1st and 2nd doses (region & age group)


```R
data <- rbind(d1, d2)
data <- merge(pivot_longer(population2, !Region, names_to = "AgeGroup", values_to = "Total"), data)
data$Count <- as.integer(data$Count)
data$Fraction <- data$Count / data$Total
data$Region <- factor(data$Region, levels = levels(regions))
data$AgeGroup <- as.factor(data$AgeGroup)
data$AgeGroup <- factor(data$AgeGroup, levels = gtools::mixedsort(levels(data$AgeGroup)))
levels(data$AgeGroup) <- agegroups2
head(data)
```

### COVID-19 cumulative vaccinations by date of vaccination and NHS region of residence (tibble `n3`)


```R
d <- as.data.frame(n3)
# Remove columns and rows with all NAs
d <- d[rowSums(is.na(d)) != ncol(d), colSums(is.na(d)) != nrow(d)]
d <- d[grep("^44", d$"Date of Vaccination", invert = FALSE), ]
head(d)
```

#### First dose (date & region)


```R
# Wide format
d3 <- d[, c(1, 3:(length(regions)+2))]
colnames(d3) <- c("Date", names(regions))
d3$Date <- as.Date(as.integer(d3$Date), origin = "1899-12-30") # Change data format
head(d3)

# Convert to long format
d3 <- pivot_longer(d3, !Date, names_to = "Region", values_to = "Count")
d3$Dose <- "First Dose"
head(d3)
```

#### Second dose (date & region)


```R
# Wide format
d4 <- d[, 11:(length(regions)+10)]
colnames(d4) <- names(regions)
d4 <- cbind(data.frame(Date = as.Date(as.integer(d$"Date of Vaccination"), origin = "1899-12-30")), d4)
head(d4)

# Convert to long format
d4 <- pivot_longer(d4, !Date, names_to = "Region", values_to = "Count")
d4$Dose <- "Second Dose"
head(d4)
```

#### Combine 1st and 2nd doses (date & region)


```R
data2 <- rbind(d3, d4)
data2$Region <- factor(data2$Region, levels = gtools::mixedsort(unique(data2$Region)))
levels(data2$Region) <- regions
data2 <- merge(data.frame(Region = population2$Region, Total = rowSums(population2[,-1])), data2)
data2$Count <- as.integer(data2$Count)
data2$Fraction <- data2$Count / data2$Total
head(data2)
```

## 7. Create plots

### ONS population estimates in England


```R
pivot_longer(population, !Region, names_to = "AgeGroup", values_to = "Count") %>% 
    mutate(AgeGroup = factor(AgeGroup, levels = unique(AgeGroup))) %>%
    ggplot(aes(AgeGroup, Count)) + geom_col(fill = "black") + 
    facet_wrap(~ Region, ncol = 4) + theme_classic(base_size = 16) +
    scale_x_discrete(labels = levels(agegroups1)) + # Change X-axis labels
    scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = "ONS population in England", x = "Age group")
```

### Cumulative COVID-19 doses by age group


```R
total <- data %>% summarise(sum(Count)) %>% pull()
total <- label_comma(accuracy = 1)(total)

data %>% group_by(AgeGroup) %>% summarise(Doses = sum(Count)) %>%
    ggplot(aes(AgeGroup, Doses)) + geom_col(fill = "black") + theme_classic(base_size = 16) +
    scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) + 
    labs(title = "Cumulative Doses from All Regions in England", subtitle = paste("Total:", total), 
         x = "Age group")
```

### Total COVID-19 doses by NHS region of residence and age group


```R
ggplot(data, aes(AgeGroup, Count)) + geom_col(fill = "black", width = 0.9) + 
    facet_wrap(~ Region, ncol = 4) + theme_classic(base_size = 16) +
    scale_y_continuous(labels = unit_format(unit = "K", scale = 1e-3)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = "Total Doses in England", x = "Age group", y = "Doses")
```

### Total COVID-19 doses by NHS region of residence and age group, coloured by 1st and 2nd dose


```R
ggplot(data, aes(AgeGroup, Count, fill = Dose)) + 
    geom_col(alpha = 0.9, position = position_dodge2(width = 0.9, padding = 0.1)) +
    facet_wrap(~ Region, ncol = 4) + theme_classic(base_size = 16) +
    scale_y_continuous(labels = unit_format(unit = "K", scale = 1e-3)) + 
    scale_fill_manual(values = c("cyan","red")) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = "Total First and Second Doses in England", x = "Age group", y = "Doses")
```

### Percentage of estimated population receiving COVID-19 vaccination, coloured by 1st and 2nd dose


```R
ggplot(data, aes(AgeGroup, Fraction, fill = Dose)) + 
    geom_col(alpha = 0.9, position = position_dodge2(width = 0.9, padding = 0.1)) + 
    facet_wrap(~ Region, ncol = 4) + theme_classic(base_size = 16) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_fill_manual(values = c("cyan","red")) + geom_hline(yintercept = 1, linetype = 2) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = "Total First and Second Doses in England", x = "Age group", y = "% Population")
```

### COVID-19 cumulative vaccinations by date of vaccination


```R
data2 %>% group_by(Date, Dose) %>% summarise(Count = sum(Count), .groups = "drop") %>%
    ggplot(aes(Date, Count, fill = Dose)) + geom_area(position = position_dodge(width = 0)) + 
    theme_classic(base_size = 16) +
    scale_x_date(date_labels = "%d %b", date_breaks = "1 week") +
    scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
    scale_fill_manual(values = c("cyan","red")) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = "Total First and Second Doses in England", x = "Date", y = "Doses")
```

### COVID-19 cumulative vaccinations by date of vaccination and NHS region of residence


```R
ggplot(data2, aes(Date, Count, color = Region)) + geom_line(size = 1) + 
    facet_wrap(~ Dose, ncol = 2) + theme_classic(base_size = 16) +
    scale_x_date(date_labels = "%d %b", date_breaks = "2 week") +
    scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
    scale_color_brewer(palette = "Dark2") + guides(color = guide_legend(override.aes = list(size = 5))) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), panel.spacing = unit(1, "lines")) +
    labs(title = "Cumulative vaccinations by date", x = "Date", y = "Doses")
```

### Percentage of estimated population receiving COVID-19 vaccination, coloured by region


```R
ggplot(data2, aes(Date, Fraction, color = Region)) + geom_line(size = 1) + 
    facet_wrap(~ Dose, ncol = 2) + theme_classic(base_size = 16) +
    scale_x_date(date_labels = "%d %b", date_breaks = "2 week") +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_color_brewer(palette = "Dark2") + guides(color = guide_legend(override.aes = list(size = 5))) +
    geom_hline(yintercept = seq(0.2, 0.8, 0.2), linetype = 2) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), panel.spacing = unit(1, "lines")) +
    labs(title = "Cumulative vaccinations by date", x = "Date", y = "% Population")
```

## 8. Session info


```R
sessionInfo()
```
