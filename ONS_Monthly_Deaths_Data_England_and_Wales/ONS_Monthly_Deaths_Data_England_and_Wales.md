# Deaths registered monthly in England and Wales

I-Hsuan Lin

University of Manchester

February 19, 2022

## 1. Introduction

This notebook shows how to use `readxl` package to retreive the *Deaths registered monthly in England and Wales* Dataset from [Office for National Statistics](https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/monthlyfiguresondeathsregisteredbyareaofusualresidence) and create various plots to show the number of deaths with `ggplot2`.

Source: [Office for National Statistics](https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/monthlyfiguresondeathsregisteredbyareaofusualresidence)

### About this dataset

> Number of deaths registered each month by area of usual residence for England and Wales, by region, county, local and unitary authority, and London borough. These are monthly provisional data covering the month before release and do not include the most up-to-date figures on deaths registered involving coronavirus (COVID-19); see our weekly deaths data.

### Important notes and usage information

> If you are looking for the latest data on deaths involving the coronavirus (COVID-19) registered in England and Wales, please see our [weekly provisional deaths dataset](https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/weeklyprovisionalfiguresondeathsregisteredinenglandandwales).

### Main points from latest release

>- The provisional number of deaths registered in England and Wales in December 2021 was 52,859; this represents an increase of 1,257 deaths in comparison with the previous month and a decrease of 3,813 deaths in comparison with the same month in 2020.
>- Moveable public holidays and the number of weekends, when register offices are closed, affect the number of registrations made in the published months and in the corresponding months in previous years.
>- Local authorities' codes and names have been updated to reflect the changes that occurred in April 2021.

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

### Set file URLs


```R
ons <- "https://www.ons.gov.uk/file?uri=/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/monthlyfiguresondeathsregisteredbyareaofusualresidence"

f2006 <- paste0(ons,"/2006/publishedoutput200tcm772274163.xls")
f2007 <- paste0(ons,"/2007/publishedoutput200tcm772274233.xls")
f2008 <- paste0(ons,"/2008/publishedoutput200tcm772274292.xls")
f2009 <- paste0(ons,"/2009/publishedoutput200tcm772274362.xls")
f2010 <- paste0(ons,"/2010/publishedoutputfeb021tcm772274383.xls")
f2011 <- paste0(ons,"/2011/publishedoutput2011finaltcm772738151.xls")
f2012 <- paste0(ons,"/2012/publishedoutput2012finaltcm773197501.xls")
f2013 <- paste0(ons,"/2013/publishedoutput2013finaltcm773717241.xls")
f2014 <- paste0(ons,"/2014/publishedoutput2014finaltcm774115982.xls")
f2015 <- paste0(ons,"/2015/publishedoutput2015final.xls")
f2016 <- paste0(ons,"/2016/publishedoutput2016final.xls")
f2017 <- paste0(ons,"/2017/publishedoutputannual2017final.xls")
f2018 <- paste0(ons,"/2018/publishedannual2018.xls")
f2019 <- paste0(ons,"/2019/annual2019publishedoutputrefresh.xls")
f2020 <- paste0(ons,"/2020/annual2020publishedoutputrefresh.xls")
f2021 <- paste0(ons,"/2021/deathsregisteredmonthlyusualareaofresidenceenglandandwales.xlsx") # Release date: 21 January 2022
```

### Download datasets


```R
years <- 2006:2021
files <- c(f2006, f2007, f2008, f2009, f2010, f2011, f2012, f2013, 
           f2014, f2015, f2016, f2017, f2018, f2019, f2020, f2021)
select <- "TOTAL REGISTRATIONS|ENGLAND, WALES AND ELSEWHERE"

data <- data.frame(matrix(NA, nrow = 0, ncol = 12), stringsAsFactors = FALSE)

for(i in 1:length(files)) {
    file <- files[i]
    filename <- paste0("Y", years[i], ".", tools::file_ext(file))
    
    if(!file.exists(filename)) {
        # Download excel file using command line tool curl
        download.file(url = file, destfile = filename, method = "curl")
    }
    
    # Pick relevant data row from the the selected sheet, and store as data.frame
    n <- read_excel(filename, sheet = paste0("Figures for ", years[i]),
                        .name_repair =  ~ make.names(.x, unique = TRUE)) %>% rename(Contents = 1) %>% 
    filter(grepl(select, X) | grepl(select, Contents)) %>% select_if(~any(grepl("^[0-9]+$", .)))
    n <- as.numeric(n)
    num <- length(n)
    if(num < 12) {
        n <- c(n, rep(NA, 12-num)) # Add NA if missing data
    }
    data <- rbind(data, n)
}
```

If you are using Jupyter Notebook, please check on your terminal that the download has completed at this stage, e.g.:

```
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100  206k    0  206k    0     0  3385k      0 --:--:-- --:--:-- --:--:-- 3385k
```

## 5. Manipulate `data.frame`

### Add row and column names


```R
names(data) = month.abb
rownames(data) = years

# Print data.frame
data
```

### Reshape to the long format


```R
d <- data %>% rownames_to_column(var = "Year") %>%     # add rowname to column
    gather(key = "Month", value = "Deaths", -Year) %>% # Convert to long format
    mutate(across(!Deaths, as.factor)) %>%             # Convert character columns to factor type
    mutate(Month = factor(Month, levels = month.abb))  # relevel month to the correct order

# Pre-2020 or not
d$Pre2020 <- ifelse(as.Date(ISOdate(years, 1, 1)) < "2020-01-01", TRUE, FALSE)
head(d)
```

### Calculate total deaths per year


```R
total <- d %>% group_by(Year, Pre2020) %>%   # Group data by Year and Pre2020
    summarise(Total = sum(Deaths, na.rm = TRUE), .groups = "drop")

total
```

### Calculate monthly mean, upper and lower bounds in pre-2020 years


```R
# Calculate mean, upper and lower bounds
stat <- d[d$Pre2020 == TRUE,] %>% group_by(Month) %>%   # Group data by Month
    summarise(Mean = mean(Deaths), SD = sd(Deaths)) %>% # Calculate month-wise mean and sd
    mutate(Upper = Mean+SD, Lower = Mean-SD)            # Calculate upper and lower bounds

stat
```

## 6. Create plots


```R
two_color <- c("cyan", "red")
two_color <- setNames(two_color, c(TRUE, FALSE))

plot_title1 <- "Deaths registered yearly in England and Wales 2006 - 2021"
plot_title2 <- "Deaths registered monthly in England and Wales 2006 - 2021"
plot_subtitle <- "Area code: K04000001, J99000001"
plot_caption <- "Dataset from Office for National Statistics"
```

### Show yearly total deaths


```R
ggplot(total, aes(Year, Total, fill = Pre2020)) + 
    geom_col(size = 1, color = NA) + 
    scale_y_continuous(labels = unit_format(unit = "K", scale = 1e-3)) +
    scale_fill_manual(values = two_color) +
    theme_classic(base_size = 18) + labs(title = plot_title1, subtitle = plot_subtitle, caption = plot_caption)
```

### Show monthly total deaths (points and pointrange)


```R
p <- ggplot(d, aes(Month, Deaths, color = Pre2020, size = Pre2020)) + 
    # Add point layer, then add pointrange layer next
    geom_point(aes(shape = Pre2020, alpha = Pre2020)) + 
    geom_pointrange(data = stat, aes(x = Month, y = Mean, ymin = Lower, ymax = Upper), 
                    color = "black", size = 2, alpha = 0.5) +
    scale_y_continuous(labels = unit_format(unit = "K", scale = 1e-3)) +
    scale_color_manual(values = two_color) + scale_alpha_manual(values = c(0.7, 1)) +
    scale_shape_manual(values = c(19, 1)) + scale_size_manual(values = c(5, 3)) +
    guides(size = guide_legend(override.aes = list(size = 5))) +
    theme_classic(base_size = 18) + labs(title = plot_title2, subtitle = plot_subtitle, caption = plot_caption)

p
```

### Show monthly total deaths (points and lines)


```R
p <- ggplot(d, aes(Month, Deaths, group = Year, color = Pre2020, size = Pre2020)) + 
    # Add line layer, then add point layer next
    geom_line() + geom_point(color = "black", size = 2) +
    scale_y_continuous(labels = unit_format(unit = "K", scale = 1e-3)) +
    scale_color_manual(values = two_color) + scale_size_manual(values = c(1.2, 0.8)) +
    theme_classic(base_size = 18) + labs(title = plot_title2, subtitle = plot_subtitle, caption = plot_caption)

p
```

### Show monthly total deaths (points, lines and smoothings)


```R
p <- ggplot(d, aes(Month, Deaths, group = Year, color = Pre2020, size = Pre2020)) + 
    # Add smoothings as the bottom layer, followed by lines and then points
    geom_smooth(data = d[grep("201", d$Year),], fill = "darkgreen", alpha = 0.1, color = "darkgreen", size = 0.6) +
    geom_line() + geom_point(color = "black", size = 2) + 
    scale_y_continuous(labels = unit_format(unit = "K", scale = 1e-3)) +
    scale_color_manual(values = two_color) + scale_size_manual(values = c(1.2, 0.8)) +
    theme_classic(base_size = 18) + labs(title = plot_title2, subtitle = plot_subtitle, caption = plot_caption)

p
```

### Save last image to file


```R
png("ONS_Monthly_Deaths_Data_England_and_Wales.png", width = 9, height = 6, units = "in", res = 150)
p
dev.off()
```

## 7. Session info


```R
sessionInfo()
```
