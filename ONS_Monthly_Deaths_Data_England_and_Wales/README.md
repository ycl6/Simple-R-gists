# Deaths registered monthly in England and Wales

I-Hsuan Lin

University of Manchester

January 08, 2022

## Introduction

This notebook shows how to use `readxl` package to retreive the *Deaths registered monthly in England and Wales* Dataset from [Office for National Statistics](https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/monthlyfiguresondeathsregisteredbyareaofusualresidence) and create various plots to show the number of deaths with `ggplot2`.

Source: [Office for National Statistics](https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/monthlyfiguresondeathsregisteredbyareaofusualresidence)

## About this dataset

> Number of deaths registered each month by area of usual residence for England and Wales, by region, county, local and unitary authority, and London borough. These are monthly provisional data covering the month before release and do not include the most up-to-date figures on deaths registered involving coronavirus (COVID-19); see our weekly deaths data.

## Important notes and usage information

> If you are looking for the latest data on deaths involving the coronavirus (COVID-19) registered in England and Wales, please see our [weekly provisional deaths dataset](https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/weeklyprovisionalfiguresondeathsregisteredinenglandandwales).

## Main points from latest release

>- In November 2021, there were 48,180 deaths registered in England, 6,511 deaths (15.6%) more than the November five-year average (2015 to 2019); there were 3,344 deaths registered in Wales, 557 deaths (20.0%) more than the November average.
>- The leading cause of death in November 2021 was dementia and Alzheimer’s disease in England (accounting for 11.8% of all deaths) and ischaemic heart diseases in Wales (accounting for 10.7% of all deaths).
>- Coronavirus (COVID-19) was the third leading cause of death in November 2021, in both England (accounting for 6.6% of all deaths) and Wales (accounting for 9.0% of all deaths).
>- Taking into account the population size and age structure, the age-standardised mortality rate (ASMR) for deaths due to COVID-19 in England increased significantly to 69.3 deaths per 100,000 people; the ASMR for deaths due to COVID-19 in Wales was 106.4 deaths per 100,000 people, which was higher than October 2021 but was not statistically significant.
>- Yorkshire and The Humber remained the English region with the highest ASMR for deaths due to COVID-19 in November 2021 (91.9 deaths per 100,000 people).
