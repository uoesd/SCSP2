library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
DA <- read_excel("SCS_BAC_and_BrAC_split_TOP.xlsx")

DA <- DA %>%
  rename(
    sex = Sex,
    beta = `Beta60 (g/kg/h)`,
    age = `Age (years)`,
    Co = `Co (g/Kg)`,
    weight = `Weight (kg)`,
    height = `Height (cm)`,
    AAC = `Amount of Alcohol Consumed (g)`,
    drinkingtime = `Drinking Time (h)`,
    maxBAC = `Maximum BAC (g/Kg)`,
    BACpeaktime = `BAC Peak time (min)`
  )








