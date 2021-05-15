# Yuhaniz Aly

### Tuberculosis ####
# Load libraries
library(HSAUR)
library(MASS)
library(Epi)
library(epitools)
library(dplyr)
library(ggplot2)
library(plotly)

# load tb dataset
tuberculosis_data <- HSAUR::BCG

# Make new columns and calculate OR for each study
tb_data <- tuberculosis_data %>%
  mutate(odds_exposed = BCGTB / NoVaccTB) %>%
  mutate(odds_unexposed = (BCGVacc - BCGTB) / (NoVacc - NoVaccTB)) %>%
  mutate(odds_ratio = odds_exposed / odds_unexposed)

# edited tb dataframe for visual
edit_tb_data <- tb_data %>%
  select(Study, odds_ratio)

# Relationship between odds ratio and year
# Scatter plot for comparison
or_year_plot <- ggplot(tb_data) +
  geom_point(mapping = aes(x = Year, y = odds_ratio, color = Study)) +
  theme(legend.position = "none") +
  labs(
    title = "Relationship between odds ratio and year",
    x = "Year",
    y = "Odds Ratio"
  )

or_year_plot_inter <- ggplotly(or_year_plot)


# Relationship between odds ratio and latitude
or_latitude_plot <- ggplot(tb_data) +
  geom_point(mapping = aes(x = Latitude, y = odds_ratio, color = Study)) +
  theme(legend.position = "none") +
  labs(
    title = "Relationship between odds ratio and latitude",
    x = "Latitude",
    y = "Odds Ratio"
  )

or_latitude_plot_inter <- ggplotly(or_latitude_plot)


# Calculate Relative Risks
tb_data_edit <- tb_data %>%
  mutate(rr_exposed = (BCGTB / BCGVacc)) %>%
  mutate(rr_unexposed = (NoVaccTB / NoVacc)) %>%
  mutate(relative_risk = rr_exposed / rr_unexposed)


# Plot for tb OR vs RR
or_vs_rr_plot <- ggplot(tb_data_edit) +
  geom_point(mapping = aes(x = odds_ratio, y = relative_risk, color = Study)) +
  theme(legend.position = "none") +
  labs(
    title = "Odds Ratio VS Relative Risk",
    x = "Odds Ratio",
    y = "Relative Risk"
  ) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(0, 1.6) +
  ylim(0, 1.6)


or_vs_rr_plot_inter <- ggplotly(or_vs_rr_plot)

### Low Birth Weight ###

# load birthweight dataset
lowbirth_df <- MASS::birthwt

# number of rows
num_rows <- nrow(lowbirth_df)

# proportion of low birth weight
proportion_lbw <- sum(as.numeric(lowbirth_df$low)) / num_rows

# proportion of smoker status
proportion_smokers <- sum(as.numeric(lowbirth_df$smoke)) / num_rows

# mean birth weight
mean_bw <- mean(lowbirth_df$bwt, na.rm = TRUE)

# keys
# - Smoking status = smoke
# - Hypertension status = ht
# - Attending 0 prenatal care visits = ftv --> renamed: zero_prenatal
# - Giving birth before age 20  = age --> renamed: under20

# created new columns in data frame
edit_low_df <- lowbirth_df %>%
  select(low, age, smoke, ht, ftv) %>%
  mutate(under20 = if_else(age < 20, 1, 0)) %>%
  mutate(zero_prenatal = if_else(ftv == 0, 1, 0))

# function for relative risk calculations
get_rr <- function(df, col) {
  a <- sum(df$low == 1 & df[col] == 1)
  b <- sum(df$low == 0 & df[col] == 1)
  c <- sum(df$low == 1 & df[col] == 0)
  d <- sum(df$low == 0 & df[col] == 0)

  rr <- (a / (a + b)) / (c / (c + d))
  return(rr)
}

rr_smoke <- get_rr(edit_low_df, "smoke")
rr_zero_prenatal <- get_rr(edit_low_df, "zero_prenatal")
rr_ht <- get_rr(edit_low_df, "ht")
rr_under20 <- get_rr(edit_low_df, "under20")

exposures <- c(
  "Smoking status", "Attending 0 prenatal care visits",
  "Hypertension status", "Giving birth before age 20"
)

rr_exposure <- c(rr_smoke, rr_zero_prenatal, rr_ht, rr_under20)

# new data frame comparing exposure of interests and relative risk
edit_rr_table <- data.frame(exposures, rr_exposure, stringsAsFactors = F)

### Endometrial Cancer ####

# load endometrial cancer dataset
data(bdendo)

endo_data <- bdendo %>%
  select(d, gall, hyp, ob)

# calculate presence of Gall bladder disease
a_gall <- sum(endo_data$d == 1 & endo_data$gall == "Yes")
b_gall <- sum(endo_data$d == 0 & endo_data$gall == "Yes")
c_gall <- sum(endo_data$d == 1 & endo_data$gall == "No")
d_gall <- sum(endo_data$d == 0 & endo_data$gall == "No")

# use epitool to find or
table_gall <- matrix(c(a_gall, b_gall, c_gall, d_gall), ncol = 2, byrow = TRUE)
colnames(table_gall) <- c("Disease", "No Disease ")
rownames(table_gall) <- c("Yes Gall", "No Gall")
table_gall <- as.table(table_gall)


table_end <- matrix(c(d_gall, c_gall, b_gall, a_gall), 2, 2)
dimnames(table_gall) <- list(
  Exposure = c("No", "Yes"),
  "Endometrial Cancer" = c("No Disease", "Disease")
)
table_gall


disease_or_gall <- epitab(table_gall, method = c("oddsratio"), verbose = F)

# pulling odds ratio for hypertension
gall_or <- disease_or_gall$tab[2, 5]


## calculate for hypertension
a_hyp <- sum(endo_data$d == 1 & endo_data$hyp == "Yes")
b_hyp <- sum(endo_data$d == 0 & endo_data$hyp == "Yes")
c_hyp <- sum(endo_data$d == 1 & endo_data$hyp == "No")
d_hyp <- sum(endo_data$d == 0 & endo_data$hyp == "No")

# use epitool to find or
table_hyp <- matrix(c(a_hyp, b_hyp, c_hyp, d_hyp), ncol = 2, byrow = TRUE)
colnames(table_hyp) <- c("Disease", "No Disease ")
rownames(table_hyp) <- c("Yes Hyp", "No Hyp")
table_hyp <- as.table(table_hyp)

table_hyp <- matrix(c(d_hyp, c_hyp, b_hyp, a_hyp), 2, 2)
dimnames(table_hyp) <- list(
  Exposure = c("No", "Yes"),
  "Endometrial Cancer" = c("No Disease", "Disease")
)
table_hyp


disease_or_hyp <- epitab(table_hyp, method = c("oddsratio"), verbose = F)

# pulling odds ratio for hypertension
hyp_or <- disease_or_hyp$tab[2, 5]


## get rid of na obesity values
edit_endo_data <- bdendo %>%
  select(d, gall, hyp, ob) %>%
  filter(!is.na(ob))

## calculate for obesity
a_ob <- sum(edit_endo_data$d == 1 & edit_endo_data$ob == "Yes")
b_ob <- sum(edit_endo_data$d == 0 & edit_endo_data$ob == "Yes")
c_ob <- sum(edit_endo_data$d == 1 & edit_endo_data$ob == "No")
d_ob <- sum(edit_endo_data$d == 0 & edit_endo_data$ob == "No")

table_ob <- matrix(c(a_ob, b_ob, c_ob, d_ob), ncol = 2, byrow = TRUE)
colnames(table_ob) <- c("Disease", "No Disease ")
rownames(table_ob) <- c("Yes Hyp", "No Hyp")
table_ob <- as.table(table_ob)

table_ob <- matrix(c(d_ob, c_ob, b_ob, a_ob), 2, 2)
dimnames(table_ob) <- list(
  Exposure = c("No", "Yes"),
  "Endometrial Cancer" = c("No Disease", "Disease")
)
table_ob

disease_or_ob <- epitab(table_ob, method = c("oddsratio"), verbose = F)

# pulling odds ratio for obesity
ob_or <- disease_or_ob$tab[2, 5]
