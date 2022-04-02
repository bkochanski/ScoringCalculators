# from https://rstudio-pubs-static.s3.amazonaws.com/376828_032c59adbc984b0ab892ce0026370352.html#113_introduction_to_woe_and_iv_and_a_complementary_r_package 


# install.packages("scorecard")
library(scorecard)
library(dplyr)
data(germancredit)
#View(germancredit)
data = germancredit %>%
  as_tibble()

#replace '.' in variable names not compatible with f_train_lasso
vars = names(data) %>%
  stringr::str_replace_all( '\\.', '_')

names(data) <- vars

# convert response factor variable to dummy variable

data = data %>%
  mutate( creditability = ifelse( creditability == 'bad', 1, 0 )
          , creditability = as.factor(creditability) )

summary(data)

Amelia::missmap(data)

iv = iv(data, y = 'creditability') %>%
  as_tibble() %>%
  mutate( info_value = round(info_value, 3) ) %>%
  arrange( desc(info_value) )

iv %>%
  knitr::kable()

bins = woebin(data, y = 'creditability')

# View(bins$status_of_existing_checking_account)
# summary(data$status_of_existing_checking_account)

bins$duration_in_month %>%
  knitr::kable()
summary(data$duration_in_month)

woebin_plot(bins$duration_in_month)

data_woe = woebin_ply( data, bins ) %>%
  as_tibble()

set.seed(1)

vars = names(data_woe)
vars = vars[ vars != 'creditability']

formula = as.formula( paste( 'creditability ~', paste( vars , collapse = '+') ) )


# install.packages("devtools")
# devtools::install_github("mtennekes/tabplot")
# devtools::install_github("drsimonj/pipelearner")
# devtools::install_github("erblast/oetteR")

lasso = oetteR::f_train_lasso( data = data_woe
                               , p = NULL
                               , formula = formula
                               , k = 50
                               , family = 'binomial'
)

rlang::last_error()




