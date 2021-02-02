
## Setup reticulate using vitrualenv
#https://support.rstudio.com/hc/en-us/articles/360023654474-installing-and-configuring-python-with-rstudio
#install.packages("reticulate")
#reticulate::py_config() ## to check correct local version is being used

## <------- Troubleshooting -------> ##

#Note: .Rprofile must have a trailing newline ('\n') it it will not be read. It should be fine for now,
#but if the python version defaults to your local machine (you can check with reticulate::py_config()),
#then make sure this is the case, then restart R.

#Shiny library setup fix:
#install.packages("Rcpp", dependencies = TRUE)

# Publishing to shinyapps.io throws an error (Error in value[[3L]](cond) :
# Error 126 occurred running /srv/connect/apps/cr_tool/python/bin/python
# Calls: local ... tryCatch -> tryCatchList -> tryCatchOne -> <Anonymous>),
#to fix:
# https://community.rstudio.com/t/problem-deploying-app-using-a-virtual-env-with-reticulate-to-run-python-code-in-app-error-virtual-environment-permission-denied/25283/15
# or
# https://github.com/rstudio/reticulate/issues/757

## <------- End Troubleshooting -------> ##

## Setup virtualenv ##
# virtualenv_create("r-reticulate")
# virtualenv_install("r-reticulate",packages = c('pandas','numpy'))
# use_virtualenv("r-reticulate", required = TRUE)

## Python Files ##
#source_python('fileParser.py')
#fileParserPy <- import("fileParser")

## set options ##

