cat("getwd():", getwd(), "\n")
cat("R home:", R.home(), "\n")
cat(".libPaths():\n")
print(.libPaths())
cat("search path:\n")
print(search())
cat("renv active:", "renv" %in% loadedNamespaces(), "\n")

library(emmeans)