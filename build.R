library(knitr)
library(markdown)

knit("download_combine_and_visualize.Rmd")
markdownToHTML("download_combine_and_visualize.md", "download_combine_and_visualize.html")
