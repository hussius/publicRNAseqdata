library(knitr)
library(markdown)
knit("combine_and_visualize.Rmd")
markdownToHTML("combine_and_visualize.md", "combine_and_visualize.html")

knit("download.Rmd")
markdownToHTML("download.md", "download.html")
