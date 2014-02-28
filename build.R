library(knitr)
library(markdown)

knit("compare_published_RPKMs.Rmd")
markdownToHTML("compare_published_RPKMs.md", "compare_published_RPKMs.html")
