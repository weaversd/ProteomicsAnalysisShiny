source("global.R")
source("R/parsers.R")
source("ui.R")
source("server.R")

shinyApp(ui = ui, server = server)