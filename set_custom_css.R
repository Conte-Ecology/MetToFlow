options(rstudio.markdownToHTML = 
             function(inputFile, outputFile) {      
                  require(markdown)
                  markdownToHTML(inputFile, outputFile, stylesheet='C:/Program Files/RStudio/resources/custom.css')   
             }
)