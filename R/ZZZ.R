#' @keywords internal


.onAttach <- function(...) {
	
	vers <- utils::packageDescription("EcoGenetics", fields = "Version")
  
	textstart<- paste("\n\n",
										"\n", "                ---------------------", "\n",
										"\r                  ","||","EcoGenetics","||",
										"\n", "                ---------------------", "\n",
										"\n", "  Spatial Analysis of Phenotypic, Genotypic and Environmental Data",
										"\n",
										"\n", "  Version",  vers, "\n\n",
										"  GitHub: https://github.com/leandroroser/EcoGenetics-devel", "\n\n",
										"  Online tutorial: https://leandroroser.github.io/EcoGenetics-Tutorial", "\n\n",
										"  Overview: help('EcoGenetics')")
	
	packageStartupMessage(textstart)
}
