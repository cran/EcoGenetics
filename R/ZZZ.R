
.onAttach <- function(...) {
	
	vers <- utils::packageDescription("EcoGenetics", fields = "Version")
  
	textstart<- paste("\n\n",
										"\n", "                #####################", "\n",
										"\r                   ","|","EcoGenetics", "|",
										"\n", "                #####################", "\n",
										"\n", "  Analysis of phenotypic, genotypic and environmental data",
										"\n",
										"\n", "  Version",  vers, "\t\t",
										" <leandroroser@ege.fcen.uba.ar>", "\n\n")
	
	packageStartupMessage(textstart)
}
