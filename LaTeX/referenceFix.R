## This script lists issues that have been pointed out and corrects them

## Currently making corrections manually

# First copy the auto-updated MHW bib file from my computer to the MHW folder
   # THis prevents my Mendeley from doing this automatically
if(file.exists("/home/rws/Documents/Project/References/MHW.bib")){
  file.copy(from = "/home/rws/Documents/Project/References/MHW.bib", 
            to = "/home/rws/MHW/LaTeX/MHW.bib", overwrite = TRUE)
}

# Then open the .bib file in R for the manual touch ups listed below
file.edit("LaTeX/MHW.bib")

## Things to correct:

# in situ -> \emph{in situ}

# Perna viridis -> \emph{Perna viridis}

# Loligo vulgaris reynaudii -> \emph{Loligo vulgaris reynaudii}

# Alosa pseudoharengus -> \emph{Alosa pseudoharengus}

# Ecklonia maxima -> \emph{Ecklonia maxima}

## When preparing the .doc file two shortcuts are taken 
  # This is a reminder how to undo them

# \cref -> \Cref

# Insert contents of MHW.bib accordingly