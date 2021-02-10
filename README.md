# Melts-Vols

moorelr/Melts-Vols is licensed under The MIT License
Copyright 2019 Lowell R. Moore

For questions or comments, email moorelr@vt.edu

This program uses a MELTS output file to draw a cumulative mineral volume figure. To use, copy the melts.out file into the folder, and run the draw plots.R script. This can be done either with the included .bat file (as long as rscript.exe is in the system path variable) or by running the .R script in RStudio or terminal.

This program was written for use in Windows, but will probably work for Linux and Mac users also (not tested!)

Workflow:
- Paste melts.out file into the working directory (same folder as the .R file)
- Update the "settings" file and choose the desired parameters for the figure.
- Run the script in the working directory.

- Note: the script will generate a "key" figure, which is basically redundant with the text labels. The colors of the text in the "key" correspond to the shaded fields in the main figure.
