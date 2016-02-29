# Ecotoner

We haven’t really published the code yet nor prepared it for sharing (though through our use of github made it openly accessible), it is actively and gradually being developed by the Schlaepfer lab, and there is no manual - we cannot give you individual support in setting up and running the code except if we agreed on a collaboration or similar agreement.

Not every part of the code has been extensively tested or is in a stable state. Similarly, not every combination of model inputs and options has been evaluated in depth and there is not guarantee for anything to work. The code comes with no warranty and no guarantees, expressed or implied, as to suitability, completeness, accuracy, and whatever other claim you would like to make.

There is no graphical user interface, help pages and available documentation may be out of date, and you will need to write your own tools to analyse outputs.

If you make use of this model, please cite appropriate references, and we would like to hear about your particular study (especially a copy of any published paper).


Some recent references
* 



### Install the package in R
```
devtools::install_git("https://github.com/dschlaep/ecotoner.git")
```

### Information about the package
```
package?ecotoner			# A few remarks
help(package = "ecotoner")	# The index page with some of the functions (which are documented to a varying degree)
```

### Using the package
I added the main code which I use to locate and measure my transects (producing the data for later analysis), as demo to the package. I also added data so that it can run  as a small contained demo ‘example’ (if the flag do.demo is set to TRUE). You could run this code directly with the demo() function, but this is probably not convenient:
```
demo("BSE-TF_EcotonesAtSlope_v10", package = "ecotoner")
```

Instead, the following command should open my main code in your text editor (at least on unix-alike systems) for easier inspection:
```
system2("open", file.path(system.file("demo", package = "ecotoner"), "BSE-TF_EcotonesAtSlope_v10.R"))
```


### For contributors only
### How to contribute
You can help us in different ways:

1. Reporting [issues](https://github.com/dschlaep/ecotoner/issues)
2. Contributing code and sending a [pull request](https://github.com/dschlaep/ecotoner/pulls)

In order to contribute to the code base of this project, you must first contact the Schlaepfer Lab. We retain any decision to accept your suggestions/contributions.
