# ToMCCA v1.0

Thanks for chosing ToMCCA for you coalescence research. ToMCCA is a small project made to facilitate Coalescence studies using a lightweight, yet powerful framework. For further information on the coalescence model used please read https://arxiv.org/abs/2302.12696 and https://arxiv.org/abs/2404.03352.pdf. Please cite the latter when you use ToMCCA in your studies.

  

For now this version is very simple, but it produces stunning results already. We will work hard in the future to improve the versatility and usability of ToMCCA as well as predefine specific studies for ease of use. If you have any suggestions for features as well as bug reports, please do not hesitate to contact me under maximilian.horst@tum.de

  

Thank you again and enjoy ToMCCA!

### Dependencies 
ToMCCA requires the **ROOT framework** (https://root.cern/install/) as well as the JSON interpreter **jq** which is preinstalled in many Linux distributions (and via brew on Mac).

## Manual

Running ToMCCA is simple! just follow these steps:

### Step 1 JSONs

ToMCCA is a versatile program that we need to configure. For this we use JSON files which allow good readability without clunky interfaces. This basic version only delivers one example JSON calles Test.json and so far there are only 5 different settings:

**Mode**: Currently there is only one mode available, called "Default". This setting can be used to run different global modes, such as specific studies, hyperthreading and more.

**NEvents**: Defines the number of events you want to simulate.

**Mult**:  Sets the mean multiplicity of the events you want to study. Note that this is given in dN/deta while experimental data often is cut on rapidity, but a transformation is built into ToMCCA. Also note that the value is **multiplied** by 10 to assure one decimal precision. 

**MultType**: Selects the type of multiplicity distribution you want for your study. For now we have 4 options: Poissonian (1), Fixed(2), Erlang(3) and ReadFromFile.  A "**Poissonian**" distribution will draw the number of particles from a Poissonian distribution with the mean=*Mult* and Variance=*Mult*. A "**Fixed**" distribution will give every event you simulate exactly the number of charged particles you set it to in the *Mult* variable (*be aware of roundoff errors!*).  The third setting "**Erlang**" uses a more realistic parameterization based on Erlang functions fitted to multiplicity distributions detemined using EPOS 3. The fourth setting "**ReadFromFile**" reads in a file that has to be provided in ln49 of ToMCCA_main.cpp.

  

### Step 2 Running

To run ToMCCA, simply execute ToMCCA.sh followed by the name of the JSON file you want to use and a short comment of the settings of this run to be added to the 
**RunPlan.txt**! This will start the process, as well as create an output folder with the current date. If you have multiple outputs from one day, it will enumerate them.
**Example**: bash ToMCCA.sh Test "Testing ToMCCA for the first time!"
Will execute ToMCCA with the setting defined in JSONs/Test.json and will write your start time and message into the RunPlan.txt

  

### Step 3 Output

The output of ToMCCA will be saved to the output folder in a file called *Output.root*. Here you can find a multitude of output histograms such as proton/deuteron spectra, multiplicity distribution or relative Momentum *q* vs distance *r* distribution of nucleons. Feel free to add your own and also to suggest further standard output histograms!
