# ToMCCA v1.0
Thanks for chosing ToMCCA for you coalescence research. ToMCCA is a small project made to facilitate Coalescence studies using a lightweight, yet powerful framework. For further information on the coalescence model used please read https://arxiv.org/abs/2302.12696

For now this version is very simple, but it produces stunning results already. We will work hard in the future to improve the versatility and usability of ToMCCA as well as predefine specific studies for ease of use. If you have any suggestions for features as well as bug reports, please do not hesitate to contact me under maximilian.horst@tum.de

Thank you again and enjoy ToMCCA!
Maximilian Horst
## Manual
Running ToMCCA is simple! just follow these steps:
### Step 1 JSONs
ToMCCA is a versatile program that we need to configure. For this we use JSON files which allow good readability without clunky interfaces. This basic version only delivers one example JSON calles Test.json and so far there are only 3 different settings: **NEvents** which tells ToMCCA the number of events you want to simulate. 
Second is **Mult**, which tells ToMCCA the mean multiplicity of the events you want to study. Note that this is given in dN/deta while experimental data often is cut on rapidity, but a transformation is built into ToMCCA. 
Lastly, **MultType** selects the type of multiplicity distribution you want for your study. For now we have 2 options: Fixed (2) and Poissonian(1). A "fixed" distribution will give every event you simulate exactly the number of charged particles you set it to in the *Mult* variable (*be aware of roundoff errors!*). A "Poissonian" distribution will draw the number of particles from a Poissonian distribution with the mean=*Mult* and Variance=*Mult*.

### Step 2 Running
To run ToMCCA, simply execute ToMCCA.sh! This will start the process, as well as create an output folder with the current date. If you have multiple outputs from one day, it will enumerate them.

### Step 3 Output
The output of ToMCCA will be saved to the output folder in a file called *Output.root*. Here you can find a multitude of output histograms such as proton/deuteron spectra, multiplicity distribution or relative Momentum *q* vs distance *r* distribution of nucleons. Feel free to add your own and also to suggest further standard output histograms!


