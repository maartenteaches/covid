# covid: Agent based model for the spread of a disease

## Description

This Agent Based Model for the spread of a disease. I have COVID-19
in mind. In its current form it extends a basic SIR model to include
a social network. 

This is a very crude model and I am not an epidemiologist, so only 
use it to play around with for your own education. This is not meant 
for making life or death decisions.


## Requirements and use

This requies [Stata](https://www.stata.com) version 15 and the 
[amb_nw](https://github.com/maartenteaches/abm_nw) class. `covid.mata`
defines the model, `nw.do` and `curves.do` run this model to create
the results shown in [this video](https://youtu.be/2KS-I74xUOM), and 
`cscript.do` the certification script that documents how I have tested
`covid.mata`.