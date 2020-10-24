# pyElectrochem
Some programs to simulate common experiments in electrochemistry : chronoamperometry, voltammetry and maybe more in a near future. Implemented in Python, with numpy, scipy and matplotlib under the hood.

Some videos are available here : 
https://www.youtube.com/playlist?list=PL_rdrSitJev4SUsyvfQLp_aizC7w2t3qi

This was made for a purely educational purpose, use it at your own risk

People interested can also have a look at Soft Potato :

https://github.com/oliverrdz/SoftPotato (Created by Oliver Rodríguez)

## Voltammetry with a Nernstian boundary condition
Use [voltammetry.py](voltammetry-nernst/voltammetry.py), some help is given with "./voltammetry.py -h" about all the arguments and default values used. You can save a movie, a csv file, a numpy version for further use.

https://www.youtube.com/watch?v=QCm5n93nyhM&list=PL_rdrSitJev4SUsyvfQLp_aizC7w2t3qi&index=2

## Voltammetry for a quasi reversible boundary condition
Use "[voltammetry-quasi-reversible.py](voltammetry-quasi-reversible/voltammetry-quasi-reversible.py)" some help is given with "voltammetry-quasi-reversible.py -h" about all the arguments and default values used. 

https://www.youtube.com/watch?v=QCm5n93nyhM&list=PL_rdrSitJev4SUsyvfQLp_aizC7w2t3qi&index=3

## Chronoamperometry
Use "[chronoamperometry.py](chronoamperometry/chronoamperometry.py)". The boundary condition imposed is to set the conventration to zero at the electrode.

https://www.youtube.com/watch?v=a9UcofmM4ns&list=PL_rdrSitJev4SUsyvfQLp_aizC7w2t3qi


## Butler Volmer (steady-state amperometry)
Use "[Butler-Volmer.py](Butler-Volmer/Butler-Volmer.py)". Requires "[widgets.py](Butler-Volmer/widgets.py)". Prefer the use of "python3 Butler-Volmer.py" over "./Butler-Volmer.py" to exectute the script 

https://www.youtube.com/watch?v=QCm5n93nyhM&list=PL_rdrSitJev4SUsyvfQLp_aizC7w2t3qi&index=4

## Diffusion controlled (steady-state amperometry)
Use "[steady-state-diffusion.py](ss-diffusion/steady-state-diffusion.py)". Requires "[widgets.py](ss-diffusion/widgets.py)". Prefer the use of "python3 steady-state-diffusion.py" over "./steady-state-diffusion.py" to exectute the script 

https://www.youtube.com/watch?v=QCm5n93nyhM&list=PL_rdrSitJev4SUsyvfQLp_aizC7w2t3qi&index=5

## Mixed control (steady-state amperometry)
Use "[ss-mixed.py](ss-mixed/ss-mixed.py)". Requires "[widgets.py](ss-mixed/widgets.py)". Prefer the use of "python3 ss-mixed.py" over "./ss-mixed.py" to exectute the script 

https://www.youtube.com/watch?v=QCm5n93nyhM&list=PL_rdrSitJev4SUsyvfQLp_aizC7w2t3qi&index=6

 
