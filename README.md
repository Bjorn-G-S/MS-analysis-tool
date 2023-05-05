# MS-analysis-tool

* [General](#general-info)
* [Purpose](#purpose)
* [Instalation](#installation)
* [How-to](#how-to)
* [Contact](#Contact)
* [License](#License)


## General

`MS_Analysiss` is a python script  used for analyzing and quantifying results experiments done in the Methane To Methanol reaction (MTM). It takes raw ASCII-type files from a MID MassSpectrometer, integrates the peaks, and calculates the yield, selectivity and productivity of the reaction based on information of the catalyst. The return data is given back as a Excell table and saved on desired file path. Currently `MS_Analysiss` is used at the section for catalysis at the University of Oslo.

## Purpose

Analyze and visualize data from MS measurments taken during MTM experiments Cu-zeolites.


## How-to
How to use the program.


1. Import the program:
```
From MS-class import MS_analysis_MTM
```
2. Define the directory of the files that are to be converted and the file with the x-vlaues:
```
obj = MS_analysis_MTM()
```
3. Run the program. the following message will the apear:
```
                   THIS CLASS IS USED FOR ANALYZING EXPERIMENTS (MS AND TEMPERATURE DATA) TAKEN AT THE M2M-rig 
                                              use the following commands in the MS_analysis:



        obj.help()                    -        Will give this list over commands, and what their function is.



        
                                               Methane to Methanol experiments:  



        obj.plot_MS()                 -        Old script for analyzing MS data from a MTM catalytic test
        
        obj.plot_MS_Area()            -        Defining the Area for the MS masses

        obj.Area_to_csv()             -        Used to convert the area peaks to .csv so that it can by analyzed outside the script. 
        
        obj.testing_resutls()         -        When Area is defined form obj.plot_MS_Area, this will do the calcualtion to get the yield, selectivity etc.
```
**** Please note that the file paths need to be changed to the desired filepaths before using the program **** 


## Dependencies

This script requires a python enviornment with the following packages, and the packaged these depend on:
```
python          (3.9.7)
pandas          (1.3.3)
numpy           (1.21.2)
matplotlib      (3.4.3)
ipywidgets      (7.6.5)
scipy           (1.7.1)
typing          (3.10.0.2)
ipympl          (0.8.5)
ipython         (7.27.0)
```

## Contact

For developer issues, please start a ticket in Github. You can also write to the dev team directly at  **b.g.solemsli@smn.uio.no**
#### Authors: 
Bjørn Gading Solemsli (@bjorngso).

## License
This script is under the MIT license scheme. 




