## Treadmill-walking-and-running-kinematic-electromyographic-and-muscle-synergies-differences

### Description
The purpose of the project, carried out during the `Neuroengineering` university course, is three-fold: it’s proposed the analysis and the comment of the results obtained comparing both treadmill walking (1.4 m/s) and running (2.8 m/s) tasks, carried out by a group of 30 young adults, kinematically, ectromyographically and focusing on the muscle synergy recruitment. 

In particular, the parameters taken into account among the investigation are: 
- stance
- swing and cadence (kinematic)
- time-normalized to gait-cycle length sEMG envelopes (electromyography) and motor modules and primitives (muscle synergies). 

The entire work, from data acquisition to data processing and evaluation, has been done in `Matlab`.

![SYNERGIES_RUNNING](https://user-images.githubusercontent.com/93050652/196741651-502b0b39-d866-4d03-9757-8edf64c18f1a.png)

![SYNERGIES_WALKING](https://user-images.githubusercontent.com/93050652/196741686-53090e53-3b10-47a0-9243-b43ea23eb81a.png)

### How to run it

You need to run just the main script, but first you have to download the library 'myMuscleSynergiesLibrary' and the starting data 'CYCLE_TIMES' and 'RAW_EMG'. Unfortunately, due to its dimentions, the file 'RAW_EMG' cannot be uploaded in GitHub, but I'm willing to share it, if you're interested.

The files 'MuscleSynergies_TW_alen' and 'MuscleSynergies_TR_alen' are also available in this repository: they contain the muscle synergies extracted for both the tasks. Using them, instead of running the extraction part of the code, you will save a lot of time and computational work.

Do you wanna know more? Please have a look at the technical [report](https://github.com/valeriopu/Treadmill-walking-and-running-kinematic-electromyographic-and-muscle-synergies-differences/blob/main/files/report.pdf).
