# visual-stimulus-controller
This code sends commands over udp to a computer that presents visual stimuli (see visual-stimulus-generator repo). It also aligns experimental information (e.g. file names) with local widefield, local ephys, or remote 2photon acquisition. Visual stimulus parameters are configured using MATLAB GUI in this folder.  It is used on a Windows 10 machine with MATLAB 2021.

Install:

1) Download all files to the same root folder.
2) Set your Matlab path to this root folder. Do not set paths to deeper folders.
3) In Stimulator.m, change 'root_controller' to be the same as this root folder.
4) If you are using widefield analysis code, change root_WFanalysis to the location of that code. 
5) Create a folder that will store the .analyzer files, which contain all the experiment information files.  We call this folder /AnalyzerFiles/.  

configureMstate.m: 

Here, you will set folder locations, ip addresses, some flags, and a few other high level defaults. For example, 'Mstate.analyzerRoot' is the folder you created for the .analyzer files above.  ip addresses are for the visual stimulus computer (see visual-stimulus-generator repo), along with the 2p computer if you happen to be doing that (Ringach's Scanbox in our case).


Schematic of connectivity:

To perform widefield imaging alone, only computer 'A' and 'B' are used.  For two-photon acquisition, a third machine ('C') is necessary.

![Slide1](https://user-images.githubusercontent.com/13107530/145627782-aee2ca0b-4889-453e-992b-73eba64d8565.jpeg)
