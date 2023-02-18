# MPC4AMV toolbox
MPC4AMVs is a comprehensive open-source toolbox, which applies advanced MPC algorithms for autonomous marine vehicles (AMVs), which are developed in the book: "Advanced Model Predictive Control for Autonomous Marine Vehicels" by Yang Shi, Chao Shen, Henglai Wei, and Kunwu Zhang. 

MPC4AMVs enables the efficient formulation and solution of control problems for nonlinear AMVs, including the dynamic positioning, path following, trajectory tracking. The m-code is primarily developed on a Linux machine using MATLAB 2020a. The code should work on any platform such as Windows and MacOS, but is developed and thus most extensively tested on Linux. In order to achieve better code reusability, and to make the code structurally clear and easier to be understood, the simulation examples are mainly coded following the Object-Oriented Programming (OOP) paradigm. Therefore, common classes and library functions (such as AUV, point2D, MPC_controller, RefGen, etc.) can be shared between chapters. The OOP coding style makes the code modular, and the readers can test and compare performances by simply replacing the controller module with different controllers. Meanwhile, to make the simulation codes extendable, several abstract classes (such as Model, Controller, Trajectory2D, etc.) are created.

The development is continued at the Applied control and information processing laboratory (ACIPL) of the University of Victoria.


# Citing MPC4AMV 
Y. Shi, C. Shen, H. Wei, and K. Zhang, Advanced Model Predictive Control for Autonomous Marine Vehicles, Springer Nature, 2023, In Press, ISBN 978-3-031-19353-8.
