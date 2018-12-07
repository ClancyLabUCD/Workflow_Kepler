# Workflow_Kepler

1) Coding Language for the cardiac models: C++

2) The post processing is done using: Python

3) Workflow is constructed in Kepler Workflow System (www.kepler-project.org)
The workflow and user manual is available in this repository.

3) How to convert your cardiac cell mode code for use with our system:

To facilitate execution of new cardian cell models, the path to C++ source code file is parametrized
in the Kepler Workflow. If user wants to use their own cardiac cell model, they can edit the Kepler 
Workflow parameter called 'sourceCode' under the category of 'SharedParameters', to point to the directory where the desired C++ source code resides.
