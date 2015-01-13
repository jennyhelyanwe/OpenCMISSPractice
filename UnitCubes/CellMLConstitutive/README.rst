.. _OpenCMISSPractice-UnitCubes-CellMLConstitutive:

=====================================
Constitutive Models Defined in CellML
=====================================

OpenCMISS has the capability to link with CellML mark-up descriptions for the 
constitutive models. This is done by mapping constitutive parameters and strain fields
as known fields to the CellML model, and mapping the second Piola-Kirchhoff stress
fields as unknowns to be mapped from the CellML model. 

A flow diagram illustrates how this linking works...

See below for a list of specific examples of defining constitutive models in CellML:

.. toctree::
   : maxdepth: 2

   mooney_rivlin
    
