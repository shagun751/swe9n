# Develpoment log for Shallow Water Equation (GWCE) model

## Version swe9n

**Details**  
continued from swe9n\_v5.34ea.sqB.tA

- Quadrilateral 9-noded element
- Simpsons's 3x3 points integration
- Lumped mass-matrix => No matrix inversion
- Predictor corrector time-stepping
- OpenMP parallelisation
  
1. [GWCE Model Notes](./gwce_model_notes.md)
1. [Initial development log](./log_swe9n_v0001.md)

-----------------------------------------------

## References
1. Kolar, R., & Westerink, J. (2000). A look back at 20 years of GWC-based shallow water models. … on Computational Methods in Water …, 2(2). Retrieved from http://www.coe.ou.edu/emgis/kolar/resources/CMWRXIIIPaper1.pdf

1. Gray, W. G., & Lynch, D. R. (1977). Time-stepping schemes for finite element tidal model computations. Advances in Water Resources, 1(2), 83–95. https://doi.org/10.1016/0309-1708(77)90026-4

1. Lynch, D. R., & Gray, W. G. (1979). A wave equation model for finite element tidal computations. Computers & Fluids, 7(3), 207–228. https://doi.org/10.1016/0045-7930(79)90037-9

1. Deltares. (2014). WES - Delft3D. Retrieved from https://content.oss.deltares.nl/delft3d/manuals/Delft3D-WES_User_Manual.pdf

1. Luettich, R., & Westerink, J. (2004). Formulation and Numerical Implementation of the 2D/3D ADCIRC Finite Element Model Version 44.XX.

1. Holland, G. J. (1980). An Analytic Model of the Wind and Pressure Profiles in Hurricanes. Monthly Weather Review, 108(8), 1212–1218. https://doi.org/10.1175/1520-0493(1980)108<1212:AAMOTW>2.0.CO;2

1. Vatvani, D., Gerritsen, H., Stelling, G. S., & Rao, A. V. R. K. (2002). Cyclone Induced Storm Surge and Flood Forecasting System for India. Solutions to Coastal Disasters ’02, 473–487. https://doi.org/10.1061/40605(258)42

1. Dresback, K. M., & Kolar, R. L. (2002). An implicit time‐marching algorithm for shallow water models based on the generalized wave continuity equation. International Journal for Numerical Methods in Fluids, 36(8), 925–945. https://doi.org/10.1002/fld.157.abs

1. Dresback, K. M., Kolar, R. L., & Dietrich, J. C. (2004). A 2D implicit time-marching algorithm for shallow water models based on the generalized wave continuity equation. International Journal for Numerical Methods in Fluids, 45(3), 253–274. https://doi.org/10.1002/fld.697

1. Carter, G. S., & Merrifield, M. A. (2007). Open boundary conditions for regional tidal simulations. Ocean Modelling, 18(3–4), 194–209. https://doi.org/10.1016/j.ocemod.2007.04.003

1. Wind Drag Based on Storm Sectors. [Link](https://ccht.ccee.ncsu.edu/wind-drag-based-on-storm-sectors/)
