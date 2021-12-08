# Theoretical_Coastal_Model
A theoretical model of the upwelling system on western coast of the United States. The model incorporates a Cell Flux Model of Phytoplankton into a basic physical model to illustrate the nutrient distribution in phytoplankton and the surrounding surface waters.


There are two different models available in this repository: Upwelling_Project_noCFM and Upwelling_Project_with-CFM. Both models contain representations for upwelling of nutrients using advection equations. The output provides nutrient concentrations of nitrate and phosphate in the surface waters and phytoplankton. The model resolves nutrient concentrations for the temporal scale of one year and spatial scale of 400 km moving away from the west coast. 

Upwelling_Project_noCFM provides a theoretical nutrient distribution in phytoplankton and surface waters without the use of the Cell Flux Model. This provides a very basic demonstration of theoretical interactions between phytoplankton and surrounding surface waters. 

Upwelling_Project_with-CFM provides a theoretical nutrient distribution in phytoplankton and surface waters with the use of the Cell Flux Model. This provides a more realistic view of the interaction between phytoplankton and nutrient availability in the surrounding surface waters. 

Parameters to be changed by user if desired:

nyr: This is the number of years to simulate. 

Nutrients in surface waters: These can be changed if the user has data to supply this information 

  Pnut: Phosphate concentration in surface waters
  
  Pnut_i:Inital surface water phosphate concentration at the coastal boundary
  
  Nnut: Nitrate concentration in surface waters
  
  Nnut_i: Inital surface water nitrate concentration at the coastal boundary
  
Kn: half-saturation constant for Nnut- This can be changed if the phytoplankton community is dominated by a specific species and the user wants to use this for that species. 

Kp:half-saturation constant for Pnut- This can be changed if the phytoplankton community is dominated by a specific species and the user wants to use this for that species. 

U0: wind forcing velocity in meters/second- If better data is available, then can be substituted in. 
