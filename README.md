# Fast/Fortran/Fantasy Map Generator (fmap)

`FMAP` is an experimental tool for generating fantasy world maps fast. The goal is to write a relatively simple tool that produces physically believable maps without the need for full-fledged Earth system models. The approach is to mathematically imitate nature where possible, and impose only a few, heavily simplified physics-informed rules, and see how believable I can make it.

The worlds are built by successively adding layers of complexity. First, tectonic plates are imitated through voronoi cells. These are tagged as oceanic or continental. Mountain ranges are added where expected. Topographic noise is added, and landscapes are modified by a simple, imposed climate.


![Voronoi cells](voronoi.png)

*Tectonic plates are imitated by voronoi cells. You can choose between Manhattan over Euclidean distance and a flat, cylindrical, torus or spherical world shape to draw the cells. The continents at the top are generated using the the torus shape (seamless texturing) and Manhattan distance. The world at the bottom is spherical, and the cells are generated using the Euclidean distance.*

![Land sea mask](landsea.png)
*You can set the % of the surface to be covered by oceans. This will make the world generator assign different densities to tectonic plates. Dark grey represent oceans, while lighter grey are land.*
