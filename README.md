This is the current research program I'm doing for doctor course. In this research, we divided the model into several groups and arrange each group local coordinates. Every group has its own rotation and solve its equations independently.
The formulation of this method is based on Linear Corotated FEM and Shape Matching for each group, and then binded together using a constraint force.
This method currently accelerates the Linear Corotated FEM (VegaFEM) by 7 times and keeps good accuracy. Our method can be implemented in various shapes.
For example, a Stanford Armadillo with right hand fixed and then dropped down by gravity.
![banner](https://github.com/MagnoliaKKK/tetfemcpp/assets/62364444/71263d49-657b-4ecb-ae39-48f4e41389c3)

We show better accuracy than an other corotated method, Operator Splitting, although they are quite fast.

<img src="https://github.com/MagnoliaKKK/tetfemcpp/assets/62364444/d7e5ac2a-74e7-4696-9f89-27745b069ecf" alt="VegaStretch" width="300"/>

<img src="https://github.com/MagnoliaKKK/tetfemcpp/assets/62364444/03cf29e7-0547-49bf-8637-0e848ff6f197" alt="LocalStretch" width="300"/>

<img src="https://github.com/MagnoliaKKK/tetfemcpp/assets/62364444/14a8c8bb-f665-46d3-b879-c06e3609b5c6" alt="OPStretch" width="300"/>



The image compares dragging a beam with a fixed left side. The beige one is VegaFEM, the blue one is our method and the grey blue is Operator Splitting.

Our method can also simulate anisotropic materials by re-formulating the stiffness matrix.

For more details, you can refer to:
[vrsj_localFEM.pdf](https://github.com/user-attachments/files/15973518/vrsj_localFEM.pdf)

Thanks for seeing my research!
