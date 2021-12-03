
<!-- PROJECT LOGO -->
<br />
<p align="center">

  <h1 align="center"><a href="https://bnicolet.com/publications/Nicolet2021Large.html">Zebrafish Traction Force Microscopy</a></h1>

<!--   <a href="https://bnicolet.com/publications/Nicolet2021Large.html">
    <img src="https://bnicolet.com/publications/images/Nicolet2021Large-teaser.jpg" alt="Logo" width="100%">
  </a> -->

  <p align="center">
    Publish info, date.
  </p>

  <p align="center">
    <a href='https://bnicolet.com/publications/Nicolet2021Large.pdf'>
      <img src='https://img.shields.io/badge/Paper-PDF-red?style=flat-square' alt='Paper PDF'>
    </a>
  </p>
</p>

<br />
<br />



<br />
<br />

This repository contains a sample implementation of our algorithm to compute traction stresses.

### Installing

The compilation is performed under Clang-12.
```bash
mkdir build
cd build
cmake ..
make -j
```

## Dependencies



## Repository structure

The `largesteps` folder contains the parameterization module made available via
`pip`. It contains:
- `geometry.py`: contains the laplacian matrix computation.
- `optimize.py`: contains the `AdamUniform` optimizer implementation
- `parameterize.py`: contains the actual parameterization code, implemented as a
  `to_differential` and `from_differential` function.
- `solvers.py`: contains the Cholesky and conjugate gradients solvers used to
  convert parameterized coordinates back to vertex coordinates.

Other functions used for the experiments are included in the `scripts` folder:
- `blender_render.py`: utility script to render meshes inside blender
- `constants.py`: contains paths to different useful folders (scenes, remesher, etc.)
- `geometry.py`: utility geometry functions (normals computation, edge length, etc.)
- `io_ply.py`: PLY mesh file loading
- `load_xml.py`: XML scene file loading
- `main.py`: contains the main optimization function
- `preamble.py`: utility scipt to a import redundant modules for the figures
- `render.py`: contains the rendering logic, using `nvdiffrast`

## License



## Citation
