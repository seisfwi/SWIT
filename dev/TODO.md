# TODOs for SWIT (release v1.1.0)

## Main tasks
- Add topography to the forward modeling （FAILED）
    - test the forward modeling with topography and its adjoint
- Add options of receiver components (Solved)
    - vx, vz, p
    - test the adjoint modeling with different receiver components
- Add Multiscale workflow
    - test the multiscale workflow
- Add source inversion workflow
    - test the source inversion workflow
- Add geometry retrieval code from su or segy files
- Add MPI support for the preprocessing and misfit calculation
- Add RTM workflow
    - test the RTM workflow
- 


## Other tasks
- Add tests for all operators
    - forward modeling benchmarks (with and without topography)
    - dot-product tests
- Add documentation
- Add examples, including
    - synthetic examples
    - field data examples
    - scripts for checking the results
- Add data processing code (MATLAB)
- Automate the installation process


## Long-term tasks
- Add GPU support by using Devito
- Add a GUI for the workflow