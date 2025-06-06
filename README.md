# Monte Carlo simulations of diffusion in a medium composed of randomly packed spheres of permeable membrane (CUDA C++)

The code implements 3d Monte Carlo simulations originally developed in [Lee, et al., Journal of Neuroscience Methods, 2021](https://doi.org/10.1016/j.jneumeth.2020.109018), demonstrating the effect of water exchange and diffusion on the time-dependent kurtosis in [Lee et al., Magnetic Resonance in Medicine 2025](https://doi.org/10.1002/mrm.30335).

* **Demo 1, randomly packed spheres, narrow pulse sequence:** We simulate water diffusion in a medium composed of randomly packed permeable spheres, identify the time to the kurtosis peak, and compare it with the theoretical prediction. We created multiple geometries with varying sphere diameter, sphere volume fraction, and membrane permeability.
* **Demo 2, randomly packed spheres, narrow pulse sequence:** We simulate water diffusion in a medium composed of randomly packed permeable spheres, identify the time to the kurtosis peak, and compare it with the theoretical prediction. We created multiple diffusion scenarios by varying intrinsic diffusivity inside and outside spheres, and membrane permeability.
* **Demo 3, randomly packed spheres, wide pulse sequence:** We simulate water diffusion in a medium composed of randomly packed permeable spheres, identify the time to the kurtosis peak, and compare it with the theoretical prediction. We created multiple geometries with varying sphere diameter, sphere volume fraction, and membrane permeability.

## References
* **Monte Carlo simulation**
  - [Fieremans, et al., NMR Biomed, 2010](https://doi.org/10.1002/nbm.1577)
  - [Novikov, et al., Nature Physics, 2011](https://doi.org/10.1038/nphys1936)
  - [Fieremans and Lee, NeuroImage 2018](https://doi.org/10.1016/j.neuroimage.2018.06.046)
  - [Lee, et al., Communications Biology 2020](https://doi.org/10.1038/s42003-020-1050-x)
  - [Lee, et al., NeuroImage 2020](https://doi.org/10.1016/j.neuroimage.2020.117228)
  - [Lee, et al., Journal of Neuroscience Methods 2021](https://doi.org/10.1016/j.jneumeth.2020.109018)
  - [Lee, et al., NMR in Biomedicine 2024](https://doi.org/10.1002/nbm.5087)
  - [Lee, et al., Magnetic Resonance in Medicine 2025](https://doi.org/10.1002/mrm.30335)

## Authors
* [Hong-Hsi Lee](https://www.martinos.org/investigator/hong-hsi-lee/)
* [Susie Y Huang](https://www.martinos.org/investigator/susie-huang/)
* [Dmitry S Novikov](http://www.diffusion-mri.com/people/dmitry-novikov)
* [Els Fieremans](http://www.diffusion-mri.com/people/els-fieremans)

## Acknowledgement
We would like to thank [Ricardo Coronado-Leija](https://scholar.google.com/citations?user=V5hykxgAAAAJ&hl=en) for the fruitful discussion of simulation implementation.

## License
This project is licensed under the [LICENSE](https://github.com/Connectome20/monte-carlo-simulation-sphere-PGSE/blob/main/LICENSE).
