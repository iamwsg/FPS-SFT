# MARS-SFT 

We propose MARS-SFT, a sparse Fourier transform for multidimensional, frequency-domain sparse signals, inspired by the idea of the Fourier projection-slice theorem. 

MARS-SFT identifies  frequencies by operating on one-dimensional slices of the discrete-time domain data, taken along specially designed lines; those lines are parametrized by slopes that are randomly generated from a set at runtime. The DFTs of the data slices  represent DFT projections  onto the  lines along which the slices were taken. When the lines' lengths and slopes are properly designed so that they allow for orthogonal and uniform frequency projections,  the multidimensional frequencies can be recovered from their projections with low sample and computational complexity. The large number of degrees of freedom in frequency projections allows MARS-SFT to be applicable to less sparse data with non-uniformly distributed frequencies. We also extend MARS-SFT into a robust version to address noisy signals that contain off-grid frequencies.  
Please see the following papers for more details.

[1] Wang, Shaogang, Vishal M. Patel, and Athina Petropulu. "FPS-SFT: A Multi-dimensional Sparse Fourier Transform Based on the Fourier Projection-slice Theorem." arXiv preprint arXiv:1711.11407 (2017).

[2] Wang, Shaogang, Vishal M. Patel, and Athina Petropulu. "Robust sparse fourier transform based on the fourier projection-slice theorem." Radar Conference (RadarConf18), 2018 IEEE. IEEE, 2018.

Please use the following entires if you want to cite those papers. 

@article{wang2017fps,
  title={FPS-SFT: A Multi-dimensional Sparse Fourier Transform Based on the Fourier Projection-slice Theorem},
  author={Wang, Shaogang and Patel, Vishal M and Petropulu, Athina},
  journal={arXiv preprint arXiv:1711.11407},
  year={2017}
}

@inproceedings{wang2018robust,
  title={Robust sparse fourier transform based on the fourier projection-slice theorem},
  author={Wang, Shaogang and Patel, Vishal M and Petropulu, Athina},
  booktitle={Radar Conference (RadarConf18), 2018 IEEE},
  pages={1427--1432},
  year={2018},
  organization={IEEE}
}




