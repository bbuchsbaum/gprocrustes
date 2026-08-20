# Reference papers

Local copies of the papers that define the package’s obligations. Publisher copyright still applies; these files are for development and citation, not redistribution.

The two Downloads copies of Bai and Bartoli (`s11263-021-01571-8.pdf` and `s11263-021-01571-8 (1).pdf`) were byte-identical. Only one is stored.

| File | Citation | Role in gprocrustes |
|---|---|---|
| [Gower_1975_generalized_procrustes.pdf](Gower_1975_generalized_procrustes.pdf) | Gower, J. C. (1975). Generalized procrustes analysis. *Psychometrika* 40:33–51. DOI [10.1007/BF02291478](https://doi.org/10.1007/BF02291478) | Symmetric consensus GPA; translation / rotation / scale; energy decomposition; historical residual sequence |
| [Goodall_1991_procrustes_methods_shape.pdf](Goodall_1991_procrustes_methods_shape.pdf) | Goodall, C. (1991). Procrustes methods in the statistical analysis of shape. *JRSS B* 53:285–339. JSTOR [2345744](https://www.jstor.org/stable/2345744) | Shape space, \(O(d)\) vs \(SO(d)\), weighted metrics, \(\Sigma_S\) vs \(\Sigma_M\), robustness, tangent geometry |
| [Ling_2024_gopp_arbitrary_adversaries.pdf](Ling_2024_gopp_arbitrary_adversaries.pdf) | Ling, S. (2024). Generalized orthogonal Procrustes problem under arbitrary adversaries. arXiv [2106.15493v3](https://arxiv.org/abs/2106.15493) | GOPP NP-hardness, spectral init, GPM, SDR, dual certificate |
| [Bai_Bartoli_2022_procrustes_deformations.pdf](Bai_Bartoli_2022_procrustes_deformations.pdf) | Bai, F. and Bartoli, A. (2022). Procrustes analysis with deformations: a closed-form solution by eigenvalue decomposition. *IJCV* 130:567–593. DOI [10.1007/s11263-021-01571-8](https://doi.org/10.1007/s11263-021-01571-8) | LBW / affine / TPS, reference-space constraints, free-translations, partial shapes, CV |
| [Pizarro_Bartoli_2011_global_gpa.pdf](Pizarro_Bartoli_2011_global_gpa.pdf) | Pizarro, D. and Bartoli, A. (2011). Global optimization for optimal generalized Procrustes analysis. *CVPR*. DOI [10.1109/CVPR.2011.5995670](https://doi.org/10.1109/CVPR.2011.5995670) | SOS / SDP global formulation for 2-D and 3-D GPA with missing data |
| [Igual_etal_2014_continuous_gpa.pdf](Igual_etal_2014_continuous_gpa.pdf) | Igual, L., Perez-Sala, X., Escalera, S., Angulo, C., and De la Torre, F. (2014). Continuous generalized Procrustes analysis. *Pattern Recognition* 47:659–671. DOI [10.1016/j.patcog.2013.08.006](https://doi.org/10.1016/j.patcog.2013.08.006) | Continuous / view-parameterized GPA; bias from discrete 3-D sampling of 2-D landmarks |

SHA-256 checksums (of the stored files):

```
2a1e96de78d38affb8ec9e31ef960328b40ccca91ceee08a9990507913737fee  Gower_1975_generalized_procrustes.pdf
240ac787bb7463027487791d3485e41f5e7d870a301fc2fea94cacd81660b6c0  Goodall_1991_procrustes_methods_shape.pdf
a384377bb4a737f5878d6446d3665e9bfbc7565e80b8751ae9445d6f4165674e  Ling_2024_gopp_arbitrary_adversaries.pdf
7a69dd3b9706b8403984dec75b03fdb9b5e84d68e5387b990632b1a5075d7c92  Bai_Bartoli_2022_procrustes_deformations.pdf
7ac75ef0613d5e9020f8667f11e95f055d0a0447ce4014b3b2c274c07ff84586  Pizarro_Bartoli_2011_global_gpa.pdf
d4d234425202fa09ae7212953db34d5e92f396939e0230184e459049914a43c9  Igual_etal_2014_continuous_gpa.pdf
```

The first four papers are the charter’s four obligations. Pizarro–Bartoli and Igual et al. are supporting references for global 2-D/3-D GPA and continuous shape models.
