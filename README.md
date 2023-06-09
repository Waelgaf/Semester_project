# Exploring methods in multi-layer clustering
 This is the git repository of the source code in _R_ of the implementation of three methods of multi-layer clustering. Additionally, there are the corresponding codes to run simulations to evaluate the performances of these methods and a real data set analysis.
 ## Code
 ### Least Square Approach
 Under the __Algorithm1__ file, the Least Square approach has to be run under _MSBM_ framework only. You can find out more about it by looking at the report attached to this repository or at the direct article related to this method [Lei et al.](https://www.stat.pitt.edu/khchen/authorversions/nettensor.pdf).
 ### TWIST
 Under the __Algorithm2__ file, the main function ``PowerInitialization`` and its relative subfunctions ``norm_vec`` and ``reg_vec`` have been imported from this [github repository](https://github.com/ChenyuzZZ73/rMultiNet). This approach has to be run under a _MMSBM_ or _MSBM_ framework. You can find out more about it by looking at the report attached to this repository or at the direct article related to this method [Jing et al.](https://arxiv.org/abs/2002.04457).
 ### ALMA
 Under the __Algorithm3__ file, you retrieve the ALMA method. Moreover, several tensor and matrix [operators](https://www.kolda.net/publication/TensorReview.pdf) have been implemented and are free to be used. This approach has to be run under_MMSBM_ framework only. You can find out more about ALMA by looking at the report to this repository or at the direct article related to this method [Fan et al.](https://arxiv.org/abs/2102.10226).
 ### Remaining code
 The remaining code concerns directly the report attached to this repository. There are implementation of:
 * Generation of networks (more details in the report)
 * Measurement of performances
 * Simulation study
 * and, real data set analysis.
 # References
 * "Jing Lei, Kehui Chen, and Brian Lynch. Consistent community detection in multi-layer network data.
Biometrika, 107(1):61–73, 2019. doi: https://doi.org/10.1093/biomet/asz068."
* "Bing-Yi Jing, Ting Li, Zhongyuan Lyu, and Dong Xia. Community detection on mixture multi-layer networks
via regularized tensor decomposition. 2020. doi: https://doi.org/10.48550/arXiv.2002.04457"
* "Xing Fan, Marianna Pensky, Feng Yu, and Teng Zhang. Alma: Alternating minimization algorithm for clus-
tering mixture multilayer network. 2021. doi: https://doi.org/10.48550/arXiv.2102.10226."

 
