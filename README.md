#  A Fast Multiscale Spatial Regularization for Sparse Hyperspectral Unmixing    #



This is the authors' implementation of [1].

The code is implemented in MATLAB and includes:  
-  demo1_MUA.m               - a demo script comparing the algorithms (DC1)  
-  demo2_MUA.m               - a demo script comparing the algorithms (DC2)  
-  demo3_MUA.m               - a demo script comparing the algorithms (DC3)  
-  demo4_supp_MUA.m          - demo script of the supplemental material  
-  demo5_supp_MUA.m          - demo script of the supplemental material  
-  demo6_supp_MUA.m          - demo script of the supplemental material  
-  sort_library_by_angle.m   - sort the signatures in the spectral library  
-  prune_library.m           - remove correlated signatures from the library  
-  tight_subplot.m           - more efficient subplot  
-  soft.m                    - soft thresholding operator function  
-  sunsal.m                  - the SUnSAL algorithm  
-  sunsal_tv.m               - the SUnSAL-TV algorithm  
-  sunsal_tv_lw_sp.m         - the S2WSU algorithm  
-  sunsal_spreg.m            - sparse unmixing at the fine spatial scale  
-  ./real_data/              - real images and spectral libraries used in the supplemental experiments  
-  ./vlfeat-0.9/             - the VLFeat toolbox (used for the SLIC superpixels alg.)  
-  ./HSI_segmentation/       - Binary partition tree HSI segmentation algorithm  
-  README                    - this file  



## IMPORTANT:
If you use this software please cite the following in any resulting
publication:

    [1] A Fast Multiscale Spatial Regularization for Sparse Hyperspectral Unmixing
        R.A. Borsoi, T. Imbiriba, J.C.M. Bermudez, C. Richard.
        IEEE Geoscience and Remote Sensing Letters, 2018.



## INSTALLING & RUNNING:
Just start MATLAB and run demo1_MUA.m, demo2_MUA.m or demo3_MUA.m.

For the demos in the supplemental material, just run demo4_supp_MUA.m, demo5_supp_MUA.m, or demo6_supp_MUA.m. The data (spectral libraries and images) for these examples is included in the "real_data" folder.

The Cuprite HS image is not included, and should be downloaded from http://www.lx.it.pt/~bioucas/code.htm.

### Important Notice About VLFeat
If you encounter problems with the "vl_slic" or "vl_setup" functions, try to download the latest version of the toolbox at http://www.vlfeat.org/install-matlab.html. 


## NOTES:
1.  The SUnSAL and SUnSAL-TV algorithms are provided by Jose Bioucas Dias 
    at http://www.lx.it.pt/~bioucas/code.htm.

2.  The S2WSU algorithm was provided by Shaoquan Zhang, and corresponds to the publication:
    S. Zhang, J. Li, H.-C. Li, C. Deng and A. Plaza. 
    Spectral-Spatial Weighted Sparse Regression for Hyperspectral Image Unmixing. 
    IEEE Transactions on Geoscience and Remote Sensing, 2018

3.  The vl_feat toolbox is included for the SLIC algorithm implementation
    http://www.vlfeat.org/

4.  The Binary Partition Tree HSI segmentation algorithm was provided by Miguel Veganzones.
    Veganzones, M. A., Tochon, G., Dalla-Mura, M., Plaza, A. J., & Chanussot, J.
    Hyperspectral image segmentation using a new spectral unmixing-based binary
    partition tree representation.
    IEEE Transactions on Image Processing, 2014.



