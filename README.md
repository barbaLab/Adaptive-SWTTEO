# Adaptive-SWTTEO
Matlab-based single-channel spike detection algorithm for extracellular neural recordings. 

The algorithm is composed of these steps: 
1. Stationary Wavelet transform (SWT) of n-th level. 
2. Mantaining only approximation coefficient at each level.
3. k-Teager Energy Operator (or k-TEO/k-NEO).  
4. Smoothing (k-SNEO).
5. Summation of all the transformed signals at each level.
6. Threshold with a windowed quantile function.
7. Peak refinement.
   
The code is inspired from the SWTTEO algorithm presented in _Florian Lieb, Hans-Georg Stark, and Christiane Thielemann. “A Stationary Wavelet Transform and a Time-Frequency Based Spike Detection Algorithm for Extracellular Recorded Data”. In: J Neural Eng 14.3 (June 2017). DOI: 10.1088/1741-2552/aa654b_. Key modifications are made from the third step till the end. Since SWTTEO had been benchmarked against multiple existing single-channel spike detection algorithms, demonstrating superior performance, our improved algorithm establishes itself as a leading solution in the field.

![AdSWTTEO Schema](./AdSWTTEO_schema1)
