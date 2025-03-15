%% Note: A more convenient toolbox will be released later; here we will first present the main code.

%% Input Data: 
%(1) GRACE Sperical Harmonic Coefficients (SHCs): each column represents
%SHCs of each month, 
load ITSG_vv1.mat %ITSG SHCs after removing mean field 2004.01-2009.12, replacing C20 and C30 terms. 

%(2) Converts geoid coefficients (gc) to mass coefficients (mc)
load ITSG_T.mat

%(3) Covariance matrices of regional mass changes propagated from those of
%SHCs after expanding 2 degrees
load QQ_yt_buffer_degree2.mat %here we take Yangzte basin as an example

%(4) Latitude (lat), longtitude (lon), row (ll), column (cc), region range (msk), converting matrix
%from GRACE SHCs to regional grids (MM).
load cs2gridM_yt_buffer_degree2.mat

%(5) GIA_new: GIA read from ICE6G-D model

%(6) time_new: 2002.04-2023.12

%% Subfunctions:
MMSE_Tik: Deriving regularization parameter based on the minimum MSE
Kalman_forward_region_rr_multistep_QQ1: Forward KF process
Kalman_back_region_rr_multistep_QQ1: Backward KF process
Ecov: Updating the signal covariance matrix of parameters; the inverse of regularization matrix.

Output data
Denoised deterministic signals (linear trends, a amd irregular signals.
