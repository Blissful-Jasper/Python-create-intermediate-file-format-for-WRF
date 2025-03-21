# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 22:41:36 2023

@author: jianpu

@blog :  https://blog.csdn.net/weixin_44237337?spm=1000.2115.3001.5343

@email: 211311040008@hhu.edu.cn

introduction : keep learning althongh walk slowly
"""

import numpy as np
import pandas as  pd
import numpy as np
from netCDF4 import Dataset
import pywinter.winter as pyw
#%% calculate

import xarray as xr

# data = xr.open_dataset("/Users/xpji/WRF/WPS/fort1201/met_em.d01.2004-06-19_00:00:00.nc")
# pres = data['PRES']
# pres = pres[0].data



#%%




nlev = 26
nlat = 181
nlon = 360
g = 9.81
# gas = 287.04
p_mean = np.zeros(nlev)
h_mean = np.zeros(nlev)
hft_mean = np.zeros(nlev)
t_mean = np.zeros(nlev)
r_mean = np.zeros(nlev)
theta_mean = np.zeros(nlev)
theta_e_mean = np.zeros(nlev)
rh_mean = np.zeros(nlev)
q_mean = np.zeros(nlev)

psfc_mean = 0.0
tsfc_mean = 0.0
rsfc_mean = 0.0
thsfc_mean = 0.0
thesfc_mean = 0.0
rhsfc_mean = 0.0
qsfc_mean = 0.0


u_file = "/Users/xpji/work_for_myself/MRG_wave_model/MRG_py_code/u3d.npy"
v_file = "/Users/xpji/work_for_myself/MRG_wave_model/MRG_py_code/v3d.npy"
h_file = "/Users/xpji/work_for_myself/MRG_wave_model/MRG_py_code/h3d.npy"
usfc_file = "/Users/xpji/work_for_myself/MRG_wave_model/MRG_py_code/u.npy"
vsfc_file = "/Users/xpji/work_for_myself/MRG_wave_model/MRG_py_code/v.npy"
hsfc_file = "/Users/xpji/work_for_myself/MRG_wave_model/MRG_py_code/h.npy"

u3d_d = np.load(u_file).transpose((2,1,0))
v3d_d = np.load(v_file).transpose((2,1,0))
h3d_d = np.load(h_file).transpose((2,1,0))

# Surface Structure
usfc2d_d = np.load(usfc_file).transpose((1,0))
vsfc2d_d = np.load(vsfc_file).transpose((1,0))
hsfc2d_d = np.load(hsfc_file).transpose((1,0))



v3d     = np.zeros((nlev,   nlat, nlon,))
u3d     = np.zeros((nlev,   nlat, nlon, ))
h3d     = np.zeros((nlev,   nlat, nlon,))
rh3d    = np.zeros((nlev,   nlat, nlon,))
t3d     = np.zeros((nlev,   nlat, nlon,))


plev = np.zeros(nlev)

h1d = np.zeros(nlev)

# pmsl = np.zeros((nlon, nlat))
# psfc = np.zeros((nlon, nlat))
# tsfc = np.zeros((nlon, nlat))

# pz = 0.0

t0 = 0.0
dtdz = 0.0
rd = 287.05
p0 = 0.0
cal_slp = 0.0

fi = '/Users/xpji/work_for_myself/MRG_wave_model/MRG_py_code/Mean_Sounding-hurricane.csv'
df = pd.read_csv(fi,header=None)

"""
df  from high to low 

100hpa,200hpa,,,,,,,1000hpa


"""
p_mean = (df.iloc[:26, 0]).values[::-1]
h_mean = (df.iloc[:26, 1]).values[::-1]
hft_mean = (df.iloc[:26, 2]).values[::-1]
t_mean = (df.iloc[:26, 3]).values[::-1]
r_mean = (df.iloc[:26, 4]).values[::-1]
theta_mean = (df.iloc[:26, 5]).values[::-1]
theta_e_mean = (df.iloc[:26, 6]).values[::-1]
rh_mean = (df.iloc[:26, 7]).values[::-1]
q_mean = (df.iloc[:26, 8]).values[::-1]
  
psfc_mean = (df.iloc[26, 0])
tsfc_mean = (df.iloc[26, 3])
rsfc_mean = (df.iloc[26, 4])
thsfc_mean = (df.iloc[26, 5])
thesfc_mean = (df.iloc[26, 6])
rhsfc_mean = (df.iloc[26, 7])
qsfc_mean = (df.iloc[26, 8])


# Construct 3D Fields

for k in range(nlev):
    print(k)
    plev[k] = p_mean[k] * 100.0
    u3d[k,:,:] = u3d_d[k,:,:,]
    v3d[k] = v3d_d[k,:,:,]
    h3d[k] = h_mean[k] / g * 9.8 + h3d_d[k,:,:,]  # ??
    rh3d[k,:,:,] = rh_mean[k]
    
def hy2t( h, p,nlev):
    """
    Isothermal P-H expression conversion from Fortran to Python.
    """
    t = np.zeros(nlev)
    g = 9.81
    gas = 287.04

    for k in range(nlev):
        # print(k)
        
        if k == 0:
            t[k] = -g / gas * (h[k+1] - h[k]) / np.log(p[k+1] / p[k])
        elif k == 25:
            t[k] = -g / gas * (h[k] - h[k - 1]) / np.log(p[k] / p[k - 1])
        else:
            t[k] = -g / gas * (h[k+1] - h[k - 1]) / np.log(p[k+1] / p[k - 1])
    return t   
# Diagnose temperature, isothermal P-H expression
for i in range(nlat):
    for j in range(nlon):
        
        h1d = h3d[ :,i, j]
        t1d = hy2t(h1d, plev, nlev) 
        t3d[:,i, j] = t1d
        

        
def cal_slp(t0, pz, hz, dtdz, g, rd):
    cal_slp = pz * (t0 / (t0 - dtdz * hz)) ** (g / (rd * dtdz))
    return cal_slp

# Calculate sea level pressure

pz = plev[0]
dtdz = (t_mean[0] - tsfc_mean) / h3d[0]
tsfc = t3d[0] - dtdz * h3d[ 0]
pmsl= cal_slp(tsfc, pz, h3d[ 0], dtdz, g, rd)


# For idealized case, set terrain = 0, so we can set psfc = pmsl
psfc = pmsl

# np.save("./u3d.npy", u3d[:, :, 0])
# np.save("./v3d.npy", v3d[:, :, 0])
# np.save("./pmsl.npy", pmsl)

np.sum(np.isnan(psfc))

################################################################################
################################################################################
################################################################################
import matplotlib.pyplot as plt
figure = plt.figure(dpi=200,figsize=(10,8))
plt.subplot(121)
plt.imshow(u3d[:,:,0])
plt.subplot(122)
plt.imshow(v3d[:,:,1])
# Read Geo-data
lat = np.load('/Users/xpji/work_for_myself/MRG_wave_model/MRG_py_code/lat.npy')
lon = np.load('/Users/xpji/work_for_myself/MRG_wave_model/MRG_py_code/lon.npy')
dlat = np.abs(lat[1]-lat[0])
dlon = np.abs(lon[1]-lon[0])

# Read 2D data: PSML\ PSFC\ TT\ UU \ VV

# Read 3D data
uu  = (u3d)
vv  = (v3d)
ght = (h3d)
tt  = (t3d)
rh  = (rh3d)
da_skintemp = xr.open_dataset('/Users/xpji/code/soil_data_for_ideal_model/download_skin_temp.nc')
skintemp =  np.full((181,360),299.7)

def get_sst_linear_decrase():
     
    # 创建一个181x360的全0数组
    sst = np.zeros((181, 360))
    
    # 设置赤道的SST为301K
    sst[90,:] = 301
    
    # 设置北极和南极的SST为293K
    sst[0,:] = 271
    sst[-1,:] = 271
    
    # 对于每个纬度值介于0和90之间的行，使用线性插值计算SST值
    for lat in range(1,90):
        # 计算该行的SST值
        sst_value = 301 - (301-271)*(lat/90)
        # 将SST值应用于该行的所有列
        sst[90+lat,:] = sst_value
        sst[90-lat,:] = sst_value
    return sst

sst      =  get_sst_linear_decrase()

plt.imshow(sst)
############################################################
# Create winter fields
print('tsfc',tsfc.shape)
print('u10m',usfc2d_d.shape)
print('v10m',vsfc2d_d.shape)
print('psfc',psfc.shape)
print('pmsl',pmsl.shape)

winter_geo  = pyw.Geo0(-90,0,dlat,dlon)
winter_t2m  = pyw.V2d('TT',tsfc)
winter_u10  = pyw.V2d('UU',usfc2d_d)
winter_v10  = pyw.V2d('VV',vsfc2d_d)
winter_psfc = pyw.V2d('PSFC',psfc)
winter_pmsl = pyw.V2d('PMSL',pmsl)
winter_sktp = pyw.V2d('SKINTEMP',skintemp)
winter_sst  = pyw.V2d('SST',sst)

winter_UU   = pyw.V3dp('UU',uu,plev)
winter_VV   = pyw.V3dp('VV',vv,plev)
winter_GHT  = pyw.V3dp('GHT',ght,plev)
winter_TT   = pyw.V3dp('TT',tt,plev)
winter_RH   = pyw.V3dp('RH',rh,plev)

da_soil = xr.open_dataset("/Users/xpji/code/soil_data_for_ideal_model/download.nc").sortby('latitude')
stl1 = da_soil['stl1'][0].data
stl2 = da_soil['stl2'][0].data
stl3 = da_soil['stl3'][0].data
stl4 = da_soil['stl4'][0].data

swvl1 = da_soil['swvl1'][0].data
swvl2 = da_soil['swvl2'][0].data
swvl3 = da_soil['swvl3'][0].data
swvl4 = da_soil['swvl4'][0].data

st = np.stack([stl1, stl2, stl3, stl4])
sm = np.stack([swvl1,swvl2, swvl3, swvl4], axis=0)

# 

slt_layer   = ['000007', '007028', '028100', '100289']
winter_sm   = pyw.Vsl('SM',sm,slt_layer)
winter_st   = pyw.Vsl('ST',st,slt_layer)

# Listing fields

total_fields = [
    winter_t2m,winter_u10,winter_v10,
    winter_psfc,winter_pmsl,
                winter_UU, winter_VV,winter_GHT,winter_TT ,winter_RH ,
                    winter_sm ,
                    winter_st ,
                   winter_sktp,winter_sst 
                ]


# Out path
path_out = '/Users/xpji/work_for_myself/MRG_wave_model/MRG_py_code/'

# Write intermediate file

pyw.cinter('FILE','2021-01-01_00',winter_geo,total_fields,path_out)





interf  =  pyw.rinter('FILE:2021-01-01_00')
print(interf['UU'].general)
print(interf['TT'].level)
print(interf['UU'].level)
print(interf['TT'].general)
print(interf['TT'].val.shape)
# print(interf['P'])

# print(pyw.rinter('/Users/xpji/MRG_wave_model/MRG_py_code/intermediate_files/FILE:1985-08-11_00').keys())
