#极射赤面投影：
        #在任意纬度lat的放大系数m=1+sin(60)/(1+sin(lat))
        #K=1
        #l_eq=a(1+(0.75)**0.5)
        #A=(l/l_eq)**(2/K)
        #lat=arcsin((1-A)/(1+A))； 
        # lat→A→l,求出l关于lat的函数       
# Lambert投影：        
        #在任意纬度lat的放大系数m=1.28303*sec(lat)*(cos(lat)/(1+sin(lat)))^0.71557
        #l_eq=(a*cos(lat_0)/K)*((1+sin(lat_0)/cos(lat_0)))**K,LAT_0=30°
        #A=(l/l_eq)**(2/K)
        #K=0.71557
        #lat=arcsin((1-A)/(1+A))；        
# Mercator投影：
        #K=0
        #m=cos(lat_0)/cos(lat)
        #lat_0=22.5°
        #l=a*cos(lat_0)*ln(cos(lat)/(1-sin(lat)))



import numpy as np
import matplotlib.pyplot as plt

# 地球参数
R = 6371.0  # 地球半径，单位：千米
Omega = 7.2921e-5  # 地球自转角速度（仅作参考）

def create_plots():
    plt.figure(figsize=(12, 8))
    
    # ========================
    # 极射赤面投影 (Stereographic)
    # ========================
    lats_stere = np.linspace(0, 89, 100)  # 避免90°奇点
    ms_stere = []
    ls_stere = []
    
    # 常数预计算
    sin60 = np.sin(np.radians(60))
    l_eq_stere = R * (1 + np.sqrt(3)/2)  # 1 + √0.75 ≈ 1.866
    
    for lat in lats_stere:
        # 放大系数
        sin_phi = np.sin(np.radians(lat))
        m = (1 + sin60) / (1 + sin_phi)
        ms_stere.append(m)
        
        # 极距
        l = l_eq_stere * np.sqrt((1 - sin_phi)/(1 + sin_phi))
        ls_stere.append(l)
    
    # ========================
    # Lambert投影
    # ========================
    lats_lcc = np.linspace(0, 89, 100)
    ms_lcc = []
    ls_lcc = []
    
    # 投影参数
    K = 0.71557
    lat_0 = 30
    cos_lat0 = np.cos(np.radians(lat_0))
    sin_lat0 = np.sin(np.radians(lat_0))
    term_l_eq = (1 + sin_lat0)/cos_lat0
    l_eq_lcc = (R * cos_lat0 / K) * (term_l_eq ** K)
    
    for lat in lats_lcc:
        sin_phi = np.sin(np.radians(lat))
        cos_phi = np.cos(np.radians(lat))
        
        # 放大系数
        term = cos_phi / (1 + sin_phi)
        m = 1.28303 * (term ** K) / cos_phi
        ms_lcc.append(m)
        
        # 极距
        A = (1 - sin_phi)/(1 + sin_phi)
        l = l_eq_lcc * (A ** (K/2))
        ls_lcc.append(l)
    
    # ========================
    # Mercator投影
    # ========================
    lats_merc = np.linspace(0, 85, 100)  # 避免极区奇点
    ms_merc = []
    ls_merc = []
    
    lat_0 = 22.5
    cos_lat0_merc = np.cos(np.radians(lat_0))
    
    for lat in lats_merc:
        # 放大系数
        m = cos_lat0_merc / np.cos(np.radians(lat))
        ms_merc.append(m)
        
        # 极距
        if lat == 0:
            l = 0.0
        else:
            lat_rad = np.radians(lat)
            l = R * cos_lat0_merc * np.log(np.tan(np.pi/4 + lat_rad/2))
        ls_merc.append(abs(l))  # 取绝对值保证正距离
    
    # ========================
    # 绘制图形
    # ========================
    ax1 = plt.subplot(2, 1, 1)
    ax2 = plt.subplot(2, 1, 2)
    
    # 放大系数曲线
    ax1.plot(lats_stere, ms_stere, label='Polar Stereographic (60°N)')
    ax1.plot(lats_lcc, ms_lcc, label='Lambert (30°N/60°N)')
    ax1.plot(lats_merc, ms_merc, label='Mercator (22.5°)')
    
    # 极距曲线
    ax2.plot(lats_stere, ls_stere, label='Polar Stereographic')
    ax2.plot(lats_lcc, ls_lcc, label='Lambert')
    ax2.plot(lats_merc, ls_merc, label='Mercator')
    
    # 图形修饰
    ax1.set_title('Map Scale Factor (m) vs Latitude')
    ax1.set_ylabel('Scale Factor')
    ax1.grid(True)
    ax1.legend()
    ax1.set_ylim(0, 10)  # 限制y轴范围
    ax1.set_ylim(0.8, 2.0)  # 限制y轴范围

    
    ax2.set_title('Polar Distance (l) vs Latitude')
    ax2.set_xlabel('Latitude (°)')
    ax2.set_ylabel('Distance (km)')
    ax2.grid(True)
    ax2.legend()
    
    plt.tight_layout()
    plt.show()

create_plots()