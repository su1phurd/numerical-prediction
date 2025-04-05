import numpy as np

# 地球参数
R = 6371.0e3  # 地球半径，单位：米
Omega = 7.2921e-5  # 地球自转角速度，单位：rad/s

def calculate_parameters(case):
    if case == 1:
        # 极射赤面投影
        d = 500e3  # 网格距转换为米
        In, Jn = -4, 8
        
        # 计算极距
        l = np.sqrt((In*d)**2 + (Jn*d)**2)
        a = R
        
        # 标准纬线60°N的辅助参数
        sin60 = np.sin(np.radians(60))
        l_eq = a * (1 + sin60)
        
        # 计算纬度
        K = 1
        A = (l / l_eq) ** (2/K)
        sin_phi = (1 - A) / (1 + A)
        phi = np.degrees(np.arcsin(sin_phi))
        
        # 放大系数
        m = (1 + sin60) / (1 + np.sin(np.radians(phi)))
        
        # 科氏参数
        f = 2 * Omega * np.sin(np.radians(phi))
        
        return l/1e3, m, f, phi

    elif case == 2:
        # Lambert投影
        d = 300e3
        In, Jn = 5, 15
        a = R
        lat_0 = 30  # 标准纬线30°N
        K = 0.71557
        
        # 计算极距
        l = np.sqrt((In*d)**2 + (Jn*d)**2)
        
        # 计算l_eq
        cos_lat0 = np.cos(np.radians(lat_0))
        sin_lat0 = np.sin(np.radians(lat_0))
        term = (1 + sin_lat0) / cos_lat0
        l_eq = (a * cos_lat0 / K) * (term ** K)
        
        # 计算纬度
        A = (l / l_eq) ** (2/K)
        sin_phi = (1 - A) / (1 + A)
        phi = np.degrees(np.arcsin(sin_phi))
        
        # 放大系数
        cos_phi = np.cos(np.radians(phi))
        sin_phi_val = sin_phi
        m = 1.28303 * (cos_phi/(1 + sin_phi_val)) ** K / cos_phi
        
        # 科氏参数
        f = 2 * Omega * sin_phi_val
        
        return l/1e3, m, f, phi

    elif case == 3:
        # Mercator投影
        d = 200e3
        Je = 3
        a = R
        lat_0 = 22.5  # 标准纬线22.5°
        
        # 计算东向坐标
        x = Je * d
        
        # 计算纬度
        A = x / (a * np.cos(np.radians(lat_0)))
        B = np.exp(2 * A)
        sin_phi = (B - 1) / (B + 1)
        phi = np.degrees(np.arcsin(sin_phi))
        
        # 放大系数
        m = np.cos(np.radians(lat_0)) / np.cos(np.radians(phi))
        
        # 计算极距（沿经线的距离）
        cos_phi = np.cos(np.radians(phi))
        sin_phi_val = sin_phi
        log_term = np.log(cos_phi / (1 - sin_phi_val))
        l = a * np.cos(np.radians(lat_0)) * log_term
        
        # 科氏参数
        f = 2 * Omega * sin_phi_val
        
        return l/1e3, m, f, phi

# 计算三种情况
cases = []
for case in [1, 2, 3]:
    l, m, f, phi = calculate_parameters(case)
    cases.append((l, m, f, phi))

# 输出结果
for i, (l, m, f, phi) in enumerate(cases, 1):
    print(f"情况 {i}:")
    print(f"  极距 l = {l:.2f} km")
    print(f"  放大系数 m = {m:.5f}")
    print(f"  Coriolis参数 f = {f:.5e} s⁻¹")
    print(f"  纬度 φ = {phi:.5f}°\n")