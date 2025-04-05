import os 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.animation import FuncAnimation

output_dir = 'output'
if not os.path.exists(output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
# 中文字体设置
# 设置字体为 SimHei（黑体）以支持中文
rcParams['font.sans-serif'] = ['SimHei']
rcParams['axes.unicode_minus'] = False


# 参数设置
dx = 400
dt = 2
num_steps = 300
m_points = 360
m_indices = np.arange(1, m_points + 1)
u0 = 20 * np.cos(np.deg2rad(3 * m_indices))

def simulate_advection(c):
    """数值求解平流方程"""
    u_history = np.zeros((num_steps + 1, m_points))
    u_history[0] = u0
    
    # 第一步用前向差分
    u_history[1] = u0 - (c * dt / (2 * dx)) * (np.roll(u0, -1) - np.roll(u0, 1))
    
    # 后续时间步用中央差分
    for n in range(1, num_steps):
        u_history[n+1] = u_history[n-1] - (c * dt / dx) * (
            np.roll(u_history[n], -1) - np.roll(u_history[n], 1))
    return u_history

# 运行模拟
c_values = [20, 210]
results = {c: simulate_advection(c) for c in c_values}

# 要观察的格点位置（Python索引从0开始）
selected_points = {
    'm60': 59,
    'm100': 99,
    'm120': 119,
    'm140': 139
}

# ====================== 热图可视化 ======================
def plot_heatmap(data, c):
    plt.figure(figsize=(12, 6))
    plt.pcolormesh(
        np.arange(m_points)+1, 
        np.arange(num_steps+1)*dt,
        data,
        shading='auto',
        cmap='RdBu_r',
        vmin=-25, 
        vmax=25
    )
    plt.colorbar(label='Velocity (m/s)')
    plt.title(f'c = {c} m/s (CFL={c*dt/dx:.2f}) - 全场演化热图')
    plt.xlabel('格点位置')
    plt.ylabel('时间 (s)')
    plt.savefig(os.path.join(output_dir, f'heatmap_c_{c}.png'))
    plt.close()

# ====================== 合并曲线图 ======================
def plot_combined_curves(data, c):
    plt.figure(figsize=(10, 5))
    time = np.arange(num_steps+1)*dt
    for label, idx in selected_points.items():
        plt.plot(time, data[:, idx], label=f'格点 {label[1:]}')
    plt.title(f'c = {c} m/s - 关键点速度时序对比')
    plt.xlabel('时间 (s)')
    plt.ylabel('速度 (m/s)')
    plt.grid(True)
    plt.legend()
    plt.ylim(-25, 25)
    plt.savefig(os.path.join(output_dir, f'combined_curves_c_{c}.png'))
    plt.close()

# ====================== 新增动画生成函数 ======================
def create_wave_animation(u_history, c, filename):
    fig = plt.figure(figsize=(10, 5))
    ax = plt.axes(xlim=(1, m_points), ylim=(-30, 30))
    line, = ax.plot([], [], lw=2)
    
    ax.set_xlabel('格点位置')
    ax.set_ylabel('速度 (m/s)')
    ax.set_title(f'c = {c} m/s - 波形东移动画 (CFL={c*dt/dx:.2f})')
    ax.grid(True)
    
    x = np.arange(1, m_points+1)
    def init():
        line.set_data(x, u_history[0])
        return line,
    
    def update(frame):
        line.set_ydata(u_history[frame])
        ax.set_title(f'c = {c} m/s - 时间: {frame*dt}s')
        return line,
    
    ani = FuncAnimation(fig, update, frames=num_steps,  # 减少帧数提升性能
                        init_func=init, blit=True, interval=30)
    ani.save(os.path.join(output_dir, filename), writer='pillow')  # 保存动画
    plt.close()


# ====================== 执行可视化 ======================
for c in c_values:
    data = results[c]
    
    # 热图
    plot_heatmap(data, c)
    #保存热图        

    # 合并曲线图
    plot_combined_curves(data, c)

    # 生成动画
    create_wave_animation(
        data, c,
        filename=f'wave_animation_c_{c}.gif'
    )

