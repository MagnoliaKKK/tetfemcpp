# 物理模拟稳定性修复

## 问题诊断

你遇到的圆环模型不稳定问题主要由以下原因导致：

1. **约束硬度过高**: `constraintHardness = 10000000000` 导致数值不稳定
2. **缺少固定锚点**: 移除所有固定点后整个系统失去约束
3. **PBD求解器参数不当**: 时间步长和阻尼设置不合理
4. **缺少数值稳定性检查**: 没有NaN检测和速度限制

## 修复方案

### 1. 参数调整
创建了稳定的参数配置：
- `constraintHardness`: 从 10000000000 降到 50000
- `timeStep`: 从 0.01 降到 0.008
- `dampingConst`: 从 20.0 改为 0.95
- `youngs`: 从 1000000 降到 100000

### 2. 数值稳定性检查
添加了以下安全机制：
- **NaN检测**: 自动检测并修复无效数值
- **速度限制**: 最大速度限制为50.0单位
- **位置限制**: 防止顶点偏离初始位置过远
- **重置功能**: 按R键重置整个模拟

### 3. 使用方法

#### 参数文件选择
有三个参数配置可选：
- `parameters.txt` - 修复后的默认配置
- `parameters_ring_stable.txt` - 专门为圆环优化的保守配置
- `parameters_armadillo.txt` - armadillo模型配置

#### 控制按键
- **R键**: 重置模拟到初始状态
- **W/A/S/D**: 原有的控制功能
- **C键**: 保存当前状态

#### 推荐设置 (圆环模型)
使用 `parameters_ring_stable.txt`:
```
youngs=50000
constraintHardness=1000  
timeStep=0.005
dampingConst=0.98
groupNumX=3
groupNumY=3  
groupNumZ=3
```

### 4. 实时监控
程序现在会输出以下警告信息：
- "NaN detected" - 发现无效数值
- "Velocity clamped" - 速度被限制
- "Extreme position detected" - 位置异常被修正

### 5. troubleshooting

如果仍然不稳定，尝试：

1. **降低约束硬度**:
   ```
   constraintHardness=100
   ```

2. **减少时间步长**:
   ```
   timeStep=0.003
   ```

3. **增加阻尼**:
   ```
   dampingConst=0.99
   ```

4. **减少分组数量**:
   ```
   groupNumX=2
   groupNumY=2
   groupNumZ=2
   ```

5. **简化网格**:
   ```
   tetgenArgs=pq2a0.05  // 更粗糙的网格
   ```

### 6. 调试技巧
- 将重力设为0开始测试: `Gravity=0.0`
- 观察控制台输出的稳定性警告
- 使用R键快速重置而不重启程序
- 逐步增加物理参数直到找到稳定值

## 技术原理

### PBD稳定性
Position-Based Dynamics需要合适的约束硬度，太高会导致：
- 数值积分误差放大
- 约束求解震荡
- 速度爆炸

### 无固定点系统
当移除所有固定点时：
- 系统失去参考框架
- 约束求解可能不收敛
- 需要更小的时间步长和更强的阻尼

这些修复应该能解决圆环模型的不稳定问题。