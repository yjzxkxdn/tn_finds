# tn_finds

面向 UTAU 的独立重采样器。保留经典位置参数和常用 flag；合成核心独立实现，运行时不加载原版 `tn_fnds`、WORLD、Python 或外部模型文件。0.4.9 的 Windows 构建静态包含 Nebula 源 F0 模型，并读写 Moresampler 的 `desc.mrq` 源 F0 轨迹。

第一次接触可以先读 [科普介绍](INTRODUCTION.md)。

## 名字

原版 `tn_fnds` 的名字来自 *tanin fundoshi*（他人のふんどし，意为“借别人的兜裆布”），指依托前人的成果。`tn_finds` 保留 **tn**（*tanin*），把 **fnds** 改成 **finds**（*find* 的第三人称单数），表示寻找、发掘并做出更好的方案。代码是独立实现，不包含原版引擎源码。

## 开发与发行

功能开发在保留完整历史的原始仓库进行；本仓库只接收整理后的发布快照，发行只从 clean 发布线（`main`）出。详见 [DEVELOPMENT.md](DEVELOPMENT.md)。

## 构建

需要 C++17 编译器。

### MinGW x64：标准 F0

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\build.ps1
```

生成 `tn_finds_v0.4.9-resampler.exe`（独立重采样）和 `tn_finds_v0.4.9-phrase.exe`（双角色整句）。两个 EXE 都是静态链接。这个直接编译路径不带 Nebula，使用窗校正自相关 F0。

### CMake：官方 Nebula 构建

0.4.9 官方 Windows EXE 使用 OpenMoresampler 的 Nebula 源 F0 后端。把 OpenMoresampler 源码放在本仓库之外，检出 `15a4b85eaac4944a0a55cae50de96fb5a6f15344`：

```powershell
git clone https://github.com/LingYi0110/OpenMoresampler.git ..\OpenMoresampler
git -C ..\OpenMoresampler checkout 15a4b85eaac4944a0a55cae50de96fb5a6f15344
```

然后配置：

```powershell
cmake -S . -B build-nebula -G Ninja -DCMAKE_BUILD_TYPE=Release `
    "-DTN_FINDS_OPENMORESAMPLER_SOURCE:PATH=$PWD\..\OpenMoresampler" -DTN_FINDS_TESTS=ON
cmake --build build-nebula
ctest --test-dir build-nebula --output-on-failure
```

目标为 `tn_finds_resampler` 和 `tn_finds`。`tn_finds` 同时可作 resampler 和 wavtool。CMake 的其他工具链也可构建标准 F0 版。

默认产物只使用 SSE2 基线，任何 x86-64 CPU 都能运行。加上 `-DTN_FINDS_NEBULA_FAST_ISA=ON` 会用 AVX2/FMA 编译 Nebula 内核（需要 2013 年之后的 CPU），源 F0 分析更快。`-DTN_FINDS_NEBULA_FAST_COLD=OFF` 和 `-DTN_FINDS_NEBULA_FAST_KERNELS=OFF` 可分别退回未优化的粗估计和内核，便于对照；这两个选项只在源码修订与已适配版本一致时生效。

## Nebula 源 F0

启用 Nebula 时，源 F0 由 OpenMoresampler 的后端分析整条录音，再按 OTO 裁切位置重新对齐。搜索范围为 70–800 Hz，帧间隔约 5 ms。相比默认自相关，它主要减少八度误判和清音误报；代价是每次单音调用都会分析整条录音，长原音会增加耗时。

`--report` 的 `f0_backend` 字段为 `nebula` 或 `standard`，可用来确认当前 EXE 的后端。

### 源 F0 缓存

Nebula 构建在 `%LOCALAPPDATA%\tn_finds\nebula-f0-v1` 缓存整条源 F0，只存轨迹，不存音频。单音入口和双角色入口共用：同一段录音换音高、OTO 裁切、时长或 flags 时不必重算源 F0，即使宿主重新生成了单音缓存。缓存按解码后的采样值和采样率识别录音，并带上估计器与模型版本；录音内容一改就自动重算，也不会写入音源目录。缓存损坏或目录不可写时直接分析。

`TN_FINDS_NEBULA_CACHE=off` 可禁用这层缓存；设为绝对路径可把缓存放到别处隔离（相对路径会被忽略）。停止渲染后可以删掉缓存目录，下一次重建。换 EXE 后仍须清一次宿主的音符缓存。

### 分析并行

Nebula 构建对较长的有声区按谱脊帧并行分析，重叠窗仍按原顺序累加，短音节不会新开线程。可用 `OMP_NUM_THREADS` 限制线程数，或在构建时设 `-DTN_FINDS_DSP_OPENMP=OFF`。普通 CMake 构建默认不启用这部分 OpenMP。

## 经典模式

```text
tn_finds input output note velocity flags offset length consonant cutoff volume modulation [!tempo [pitchbend]]
```

```powershell
.\tn_finds_v0.4.9-resampler.exe in.wav out.wav A3 100 eTnr15g-15 0 900 100 0 100 0 "!120"
```

输出单声道 PCM16，采样率与输入相同。`consonant` 是原音设定的固定区；固定区输出长度乘 `2^(1-velocity/100)`。`cutoff` 为正时从末尾裁除，为负时从 `offset` 向后取 `abs(cutoff)` 毫秒。长度四舍五入到样本。

默认允许循环；`e` 关闭循环，改用线性时间映射。`modulation` 保留源 F0 的相对变化比例。音高弯曲使用 UTAU 的 12-bit Base64，`#n#` 重复上一个值 n 次，间隔为 1/96 拍。输出峰值先归一到 0.5，再乘音量。

| Flag | 行为 |
|---|---|
| `g` | 共振峰比例 `10^(-g/200)`。正值下移，`g-15.0515` 约上移 3 半音。 |
| `t` | 调音偏移，单位为 10 音分；不改共振峰。 |
| `b` | `0-100` 增强无声部分；`b10` 约 +5%，`b100` 约 +500%。 |
| `A` / `O` | 随局部音高变化的轻微音量调整 / 浊音亮度。 |
| `Tnr` | `0-100`，高频气声退回系数，默认 `Tnr75`。`Tnr0` 使用原始五周期平均。 |
| `Tns` / `Tnc` | `0-100`，把固定区之后 / 固定区内的残差替换为平滑频谱噪声。默认 `0` 保留录音残差；中间值按功率混合。 |
| `W` | 负值强制无声；`50-1000` 覆盖源 F0；`0` 或小于 50 的正值仍检测。 |
| `e` / `Me` | `e` 或 `e1` 关闭循环，使用线性伸缩；`Me` 或 `Me1` 优先尝试安全的稳态循环；`e0` / `Me0` 保持默认。同时出现时最后一个生效。 |
| `Mt` / `Mb` / `Mo` / `Mr` / `ME` | Moresampler 风格的松紧、呼吸度、开口度、共振和共振峰强度，范围 `-100..100`，默认 `0`。 |

`Mt` 从周期谐波估计声门参数 `Rd`，按非线性映射修改；`Mb` 缩放原录音分离出的气声并渐隐周期声；`Mo` 施加以目标基频为中心的频谱幅度曲线；`Mr` 调整 3 kHz 附近；`ME` 调整已有谱峰与谱谷的对比。`Mt0Mb0Mo0Mr0ME0` 与省略这些 flag 逐字节一致。

`e0/e1` 和 `Me0/Me1` 是为 OpenUtau 的 Options 表情准备的。OpenUtau 只拼接选项文本，不会自动加简称；选项值应写 `e0,e1`（或 `e,`），不要写 `0,1`，否则会得到 `...Mt00` / `...Mt01`，数字粘到前一个数值 flag 上。

## 独立模式

```text
tn_finds --render input.wav output.wav [options]
tn_finds --estimate-f0 input.wav track.csv [--no-mrq]
tn_finds --write-mrq input.wav | voicebank-directory [desc.mrq]
tn_finds --bind-mrq original.wav | voicebank-directory        (Windows)
```

默认线性伸缩、float32 WAV、不归一化。零效果、等长、音量 100 时走样本直通，不检测 F0、不合成。

| 选项 | 行为 |
|---|---|
| `--pitch-st N` | 相对音高，半音，默认 0。 |
| `--target-hz N` | 目标音高，0 保留源轮廓。 |
| `--formant-st N` | 共振峰移动，半音；正数上移。 |
| `--length-ms N` | 输出长度，-1 保持选中输入长度。 |
| `--offset-ms N` / `--cutoff-ms N` | 输入起点 / 末端裁除；负 cutoff 从 offset 向后取。 |
| `--consonant-ms N` / `--velocity N` | 固定区长度和辅音速度。 |
| `--mask-return N` | 气声退回系数 `0..1`，默认 0.75；`--harmonic-amount` 是别名。 |
| `--residual-mix N` / `--fixed-residual-mix N` | 对应 `Tns` / `Tnc`，范围 `0..100`。 |
| `--loop` / `--force-loop` | 允许循环 / 优先短稳态循环。 |
| `--tension N` 等 | 对应 `Mt`、`Mb`、`Mo`、`Mr`、`ME`。 |
| `--f0 TRACK.csv` | 源 F0 覆盖，绝对输入秒数；0 表示不提供周期提示。 |
| `--no-mrq` | 忽略本来可用的 MRQ；`TN_FINDS_MRQ=off` 对经典调用同样生效。 |
| `--write-mrq` / `--bind-mrq` | 生成 / 登记 MRQ 条目，见上文。 |
| `--volume N`、`--pcm16`、`--peak-normalize`、`--report` | 音量、输出格式、归一化和 JSON 计时报告。 |

`--estimate-f0` 输出源 F0 轨迹，零 Hz 表示清音。`--f0` 格式：

```csv
time_s,f0_hz
0.000,0
0.115,0
0.120,220.1
0.125,221.0
```

时间是完整原文件的秒数，即使设置了 `--offset-ms` 也不用重减。时间递增；首尾外侧沿用最近端点。正值是候选周期，0 表示没有提示，由音频相干度决定提取多少。

## MRQ 源 F0

经典调用、双角色音符准备和独立渲染默认读取输入原音同目录的 `desc.mrq`，按**完整文件名（含扩展名）**匹配条目。可用的 MRQ v1/v2 以自己的采样间隔和 Hz 为准，`0` 不当作“请重新检测”；只有读不到匹配项或条目无效时才走自动后端（Nebula 或默认估计器）。它只提供源 F0，谐波分离、残差、目标音高和整句相位仍由 tn_finds 计算。

`--f0` 和强制的 `W` / `--forced-hz` 覆盖默认选择；`--no-mrq` 单次关闭，`TN_FINDS_MRQ=off` 在该进程及其子进程关闭。`--estimate-f0 输入.wav 输出.csv` 导出实际选中的轨迹，加 `--no-mrq` 才强制导出自动后端的结果。独立渲染加 `--report` 时，`f0_backend` 为 `mrq` 表示命中。

**生成 MRQ。** `--write-mrq` 先测量录音，再按 Moresampler 的 v2 规则写入描述文件：

```powershell
& .\tn_finds.exe --write-mrq 'D:\Voicebanks\MyBank\F4\a.wav'   # 写入同目录 desc.mrq
& .\tn_finds.exe --write-mrq 'D:\Voicebanks\MyBank'            # 递归整库，逐目录写各自的 desc.mrq
& .\tn_finds.exe --write-mrq a.wav D:\cache\desc.mrq           # 指定目标文件
```

记录包含完整文件名（UTF-16LE）、`nf0` / `fs` / `nhop`、float32 轨迹、时间戳和 `modified=524`（表示工具生成、非手工修改）。同名条目放得下时原地重写，放不下时把旧名字清零（墓碑）并在末尾追加，条目计数加一；这与 Moresampler 的 `mrq_write_entry` 逐条对应，因此两边的文件可以互相读取。写入期间独占锁住版本字段，Moresampler 读到正在写的文件会跳过，而不是读半条记录。

只允许 Nebula 后端写 MRQ：普通后端的帧间隔（44.1 kHz 下 221 采样）和清浊判定与 Moresampler 不同，写出的文件会被别的重采样器当作 Nebula 结果使用，所以那种构建直接拒绝并说明原因，不改动任何文件。

**OpenUtau 外部引擎。** 宿主可能把源录音改名为 `src-….wav`。缓存目录里已有同名条目时可以直接读；若 MRQ 只在原声库，给原 WAV 或音库目录登记一次：

```powershell
& .\tn_finds.exe --bind-mrq 'D:\Voicebanks\MyBank'
```

用 OpenUtau 实际配置的原音路径（保留路径拼写）。登记只扫描指定文件或目录，索引存在 `%LOCALAPPDATA%\tn_finds\mrq-bindings-v1`，不改声库或宿主。匹配路径指纹且音频内容相同后直接读取原声库当前的条目，优先级高于临时目录内的旧同名条目。原音换过或搬过要重新登记；只改 F0 不必重新登记。切换引擎或修订 MRQ 后清一次宿主音符缓存。

## 双角色整句

`tn_finds` 可以同时担任 resampler 和 wavtool：

- **UTAU**：项目属性里 Tool1 和 Tool2 都指向 `tn_finds_v0.4.9-phrase.exe`，使用标准单线程批处理。
- **OpenUtau**：把同一份 EXE 放进 `Resamplers` 和 `Wavtools`，CLASSIC 轨道的两项都选它。

resampler 阶段只写目标 F0、谐波参数和录音残差；wavtool 收齐音符后，用一条共享的 F0/相位时钟合成整句。切换版本或工具组合后清理一次音符缓存。

标准 Windows 批处理会自动识别。手动调用：

```text
tn_finds --prepare-note ...        # 经典重采样参数
tn_finds --phrase ... LAST_NOTE    # 经典拼接参数，最后一次追加 LAST_NOTE
tn_finds --wavtool ...             # 普通 PCM 包络拼接
```

未识别的普通 wavtool 调用走 PCM 包络相加。

## 分析缓存

CMake 的 Windows 构建在 `%LOCALAPPDATA%\tn_finds\oto-analysis-v1` 缓存 OTO 实际分析区间的局部结果：时钟、源脉冲、谐波分量、残差修正和所需的声门 `Rd`。单音入口和双角色入口共用。键按解码采样、采样率、分析区间（包含短尾所需上下文）、源 F0 与后端、`Tnr` 和分析选项校验。只改目标音高或曲线、区间不变时不用再分析原音；改 OTO 或改变短尾上下文会按新范围重算。`Tnc` / `Tns` 的噪声生成、时间映射、峰值处理和目标相位仍按当前音符重新计算。双角色 wavtool 仍按原顺序快照音符参数，用同一个相位时钟合成整句，不拼接局部 WAV。`--report` 的 `analysis_cache_hit` 表示本次复用了这层分析。

`TN_FINDS_ANALYSIS_CACHE=off` 可禁用这层缓存；设为绝对路径可把它放到别处隔离。缓存包含源音频的 DSP 中间分量，44.1 kHz 下约每秒分析区间占 0.34 MiB，不写入音库；单条结果超过 16 MiB 不会缓存。目前没有自动清理，停止渲染后可删掉缓存目录；缓存异常或不可用时回退到原计算。要强制重做自动分析，同时关闭 MRQ、这层缓存和 `TN_FINDS_NEBULA_CACHE`；`--estimate-f0` 不经过这层缓存。

## 实现

- **F0**：默认在约 6.3 kHz 的带限视图上用 Hann 窗校正自相关选周期，再用局部相关排除邻近元音造成的误报，搜索约 55–1300 Hz。Nebula 构建改为整条录音的 70–800 Hz 模型轨迹；同目录的 `desc.mrq` 可用时优先采用它给出的轨迹。
- **谐波与气声**：用对齐的五周期估计中高频谐波，再用窄谱线软 mask 把不够周期的部分退回残差。低频按前后周期的一致性分离。
- **变调与共振峰**：谐波幅度来自源频谱包络，相位使用连续的目标准时钟；共振峰包络分别作用于谐波和浊音气声。
- **时间映射**：短窗 WSOLA 伸缩原音、原谐波和残差修正；长音只循环有实测周期支撑的稳定段，没有合适片段就使用线性映射。
- **残差**：`Tns` / `Tnc` 为 0 时保留录音残差；大于 0 时按源时间取平滑谱，用独立噪声重建，避免反复拼接形成固定嗡声。
- **双角色**：音符参数只准备一次，整句只合成一次，重叠区不重新估计 F0。

## 测试

需要 NumPy 和 SoundFile：

```powershell
python -B tests/test_audio.py --exe .\tn_finds_v0.4.9-resampler.exe
python -B tests/test_harmonic_flag.py --exe .\tn_finds_v0.4.9-resampler.exe
python -B tests/test_voice_flags.py --exe .\tn_finds_v0.4.9-resampler.exe
python -B tests/test_wavtool.py --exe .\tn_finds_v0.4.9-phrase.exe --resampler .\tn_finds_v0.4.9-resampler.exe
python -B tests/test_phrase.py --exe .\tn_finds_v0.4.9-phrase.exe
python -B tests/test_nebula_f0.py .\tn_finds_v0.4.9-resampler.exe
python -B tests/test_mrq.py .\tn_finds_v0.4.9-resampler.exe .\tn_finds_v0.4.9-phrase.exe
python -B tests/test_nebula_cache.py .\tn_finds_v0.4.9-resampler.exe
python -B tests/test_source_analysis_cache.py .\tn_finds_v0.4.9-resampler.exe .\tn_finds_v0.4.9-phrase.exe
```

CMake 构建可加 `-DTN_FINDS_TESTS=ON` 后运行 `ctest`，会执行同一组测试和 FFT 回归；Nebula 与缓存测试只在对应后端可用或 Windows 上注册。`build.ps1` 的直编产物不含 Nebula 后端和两层缓存，用它跑 `test_source_analysis_cache.py` 会提示跳过，MRQ 也只能读取。

## 限制

严重气声、吼腔、多声部、极短低音或极端移调可能误判 F0 或出现周期泄漏。Nebula 的搜索上限为 800 Hz，且每次分析整条录音；高音或长原音需要留意。共振峰位移是近似，不能保证每个峰严格到位；与变调或伸缩同时使用时不保证原谐波相位不变。兼容参数不等价于原版音色。

两层缓存只按当前设置复用结果，不做容量管理：长期反复渲染同一批长原音会让 `%LOCALAPPDATA%\tn_finds` 变大，删掉目录即可清空。MRQ 写入只支持 Nebula 后端构建；普通构建读 MRQ，但会拒绝写入。

## 许可与第三方

[GPL-3.0-only](LICENSE)。感谢 tn_fnds、WORLD/EFB-GW、world4utau 及历代 UTAU 工具作者。设计延续 TD-PSOLA 的基音同步读写思想；没有将旧引擎代码编入本仓库。

Nebula 后端来自 [OpenMoresampler](https://github.com/LingYi0110/OpenMoresampler/tree/15a4b85eaac4944a0a55cae50de96fb5a6f15344)（GPL-3.0-only），算法与模型来自 [Sleepwalking/Nebula](https://github.com/Sleepwalking/nebula)（GPL-3.0）。官方 Windows EXE 静态链接该后端；源码构建不包含这些文件，需按上面的 CMake 步骤指向本地 OpenMoresampler。
