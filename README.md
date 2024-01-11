# testSpectrogram

<table>
<thead>
  <tr>
    <th><a href="https://ieeexplore.ieee.org/abstract/document/10233699">SDP3 (IOTJ)</a></th>
    <th><a href="https://ieeexplore.ieee.org/abstract/document/9593198">PBAH (SPAWC)</a></th>
    <th><a href="https://arxiv.org/pdf/2311.07169.pdf">CASTER (arxiv preprint)</a></th>
    <th><a href="https://www.bilibili.com/video/BV14G411y7nn/?spm_id_from=333.999.0.0&vd_source=acf5c2e5837a698024101aaf6bf0d161">Video demo on bilibili (chinese)</a></th>
  </tr>
</thead>
<tbody>
  <tr>
    <td colspan="2"><img src="./README.assets/DAHC.png"></td>
    <td><img src="./README.assets/channel_simulation.png"></td>
    <td><a href="https://www.bilibili.com/video/BV14G411y7nn/?spm_id_from=333.999.0.0&vd_source=acf5c2e5837a698024101aaf6bf0d161"><img src="./README.assets/image-20231116131857577.png"></a></td>
  </tr>
</tbody>
</table>

## Introduction

[testSpectrogram](https://github.com/rzy0901/testSpectrogram) is an open-source platform for wireless channel simulation, human/hand pose extraction, gesture spectrogram generation, and real-time gesture recognition based on millimeter-wave passive sensing and communication systems.

Detailed documentation, source codes, and paper information will be coming soon (maybe after paper acceptance).

## Codes Overview

### Codes for "SDP3" and "PBAH" papers

+ [Micro_Doppler_Radar_Simulator](./Micro_Doppler_Radar_Simulator)
  + Data driven hybrid channel model simulation using a Boulic Human walking model.
+ [testZED](https://github.com/rzy0901/testZED) and [zed_pose](./zed_pose)
  + Simple Mocap-based channel simulation example.
  + Camera coordinate 3D human keypoints extraction based on the depth camera ZED 2i, using [zed-sdk](https://github.com/stereolabs/zed-sdk).

### Codes for "CASTER" paper

+ [mediapipe_spectrogram](https://github.com/rzy0901/mediapipe_spectrogram) (Not open-sourced yet)
  + Primitive-based wireless channel simulation for hand gesture recognition.
  + Camera coordinate 3D hand keypoints extraction based on a monocular camera, using [mediapipe](https://github.com/google/mediapipe) and [opencv](https://github.com/opencv/opencv).
+ [mediapipe_spectrogram_classification](https://github.com/Jcq242818/CASTER/tree/retraining2) (Not open-sourced yet) and [RxRealTime_GUI_rzy](https://github.com/rzy0901/RxRealTime_GUI_rzy/tree/master)
  + "Simulation-to-reality" hand gesture recognition based on ResNet18.
  + Real-time gesture recognition based on millimeter-wave passive sensing and communication systems, using a model trained by a simulated dataset.

## Cite this repository

```
@inproceedings{li2021wireless,
  title={Wireless sensing with deep spectrogram network and primitive based autoregressive hybrid channel model},
  author={Li, Guoliang and Wang, Shuai and Li, Jie and Wang, Rui and Peng, Xiaohui and Han, Tony Xiao},
  booktitle={2021 IEEE 22nd International Workshop on Signal Processing Advances in Wireless Communications (SPAWC)},
  pages={481--485},
  year={2021},
  organization={IEEE}
}
```

```
@article{li2023integrated,
  title={Integrated Sensing and Communication from Learning Perspective: An SDP3 Approach},
  author={Li, Guoliang and Wang, Shuai and Li, Jie and Wang, Rui and Liu, Fan and Peng, Xiaohui and Han, Tony Xiao and Xu, Chengzhong},
  journal={IEEE Internet of Things Journal},
  year={2023},
  publisher={IEEE}
}
```

```
@misc{ren2023caster,
      title={CASTER: A Computer-Vision-Assisted Wireless Channel Simulator for Gesture Recognition}, 
      author={Zhenyu Ren and Guoliang Li and Chenqing Ji and Chao Yu and Shuai Wang and Rui Wang},
      year={2023},
      eprint={2311.07169},
      archivePrefix={arXiv},
      primaryClass={eess.SP}
}
```

## Authors

[Zhenyu Ren](https://github.com/rzy0901)

[Guoliang Li](https://github.com/GuoliangLI1998)

[Shuai Wang](https://github.com/bearswang)

[Chenqing Ji](https://github.com/Jcq242818)

[Chao Yu](https://github.com/Ychao12032212)

## Acknowledgement

This series of work is under supervision of Prof. Rui Wang and Prof. Shuai Wang.
