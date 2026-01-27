# Installation
## snap VS systemd

这是两种安装方式.由于systemd安装更自由，docker源设置更方便，加上目前网上snap debug教程极少，强烈推荐systemd安装，基本不考虑snap安装
tutorial: https://zhuanlan.zhihu.com/p/1924864194345431337

|特性	|Snap 安装| Systemd 安装（apt/docker官方源）|
|-|-|-|
|安装方式|	sudo snap install docker | sudo apt install docker-ce|
|管理方式|	Snap 系统管理 | systemd 服务管理|
|权限模型|	严格的沙箱限制 | 传统的 Linux 权限|
|隔离性	|高度隔离（沙箱）	|较少隔离|
|更新方式 |自动更新 | 手动或通过包管理器|
|文件位置 | /snap/docker/current/ | /usr/bin/docker, /var/lib/docker/|
|性能|可能有轻微开销 | 原生性能|

## 镜像源推荐

```sh
{
    "registry-mirrors": [
     "https://hub.rat.dev",
     "https://docker.unsee.tech/",
     "https://docker.anye.in",
     "https://hub.geekery.cn",
     "https://dockerpull.org",
     "https://docker.1panel.live"
  ]
}
```


# set and check user permission
```bash
# 1. 将用户添加到docker组（只需要做一次）
sudo usermod -aG docker $USER

# 2. 重新登录或刷新组权限
newgrp docker

# 3. 验证不需要sudo就能用docker
docker ps
```