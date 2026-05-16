#!/usr/bin/env python3
"""
显式构建或检查 workflow family images。

示例：
  python scripts/bootstrap_env_images.py
  python scripts/bootstrap_env_images.py --profile py-scverse-v1
  python scripts/bootstrap_env_images.py --check-only
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.utils.execution_environment import image_catalog_root, load_image_catalog


def _run(cmd: list[str], *, cwd: Path | None = None) -> None:
    subprocess.run(cmd, check=True, env=os.environ.copy(), cwd=str(cwd) if cwd else None)


def _capture(cmd: list[str], *, cwd: Path | None = None) -> str:
    result = subprocess.run(
        cmd,
        check=True,
        env=os.environ.copy(),
        cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    return result.stdout


def _image_exists(tag: str) -> bool:
    cmd = ["docker", "image", "inspect", tag]
    result = subprocess.run(cmd, env=os.environ.copy(), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return result.returncode == 0


def _select_shell(runtime: str) -> str:
    # 当前 catalog 使用的基础镜像都应具备 /bin/sh，统一使用 sh -lc 可减少跨镜像差异。
    if runtime in {"python", "r", "mixed"}:
        return "sh"
    return "sh"


def _run_smoke(entry) -> None:
    if not entry.smoke_test:
        print("  smoke: (skip, no smoke_test)")
        return
    shell = _select_shell(entry.runtime)
    cmd = [
        "docker",
        "run",
        "--rm",
        entry.image_tag,
        shell,
        "-lc",
        entry.smoke_test,
    ]
    print(f"  smoke start: {' '.join(cmd)}")
    _run(cmd)
    print("  smoke done")


def _build_image(entry, root: Path) -> None:
    dockerfile_path = root / entry.dockerfile
    context_path = (root / entry.build_context).resolve()
    if not dockerfile_path.is_file():
        raise FileNotFoundError(f"dockerfile 不存在: {dockerfile_path}")

    cmd = [
        "docker",
        "build",
        "-t",
        entry.image_tag,
        "-f",
        str(dockerfile_path),
        str(context_path),
    ]
    print(f"  build start: {' '.join(cmd)}")
    _run(cmd)
    print("  build done")


def _python_package_versions(entry) -> str:
    cmd = [
        "docker",
        "run",
        "--rm",
        entry.image_tag,
        "sh",
        "-lc",
        "python -m pip freeze --all | sort",
    ]
    return _capture(cmd).strip()


def _r_package_versions(entry) -> str:
    cmd = [
        "docker",
        "run",
        "--rm",
        entry.image_tag,
        "sh",
        "-lc",
        (
            "Rscript -e \"ip <- installed.packages()[, c('Package', 'Version')]; "
            "ip <- ip[order(ip[,1]), , drop=FALSE]; "
            "cat(paste(ip[,1], ip[,2], sep='=='), sep='\\n')\""
        ),
    ]
    return _capture(cmd).strip()


def _capture_package_versions(entry) -> str:
    runtime = (entry.runtime or "").strip()
    if runtime == "python":
        return _python_package_versions(entry)
    if runtime == "r":
        return _r_package_versions(entry)
    if runtime == "mixed":
        py = _python_package_versions(entry)
        r = _r_package_versions(entry)
        return f"[python]\n{py}\n\n[r]\n{r}".strip()
    return ""


def _write_package_versions(entry) -> None:
    relpath = getattr(entry, "package_versions_path", None) or ""
    if not relpath:
        print("  package versions: (skip, no package_versions_path)")
        return
    content = _capture_package_versions(entry)
    if not content:
        print("  package versions: (skip, empty capture)")
        return
    output_path = (image_catalog_root() / relpath).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(content + "\n", encoding="utf-8")
    print(f"  package versions written: {output_path}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", help="仅处理单个 env_profile")
    parser.add_argument("--check-only", action="store_true", help="只打印 catalog，不执行 docker build")
    parser.add_argument("--smoke-only", action="store_true", help="只运行 smoke test，不执行 docker build")
    parser.add_argument("--write-package-versions-only", action="store_true", help="仅从现有镜像导出包版本清单，不执行 docker build")
    parser.add_argument("--skip-smoke", action="store_true", help="构建后跳过 smoke test")
    args = parser.parse_args()

    if args.check_only and (args.smoke_only or args.write_package_versions_only or args.skip_smoke):
        parser.error("--check-only 不能与 --smoke-only / --write-package-versions-only / --skip-smoke 同时使用")
    if args.smoke_only and args.skip_smoke:
        parser.error("--smoke-only 与 --skip-smoke 不能同时使用")
    if args.smoke_only and args.write_package_versions_only:
        parser.error("--smoke-only 与 --write-package-versions-only 不能同时使用")
    if args.write_package_versions_only and args.skip_smoke:
        parser.error("--write-package-versions-only 与 --skip-smoke 不能同时使用")

    catalog = load_image_catalog()
    if not catalog:
        print("image catalog 为空")
        return 1

    root = image_catalog_root()
    if args.profile:
        if args.profile not in catalog:
            raise KeyError(f"未知 env_profile: {args.profile}")
        entries = [catalog[args.profile]]
    else:
        entries = list(catalog.values())

    for entry in entries:
        print(f"[IMAGE] profile={entry.env_profile} tag={entry.image_tag} runtime={entry.runtime}")
        print(f"  dockerfile={entry.dockerfile}")
        print(f"  smoke_test={entry.smoke_test}")
        if args.check_only:
            continue
        if args.smoke_only:
            if not _image_exists(entry.image_tag):
                raise RuntimeError(f"镜像不存在，无法仅执行 smoke test: {entry.image_tag}")
            _run_smoke(entry)
            continue
        if args.write_package_versions_only:
            if not _image_exists(entry.image_tag):
                raise RuntimeError(f"镜像不存在，无法导出包版本清单: {entry.image_tag}")
            _write_package_versions(entry)
            continue

        _build_image(entry, root)
        if not args.skip_smoke:
            _run_smoke(entry)
        _write_package_versions(entry)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
