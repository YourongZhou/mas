This directory stores per-profile package-version exports captured from built Docker images.

Current convention:
- `*.pip-freeze.txt`: Python environments captured with `python -m pip freeze --all | sort`
- `*.installed-packages.txt`: R environments captured from `installed.packages()`
- `*.packages.txt`: mixed environments that include both Python and R sections

Regenerate these files with:

```bash
python mas_2/scripts/bootstrap_env_images.py
```
