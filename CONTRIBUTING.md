# Contributing to diffusion_sicm

## Development Setup

```bash
git clone https://github.com/zh-ga/Diffusion.git
cd Diffusion
python -m venv .venv
source .venv/bin/activate
pip install -e ".[fdm]"
pip install pytest
```

## Branching Strategy

- `main` -- stable, production-ready code
- `feature/xxx` -- new features
- `fix/xxx` -- bug fixes

## Pull Request Process

1. Ensure all tests pass: `pytest tests/ -v`
2. Keep PRs focused on a single concern
3. Squash merge is preferred
