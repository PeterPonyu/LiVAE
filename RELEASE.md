# Release Process

This document describes the automated and manual release process for LiVAE.

## Automated PyPI Publishing

### How It Works

When you create a new GitHub release, the package is **automatically published to PyPI** via GitHub Actions.

**Workflow:** `.github/workflows/publish.yml`

**Trigger:** When a new release is published on GitHub

**Steps:**
1. Checks out the code
2. Verifies the git tag version matches `livae.__version__`
3. Builds wheel and sdist distributions
4. Validates distributions with twine
5. Uploads to PyPI using `PYPI_API_TOKEN` secret

### Prerequisites

✅ **PYPI_API_TOKEN** secret is configured in repository settings  
✅ Version in `livae/__init__.py` matches the release tag  
✅ Version in `pyproject.toml` matches the release tag

## Release Checklist

### 1. Update Version

Update version in **both** files:

```python
# livae/__init__.py
__version__ = '0.2.2'  # Update this
```

```toml
# pyproject.toml
[project]
version = "0.2.2"  # Update this
```

### 2. Update CHANGELOG (Optional but Recommended)

Document changes in `CHANGELOG.md` or release notes.

### 3. Run Tests Locally

```bash
pytest -q
```

### 4. Commit and Push

```bash
git add livae/__init__.py pyproject.toml
git commit -m "chore: bump version to 0.2.2"
git push origin main
```

### 5. Create GitHub Release

**Option A: Using GitHub CLI**

```bash
# Build distributions locally (optional, for release assets)
python -m build

# Create release
gh release create v0.2.2 \
  dist/livae-0.2.2.tar.gz \
  dist/livae-0.2.2-py3-none-any.whl \
  --title "v0.2.2: Release Title" \
  --notes "Release notes here"
```

**Option B: Using GitHub Web Interface**

1. Go to https://github.com/PeterPonyu/LiVAE/releases/new
2. Choose tag: `v0.2.2` (must match package version)
3. Set title: `v0.2.2: Release Title`
4. Add release notes
5. Click "Publish release"

### 6. Automated PyPI Publication

Once the release is published:
- GitHub Actions workflow triggers automatically
- Package is built and published to PyPI
- Check workflow status: https://github.com/PeterPonyu/LiVAE/actions

### 7. Verify Publication

```bash
# Wait ~2 minutes for PyPI to index
pip install --upgrade livae

# Verify version
python -c "import livae; print(livae.__version__)"
```

**PyPI Package:** https://pypi.org/project/livae/

## Manual PyPI Publishing (Fallback)

If automated publishing fails, you can publish manually:

```bash
# Build distributions
python -m build

# Check distributions
python -m twine check dist/*

# Upload to PyPI
python -m twine upload dist/*
# Username: __token__
# Password: [PYPI_API_TOKEN value]
```

## Version Numbering

Follow [Semantic Versioning](https://semver.org/):

- **MAJOR** (0.x.0): Breaking changes
- **MINOR** (0.x.0): New features, backward compatible
- **PATCH** (0.0.x): Bug fixes, backward compatible

## Troubleshooting

### Version Mismatch Error

```
Error: Tag version (0.2.2) does not match package version (0.2.1)
```

**Solution:** Update `livae/__init__.py` and `pyproject.toml` to match the tag.

### PyPI Token Invalid

**Solution:** Regenerate token at https://pypi.org/manage/account/token/ and update `PYPI_API_TOKEN` secret:

```bash
gh secret set PYPI_API_TOKEN --repo PeterPonyu/LiVAE
# Paste token when prompted
```

### Build Fails

**Solution:** Ensure all dependencies in `pyproject.toml` are correct and test locally:

```bash
pip install build
python -m build
```

## Current Release

**Latest Version:** v0.2.1  
**PyPI:** https://pypi.org/project/livae/0.2.1/  
**GitHub:** https://github.com/PeterPonyu/LiVAE/releases/tag/v0.2.1
