# OpenMM-KCX Docker Image

OpenMM Docker image with support for KCX (N6-carboxylysine / carboxylated lysine) residue parameterization for AMBER ff19SB force field.

## Quick Setup

### 1. Create GitHub Repository

1. Go to [github.com/new](https://github.com/new)
2. Name it `openmm-kcx`
3. Make it **Public** (required for DockerHub)
4. Don't initialize with README

### 2. Add DockerHub Secrets

1. Create a DockerHub account at [hub.docker.com](https://hub.docker.com) if you don't have one
2. Generate an access token: DockerHub → Account Settings → Security → New Access Token
3. In your GitHub repo: Settings → Secrets and variables → Actions → New repository secret
   - `DOCKERHUB_USERNAME`: Your DockerHub username
   - `DOCKERHUB_TOKEN`: Your DockerHub access token

### 3. Push Code to GitHub

```bash
cd "c:\Python32\Coding Projects\Hydantoinase\Custom_OpenMM_KCX"
git init
git add .
git commit -m "Initial commit: OpenMM-KCX Docker image"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/openmm-kcx.git
git push -u origin main
```

### 4. Wait for Build

GitHub Actions will automatically build and push the image. Check progress at:
`https://github.com/YOUR_USERNAME/openmm-kcx/actions`

### 5. Deploy to Tamarind

Once the image is pushed to DockerHub, deploy to Tamarind:

```python
import requests

api_key = "your-tamarind-api-key"
headers = {'x-api-key': api_key}

params = {
    "modelName": "openmm-kcx",
    "dockerImage": "YOUR_DOCKERHUB_USERNAME/openmm-kcx:latest",
    "description": "OpenMM with KCX (carboxylated lysine) force field support"
}

response = requests.post(
    "https://app.tamarind.bio/api/deploy-model",
    headers=headers,
    json=params
)
print(response.text)
```

## Files

- `Dockerfile` - Docker build configuration
- `kcx.frcmod` - AMBER force field modifications for KCX
- `kcx.lib` - AMBER library file for KCX residue
- `.github/workflows/docker-build.yml` - GitHub Actions workflow

## Usage

The KCX parameters are located at `/app/kcx.frcmod` and `/app/kcx.lib` inside the container.

To use in tleap:
```
loadamberparams /app/kcx.frcmod
loadoff /app/kcx.lib
```
