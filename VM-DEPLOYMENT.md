# Ubuntu VM deployment

The app binds to `127.0.0.1:6875` on the VM. The VM's existing nginx and Certbot configuration publishes that private upstream at `https://rapid-2905.vm.duke.edu/`. Uploaded files and temporary session data remain inside the container or its temporary filesystem. Project and table downloads are returned through the user's browser.

## Upload from Windows PowerShell

```powershell
scp .\proteomics-data-workup-vm.tar.gz rapiduser@rapid-2905.vm.duke.edu:~/
ssh rapiduser@rapid-2905.vm.duke.edu
```

## Install on the VM

```bash
mkdir -p ~/proteomics-data-workup
tar -xzf ~/proteomics-data-workup-vm.tar.gz -C ~/proteomics-data-workup
cd ~/proteomics-data-workup
docker compose up -d --build
docker compose ps
```

Open `https://rapid-2905.vm.duke.edu/` after the container becomes healthy. Nginx owns public ports 80 and 443 and proxies to `http://127.0.0.1:6875`; the Shiny port is not exposed publicly.

## View logs

```bash
cd ~/proteomics-data-workup
docker compose logs --tail=100 -f
```

## Update later

Upload a new archive, extract it into the same directory, and rebuild:

```bash
cd ~/proteomics-data-workup
docker compose up -d --build
```
