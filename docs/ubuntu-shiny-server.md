# Ubuntu 24.04 Shiny Server Deployment

These notes describe a straightforward VM deployment for the Proteomics Data Workup Shiny app.

## Recommended VM

- Ubuntu 24.04 LTS
- 4 CPU
- 40 GB RAM
- 200 GB disk

## Install System Packages

```bash
sudo apt update
sudo apt install -y build-essential git gdebi-core libcurl4-openssl-dev libssl-dev libxml2-dev \
  libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
  pandoc r-base r-base-dev
```

## Install Shiny Server

Check Posit's download page for the current Shiny Server package URL. A typical flow is:

```bash
wget https://download3.rstudio.org/ubuntu-20.04/x86_64/shiny-server-1.5.23.1030-amd64.deb
sudo gdebi shiny-server-1.5.23.1030-amd64.deb
```

If the exact version changes, use the current Ubuntu/Debian `.deb` package from Posit.

## Deploy the App

```bash
sudo mkdir -p /srv/shiny-server/proteomics-data-workup
sudo chown -R $USER:$USER /srv/shiny-server/proteomics-data-workup
cd /srv/shiny-server

git clone https://github.com/mwfoster/proteomics-data-workup.git proteomics-data-workup
cd proteomics-data-workup
Rscript install_packages.R
```

If the repository is private, use SSH deploy keys or GitHub's HTTPS authentication according to your institution's policy.

## Test Locally on the VM

```bash
cd /srv/shiny-server/proteomics-data-workup
Rscript run_app.R
```

Then browse to:

```text
http://VM_HOSTNAME_OR_IP:6875/
```

For Shiny Server's default hosting, browse to:

```text
http://VM_HOSTNAME_OR_IP:3838/proteomics-data-workup/
```

## Update the App

```bash
cd /srv/shiny-server/proteomics-data-workup
git pull
Rscript install_packages.R
sudo systemctl restart shiny-server
```

## Data Policy

Do not store private study uploads in the Git repository. Users should upload data through the Shiny app or use approved internal storage. The app's project ZIP export can include study files, so treat exported ZIPs as sensitive data.

## Troubleshooting

View Shiny Server logs:

```bash
sudo journalctl -u shiny-server --no-pager -n 200
ls -lh /var/log/shiny-server/
```

Check R package install issues by running:

```bash
cd /srv/shiny-server/proteomics-data-workup
Rscript install_packages.R
```
