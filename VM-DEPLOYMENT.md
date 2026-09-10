# Ubuntu VM deployment

The app binds to `127.0.0.1:6875` on the VM. Nginx can publish that private upstream through a public HTTPS domain. Uploaded files and temporary session data remain inside the container or its temporary filesystem. Project, table, and figure downloads are returned through the user's browser.

Replace `your-vm.example.org` and `your-user` below with values for your VM.

## Install prerequisites

Install Docker Engine, the Docker Compose plugin, nginx, and Certbot using the instructions for your Ubuntu release. Confirm that your domain resolves to the VM and that inbound TCP ports 80 and 443 are allowed.

## Install the app

```bash
git clone https://github.com/mwfoster/proteomics-data-workup.git
cd ~/proteomics-data-workup
docker compose up -d --build
docker compose ps
```

Wait until the container reports `healthy`. The Compose configuration exposes the Shiny app only at `http://127.0.0.1:6875` on the VM.

## Configure nginx and HTTPS

Create `/etc/nginx/sites-available/proteomics-data-workup`:

```nginx
server {
    listen 80;
    server_name your-vm.example.org;

    client_max_body_size 1024M;
    client_body_timeout 600s;

    location / {
        proxy_pass http://127.0.0.1:6875;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_request_buffering off;
        proxy_read_timeout 600s;
        proxy_send_timeout 600s;
    }
}
```

Enable the site and request a certificate:

```bash
sudo ln -s /etc/nginx/sites-available/proteomics-data-workup /etc/nginx/sites-enabled/proteomics-data-workup
sudo nginx -t
sudo systemctl reload nginx
sudo certbot --nginx -d your-vm.example.org
```

Open `https://your-vm.example.org/`. Nginx owns public ports 80 and 443; the Shiny port remains private.

## View logs

```bash
cd ~/proteomics-data-workup
docker compose logs --tail=100 -f
```

## Update later

```bash
cd ~/proteomics-data-workup
git pull --ff-only
docker compose up -d --build
```

## Stop the app

```bash
cd ~/proteomics-data-workup
docker compose down
```
