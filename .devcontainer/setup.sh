pip install -e .[dev]
pip install -e BEERS_UTILS
pip install -e CAMPAREE

echo 'export PYTHONDONTWRITEBYTECODE=1' >> ~/.bashrc

ssh-keyscan -t rsa github.com >> ~/.ssh/known_hosts
