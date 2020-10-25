# Amazon Linux 2

# update
sudo yum update

# git
sudo yum install git

# git-lfs: https://github.com/git-lfs/git-lfs/blob/master/INSTALLING.md
curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.rpm.sh | sudo bash
sudo yum install git-lfs

# Docker: https://docs.aws.amazon.com/AmazonECS/latest/developerguide/docker-basics.html
sudo amazon-linux-extras install docker
sudo service docker start
sudo usermod -a -G docker ec2-user
# Re-Login or Re-Boot

# Docker-Compose: https://docs.docker.com/compose/install/
sudo curl -L "https://github.com/docker/compose/releases/download/1.27.4/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
sudo chmod +x /usr/local/bin/docker-compose

# Clone repo
git clone https://git.rwth-aachen.de/lamBOO/fenicsR13.git

# Start Docker Environment
chmod -R 777 fenicsR13/
cd fenicsR13
docker-compose run fenicsr13_debug
