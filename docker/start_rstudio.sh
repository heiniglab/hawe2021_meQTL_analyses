#!/bin/bash

# generate a new password to be used for the current user
password=$(openssl rand -base64 20)
export RSTUDIO_PASSWORD=${password}

# port is the first parameter, must not be empty
port=${1}
if [[ -z ${port} ]]; then
  echo "You must specify a port. Suggested range: 8000-9000"
  exit 1
fi

echo "Password for this session is:"
echo ${password}

# this path has to match the one in the dockerfile!
RSTUDIO_AUTH="/bin/rstudio_light_auth"

# create a user-specififc session cookie file
mkdir -p /localscratch/${USER}/
COOKIE_FILE=/localscratch/${USER}/rstudio-session-cookie
test -f ${COOKIE_FILE} || echo `uuidgen` > ${COOKIE_FILE}

echo ""
echo "Running RStudio server at http://${SLURMD_NODENAME}:${port}"

# run rstudio
/usr/lib/rstudio-server/bin/rserver \
  --www-port=${port} \
  --auth-none=0 \
  --auth-pam-helper-path=${RSTUDIO_AUTH} \
  --auth-encrypt-password=0 \
  --secure-cookie-key-file ${COOKIE_FILE}

